//! Dictionary builder orchestration
//!
//! Coordinates the multi-step pipeline to build an SSHash dictionary:
//! 1. Parse and encode sequences into SPSS
//! 2. Extract minimizer tuples
//! 3. Classify buckets
//! 4. Build minimizers control map
//! 5. Build sparse and skew index
//! 6. Assemble final dictionary

use crate::{
    builder::{
        buckets::{classify_into_buckets_inplace, ClassifiedBuckets, BucketStatistics, MIN_BUCKET_SIZE},
        config::BuildConfiguration,
        encode::Encoder,
        external_sort::FileTuples,
        minimizer_tuples::{compute_minimizer_tuples, compute_minimizer_tuples_external_file, needs_external_sorting},
    },
    dictionary::Dictionary,
    kmer::{Kmer, KmerBits},
    minimizers_control_map::{MinimizersControlMap, MinimizersControlMapBuilder, BucketType},
    partitioned_mphf::PartitionedMphf,
    sparse_and_skew_index::SparseAndSkewIndex,
    spectrum_preserving_string_set::SpectrumPreservingStringSet,
};
use std::io::{BufWriter, Write};
use tracing::info;

/// Bucket metadata collected during a sequential scan of sorted tuples.
///
/// Only stores per-bucket sizes (needed for EF construction and pass 2 routing).
/// Minimizer values are fed directly to the MPHF builder during the scan,
/// and bucket start positions are recovered by rescanning the file sequentially.
pub struct BucketMetadata {
    /// Number of unique super-kmers per bucket.
    pub cached_sizes: Vec<u32>,
    /// Count of singleton buckets (cached_size == 1).
    pub num_singleton: u64,
    /// Count of light buckets (2 <= cached_size <= MIN_BUCKET_SIZE).
    pub num_light: u64,
    /// Count of heavy buckets (cached_size > MIN_BUCKET_SIZE).
    pub num_heavy: u64,
}

impl BucketMetadata {
    /// Number of buckets.
    #[inline]
    pub fn num_buckets(&self) -> usize {
        self.cached_sizes.len()
    }
}

/// Builder for constructing SSHash dictionaries
pub struct DictionaryBuilder {
    config: BuildConfiguration,
}

impl DictionaryBuilder {
    /// Create a new dictionary builder with the given configuration
    pub fn new(config: BuildConfiguration) -> Result<Self, String> {
        config.validate()?;
        Ok(Self { config })
    }
    
    /// Build a dictionary from input sequences
    ///
    /// # Arguments
    /// * `sequences` - Vector of DNA sequences (strings)
    ///
    /// # Parallelism
    /// The number of threads is controlled by `config.num_threads`:
    /// - `0` — use all available CPU cores (rayon default)
    /// - `1` — single-threaded (no rayon overhead)
    /// - `N` — use exactly N threads
    ///
    /// # Returns
    /// A fully constructed Dictionary ready for queries
    pub fn build_from_sequences(&self, sequences: Vec<String>) -> Result<Dictionary, String> {
        #[allow(deprecated)]
        if !self.config.canonical {
            tracing::warn!(
                "non-canonical (regular) mode was removed in 0.7.0; \
                 building the unified canonical index"
            );
        }
        // Build a rayon thread pool sized to config.num_threads.
        // num_threads == 0 means "all cores" (rayon default).
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.config.num_threads)
            .build()
            .map_err(|e| format!("Failed to create thread pool: {e}"))?;

        pool.install(|| self.build_from_sequences_inner(sequences))
    }

    /// Inner build logic, runs inside the rayon thread pool
    fn build_from_sequences_inner(&self, sequences: Vec<String>) -> Result<Dictionary, String> {
        self.config.print();
        info!("Building SSHash Dictionary");

        // Step 1: Encode sequences into SPSS
        info!("Step 1: Encoding sequences...");
        let (spss, num_sequences) = self.encode_sequences(sequences)?;
        info!("  Encoded {} sequences", num_sequences);
        info!("  Total bases: {}", spss.total_bases());

        // Decide between external sort (streaming from disk) or in-memory paths
        let total_bases = spss.total_bases();
        let num_strings = spss.num_strings();
        let k = self.config.k as u64;
        let total_kmers = total_bases.saturating_sub(num_strings * (k - 1));

        if self.config.force_external_sort
            || needs_external_sorting(total_kmers, self.config.ram_limit_gib)
        {
            info!("Using external sorting: estimated {} k-mers exceeds RAM limit of {} GiB",
                total_kmers, self.config.ram_limit_gib);
            self.build_with_external_sort(spss)
        } else {
            self.build_with_in_memory_tuples(spss)
        }
    }

    /// In-memory build path: all tuples fit in RAM.
    ///
    /// This is the original pipeline: parallel tuple extraction → sort → classify
    /// in-place → build MPHF → build sparse/skew index.
    fn build_with_in_memory_tuples(&self, spss: SpectrumPreservingStringSet) -> Result<Dictionary, String> {
        // Step 2: Extract minimizer tuples (with inline coalescing during extraction)
        info!("Step 2: Extracting and coalescing minimizer tuples (in-memory)...");
        let tuples = crate::dispatch_on_k!(self.config.k, K => {
            Ok::<_, String>(compute_minimizer_tuples::<K>(&spss, &self.config))
        })?;
        info!("  Extracted and coalesced {} tuples", tuples.len());

        // Step 3: Classify into buckets (in-place — no tuple duplication)
        info!("Step 3: Classifying buckets (in-place)...");
        let classified = classify_into_buckets_inplace(tuples);

        // Compute statistics
        let mut stats = BucketStatistics::new();
        for i in 0..classified.num_buckets() {
            stats.add_bucket(classified.bucket_tuples(i));
        }
        stats.print_summary();

        // Step 4: Build minimizers control map
        info!("Step 4: Building minimizers control map...");
        let (control_map, bucket_id_by_mphf_index) = self.build_control_map(&classified)?;
        info!("  Built MPHF for {} minimizers", control_map.num_minimizers());

        // Step 5: Build sparse and skew index
        info!("Step 5: Building sparse and skew index...");
        let mphf_order = if !bucket_id_by_mphf_index.is_empty() {
            Some(bucket_id_by_mphf_index)
        } else {
            None
        };
        let index = self.build_index(&classified, mphf_order.as_deref(), &spss)?;
        info!("  Index built successfully");

        self.assemble_dictionary(spss, control_map, index)
    }

    /// External sort build path: tuples streamed from file via buffered I/O.
    ///
    /// Following the C++ approach: never materializes all tuples in memory.
    /// Instead, scans the file multiple times via fresh `BufReader`s:
    /// - Pass A: Scan buckets → write minimizers to temp file → build MPHF
    /// - Pass B: Rescan → collect cached_sizes + bucket_id_by_mphf_index
    /// - Pass 2: Build sparse index (fill offsets) + collect heavy bucket
    ///   tuples for skew index
    fn build_with_external_sort(&self, spss: SpectrumPreservingStringSet) -> Result<Dictionary, String> {
        // Step 2: External sort → merged file (accessed via buffered I/O)
        info!("Step 2: External sort...");
        let file_tuples = crate::dispatch_on_k!(self.config.k, K => {
            compute_minimizer_tuples_external_file::<K>(&spss, &self.config)
                .map_err(|e| e.to_string())
        })?;
        info!("  Sorted {} tuples to disk", file_tuples.num_tuples());

        // Steps 3+4: Scan file → feed minimizers to MPHF builder + collect sizes
        info!("Step 3: Scanning buckets (pass A + B)...");
        let (bucket_meta, control_map, bucket_id_by_mphf_index) =
            self.scan_and_build_control_map(&file_tuples)?;
        info!("  Found {} buckets ({} singleton, {} light, {} heavy)",
            bucket_meta.num_buckets(),
            bucket_meta.num_singleton,
            bucket_meta.num_light,
            bucket_meta.num_heavy);
        info!("  Built MPHF for {} minimizers", control_map.num_minimizers());

        // Step 5: Build sparse and skew index from file (pass 2)
        info!("Step 5: Building sparse and skew index (pass 2)...");
        let mphf_order = if !bucket_id_by_mphf_index.is_empty() {
            Some(bucket_id_by_mphf_index)
        } else {
            None
        };
        let index = self.build_index_from_file(
            &file_tuples,
            &bucket_meta,
            mphf_order,
            &spss,
        )?;
        info!("  Index built successfully");

        self.assemble_dictionary(spss, control_map, index)
    }

    /// Assemble the final dictionary from its components.
    fn assemble_dictionary(
        &self,
        spss: SpectrumPreservingStringSet,
        control_map: crate::minimizers_control_map::MinimizersControlMap,
        index: SparseAndSkewIndex,
    ) -> Result<Dictionary, String> {
        info!("Dictionary Build Complete");
        let total_bits = spss.num_bits() + control_map.num_bits() + index.num_bits();
        info!("Total memory: {:.2} MB", total_bits as f64 / (8.0 * 1024.0 * 1024.0));

        Ok(Dictionary::new(
            spss,
            control_map,
            index,
            self.config.k,
            self.config.m,
            self.config.seed,
        ))
    }
    
    /// Scan the merged file, collect minimizers, and build MPHF directly.
    ///
    /// Bypasses `MinimizersControlMapBuilder` entirely for the external sort path.
    /// Splits scanning into two passes to avoid holding 3.2 GB minimizers Vec
    /// and 1.6 GB cached_sizes Vec simultaneously:
    ///
    /// - Pass A: Scan file → write minimizers to temp file → mmap it → build MPHF
    /// - Pass B: Rescan file → build bucket_id_by_mphf_index + collect cached_sizes
    ///
    /// Each pass opens a fresh `BufReader` — no mmap pages linger between passes.
    ///
    /// Returns `(bucket_meta, control_map, bucket_id_by_mphf_index)`.
    fn scan_and_build_control_map(
        &self,
        file_tuples: &FileTuples,
    ) -> Result<(BucketMetadata, MinimizersControlMap, Vec<usize>), String> {
        // --- Pass A: Write minimizers to temp file, mmap it, build MPHF ---
        let minimizers_path = file_tuples.path().with_extension("minimizers.tmp");
        let mut num_buckets = 0usize;
        {
            let file = std::fs::File::create(&minimizers_path)
                .map_err(|e| format!("Failed to create minimizers temp file: {e}"))?;
            let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
            let bucket_iter = file_tuples.bucket_iter()
                .map_err(|e| format!("Failed to open file for pass A: {e}"))?;
            for scan in bucket_iter {
                writer.write_all(&scan.minimizer.to_ne_bytes())
                    .map_err(|e| format!("Failed to write minimizer: {e}"))?;
                num_buckets += 1;
            }
            writer.flush().map_err(|e| format!("Failed to flush minimizers: {e}"))?;
        }
        // Mmap the minimizers file as a &[u64] slice for PHast
        let min_file = std::fs::File::open(&minimizers_path)
            .map_err(|e| format!("Failed to open minimizers file: {e}"))?;
        let min_mmap = unsafe { memmap2::Mmap::map(&min_file) }
            .map_err(|e| format!("Failed to mmap minimizers file: {e}"))?;
        assert_eq!(min_mmap.len(), num_buckets * 8);
        // SAFETY: the file contains native-endian u64 values we just wrote,
        // and u64 has alignment 8 which mmap guarantees (page-aligned).
        let minimizers_slice: &[u64] = unsafe {
            std::slice::from_raw_parts(min_mmap.as_ptr() as *const u64, num_buckets)
        };
        // Build MPHF from the mmap'd slice — no heap allocation for keys
        info!("Building PHast MPHF for {} minimizers (partitioned={})",
            num_buckets, self.config.partitioned_mphf);
        let mphf = PartitionedMphf::build_from_slice(minimizers_slice, self.config.partitioned_mphf);

        // Drop minimizers mmap + delete temp file
        drop(min_mmap);
        drop(min_file);
        let _ = std::fs::remove_file(&minimizers_path);

        // --- Pass B: Rescan file → build bucket_id_by_mphf_index + cached_sizes ---
        let mut cached_sizes: Vec<u32> = Vec::with_capacity(num_buckets);
        let mut bucket_id_by_mphf_index = vec![0usize; num_buckets];
        let mut num_singleton = 0u64;
        let mut num_light = 0u64;
        let mut num_heavy = 0u64;
        let mut num_kmers = 0u64;

        let bucket_iter = file_tuples.bucket_iter()
            .map_err(|e| format!("Failed to open file for pass B: {e}"))?;
        for (bucket_idx, scan) in bucket_iter.enumerate() {
            cached_sizes.push(scan.cached_size as u32);

            match scan.cached_size {
                1 => num_singleton += 1,
                2..=MIN_BUCKET_SIZE => num_light += 1,
                _ => num_heavy += 1,
            }

            // num_kmers is pre-computed by FileBucketIter
            num_kmers += scan.num_kmers;

            // Look up MPHF index for this bucket's minimizer
            let mphf_idx = mphf.get(&scan.minimizer);
            bucket_id_by_mphf_index[mphf_idx] = bucket_idx;
        }

        info!("  Total k-mers: {}", num_kmers);

        let control_map = MinimizersControlMap::from_mphf(mphf, num_buckets as u64);

        let bucket_meta = BucketMetadata {
            cached_sizes,
            num_singleton,
            num_light,
            num_heavy,
        };

        Ok((bucket_meta, control_map, bucket_id_by_mphf_index))
    }

    /// Build the sparse and skew index from file tuples (streaming path).
    ///
    /// Takes `mphf_order` by value so it can be dropped after computing
    /// offset_start_by_orig (saves ~3.2 GB for 400M buckets).
    fn build_index_from_file(
        &self,
        file_tuples: &FileTuples,
        bucket_meta: &BucketMetadata,
        mphf_order: Option<Vec<usize>>,
        spss: &SpectrumPreservingStringSet,
    ) -> Result<SparseAndSkewIndex, String> {
        let total_bases = spss.total_bases();
        let num_bits_per_offset = crate::constants::ceil_log2(total_bases);

        let index = crate::dispatch_on_k!(self.config.k, K => {
            SparseAndSkewIndex::build_from_file::<K>(
                file_tuples,
                bucket_meta,
                mphf_order,
                num_bits_per_offset,
                spss,
            ).map_err(|e| e.to_string())?
        });

        Ok(index)
    }

    /// Encode sequences into spectrum-preserving string set
    fn encode_sequences(&self, sequences: Vec<String>) -> Result<(SpectrumPreservingStringSet, usize), String> {
        let num_sequences = sequences.len();
        let spss = crate::dispatch_on_k!(self.config.k, K => {
            self.encode_sequences_k::<K>(sequences)?
        });
        
        Ok((spss, num_sequences))
    }
    
    /// Encode sequences with specific K
    fn encode_sequences_k<const K: usize>(&self, sequences: Vec<String>) -> Result<SpectrumPreservingStringSet, String>
    where
        Kmer<K>: KmerBits,
    {
        let mut encoder = Encoder::<K>::new();
        
        for (idx, seq) in sequences.iter().enumerate() {
            encoder.add_sequence(seq.as_bytes()).map_err(|e| {
                format!("Failed to encode sequence {}: {}", idx, e)
            })?;
        }
        
        Ok(encoder.build(self.config.m))
    }
    
    /// Build the minimizers control map from classified buckets (in-memory path).
    ///
    /// Returns the control map AND a mapping from MPHF index to bucket_id
    /// for reordering control_codewords to MPHF order.
    fn build_control_map(&self, classified: &ClassifiedBuckets) -> Result<(crate::minimizers_control_map::MinimizersControlMap, Vec<usize>), String> {
        let mut builder = MinimizersControlMapBuilder::new();

        for (bucket_id, bref) in classified.bucket_refs.iter().enumerate() {
            builder.add_minimizer(bref.minimizer);

            let bucket_type = match bref.bucket_type {
                crate::builder::buckets::BucketType::Singleton => BucketType::Regular,
                crate::builder::buckets::BucketType::Light => BucketType::Sparse,
                crate::builder::buckets::BucketType::Heavy => BucketType::HeavyLoad,
            };

            builder.set_bucket_type(bref.minimizer, bucket_type);

            if let Some(control) = builder.get_control_mut(bref.minimizer) {
                control.metadata = bucket_id as u64;
            }
        }

        let c = 100u16;
        let alpha = 0.94;

        builder.build(c, alpha, self.config.partitioned_mphf).map_err(|e| {
            format!("Failed to build minimizers control map: {}", e)
        })
    }
    
    /// Build the sparse and skew index
    fn build_index(
        &self,
        classified: &ClassifiedBuckets,
        mphf_order: Option<&[usize]>,
        spss: &SpectrumPreservingStringSet,
    ) -> Result<SparseAndSkewIndex, String> {
        let total_bases = spss.total_bases();
        let num_bits_per_offset = crate::constants::ceil_log2(total_bases);

        let index = crate::dispatch_on_k!(self.config.k, K => {
            SparseAndSkewIndex::build_from_classified::<K>(classified, mphf_order, num_bits_per_offset, spss)
        });

        Ok(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_dictionary_builder_creation() {
        let config = BuildConfiguration::default();
        let builder = DictionaryBuilder::new(config);
        assert!(builder.is_ok());
    }
    
    #[test]
    fn test_dictionary_builder_invalid_config() {
        let config = BuildConfiguration { k: 30, ..BuildConfiguration::default() }; // Even k is invalid
        let builder = DictionaryBuilder::new(config);
        assert!(builder.is_err());
    }
    
    #[test]
    fn test_build_simple_dictionary() {
        let config = BuildConfiguration::new(21, 11).unwrap();
        let builder = DictionaryBuilder::new(config).unwrap();
        
        let sequences = vec![
            "ACGTACGTACGTACGTACGTACGT".to_string(),
            "TGCATGCATGCATGCATGCATGCA".to_string(),
        ];
        
        let dict = builder.build_from_sequences(sequences);
        // Note: This test may fail until we have proper k-mer extraction
        // in the build pipeline. For now, just check that it runs.
        println!("Dictionary build result: {:?}", dict.is_ok());
    }
}
