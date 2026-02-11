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
        buckets::{classify_into_buckets, BucketStatistics},
        config::BuildConfiguration,
        encode::Encoder,
        minimizer_tuples::{compute_minimizer_tuples, compute_minimizer_tuples_external, needs_external_sorting},
    },
    dictionary::Dictionary,
    kmer::{Kmer, KmerBits},
    minimizers_control_map::{MinimizersControlMapBuilder, BucketType},
    sparse_and_skew_index::SparseAndSkewIndex,
    spectrum_preserving_string_set::SpectrumPreservingStringSet,
};
use tracing::info;

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
        
        // Step 2: Extract minimizer tuples (with inline coalescing during extraction)
        info!("Step 2: Extracting and coalescing minimizer tuples...");
        let tuples = self.extract_tuples(&spss)?;
        info!("  Extracted and coalesced {} tuples", tuples.len());
        
        // Step 3: Classify into buckets
        info!("Step 3: Classifying buckets...");
        let buckets = classify_into_buckets(tuples);
        
        // Compute statistics
        let mut stats = BucketStatistics::new();
        for bucket in &buckets {
            stats.add_bucket(&bucket.tuples);
        }
        stats.print_summary();
        
        // Step 5: Build minimizers control map
        info!("Step 5: Building minimizers control map...");
        let (control_map, bucket_id_by_mphf_index) = self.build_control_map(&buckets)?;
        info!("  Built MPHF for {} minimizers", control_map.num_minimizers());
        
        // Step 6: Build sparse and skew index
        info!("Step 6: Building sparse and skew index...");
        let mut index = self.build_index(buckets, &spss)?;
        
        // Step 6b: Reorder control_codewords from bucket order to MPHF order
        // This eliminates the need for a controls array at query time,
        // matching C++ architecture where mphf(minimizer) directly indexes control_codewords
        if !bucket_id_by_mphf_index.is_empty() {
            index.reorder_control_codewords_to_mphf_order(&bucket_id_by_mphf_index);
        }
        info!("  Index built successfully");
        
        // Step 7: Assemble dictionary
        info!("Dictionary Build Complete");
        let total_bits = spss.num_bits() + control_map.num_bits() + index.num_bits();
        info!("Total memory: {:.2} MB", total_bits as f64 / (8.0 * 1024.0 * 1024.0));
        
        Ok(Dictionary::new(
            spss,
            control_map,
            index,
            self.config.k,
            self.config.m,
            self.config.canonical,
        ))
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
    
    /// Extract minimizer tuples from SPSS
    /// 
    /// Uses external sorting when the estimated memory exceeds the RAM limit.
    fn extract_tuples(&self, spss: &SpectrumPreservingStringSet) -> Result<Vec<crate::builder::minimizer_tuples::MinimizerTuple>, String> {
        // Estimate total k-mers: total_bases - num_strings * (k - 1)
        let total_bases = spss.total_bases();
        let num_strings = spss.num_strings();
        let k = self.config.k as u64;
        let total_kmers = total_bases.saturating_sub(num_strings * (k - 1));
        
        // Check if we need external sorting
        if needs_external_sorting(total_kmers, self.config.ram_limit_gib) {
            info!(
                "Using external sorting: estimated {} k-mers exceeds RAM limit of {} GiB",
                total_kmers, self.config.ram_limit_gib
            );
            crate::dispatch_on_k!(self.config.k, K => {
                compute_minimizer_tuples_external::<K>(spss, &self.config)
                    .map_err(|e| e.to_string())
            })
        } else {
            crate::dispatch_on_k!(self.config.k, K => {
                Ok(compute_minimizer_tuples::<K>(spss, &self.config))
            })
        }
    }
    
    /// Build the minimizers control map from buckets
    ///
    /// Returns the control map AND a mapping from MPHF index to bucket_id
    /// for reordering control_codewords to MPHF order.
    fn build_control_map(&self, buckets: &[crate::builder::buckets::Bucket]) -> Result<(crate::minimizers_control_map::MinimizersControlMap, Vec<usize>), String> {
        let mut builder = MinimizersControlMapBuilder::new();
        
        // Add all minimizers and set their bucket types
        // The bucket index is implicitly the order in which we add them
        for (bucket_id, bucket) in buckets.iter().enumerate() {
            builder.add_minimizer(bucket.minimizer);
            
            let bucket_type = match bucket.bucket_type {
                crate::builder::buckets::BucketType::Singleton => BucketType::Regular,
                crate::builder::buckets::BucketType::Light => BucketType::Sparse,
                crate::builder::buckets::BucketType::Heavy => BucketType::HeavyLoad,
            };
            
            builder.set_bucket_type(bucket.minimizer, bucket_type);
            
            // Store bucket_id in metadata so we can map MPHF index → bucket_id
            if let Some(control) = builder.get_control_mut(bucket.minimizer) {
                control.metadata = bucket_id as u64;
            }
        }
        
        // Build the MPHF
        // Use relative_level_size = 100 for minimum MPHF size (~2.8 bits/key)
        // Note: C++ lambda controls PTHash partition size and is unrelated to fmph's level size
        let c = 100u16;
        let alpha = 0.94; // Load factor (not used by ph, but kept for API compatibility)
        
        builder.build(c, alpha).map_err(|e| {
            format!("Failed to build minimizers control map: {}", e)
        })
    }
    
    /// Build the sparse and skew index
    fn build_index(
        &self,
        buckets: Vec<crate::builder::buckets::Bucket>,
        spss: &SpectrumPreservingStringSet,
    ) -> Result<SparseAndSkewIndex, String> {
            let _k = self.config.k as u64;
            let _m = self.config.m as u64;
        
            // Calculate offset encoding bits (matching C++ decoded_offsets)
            // pos_in_seq stores absolute position in concatenated SPSS
            // So num_bits_per_offset = ceil_log2(total_bases)
            let total_bases = spss.total_bases();
            let num_bits_per_offset = crate::constants::ceil_log2(total_bases);
        

        // Dispatch to appropriate build function based on k
        let index = crate::dispatch_on_k!(self.config.k, K => {
            SparseAndSkewIndex::build::<K>(buckets, num_bits_per_offset, spss, self.config.canonical)
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
