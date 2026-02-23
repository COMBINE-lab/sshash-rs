//! Sparse and Skew Index for k-mer lookup
//!
//! This module implements the SSHash lookup index matching the C++ sshash
//! architecture:
//! - `num_super_kmers_before_bucket`: Elias-Fano monotone sequence for bucket
//!   boundary routing (n = num_buckets + 1)
//! - `offsets`: flat BitFieldVec storing one pos_in_seq per super-kmer
//! - `SkewIndex`: partitioned MPHFs for heavy buckets (size > MIN_BUCKET_SIZE)
//!
//! Lookup path:
//!   bucket_id = mphf(minimizer)
//!   (begin, end) = locate_bucket(bucket_id)
//!   n = end - begin
//!   if ceil_log2(n) > MIN_L: skew_index.lookup(kmer, ceil_log2(n)) -> within-bucket pos -> offsets[begin + pos]
//!   else:                     linear scan offsets[begin..end]

use crate::builder::buckets::{Bucket, ClassifiedBuckets, MIN_BUCKET_SIZE};
use crate::builder::dictionary_builder::BucketMetadata;
use crate::builder::external_sort::FileTuples;
use crate::builder::minimizer_tuples::MinimizerTuple;
use crate::constants::{ceil_log2, MIN_L, MAX_L, INVALID_UINT64};
use crate::kmer::{Kmer, KmerBits};
use crate::partitioned_mphf::PartitionedMphf;
use crate::offsets::EliasFanoOffsets;
use std::hash::Hash;
use sux::bits::bit_field_vec::BitFieldVec;
use tracing::info;
use value_traits::slices::{SliceByValue, SliceByValueMut};

/// Sparse and skew index data structure (C++ compatible architecture)
///
/// This index enables O(1) lookup of k-mer positions by:
/// 1. Querying the minimizers control map to get the bucket_id
/// 2. Using locate_bucket(bucket_id) to find the range [begin, end) in offsets
/// 3. Size-based routing: linear scan for small buckets, MPHF for heavy buckets
pub struct SparseAndSkewIndex {
    /// Monotone Elias-Fano sequence: extras[i] = cumulative (bucket_size - 1)
    /// for buckets 0..i. Used to locate bucket boundaries:
    ///   begin = ef.access(id) + id
    ///   end   = ef.access(id+1) + id + 1
    pub num_super_kmers_before_bucket: EliasFanoOffsets,

    /// One entry per super-kmer: absolute offset (pos_in_seq) into SPSS.
    /// Entries are ordered by bucket (in MPHF order), then by position within
    /// each bucket.
    pub offsets: BitFieldVec<usize>,

    /// Skew index for heavy buckets (size > MIN_BUCKET_SIZE)
    pub skew_index: SkewIndex,
}

impl SparseAndSkewIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self {
            num_super_kmers_before_bucket: EliasFanoOffsets::from_vec(&[0]),
            offsets: BitFieldVec::new(1, 0),
            skew_index: SkewIndex::new(),
        }
    }

    /// Locate a bucket's range in the offsets array.
    ///
    /// Returns `(begin, end)` such that `offsets[begin..end]` contains
    /// all super-kmer positions for the given bucket.
    #[inline]
    pub fn locate_bucket(&self, bucket_id: usize) -> (usize, usize) {
        let begin = self.num_super_kmers_before_bucket.access(bucket_id) as usize + bucket_id;
        let end = self.num_super_kmers_before_bucket.access(bucket_id + 1) as usize + bucket_id + 1;
        (begin, end)
    }

    /// Build the sparse and skew index from owned `Bucket` objects (legacy API).
    ///
    /// Prefer [`build_from_classified`] for production use — it avoids
    /// duplicating the tuple payload.
    pub fn build<const K: usize>(
        buckets: &[Bucket],
        mphf_order: Option<&[usize]>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let num_buckets = buckets.len();
        let t0 = std::time::Instant::now();

        let mut extras = Vec::with_capacity(num_buckets + 1);
        extras.push(0u64);
        let mut total_super_kmers: usize = 0;
        let mut max_bucket_size: usize = 0;
        let mut offset_start_by_orig: Vec<u32> = vec![0; num_buckets];
        for i in 0..num_buckets {
            let orig_idx = mphf_order.map_or(i, |order| order[i]);
            let size = buckets[orig_idx].size();
            if size > max_bucket_size { max_bucket_size = size; }
            offset_start_by_orig[orig_idx] = total_super_kmers as u32;
            total_super_kmers += size;
            let prev = *extras.last().unwrap();
            extras.push(prev + (size as u64).saturating_sub(1));
        }
        let num_super_kmers_before_bucket = EliasFanoOffsets::from_vec(&extras);
        drop(extras);
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size);
        let t1 = std::time::Instant::now();

        let mut offsets = BitFieldVec::new(num_bits_per_offset.max(1), total_super_kmers);
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = Vec::new();
        for orig_idx in 0..num_buckets {
            let bucket = &buckets[orig_idx];
            let size = bucket.size();
            let mut idx = offset_start_by_orig[orig_idx] as usize;
            let mut prev_pos_in_seq = INVALID_UINT64;
            for tuple in &bucket.tuples {
                if tuple.pos_in_seq != prev_pos_in_seq {
                    offsets.set_value(idx, tuple.pos_in_seq as usize);
                    idx += 1;
                    prev_pos_in_seq = tuple.pos_in_seq;
                }
            }
            if size > MIN_BUCKET_SIZE { heavy_buckets_for_skew.push((orig_idx, size)); }
        }
        drop(offset_start_by_orig);
        info!("  Pass 2 (fill offsets): {:.2}s — {} heavy buckets collected",
            t1.elapsed().as_secs_f64(), heavy_buckets_for_skew.len());
        let t2 = std::time::Instant::now();

        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);
        let heavy_refs: Vec<(usize, &Bucket)> = heavy_buckets_for_skew.iter()
            .map(|&(orig_idx, _size)| (orig_idx, &buckets[orig_idx]))
            .collect();
        let skew_index = SkewIndex::build::<K>(&heavy_refs, max_bucket_size, spss, canonical);
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Self { num_super_kmers_before_bucket, offsets, skew_index }
    }

    /// Build the sparse and skew index from in-place classified buckets.
    ///
    /// This is the production path — avoids tuple duplication.
    pub fn build_from_classified<const K: usize>(
        classified: &ClassifiedBuckets,
        mphf_order: Option<&[usize]>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let num_buckets = classified.num_buckets();
        let t0 = std::time::Instant::now();

        // --- Pass 1: Build EF + compute per-bucket offset positions ---
        let mut extras = Vec::with_capacity(num_buckets + 1);
        extras.push(0u64);
        let mut total_super_kmers: usize = 0;
        let mut max_bucket_size: usize = 0;
        let mut offset_start_by_orig: Vec<u32> = vec![0; num_buckets];
        for i in 0..num_buckets {
            let orig_idx = mphf_order.map_or(i, |order| order[i]);
            let size = classified.bucket_refs[orig_idx].cached_size;
            if size > max_bucket_size {
                max_bucket_size = size;
            }
            offset_start_by_orig[orig_idx] = total_super_kmers as u32;
            total_super_kmers += size;
            let prev = *extras.last().unwrap();
            extras.push(prev + (size as u64).saturating_sub(1));
        }
        let num_super_kmers_before_bucket = EliasFanoOffsets::from_vec(&extras);
        drop(extras);
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size);
        let t1 = std::time::Instant::now();

        // --- Pass 2: Fill flat offsets array ---
        let mut offsets = BitFieldVec::new(num_bits_per_offset.max(1), total_super_kmers);
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = Vec::new();
        #[allow(clippy::needless_range_loop)] // orig_idx indexes multiple arrays + method calls
        for orig_idx in 0..num_buckets {
            let bref = &classified.bucket_refs[orig_idx];
            let size = bref.cached_size;
            let mut idx = offset_start_by_orig[orig_idx] as usize;
            let mut prev_pos_in_seq = INVALID_UINT64;
            for tuple in classified.bucket_tuples(orig_idx) {
                if tuple.pos_in_seq != prev_pos_in_seq {
                    offsets.set_value(idx, tuple.pos_in_seq as usize);
                    idx += 1;
                    prev_pos_in_seq = tuple.pos_in_seq;
                }
            }

            if size > MIN_BUCKET_SIZE {
                heavy_buckets_for_skew.push((orig_idx, size));
            }
        }
        drop(offset_start_by_orig);

        info!("  Pass 2 (fill offsets): {:.2}s — {} heavy buckets collected",
            t1.elapsed().as_secs_f64(), heavy_buckets_for_skew.len());
        let t2 = std::time::Instant::now();

        // Sort heavy buckets by size for partition processing
        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);

        // --- Pass 3: Build skew index for heavy buckets ---
        let skew_index = SkewIndex::build_from_classified::<K>(
            classified,
            &heavy_buckets_for_skew,
            max_bucket_size,
            spss,
            canonical,
        );
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Self {
            num_super_kmers_before_bucket,
            offsets,
            skew_index,
        }
    }

    /// Build the sparse and skew index from file tuples (streaming path).
    ///
    /// This avoids materializing all tuples in memory. Instead it:
    /// - Pass 1: Uses `cached_sizes` + `mphf_order` to build EF and compute
    ///   `offset_start_by_orig`. Then drops `mphf_order` (saves ~3.2 GB).
    /// - Pass 2: Reads file sequentially via `SequentialTupleReader` (no mmap).
    ///   Fills offsets and collects heavy bucket tuples during the forward scan.
    /// - Pass 3: Builds skew index from materialized heavy buckets.
    pub fn build_from_file<const K: usize>(
        file_tuples: &FileTuples,
        bucket_meta: &BucketMetadata,
        mphf_order: Option<Vec<usize>>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Result<Self, std::io::Error>
    where
        Kmer<K>: KmerBits,
    {
        let num_buckets = bucket_meta.num_buckets();
        let t0 = std::time::Instant::now();

        // --- Pass 1: Build EF + compute per-bucket offset positions ---
        // Then drop mphf_order to free ~3.2 GB.
        let mut extras = Vec::with_capacity(num_buckets + 1);
        extras.push(0u64);
        let mut total_super_kmers: usize = 0;
        let mut max_bucket_size: usize = 0;
        let mut offset_start_by_orig: Vec<u32> = vec![0; num_buckets];
        for i in 0..num_buckets {
            let orig_idx = mphf_order.as_ref().map_or(i, |order| order[i]);
            let size = bucket_meta.cached_sizes[orig_idx] as usize;
            if size > max_bucket_size {
                max_bucket_size = size;
            }
            offset_start_by_orig[orig_idx] = total_super_kmers as u32;
            total_super_kmers += size;
            let prev = *extras.last().unwrap();
            extras.push(prev + (size as u64).saturating_sub(1));
        }
        let num_super_kmers_before_bucket = EliasFanoOffsets::from_vec(&extras);
        drop(extras);
        // Free mphf_order — no longer needed after offset_start_by_orig is computed
        drop(mphf_order);
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size);
        let t1 = std::time::Instant::now();

        // --- Pass 2: Sequential file scan → fill offsets + collect heavy tuples ---
        // Uses a SequentialTupleReader (BufReader) — no mmap pages kept resident.
        let mut offsets = BitFieldVec::new(num_bits_per_offset.max(1), total_super_kmers);
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = Vec::new();
        let mut heavy_tuples: Vec<Vec<MinimizerTuple>> = Vec::new();
        let mut bucket_idx = 0usize;

        let mut reader = file_tuples.sequential_reader()?;
        let mut next_tuple = reader.read_next()?;

        while let Some(first) = next_tuple {
            let minimizer = first.minimizer;
            let size = bucket_meta.cached_sizes[bucket_idx] as usize;
            let is_heavy = size > MIN_BUCKET_SIZE;

            let mut write_idx = offset_start_by_orig[bucket_idx] as usize;
            let mut prev_pos_in_seq = INVALID_UINT64;

            if is_heavy {
                // Collect heavy bucket tuples during the forward scan
                let mut bucket_tuples = Vec::new();

                // Process the first tuple
                let first_pos = first.pos_in_seq;
                if first_pos != prev_pos_in_seq {
                    offsets.set_value(write_idx, first_pos as usize);
                    write_idx += 1;
                    prev_pos_in_seq = first_pos;
                }
                bucket_tuples.push(first.to_internal());

                // Read remaining tuples in this bucket
                next_tuple = reader.read_next()?;
                while let Some(t) = next_tuple {
                    if t.minimizer != minimizer { break; }
                    let t_pos = t.pos_in_seq;
                    if t_pos != prev_pos_in_seq {
                        offsets.set_value(write_idx, t_pos as usize);
                        write_idx += 1;
                        prev_pos_in_seq = t_pos;
                    }
                    bucket_tuples.push(t.to_internal());
                    next_tuple = reader.read_next()?;
                }

                heavy_buckets_for_skew.push((heavy_tuples.len(), size));
                heavy_tuples.push(bucket_tuples);
            } else {
                // Stream non-heavy bucket tuples without materialization
                let first_pos = first.pos_in_seq;
                if first_pos != prev_pos_in_seq {
                    offsets.set_value(write_idx, first_pos as usize);
                    write_idx += 1;
                    prev_pos_in_seq = first_pos;
                }

                next_tuple = reader.read_next()?;
                while let Some(ext) = next_tuple {
                    if ext.minimizer != minimizer { break; }
                    let ext_pos = ext.pos_in_seq;
                    if ext_pos != prev_pos_in_seq {
                        offsets.set_value(write_idx, ext_pos as usize);
                        write_idx += 1;
                        prev_pos_in_seq = ext_pos;
                    }
                    next_tuple = reader.read_next()?;
                }
            }

            bucket_idx += 1;
        }
        drop(offset_start_by_orig);

        info!("  Pass 2 (fill offsets from file): {:.2}s — {} heavy buckets materialized",
            t1.elapsed().as_secs_f64(), heavy_tuples.len());
        let t2 = std::time::Instant::now();

        // Sort heavy buckets by size for partition processing
        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);

        // --- Pass 3: Build skew index from materialized heavy bucket tuples ---
        let tuples_iter = heavy_buckets_for_skew.iter()
            .map(|&(heavy_idx, _size)| {
                let tuples = &heavy_tuples[heavy_idx];
                // Compute cached_size from tuples
                let mut bucket_size = 0usize;
                let mut prev = INVALID_UINT64;
                for t in tuples.iter() {
                    if t.pos_in_seq != prev {
                        bucket_size += 1;
                        prev = t.pos_in_seq;
                    }
                }
                (bucket_size, tuples.as_slice())
            });
        let skew_index = SkewIndex::build_inner::<K>(tuples_iter, max_bucket_size, spss, canonical);
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Ok(Self {
            num_super_kmers_before_bucket,
            offsets,
            skew_index,
        })
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        let ef_bits = self.num_super_kmers_before_bucket.num_bits();
        let offsets_bits = (self.offsets.len() * self.offsets.bit_width()) as u64;
        let skew_bits = self.skew_index.num_bits();
        ef_bits + offsets_bits + skew_bits
    }

    /// Byte size of EF sequence
    pub fn ef_bytes(&self) -> usize {
        self.num_super_kmers_before_bucket.num_bytes() as usize
    }

    /// Byte size of offsets (BitFieldVec backing data)
    pub fn offsets_bytes(&self) -> usize {
        std::mem::size_of_val(self.offsets.as_slice())
    }

    /// Byte size of skew index data (positions only, excluding MPHFs)
    pub fn skew_index_bytes(&self) -> usize {
        self.skew_index.data_bytes()
    }

    /// Byte size of skew index MPHFs (serialized)
    pub fn skew_mphf_bytes(&self) -> usize {
        self.skew_index.mphf_bytes()
    }

    /// Number of buckets stored
    pub fn num_buckets(&self) -> usize {
        // EF has num_buckets + 1 entries
        self.num_super_kmers_before_bucket.len().saturating_sub(1)
    }

    /// Total number of super-kmers (offsets entries)
    pub fn num_offsets(&self) -> usize {
        self.offsets.len()
    }
}

impl SparseAndSkewIndex {
    /// Serialize without MPHF (for splitting MPHF across separate container)
    pub fn serialize_without_mphf<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> std::io::Result<()> {
        use epserde::ser::Serialize;

        // Serialize EF sequence
        self.num_super_kmers_before_bucket.write_to(writer)?;

        // Serialize offsets
        // SAFETY: BitFieldVec serialization preserves data layout
        unsafe { self.offsets
            .serialize(writer)
            .map_err(std::io::Error::other)? };

        // Serialize skew index (excluding MPHF)
        self.skew_index.serialize_without_mphf(writer)?;

        Ok(())
    }

    /// Deserialize without MPHF (MPHFs loaded separately from container)
    pub fn deserialize_without_mphf<R: std::io::Read>(
        reader: &mut R,
    ) -> std::io::Result<Self> {
        use epserde::deser::Deserialize;

        // Deserialize EF sequence
        let num_super_kmers_before_bucket = EliasFanoOffsets::read_from(reader)?;

        // Deserialize offsets
        // SAFETY: Reading data serialized by our own serialize method
        let offsets = unsafe { BitFieldVec::deserialize_full(reader) }
            .map_err(std::io::Error::other)?;

        // Deserialize skew index (excluding MPHF)
        let skew_index = SkewIndex::deserialize_without_mphf(reader)?;

        Ok(Self {
            num_super_kmers_before_bucket,
            offsets,
            skew_index,
        })
    }
}

impl Default for SparseAndSkewIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Skew index for heavy buckets (size > MIN_BUCKET_SIZE)
///
/// Heavy buckets are partitioned by size into logarithmic partitions:
/// - Partition 0: buckets of size (64, 128]
/// - Partition 1: buckets of size (128, 256]
/// - Partition 2: buckets of size (256, 512]
/// - etc.
///
/// For each partition, we build a separate MPHF over all k-mers in that partition.
/// `positions[partition][mphf_idx]` stores the within-bucket super-kmer index.
/// The caller resolves to an absolute position via `offsets[begin + within_pos]`.
pub struct SkewIndex {
    /// Vector of MPHFs, one per partition
    /// Each MPHF maps k-mers to positions within the partition
    pub mphfs: Vec<Option<PartitionedMphf>>,

    /// Vector of position arrays, one per partition.
    /// `positions[i]` stores the within-bucket super-kmer index for each k-mer.
    pub positions: Vec<BitFieldVec<usize>>,
}

impl SkewIndex {
    /// Serialize without MPHF (for splitting MPHF across separate container)
    pub fn serialize_without_mphf<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> std::io::Result<()> {
        use epserde::ser::Serialize;

        // Serialize number of partitions
        let num_partitions = self.mphfs.len() as u32;
        writer.write_all(&num_partitions.to_le_bytes())?;

        // Serialize positions (one BitFieldVec per partition)
        writer.write_all(&(self.positions.len() as u32).to_le_bytes())?;
        for pos_vec in &self.positions {
            // SAFETY: BitFieldVec serialization preserves data layout
            unsafe { pos_vec
                .serialize(writer)
                .map_err(std::io::Error::other)? };
        }

        Ok(())
    }

    /// Deserialize without MPHF (MPHFs loaded separately)
    pub fn deserialize_without_mphf<R: std::io::Read>(
        reader: &mut R,
    ) -> std::io::Result<Self> {
        use epserde::deser::Deserialize;

        // Deserialize number of partitions
        let mut num_parts_bytes = [0u8; 4];
        reader.read_exact(&mut num_parts_bytes)?;
        let num_partitions = u32::from_le_bytes(num_parts_bytes) as usize;

        // Initialize with empty MPHFs (will be filled later)
        let mphfs = (0..num_partitions).map(|_| None).collect();

        // Deserialize positions
        let mut pos_bytes = [0u8; 4];
        reader.read_exact(&mut pos_bytes)?;
        let num_pos_vecs = u32::from_le_bytes(pos_bytes) as usize;

        let mut positions = Vec::with_capacity(num_pos_vecs);
        for _ in 0..num_pos_vecs {
            // SAFETY: Reading data serialized by our own serialize method
            let pos_vec = unsafe { BitFieldVec::deserialize_full(reader) }
                .map_err(std::io::Error::other)?;
            positions.push(pos_vec);
        }

        Ok(Self {
            mphfs,
            positions,
        })
    }

    /// Build the skew index from heavy buckets (legacy API for tests).
    ///
    /// See [`build_from_classified`] for the production path.
    pub fn build<const K: usize>(
        heavy_buckets: &[(usize, &Bucket)],
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        // Adapt: extract tuple slices from the owned Bucket refs
        let tuples_iter = heavy_buckets.iter()
            .map(|&(_orig_idx, bucket)| (bucket.size(), bucket.tuples.as_slice()));
        Self::build_inner::<K>(tuples_iter, max_bucket_size, spss, canonical)
    }

    /// Build the skew index from [`ClassifiedBuckets`] heavy bucket refs.
    pub fn build_from_classified<const K: usize>(
        classified: &ClassifiedBuckets,
        heavy_buckets_for_skew: &[(usize, usize)], // (orig_idx, size), sorted by size
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let tuples_iter = heavy_buckets_for_skew.iter()
            .map(|&(orig_idx, _size)| {
                let bref = &classified.bucket_refs[orig_idx];
                (bref.cached_size, classified.bucket_tuples(orig_idx))
            });
        Self::build_inner::<K>(tuples_iter, max_bucket_size, spss, canonical)
    }

    /// Core skew index builder.
    ///
    /// `heavy_buckets` yields `(bucket_size, tuple_slice)` pairs sorted by size.
    fn build_inner<'a, const K: usize>(
        heavy_buckets: impl Iterator<Item = (usize, &'a [MinimizerTuple])>,
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let min_size = MIN_BUCKET_SIZE;

        // Calculate number of partitions
        let log2_max_bucket_size = ceil_log2(max_bucket_size as u64);
        let num_partitions = if max_bucket_size <= min_size {
            0
        } else if max_bucket_size < (1 << MAX_L) {
            log2_max_bucket_size - MIN_L
        } else {
            MAX_L - MIN_L + 1
        };

        if num_partitions == 0 {
            return Self::new();
        }

        info!("    SkewIndex: max_bucket_size={}, num_partitions={}",
            max_bucket_size, num_partitions);

        let mut mphfs = Vec::with_capacity(num_partitions);
        let mut positions = Vec::with_capacity(num_partitions);

        let mut partition_id = 0;
        let mut lower = min_size;
        let mut upper = 2 * lower;
        let mut num_bits_per_pos = MIN_L + 1;

        let mut partition_kmers: Vec<<Kmer<K> as KmerBits>::Storage> = Vec::new();
        let mut partition_positions: Vec<usize> = Vec::new();

        for (bucket_size, bucket_tuples) in heavy_buckets {
            // Finalize current partition(s) if this bucket is too large
            while bucket_size > upper {
                if !partition_kmers.is_empty() {
                    let mphf_opt = build_partition_mphf(&partition_kmers);
                    let reordered_positions = if let Some(ref mphf) = mphf_opt {
                        let mut rp = BitFieldVec::new(num_bits_per_pos, partition_positions.len());
                        for (i, &kmer_value) in partition_kmers.iter().enumerate() {
                            let pos = mphf.get(&kmer_value);
                            rp.set_value(pos, partition_positions[i]);
                        }
                        rp
                    } else {
                        BitFieldVec::new(num_bits_per_pos, 0)
                    };
                    mphfs.push(mphf_opt);
                    positions.push(reordered_positions);
                    partition_kmers.clear();
                    partition_positions.clear();
                }

                lower = upper;
                upper = 2 * lower;
                num_bits_per_pos += 1;
                partition_id += 1;

                if partition_id == num_partitions - 1 {
                    upper = max_bucket_size;
                    num_bits_per_pos = log2_max_bucket_size;
                }
            }

            // Add k-mers from this bucket to current partition
            let mut pos_in_bucket = 0u32;
            let mut prev_pos_in_seq = INVALID_UINT64;

            for tuple in bucket_tuples {
                if tuple.pos_in_seq != prev_pos_in_seq {
                    prev_pos_in_seq = tuple.pos_in_seq;
                    pos_in_bucket += 1;
                }

                debug_assert!(tuple.pos_in_seq >= tuple.pos_in_kmer as u64);
                let starting_kmer_pos =
                    (tuple.pos_in_seq - tuple.pos_in_kmer as u64) as usize;

                for kmer_offset in 0..tuple.num_kmers_in_super_kmer {
                    let kmer: Kmer<K> = spss.decode_kmer_at(starting_kmer_pos + kmer_offset as usize);
                    let mut kmer_value = kmer.bits();

                    if canonical {
                        let rc = kmer.reverse_complement();
                        let rc_value = rc.bits();
                        if rc_value < kmer_value {
                            kmer_value = rc_value;
                        }
                    }

                    partition_kmers.push(kmer_value);
                    partition_positions.push((pos_in_bucket - 1) as usize);
                }
            }
        }

        // Finalize last partition
        if !partition_kmers.is_empty() {
            let mphf_opt = build_partition_mphf(&partition_kmers);
            let reordered_positions = if let Some(ref mphf) = mphf_opt {
                let mut rp = BitFieldVec::new(num_bits_per_pos, partition_positions.len());
                for (i, &kmer_value) in partition_kmers.iter().enumerate() {
                    let pos = mphf.get(&kmer_value);
                    rp.set_value(pos, partition_positions[i]);
                }
                rp
            } else {
                BitFieldVec::new(num_bits_per_pos, 0)
            };
            mphfs.push(mphf_opt);
            positions.push(reordered_positions);
        }

        Self { mphfs, positions }
    }

    /// Lookup a k-mer in the skew index.
    ///
    /// Returns the within-bucket super-kmer index. The caller translates this
    /// to an absolute position via `offsets[begin + result]`.
    ///
    /// # Arguments
    /// * `kmer_value` - The k-mer value (same type used to build the MPHF)
    /// * `log2_bucket_size` - ceil_log2 of the bucket size (determines partition)
    ///
    /// # Returns
    /// The within-bucket super-kmer index, or INVALID_UINT64 if not found
    pub fn lookup<T: Hash>(&self, kmer_value: &T, log2_bucket_size: usize) -> u64 {
        if log2_bucket_size <= MIN_L {
            return INVALID_UINT64;
        }

        let partition_id = (log2_bucket_size - (MIN_L + 1))
            .min(self.mphfs.len().saturating_sub(1));

        if partition_id >= self.mphfs.len() {
            return INVALID_UINT64;
        }

        let mphf = match &self.mphfs[partition_id] {
            Some(m) => m,
            None => return INVALID_UINT64,
        };
        let pos = mphf.get(kmer_value);

        if pos >= self.positions[partition_id].len() {
            return INVALID_UINT64;
        }

        self.positions[partition_id].index_value(pos) as u64
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        // MPHF memory is estimated at ~4 bits per key
        let total_keys: usize = self.positions.iter().map(|p| p.len()).sum();
        let mphf_bits = (total_keys as f64 * 4.0) as u64;

        let pos_bits: u64 = self.positions.iter()
            .map(|p| (p.len() * p.bit_width()) as u64)
            .sum();

        mphf_bits + pos_bits
    }

    /// Byte size of data (positions only), excluding MPHFs
    pub fn data_bytes(&self) -> usize {
        self.positions.iter()
            .map(|p| std::mem::size_of_val(p.as_slice()))
            .sum()
    }

    /// Byte size of serialized MPHFs
    pub fn mphf_bytes(&self) -> usize {
        self.mphfs.iter()
            .filter_map(|opt| opt.as_ref())
            .map(|pmphf| pmphf.write_bytes())
            .sum()
    }

    /// Create a new empty skew index
    pub fn new() -> Self {
        Self {
            mphfs: Vec::new(),
            positions: Vec::new(),
        }
    }

    /// Get the number of MPHF partitions
    pub fn num_partitions(&self) -> usize {
        self.mphfs.len()
    }

    /// Get a reference to the MPHFs vector
    pub fn mphfs_ref(&self) -> &[Option<PartitionedMphf>] {
        &self.mphfs
    }

    /// Set the MPHFs vector (used during deserialization)
    pub fn set_mphfs(&mut self, mphfs: Vec<Option<PartitionedMphf>>) {
        self.mphfs = mphfs;
    }
}

/// Build a PartitionedMphf for a skew index partition's k-mers.
///
/// Generic over the key type (`u64` for K <= 31, `u128` for K > 31).
/// Skew index partitions are typically small enough that a single (non-partitioned)
/// MPHF is used — partitioning is only relevant for the main minimizer MPHF.
fn build_partition_mphf<T: Hash + Clone + Sync + Send>(kmers: &[T]) -> Option<PartitionedMphf> {
    if kmers.is_empty() {
        return None;
    }

    // Skew index MPHFs are always monolithic (single partition) since they cover
    // a single size-class bucket, not the full minimizer set.
    Some(PartitionedMphf::build_from_slice(kmers, false))
}

impl Default for SkewIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::minimizer_tuples::MinimizerTuple;

    #[test]
    fn test_sparse_index_creation() {
        let index = SparseAndSkewIndex::new();
        assert_eq!(index.num_buckets(), 0);
        assert_eq!(index.num_offsets(), 0);
    }

    #[test]
    fn test_sparse_index_build_singleton() {
        use crate::builder::encode::Encoder;

        // Create SPSS for k=31
        let mut encoder = Encoder::<31>::new();
        encoder.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        let spss = encoder.build(13);

        // Create a single singleton bucket
        let bucket = Bucket::new(100, vec![MinimizerTuple::new(100, 50, 0)]);
        let buckets = vec![bucket];

        let index = SparseAndSkewIndex::build::<31>(&buckets, None, 32, &spss, false);

        // Should have 1 bucket
        assert_eq!(index.num_buckets(), 1);

        // locate_bucket should give (0, 1) for singleton
        let (begin, end) = index.locate_bucket(0);
        assert_eq!((begin, end), (0, 1));

        // offsets[0] should be 50 (the pos_in_seq)
        assert_eq!(index.offsets.index_value(0), 50);
    }

    #[test]
    fn test_sparse_index_build_light() {
        use crate::builder::encode::Encoder;

        // Create SPSS
        let mut encoder = Encoder::<31>::new();
        encoder.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        let spss = encoder.build(13);

        // Create a light bucket (size = 3)
        let bucket = Bucket::new(200, vec![
            MinimizerTuple::new(200, 100, 0),
            MinimizerTuple::new(200, 101, 0),
            MinimizerTuple::new(200, 102, 0),
        ]);
        let buckets = vec![bucket];

        let index = SparseAndSkewIndex::build::<31>(&buckets, None, 32, &spss, false);

        // Should have 1 bucket with 3 offsets
        assert_eq!(index.num_buckets(), 1);
        assert_eq!(index.num_offsets(), 3);

        // locate_bucket should give (0, 3) for size-3 bucket
        let (begin, end) = index.locate_bucket(0);
        assert_eq!((begin, end), (0, 3));

        // Check positions
        assert_eq!(index.offsets.index_value(0), 100);
        assert_eq!(index.offsets.index_value(1), 101);
        assert_eq!(index.offsets.index_value(2), 102);
    }

    #[test]
    fn test_sparse_index_build_mixed() {
        use crate::builder::encode::Encoder;

        // Create SPSS
        let mut encoder = Encoder::<31>::new();
        encoder.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        let spss = encoder.build(13);

        // Create mixed buckets
        let buckets = vec![
            Bucket::new(100, vec![MinimizerTuple::new(100, 50, 0)]), // Singleton
            Bucket::new(200, vec![
                MinimizerTuple::new(200, 100, 0),
                MinimizerTuple::new(200, 101, 0),
            ]), // Light (size=2)
            Bucket::new(300, vec![MinimizerTuple::new(300, 200, 0)]), // Singleton
        ];

        let index = SparseAndSkewIndex::build::<31>(&buckets, None, 32, &spss, false);

        assert_eq!(index.num_buckets(), 3);
        assert_eq!(index.num_offsets(), 4); // 1 + 2 + 1

        // Bucket 0: singleton at offset 50
        let (b0, e0) = index.locate_bucket(0);
        assert_eq!((b0, e0), (0, 1));
        assert_eq!(index.offsets.index_value(0), 50);

        // Bucket 1: light with 2 offsets
        let (b1, e1) = index.locate_bucket(1);
        assert_eq!((b1, e1), (1, 3));
        assert_eq!(index.offsets.index_value(1), 100);
        assert_eq!(index.offsets.index_value(2), 101);

        // Bucket 2: singleton at offset 200
        let (b2, e2) = index.locate_bucket(2);
        assert_eq!((b2, e2), (3, 4));
        assert_eq!(index.offsets.index_value(3), 200);
    }
}
