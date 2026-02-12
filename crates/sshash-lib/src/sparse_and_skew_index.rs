//! Sparse and Skew Index for k-mer lookup
//!
//! This module implements the SSHash lookup index consisting of:
//! - Control codewords (singleton bucket lookup)
//! - Sparse index (light buckets with linear scan)
//! - Skew index (heavy buckets with MPHF)

use crate::builder::buckets::{Bucket, BucketType, MIN_BUCKET_SIZE};
use crate::constants::{MIN_L, MAX_L, INVALID_UINT64};
use crate::kmer::{Kmer, KmerBits};
use crate::mphf_config::{Mphf, build_mphf_from_slice};
use std::hash::Hash;
use sux::bits::bit_field_vec::BitFieldVec;
use value_traits::slices::{SliceByValue, SliceByValueMut};

/// Sparse and skew index data structure
/// 
/// This index enables O(1) lookup of k-mer positions by:
/// 1. Querying the minimizers control map to get the bucket type
/// 2. Using bucket-specific lookup strategy (singleton/light/heavy)
pub struct SparseAndSkewIndex {
    /// Control codewords for each minimizer (bucket type + pointer)
    /// Size: num_minimizers elements, each (log(N) + 1) bits
    /// LSB encodes bucket type:
    ///   - 0: Singleton bucket (rest is position)
    ///   - 01: Light/mid-load bucket (rest is bucket_id and size)
    ///   - 11: Heavy bucket (rest is partition_id and offset)
    pub control_codewords: BitFieldVec<usize>,
    
    /// Starting positions for light buckets of each size
    /// Index by bucket_size, value is offset into mid_load_buckets
    pub begin_buckets_of_size: Vec<u32>,
    
    /// Positions for all light/mid-load buckets
    /// Linearized storage of all (pos_in_seq) values for buckets of size 2..MIN_BUCKET_SIZE
    pub mid_load_buckets: BitFieldVec<usize>,
    
    /// Skew index for heavy buckets (size > MIN_BUCKET_SIZE)
    pub skew_index: SkewIndex,
}

impl SparseAndSkewIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self {
            control_codewords: BitFieldVec::new(1, 0),
            begin_buckets_of_size: Vec::new(),
            mid_load_buckets: BitFieldVec::new(1, 0),
            skew_index: SkewIndex::new(),
        }
    }
    
    /// Build the sparse and skew index from classified buckets
    ///
    /// # Arguments
    /// * `buckets` - Classified buckets (sorted by minimizer)
    /// * `num_bits_per_offset` - Number of bits needed to represent max offset
    /// * `spss` - The SPSS for decoding k-mers (needed for skew index)
    /// * `k` - K-mer size
    ///
    /// # Returns
    /// The built sparse and skew index
    pub fn build<const K: usize>(
        buckets: Vec<Bucket>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let num_minimizers = buckets.len();
        
        // Initialize control codewords (one per minimizer)
        let mut control_codewords_tmp = vec![0u64; num_minimizers];
        
        // Initialize begin_buckets_of_size (indexed by bucket size 2..=MIN_BUCKET_SIZE)
        let mut begin_buckets_of_size = vec![0u32; MIN_BUCKET_SIZE + 1];
        
        // Separate heavy buckets for skew index (store original index for O(1) codeword update)
        let mut heavy_buckets: Vec<(usize, Bucket)> = Vec::new();
        let mut max_bucket_size = 0;
        
        // Track bucket counts by size using fixed-size array (avoids HashMap overhead)
        let mut bucket_counts_by_size = vec![0usize; MIN_BUCKET_SIZE + 1];
        
        // First pass: count buckets by size and collect heavy buckets
        for (idx, bucket) in buckets.iter().enumerate() {
            let size = bucket.size();
            if size > max_bucket_size {
                max_bucket_size = size;
            }
            
            if size > 1 && size <= MIN_BUCKET_SIZE {
                bucket_counts_by_size[size] += 1;
            } else if size > MIN_BUCKET_SIZE {
                heavy_buckets.push((idx, bucket.clone()));
            }
        }
        
        // Sort heavy buckets by size for partition processing
        heavy_buckets.sort_by_key(|(_idx, b)| b.size());
        
        // Setup begin_buckets_of_size (cumulative offsets)
        let mut cumulative_offset = 0u64;
        for (size, slot) in begin_buckets_of_size.iter_mut().enumerate().take(MIN_BUCKET_SIZE + 1).skip(2) {
            *slot = cumulative_offset as u32;
            let count = bucket_counts_by_size[size];
            cumulative_offset += (count * size) as u64;
        }
        
        // Pre-allocate mid_load_buckets with correct total size.
        // Entries must be grouped by bucket size to match begin_buckets_of_size layout.
        let total_mid_positions: usize = bucket_counts_by_size.iter().enumerate()
            .skip(2)
            .map(|(size, &count)| size * count)
            .sum();
        let mut mid_load_buckets = BitFieldVec::new(num_bits_per_offset, total_mid_positions);
        
        // Track current bucket_id for each size using fixed-size array
        let mut bucket_id_by_size = vec![0usize; MIN_BUCKET_SIZE + 1];
        
        // Second pass: build control codewords and populate mid_load_buckets
        for (minimizer_idx, bucket) in buckets.iter().enumerate() {
            let size = bucket.size();
            
            match bucket.bucket_type {
                BucketType::Singleton => {
                    // Singleton: encode position with LSB = 0 (last bit must be 0)
                    assert_eq!(size, 1);
                    let code = bucket.tuples[0].pos_in_seq << 1; // LSB = 0
                    control_codewords_tmp[minimizer_idx] = code;
                }
                
                BucketType::Light => {
                    // Light/mid-load: encode bucket_id and size with LSB = 01
                    assert!(size > 1 && size <= MIN_BUCKET_SIZE);
                    
                    let bucket_id = bucket_id_by_size[size];
                    bucket_id_by_size[size] += 1;
                    
                    // Encode: | bucket_id (variable) | size-2 (MIN_L bits) | 01 (2 bits) |
                    let size_code = (size - 2) as u64; // sizes 2..=64 -> 0..=62
                    let p = (bucket_id as u64) << MIN_L | size_code;
                    let code = (p << 2) | 0b01; // LSB 2 bits = 01
                    control_codewords_tmp[minimizer_idx] = code;
                    
                    // Place positions at the correct offset in mid_load_buckets.
                    // Entries are grouped by bucket size to match begin_buckets_of_size layout,
                    // mirroring the C++ code which sorts buckets by size before writing.
                    let begin = begin_buckets_of_size[size] as usize;
                    let bucket_start = begin + bucket_id * size;
                    let mut j = 0;
                    let mut prev_pos_in_seq = INVALID_UINT64;
                    for tuple in &bucket.tuples {
                        if tuple.pos_in_seq != prev_pos_in_seq {
                            mid_load_buckets.set_value(bucket_start + j, tuple.pos_in_seq as usize);
                            j += 1;
                            prev_pos_in_seq = tuple.pos_in_seq;
                        }
                    }
                }
                
                BucketType::Heavy => {
                    // Heavy: codewords are set below after skew index build
                }
            }
        }
        
        // Build skew index for heavy buckets
        let skew_index = SkewIndex::build::<K>(
            &heavy_buckets,
            max_bucket_size,
            num_bits_per_offset,
            spss,
            canonical,
        );
        
        // Now update control codewords for heavy buckets with actual partition info
        // We need to track the offset into heavy_load_buckets for each bucket
        let mut heavy_offset_tracker = 0u64;
        let mut partition_id = 0usize;
        let mut lower = MIN_BUCKET_SIZE;
        let mut upper = 2 * lower;
        
        for &(original_idx, ref bucket) in &heavy_buckets {
            let size = bucket.size();
            
            // Advance partition if needed
            while size > upper {
                lower = upper;
                upper = 2 * lower;
                partition_id += 1;
                
                let log2_max = (max_bucket_size as f64).log2().ceil() as usize;
                let num_partitions = if max_bucket_size < (1 << MAX_L) {
                    log2_max - MIN_L
                } else {
                    MAX_L - MIN_L + 1
                };
                
                if partition_id == num_partitions - 1 {
                    upper = max_bucket_size;
                }
            }
            
            // Use stored original index directly — O(1) instead of O(N) linear scan
            let p = (heavy_offset_tracker << 3) | (partition_id as u64);
            let code = (p << 2) | 0b11; // LSB 2 bits = 11
            control_codewords_tmp[original_idx] = code;
            
            // Count positions in this bucket
            heavy_offset_tracker += bucket.size() as u64;
        }
        
        // Use num_bits_per_offset + 1, matching C++ exactly.
        // Singleton codes use num_bits_per_offset + 1 bits (pos_in_seq << 1),
        // and light/heavy codes are bounded by the same width.
        let control_codeword_bits = num_bits_per_offset + 1;
        let mut control_codewords = BitFieldVec::new(control_codeword_bits, num_minimizers);
        for (i, &codeword) in control_codewords_tmp.iter().enumerate() {
            control_codewords.set_value(i, codeword as usize);
        }

        Self {
            control_codewords,
            begin_buckets_of_size,
            mid_load_buckets,
            skew_index,
        }
    }
    
    /// Look up a k-mer using its control codeword
    ///
    /// # Arguments
    /// * `kmer_value` - The k-mer value
    /// * `control_code` - The control codeword from the minimizers control map
    ///
    /// # Returns
    /// The position in the sequence, or INVALID_UINT64 if not found
    pub fn lookup<T: std::hash::Hash>(&self, kmer_value: &T, control_code: u64) -> u64 {
        // Decode the control codeword based on LSB
        // Singleton: last bit = 0 (LSB = 0b0x, where x can be 0 or 1)
        // Light: last 2 bits = 0b01
        // Heavy: last 2 bits = 0b11
        
        if (control_code & 1) == 0 {
            // Singleton: last bit is 0
            // code >> 1 gives the position
            control_code >> 1
        } else if (control_code & 0b11) == 0b01 {
            // Light bucket: extract bucket_id and size
            let p = control_code >> 2;
            let size_code = p & ((1 << MIN_L) - 1);
            let bucket_id = p >> MIN_L;
            let size = size_code + 2; // sizes 0..=62 -> 2..=64
            
            // Get the start offset for buckets of this size
            let begin = self.begin_buckets_of_size[size as usize] as u64;
            let bucket_offset = bucket_id * size;
            
            // Linear search within the bucket
            // TODO: This is a simplified version; the C++ version may use
            // different strategies depending on bucket size
            let start = (begin + bucket_offset) as usize;
            let end = start + (size as usize);
            
            if end <= self.mid_load_buckets.len() {
                // Return first position (placeholder - needs proper k-mer matching)
                self.mid_load_buckets.index_value(start) as u64
            } else {
                INVALID_UINT64
            }
        } else if (control_code & 0b11) == 0b11 {
            // Heavy bucket: use skew index
            self.skew_index.lookup(kmer_value, control_code)
        } else {
            INVALID_UINT64 // Invalid LSB pattern
        }
    }
    
    /// Reorder control_codewords from bucket order to MPHF order.
    ///
    /// After this, `control_codewords[mphf_index]` gives the control codeword
    /// for the minimizer that hashes to `mphf_index`, matching C++ architecture
    /// where `lookup(minimizer) = control_codewords.access(mphf(minimizer))`.
    ///
    /// # Arguments
    /// * `bucket_id_by_mphf_index` - Mapping from MPHF index to bucket_id.
    ///   `bucket_id_by_mphf_index[mphf_index] = bucket_id` where bucket_id
    ///   is the index into the original (pre-reorder) control_codewords.
    pub fn reorder_control_codewords_to_mphf_order(&mut self, bucket_id_by_mphf_index: &[usize]) {
        let n = self.control_codewords.len();
        assert_eq!(n, bucket_id_by_mphf_index.len());
        
        let bit_width = self.control_codewords.bit_width();
        let mut reordered = BitFieldVec::new(bit_width, n);
        
        for (mphf_index, &bucket_id) in bucket_id_by_mphf_index.iter().enumerate() {
            let codeword = self.control_codewords.index_value(bucket_id);
            reordered.set_value(mphf_index, codeword);
        }
        
        self.control_codewords = reordered;
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        let control_bits = (self.control_codewords.len() * self.control_codewords.bit_width()) as u64;
        let begin_bits = (self.begin_buckets_of_size.len() * 32) as u64;
        let mid_bits = (self.mid_load_buckets.len() * self.mid_load_buckets.bit_width()) as u64;
        let skew_bits = self.skew_index.num_bits();
        
        control_bits + begin_bits + mid_bits + skew_bits
    }

    /// Byte size of control codewords (BitFieldVec backing data)
    pub fn control_codewords_bytes(&self) -> usize {
        std::mem::size_of_val(self.control_codewords.as_slice())
    }

    /// Byte size of mid_load_buckets (BitFieldVec backing data)
    pub fn mid_load_buckets_bytes(&self) -> usize {
        std::mem::size_of_val(self.mid_load_buckets.as_slice())
    }

    /// Byte size of `begin_buckets_of_size` (`Vec<u32>`)
    pub fn begin_buckets_of_size_bytes(&self) -> usize {
        self.begin_buckets_of_size.len() * 4
    }

    /// Byte size of skew index data (positions + heavy_load_buckets, excluding MPHFs)
    pub fn skew_index_bytes(&self) -> usize {
        self.skew_index.data_bytes()
    }

    /// Byte size of skew index MPHFs (serialized)
    pub fn skew_mphf_bytes(&self) -> usize {
        self.skew_index.mphf_bytes()
    }
}

impl SparseAndSkewIndex {
    /// Serialize without MPHF (for splitting MPHF across separate container)
    pub fn serialize_without_mphf<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> std::io::Result<()> {
        use epserde::ser::Serialize;
        
        // Serialize control codewords
        // SAFETY: BitFieldVec serialization preserves data layout
        unsafe { self.control_codewords
            .serialize(writer)
            .map_err(std::io::Error::other)? };
        
        // Serialize begin_buckets_of_size manually
        let num_begin_entries = self.begin_buckets_of_size.len() as u32;
        writer.write_all(&num_begin_entries.to_le_bytes())?;
        for &entry in &self.begin_buckets_of_size {
            writer.write_all(&entry.to_le_bytes())?;
        }
        
        // Serialize mid_load_buckets
        // SAFETY: BitFieldVec serialization preserves data layout
        unsafe { self.mid_load_buckets
            .serialize(writer)
            .map_err(std::io::Error::other)? };
        
        // Serialize skew index (excluding MPHF)
        self.skew_index.serialize_without_mphf(writer)?;
        
        Ok(())
    }
    
    /// Deserialize without MPHF (MPHFs loaded separately from container)
    ///
    /// Returns just the SparseAndSkewIndex (MinimizersControlMap is built elsewhere)
    pub fn deserialize_without_mphf<R: std::io::Read>(
        reader: &mut R,
    ) -> std::io::Result<Self> {
        use epserde::deser::Deserialize;
        
        // Deserialize control codewords
        // SAFETY: Reading data serialized by our own serialize method
        let control_codewords = unsafe { BitFieldVec::deserialize_full(reader) }
            .map_err(std::io::Error::other)?;
        
        // Deserialize begin_buckets_of_size manually
        let mut num_begin_bytes = [0u8; 4];
        reader.read_exact(&mut num_begin_bytes)?;
        let num_begin_entries = u32::from_le_bytes(num_begin_bytes) as usize;
        
        let mut begin_buckets_of_size = Vec::with_capacity(num_begin_entries);
        for _ in 0..num_begin_entries {
            let mut entry_bytes = [0u8; 4];
            reader.read_exact(&mut entry_bytes)?;
            begin_buckets_of_size.push(u32::from_le_bytes(entry_bytes));
        }
        
        // Deserialize mid_load_buckets
        // SAFETY: Reading data serialized by our own serialize method
        let mid_load_buckets = unsafe { BitFieldVec::deserialize_full(reader) }
            .map_err(std::io::Error::other)?;
        
        // Deserialize skew index (excluding MPHF)
        let skew_index = SkewIndex::deserialize_without_mphf(reader)?;
        
        let index = Self {
            control_codewords,
            begin_buckets_of_size,
            mid_load_buckets,
            skew_index,
        };
        
        Ok(index)
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
pub struct SkewIndex {
    /// Vector of MPHFs, one per partition
    /// Each MPHF maps k-mers to positions within the partition
    pub mphfs: Vec<Option<Mphf>>,
    
    /// Vector of position arrays, one per partition.
    /// `positions[i]` stores the relative position within each bucket.
    /// Stored as compact vectors with variable bits per position.
    pub positions: Vec<BitFieldVec<usize>>,
    
    /// All positions from heavy buckets (linearized)
    pub heavy_load_buckets: BitFieldVec<usize>,
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
        
        // Serialize heavy_load_buckets
        // SAFETY: BitFieldVec serialization preserves data layout
        unsafe { self.heavy_load_buckets
            .serialize(writer)
            .map_err(std::io::Error::other)? };
        
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
        
        // Deserialize heavy_load_buckets
        // SAFETY: Reading data serialized by our own serialize method
        let heavy_load_buckets = unsafe { BitFieldVec::deserialize_full(reader) }
            .map_err(std::io::Error::other)?;
        
        Ok(Self {
            mphfs,
            positions,
            heavy_load_buckets,
        })
    }

    
    /// Build the skew index from heavy buckets
    ///
    /// # Arguments
    /// * `heavy_buckets` - Buckets with size > MIN_BUCKET_SIZE, sorted by size
    /// * `max_bucket_size` - Maximum bucket size
    /// * `num_bits_per_offset` - Number of bits needed to represent max offset
    /// * `spss` - The SPSS for decoding k-mers
    /// * `k` - K-mer size
    ///
    /// # Returns
    /// The built skew index with partitioned MPHFs
    pub fn build<const K: usize>(
        heavy_buckets: &[(usize, Bucket)],
        max_bucket_size: usize,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
        canonical: bool,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        if heavy_buckets.is_empty() {
            return Self::new();
        }
        
        let min_size = MIN_BUCKET_SIZE;
        
        // Calculate number of partitions
        let log2_max_bucket_size = (max_bucket_size as f64).log2().ceil() as usize;
        let num_partitions = if max_bucket_size < min_size {
            0
        } else if max_bucket_size < (1 << MAX_L) {
            log2_max_bucket_size - MIN_L
        } else {
            MAX_L - MIN_L + 1
        };
        
        if num_partitions == 0 {
            return Self::new();
        }
        
        // Initialize storage
        let mut mphfs = Vec::with_capacity(num_partitions);
        let mut positions = Vec::with_capacity(num_partitions);
        let mut heavy_load_buckets = BitFieldVec::with_capacity(num_bits_per_offset, 0);
        
        // Process heavy buckets by partition
        let mut partition_id = 0;
        let mut lower = min_size;
        let mut upper = 2 * lower;
        let mut num_bits_per_pos = MIN_L + 1;
        
        // Temporary storage for current partition.
        // Uses KmerBits::Storage (u64 for K ≤ 31, u128 for K > 31) so the
        // MPHF is built over the native k-mer representation.
        let mut partition_kmers: Vec<<Kmer<K> as KmerBits>::Storage> = Vec::new();
        let mut partition_positions: Vec<usize> = Vec::new();
        
        for (_orig_idx, bucket) in heavy_buckets {
            let bucket_size = bucket.size();
            
            // Check if we need to finalize current partition and start new one
            while bucket_size > upper {
                // Build MPHF for current partition if it has k-mers
                if !partition_kmers.is_empty() {
                    let mphf_opt = build_partition_mphf(&partition_kmers);
                    
                    // Reorder positions according to MPHF
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
                
                // Move to next partition
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
            
            for tuple in &bucket.tuples {
                if tuple.pos_in_seq != prev_pos_in_seq {
                    prev_pos_in_seq = tuple.pos_in_seq;
                    // Store minimizer position once per unique position
                    heavy_load_buckets.push(tuple.pos_in_seq as usize);
                    pos_in_bucket += 1;
                }
                
                // Extract all k-mers from this super-k-mer using decode_kmer_at.
                // Compute starting position once, then offset per kmer — no binary search.
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
        
        Self {
            mphfs,
            positions,
            heavy_load_buckets,
        }
    }
    
    /// Lookup a k-mer in the skew index
    ///
    /// Generic over the key type: `u64` for K ≤ 31, `u128` for K > 31.
    /// The key type must match the one used during [`build`].
    ///
    /// # Arguments
    /// * `kmer_value` - The k-mer value (same type used to build the MPHF)
    /// * `code` - Control codeword from sparse index
    ///
    /// # Returns
    /// The position in the sequence
    pub fn lookup<T: Hash>(&self, kmer_value: &T, code: u64) -> u64 {
        let code = code >> 2;
        let partition_id = (code & 7) as usize;
        let begin = code >> 3;
        
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
        
        let pos_in_bucket = self.positions[partition_id].index_value(pos) as u64;
        
        self.heavy_load_buckets.index_value((begin + pos_in_bucket) as usize) as u64
    }
    
    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        // MPHF memory is estimated at ~4 bits per key
        let total_keys: usize = self.positions.iter().map(|p| p.len()).sum();
        let mphf_bits = (total_keys as f64 * 4.0) as u64;
        
        let pos_bits: u64 = self.positions.iter()
            .map(|p| (p.len() * p.bit_width()) as u64)
            .sum();
        let heavy_bits = (self.heavy_load_buckets.len() * self.heavy_load_buckets.bit_width()) as u64;
        
        mphf_bits + pos_bits + heavy_bits
    }

    /// Byte size of data (positions + heavy_load_buckets), excluding MPHFs
    pub fn data_bytes(&self) -> usize {
        let pos_bytes: usize = self.positions.iter()
            .map(|p| std::mem::size_of_val(p.as_slice()))
            .sum();
        let heavy_bytes = std::mem::size_of_val(self.heavy_load_buckets.as_slice());
        pos_bytes + heavy_bytes
    }

    /// Byte size of serialized MPHFs
    pub fn mphf_bytes(&self) -> usize {
        self.mphfs.iter()
            .filter_map(|opt| opt.as_ref())
            .map(|mphf| mphf.write_bytes())
            .sum()
    }

    /// Create a new empty skew index
    pub fn new() -> Self {
        Self {
            mphfs: Vec::new(),
            positions: Vec::new(),
            heavy_load_buckets: BitFieldVec::new(1, 0),
        }
    }

    /// Get the number of MPHF partitions
    pub fn num_partitions(&self) -> usize {
        self.mphfs.len()
    }

    /// Get a reference to the MPHFs vector
    pub fn mphfs_ref(&self) -> &[Option<Mphf>] {
        &self.mphfs
    }

    /// Set the MPHFs vector (used during deserialization)
    pub fn set_mphfs(&mut self, mphfs: Vec<Option<Mphf>>) {
        self.mphfs = mphfs;
    }
}

/// Build an MPHF for a partition's k-mers.
///
/// Generic over the key type (`u64` for K ≤ 31, `u128` for K > 31).
fn build_partition_mphf<T: Hash + Clone>(kmers: &[T]) -> Option<Mphf> {
    if kmers.is_empty() {
        return None;
    }
    
    // Build MPHF using PHast with ahash
    Some(build_mphf_from_slice(kmers))
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
        assert_eq!(index.control_codewords.len(), 0);
        assert_eq!(index.mid_load_buckets.len(), 0);
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
        
        let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
        
        // Should have 1 control codeword
        assert_eq!(index.control_codewords.len(), 1);
        
        // Decode: position should be 50, LSB should be 0
        let code = index.control_codewords.index_value(0) as u64;
        assert_eq!(code & 1, 0); // Last bit = 0 for singleton
        assert_eq!(code >> 1, 50); // Position = 50
        
        // Mid-load buckets should be empty
        assert_eq!(index.mid_load_buckets.len(), 0);
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
        
        let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
        
        // Should have 1 control codeword
        assert_eq!(index.control_codewords.len(), 1);
        
        // Decode: LSB 2 bits should be 01 for light bucket
        let code = index.control_codewords.index_value(0) as u64;
        assert_eq!(code & 0b11, 0b01); // LSB 2 bits = 01
        
        // Mid-load buckets should have 3 positions
        assert_eq!(index.mid_load_buckets.len(), 3);
        assert_eq!(index.mid_load_buckets.index_value(0), 100);
        assert_eq!(index.mid_load_buckets.index_value(1), 101);
        assert_eq!(index.mid_load_buckets.index_value(2), 102);
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
        
        let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
        
        assert_eq!(index.control_codewords.len(), 3);
        
        // First bucket: singleton
        assert_eq!(index.control_codewords.index_value(0) & 1, 0);
        
        // Second bucket: light
        assert_eq!(index.control_codewords.index_value(1) & 0b11, 0b01);
        
        // Third bucket: singleton
        assert_eq!(index.control_codewords.index_value(2) & 1, 0);
        
        // Mid-load should have 2 positions (from second bucket)
        assert_eq!(index.mid_load_buckets.len(), 2);
    }
}
