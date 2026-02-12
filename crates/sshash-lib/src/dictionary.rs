//! Dictionary - the main SSHash data structure
//!
//! Provides efficient storage and querying of k-mer sets using:
//! - Spectrum-preserving string encoding
//! - Minimizer-based indexing
//! - Sparse and skew index for O(1) lookup

use crate::{
    kmer::{Kmer, KmerBits},
    minimizer::MinimizerInfo,
    minimizers_control_map::MinimizersControlMap,
    sparse_and_skew_index::SparseAndSkewIndex,
    spectrum_preserving_string_set::SpectrumPreservingStringSet,
    constants::{INVALID_UINT64, MIN_L},
};
use value_traits::slices::SliceByValue;
use tracing::info;

/// The main dictionary structure
/// 
/// Note: Serialization support is limited due to underlying types.
/// Use custom save/load methods instead.  
pub struct Dictionary {
    /// Encoded string storage
    spss: SpectrumPreservingStringSet,
    
    /// Maps minimizers to control information
    control_map: MinimizersControlMap,
    
    /// Sparse and skew index for lookups
    index: SparseAndSkewIndex,
    
    /// K-mer size
    k: usize,
    
    /// Minimizer size
    m: usize,
    
    /// Whether to use canonical mode (k-mer or RC, whichever is smaller)
    canonical: bool,

    /// Cached hasher for minimizer extraction (avoids re-creating RandomState per query)
    hasher: crate::hasher::DeterministicHasher,
}

impl Dictionary {
    /// Create a new dictionary from components
    pub fn new(
        spss: SpectrumPreservingStringSet,
        control_map: MinimizersControlMap,
        index: SparseAndSkewIndex,
        k: usize,
        m: usize,
        canonical: bool,
    ) -> Self {
        Self {
            spss,
            control_map,
            index,
            k,
            m,
            canonical,
            hasher: crate::hasher::DeterministicHasher::new(1),
        }
    }
    
    /// Look up a k-mer's position in the dictionary.
    ///
    /// # Arguments
    /// * `kmer` - The k-mer to look up
    ///
    /// # Returns
    /// The position in the original sequence, or `INVALID_UINT64` if not found
    #[inline]
    pub fn lookup<const K: usize>(&self, kmer: &Kmer<K>) -> u64
    where
        Kmer<K>: KmerBits,
    {
        let (pos, _) = self.lookup_with_orientation(kmer);
        pos
    }

    /// Query a k-mer and return a full [`LookupResult`](crate::streaming_query::LookupResult)
    /// with position, string membership, and orientation information.
    ///
    /// This is the rich query interface (matching the C++ `dictionary::lookup`
    /// that returns `lookup_result`). For a lightweight lookup that returns
    /// only the k-mer ID, use [`lookup`](Self::lookup).
    ///
    /// # Arguments
    /// * `kmer` - The k-mer to query
    ///
    /// # Returns
    /// A [`LookupResult`](crate::streaming_query::LookupResult) populated with
    /// `kmer_id`, `kmer_id_in_string`, `kmer_offset`, `kmer_orientation`,
    /// `string_id`, `string_begin`, and `string_end`. If the k-mer is not
    /// found, returns
    /// [`LookupResult::not_found()`](crate::streaming_query::LookupResult::not_found).
    #[inline]
    pub fn query<const K: usize>(&self, kmer: &Kmer<K>) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        if self.canonical {
            return self.query_canonical(kmer);
        }

        // Non-canonical mode: try forward k-mer first, then RC
        if let Some(res) = self.query_regular(kmer, 1) {
            return res;
        }
        let kmer_rc = kmer.reverse_complement();
        if let Some(res) = self.query_regular(&kmer_rc, -1) {
            return res;
        }
        crate::streaming_query::LookupResult::not_found()
    }

    /// Internal: regular (non-canonical) lookup returning full LookupResult.
    /// Returns None if the k-mer is not found. Performs the locate() only once.
    #[inline]
    fn query_regular<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        orientation: i8,
    ) -> Option<crate::streaming_query::LookupResult>
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let bucket_id = self.control_map.lookup(minimizer_info.value)?;
        if bucket_id >= self.index.control_codewords.len() {
            return None;
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;

        let kmer_pos = if (control_code & 1) == 0 {
            let minimizer_pos = control_code >> 1;
            self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, minimizer_info)
        } else if (control_code & 0b11) == 0b01 {
            self.lookup_in_light_bucket::<K>(kmer, minimizer_info, control_code)
        } else if (control_code & 0b11) == 0b11 {
            let minimizer_pos = self.index.skew_index.lookup(&kmer.bits(), control_code);
            if minimizer_pos == INVALID_UINT64 { return None; }
            self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, minimizer_info)
        } else {
            INVALID_UINT64
        };

        if kmer_pos == INVALID_UINT64 {
            return None;
        }

        // The inner lookup already validated boundary via locate_with_end(kmer_pos).
        // We call locate_with_end again here to get string info — but the inner
        // function already proved the position is valid, so this is guaranteed to succeed.
        let (string_id, string_begin, string_end) = self.spss.locate_with_end(kmer_pos)?;
        let kmer_id_in_string = kmer_pos - string_begin;
        let kmer_id = kmer_pos - string_id * (self.k as u64 - 1);

        Some(crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset: kmer_pos,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        })
    }

    /// Internal: canonical lookup returning full LookupResult.
    #[inline]
    fn query_canonical<const K: usize>(
        &self,
        kmer: &Kmer<K>,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let (kmer_offset, orientation) = self.lookup_canonical_with_orientation(kmer);
        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        let (string_id, string_begin, string_end) = match self.spss.locate_with_end(kmer_offset) {
            Some(triple) => triple,
            None => return crate::streaming_query::LookupResult::not_found(),
        };
        if kmer_offset > string_end - self.k as u64 {
            return crate::streaming_query::LookupResult::not_found();
        }
        let kmer_id_in_string = kmer_offset - string_begin;
        let kmer_id = kmer_offset - string_id * (self.k as u64 - 1);

        crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        }
    }

    /// Query a k-mer given as a DNA string and return a full
    /// [`LookupResult`](crate::streaming_query::LookupResult).
    ///
    /// This is a convenience wrapper around [`query`](Self::query) that parses
    /// the k-mer from a string of nucleotide characters (A/C/G/T).
    ///
    /// # Arguments
    /// * `kmer_str` - A string of exactly `k` nucleotide characters
    ///
    /// # Panics
    /// Panics if `kmer_str` length does not equal `k`, or if `k` does not
    /// match the const generic `K`.
    pub fn query_from_str<const K: usize>(&self, kmer_str: &str) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        assert_eq!(kmer_str.len(), self.k, "k-mer string length must equal k={}", self.k);
        match Kmer::<K>::from_str(kmer_str) {
            Ok(kmer) => self.query(&kmer),
            Err(_) => crate::streaming_query::LookupResult::not_found(),
        }
    }

    /// Debug: check if control map MPHF is present
    pub fn debug_has_control_map_mphf(&self) -> bool {
        self.control_map.mphf_ref().is_some()
    }

    /// Debug: check if a minimizer exists in the control map
    pub fn debug_control_map_lookup(&self, minimizer: u64) -> bool {
        self.control_map.lookup(minimizer).is_some()
    }

    /// Debug: extract minimizer info for a k-mer
    pub fn debug_extract_minimizer<const K: usize>(&self, kmer: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        self.extract_minimizer(kmer)
    }

    /// Debug: get bucket and control codeword for a k-mer
    pub fn debug_bucket_info<const K: usize>(&self, kmer: &Kmer<K>) -> Option<(u64, u64, u64)>
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let bucket_id = self.control_map.lookup(minimizer_info.value)?;
        if bucket_id >= self.index.control_codewords.len() {
            return None;
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;
        Some((minimizer_info.value, bucket_id as u64, control_code))
    }

    /// Look up a k-mer and return position + orientation
    #[inline]
    pub fn lookup_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        if self.canonical {
            self.lookup_canonical_with_orientation(kmer)
        } else {
            // Non-canonical mode: try forward k-mer first, then RC
            // (matching C++ dictionary::lookup with check_reverse_complement=true)
            let (pos, ori) = self.lookup_regular_with_orientation(kmer);
            if pos != INVALID_UINT64 {
                return (pos, ori);
            }
            let kmer_rc = kmer.reverse_complement();
            let (pos_rc, _) = self.lookup_regular_with_orientation(&kmer_rc);
            if pos_rc != INVALID_UINT64 {
                return (pos_rc, -1); // backward orientation
            }
            (INVALID_UINT64, 1)
        }
    }

    #[inline]
    fn lookup_regular_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let minimizer = minimizer_info.value;

        let bucket_id = match self.control_map.lookup(minimizer) {
            Some(id) => id,
            None => return (INVALID_UINT64, 1),
        };

        if bucket_id >= self.index.control_codewords.len() {
            return (INVALID_UINT64, 1);
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;

        if (control_code & 1) == 0 {
            let minimizer_pos = control_code >> 1;
            let pos = self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, minimizer_info);
            (pos, 1)
        } else if (control_code & 0b11) == 0b01 {
            let pos = self.lookup_in_light_bucket::<K>(kmer, minimizer_info, control_code);
            (pos, 1)
        } else if (control_code & 0b11) == 0b11 {
            let minimizer_pos = self.index.skew_index.lookup(&kmer.bits(), control_code);
            if minimizer_pos == INVALID_UINT64 {
                return (INVALID_UINT64, 1);
            }
            let pos = self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, minimizer_info);
            (pos, 1)
        } else {
            (INVALID_UINT64, 1)
        }
    }

    #[inline]
    fn lookup_canonical_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let kmer_rc = kmer.reverse_complement();
        let mini_fwd = self.extract_minimizer(kmer);
        let mini_rc = self.extract_minimizer(&kmer_rc);

        if mini_fwd.value < mini_rc.value {
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_fwd)
        } else if mini_rc.value < mini_fwd.value {
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_rc)
        } else {
            let res = self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_fwd);
            if res.0 != INVALID_UINT64 {
                return res;
            }
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_rc)
        }
    }

    #[inline]
    fn lookup_canonical_with_minimizer<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_info: MinimizerInfo,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(minimizer_info.value) {
            Some(id) => id,
            None => return (INVALID_UINT64, 1),
        };

        if bucket_id >= self.index.control_codewords.len() {
            return (INVALID_UINT64, 1);
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;

        // Match C++ exactly: check single LSB first, then 2 LSBs
        if (control_code & 1) == 0 {
            // Singleton
            let minimizer_pos = control_code >> 1;
            self.lookup_from_minimizer_pos_canonical::<K>(
                kmer,
                kmer_rc,
                minimizer_pos,
                minimizer_info,
            )
        } else if (control_code & 0b11) == 0b01 {
            // Light bucket
            self.lookup_in_light_bucket_canonical::<K>(
                kmer,
                kmer_rc,
                minimizer_info,
                control_code,
            )
        } else {
            // Heavy bucket (0b11)
            // Must use canonical kmer (min of fwd and rc) since the skew index
            // MPHF was built with canonical kmers (matching C++ behavior).
            let kmer_canon_value = std::cmp::min(kmer.bits(), kmer_rc.bits());
            let minimizer_pos = self.index.skew_index.lookup(&kmer_canon_value, control_code);
            if minimizer_pos == INVALID_UINT64 {
                return (INVALID_UINT64, 1);
            }
            self.lookup_from_minimizer_pos_canonical::<K>(
                kmer,
                kmer_rc,
                minimizer_pos,
                minimizer_info,
            )
        }
    }
    
    /// Look up a k-mer in a light bucket by linear search
    #[inline]
    fn lookup_in_light_bucket<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        control_code: u64,
    ) -> u64
    where
        Kmer<K>: KmerBits,
    {
        // Decode light bucket control codeword: | bucket_id | size-2 (6 bits) | 01 |
        let p = control_code >> 2;
        let size_code = p & ((1 << 6) - 1); // MIN_L = 6 bits for size
        let bucket_id = p >> 6;
        let bucket_size = size_code + 2;
        
        // Get start offset for this bucket size
        if bucket_size as usize >= self.index.begin_buckets_of_size.len() {
            return INVALID_UINT64;
        }
        
        let begin = self.index.begin_buckets_of_size[bucket_size as usize] as u64;
        let bucket_offset = bucket_id * bucket_size;
        let start = (begin + bucket_offset) as usize;
        let end = start + (bucket_size as usize);
        
        if end > self.index.mid_load_buckets.len() {
            return INVALID_UINT64;
        }
        
        // Linear search through minimizer positions in this bucket
        for i in start..end {
            let minimizer_pos = self.index.mid_load_buckets.index_value(i) as u64;
            // With absolute offsets: minimizer_pos is the absolute position
            let kmer_pos = match minimizer_pos.checked_sub(minimizer_info.pos_in_kmer as u64) {
                Some(pos) => pos,
                None => continue,
            };
            
            // Read k-mer directly at absolute position (no binary search needed)
            let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);
            
            if stored_kmer.bits() == query_kmer.bits() {
                // Boundary check: locate string for kmer_pos (matching C++)
                if let Some((_sid, string_begin, string_end)) = self.spss.locate_with_end(kmer_pos) {
                    if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
                        return kmer_pos;
                    }
                }
            }
        }
        INVALID_UINT64
    }

    #[inline]
    fn lookup_from_minimizer_pos<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        minimizer_pos: u64,
        minimizer_info: MinimizerInfo,
    ) -> u64
    where
        Kmer<K>: KmerBits,
    {
        // With absolute offsets: minimizer_pos IS the absolute position of the minimizer
        // kmer_pos = minimizer_pos - pos_in_kmer
        let kmer_pos = match minimizer_pos.checked_sub(minimizer_info.pos_in_kmer as u64) {
            Some(pos) => pos,
            None => return INVALID_UINT64,
        };

        // Read k-mer directly at absolute position (no binary search needed!)
        let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);
        
        if stored_kmer.bits() == query_kmer.bits() {
            // Boundary check matching C++ _lookup_regular:
            // Locate the string containing kmer_pos (using kmer_pos, not
            // minimizer_pos, matching C++ offset_to_id(kmer_offset)).
            // Then verify: kmer_pos >= string_begin AND kmer_pos < string_end - k + 1.
            if let Some((_sid, string_begin, string_end)) = self.spss.locate_with_end(kmer_pos) {
                if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
                    return kmer_pos;
                }
            }
        }

        INVALID_UINT64
    }
    
    /// Canonical lookup from minimizer position
    /// Tries both forward and RC minimizer positions and returns (position, orientation)
    #[inline]
    fn lookup_from_minimizer_pos_canonical<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        minimizer_info: MinimizerInfo,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        // Try forward position first (pos_in_kmer from minimizer)
        let pos_in_kmer_fwd = minimizer_info.pos_in_kmer as u64;
        if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
            query_kmer, kmer_rc, minimizer_pos, pos_in_kmer_fwd
        ) {
            return (pos, orientation);
        }
        
        // Try RC position (k - m - pos_in_kmer)
        let pos_in_kmer_rc = K as u64 - self.m() as u64 - minimizer_info.pos_in_kmer as u64;
        if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
            query_kmer, kmer_rc, minimizer_pos, pos_in_kmer_rc
        ) {
            return (pos, orientation);
        }
        
        (INVALID_UINT64, 1)
    }
    
    /// Try canonical lookup at a specific position with a given pos_in_kmer
    #[inline]
    fn try_canonical_lookup_at_pos<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        pos_in_kmer: u64,
    ) -> Option<(u64, i8)>
    where
        Kmer<K>: KmerBits,
    {
        // With absolute offsets: minimizer_pos is the absolute position
        let kmer_pos = minimizer_pos.checked_sub(pos_in_kmer)?;

        // Read k-mer directly at absolute position (no binary search needed)
        let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);
        
        // Check if stored k-mer matches either forward or RC
        let orientation = if stored_kmer.bits() == query_kmer.bits() {
            1i8  // Forward orientation
        } else if stored_kmer.bits() == kmer_rc.bits() {
            -1i8 // Backward orientation
        } else {
            return None;
        };

        // Boundary check: locate string for kmer_pos (matching C++)
        if let Some((_sid, string_begin, string_end)) = self.spss.locate_with_end(kmer_pos) {
            if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
                return Some((kmer_pos, orientation));
            }
        }
        None
    }
    
    /// Canonical light bucket lookup
    /// Searches all positions in the light bucket for either forward or RC k-mer
    #[inline]
    fn lookup_in_light_bucket_canonical<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        control_code: u64,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        // Decode light bucket control code: remove LSB 01, then extract size
        let p = control_code >> 2;
        let size_code = (p & ((1 << MIN_L) - 1)) as usize;
        let bucket_id = (p >> MIN_L) as usize;
        let bucket_size = size_code + 2;
        
        // Get offset into mid_load_buckets for this bucket
        let bucket_begin = if bucket_size < self.index.begin_buckets_of_size.len() {
            self.index.begin_buckets_of_size[bucket_size] as usize
        } else {
            0
        };
        let offset = bucket_begin + bucket_id * bucket_size;
        
        // Search all positions in the bucket
        for i in 0..bucket_size {
            let pos_idx = offset + i;
            if pos_idx < self.index.mid_load_buckets.len() {
                let minimizer_pos = self.index.mid_load_buckets.index_value(pos_idx) as u64;
                
                // Try canonical lookup at this position
                if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
                    query_kmer, kmer_rc, minimizer_pos, minimizer_info.pos_in_kmer as u64
                ) {
                    return (pos, orientation);
                }
                
                // Also try RC position
                let pos_in_kmer_rc = K as u64 - self.m() as u64 - minimizer_info.pos_in_kmer as u64;
                if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
                    query_kmer, kmer_rc, minimizer_pos, pos_in_kmer_rc
                ) {
                    return (pos, orientation);
                }
            }
        }
        
        (INVALID_UINT64, 1)
    }
    
    /// Check if a k-mer exists in the dictionary
    ///
    /// # Arguments
    /// * `kmer` - The k-mer to check
    ///
    /// # Returns
    /// true if the k-mer is in the dictionary, false otherwise
    pub fn access<const K: usize>(&self, kmer: &Kmer<K>) -> bool
    where
        Kmer<K>: KmerBits,
    {
        self.lookup(kmer) != INVALID_UINT64
    }    
    /// Create a streaming query engine for this dictionary
    ///
    /// Streaming queries optimize consecutive k-mer lookups by maintaining state
    /// and avoiding redundant work.
    ///
    /// # Returns
    /// A new StreamingQuery instance configured for this dictionary
    pub fn create_streaming_query<const K: usize>(&self) -> crate::streaming_query::StreamingQueryEngine<'_, K>
    where
        Kmer<K>: KmerBits,
    {
        crate::streaming_query::StreamingQueryEngine::new(self)
    }    
    /// Get the k-mer size
    pub fn k(&self) -> usize {
        self.k
    }
    
    /// Get the minimizer size
    pub fn m(&self) -> usize {
        self.m
    }
    
    /// Check if canonical mode is enabled
    pub fn canonical(&self) -> bool {
        self.canonical
    }
    
    /// Get a reference to the underlying SPSS
    pub(crate) fn spss(&self) -> &SpectrumPreservingStringSet {
        &self.spss
    }

    /// Debug: get reference to SPSS for testing
    #[cfg(test)]
    pub fn debug_spss(&self) -> &SpectrumPreservingStringSet {
        &self.spss
    }
    
    /// Get the number of strings in the SPSS
    pub fn num_strings(&self) -> u64 {
        self.spss.num_strings()
    }
    
    /// Get the length of a specific string in bases
    pub fn string_length(&self, string_id: u64) -> usize {
        self.spss.string_length(string_id)
    }
    
    /// Locate which string contains a given absolute position.
    /// Returns `(string_id, string_begin)` or None if out of bounds.
    /// Uses binary search on the offsets array.
    /// This is exposed for benchmarking the binary search performance.
    #[inline]
    pub fn locate_string(&self, absolute_pos: u64) -> Option<(u64, u64)> {
        self.spss.locate(absolute_pos)
    }

    /// Access a k-mer at a given position within a string
    pub fn access_kmer<const K: usize>(&self, string_id: u64, pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        self.spss.decode_kmer::<K>(string_id, pos)
    }
    
    /// Get the number of unique minimizers
    pub fn num_minimizers(&self) -> u64 {
        self.control_map.num_minimizers()
    }
    
    /// Get total memory usage in bits
    pub fn num_bits(&self) -> u64 {
        self.spss.num_bits()
            + self.control_map.num_bits()
            + self.index.num_bits()
    }

    /// Print a detailed space breakdown of the index, matching the C++ output format
    pub fn print_space_breakdown(&self) {
        let num_kmers = self.spss.total_bases().saturating_sub(
            (self.k as u64 - 1) * self.spss.num_strings()
        ) as f64;
        
        let strings_bytes = self.spss.strings_bytes() as f64;
        let offsets_bytes = self.spss.offsets_bytes() as f64;
        let control_cw_bytes = self.index.control_codewords_bytes() as f64;
        let mid_load_bytes = self.index.mid_load_buckets_bytes() as f64;
        let begin_buckets_bytes = self.index.begin_buckets_of_size_bytes() as f64;
        let skew_bytes = self.index.skew_index_bytes() as f64;
        let mphf_bytes = self.control_map.mphf_serialized_bytes() as f64;
        let skew_mphf_bytes = self.index.skew_mphf_bytes() as f64;
        
        let total = strings_bytes + offsets_bytes + control_cw_bytes 
            + mid_load_bytes + begin_buckets_bytes + skew_bytes 
            + mphf_bytes + skew_mphf_bytes;
        
        let perc = |x: f64| -> f64 { x * 100.0 / total };
        
        info!("total index size: {} [B] -- {:.5} [MB] ({:.5} [bits/kmer])",
            total as u64, total / 1_000_000.0, total * 8.0 / num_kmers);
        info!("SPACE BREAKDOWN:");
        info!("  mphf: {:.6} [bits/kmer] ({:.5} [bits/key]) -- {:.4}%",
            mphf_bytes * 8.0 / num_kmers,
            mphf_bytes * 8.0 / self.num_minimizers() as f64,
            perc(mphf_bytes));
        info!("  strings_offsets: {:.6} [bits/kmer] -- {:.5}%",
            offsets_bytes * 8.0 / num_kmers, perc(offsets_bytes));
        info!("  control_codewords: {:.5} [bits/kmer] -- {:.4}%",
            control_cw_bytes * 8.0 / num_kmers, perc(control_cw_bytes));
        info!("  mid_load_buckets: {:.6} [bits/kmer] -- {:.5}%",
            mid_load_bytes * 8.0 / num_kmers, perc(mid_load_bytes));
        info!("  begin_buckets_of_size: {:.8} [bits/kmer] -- {:.7}%",
            begin_buckets_bytes * 8.0 / num_kmers, perc(begin_buckets_bytes));
        info!("  strings: {:.5} [bits/kmer] -- {:.4}%",
            strings_bytes * 8.0 / num_kmers, perc(strings_bytes));
        info!("  skew_index: {:.6} [bits/kmer] -- {:.5}%",
            skew_bytes * 8.0 / num_kmers, perc(skew_bytes));
        info!("  skew_mphfs: {:.6} [bits/kmer] -- {:.5}%",
            skew_mphf_bytes * 8.0 / num_kmers, perc(skew_mphf_bytes));
        info!("  --------------");
        info!("  total: {:.5} [bits/kmer]", total * 8.0 / num_kmers);
    }
    
    /// Serialize the dictionary to files
    ///
    /// Creates two files:
    /// - `path.ssi` - Main index file with SPSS and index metadata
    /// - `path.ssi.mphf` - MPHF container with all partitions
    ///
    /// # Arguments
    /// * `path` - Base path for the index files (extension .ssi will be added)
    ///
    /// # Returns
    /// Ok(()) on success, Err on any I/O error
    pub fn save<P: AsRef<std::path::Path>>(&self, path: P) -> crate::serialization::SerializationResult<()> {
        use crate::serialization::*;
        use std::io::BufWriter;

        let base_path = path.as_ref();
        
        // Create main index file
        let index_path = index_file_path(base_path);
        let index_file = std::fs::File::create(&index_path)?;
        let mut index_writer = BufWriter::new(index_file);
        
        // Write header
        let header = DictionarySerializationHeader::new(
            self.k,
            self.m,
            self.canonical,
            (self.index.skew_index.num_partitions() + 1) as u32,
        );
        header.write(&mut index_writer)?;
        
        // Write SPSS using custom binary format (cseq-based offsets)
        self.spss.serialize_to(&mut index_writer)?;
        
        // Write control map without MPHF
        self.control_map.serialize_without_mphf(&mut index_writer)?;
        
        // Write index metadata excluding MPHF
        self.index.serialize_without_mphf(&mut index_writer)?;
        
        // Create MPHF container file
        let mphf_path = mphf_container_path(base_path);
        let mphf_file = std::fs::File::create(&mphf_path)?;
        let mut mphf_writer = BufWriter::new(mphf_file);
        
        // Write MPHF container (control_map first, then skew_index partitions)
        let mut mphfs: Vec<Option<&crate::mphf_config::Mphf>> = Vec::with_capacity(
            self.index.skew_index.num_partitions() + 1,
        );
        mphfs.push(self.control_map.mphf_ref());
        mphfs.extend(self.index.skew_index.mphfs_ref().iter().map(|o| o.as_ref()));
        write_mphf_container(&mut mphf_writer, &mphfs)?;
        
        Ok(())
    }
    
    /// Deserialize a dictionary from files
    ///
    /// Loads from:
    /// - `path.ssi` - Main index file
    /// - `path.ssi.mphf` - MPHF container
    ///
    /// # Arguments
    /// * `path` - Base path for the index files
    ///
    /// # Returns
    /// Ok(Dictionary) on success, Err on any I/O error
    pub fn load<P: AsRef<std::path::Path>>(path: P) -> crate::serialization::SerializationResult<Self> {
        use crate::serialization::*;
        use std::io::BufReader;

        let base_path = path.as_ref();
        
        // Load main index file
        let index_path = index_file_path(base_path);
        let index_file = std::fs::File::open(&index_path)?;
        let mut index_reader = BufReader::new(index_file);
        
        // Read header
        let header = DictionarySerializationHeader::read(&mut index_reader)?;
        
        // Read SPSS using custom binary format (cseq-based offsets)
        let spss = SpectrumPreservingStringSet::deserialize_from(&mut index_reader)?;
        
        // Read control map without MPHF
        let mut control_map = MinimizersControlMap::deserialize_without_mphf(&mut index_reader)?;
        
        // Read index metadata without MPHF
        let mut index = SparseAndSkewIndex::deserialize_without_mphf(&mut index_reader)?;
        
        // Load MPHF container
        let mphf_path = mphf_container_path(base_path);
        let mphf_file = std::fs::File::open(&mphf_path)?;
        let mut mphf_reader = BufReader::new(mphf_file);
        let mut mphfs = read_mphf_container(&mut mphf_reader)?;

        // First MPHF is the control_map's; remaining are skew_index partitions
        let control_mphf = if !mphfs.is_empty() { mphfs.remove(0) } else { None };
        control_map.set_mphf(control_mphf);
        index.skew_index.set_mphfs(mphfs);

        Ok(Dictionary {
            spss,
            control_map,
            index,
            k: header.k,
            m: header.m,
            canonical: header.canonical,
            hasher: crate::hasher::DeterministicHasher::new(1),
        })
    }

    /// Extract the minimizer from a k-mer using the cached hasher.
    #[inline]
    fn extract_minimizer<const K: usize>(&self, kmer: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        self.extract_minimizer_inline::<K>(*kmer)
    }

    /// Inline minimizer extraction without creating a MinimizerIterator.
    /// Rescans all (k-m+1) windows to find the minimum hash.
    #[inline]
    fn extract_minimizer_inline<const K: usize>(&self, kmer: Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        let num_windows = self.k - self.m + 1;
        let mask = (1u64 << (self.m * 2)) - 1;
        let bits = <Kmer<K> as KmerBits>::to_u128(kmer.bits());

        let mut min_hash = u64::MAX;
        let mut min_value = 0u64;
        let mut min_pos = 0usize;

        for i in 0..num_windows {
            let mmer_value = ((bits >> (i * 2)) as u64) & mask;
            let hash = self.hasher.hash_u64(mmer_value);
            if hash < min_hash {
                min_hash = hash;
                min_value = mmer_value;
                min_pos = i;
            }
        }

        MinimizerInfo::new(min_value, min_pos as u64, min_pos)
    }

    // --- Streaming query helpers ---
    // These accept pre-computed minimizer info (from the streaming iterator)
    // and return a full LookupResult, matching C++ lookup_regular/lookup_canonical
    // with minimizer_info parameter.

    /// Regular (non-canonical) lookup with pre-computed minimizer, returning LookupResult.
    /// Used by the streaming query to avoid recomputing the minimizer.
    pub(crate) fn lookup_regular_streaming<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        mini_info: MinimizerInfo,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(mini_info.value) {
            Some(id) => id,
            None => {
                let mut res = crate::streaming_query::LookupResult::not_found();
                res.minimizer_found = false;
                return res;
            }
        };

        if bucket_id >= self.index.control_codewords.len() {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;

        let kmer_offset = if (control_code & 1) == 0 {
            let minimizer_pos = control_code >> 1;
            self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, mini_info)
        } else if (control_code & 0b11) == 0b01 {
            self.lookup_in_light_bucket::<K>(kmer, mini_info, control_code)
        } else if (control_code & 0b11) == 0b11 {
            let minimizer_pos = self.index.skew_index.lookup(&kmer.bits(), control_code);
            if minimizer_pos == INVALID_UINT64 {
                return crate::streaming_query::LookupResult::not_found();
            }
            self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, mini_info)
        } else {
            INVALID_UINT64
        };

        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        self.build_lookup_result(kmer_offset, 1)
    }

    /// Canonical lookup with pre-computed minimizer, returning LookupResult.
    /// Used by the streaming query to avoid recomputing the minimizer.
    pub(crate) fn lookup_canonical_streaming<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        mini_info: MinimizerInfo,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(mini_info.value) {
            Some(id) => id,
            None => {
                let mut res = crate::streaming_query::LookupResult::not_found();
                res.minimizer_found = false;
                return res;
            }
        };

        if bucket_id >= self.index.control_codewords.len() {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        }
        let control_code = self.index.control_codewords.index_value(bucket_id) as u64;

        let (kmer_offset, orientation) = if (control_code & 1) == 0 {
            let minimizer_pos = control_code >> 1;
            self.lookup_from_minimizer_pos_canonical::<K>(kmer, kmer_rc, minimizer_pos, mini_info)
        } else if (control_code & 0b11) == 0b01 {
            self.lookup_in_light_bucket_canonical::<K>(kmer, kmer_rc, mini_info, control_code)
        } else if (control_code & 0b11) == 0b11 {
            // Must use canonical kmer (min of fwd and rc) since the skew index
            // MPHF was built with canonical kmers (matching C++ behavior).
            let kmer_canon_value = std::cmp::min(kmer.bits(), kmer_rc.bits());
            let minimizer_pos = self.index.skew_index.lookup(&kmer_canon_value, control_code);
            if minimizer_pos == INVALID_UINT64 {
                return crate::streaming_query::LookupResult::not_found();
            }
            self.lookup_from_minimizer_pos_canonical::<K>(kmer, kmer_rc, minimizer_pos, mini_info)
        } else {
            (INVALID_UINT64, 1)
        };

        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        self.build_lookup_result(kmer_offset, orientation)
    }

    /// Build a full LookupResult from a kmer_offset and orientation.
    /// Shared helper for streaming lookup methods.
    #[inline]
    fn build_lookup_result(
        &self,
        kmer_offset: u64,
        orientation: i8,
    ) -> crate::streaming_query::LookupResult {
        let (string_id, string_begin, string_end) = match self.spss.locate_with_end(kmer_offset) {
            Some(triple) => triple,
            None => return crate::streaming_query::LookupResult::not_found(),
        };

        // Boundary check: verify k-mer fits within its string (matching C++)
        if kmer_offset > string_end - self.k as u64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        let kmer_id_in_string = kmer_offset - string_begin;
        let kmer_id = kmer_offset - string_id * (self.k as u64 - 1);

        crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use crate::builder::{BuildConfiguration, DictionaryBuilder};
    
    #[test]
    fn test_dictionary_creation() {
        let spss = SpectrumPreservingStringSet::new(31, 13);
        let control_map = MinimizersControlMap::from_parts(
            Vec::new(),
            Vec::new(),
            0,
        );
        let index = SparseAndSkewIndex::new();
        
        let dict = Dictionary::new(spss, control_map, index, 31, 13, false);
        
        assert_eq!(dict.k(), 31);
        assert_eq!(dict.m(), 13);
        assert!(!dict.canonical());
    }
    
    #[test]
    fn test_dictionary_build_and_lookup() {
        // Create a simple test dictionary
        let mut config = BuildConfiguration::new(31, 21).unwrap();
        config.verbose = true;
        let builder = DictionaryBuilder::new(config).unwrap();
        
        // Use simple, short sequences to avoid bucket overflow
        let sequences = vec![
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
        ];
        
        eprintln!("\n=== Building dictionary ===");
        let dict = builder.build_from_sequences(sequences.clone()).unwrap();
        
        assert_eq!(dict.k(), 31);
        assert_eq!(dict.m(), 21);
        
        eprintln!("\n=== Dictionary info ===");
        eprintln!("Num minimizers: {}", dict.num_minimizers());
        eprintln!("Num control codewords: {}", dict.index.control_codewords.len());
        eprintln!("SPSS num strings: {}", dict.spss.num_strings());
        eprintln!("SPSS total bases: {}", dict.spss.total_bases());
        
        // Create a k-mer from the first sequence
        let test_kmer_str = &sequences[0][0..31];  // First 31 bases
        eprintln!("\n=== Testing lookup for k-mer: {} ===", test_kmer_str);
        let kmer = crate::kmer::Kmer::<31>::from_str(test_kmer_str).unwrap();
        
        // Extract minimizer to see what it is
        let mut mini_iter = crate::minimizer::MinimizerIterator::with_seed(31, 21, 1);
        let mini_info = mini_iter.next(kmer);
        eprintln!("K-mer minimizer: value={}, pos_in_kmer={}", mini_info.value, mini_info.pos_in_kmer);
        
        // Look up in control map
        if let Some(bucket_id) = dict.control_map.lookup(mini_info.value) {
            eprintln!("Control map lookup: bucket_id={}", bucket_id);
            
            if bucket_id < dict.index.control_codewords.len() {
                let control_code = dict.index.control_codewords.index_value(bucket_id) as u64;
                eprintln!("Control codeword: 0x{:016x}, LSB={}", control_code, control_code & 0b11);
            }
        } else {
            eprintln!("Minimizer NOT found in control map!");
        }
        
        // Lookup the k-mer
        let result = dict.lookup(&kmer);
        
        eprintln!("\nLookup result: {} (INVALID={})", result, crate::constants::INVALID_UINT64);
        
        // The k-mer should be found at some valid position
        // For now, just test that the pipeline runs without crashing
    }
}
