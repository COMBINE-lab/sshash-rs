//! Streaming query for efficient k-mer lookups
//!
//! This module implements streaming queries, which optimize lookup performance
//! when querying consecutive k-mers (sliding window over a sequence).
//!
//! Key optimizations:
//! - Incremental k-mer updates (drop first base, add last base)
//! - Reuse minimizer state across adjacent k-mers
//! - Extend within the same string when possible (avoiding MPHF lookups)
//! - Skip searches when minimizer unchanged and previous lookup failed

use crate::kmer::{Kmer, KmerBits};
use crate::minimizer::{MinimizerInfo, MinimizerIterator, MinimizerIteratorRc};

/// Result of a k-mer lookup
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LookupResult {
    /// Absolute k-mer ID (global across all strings)
    pub kmer_id: u64,
    /// Relative k-mer ID within the string (0 <= kmer_id_in_string < string_size)
    pub kmer_id_in_string: u64,
    /// Bit offset into the string data
    pub kmer_offset: u64,
    /// Orientation: +1 for forward, -1 for reverse complement
    pub kmer_orientation: i64,
    
    /// String ID containing this k-mer
    pub string_id: u64,
    /// Start position of the string (in bases)
    pub string_begin: u64,
    /// End position of the string (in bases)
    pub string_end: u64,
    
    /// Whether the minimizer was found in the index
    pub minimizer_found: bool,
}

impl LookupResult {
    /// Create a new lookup result indicating "not found"
    pub fn not_found() -> Self {
        Self {
            kmer_id: u64::MAX,
            kmer_id_in_string: u64::MAX,
            kmer_offset: u64::MAX,
            kmer_orientation: 1, // Forward by default
            string_id: u64::MAX,
            string_begin: u64::MAX,
            string_end: u64::MAX,
            minimizer_found: true,
        }
    }

    /// Check if this result represents a found k-mer
    #[inline]
    pub fn is_found(&self) -> bool {
        self.kmer_id != u64::MAX
    }

    /// Get the string length
    #[inline]
    pub fn string_length(&self) -> u64 {
        if self.is_found() {
            self.string_end - self.string_begin
        } else {
            0
        }
    }
}

impl Default for LookupResult {
    fn default() -> Self {
        Self::not_found()
    }
}

/// Streaming query engine for efficient consecutive k-mer lookups
///
/// This struct maintains state across multiple lookups to optimize
/// queries for sliding windows over sequences.
///
/// # Example
/// ```no_run
/// use sshash_lib::streaming_query::StreamingQuery;
/// // Assuming we have a dictionary...
/// // let mut query = StreamingQuery::new(&dict, true); // canonical=true
/// // 
/// // Process consecutive k-mers efficiently
/// // let result1 = query.lookup("ACGTACGTACGTACGTACGTACGTACGTACG");
/// // let result2 = query.lookup("CGTACGTACGTACGTACGTACGTACGTACGT"); // Sliding by 1
/// ```
pub struct StreamingQuery<const K: usize>
where
    Kmer<K>: KmerBits,
{
    _k: usize,
    _m: usize,
    _canonical: bool,

    // K-mer state (zero-initialized; `start` flag tracks validity)
    start: bool,
    kmer: Kmer<K>,
    kmer_rc: Kmer<K>,

    // Minimizer state
    minimizer_it: MinimizerIterator,
    minimizer_it_rc: MinimizerIteratorRc,
    curr_mini_info: MinimizerInfo,
    prev_mini_info: MinimizerInfo,
    curr_mini_info_rc: MinimizerInfo,
    prev_mini_info_rc: MinimizerInfo,

    // String extension state
    remaining_string_bases: u64,
    // Buffered SPSS reader: bit_pos is the current read position in the strings
    // (in 2-bit units). buf holds data at bit_pos; buf_avail tracks valid bits.
    buf: u64,
    buf_avail: u32,
    buf_bit_pos: u64, // current position in strings data (2-bit units)

    // Result state
    result: LookupResult,

    // Performance counters
    num_searches: u64,
    num_extensions: u64,
    num_invalid: u64,
    num_negative: u64,
}

impl<const K: usize> StreamingQuery<K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new streaming query engine
    ///
    /// # Arguments
    /// * `k` - K-mer size
    /// * `m` - Minimizer size
    /// * `canonical` - Whether to use canonical k-mers (min of forward/RC)
    pub fn new(k: usize, m: usize, canonical: bool) -> Self {
        assert_eq!(k, K, "k parameter must match const generic K");

        let dummy_mini = MinimizerInfo::new(u64::MAX, 0, 0);

        Self {
            _k: k,
            _m: m,
            _canonical: canonical,
            start: true,
            kmer: Kmer::new(<Kmer<K> as KmerBits>::from_u64(0)),
            kmer_rc: Kmer::new(<Kmer<K> as KmerBits>::from_u64(0)),
            minimizer_it: MinimizerIterator::with_seed(k, m, 1),
            minimizer_it_rc: MinimizerIteratorRc::with_seed(k, m, 1),
            curr_mini_info: dummy_mini,
            prev_mini_info: dummy_mini,
            curr_mini_info_rc: dummy_mini,
            prev_mini_info_rc: dummy_mini,
            remaining_string_bases: 0,
            buf: 0,
            buf_avail: 0,
            buf_bit_pos: 0,
            result: LookupResult::not_found(),
            num_searches: 0,
            num_extensions: 0,
            num_invalid: 0,
            num_negative: 0,
        }
    }

    /// Reset the query state (call this when starting a new sequence)
    pub fn reset(&mut self) {
        self.start = true;
        self.remaining_string_bases = 0;
        self.result = LookupResult::not_found();
        self.minimizer_it.set_position(0);
        self.minimizer_it_rc.set_position(0);
    }

    /// Perform a streaming lookup for a k-mer
    ///
    /// This is the main entry point for queries. For optimal performance,
    /// call this with consecutive k-mers (sliding by 1 base at a time).
    ///
    /// # Arguments
    /// * `kmer_str` - DNA string of length K
    ///
    /// # Returns
    /// A LookupResult indicating whether the k-mer was found and its location
    pub fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        // MVP version without Dictionary integration (always seeds)
        self.lookup_internal(kmer_bytes, None)
    }

    /// Perform a streaming lookup with dictionary integration.
    ///
    /// Accepts a `&Dictionary` at call time rather than storing a reference,
    /// so callers can manage the dictionary's lifetime independently (e.g. via `Arc`).
    #[inline(always)]
    pub fn lookup_with_dict(&mut self, kmer_bytes: &[u8], dict: &crate::dictionary::Dictionary) -> LookupResult {
        // 1. Validation + 2. Compute k-mer and reverse complement
        if self.start {
            if !self.is_valid_kmer_bytes(kmer_bytes) {
                self.num_invalid += 1;
                self.reset();
                return self.result;
            }
            self.kmer = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
            self.kmer_rc = self.kmer.reverse_complement();
        } else {
            let last_base = unsafe { *kmer_bytes.get_unchecked(K - 1) };
            if !self.is_valid_base(last_base) {
                self.num_invalid += 1;
                self.reset();
                return self.result;
            }
            let encoded = crate::encoding::encode_base_unchecked(last_base);
            self.kmer = self.kmer.roll_right_base(encoded);
            let complement = crate::encoding::complement_base(encoded);
            self.kmer_rc = self.kmer_rc.append_base(complement);
        }

        self.curr_mini_info = self.minimizer_it.next(self.kmer);
        self.curr_mini_info_rc = self.minimizer_it_rc.next(self.kmer_rc);

        // 3. Compute result (either extend or search)
        if self.remaining_string_bases == 0 {
            self.seed_with_dict(dict);
        } else {
            self.try_extend(dict);
        }

        // 4. Update state
        self.prev_mini_info = self.curr_mini_info;
        self.prev_mini_info_rc = self.curr_mini_info_rc;
        self.start = false;

        self.result
    }

    #[inline(always)]
    fn lookup_internal(&mut self, kmer_bytes: &[u8], dict_opt: Option<&crate::dictionary::Dictionary>) -> LookupResult {
        // 1. Validation
        let is_valid = if self.start {
            self.is_valid_kmer_bytes(kmer_bytes)
        } else {
            let last_base = unsafe { *kmer_bytes.get_unchecked(K - 1) };
            self.is_valid_base(last_base)
        };

        if !is_valid {
            self.num_invalid += 1;
            self.reset();
            return self.result;
        }

        // 2. Compute k-mer and reverse complement, update minimizers
        if self.start {
            self.kmer = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
            self.kmer_rc = self.kmer.reverse_complement();
        } else {
            let last_base = unsafe { *kmer_bytes.get_unchecked(K - 1) };
            let encoded = crate::encoding::encode_base_unchecked(last_base);
            self.kmer = self.kmer.roll_right_base(encoded);
            let complement = crate::encoding::complement_base(encoded);
            self.kmer_rc = self.kmer_rc.append_base(complement);
        }

        self.curr_mini_info = self.minimizer_it.next(self.kmer);
        self.curr_mini_info_rc = self.minimizer_it_rc.next(self.kmer_rc);

        // 3. Compute result (either extend or search)
        if self.remaining_string_bases == 0 {
            self.seed(dict_opt);
        } else {
            if let Some(dict) = dict_opt {
                self.try_extend(dict);
            } else {
                self.seed(dict_opt);
            }
        }

        // 4. Update state
        self.prev_mini_info = self.curr_mini_info;
        self.prev_mini_info_rc = self.curr_mini_info_rc;
        self.start = false;

        self.result
    }

    /// Validate a full k-mer byte slice
    fn is_valid_kmer_bytes(&self, bytes: &[u8]) -> bool {
        if bytes.len() != K {
            return false;
        }
        for &b in bytes {
            if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
                return false;
            }
        }
        true
    }

    /// Validate a single base
    fn is_valid_base(&self, b: u8) -> bool {
        matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
    }

    /// Perform a full search (seed) for the current k-mer
    ///
    /// This is called when we can't extend within the current string.
    fn seed(&mut self, dict_opt: Option<&crate::dictionary::Dictionary>) {
        self.remaining_string_bases = 0;

        // Optimization: if minimizer unchanged and previous was not found, skip
        if !self.start
            && self.curr_mini_info.value == self.prev_mini_info.value
            && self.curr_mini_info_rc.value == self.prev_mini_info_rc.value
            && !self.result.minimizer_found
        {
            assert_eq!(self.result.kmer_id, u64::MAX);
            self.num_negative += 1;
            return;
        }

        if let Some(dict) = dict_opt {
            let kmer = self.kmer;
            let kmer_rc = self.kmer_rc;
            let mini_fwd = self.curr_mini_info;
            let mini_rc = self.curr_mini_info_rc;

            if self._canonical {
                if mini_fwd.value < mini_rc.value {
                    self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_fwd);
                } else if mini_rc.value < mini_fwd.value {
                    self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_rc);
                } else {
                    self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_fwd);
                    if self.result.kmer_id == u64::MAX {
                        self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_rc);
                    }
                }
            } else {
                self.result = dict.lookup_regular_streaming::<K>(&kmer, mini_fwd);
                let minimizer_found = self.result.minimizer_found;
                if self.result.kmer_id == u64::MAX {
                    assert_eq!(self.result.kmer_orientation, 1);
                    self.result = dict.lookup_regular_streaming::<K>(&kmer_rc, mini_rc);
                    self.result.kmer_orientation = -1;
                    let minimizer_rc_found = self.result.minimizer_found;
                    self.result.minimizer_found = minimizer_rc_found || minimizer_found;
                }
            }

            if self.result.kmer_id == u64::MAX {
                self.num_negative += 1;
                return;
            }

            assert!(self.result.minimizer_found);
            self.num_searches += 1;

            let string_size = self.result.string_end - self.result.string_begin;
            if self.result.kmer_orientation > 0 {
                self.remaining_string_bases =
                    (string_size - K as u64) - self.result.kmer_id_in_string;
            } else {
                self.remaining_string_bases = self.result.kmer_id_in_string;
            }
            self.buf_bit_pos = self.result.string_begin + self.result.kmer_id_in_string;
            self.buf_avail = 0;
        } else {
            self.result = LookupResult::not_found();
            self.num_negative += 1;
        }
    }

    /// Perform a full search (seed) for the current k-mer — non-Option hot path
    #[inline(always)]
    fn seed_with_dict(&mut self, dict: &crate::dictionary::Dictionary) {
        self.remaining_string_bases = 0;

        if !self.start
            && self.curr_mini_info.value == self.prev_mini_info.value
            && self.curr_mini_info_rc.value == self.prev_mini_info_rc.value
            && !self.result.minimizer_found
        {
            debug_assert_eq!(self.result.kmer_id, u64::MAX);
            self.num_negative += 1;
            return;
        }

        let kmer = self.kmer;
        let kmer_rc = self.kmer_rc;
        let mini_fwd = self.curr_mini_info;
        let mini_rc = self.curr_mini_info_rc;

        if self._canonical {
            if mini_fwd.value < mini_rc.value {
                self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_fwd);
            } else if mini_rc.value < mini_fwd.value {
                self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_rc);
            } else {
                self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_fwd);
                if self.result.kmer_id == u64::MAX {
                    self.result = dict.lookup_canonical_streaming::<K>(&kmer, &kmer_rc, mini_rc);
                }
            }
        } else {
            self.result = dict.lookup_regular_streaming::<K>(&kmer, mini_fwd);
            let minimizer_found = self.result.minimizer_found;
            if self.result.kmer_id == u64::MAX {
                debug_assert_eq!(self.result.kmer_orientation, 1);
                self.result = dict.lookup_regular_streaming::<K>(&kmer_rc, mini_rc);
                self.result.kmer_orientation = -1;
                let minimizer_rc_found = self.result.minimizer_found;
                self.result.minimizer_found = minimizer_rc_found || minimizer_found;
            }
        }

        if self.result.kmer_id == u64::MAX {
            self.num_negative += 1;
            return;
        }

        debug_assert!(self.result.minimizer_found);
        self.num_searches += 1;

        let string_size = self.result.string_end - self.result.string_begin;
        if self.result.kmer_orientation > 0 {
            self.remaining_string_bases =
                (string_size - K as u64) - self.result.kmer_id_in_string;
        } else {
            self.remaining_string_bases = self.result.kmer_id_in_string;
        }
        self.buf_bit_pos = self.result.string_begin + self.result.kmer_id_in_string;
        self.buf_avail = 0;
    }

    /// Try to extend within the current string using buffered SPSS reader.
    #[inline(always)]
    fn try_extend(&mut self, dict: &crate::dictionary::Dictionary) {
        let kmer_bits = (K * 2) as u32;

        let expected_bits = if self.result.kmer_orientation > 0 {
            self.buf_bit_pos += 1;
            if self.buf_avail >= kmer_bits + 2 {
                self.buf >>= 2;
                self.buf_avail -= 2;
                self.buf & ((1u64 << kmer_bits) - 1)
            } else {
                self.load_kmer_at(dict, self.buf_bit_pos as usize)
            }
        } else {
            self.buf_bit_pos -= 1;
            self.load_kmer_at(dict, self.buf_bit_pos as usize)
        };

        let kmer_u64 = <Kmer<K> as KmerBits>::to_u64(self.kmer.bits());
        let kmer_rc_u64 = <Kmer<K> as KmerBits>::to_u64(self.kmer_rc.bits());
        if expected_bits == kmer_u64 || expected_bits == kmer_rc_u64 {
            self.num_extensions += 1;
            let ori = self.result.kmer_orientation;
            self.result.kmer_id = (self.result.kmer_id as i64 + ori) as u64;
            self.result.kmer_id_in_string =
                (self.result.kmer_id_in_string as i64 + ori) as u64;
            self.remaining_string_bases -= 1;
            return;
        }

        self.seed_with_dict(dict);
    }

    /// Load a kmer from the SPSS strings at the given absolute base position.
    /// Also refills the internal buffer for subsequent shifts.
    #[inline(always)]
    fn load_kmer_at(&mut self, dict: &crate::dictionary::Dictionary, abs_pos: usize) -> u64 {
        let byte_offset = abs_pos / 4;
        let bit_shift = (abs_pos % 4) * 2;
        let strings = dict.spss().strings_data();
        let raw = unsafe {
            std::ptr::read_unaligned(strings.as_ptr().add(byte_offset) as *const u64)
        };
        let mut result = raw >> bit_shift;
        if bit_shift != 0 {
            let extra = unsafe { *strings.as_ptr().add(byte_offset + 8) };
            result |= (extra as u64) << (64 - bit_shift);
        }
        // Buffer always has 64 valid bits after refill (matching C++ get_word64 behavior)
        self.buf = result;
        self.buf_avail = 64;
        let kmer_bits = (K * 2) as u32;
        let mask = (1u64 << kmer_bits) - 1;
        result & mask
    }

    /// Get the number of full searches performed
    pub fn num_searches(&self) -> u64 {
        self.num_searches
    }

    /// Get the number of extensions (no search needed)
    pub fn num_extensions(&self) -> u64 {
        self.num_extensions
    }

    /// Get the number of positive lookups (found)
    pub fn num_positive_lookups(&self) -> u64 {
        self.num_searches + self.num_extensions
    }

    /// Get the number of negative lookups (not found)
    pub fn num_negative_lookups(&self) -> u64 {
        self.num_negative
    }

    /// Get the number of invalid lookups (malformed input)
    pub fn num_invalid_lookups(&self) -> u64 {
        self.num_invalid
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup_result_creation() {
        let result = LookupResult::not_found();
        assert!(!result.is_found());
        assert_eq!(result.kmer_id, u64::MAX);
    }

    #[test]
    fn test_lookup_result_string_length() {
        let mut result = LookupResult::not_found();
        result.string_begin = 100;
        result.string_end = 200;
        result.kmer_id = 42; // Mark as found
        
        assert_eq!(result.string_length(), 100);
    }

    #[test]
    fn test_streaming_query_creation() {
        let query: StreamingQuery<31> = StreamingQuery::new(31, 13, true);
        assert_eq!(query._k, 31);
        assert_eq!(query._m, 13);
        assert!(query._canonical);
        assert_eq!(query.num_searches(), 0);
    }

    #[test]
    fn test_streaming_query_reset() {
        let mut query: StreamingQuery<31> = StreamingQuery::new(31, 13, false);
        query.num_searches = 10;
        query.num_extensions = 5;
        
        query.reset();
        
        assert!(query.start);
        assert_eq!(query.remaining_string_bases, 0);
    }

    #[test]
    fn test_streaming_query_validation() {
        let query: StreamingQuery<31> = StreamingQuery::new(31, 13, true);
        
        assert!(query.is_valid_kmer_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACG")); // 31 bases
        assert!(!query.is_valid_kmer_bytes(b"ACGT")); // Too short
        assert!(!query.is_valid_kmer_bytes(b"ACGTACGTACGTACGTACGTACGTACGTACGN")); // Invalid base
        
        assert!(query.is_valid_base(b'A'));
        assert!(query.is_valid_base(b'a'));
        assert!(!query.is_valid_base(b'N'));
    }

    #[test]
    fn test_streaming_query_lookup_invalid() {
        let mut query: StreamingQuery<15> = StreamingQuery::new(15, 7, true);
        
        // Invalid: too short
        let result = query.lookup(b"ACGT");
        assert!(!result.is_found());
        assert_eq!(query.num_invalid_lookups(), 1);

        // Invalid: has 'N'
        query.reset();
        let result = query.lookup(b"ACGTACGTACGTACN");
        assert!(!result.is_found());
        assert_eq!(query.num_invalid_lookups(), 2);
    }

    #[test]
    fn test_streaming_query_incremental_update() {
        let mut query: StreamingQuery<9> = StreamingQuery::new(9, 5, false);

        // First lookup
        let _result1 = query.lookup(b"ACGTACGTA");
        assert!(!query.start); // No longer in start state

        // Second lookup (sliding by 1)
        let _result2 = query.lookup(b"CGTACGTAC");
        
        // Even though lookups fail (no dictionary), state should update
        assert!(!query.start);
    }
}

/// Streaming query engine integrated with Dictionary
///
/// This provides the full streaming query functionality by connecting
/// to a Dictionary instance for actual k-mer lookups.
pub struct StreamingQueryEngine<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    dict: &'a crate::dictionary::Dictionary,
    query: StreamingQuery<K>,
}

impl<'a, const K: usize> StreamingQueryEngine<'a, K>
where
    Kmer<K>: KmerBits,
{
    /// Create a new streaming query engine for a dictionary
    pub fn new(dict: &'a crate::dictionary::Dictionary) -> Self {
        let canonical = dict.canonical();
        Self {
            dict,
            query: StreamingQuery::new(dict.k(), dict.m(), canonical),
        }
    }
    
    /// Reset the query state
    pub fn reset(&mut self) {
        self.query.reset();
    }
    
    /// Perform a streaming lookup
    #[inline(always)]
    pub fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        self.query.lookup_with_dict(kmer_bytes, self.dict)
    }
    
    /// Get the number of full searches performed
    pub fn num_searches(&self) -> u64 {
        self.query.num_searches()
    }
    
    /// Get the number of extensions (no search needed)
    pub fn num_extensions(&self) -> u64 {
        self.query.num_extensions()
    }
    
    /// Get statistics
    pub fn stats(&self) -> StreamingQueryStats {
        StreamingQueryStats {
            num_searches: self.query.num_searches(),
            num_extensions: self.query.num_extensions(),
            num_invalid: self.query.num_invalid_lookups(),
            num_negative: self.query.num_negative_lookups(),
        }
    }
}

/// Statistics from streaming queries
#[derive(Debug, Clone)]
pub struct StreamingQueryStats {
    /// Number of full MPHF lookups performed
    pub num_searches: u64,
    /// Number of k-mers resolved by extending from a previous result
    pub num_extensions: u64,
    /// Number of lookups that failed validation (hash collision)
    pub num_invalid: u64,
    /// Number of k-mers not found in the dictionary
    pub num_negative: u64,
}
