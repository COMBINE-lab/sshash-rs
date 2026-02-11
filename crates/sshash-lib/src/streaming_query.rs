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
use crate::minimizer::{MinimizerInfo, MinimizerIterator};
use crate::encoding::encode_base;

/// Result of a k-mer lookup
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LookupResult {
    /// Absolute k-mer ID (global across all strings)
    pub kmer_id: u64,
    /// Relative k-mer ID within the string (0 <= kmer_id_in_string < string_size)
    pub kmer_id_in_string: u64,
    /// Bit offset into the string data
    pub kmer_offset: u64,
    /// Orientation: +1 for forward, -1 for reverse complement
    pub kmer_orientation: i8,
    
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
    k: usize,
    _m: usize, // Will be used in full Dictionary lookup
    _canonical: bool, // Will be used in full Dictionary lookup
    
    // K-mer state
    start: bool,
    kmer: Option<Kmer<K>>,
    kmer_rc: Option<Kmer<K>>,
    
    // Minimizer state
    minimizer_it: MinimizerIterator,
    minimizer_it_rc: MinimizerIterator,
    curr_mini_info: MinimizerInfo,
    prev_mini_info: MinimizerInfo,
    curr_mini_info_rc: MinimizerInfo,
    prev_mini_info_rc: MinimizerInfo,
    
    // String extension state
    remaining_string_bases: u64,
    
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
            k,
            _m: m,
            _canonical: canonical,
            start: true,
            kmer: None,
            kmer_rc: None,
            minimizer_it: MinimizerIterator::with_seed(k, m, 1),
            minimizer_it_rc: MinimizerIterator::with_seed(k, m, 1),
            curr_mini_info: dummy_mini,
            prev_mini_info: dummy_mini,
            curr_mini_info_rc: dummy_mini,
            prev_mini_info_rc: dummy_mini,
            remaining_string_bases: 0,
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
    pub fn lookup(&mut self, kmer_str: &str) -> LookupResult {
        // MVP version without Dictionary integration (always seeds)
        self.lookup_internal(kmer_str, None)
    }
    
    /// Perform a streaming lookup with dictionary integration
    ///
    /// Internal method used by StreamingQueryEngine
    pub(crate) fn lookup_with_dict(&mut self, kmer_str: &str, dict: &crate::dictionary::Dictionary) -> LookupResult {
        self.lookup_internal(kmer_str, Some(dict))
    }
    
    fn lookup_internal(&mut self, kmer_str: &str, dict_opt: Option<&crate::dictionary::Dictionary>) -> LookupResult {
        // 1. Validation
        let is_valid = if self.start {
            self.is_valid_kmer(kmer_str)
        } else {
            self.is_valid_base(kmer_str.as_bytes()[self.k - 1])
        };

        if !is_valid {
            self.num_invalid += 1;
            self.reset();
            return self.result.clone();
        }

        // 2. Compute k-mer and reverse complement, update minimizers
        if self.start {
            // First k-mer: parse from scratch
            match Kmer::from_str(kmer_str) {
                Ok(km) => {
                    self.kmer = Some(km);
                    let rc = km.reverse_complement();
                    self.kmer_rc = Some(rc);
                    
                    self.curr_mini_info = self.minimizer_it.next(km);
                    self.curr_mini_info_rc = self.minimizer_it_rc.next(rc);
                }
                Err(_) => {
                    self.num_invalid += 1;
                    self.reset();
                    return self.result.clone();
                }
            }
        } else {
            // Update incrementally: drop first base, add new last base
            if let Some(mut km) = self.kmer {
                // Drop first base (shift left)
                for i in 0..(self.k - 1) {
                    let base = km.get_base(i + 1);
                    km.set_base(i, base);
                }
                
                // Add new last base
                let new_base = kmer_str.as_bytes()[self.k - 1];
                if let Ok(encoded) = encode_base(new_base) {
                    km.set_base(self.k - 1, encoded);
                    self.kmer = Some(km);
                    
                    // Update RC: pad (shift right), set first base to complement
                    if let Some(mut km_rc) = self.kmer_rc {
                        for i in (1..self.k).rev() {
                            let base = km_rc.get_base(i - 1);
                            km_rc.set_base(i, base);
                        }
                        
                        // Complement of new base at position 0
                        let complement = crate::encoding::complement_base(encoded);
                        km_rc.set_base(0, complement);
                        self.kmer_rc = Some(km_rc);
                        
                        self.curr_mini_info = self.minimizer_it.next(km);
                        self.curr_mini_info_rc = self.minimizer_it_rc.next(km_rc);
                    }
                }
            }
        }

        // 3. Compute result (either extend or search)
        if self.remaining_string_bases == 0 {
            self.seed(dict_opt);
        } else {
            // Try to extend within current string
            if let Some(dict) = dict_opt {
                self.try_extend(dict);
            } else {
                // No dictionary, can't extend
                self.seed(dict_opt);
            }
        }

        // 4. Update state
        self.prev_mini_info = self.curr_mini_info;
        self.prev_mini_info_rc = self.curr_mini_info_rc;
        self.start = false;

        self.result.clone()
    }

    /// Validate a full k-mer string
    fn is_valid_kmer(&self, s: &str) -> bool {
        if s.len() != self.k {
            return false;
        }
        for &b in s.as_bytes() {
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

        if let (Some(dict), Some(kmer)) = (dict_opt, self.kmer) {
            if self._canonical {
                // Canonical mode: matching C++ lookup_canonical logic in seed
                let kmer_rc = kmer.reverse_complement();
                let mini_fwd = self.curr_mini_info;
                let mini_rc = self.curr_mini_info_rc;

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
                // Regular mode: try forward, then RC with backward orientation
                // Matches C++ streaming_query::seed() for non-canonical
                self.result = dict.lookup_regular_streaming::<K>(&kmer, self.curr_mini_info);
                let minimizer_found = self.result.minimizer_found;
                if self.result.kmer_id == u64::MAX {
                    assert_eq!(self.result.kmer_orientation, 1); // forward
                    let kmer_rc = kmer.reverse_complement();
                    self.result = dict.lookup_regular_streaming::<K>(&kmer_rc, self.curr_mini_info_rc);
                    self.result.kmer_orientation = -1; // backward
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

            // Calculate remaining bases for extension, matching C++ exactly:
            //   forward:  (string_end - string_begin - k) - kmer_id_in_string
            //   backward: kmer_id_in_string
            let string_size = self.result.string_end - self.result.string_begin;
            if self.result.kmer_orientation > 0 {
                self.remaining_string_bases =
                    (string_size - self.k as u64) - self.result.kmer_id_in_string;
            } else {
                self.remaining_string_bases = self.result.kmer_id_in_string;
            }
        } else {
            // No dictionary available
            self.result = LookupResult::not_found();
            self.num_negative += 1;
        }
    }
    
    /// Try to extend within the current string
    ///
    /// Matches C++ streaming_query extension logic:
    /// - Read the expected next k-mer from the string data
    /// - If it matches the current k-mer (or its RC), update result fields
    fn try_extend(&mut self, dict: &crate::dictionary::Dictionary) {
        if let (Some(kmer), Some(kmer_rc)) = (self.kmer, self.kmer_rc) {
            // Compute the absolute position of the expected next k-mer
            // C++: kmer_offset = 2 * (kmer_id + string_id * (k-1))
            // The absolute base position in the concatenated strings
            let abs_pos = self.result.kmer_id_in_string as usize
                + self.result.string_begin as usize;

            let next_abs_pos = if self.result.kmer_orientation > 0 {
                abs_pos + 1
            } else {
                abs_pos.wrapping_sub(1)
            };

            // Read expected k-mer from string data at the next position
            let expected_kmer: Kmer<K> = dict.spss().decode_kmer_at(next_abs_pos);

            if expected_kmer.as_u64() == kmer.as_u64()
                || expected_kmer.as_u64() == kmer_rc.as_u64()
            {
                // Successfully extended!
                self.num_extensions += 1;
                let delta = self.result.kmer_orientation as i64;
                self.result.kmer_id = (self.result.kmer_id as i64 + delta) as u64;
                self.result.kmer_id_in_string =
                    (self.result.kmer_id_in_string as i64 + delta) as u64;
                self.result.kmer_offset =
                    (self.result.kmer_offset as i64 + delta) as u64;
                self.remaining_string_bases -= 1;
                return;
            }
        }
        
        // Extension failed, do a full search
        self.seed(Some(dict));
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
        assert_eq!(query.k, 31);
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
        
        assert!(query.is_valid_kmer("ACGTACGTACGTACGTACGTACGTACGTACG")); // 31 bases
        assert!(!query.is_valid_kmer("ACGT")); // Too short
        assert!(!query.is_valid_kmer("ACGTACGTACGTACGTACGTACGTACGTACGN")); // Invalid base
        
        assert!(query.is_valid_base(b'A'));
        assert!(query.is_valid_base(b'a'));
        assert!(!query.is_valid_base(b'N'));
    }

    #[test]
    fn test_streaming_query_lookup_invalid() {
        let mut query: StreamingQuery<15> = StreamingQuery::new(15, 7, true);
        
        // Invalid: too short
        let result = query.lookup("ACGT");
        assert!(!result.is_found());
        assert_eq!(query.num_invalid_lookups(), 1);
        
        // Invalid: has 'N'
        query.reset();
        let result = query.lookup("ACGTACGTACGTACN");
        assert!(!result.is_found());
        assert_eq!(query.num_invalid_lookups(), 2);
    }

    #[test]
    fn test_streaming_query_incremental_update() {
        let mut query: StreamingQuery<9> = StreamingQuery::new(9, 5, false);
        
        // First lookup
        let _result1 = query.lookup("ACGTACGTA");
        assert!(!query.start); // No longer in start state
        
        // Second lookup (sliding by 1)
        let _result2 = query.lookup("CGTACGTAC");
        
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
    pub fn lookup(&mut self, kmer_str: &str) -> LookupResult {
        // Perform streaming lookup with dictionary integration
        self.query.lookup_with_dict(kmer_str, self.dict)
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
