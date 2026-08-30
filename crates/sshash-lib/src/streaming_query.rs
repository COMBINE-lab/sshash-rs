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
    m: usize,

    // K-mer state (zero-initialized; `start` flag tracks validity)
    start: bool,
    kmer: Kmer<K>,
    kmer_rc: Kmer<K>,

    // Minimizer state: one iterator, one order, both strands (the unified
    // canonical scheme is mirror-equivariant, so no RC iterator is needed).
    minimizer_it: MinimizerIterator,
    curr_mini_info: MinimizerInfo,

    // String extension state
    remaining_string_bases: u64,
    buf_bit_pos: u64, // current k-mer position in strings data (base units)
    buf: u64,         // buffered bases ahead of the current new-base position
    buf_avail: u32,   // number of valid bits remaining in buf

    // Result state
    result: LookupResult,

    /// Last real-seed state, for the same-minimizer memos of `seed_with_dict`:
    /// the minimizer of the last dictionary seed (with, for the singleton
    /// shortcut, the stream position of its occurrence) and whether a
    /// positive match anchors it. The sentinel value `u64::MAX` can never be
    /// a real minimizer (m <= 31 means values fit 62 bits).
    last_seed_mini_info: MinimizerInfo,
    last_seed_positive: bool,

    // Performance counters
    num_searches: u64,
    num_extensions: u64,
    num_invalid: u64,
    num_negative: u64,
    num_skipped_singleton_lookups: u64,
    num_bucket_cache_hits: u64,

    /// Memoized bucket bounds + decoded locate set for the current minimizer
    /// (see `BucketCache`).
    bucket_cache: crate::dictionary::BucketCache,
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
    /// * `_canonical` - Ignored since 0.7.0: the index is always canonical.
    ///   Kept so existing callers compile.
    pub fn new(k: usize, m: usize, _canonical: bool) -> Self {
        Self::with_seed(k, m, 1)
    }

    /// Create a new streaming query engine with an explicit hasher seed
    /// (must match the seed the queried dictionary was built with).
    pub fn with_seed(k: usize, m: usize, seed: u64) -> Self {
        assert_eq!(k, K, "k parameter must match const generic K");

        let dummy_mini = MinimizerInfo::new(u64::MAX, 0, 0);

        Self {
            _k: k,
            m,
            start: true,
            kmer: Kmer::new(<Kmer<K> as KmerBits>::from_u64(0)),
            kmer_rc: Kmer::new(<Kmer<K> as KmerBits>::from_u64(0)),
            minimizer_it: MinimizerIterator::with_seed(k, m, seed),
            curr_mini_info: dummy_mini,
            remaining_string_bases: 0,
            buf_bit_pos: 0,
            buf: 0,
            buf_avail: 0,
            result: LookupResult::not_found(),
            last_seed_mini_info: dummy_mini,
            last_seed_positive: false,
            num_searches: 0,
            num_extensions: 0,
            num_invalid: 0,
            num_negative: 0,
            num_skipped_singleton_lookups: 0,
            num_bucket_cache_hits: 0,
            bucket_cache: Default::default(),
        }
    }

    /// Reset the query state (call this when starting a new sequence)
    pub fn reset(&mut self) {
        self.start = true;
        self.remaining_string_bases = 0;
        self.result = LookupResult::not_found();
        self.minimizer_it.set_position(0);
        // The stream coordinate frame restarts, so the position-anchored
        // singleton memo must not survive; the bucket cache's contents are
        // pure functions of the minimizer value but are dropped with it.
        self.last_seed_mini_info = MinimizerInfo::new(u64::MAX, 0, 0);
        self.last_seed_positive = false;
        self.bucket_cache.reset();
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

        self.curr_mini_info = self.minimizer_it.next(&self.kmer, &self.kmer_rc);

        // 3. Compute result (either extend or search)
        if self.remaining_string_bases == 0 {
            self.seed_with_dict(dict);
        } else {
            self.try_extend(dict);
        }

        // 4. Update state
        self.start = false;

        // Ground-truth check (debug builds): every streaming result — memo
        // answers included — must equal the corresponding point lookup on the
        // fields C++ `equal_lookup_result` compares.
        #[cfg(debug_assertions)]
        {
            let point = dict.query::<K>(&self.kmer);
            debug_assert_eq!(self.result.kmer_id, point.kmer_id, "streaming != point (kmer_id)");
            if self.result.is_found() {
                debug_assert_eq!(self.result.kmer_id_in_string, point.kmer_id_in_string);
                debug_assert_eq!(self.result.kmer_orientation, point.kmer_orientation);
                debug_assert_eq!(self.result.string_id, point.string_id);
                debug_assert_eq!(self.result.string_begin, point.string_begin);
                debug_assert_eq!(self.result.string_end, point.string_end);
            }
        }

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

        self.curr_mini_info = self.minimizer_it.next(&self.kmer, &self.kmer_rc);

        // 3. Compute result (either extend or search)
        match dict_opt {
            Some(dict) => {
                if self.remaining_string_bases == 0 {
                    self.seed_with_dict(dict);
                } else {
                    self.try_extend(dict);
                }
            }
            None => {
                // MVP path without a dictionary: nothing to search.
                self.remaining_string_bases = 0;
                self.result = LookupResult::not_found();
                self.num_negative += 1;
            }
        }

        // 4. Update state
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

    /// Perform a full search (seed) for the current k-mer.
    ///
    /// Called when we can't extend within the current string. Outlined: it
    /// inlines into `lookup` at two call sites, and its body (grown by the
    /// same-minimizer memos) measurably slows down streams that rarely seed
    /// when left inline (observed in C++, commit c22c897).
    ///
    /// When the seeded k-mer carries the same minimizer as the last seed,
    /// three memos answer without (or with reduced) dictionary work:
    ///
    /// (a) the minimizer was absent at the last seed: any k-mer having the
    ///     same minimizer is surely absent as well — negative for free;
    /// (b) the last seed matched via a singleton bucket at stream position L
    ///     and the current k-mer still carries the very same minimizer
    ///     occurrence (same stream position: the scheme is forward, so a
    ///     sampled position, once abandoned, is never re-selected). The only
    ///     locus the k-mer could occupy is the one the just-failed extension
    ///     step already rejected, or it lies beyond the string boundary; the
    ///     mirror locus hosts the reverse complement of the occurrence, so
    ///     it is excluded whenever the minimizer is not its own reverse
    ///     complement (self-RC minimizers, possible for even m only, fall
    ///     through) — negative with no text access at all;
    /// (c) the decoded locate set of the bucket is cached: verify the
    ///     candidate loci directly against the text, with no MPHF access and
    ///     no offset decoding.
    ///
    /// Heavy buckets fall through to the full lookup. Error-free positive
    /// streams never enter the memos and are unaffected.
    #[inline(never)]
    fn seed_with_dict(&mut self, dict: &crate::dictionary::Dictionary) {
        self.remaining_string_bases = 0;

        let kmer = self.kmer;
        let kmer_rc = self.kmer_rc;
        let curr = self.curr_mini_info;

        if curr.value == self.last_seed_mini_info.value {
            // (a) minimizer absent at the last seed.
            if !self.result.minimizer_found {
                debug_assert_eq!(self.result.kmer_id, u64::MAX);
                self.num_negative += 1;
                return;
            }

            // (b) singleton bucket, same minimizer occurrence, positive anchor.
            if self.last_seed_positive
                && self.bucket_cache.size == 1
                && curr.position == self.last_seed_mini_info.position
                && !crate::minimizer::is_self_rc_mmer(curr.value, self.m)
            {
                self.result = LookupResult::not_found();
                self.num_negative += 1;
                self.num_skipped_singleton_lookups += 1;
                return;
            }

            // (c) cached locate set: verify directly against the text.
            if self.bucket_cache.size != 0 {
                let n = self.bucket_cache.size;
                let positions = self.bucket_cache.positions;
                self.result = dict.lookup_at_positions(&positions[..n], &kmer, &kmer_rc, curr);
                self.num_bucket_cache_hits += 1;
                if self.result.kmer_id == u64::MAX {
                    self.num_negative += 1;
                    return;
                }
                // Keep the singleton shortcut anchored to the latest match.
                self.last_seed_mini_info = curr;
                self.last_seed_positive = true;
                self.num_searches += 1;
                self.begin_extension();
                return;
            }

            // Heavy bucket: fall through to the full lookup.
        }

        self.result =
            dict.lookup_with_minimizer(&kmer, &kmer_rc, curr, Some(&mut self.bucket_cache));
        self.last_seed_mini_info = curr;
        self.last_seed_positive = self.result.kmer_id != u64::MAX;

        if self.result.kmer_id == u64::MAX {
            self.num_negative += 1;
            return;
        }

        debug_assert!(self.result.minimizer_found);
        self.num_searches += 1;
        self.begin_extension();
    }

    /// Set up single-base extension state from a fresh positive result.
    #[inline(always)]
    fn begin_extension(&mut self) {
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

    /// Try to extend within the current string using single-base comparison
    /// with buffered SPSS reads.
    ///
    /// Combines two optimizations:
    /// 1. Buffer amortization: 64-bit reads cover 32 bases, shift-extracted on
    ///    successive forward extensions (avoids memory access every other call).
    /// 2. Single-base comparison: only check the one new base instead of loading
    ///    and comparing a full k-mer. Correct because consecutive k-mers in a
    ///    string share K-1 bases.
    #[inline(always)]
    fn try_extend(&mut self, dict: &crate::dictionary::Dictionary) {
        if self.result.kmer_orientation > 0 {
            self.buf_bit_pos += 1;
            let new_base_pos = self.buf_bit_pos as usize + K - 1;
            let spss_base = if self.buf_avail >= 2 {
                self.buf_avail -= 2;
                let b = self.buf & 3;
                self.buf >>= 2;
                b
            } else {
                self.load_forward_base(dict, new_base_pos)
            };
            // The newest base sits at bits [2(K-1), 2K): extract it in storage
            // width. Going through to_u64 truncates u128 storage and the
            // 2*(K-1) shift overflows u64 for K >= 33.
            if spss_base == u64::from(self.kmer.get_base(K - 1)) {
                self.num_extensions += 1;
                self.result.kmer_id += 1;
                self.result.kmer_id_in_string += 1;
                self.remaining_string_bases -= 1;
                return;
            }
        } else {
            self.buf_bit_pos -= 1;
            let spss_base = Self::get_base(dict, self.buf_bit_pos as usize);
            if spss_base == u64::from(self.kmer_rc.get_base(0)) {
                self.num_extensions += 1;
                self.result.kmer_id -= 1;
                self.result.kmer_id_in_string -= 1;
                self.remaining_string_bases -= 1;
                return;
            }
        }
        self.seed_with_dict(dict);
    }

    /// Load a single base from SPSS at `abs_pos` and refill the forward buffer
    /// with subsequent bases for shift-extraction on the next calls.
    #[inline(always)]
    fn load_forward_base(&mut self, dict: &crate::dictionary::Dictionary, abs_pos: usize) -> u64 {
        let strings = dict.spss().strings_data();
        let byte_offset = abs_pos / 4;
        let bit_shift = (abs_pos % 4) * 2;
        let raw = unsafe {
            std::ptr::read_unaligned(strings.as_ptr().add(byte_offset) as *const u64)
        };
        let word = raw >> bit_shift;
        let base = word & 3;
        self.buf = word >> 2;
        self.buf_avail = if bit_shift == 0 { 62 } else { 64 - bit_shift as u32 - 2 };
        base
    }

    /// Extract a single 2-bit base from the SPSS at the given absolute position.
    #[inline(always)]
    fn get_base(dict: &crate::dictionary::Dictionary, abs_pos: usize) -> u64 {
        let strings = dict.spss().strings_data();
        let byte_offset = abs_pos / 4;
        let bit_shift = (abs_pos % 4) * 2;
        let byte = unsafe { *strings.as_ptr().add(byte_offset) };
        ((byte >> bit_shift) & 3) as u64
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

    /// Negatives answered by the singleton same-occurrence memo (no text access).
    pub fn num_skipped_singleton_lookups(&self) -> u64 {
        self.num_skipped_singleton_lookups
    }

    /// Seeds answered from the cached locate set (no MPHF/offset decoding).
    pub fn num_bucket_cache_hits(&self) -> u64 {
        self.num_bucket_cache_hits
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
        assert_eq!(query.m, 13);
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
        Self {
            dict,
            query: StreamingQuery::with_seed(dict.k(), dict.m(), dict.seed()),
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
    
    /// Negatives answered by the singleton same-occurrence memo.
    pub fn num_skipped_singleton_lookups(&self) -> u64 {
        self.query.num_skipped_singleton_lookups()
    }

    /// Seeds answered from the cached locate set.
    pub fn num_bucket_cache_hits(&self) -> u64 {
        self.query.num_bucket_cache_hits()
    }

    /// Get statistics
    pub fn stats(&self) -> StreamingQueryStats {
        StreamingQueryStats {
            num_searches: self.query.num_searches(),
            num_extensions: self.query.num_extensions(),
            num_invalid: self.query.num_invalid_lookups(),
            num_negative: self.query.num_negative_lookups(),
            num_skipped_singleton_lookups: self.query.num_skipped_singleton_lookups(),
            num_bucket_cache_hits: self.query.num_bucket_cache_hits(),
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
    /// Negatives answered by the singleton same-occurrence memo
    pub num_skipped_singleton_lookups: u64,
    /// Seeds answered from the cached locate set (no MPHF/offset decoding)
    pub num_bucket_cache_hits: u64,
}
