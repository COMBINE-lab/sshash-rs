//! Experimental split-phase streaming queries for a single-threaded read scheduler.
//! The ordinary streaming path is unchanged. These preparation/memo/extension
//! routines mirror that path and are checked against it field-for-field. This
//! duplication isolates the scheduling experiment from baseline code generation.
use super::*;
use crate::Dictionary;

/// One independent read's streaming state, with dictionary seeds deferred.
///
/// `prepare` handles rolling minimizers, extensions, and existing negative,
/// singleton, and decoded-bucket memos. `None` suspends at a genuine dictionary
/// seed; resolve it with `complete_scalar` or `Dictionary::complete_seed_batch`.
/// Consume the returned result before advancing this read. Different reads may
/// advance independently. No worker threads, allocation per query, or reordering
/// within a read is required. Reset between reads and before nonconsecutive
/// windows, exactly as for the original streaming engine.
///
/// A query is tied to its dictionary for its entire lifetime, preventing seeds
/// from accidentally being completed against a different index. The dictionary
/// must remain immutable, including its derived filter, during this lifetime.
pub struct DeferredQuery<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    pub(crate) dict: &'a Dictionary,
    state: StreamingQuery<K>,
    pending: bool,
}
impl<'a, const K: usize> DeferredQuery<'a, K>
where
    Kmer<K>: KmerBits,
{
    /// Allocate one stream state, with the index's minimizer hash seed.
    pub fn new(dict: &'a Dictionary) -> Self {
        assert_eq!(K, dict.k(), "query/index k mismatch");
        Self {
            dict,
            state: StreamingQuery::with_seed(K, dict.m(), dict.seed()),
            pending: false,
        }
    }
    /// Discard this read's state, including any outstanding seed.
    pub fn reset(&mut self) {
        self.state.reset();
        self.pending = false;
    }
    /// Advance one window. Panics if a previous seed is still pending.
    /// Malformed lengths/bases produce the ordinary invalid-query result.
    #[inline(always)]
    pub fn prepare(&mut self, bytes: &[u8]) -> Option<LookupResult> {
        assert!(!self.pending, "complete the pending seed before advancing");
        if bytes.len() != K {
            self.state.num_invalid += 1;
            self.reset();
            return Some(self.state.result);
        }
        let result = self.state.prepare_deferred(bytes, self.dict);
        self.pending = result.is_none();
        result
    }
    /// Whether this read is suspended on a dictionary seed.
    #[inline(always)]
    pub fn is_pending(&self) -> bool {
        self.pending
    }
    /// Resolve one pending seed through the original scalar dictionary path.
    /// This is also the scheduling-only ablation for the batched resolver.
    pub fn complete_scalar(&mut self) -> LookupResult {
        assert!(self.pending, "no pending seed");
        let s = &mut self.state;
        let result = self.dict.lookup_with_minimizer(
            &s.kmer,
            &s.kmer_rc,
            s.curr_mini_info,
            Some(&mut s.bucket_cache),
        );
        self.finish(result);
        result
    }
    /// Complete a seed using an independent selected-minimizer absence filter.
    ///
    /// This research path preserves membership and complete positive results,
    /// but may strengthen a negative `minimizer_found` hint from true to false
    /// when the legacy MPHF would route an absent key into a heavy bucket. It
    /// therefore does **not** promise field-identical negative results. Ordinary
    /// scalar and batched completion retain their original full-field contract.
    ///
    /// For correct membership, `may_contain` must return true for every canonical
    /// minimizer selected anywhere in this dictionary. False positives are safe;
    /// false negatives suppress real hits, including subsequent windows sharing
    /// that minimizer. Absence of an arbitrary suffix q-mer is NOT such a proof.
    /// The predicate is called once, only for an outstanding uncached seed.
    ///
    /// A certified rejection installs the same negative-minimizer memo as an
    /// exact missing-minimizer lookup, retaining rolling state across windows.
    /// It invalidates prior bucket bounds/anchors and never installs a positive offset.
    /// Panics if no seed is pending. Available only with `experimental-query`.
    #[inline]
    pub fn complete_with_minimizer_filter(
        &mut self,
        may_contain: impl FnOnce(u64) -> bool,
    ) -> LookupResult {
        assert!(self.pending, "no pending seed");
        if may_contain(self.state.curr_mini_info.value) {
            return self.complete_scalar();
        }
        self.state.bucket_cache.reset();
        let mut result = LookupResult::not_found();
        result.minimizer_found = false;
        self.finish(result);
        result
    }
    /// Direct streaming control for independent minimizer filtering.
    ///
    /// This runs a complete window without a prepare/complete suspension or an
    /// Option-valued intermediate result. Rolling/extension behavior mirrors
    /// the original stream; filtering occurs only in the outlined seed path.
    /// It has the same no-false-negative predicate requirement and strengthened
    /// negative-hint contract as [`Self::complete_with_minimizer_filter`].
    /// Reset before discontinuities; panics if a deferred seed is pending.
    /// Malformed lengths/bases return the ordinary invalid-query result.
    #[inline(always)]
    pub fn lookup_filtered(
        &mut self,
        bytes: &[u8],
        may_contain: impl FnOnce(u64) -> bool,
    ) -> LookupResult {
        assert!(!self.pending, "complete the pending seed before advancing");
        if bytes.len() != K {
            self.state.num_invalid += 1;
            self.reset();
            return self.state.result;
        }
        self.state
            .lookup_filtered_internal(bytes, self.dict, may_contain)
    }
    /// Most recent completed lookup. Panics while a seed is outstanding.
    #[inline(always)]
    pub fn result(&self) -> LookupResult {
        assert!(!self.pending);
        self.state.result
    }
    /// Read-only encoded pending seed and its rolling minimizer.
    /// Used by diagnostics and research filters; does not complete or advance it.
    pub fn pending_seed(&self) -> Option<(Kmer<K>, MinimizerInfo)> {
        self.pending
            .then_some((self.state.kmer, self.state.curr_mini_info))
    }
    #[inline(always)]
    pub(crate) fn parts(&self) -> (Kmer<K>, Kmer<K>, MinimizerInfo) {
        (
            self.state.kmer,
            self.state.kmer_rc,
            self.state.curr_mini_info,
        )
    }
    #[inline(always)]
    pub(crate) fn cache_mut(&mut self) -> &mut crate::dictionary::BucketCache {
        &mut self.state.bucket_cache
    }
    pub(crate) fn finish(&mut self, result: LookupResult) {
        let s = &mut self.state;
        s.result = result;
        s.last_seed_mini_info = s.curr_mini_info;
        s.last_seed_positive = result.is_found();
        if result.is_found() {
            s.num_searches += 1;
            s.begin_extension();
        } else {
            s.num_negative += 1;
        }
        self.pending = false;
    }
}
impl<const K: usize> StreamingQuery<K>
where
    Kmer<K>: KmerBits,
{
    /// Research copy of the original hot path: filter work stays outlined.
    #[inline(always)]
    fn lookup_filtered_internal(
        &mut self,
        kmer_bytes: &[u8],
        dict: &Dictionary,
        may_contain: impl FnOnce(u64) -> bool,
    ) -> LookupResult {
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

        if self.remaining_string_bases == 0 || !self.try_extend_deferred(dict) {
            self.seed_with_filter(dict, may_contain);
        }

        // 4. Update state
        self.start = false;

        self.result
    }
    /// Original memo rules with one predicate before uncached MPHF resolution.
    #[inline(never)]
    fn seed_with_filter(&mut self, dict: &Dictionary, may_contain: impl FnOnce(u64) -> bool) {
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

        // Repeated heavy buckets already have resolved bounds; do not add a
        // filter probe there. A new minimizer alone reaches the predicate.
        self.result = if curr.value == self.last_seed_mini_info.value || may_contain(curr.value) {
            dict.lookup_with_minimizer(&kmer, &kmer_rc, curr, Some(&mut self.bucket_cache))
        } else {
            self.bucket_cache.reset();
            let mut result = LookupResult::not_found();
            result.minimizer_found = false;
            result
        };
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

    /// Shared-state preparation; no MPHF work occurs before suspension.
    #[inline(always)]
    fn prepare_deferred(&mut self, kmer_bytes: &[u8], dict: &Dictionary) -> Option<LookupResult> {
        if self.start {
            if !self.is_valid_kmer_bytes(kmer_bytes) {
                self.num_invalid += 1;
                self.reset();
                return Some(self.result);
            }
            self.kmer = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
            self.kmer_rc = self.kmer.reverse_complement();
        } else {
            let last_base = unsafe { *kmer_bytes.get_unchecked(K - 1) };
            if !self.is_valid_base(last_base) {
                self.num_invalid += 1;
                self.reset();
                return Some(self.result);
            }
            let encoded = crate::encoding::encode_base_unchecked(last_base);
            self.kmer = self.kmer.roll_right_base(encoded);
            let complement = crate::encoding::complement_base(encoded);
            self.kmer_rc = self.kmer_rc.append_base(complement);
        }

        self.curr_mini_info = self.minimizer_it.next(&self.kmer, &self.kmer_rc);

        self.start = false;
        if self.remaining_string_bases != 0 && self.try_extend_deferred(dict) {
            return Some(self.result);
        }
        if self.try_local_seed(dict) {
            return Some(self.result);
        }
        None
    }
    /// Preserve all three existing memo rules, including conservative miss hints.
    #[inline(never)]
    fn try_local_seed(&mut self, dict: &Dictionary) -> bool {
        self.remaining_string_bases = 0;
        let kmer = self.kmer;
        let kmer_rc = self.kmer_rc;
        let curr = self.curr_mini_info;

        if curr.value == self.last_seed_mini_info.value {
            // (a) minimizer absent at the last seed.
            if !self.result.minimizer_found {
                debug_assert_eq!(self.result.kmer_id, u64::MAX);
                self.num_negative += 1;
                return true;
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
                return true;
            }

            // (c) cached locate set: verify directly against the text.
            if self.bucket_cache.size != 0 {
                let n = self.bucket_cache.size;
                let positions = self.bucket_cache.positions;
                self.result = dict.lookup_at_positions(&positions[..n], &kmer, &kmer_rc, curr);
                self.num_bucket_cache_hits += 1;
                if self.result.kmer_id == u64::MAX {
                    self.num_negative += 1;
                    return true;
                }
                // Keep the singleton shortcut anchored to the latest match.
                self.last_seed_mini_info = curr;
                self.last_seed_positive = true;
                self.num_searches += 1;
                self.begin_extension();
                return true;
            }

            // A repeated heavy bucket already has cached bounds. Preserve that
            // reuse locally; only uncached resolutions belong in the batch.
            self.seed_with_dict(dict);
            return true;
        }

        false
    }
    /// Original extension logic, including corrected absolute offsets; failure
    /// returns to the scheduler instead of immediately issuing a dictionary seed.
    #[inline(always)]
    fn try_extend_deferred(&mut self, dict: &Dictionary) -> bool {
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
                self.result.kmer_offset += 1;
                self.remaining_string_bases -= 1;
                return true;
            }
        } else {
            self.buf_bit_pos -= 1;
            let spss_base = Self::get_base(dict, self.buf_bit_pos as usize);
            if spss_base == u64::from(self.kmer_rc.get_base(0)) {
                self.num_extensions += 1;
                self.result.kmer_id -= 1;
                self.result.kmer_id_in_string -= 1;
                self.result.kmer_offset -= 1;
                self.remaining_string_bases -= 1;
                return true;
            }
        }
        false
    }
}
