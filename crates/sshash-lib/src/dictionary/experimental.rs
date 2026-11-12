//! Query experiments enabled by the `experimental-query` Cargo feature.
//!
//! Both paths ultimately use `Dictionary::lookup_in_bucket`, preserving candidate
//! order, orientation, IDs, and the conservative `minimizer_found` hint. The
//! fingerprint table is derived from the loaded index; neither experiment changes
//! index construction or the `.ssi` / `.ssi.mphf` representation.
//!
//! This module deliberately does not replace the streaming engine with point
//! batches: its rolling minimizer, extension, and same-minimizer caches remain
//! valuable. A caller must batch independent seeds rather than reorder dependent
//! windows within a read.
use super::*;
use crate::streaming_query::LookupResult;
use sux::utils::prefetch_index;

#[inline]
/// Seven-bit tag of a canonical minimizer, independent of its MPHF address.
///
/// The SplitMix64 finalizer mixes all bases before masking. This is a deterministic
/// non-cryptographic mixer, not a guaranteed-uniform random fingerprint: the
/// often-useful 1/128 collision estimate is a model, not a correctness assumption.
/// Values 0..=127 are tags; 128 is reserved for unconditional bypass.
pub(super) fn fingerprint(value: u64) -> u8 {
    // Mix all bases; using low bits directly biases selected minimizers.
    let mut x = value;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    ((x ^ (x >> 31)) & 127) as u8
}

impl Dictionary {
    /// Enable a derived minimizer-absence filter for point, batch, and streaming queries.
    ///
    /// Available only with the `experimental-query` feature. New and loaded
    /// dictionaries start with this filter disabled. Call this before borrowing
    /// the dictionary for queries; subsequent queries only read the table.
    ///
    /// # Representation and exactness
    ///
    /// An MPHF-indexed byte array stores a seven-bit tag of each light bucket's
    /// canonical minimizer. A tag mismatch proves that the selected minimizer
    /// is absent from the minimizer map, and skips bucket-bound/anchor/text reads.
    /// A match is only a *possible* hit and still receives full verification.
    /// Fingerprint collisions therefore cannot create incorrect lookup results.
    ///
    /// Heavy buckets (more than 64 anchors) and empty buckets use a bypass marker.
    /// In particular, heavy-bucket misses must keep the original conservative
    /// `minimizer_found = true` hint. This filter preserves every result field,
    /// not just membership. It does not prove that the minimizer's bases occur
    /// nowhere in the reference: some occurrences may never be selected minima.
    ///
    /// # Cost and lifecycle
    ///
    /// Construction is sequential, visits every bucket, and reads one anchor and
    /// its sequence for each nonempty light bucket. Retained payload is one byte
    /// per bucket and is included in [`Self::num_bits`]. The return value is that
    /// payload in bytes, excluding allocator/struct overhead. Calling this again
    /// rebuilds the table; the old allocation remains live until replacement.
    ///
    /// The table is not serialized: saving with the filter enabled produces the
    /// same index files, and loading them leaves it disabled. Disable it with
    /// [`Self::disable_bucket_fingerprints`]. The fixed seven-bit width is the
    /// measured baseline, not a configurable or proven-optimal layout.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// # #[cfg(feature = "experimental-query")]
    /// # {
    /// use sshash_lib::{Dictionary, Kmer};
    /// let mut dict = Dictionary::load("index")?;
    /// assert_eq!(dict.k(), 31);
    /// let key = Kmer::<31>::from_str("ACGTACGTACGTACGTACGTACGTACGTACG")?;
    /// let expected = dict.query(&key);
    /// let bytes = dict.enable_bucket_fingerprints();
    /// assert_eq!(bytes as u64, dict.num_minimizers()); // builder-produced index
    /// assert_eq!(dict.query(&key), expected);
    /// dict.disable_bucket_fingerprints();
    /// # }
    /// # Ok(())
    /// # }
    /// ```
    pub fn enable_bucket_fingerprints(&mut self) -> usize {
        let tags: Vec<u8> = (0..self.index.num_buckets())
            .map(|id| {
                let (begin, end) = self.index.locate_bucket(id);
                if begin == end || ceil_log2((end - begin) as u64) > MIN_L {
                    return 128;
                }
                let pos = self.index.offsets.index_value(begin);
                let value = self.spss.decode_mmer_at(pos, self.m);
                fingerprint(value.min(crate::minimizer::reverse_complement_mmer(value, self.m)))
            })
            .collect();
        let bytes = tags.len();
        self.bucket_fingerprints = Some(tags.into_boxed_slice());
        bytes
    }

    /// Disable the derived filter and release its allocation.
    ///
    /// Available only with `experimental-query`. This is an idempotent operation:
    /// calling it on an already-disabled dictionary has no effect. It changes
    /// neither the dictionary's contents nor its serialized representation.
    pub fn disable_bucket_fingerprints(&mut self) {
        self.bucket_fingerprints = None;
    }

    /// Query `B` independent encoded k-mers, returning exact results in input order.
    ///
    /// Available only with `experimental-query`. Each result equals an individual
    /// [`Self::query`] of the corresponding key, including negative-result hints.
    /// Duplicates and both orientations are allowed; `B = 0` returns an empty
    /// array. All `B` slots are queried, so a caller padding a final partial batch
    /// must discard the padded results. Intermediates occupy `O(B)` stack space;
    /// there is no per-call heap allocation. Large `B` can exhaust the stack.
    ///
    /// Queries progress through separate minimizer, MPHF, bounds, and verification
    /// stages. On a single-partition MPHF, PHastR prefetches first-level seeds;
    /// multiple partitions retain scalar routing. `PREFETCH = true` additionally
    /// hints Elias–Fano select inventories/low bits, packed anchor offsets, and
    /// the first light-bucket sequence region. Heavy buckets and later candidates
    /// retain ordinary verification; this is not a fully pipelined locator.
    ///
    /// `PREFETCH = false` is an ablation of those *downstream* hints, not of
    /// PHastR's batching/prefetching. Hints use sux's architecture-dependent
    /// implementation and may be no-ops on unsupported targets. They do not
    /// affect correctness. Batch sizes 16–64 are useful starting points to
    /// measure, not a universal optimum. This API does not reuse the streaming
    /// engine's rolling minimizer, extension, or memoization state.
    ///
    /// # Panics
    ///
    /// Panics if `K != self.k()`, including for an empty batch. As for scalar
    /// queries, keys and dictionary must satisfy the ordinary k-mer/index
    /// invariants; the `KmerBits` bound permits odd `K` in 3..=63.
    ///
    /// # Example
    ///
    /// ```no_run
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// # #[cfg(feature = "experimental-query")]
    /// # {
    /// use sshash_lib::{Dictionary, Kmer};
    /// let dict = Dictionary::load("index")?;
    /// let key = Kmer::<31>::from_str("ACGTACGTACGTACGTACGTACGTACGTACG")?;
    /// let keys = [key; 32]; // in practice, seeds from independent reads
    /// let results = dict.query_batch::<31, 32, true>(&keys);
    /// for (key, result) in keys.iter().zip(results) {
    ///     assert_eq!(dict.query(key), result);
    /// }
    /// # }
    /// # Ok(())
    /// # }
    /// ```
    pub fn query_batch<const K: usize, const B: usize, const PREFETCH: bool>(
        &self,
        kmers: &[Kmer<K>; B],
    ) -> [LookupResult; B]
    where
        Kmer<K>: KmerBits,
    {
        assert_eq!(K, self.k);
        let rcs = kmers.map(|k| k.reverse_complement());
        let minis: [_; B] = std::array::from_fn(|i| self.extract_minimizer(&kmers[i], &rcs[i]));
        let keys = minis.map(|m| m.value);
        let ids = self.control_map.mphf_ref().map(|m| m.get_batch(&keys));
        let Some(ids) = ids else {
            return std::array::from_fn(|i| self.query(&kmers[i]));
        };
        let active: [_; B] = std::array::from_fn(|i| {
            ids[i] < self.control_map.num_minimizers() as usize
                && self.bucket_may_contain(ids[i], keys[i])
        });
        if PREFETCH {
            for i in 0..B {
                if active[i] {
                    self.index
                        .num_super_kmers_before_bucket
                        .prefetch_access(ids[i]);
                }
            }
        }
        let bounds: [_; B] = std::array::from_fn(|i| {
            if !active[i] {
                return None;
            }
            let (begin, end) = self.index.locate_bucket(ids[i]);
            if PREFETCH && begin < end {
                let word = begin * self.index.offsets.bit_width() / 64;
                prefetch_index(self.index.offsets.as_slice(), word);
                prefetch_index(self.index.offsets.as_slice(), word + 1);
            }
            Some((begin, end))
        });
        if PREFETCH {
            for &(begin, end) in bounds.iter().flatten() {
                if begin == end || ceil_log2((end - begin) as u64) > MIN_L {
                    continue;
                }
                let pos = self.index.offsets.index_value(begin);
                let text = self.spss.strings_data();
                prefetch_index(text, pos / 4);
                // The two candidate starts are at most K bases before the
                // anchor. Cover both ends even across a cache-line boundary.
                prefetch_index(text, pos.saturating_sub(K) / 4);
                prefetch_index(text, pos / 4 + 16);
            }
        }
        std::array::from_fn(|i| {
            self.lookup_in_bucket(&kmers[i], &rcs[i], minis[i], bounds[i], None)
        })
    }
}
