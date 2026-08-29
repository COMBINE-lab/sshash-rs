//! Register-blocked Bloom filter used as a negative prefilter for
//! [`TinyDictionary`](crate::TinyDictionary) lookups.
//!
//! # Why
//!
//! For a probe-panel reference (10x Flex and friends) the mapped read is far
//! longer than the region that can possibly match: with a 50 bp probe set,
//! k=23 and a 150 bp read, ~79% of k-mer positions are *guaranteed* absent
//! from the dictionary (measured, not estimated). Every one of those currently
//! pays a full hashbrown probe.
//!
//! Profiling that workload (AMD EPYC 9575F, 20M read pairs) attributes ~35% of
//! all cycles to the probe itself — `_mm_movemask_epi8` 22%, hashbrown's
//! `likely` 9%, `lowest_set_bit` 3.9% — at IPC 2.38 with only 5.7% backend
//! stalls. The probe is expensive because it *executes* a lot of work, not
//! because it stalls on memory, so the win here comes from replacing ~40-60
//! instructions with ~12 and one load, not from avoiding cache misses.
//!
//! # Correctness
//!
//! A Bloom filter has no false negatives, so [`BlockedBloom::might_contain`]
//! returning `false` is proof the key is absent. False positives simply fall
//! through to the real hashbrown probe, which then answers authoritatively.
//! Mapping results are therefore bit-identical with the filter enabled or
//! disabled — the filter can only make a miss cheaper, never change an answer.
//!
//! The filter is derived from the key set at construction/load time and is not
//! serialized, so `.tdct` files stay byte-compatible in both directions.

/// Bits of filter per key.
///
/// 8 bits/key with 4 bits set per 64-bit block puts the false-positive rate
/// under ~1%, so essentially every one of the ~79% guaranteed-absent lookups is
/// answered from the filter alone. For a 1.5M-k-mer probe panel that is a 2 MiB
/// filter; 16 bits/key (4 MiB, ~0.24% FPR) measured very slightly *slower* on
/// the Flex workload — the lower FPR does not pay for the extra cache
/// footprint, so the smaller filter wins.
const BITS_PER_KEY: usize = 8;

/// Register-blocked Bloom filter: one 64-bit block per key, all `k` bits of a
/// key land in the same block, so a query touches exactly one cache line.
#[derive(Debug, Clone)]
pub(crate) struct BlockedBloom {
    /// Power-of-two number of 64-bit blocks.
    blocks: Box<[u64]>,
    /// `blocks.len() - 1`, for index masking.
    block_mask: u64,
}

/// Mix a canonical k-mer's bits into a well-distributed 64-bit hash.
///
/// Keys are 2-bit-packed k-mers — highly structured, with the low bits varying
/// fastest — so they cannot be used as a hash directly. One 64x64->128 multiply
/// folded down is enough to spread them, and costs a single `mulx`.
#[inline(always)]
fn mix(key: u64) -> u64 {
    const GOLDEN: u128 = 0x9E37_79B9_7F4A_7C15;
    let p = (key as u128).wrapping_mul(GOLDEN);
    ((p >> 64) as u64) ^ (p as u64)
}

/// The 4 bits this key sets within its block, as a single mask.
#[inline(always)]
fn block_bits(h: u64) -> u64 {
    // 4 disjoint 6-bit fields from the low 24 bits select 4 bit positions.
    (1u64 << (h & 63))
        | (1u64 << ((h >> 6) & 63))
        | (1u64 << ((h >> 12) & 63))
        | (1u64 << ((h >> 18) & 63))
}

impl BlockedBloom {
    /// Build a filter over `keys`. `n_hint` sizes the table; passing the exact
    /// key count is best but any over-estimate is safe.
    pub(crate) fn build<I: IntoIterator<Item = u64>>(keys: I, n_hint: usize) -> Self {
        let want_blocks = (n_hint * BITS_PER_KEY).div_ceil(64).max(1);
        let num_blocks = want_blocks.next_power_of_two();
        let mut blocks = vec![0u64; num_blocks].into_boxed_slice();
        let block_mask = (num_blocks - 1) as u64;

        for key in keys {
            let h = mix(key);
            let idx = ((h >> 32) & block_mask) as usize;
            blocks[idx] |= block_bits(h);
        }

        Self { blocks, block_mask }
    }

    /// `false` means the key is definitely absent. `true` means "probably
    /// present — go ask the real table".
    #[inline(always)]
    pub(crate) fn might_contain(&self, key: u64) -> bool {
        let h = mix(key);
        let idx = ((h >> 32) & self.block_mask) as usize;
        // SAFETY-free: `block_mask` is `len - 1` for a power-of-two `len`, so
        // the index is always in bounds and the check folds away.
        let word = self.blocks[idx];
        let m = block_bits(h);
        (word & m) == m
    }

    /// Size of the filter in bytes (for logging / diagnostics).
    #[allow(dead_code)]
    pub(crate) fn size_bytes(&self) -> usize {
        self.blocks.len() * 8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_false_negatives() {
        // Every inserted key must be reported present. This is the property the
        // dictionary's correctness depends on.
        let keys: Vec<u64> = (0..200_000u64).map(|i| i.wrapping_mul(0x1234_5678_9ABC_DEF)).collect();
        let f = BlockedBloom::build(keys.iter().copied(), keys.len());
        for &k in &keys {
            assert!(f.might_contain(k), "false negative for key {k:#x}");
        }
    }

    #[test]
    fn kmer_like_keys_have_no_false_negatives() {
        // Structured, 2-bit-packed-looking keys: the realistic case.
        let keys: Vec<u64> = (0..100_000u64).map(|i| i * 4 + (i % 3)).collect();
        let f = BlockedBloom::build(keys.iter().copied(), keys.len());
        for &k in &keys {
            assert!(f.might_contain(k), "false negative for key {k:#x}");
        }
    }

    #[test]
    fn false_positive_rate_is_low() {
        let n = 1_500_000usize;
        let keys: Vec<u64> = (0..n as u64).map(|i| mix(i ^ 0xDEAD_BEEF)).collect();
        let f = BlockedBloom::build(keys.iter().copied(), n);

        let mut fp = 0;
        let trials = 200_000u64;
        for i in 0..trials {
            // Keys guaranteed not in the set.
            let probe = mix(i ^ 0xFEED_FACE) | 1 << 63;
            if f.might_contain(probe) {
                fp += 1;
            }
        }
        let rate = fp as f64 / trials as f64;
        assert!(rate < 0.02, "false positive rate too high: {rate}");
    }

    #[test]
    fn empty_is_safe() {
        let f = BlockedBloom::build(std::iter::empty(), 0);
        assert!(!f.might_contain(12345));
    }
}
