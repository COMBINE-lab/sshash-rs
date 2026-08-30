//! Minimizer extraction and iteration — the unified canonical scheme.
//!
//! Port of C++ SSHash v6.0.0 (commits 6c53ae9 + b1d0706): the minimizer of a
//! k-mer x is the locus i in [0, k-m] minimizing `h(kappa(i))`, where
//! `kappa(i) = min(mmer_i, rc(mmer_i))` is the canonical m-mer at locus i.
//!
//! Because kappa is invariant under reverse-complementation and the
//! kappa-sequence of rc(x) is the reversal of that of x, rc(x) selects the
//! mirrored locus and literally the same m-mer: x and rc(x) always share a
//! bucket (one probe per lookup) while the selection remains an argmin under a
//! single random order, keeping the plain forward density 2/(k-m+2). This is
//! what lets one indexing modality replace the old regular/canonical pair.
//!
//! Ties (two loci with the same minimum hash) must be broken in a
//! mirror-equivariant way or x and rc(x) would land in different buckets — a
//! correctness failure, not a density one. No leftmost rule can work: the
//! mirror reverses the order of tied loci. We use the centre-closest rule of
//! [Cologni & Pibiri, "Canonical Schemes and Minimizers", Proposition 24]:
//! among tied loci take the one closest to the window centre (k-m)/2, and
//! between two equally close (then mirror images i and k-m-i) take the smaller
//! index if x <= rc(x) and the larger otherwise. The mirror fixes the centre,
//! so the rule is equivariant; the distance of a fixed locus to the centre is
//! monotone as the window slides right, so it is also *forward*: the sampled
//! position never decreases, and the number of super-kmers equals the number
//! of minimizer positions exactly (checked at tuple-merge time by the
//! builder).

use crate::hasher::DeterministicHasher;
use crate::kmer::{Kmer, KmerBits};

/// Reverse complement of an m-mer packed in the low `2m` bits of a word.
///
/// This is the u64 body of [`Kmer::reverse_complement`] specialized to a value
/// known to fit one word, regardless of how wide the k-mer type is: complement
/// (A<->T is 00<->10, C<->G is 01<->11, i.e. XOR with 0b10 per base), reverse
/// the 2-bit pairs, then right-align.
#[inline(always)]
pub fn reverse_complement_mmer(mmer: u64, m: usize) -> u64 {
    debug_assert!(m >= 1 && 2 * m < 64, "m must be in [1, 31]");
    let mut x = mmer ^ 0xAAAA_AAAA_AAAA_AAAAu64;
    x = ((x >> 2) & 0x3333_3333_3333_3333u64) | ((x & 0x3333_3333_3333_3333u64) << 2);
    x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0Fu64) | ((x & 0x0F0F_0F0F_0F0F_0F0Fu64) << 4);
    x = x.swap_bytes();
    x >> (64 - 2 * m)
}

/// True if the m-mer equals its own reverse complement (possible only for
/// even m). The streaming same-minimizer memo must fall through to a full
/// verification for such minimizers.
#[inline(always)]
pub fn is_self_rc_mmer(mmer: u64, m: usize) -> bool {
    reverse_complement_mmer(mmer, m) == mmer
}

/// The canonical m-mer at locus `i` of a k-mer, given the k-mer's reverse
/// complement: `kappa(i) = min(x_i, rc(x_i))`.
///
/// Uses `rc(x_i) = rc(x)_{K-m-i}`: the reverse complements of all k-m+1 loci
/// are windows of the single reverse complement of the whole k-mer, so the
/// k-mer is reverse-complemented once instead of once per locus (the caller
/// usually has rc(x) already — a lookup needs it anyway for orientation).
#[inline(always)]
pub fn canonical_mmer_at<const K: usize>(
    kmer: &Kmer<K>,
    kmer_rc: &Kmer<K>,
    m: usize,
    i: usize,
) -> u64
where
    Kmer<K>: KmerBits,
{
    debug_assert!(i + m <= K);
    let mask = (1u64 << (2 * m)) - 1;
    let (fwd, rc) = if K <= 31 {
        let bits = <Kmer<K> as KmerBits>::to_u64(kmer.bits());
        let rc_bits = <Kmer<K> as KmerBits>::to_u64(kmer_rc.bits());
        ((bits >> (2 * i)) & mask, (rc_bits >> (2 * (K - m - i))) & mask)
    } else {
        let bits = <Kmer<K> as KmerBits>::to_u128(kmer.bits());
        let rc_bits = <Kmer<K> as KmerBits>::to_u128(kmer_rc.bits());
        (
            ((bits >> (2 * i)) as u64) & mask,
            ((rc_bits >> (2 * (K - m - i))) as u64) & mask,
        )
    };
    debug_assert_eq!(rc, reverse_complement_mmer(fwd, m));
    fwd.min(rc)
}

/// Information about a minimizer within a k-mer
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinimizerInfo {
    /// The minimizer value: the canonical m-mer at the selected locus
    pub value: u64,
    /// The absolute position of the selected locus in the full string
    /// (not within the k-mer)
    pub position: u64,
    /// The selected locus within the k-mer window, in [0, k-m]
    /// (0 = leftmost m-mer of the window)
    pub pos_in_kmer: usize,
}

impl MinimizerInfo {
    /// Create a new minimizer info
    pub fn new(value: u64, position: u64, pos_in_kmer: usize) -> Self {
        Self {
            value,
            position,
            pos_in_kmer,
        }
    }
}

/// Centre-closest tie resolution (Proposition 24), shared by the batch
/// [`compute_minimizer`] and the incremental [`MinimizerIterator`].
///
/// Scans loci `[lo_locus, hi_locus]` for hashes equal to `min_hash`, takes the
/// tied locus closest to the window centre via the doubled distance
/// `|2i - (K-m)|` (avoiding the half-integer centre of an odd window); between
/// two equally close loci — mirror images — the smaller index wins when
/// `x <= rc(x)`, the larger otherwise.
///
/// Deliberately not inlined: ties fire on ~1e-4 of windows, and letting this
/// loop inline into the lookup path costs measurably more than the call
/// (observed in C++, commit b1d0706).
#[cold]
#[inline(never)]
fn resolve_tie<const K: usize>(
    kmer: &Kmer<K>,
    kmer_rc: &Kmer<K>,
    m: usize,
    hasher: &DeterministicHasher,
    min_hash: u64,
    min_value: u64,
    lo_locus: usize,
    hi_locus: usize,
) -> usize
where
    Kmer<K>: KmerBits,
{
    let two_c = K - m;
    let mut best_dist = u64::MAX;
    let mut lo = lo_locus;
    let mut hi = lo_locus;
    for i in lo_locus..=hi_locus {
        let value = canonical_mmer_at(kmer, kmer_rc, m, i);
        if hasher.hash(value) == min_hash {
            // Our multiply-XOR mixer is a bijection on u64, so equal hashes
            // imply equal kappa values; the scheme would still be correct with
            // a non-bijective hasher (C++'s is), but this documents ours.
            debug_assert_eq!(value, min_value);
            let dist = (2 * i as u64).abs_diff(two_c as u64);
            if dist < best_dist {
                best_dist = dist;
                lo = i;
                hi = i;
            } else if dist == best_dist {
                hi = i;
            }
        }
    }
    // lo and hi are the two equally-closest loci (mirror images), or one locus.
    if lo == hi || !(kmer_rc.bits() < kmer.bits()) {
        lo
    } else {
        hi
    }
}

/// Compute the minimizer of a single k-mer: the locus minimizing
/// `h(kappa(i))`, with the centre-closest tie-break.
///
/// The returned `position` is the locus within the k-mer (identical to
/// `pos_in_kmer`); callers anchoring in a sequence add the window's absolute
/// begin themselves.
pub fn compute_minimizer<const K: usize>(
    kmer: &Kmer<K>,
    kmer_rc: &Kmer<K>,
    m: usize,
    hasher: &DeterministicHasher,
) -> MinimizerInfo
where
    Kmer<K>: KmerBits,
{
    debug_assert!(m >= 1 && m <= K && 2 * m < 64);
    let mask = (1u64 << (2 * m)) - 1;

    // Hot loop, structured like C++'s: the forward m-mer comes from a running
    // window shifted by a constant 2 per locus; only the RC extraction pays a
    // variable shift (`rc(x_i) = rc(x)_{K-m-i}`). The first locus is peeled
    // off so that `min_hash` starts out at a real hash value: initializing it
    // to u64::MAX would make an actual hash of u64::MAX register as a tie
    // rather than as the minimum.
    macro_rules! scan {
        ($bits:expr, $rc_bits:expr) => {{
            let bits = $bits;
            let rc_bits = $rc_bits;
            let mut window = bits;
            let kappa = |window, i: usize| {
                let fwd = (window as u64) & mask;
                let rc = (((rc_bits) >> (2 * (K - m - i))) as u64) & mask;
                debug_assert_eq!(rc, reverse_complement_mmer(fwd, m));
                fwd.min(rc)
            };

            let mut min_value = kappa(window, 0);
            let mut min_hash = hasher.hash(min_value);
            let mut leftmost = 0usize;
            let mut rightmost = 0usize;
            let mut tie = false;

            for i in 1..=(K - m) {
                window >>= 2;
                let value = kappa(window, i);
                let hash = hasher.hash(value);
                if hash < min_hash {
                    min_hash = hash;
                    min_value = value;
                    leftmost = i;
                    rightmost = i;
                    tie = false;
                } else if hash == min_hash {
                    rightmost = i;
                    tie = true;
                }
            }
            (min_value, min_hash, leftmost, rightmost, tie)
        }};
    }

    let (min_value, min_hash, leftmost, rightmost, tie) = if K <= 31 {
        scan!(
            <Kmer<K> as KmerBits>::to_u64(kmer.bits()),
            <Kmer<K> as KmerBits>::to_u64(kmer_rc.bits())
        )
    } else {
        scan!(
            <Kmer<K> as KmerBits>::to_u128(kmer.bits()),
            <Kmer<K> as KmerBits>::to_u128(kmer_rc.bits())
        )
    };

    let chosen = if tie {
        resolve_tie(kmer, kmer_rc, m, hasher, min_hash, min_value, leftmost, rightmost)
    } else {
        leftmost
    };
    MinimizerInfo::new(min_value, chosen as u64, chosen)
}

/// Iterator that computes the minimizer of each k-mer of a sequence, sliding
/// one character at a time ("re-scan" method).
///
/// Computes exactly the same thing as [`compute_minimizer`] incrementally,
/// which a debug assertion in [`Self::next`] checks on every call.
///
/// The only extra state compared to a plain forward minimizer is `num_mins`,
/// the number of loci of the current window attaining the minimum hash: a tie
/// cannot be broken by position alone without breaking mirror-equivariance,
/// so when there is one we have to look at the k-mer's own orientation.
pub struct MinimizerIterator {
    k: usize,
    m: usize,
    window_size: usize, // k - m
    position: u64,
    min_pos_in_kmer: usize,
    min_value: u64,
    min_position: u64,
    min_hash: u64,
    num_mins: u64,
    hasher: DeterministicHasher,
}

impl MinimizerIterator {
    /// Create a new minimizer iterator with a seeded hasher
    ///
    /// # Arguments
    /// * `k` - k-mer size (must be >= m)
    /// * `m` - minimizer size (must be <= k and <= 31)
    /// * `seed` - seed for the hash function (for deterministic results)
    pub fn with_seed(k: usize, m: usize, seed: u64) -> Self {
        assert!(k > 0 && m <= k, "k must be > 0 and m <= k");
        assert!(2 * m < 64, "m must be <= 31 (minimizers are u64)");
        let hasher = DeterministicHasher::new(seed);
        let mut iter = Self {
            k,
            m,
            window_size: k - m,
            position: 0,
            min_pos_in_kmer: 0,
            min_value: 0,
            min_position: 0,
            min_hash: u64::MAX,
            num_mins: 0,
            hasher,
        };
        iter.reset();
        iter
    }

    /// Create a new minimizer iterator (defaults to seed=1)
    pub fn new(k: usize, m: usize, _position: u64) -> Self {
        Self::with_seed(k, m, 1)
    }

    /// Set a new starting position
    pub fn set_position(&mut self, position: u64) {
        self.position = position;
        self.reset();
    }

    /// Reset internal state (called when position changes)
    pub fn reset(&mut self) {
        // min_pos_in_kmer == 0 triggers a rescan on the first next(), with
        // begin = min_position + 1 = position (wrapping, like C++ unsigned).
        self.min_pos_in_kmer = 0;
        self.min_position = self.position.wrapping_sub(1);
        self.num_mins = 0;
    }

    /// Compute the minimizer of the next k-mer of the sequence.
    ///
    /// `kmer_rc` must be the reverse complement of `kmer`; taking it as a
    /// parameter (unlike C++, which recomputes it) lets the streaming query
    /// reuse the RC it already maintains incrementally.
    #[inline]
    pub fn next<const K: usize>(&mut self, kmer: &Kmer<K>, kmer_rc: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        debug_assert_eq!(K, self.k, "k-mer size must match iterator k");
        debug_assert_eq!(kmer.reverse_complement(), *kmer_rc);

        if self.min_pos_in_kmer == 0 {
            // The minimum leaves the window: re-scan to compute the new one.
            self.position = self.min_position.wrapping_add(1);
            self.rescan(kmer, kmer_rc);
        } else {
            self.position += 1;
            let value = canonical_mmer_at(kmer, kmer_rc, self.m, self.window_size);
            let hash = self.hasher.hash(value);
            if hash < self.min_hash {
                self.min_hash = hash;
                self.min_value = value;
                self.min_position = self.position;
                self.min_pos_in_kmer = self.window_size;
                self.num_mins = 1;
            } else {
                // Only the leftmost minimum can leave the window without a
                // re-scan, so the count stays valid across the slide.
                self.min_pos_in_kmer -= 1;
                if hash == self.min_hash {
                    self.num_mins += 1;
                }
            }
        }

        let mut info = MinimizerInfo::new(self.min_value, self.min_position, self.min_pos_in_kmer);
        if self.num_mins > 1 {
            self.break_tie(kmer, kmer_rc, &mut info);
        }

        debug_assert_eq!(
            (info.value, info.pos_in_kmer),
            {
                let reference = compute_minimizer(kmer, kmer_rc, self.m, &self.hasher);
                (reference.value, reference.pos_in_kmer)
            },
            "incremental minimizer disagrees with the batch reference"
        );

        info
    }

    /// Rescan the window to find the minimum over all loci.
    fn rescan<const K: usize>(&mut self, kmer: &Kmer<K>, kmer_rc: &Kmer<K>)
    where
        Kmer<K>: KmerBits,
    {
        let begin = self.position;

        // First locus peeled off the loop: see `compute_minimizer`.
        self.min_value = canonical_mmer_at(kmer, kmer_rc, self.m, 0);
        self.min_hash = self.hasher.hash(self.min_value);
        self.min_pos_in_kmer = 0;
        self.num_mins = 1;

        for i in 1..=self.window_size {
            let value = canonical_mmer_at(kmer, kmer_rc, self.m, i);
            let hash = self.hasher.hash(value);
            if hash < self.min_hash {
                // leftmost anchoring; ties resolved by the caller
                self.min_hash = hash;
                self.min_value = value;
                self.min_pos_in_kmer = i;
                self.num_mins = 1;
            } else if hash == self.min_hash {
                self.num_mins += 1;
            }
        }

        self.position = begin + self.window_size as u64;
        self.min_position = begin + self.min_pos_in_kmer as u64;
    }

    /// Apply the centre-closest tie-break to the *returned* info only: the
    /// iterator state stays leftmost-anchored (tied loci cannot precede the
    /// leftmost minimum, so the rescan trigger is unaffected).
    fn break_tie<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        info: &mut MinimizerInfo,
    ) where
        Kmer<K>: KmerBits,
    {
        debug_assert!(self.num_mins > 1);
        let window_begin = self.min_position - self.min_pos_in_kmer as u64;
        let chosen = resolve_tie(
            kmer,
            kmer_rc,
            self.m,
            &self.hasher,
            self.min_hash,
            self.min_value,
            self.min_pos_in_kmer,
            self.window_size,
        );
        *info = MinimizerInfo::new(self.min_value, window_begin + chosen as u64, chosen);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Kmer;

    #[test]
    fn test_reverse_complement_mmer_matches_kmer_rc() {
        // Compare against Kmer::reverse_complement at matching widths.
        let mmer13 = Kmer::<13>::from_str("ACGTTGCAATCGA").unwrap();
        let rc = mmer13.reverse_complement();
        assert_eq!(
            reverse_complement_mmer(<Kmer<13> as KmerBits>::to_u64(mmer13.bits()), 13),
            <Kmer<13> as KmerBits>::to_u64(rc.bits())
        );
    }

    #[test]
    fn test_self_rc_mmer() {
        // "ACGT" is its own reverse complement (even m only).
        let mmer = Kmer::<5>::from_str("ACGTA").unwrap(); // container only
        let _ = mmer;
        let acgt = {
            let mut bits = 0u64;
            for (i, b) in "ACGT".bytes().enumerate() {
                bits |= ((((b >> 1) & 3) as u64) << (i * 2)) as u64;
            }
            bits
        };
        assert!(is_self_rc_mmer(acgt, 4));
        let acgg = acgt & !(3 << 6) | (((b'G' >> 1) & 3) as u64) << 6;
        assert!(!is_self_rc_mmer(acgg, 4));
    }

    #[test]
    fn test_mirror_equivariance_basic() {
        // The scheme's defining property: x and rc(x) select the same value at
        // mirrored loci.
        let hasher = DeterministicHasher::new(1);
        let seqs = [
            "ACGTACGTACGTACGTACGTACGTACGTACG",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA",
            "ACGTTGCAATCGATCGGCTAGCTAGCTAGCA",
        ];
        for m in [3usize, 4, 7, 12, 13, 16, 19] {
            for s in seqs {
                let kmer = Kmer::<31>::from_str(s).unwrap();
                let rc = kmer.reverse_complement();
                let a = compute_minimizer(&kmer, &rc, m, &hasher);
                let b = compute_minimizer(&rc, &kmer, m, &hasher);
                assert_eq!(a.value, b.value, "value differs, m={m} s={s}");
                assert_eq!(
                    a.pos_in_kmer + b.pos_in_kmer,
                    31 - m,
                    "loci are not mirror images, m={m} s={s}"
                );
                assert_eq!(a.value, canonical_mmer_at(&kmer, &rc, m, a.pos_in_kmer));
            }
        }
    }

    #[test]
    fn test_minimizer_iterator_basic() {
        let kmer: Kmer<15> = Kmer::from_str("ACGTACGTACGTACG").unwrap();
        let rc = kmer.reverse_complement();

        let mut iter = MinimizerIterator::new(15, 7, 0);
        let mini = iter.next(&kmer, &rc);

        assert!(mini.pos_in_kmer <= 15 - 7);
        assert_eq!(mini.position, mini.pos_in_kmer as u64);
    }

    #[test]
    fn test_minimizer_consistency_fresh_vs_sequential() {
        // A fresh iterator gives the same minimizer as processing k-mers
        // sequentially through a string (and both match the batch reference,
        // via the debug assertion inside next()).
        let seq = "ATTTTCAGGATGTTTTCAGGTTCATCATCTCCCTTCTTTGCAGGATAGTAGATAAGATCGCTCATCAACGGATGTTGTGT";
        let k = 31usize;
        let m = 19usize;
        let seed = 1u64;
        let num_kmers = seq.len() - k + 1;

        let mut seq_iter = MinimizerIterator::with_seed(k, m, seed);
        seq_iter.set_position(0);

        for kmer_pos in 0..num_kmers {
            let kmer_str = &seq[kmer_pos..kmer_pos + k];
            let kmer = Kmer::<31>::from_str(kmer_str).unwrap();
            let rc = kmer.reverse_complement();

            let seq_mini = seq_iter.next(&kmer, &rc);

            let mut fresh_iter = MinimizerIterator::with_seed(k, m, seed);
            let fresh_mini = fresh_iter.next(&kmer, &rc);

            assert_eq!(seq_mini.value, fresh_mini.value, "value mismatch at {kmer_pos}");
            assert_eq!(
                seq_mini.pos_in_kmer, fresh_mini.pos_in_kmer,
                "pos mismatch at {kmer_pos}"
            );
            // Forwardness: absolute sampled position never decreases and
            // matches window anchoring.
            assert_eq!(seq_mini.position, kmer_pos as u64 + seq_mini.pos_in_kmer as u64);
        }
    }

    #[test]
    fn test_forwardness_on_random_sequence() {
        // The sampled absolute position must never decrease as the window
        // slides (Proposition 24) — this is what makes one-tuple-per-position
        // sound in the builder. Small m makes ties frequent.
        let mut state = 0xDEADBEEFu64;
        let mut next_base = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 33) as usize % 4]
        };
        let seq: Vec<u8> = (0..4000).map(|_| next_base()).collect();
        let k = 21usize;
        for m in [3usize, 4, 5] {
            let mut iter = MinimizerIterator::with_seed(k, m, 1);
            iter.set_position(0);
            let mut prev_pos = 0u64;
            for w in seq.windows(k) {
                let kmer = Kmer::<21>::from_ascii_unchecked(w);
                let rc = kmer.reverse_complement();
                let mini = iter.next(&kmer, &rc);
                assert!(
                    mini.position >= prev_pos,
                    "sampled position decreased: {} -> {} (m={m})",
                    prev_pos,
                    mini.position
                );
                prev_pos = mini.position;
            }
        }
    }
}
