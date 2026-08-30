//! Property tests for the unified canonical minimizer scheme.
//!
//! Port of C++ sshash `test/test_minimizer.cpp`: an independent brute-force
//! reference (re-derived from the spec, sharing no code with the production
//! implementation beyond the hasher) is checked against `compute_minimizer`
//! and the incremental `MinimizerIterator` for:
//!
//! - mirror-equivariance: x and rc(x) select the same m-mer at mirrored loci;
//! - incremental == batch == brute-force on every window of a sequence;
//! - forwardness: the sampled absolute position never decreases
//!   (Proposition 24 — what makes one-tuple-per-position sound);
//! - the tie-break is actually exercised (small m makes ties frequent), so a
//!   hashing change cannot silently make these tests vacuous.
//!
//! K is odd-only in [3, 63] (crate constraint); m runs over every value in
//! [1, min(K-1, 31)] including even m (salmon ships m=16).

use sshash_lib::hasher::DeterministicHasher;
use sshash_lib::kmer::{Kmer, KmerBits};
use sshash_lib::minimizer::{MinimizerIterator, canonical_mmer_at, compute_minimizer};

// ---------------------------------------------------------------------------
// Independent reference implementation (from the spec)
// ---------------------------------------------------------------------------

fn encode(b: u8) -> u64 {
    ((b >> 1) & 3) as u64
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => unreachable!("non-ACGT base"),
    }
}

/// Packed numeric value of an ASCII sequence (base i at bits 2i), as u128 so
/// it covers K up to 63.
fn ref_pack(bytes: &[u8]) -> u128 {
    bytes
        .iter()
        .enumerate()
        .fold(0u128, |acc, (i, &b)| acc | ((encode(b) as u128) << (2 * i)))
}

fn ref_revcomp(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().rev().map(|&b| comp(b)).collect()
}

/// kappa: the canonical m-mer value at one locus.
fn ref_canonical_mmer(mmer: &[u8]) -> u64 {
    let fwd = ref_pack(mmer) as u64;
    let rc = ref_pack(&ref_revcomp(mmer)) as u64;
    fwd.min(rc)
}

/// The spec, written directly: argmin over loci of h(kappa(i)); ties broken
/// centre-closest, and between two equally close (mirror-image) loci the
/// smaller index iff x <= rc(x) under the packed numeric order.
fn ref_compute_minimizer(window: &[u8], m: usize, hasher: &DeterministicHasher) -> (u64, usize, bool) {
    let k = window.len();
    let vals: Vec<u64> = (0..=k - m).map(|i| ref_canonical_mmer(&window[i..i + m])).collect();
    let hashes: Vec<u64> = vals.iter().map(|&v| hasher.hash(v)).collect();
    let min_hash = *hashes.iter().min().unwrap();
    let tied: Vec<usize> = (0..hashes.len()).filter(|&i| hashes[i] == min_hash).collect();
    let had_tie = tied.len() > 1;

    let two_c = (k - m) as i64;
    let dist = |i: usize| (2 * i as i64 - two_c).abs();
    let best = tied.iter().map(|&i| dist(i)).min().unwrap();
    let closest: Vec<usize> = tied.into_iter().filter(|&i| dist(i) == best).collect();

    let chosen = if closest.len() == 1 {
        closest[0]
    } else {
        assert_eq!(closest.len(), 2, "at most two loci can be equally centre-close");
        assert_eq!(closest[0] + closest[1], k - m, "equally close loci are mirror images");
        let fwd = ref_pack(window);
        let rc = ref_pack(&ref_revcomp(window));
        if fwd <= rc { closest[0] } else { closest[1] }
    };
    (vals[chosen], chosen, had_tie)
}

// ---------------------------------------------------------------------------
// Deterministic sequence generation
// ---------------------------------------------------------------------------

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() % 4) as usize]
    }
    fn seq(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.base()).collect()
    }
}

// ---------------------------------------------------------------------------
// Generic property checks
// ---------------------------------------------------------------------------

fn check_k<const K: usize>(num_kmers: usize, seq_len: usize)
where
    Kmer<K>: KmerBits,
{
    let hasher = DeterministicHasher::new(1);
    let max_m = (K - 1).min(31);
    let mut total_ties_small_m = 0u64;

    for m in 1..=max_m {
        let mut rng = Lcg(0x9E37_79B9 ^ (K as u64) << 8 ^ m as u64);

        // --- Mirror-equivariance + batch == reference, on random k-mers ---
        for _ in 0..num_kmers {
            let bytes = rng.seq(K);
            let kmer = Kmer::<K>::from_ascii_unchecked(&bytes);
            let rc = kmer.reverse_complement();

            let a = compute_minimizer(&kmer, &rc, m, &hasher);
            let b = compute_minimizer(&rc, &kmer, m, &hasher);
            assert_eq!(a.value, b.value, "mirror value K={K} m={m}");
            assert_eq!(a.pos_in_kmer + b.pos_in_kmer, K - m, "mirror locus K={K} m={m}");
            assert_eq!(
                a.value,
                canonical_mmer_at(&kmer, &rc, m, a.pos_in_kmer),
                "value must be kappa at the reported locus K={K} m={m}"
            );

            let (ref_val, ref_pos, had_tie) = ref_compute_minimizer(&bytes, m, &hasher);
            assert_eq!(a.value, ref_val, "batch != reference (value) K={K} m={m}");
            assert_eq!(a.pos_in_kmer, ref_pos, "batch != reference (locus) K={K} m={m}");
            if m <= 4 && had_tie {
                total_ties_small_m += 1;
            }
        }

        // --- Incremental == reference + forwardness, along a sequence ---
        let seq = rng.seq(seq_len);
        let mut iter = MinimizerIterator::with_seed(K, m, 1);
        iter.set_position(0);
        let mut prev_position = 0u64;
        for (start, w) in seq.windows(K).enumerate() {
            let kmer = Kmer::<K>::from_ascii_unchecked(w);
            let rc = kmer.reverse_complement();
            let mini = iter.next(&kmer, &rc);

            let (ref_val, ref_pos, _) = ref_compute_minimizer(w, m, &hasher);
            assert_eq!(mini.value, ref_val, "iterator != reference (value) K={K} m={m} start={start}");
            assert_eq!(mini.pos_in_kmer, ref_pos, "iterator != reference (locus) K={K} m={m} start={start}");
            assert_eq!(
                mini.position,
                start as u64 + mini.pos_in_kmer as u64,
                "absolute position anchoring K={K} m={m} start={start}"
            );
            assert!(
                mini.position >= prev_position,
                "forwardness violated: {} -> {} (K={K} m={m} start={start})",
                prev_position,
                mini.position
            );
            prev_position = mini.position;
        }
    }

    // The tie-break must actually fire for small m, or the tests above prove
    // nothing about it.
    assert!(
        total_ties_small_m > 0,
        "no ties encountered at m <= 4 for K={K}: the tie-break is untested"
    );
}

// One test per K so they run in parallel; the reduced grid keeps the default
// `cargo test` fast, and `--ignored` runs the heavy counts.

#[test]
fn properties_k5() {
    check_k::<5>(2_000, 2_000);
}

#[test]
fn properties_k15() {
    check_k::<15>(1_000, 2_000);
}

#[test]
fn properties_k21() {
    check_k::<21>(1_000, 2_000);
}

#[test]
fn properties_k31() {
    check_k::<31>(1_000, 2_000);
}

#[test]
fn properties_k33() {
    check_k::<33>(500, 1_000);
}

#[test]
fn properties_k47() {
    check_k::<47>(500, 1_000);
}

#[test]
fn properties_k63() {
    check_k::<63>(500, 1_000);
}

#[test]
#[ignore = "heavy grid; run with cargo test --release -- --ignored"]
fn properties_heavy() {
    check_k::<5>(20_000, 20_000);
    check_k::<15>(20_000, 20_000);
    check_k::<21>(20_000, 20_000);
    check_k::<31>(20_000, 20_000);
    check_k::<33>(10_000, 10_000);
    check_k::<47>(10_000, 10_000);
    check_k::<63>(10_000, 10_000);
}
