//! MPHF (Minimal Perfect Hash Function) type configuration
//!
//! Central module for MPHF type aliases and helpers used across the crate.
//! We use PHast (Perfect Hashing made fast) with rapidhash instead of the default
//! SipHash hasher for faster hash evaluations during both construction and query.
//!
//! rapidhash is CPU-feature independent (unlike ahash which switches algorithm
//! when AES-NI is available via target-cpu=native, breaking serialized indices).
//!
//! PHast offers very fast evaluation and size below 2 bits/key.

use ph::seedable_hash::BuildRapidHash;
use ph::phast;
use ph::seeds::Bits8;
use std::hash::Hash;
use std::io;

/// The seeded hasher used inside our MPHF functions.
///
/// Uses rapidhash via the `ph` crate's `BuildRapidHash` adapter, which is
/// CPU-feature independent (no AES-NI divergence). This ensures indices built
/// with or without `target-cpu=native` produce the same MPHF evaluations.
pub type MphfHasher = BuildRapidHash;

/// Our MPHF type — PHast with rapidhash instead of default SipHash.
///
/// Type parameters:
/// - `Bits8`: 8 bits per seed (default, recommended)
/// - `phast::SeedOnly`: Regular PHast variant (not PHast+)
/// - `phast::DefaultCompressedArray`: Default compressed array implementation
/// - `MphfHasher`: rapidhash-based seeded hasher (faster than default SipHash)
pub type Mphf = phast::Function<Bits8, phast::SeedOnly, phast::DefaultCompressedArray, MphfHasher>;

/// Create the deterministic MPHF hasher.
///
/// Must use the same hasher at both build and load time to ensure
/// the serialized MPHF produces correct results after deserialization.
pub fn mphf_hasher() -> MphfHasher {
    BuildRapidHash
}

/// Create PHast parameters with default settings (Bits8, optimal bucket size).
pub fn mphf_params() -> phast::Params<Bits8> {
    phast::Params::new(Bits8, phast::bits_per_seed_to_100_bucket_size(8))
}

/// Build an MPHF from a slice of keys (single-threaded).
///
/// Generic over the key type — works with `u64` (minimizers, k ≤ 31 k-mers)
/// or `u128` (k > 31 k-mers). The MPHF itself is key-type-erased after
/// construction; only the build and query calls need to agree on the type.
pub fn build_mphf_from_slice<T: Hash + Clone>(keys: &[T]) -> Mphf {
    Mphf::with_slice_p_hash_sc(keys, &mphf_params(), mphf_hasher(), phast::SeedOnly)
}

/// Build an MPHF from an owned Vec of keys (single-threaded, avoids cloning).
///
/// Generic over the key type — see [`build_mphf_from_slice`] for details.
pub fn build_mphf_from_vec<T: Hash>(keys: Vec<T>) -> Mphf {
    Mphf::with_vec_p_hash_sc(keys, &mphf_params(), mphf_hasher(), phast::SeedOnly)
}

/// Build an MPHF from a slice of keys using multiple threads.
///
/// Falls back to single-threaded construction when `threads == 1`.
/// Uses rayon internally (PHast's `_mt` variant).
pub fn build_mphf_from_slice_mt<T: Hash + Sync + Send + Clone>(keys: &[T], threads: usize) -> Mphf {
    Mphf::with_slice_p_threads_hash_sc(keys, &mphf_params(), threads, mphf_hasher(), phast::SeedOnly)
}

/// Read (deserialize) an MPHF from a reader.
///
/// Uses the same deterministic rapidhash hasher and SeedOnly chooser
/// as used during construction, ensuring correct round-trip behavior.
pub fn read_mphf(reader: &mut dyn io::Read) -> io::Result<Mphf> {
    Mphf::read_with_hasher_sc(reader, mphf_hasher(), phast::SeedOnly)
}
