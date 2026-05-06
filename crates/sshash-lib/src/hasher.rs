//! Deterministic hasher for minimizers.
//!
//! Uses a fast multiply-XOR mixer (matching the C++ sshash mixer_64 approach)
//! for the hot-path minimizer computation. The magic constant is derived from
//! the seed via rapidhash for determinism and portability.

use rapidhash::fast::RapidHasher;
use std::hash::Hasher;

const MIXER_MULTIPLIER: u64 = 0x517cc1b727220a95;

/// A deterministic hasher with a seeded state.
///
/// Uses a fast multiply-XOR mixer for the minimizer hot path.
/// The magic constant is derived from the seed via rapidhash at construction
/// time, ensuring CPU-feature independence and portability.
#[derive(Clone)]
pub struct DeterministicHasher {
    seed: u64,
    magic: u64,
}

impl DeterministicHasher {
    /// Create a new deterministic hasher with the given seed
    pub fn new(seed: u64) -> Self {
        let mut h = RapidHasher::new(0);
        h.write_u64(seed);
        let magic = h.finish();
        Self { seed, magic }
    }

    /// Hash a u64 value using the fast multiply-XOR mixer.
    /// This matches C++ sshash's mixer_64: `(x * MULTIPLIER) ^ magic`
    #[inline(always)]
    pub fn hash_u64(&self, value: u64) -> u64 {
        value.wrapping_mul(MIXER_MULTIPLIER) ^ self.magic
    }

    /// Hash a u64 using the seeded hasher (for compatibility with iterator API)
    #[inline(always)]
    pub fn hash(&self, value: u64) -> u64 {
        self.hash_u64(value)
    }

    /// Get the seed value
    pub fn seed(&self) -> u64 {
        self.seed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_deterministic_hashing() {
        let hasher1 = DeterministicHasher::new(42);
        let hasher2 = DeterministicHasher::new(42);
        let hasher3 = DeterministicHasher::new(43);

        let value = 0x123456789abcdef0u64;

        // Same seed should produce same hash
        assert_eq!(hasher1.hash(value), hasher2.hash(value));

        // Different seed should produce different hash
        assert_ne!(hasher1.hash(value), hasher3.hash(value));
    }

    #[test]
    fn test_different_values_produce_different_hashes() {
        let hasher = DeterministicHasher::new(1);

        let hash1 = hasher.hash(100);
        let hash2 = hasher.hash(101);

        // Different inputs should produce different outputs
        assert_ne!(hash1, hash2);
    }
}
