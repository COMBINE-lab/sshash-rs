//! Deterministic hasher for minimizers using ahash.
//!
//! This uses AHasher with explicit seeds to provide deterministic hashing.
//! We can swap the hash implementation later without changing callers.

use ahash::RandomState;
use std::hash::{BuildHasher, Hasher};

/// A deterministic hasher with a seeded state
#[derive(Clone)]
pub struct DeterministicHasher {
    seed: u64,
    state: RandomState,
}

impl DeterministicHasher {
    /// Create a new deterministic hasher with the given seed
    pub fn new(seed: u64) -> Self {
        let state = RandomState::with_seeds(seed, !seed, seed, !seed);
        Self { seed, state }
    }

    /// Hash a u64 value using a seeded AHasher
    #[inline]
    pub fn hash_u64(&self, value: u64) -> u64 {
        let mut hasher = self.state.build_hasher();
        hasher.write_u64(value);
        hasher.finish()
    }

    /// Hash a u64 using the seeded hasher (for compatibility with iterator API)
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
