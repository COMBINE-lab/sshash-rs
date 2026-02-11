//! Constants and configuration for SSHash
//!
//! This module defines compile-time and runtime constants used throughout
//! the library, including valid k-mer sizes and algorithm parameters.

/// Invalid value sentinel (matching C++ behavior)
pub const INVALID_UINT64: u64 = u64::MAX;

/// Default seed for hash functions
pub const DEFAULT_SEED: u64 = 1;

/// Default RAM limit in GiB for construction
pub const DEFAULT_RAM_LIMIT_GIB: usize = 8;

/// For PTHash/PHast MPHF construction
pub const DEFAULT_LAMBDA: f64 = 5.0;
/// Average partition size for MPHF construction
pub const AVG_PARTITION_SIZE: usize = 3_000_000;

/// Bucket size range parameters (for sparse and skew index)
pub const MIN_L: usize = 6;
/// Maximum bucket size parameter
pub const MAX_L: usize = 13;

/// Orientation constants
pub const FORWARD_ORIENTATION: i8 = 1;
/// Backward (reverse complement) orientation
pub const BACKWARD_ORIENTATION: i8 = -1;

/// Version number
pub const VERSION: (u8, u8, u8) = (0, 1, 0);

/// All valid k-mer sizes (odd numbers from 3 to 63)
/// This is the single source of truth for supported K values
pub const VALID_K_VALUES: &[usize] = &[
    3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49,
    51, 53, 55, 57, 59, 61, 63,
];

/// Check if a k-mer size is valid
#[inline]
pub const fn is_valid_k(k: usize) -> bool {
    k >= 3 && k <= 63 && k % 2 == 1
}

/// Maximum k-mer size supported
pub const MAX_K: usize = 63;

/// Minimum k-mer size supported
pub const MIN_K: usize = 3;

/// Compute ceil(log2(x)), matching C++ std::ceil(std::log2(x)).
///
/// Returns 0 for x <= 1, and the minimum number of bits needed to
/// represent values in [0, x) for x >= 2.
#[inline]
pub const fn ceil_log2(x: u64) -> usize {
    if x <= 1 {
        0
    } else {
        64 - (x - 1).leading_zeros() as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_k_values() {
        // Verify all values in VALID_K_VALUES are actually valid
        for &k in VALID_K_VALUES {
            assert!(is_valid_k(k), "k={} should be valid", k);
        }

        // Verify count (31 odd numbers from 3 to 63)
        assert_eq!(VALID_K_VALUES.len(), 31);

        // Verify no duplicates
        let mut sorted = VALID_K_VALUES.to_vec();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), VALID_K_VALUES.len());
    }

    #[test]
    fn test_is_valid_k() {
        // Valid cases
        assert!(is_valid_k(3));
        assert!(is_valid_k(31));
        assert!(is_valid_k(63));
        assert!(is_valid_k(21));

        // Invalid cases (even)
        assert!(!is_valid_k(2));
        assert!(!is_valid_k(4));
        assert!(!is_valid_k(30));
        assert!(!is_valid_k(32));

        // Invalid cases (out of range)
        assert!(!is_valid_k(1));
        assert!(!is_valid_k(65));
        assert!(!is_valid_k(100));
    }
}
