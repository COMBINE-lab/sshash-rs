// SSHash: Sparse and Skew Hashing of k-mers
//
// A Rust implementation of the SSHash compressed k-mer dictionary,
// providing efficient storage and querying of k-mer sets.
//
// Note: This branch contains the *initial* Claude-driven rewrite, which adopted
// a slightly different structure for encoding light/mid/heavy buckets than the C++
// version. The main branch now contains the C++-equivalent structure, but we maintain
// this branch for posterity and potential comparison.

#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

pub mod builder;
pub mod constants;
pub mod dictionary;
pub mod encoding;
pub mod hasher;
pub mod kmer;
pub mod minimizer;
pub mod minimizers_control_map;
pub mod mphf_config;
pub mod offsets;
pub mod serialization;
pub mod sparse_and_skew_index;
pub mod spectrum_preserving_string_set;
pub mod streaming_query;

// Re-export common types at crate root
pub use builder::{BuildConfiguration, CfSegData, DictionaryBuilder, parse_cf_seg};
pub use dictionary::Dictionary;
pub use kmer::{Kmer, Kmer21, Kmer31, Kmer63, KmerBits};
pub use minimizer::{MinimizerInfo, MinimizerIterator};
pub use minimizers_control_map::{BucketType, MinimizersControlMap, MinimizersControlMapBuilder};
pub use streaming_query::{LookupResult, StreamingQuery};

/// Dispatch to the correct const generic `K` based on a runtime `k` value.
///
/// The const generic `K` determines both the storage type (`u64` for K ≤ 31,
/// `u128` for K > 31) **and** the k-mer length used in operations like
/// `reverse_complement()`, `from_str()`, and `decode_kmer()`.  Therefore `K`
/// must always equal the actual k value, not merely the maximum for the
/// storage class.
///
/// `k` must be an odd value in \[3, 63\]; these are the only values for which
/// [`KmerBits`] is implemented.
///
/// # Usage
///
/// ```ignore
/// use sshash_lib::dispatch_on_k;
///
/// dispatch_on_k!(k, K => {
///     let kmer = Kmer::<K>::from_string(s)?;
///     dict.lookup::<K>(&kmer)
/// })
/// ```
///
/// # Panics
///
/// Panics at runtime if `k` is even or outside the \[3, 63\] range.
#[macro_export]
macro_rules! dispatch_on_k {
    ($k:expr, $K:ident => $body:expr) => {{
        match $k {
            3 => {
                const $K: usize = 3;
                $body
            }
            5 => {
                const $K: usize = 5;
                $body
            }
            7 => {
                const $K: usize = 7;
                $body
            }
            9 => {
                const $K: usize = 9;
                $body
            }
            11 => {
                const $K: usize = 11;
                $body
            }
            13 => {
                const $K: usize = 13;
                $body
            }
            15 => {
                const $K: usize = 15;
                $body
            }
            17 => {
                const $K: usize = 17;
                $body
            }
            19 => {
                const $K: usize = 19;
                $body
            }
            21 => {
                const $K: usize = 21;
                $body
            }
            23 => {
                const $K: usize = 23;
                $body
            }
            25 => {
                const $K: usize = 25;
                $body
            }
            27 => {
                const $K: usize = 27;
                $body
            }
            29 => {
                const $K: usize = 29;
                $body
            }
            31 => {
                const $K: usize = 31;
                $body
            }
            33 => {
                const $K: usize = 33;
                $body
            }
            35 => {
                const $K: usize = 35;
                $body
            }
            37 => {
                const $K: usize = 37;
                $body
            }
            39 => {
                const $K: usize = 39;
                $body
            }
            41 => {
                const $K: usize = 41;
                $body
            }
            43 => {
                const $K: usize = 43;
                $body
            }
            45 => {
                const $K: usize = 45;
                $body
            }
            47 => {
                const $K: usize = 47;
                $body
            }
            49 => {
                const $K: usize = 49;
                $body
            }
            51 => {
                const $K: usize = 51;
                $body
            }
            53 => {
                const $K: usize = 53;
                $body
            }
            55 => {
                const $K: usize = 55;
                $body
            }
            57 => {
                const $K: usize = 57;
                $body
            }
            59 => {
                const $K: usize = 59;
                $body
            }
            61 => {
                const $K: usize = 61;
                $body
            }
            63 => {
                const $K: usize = 63;
                $body
            }
            other => panic!(
                "Unsupported k value: {}. k must be an odd value between 3 and 63.",
                other,
            ),
        }
    }};
}

/// Version information
pub fn version() -> (u8, u8, u8) {
    constants::VERSION
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_version() {
        let (major, minor, patch) = version();
        assert_eq!(major, 0);
        assert_eq!(minor, 1);
        assert_eq!(patch, 0);
    }
}
