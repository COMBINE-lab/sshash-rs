// SSHash: Sparse and Skew Hashing of k-mers
//
// A Rust implementation of the SSHash compressed k-mer dictionary,
// providing efficient storage and querying of k-mer sets.

#![doc = include_str!("../README.md")]
#![warn(missing_docs)]

pub mod constants;
pub mod encoding;
pub mod hasher;
pub mod kmer;
pub mod minimizer;
pub mod mphf_config;
pub mod offsets;
pub mod minimizers_control_map;
pub mod spectrum_preserving_string_set;
pub mod sparse_and_skew_index;
pub mod streaming_query;
pub mod dictionary;
pub mod builder;
pub mod serialization;

// Re-export common types at crate root 
pub use kmer::{Kmer, Kmer21, Kmer31, Kmer63, KmerBits};
pub use minimizer::{MinimizerInfo, MinimizerIterator};
pub use minimizers_control_map::{MinimizersControlMap, MinimizersControlMapBuilder, BucketType};
pub use streaming_query::{LookupResult, StreamingQuery};
pub use dictionary::Dictionary;
pub use builder::{BuildConfiguration, DictionaryBuilder};

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
