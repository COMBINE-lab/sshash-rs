//! Builder module for constructing SSHash dictionaries
//!
//! This module implements the multi-step build pipeline:
//! 1. Parse and encode input strings (FASTA/FASTQ)
//! 2. Compute minimizer tuples for each k-mer
//! 3. Merge and sort minimizer tuples
//! 4. Build minimizers control map (MPHF)
//! 5. Hash minimizers with MPHF IDs
//! 6. Build sparse and skew index
//! 7. Finalize dictionary structure

pub mod config;
pub mod parse;
pub mod encode;
pub mod minimizer_tuples;
pub mod buckets;
pub mod dictionary_builder;
pub mod external_sort;

pub use config::BuildConfiguration;
pub use minimizer_tuples::MinimizerTuple;
pub use dictionary_builder::DictionaryBuilder;

// Re-export for convenience
pub use crate::kmer::Kmer;
