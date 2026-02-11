//! Integration tests for the build pipeline
//!
//! These tests exercise the full build pipeline from parsing to index construction.

use sshash_lib::builder::{
    config::BuildConfiguration,
    encode::Encoder,
    minimizer_tuples::compute_minimizer_tuples,
    buckets::{classify_into_buckets, BucketStatistics, BucketType},
};
use sshash_lib::sparse_and_skew_index::SparseAndSkewIndex;
use value_traits::slices::SliceByValue;

#[test]
fn test_end_to_end_pipeline_simple() {
    // Step 1: Setup configuration
    let config = BuildConfiguration::new(31, 13).unwrap();
    
    // Step 2: Parse and encode sequences
    let mut encoder = Encoder::<31>::new();
    let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 40 bases -> 10 k-mers
    encoder.add_sequence(sequence).unwrap();
    
    let spss = encoder.build(config.m);
    
    // Verify SPSS
    assert_eq!(spss.k(), 31);
    assert_eq!(spss.num_strings(), 1);
    
    // Step 3: Compute minimizer tuples
    let tuples = compute_minimizer_tuples::<31>(&spss, &config);
    
    // Should have tuples for each k-mer
    assert!(!tuples.is_empty());
    assert!(tuples.len() <= 10); // May be coalesced
    
    // Step 4: Classify into buckets
    let buckets = classify_into_buckets(tuples);
    
    assert!(!buckets.is_empty());
    
    // Step 5: Compute statistics
    let mut stats = BucketStatistics::new();
    for bucket in &buckets {
        stats.add_bucket(&bucket.tuples);
    }
    
    assert_eq!(stats.num_buckets, buckets.len() as u64);
    
    // Step 6: Build sparse and skew index
    let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
    
    assert_eq!(index.control_codewords.len(), stats.num_buckets as usize);
}

#[test]
fn test_end_to_end_multiple_sequences() {
    // Step 1: Setup configuration
    let config = BuildConfiguration::new(31, 13).unwrap();
    
    // Step 2: Parse and encode multiple sequences
    let mut encoder = Encoder::<31>::new();
    
    // Add several sequences
    let sequences = vec![
        b"ACGTACGTACGTACGTACGTACGTACGTACG".as_slice(), // Exactly 31 bases
        b"TGCATGCATGCATGCATGCATGCATGCATGC".as_slice(),
        b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".as_slice(),
        b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".as_slice(),
    ];
    
    for seq in &sequences {
        encoder.add_sequence(seq).unwrap();
    }
    
    let spss = encoder.build(config.m);
    assert_eq!(spss.num_strings(), 4);
    
    // Step 3: Compute minimizer tuples
    let tuples = compute_minimizer_tuples::<31>(&spss, &config);
    
    // Step 4: Classify and build index
    let buckets = classify_into_buckets(tuples);
    let mut stats = BucketStatistics::new();
    for bucket in &buckets {
        stats.add_bucket(&bucket.tuples);
    }
    
    // Should have at least 1 bucket
    assert!(stats.num_buckets > 0);
    
    // Build index
    let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
    assert_eq!(index.control_codewords.len(), stats.num_buckets as usize);
}

#[test]
fn test_bucket_type_distribution() {
    // Create a configuration
    let config = BuildConfiguration::new(31, 13).unwrap();
    
    // Create sequences that will generate different bucket types
    let mut encoder = Encoder::<31>::new();
    
    // Sequence with repetitive k-mers (will create larger buckets)
    let long_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    encoder.add_sequence(long_seq).unwrap();
    
    let spss = encoder.build(config.m);
    let tuples = compute_minimizer_tuples::<31>(&spss, &config);
    let buckets = classify_into_buckets(tuples);
    
    // Count bucket types
    let mut singleton_count = 0;
    let mut light_count = 0;
    let mut heavy_count = 0;
    
    for bucket in &buckets {
        match bucket.bucket_type() {
            BucketType::Singleton => singleton_count += 1,
            BucketType::Light => light_count += 1,
            BucketType::Heavy => heavy_count += 1,
        }
    }
    
    // Should have at least some buckets
    let total = singleton_count + light_count + heavy_count;
    assert!(total > 0);
    assert_eq!(total, buckets.len());
}

#[test]
fn test_canonical_mode_pipeline() {
    // Test with canonical mode enabled
    let mut config = BuildConfiguration::new(31, 13).unwrap();
    config.canonical = true;
    
    let mut encoder = Encoder::<31>::new();
    let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACG";
    encoder.add_sequence(sequence).unwrap();
    
    let spss = encoder.build(config.m);
    let tuples = compute_minimizer_tuples::<31>(&spss, &config);
    
    // In canonical mode, should still get valid tuples
    assert!(!tuples.is_empty());
    
    // All pos_in_kmer values should be valid (0 to k-m)
    for tuple in &tuples {
        assert!(tuple.pos_in_kmer <= (31 - 13) as u8);
    }
    
    // Build index in canonical mode
    let buckets = classify_into_buckets(tuples);
    let index = SparseAndSkewIndex::build::<31>(buckets, 32, &spss, false);
    
    assert!(index.control_codewords.len() > 0);
}
