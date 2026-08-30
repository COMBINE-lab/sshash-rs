//! Bucket classification and statistics
//!
//! This module classifies minimizer buckets into singleton, light (sparse),
//! and heavy (skew) categories based on bucket size.

use crate::builder::minimizer_tuples::MinimizerTuple;
use crate::constants::MIN_L;
use tracing::info;

/// Bucket size threshold between light and heavy buckets
pub const MIN_BUCKET_SIZE: usize = 1 << MIN_L; // 2^6 = 64

/// Statistics about bucket sizes and distribution
#[derive(Debug, Clone)]
pub struct BucketStatistics {
    /// Total number of buckets (unique minimizers)
    pub num_buckets: u64,
    
    /// Total number of k-mers
    pub num_kmers: u64,
    
    /// Total number of minimizer positions
    pub num_positions: u64,
    
    /// Number of singleton buckets (size == 1)
    pub num_singleton_buckets: u64,
    
    /// Number of light buckets (1 < size <= MIN_BUCKET_SIZE)
    pub num_light_buckets: u64,
    
    /// Number of heavy buckets (size > MIN_BUCKET_SIZE)
    pub num_heavy_buckets: u64,
    
    /// Maximum observed bucket size
    pub max_bucket_size: usize,
    
    /// Total positions in light buckets
    pub num_positions_in_light: u64,
    
    /// Total positions in heavy buckets
    pub num_positions_in_heavy: u64,
    
    /// Total super-k-mers in buckets larger than 1
    pub num_super_kmers_in_non_singleton: u64,
}

impl BucketStatistics {
    /// Create a new statistics tracker
    pub fn new() -> Self {
        Self {
            num_buckets: 0,
            num_kmers: 0,
            num_positions: 0,
            num_singleton_buckets: 0,
            num_light_buckets: 0,
            num_heavy_buckets: 0,
            max_bucket_size: 0,
            num_positions_in_light: 0,
            num_positions_in_heavy: 0,
            num_super_kmers_in_non_singleton: 0,
        }
    }
    
    /// Record statistics for a bucket
    pub fn add_bucket(&mut self, bucket: &[MinimizerTuple]) {
        self.num_buckets += 1;
        // Forward scheme: one tuple per minimizer position (enforced during
        // classification/merge), so the size is just the tuple count.
        debug_assert!(bucket.windows(2).all(|p| p[0].pos_in_seq != p[1].pos_in_seq));
        let bucket_size: usize = bucket.len();

        if bucket_size > self.max_bucket_size {
            self.max_bucket_size = bucket_size;
        }
        
        // Count total k-mers in this bucket
        let kmers_in_bucket: u64 = bucket.iter()
            .map(|mt| mt.num_kmers_in_super_kmer as u64)
            .sum();
        self.num_kmers += kmers_in_bucket;
        
        // Classify bucket
        match bucket_size {
            1 => {
                self.num_singleton_buckets += 1;
                self.num_positions += 1;
            }
            2..=MIN_BUCKET_SIZE => {
                self.num_light_buckets += 1;
                self.num_positions_in_light += bucket_size as u64;
                self.num_positions += bucket_size as u64;
                self.num_super_kmers_in_non_singleton += bucket.len() as u64;
            }
            _ => {
                self.num_heavy_buckets += 1;
                self.num_positions_in_heavy += bucket_size as u64;
                self.num_positions += bucket_size as u64;
                self.num_super_kmers_in_non_singleton += bucket.len() as u64;
            }
        }
    }
    
    /// Log statistics summary via tracing
    pub fn print_summary(&self) {
        info!("Bucket Statistics:");
        info!("  Total buckets: {}", self.num_buckets);
        info!("  Total k-mers: {}", self.num_kmers);
        info!("  Total positions: {}", self.num_positions);
        
        info!("  Singleton buckets: {} ({:.2}%)", 
                 self.num_singleton_buckets,
                 (self.num_singleton_buckets as f64 * 100.0) / self.num_buckets as f64);
        
        info!("  Light buckets (2-{}): {} ({:.2}%)", 
                 MIN_BUCKET_SIZE,
                 self.num_light_buckets,
                 (self.num_light_buckets as f64 * 100.0) / self.num_buckets as f64);
        
        info!("  Heavy buckets (>{}): {} ({:.2}%)", 
                 MIN_BUCKET_SIZE,
                 self.num_heavy_buckets,
                 (self.num_heavy_buckets as f64 * 100.0) / self.num_buckets as f64);
        
        info!("  Max bucket size: {}", self.max_bucket_size);
        info!("  Positions in light buckets: {} ({:.2}%)",
                 self.num_positions_in_light,
                 (self.num_positions_in_light as f64 * 100.0) / self.num_positions as f64);
        info!("  Positions in heavy buckets: {} ({:.2}%)",
                 self.num_positions_in_heavy,
                 (self.num_positions_in_heavy as f64 * 100.0) / self.num_positions as f64);
    }
}

impl Default for BucketStatistics {
    fn default() -> Self {
        Self::new()
    }
}

/// Bucket type classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BucketType {
    /// Singleton bucket (size == 1)
    Singleton,
    /// Light bucket (1 < size <= MIN_BUCKET_SIZE) - goes into sparse index
    Light,
    /// Heavy bucket (size > MIN_BUCKET_SIZE) - goes into skew index
    Heavy,
}

impl BucketType {
    /// Classify a bucket by its size
    pub fn from_bucket_size(size: usize) -> Self {
        match size {
            1 => BucketType::Singleton,
            2..=MIN_BUCKET_SIZE => BucketType::Light,
            _ => BucketType::Heavy,
        }
    }
}

/// A bucket of minimizer tuples sharing the same minimizer value
#[derive(Debug, Clone)]
pub struct Bucket {
    /// The minimizer value for this bucket
    pub minimizer: u64,

    /// All tuples with this minimizer
    pub tuples: Vec<MinimizerTuple>,

    /// Cached number of unique super-kmers (distinct pos_in_seq values)
    pub cached_size: usize,

    /// Bucket type classification
    pub bucket_type: BucketType,
}

impl Bucket {
    /// Create a new bucket
    pub fn new(minimizer: u64, tuples: Vec<MinimizerTuple>) -> Self {
        let (bucket_size, bucket_type) = compute_bucket_size_from_slice(&tuples);
        Self {
            minimizer,
            tuples,
            cached_size: bucket_size,
            bucket_type,
        }
    }

    /// Get the bucket size (number of unique super-kmers)
    #[inline]
    pub fn size(&self) -> usize {
        self.cached_size
    }

    /// Get the bucket type
    #[inline]
    pub fn bucket_type(&self) -> BucketType {
        self.bucket_type
    }
}

/// Classify sorted minimizer tuples into buckets by minimizer value
///
/// Takes a sorted vector of MinimizerTuples and groups them by minimizer,
/// creating one bucket per unique minimizer value.
///
/// # Arguments
/// * `tuples` - Sorted minimizer tuples (sorted by minimizer, then pos_in_seq)
///
/// # Returns
/// A vector of buckets, one per unique minimizer
pub fn classify_into_buckets(tuples: Vec<MinimizerTuple>) -> Vec<Bucket> {
    if tuples.is_empty() {
        return Vec::new();
    }

    let mut buckets = Vec::new();
    let mut current_minimizer = tuples[0].minimizer;
    let mut current_bucket_tuples = Vec::new();

    for tuple in tuples {
        if tuple.minimizer != current_minimizer {
            // Start new bucket
            buckets.push(Bucket::new(current_minimizer, current_bucket_tuples));
            current_minimizer = tuple.minimizer;
            current_bucket_tuples = Vec::new();
        }
        current_bucket_tuples.push(tuple);
    }

    // Push the last bucket
    if !current_bucket_tuples.is_empty() {
        buckets.push(Bucket::new(current_minimizer, current_bucket_tuples));
    }

    buckets
}

// ---------------------------------------------------------------------------
// In-place (zero-copy) bucket classification
// ---------------------------------------------------------------------------

/// Lightweight bucket descriptor referencing a range in a sorted tuples array.
///
/// Unlike [`Bucket`], this does not own the tuples — it references them by
/// index into the parent [`ClassifiedBuckets::tuples`] vector.
#[derive(Debug, Clone, Copy)]
pub struct BucketRef {
    /// The minimizer value for this bucket.
    pub minimizer: u64,
    /// Start index (inclusive) in the tuples array.
    pub start: usize,
    /// Number of tuples in this bucket.
    pub len: usize,
    /// Number of unique super-kmers (distinct `pos_in_seq` values).
    pub cached_size: usize,
    /// Bucket type classification.
    pub bucket_type: BucketType,
}

/// In-place classified buckets: owns a sorted tuple array plus lightweight
/// bucket descriptors.
///
/// This avoids the memory duplication of [`classify_into_buckets`] which
/// moves every tuple into per-bucket `Vec`s (briefly doubling memory for
/// the tuple payload).
pub struct ClassifiedBuckets {
    /// The original sorted tuples — never moved or copied.
    pub tuples: Vec<MinimizerTuple>,
    /// One descriptor per unique minimizer, in the same order as the tuples.
    pub bucket_refs: Vec<BucketRef>,
}

impl ClassifiedBuckets {
    /// Number of buckets (unique minimizers).
    #[inline]
    pub fn num_buckets(&self) -> usize {
        self.bucket_refs.len()
    }

    /// Get the tuple slice for bucket at index `idx`.
    #[inline]
    pub fn bucket_tuples(&self, idx: usize) -> &[MinimizerTuple] {
        let bref = &self.bucket_refs[idx];
        &self.tuples[bref.start..bref.start + bref.len]
    }
}

/// Classify sorted minimizer tuples into buckets **in place**, without
/// copying or moving the tuple data.
///
/// Returns a [`ClassifiedBuckets`] that owns the original tuples vector
/// and provides lightweight [`BucketRef`] descriptors.
pub fn classify_into_buckets_inplace(tuples: Vec<MinimizerTuple>) -> ClassifiedBuckets {
    if tuples.is_empty() {
        return ClassifiedBuckets {
            tuples,
            bucket_refs: Vec::new(),
        };
    }

    let mut bucket_refs = Vec::new();
    let mut start = 0usize;
    let mut current_minimizer = tuples[0].minimizer;

    for i in 1..tuples.len() {
        if tuples[i].minimizer != current_minimizer {
            let len = i - start;
            let (cached_size, bucket_type) =
                compute_bucket_size_from_slice(&tuples[start..i]);
            bucket_refs.push(BucketRef {
                minimizer: current_minimizer,
                start,
                len,
                cached_size,
                bucket_type,
            });
            current_minimizer = tuples[i].minimizer;
            start = i;
        }
    }

    // Last bucket
    let len = tuples.len() - start;
    let (cached_size, bucket_type) =
        compute_bucket_size_from_slice(&tuples[start..]);
    bucket_refs.push(BucketRef {
        minimizer: current_minimizer,
        start,
        len,
        cached_size,
        bucket_type,
    });

    ClassifiedBuckets {
        tuples,
        bucket_refs,
    }
}

/// Size and classify a bucket from a (sorted) tuple slice.
///
/// The minimizer scheme is forward, so a position is never abandoned and
/// later re-selected: each tuple carries a distinct `pos_in_seq` and the
/// number of super-kmers equals the number of tuples. This function is the
/// in-memory pipeline's enforcement point for that invariant: everything
/// downstream sizes the index by the tuple count, so a violation must stop
/// the build rather than corrupt it.
fn compute_bucket_size_from_slice(tuples: &[MinimizerTuple]) -> (usize, BucketType) {
    for pair in tuples.windows(2) {
        assert_ne!(
            pair[0].pos_in_seq, pair[1].pos_in_seq,
            "the minimizer scheme is not forward: (minimizer {:#x}, pos_in_seq {}) \
             occurs in more than one super-kmer",
            pair[0].minimizer, pair[0].pos_in_seq
        );
    }
    (tuples.len(), BucketType::from_bucket_size(tuples.len()))
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_bucket_type_classification() {
        assert_eq!(BucketType::from_bucket_size(1), BucketType::Singleton);
        assert_eq!(BucketType::from_bucket_size(2), BucketType::Light);
        assert_eq!(BucketType::from_bucket_size(64), BucketType::Light);
        assert_eq!(BucketType::from_bucket_size(65), BucketType::Heavy);
        assert_eq!(BucketType::from_bucket_size(1000), BucketType::Heavy);
    }
    
    #[test]
    fn test_classify_into_buckets_empty() {
        let tuples = Vec::new();
        let buckets = classify_into_buckets(tuples);
        assert_eq!(buckets.len(), 0);
    }
    
    #[test]
    fn test_classify_into_buckets_single() {
        let tuples = vec![
            MinimizerTuple::new(100, 50, 0),
        ];
        let buckets = classify_into_buckets(tuples);
        assert_eq!(buckets.len(), 1);
        assert_eq!(buckets[0].minimizer, 100);
        assert_eq!(buckets[0].size(), 1);
        assert_eq!(buckets[0].bucket_type(), BucketType::Singleton);
    }
    
    #[test]
    fn test_classify_into_buckets_multiple() {
        let tuples = vec![
            MinimizerTuple::new(100, 50, 0),
            MinimizerTuple::new(100, 51, 0),
            MinimizerTuple::new(200, 100, 0),
            MinimizerTuple::new(300, 150, 0),
            MinimizerTuple::new(300, 151, 0),
            MinimizerTuple::new(300, 152, 0),
        ];
        
        let buckets = classify_into_buckets(tuples);
        assert_eq!(buckets.len(), 3);
        
        assert_eq!(buckets[0].minimizer, 100);
        assert_eq!(buckets[0].size(), 2);
        assert_eq!(buckets[0].bucket_type(), BucketType::Light);
        
        assert_eq!(buckets[1].minimizer, 200);
        assert_eq!(buckets[1].size(), 1);
        assert_eq!(buckets[1].bucket_type(), BucketType::Singleton);
        
        assert_eq!(buckets[2].minimizer, 300);
        assert_eq!(buckets[2].size(), 3);
        assert_eq!(buckets[2].bucket_type(), BucketType::Light);
    }
    
    #[test]
    fn test_bucket_statistics() {
        let mut stats = BucketStatistics::new();
        
        // Add a singleton bucket
        let bucket1 = vec![MinimizerTuple::new(100, 50, 0)];
        stats.add_bucket(&bucket1);
        
        // Add a light bucket
        let bucket2 = vec![
            MinimizerTuple::new(200, 100, 0),
            MinimizerTuple::new(200, 101, 0),
        ];
        stats.add_bucket(&bucket2);
        
        // Add a heavy bucket (65 tuples)
        let mut bucket3 = Vec::new();
        for i in 0..65 {
            bucket3.push(MinimizerTuple::new(300, 200 + i, 0));
        }
        stats.add_bucket(&bucket3);
        
        assert_eq!(stats.num_buckets, 3);
        assert_eq!(stats.num_singleton_buckets, 1);
        assert_eq!(stats.num_light_buckets, 1);
        assert_eq!(stats.num_heavy_buckets, 1);
        assert_eq!(stats.max_bucket_size, 65);
    }
}
