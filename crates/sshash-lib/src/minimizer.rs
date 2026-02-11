//! Minimizer extraction and iteration
//!
//! This module implements the core algorithm for extracting minimizers from k-mers.
//! A minimizer is the smallest m-mer (subsequence of length m) within a k-mer,
//! where "smallest" is determined by hashing.

use crate::kmer::Kmer;
use crate::hasher::DeterministicHasher;

/// Information about a minimizer within a k-mer
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinimizerInfo {
    /// The minimizer value (encoded as bits)
    pub value: u64,
    /// The absolute position in the full string (not within the k-mer)
    pub position: u64,
    /// The position of the minimizer within the current k-mer (0 = rightmost)
    pub pos_in_kmer: usize,
}

impl MinimizerInfo {
    /// Create a new minimizer info
    pub fn new(value: u64, position: u64, pos_in_kmer: usize) -> Self {
        Self {
            value,
            position,
            pos_in_kmer,
        }
    }
}

/// Iterator that extracts minimizers from k-mers using the "re-scan" method
///
/// This is an efficient minimizer extraction algorithm that:
/// - Rescans the window when the minimum falls out
/// - Performs single comparisons for sliding updates
/// - Tracks the position of the minimum within the k-mer window
///
/// # Example
/// ```
/// use sshash_lib::kmer::Kmer;
/// use sshash_lib::minimizer::MinimizerIterator;
///
/// // Create iterator for k-mer size 31, minimizer size 13
/// let mut iter = MinimizerIterator::new(31, 13, 0);
///
/// // Process a k-mer and extract minimizer
/// let kmer: Kmer<31> = Kmer::from_str("ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
/// let mini = iter.next(kmer);
/// println!("Minimizer at position {} in k-mer, global position {}", 
///     mini.pos_in_kmer, mini.position);
/// ```
pub struct MinimizerIterator {
    k: usize,
    m: usize,
    position: u64,
    min_pos_in_kmer: usize,
    min_value: u64,
    min_position: u64,
    min_hash: u64,
    hasher: DeterministicHasher,
}

impl MinimizerIterator {
    /// Create a new minimizer iterator with a seeded hasher
    ///
    /// # Arguments
    /// * `k` - k-mer size (must be >= m)
    /// * `m` - minimizer size (must be <= k)
    /// * `seed` - seed for the hash function (for deterministic results)
    pub fn with_seed(k: usize, m: usize, seed: u64) -> Self {
        assert!(k > 0 && m <= k, "k must be > 0 and m <= k");
        let hasher = DeterministicHasher::new(seed);
        let mut iter = Self {
            k,
            m,
            position: 0,
            min_pos_in_kmer: 0,
            min_value: 0,
            min_position: 0,
            min_hash: u64::MAX,
            hasher,
        };
        iter.reset();
        iter
    }

    /// Create a new minimizer iterator (defaults to seed=1)
    ///
    /// # Arguments
    /// * `k` - k-mer size (must be >= m)
    /// * `m` - minimizer size (must be <= k)
    /// * `_position` - (deprecated parameter for API compatibility)
    pub fn new(k: usize, m: usize, _position: u64) -> Self {
        Self::with_seed(k, m, 1)
    }

    /// Set a new starting position
    pub fn set_position(&mut self, position: u64) {
        self.position = position;
        self.reset();
    }

    /// Reset internal state (called when position changes)
    fn reset(&mut self) {
        // Match C++ behavior: set min_pos_in_kmer to trigger rescan on first next()
        // and set min_position to position - 1 (wrapping, matching C++ unsigned arithmetic)
        // so that next() will compute begin = min_position + 1 = position
        self.min_pos_in_kmer = 0;
        self.min_position = self.position.wrapping_sub(1);
        self.min_hash = u64::MAX;
    }

    /// Compute hash of a u64 value using the seeded hasher
    fn hash_u64(&self, value: u64) -> u64 {
        self.hasher.hash_u64(value)
    }

    /// Extract the next minimizer from a k-mer
    ///
    /// For sliding through a sequence, call this repeatedly with consecutive k-mers
    /// to efficiently compute minimizers at each position.
    pub fn next<const K: usize>(&mut self, kmer: Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: crate::kmer::KmerBits,
    {
        assert_eq!(K, self.k, "k-mer size must match iterator k");

        if self.min_pos_in_kmer == 0 {
            // Minimum fell out of window: rescan to find new minimum
            // This matches C++ logic in minimizer_iterator.hpp line 35-36
            // Use wrapping_add because min_position starts at u64::MAX after reset
            self.position = self.min_position.wrapping_add(1);
            self.rescan(kmer);
        } else {
            // Sliding: check new m-mer at the end of the k-mer
            self.position += 1;

            // Extract the rightmost m-mer
            let mmer_value = self.extract_mmer(kmer, self.k - self.m);
            let hash = self.hash_u64(mmer_value);

            if hash < self.min_hash {
                // Found a new minimum (always at position k-m in new k-mer)
                self.min_hash = hash;
                self.min_value = mmer_value;
                self.min_position = self.position;
                self.min_pos_in_kmer = self.k - self.m;
            } else {
                // Minimum is still in window, just shifted left
                self.min_pos_in_kmer -= 1;
            }
        }

        MinimizerInfo::new(self.min_value, self.min_position, self.min_pos_in_kmer)
    }

    /// Rescan the k-mer window to find all m-mers and their minimum
    /// This matches the C++ implementation in minimizer_iterator.hpp lines 61-73
    fn rescan<const K: usize>(&mut self, kmer: Kmer<K>)
    where
        Kmer<K>: crate::kmer::KmerBits,
    {
        self.min_hash = u64::MAX;
        self.min_value = 0;
        self.min_pos_in_kmer = 0;

        let begin = self.position;
        
        // Try each m-mer position in the k-mer window
        for i in 0..=(self.k - self.m) {
            let mmer_value = self.extract_mmer(kmer, i);
            let hash = self.hash_u64(mmer_value);

            if hash < self.min_hash {
                // Leftmost minimum wins ties (matching C++ behavior)
                self.min_hash = hash;
                self.min_value = mmer_value;
                self.min_pos_in_kmer = i;
            }
        }
        
        // Set position to represent the position of the rightmost m-mer after rescan
        // This matches C++ logic: m_position = begin + (k - m) at end of rescan
        self.position = begin + (self.k - self.m) as u64;
        self.min_position = begin + self.min_pos_in_kmer as u64;
    }

    /// Extract m bases starting at position `start_pos` from a k-mer
    /// Returns the m-mer as a u64
    #[inline]
    fn extract_mmer<const K: usize>(&self, kmer: Kmer<K>, start_pos: usize) -> u64
    where
        Kmer<K>: crate::kmer::KmerBits,
    {
        // Extract m bases as a single shift+mask on the native storage type
        let shift = start_pos * 2;
        let mask = (1u64 << (self.m * 2)) - 1;
        let bits_u128 = <Kmer<K> as crate::kmer::KmerBits>::to_u128(kmer.bits());
        ((bits_u128 >> shift) as u64) & mask
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::Kmer;

    #[test]
    fn test_minimizer_iterator_basic() {
        // Create a simple k-mer: ACGTACGTACGTACGT (K=15, so we'll use K=15)
        let kmer: Kmer<15> = Kmer::from_str("ACGTACGTACGTACG").unwrap();

        let mut iter = MinimizerIterator::new(15, 7, 0);
        let mini = iter.next(kmer);

        // Verify it found a minimizer with valid properties
        assert!(mini.pos_in_kmer <= 15 - 7); // pos_in_kmer is 0..=(k-m)
        // position = pos_in_kmer for the first k-mer (min_position = begin + pos_in_kmer, begin=0)
        assert_eq!(mini.position, mini.pos_in_kmer as u64);
    }

    #[test]
    fn test_minimizer_iterator_position() {
        // Note: This test verifies minimizer iteration works correctly.
        // The exact position values depend on hash values of m-mers.
        let kmer1: Kmer<9> = Kmer::from_str("ACGTACGTA").unwrap();
        let kmer2: Kmer<9> = Kmer::from_str("CGTACGTAC").unwrap();

        let mut iter = MinimizerIterator::new(9, 5, 0);
        let mini1 = iter.next(kmer1);
        let mini2 = iter.next(kmer2);

        // Verify valid results were returned
        assert!(mini1.pos_in_kmer < 9);
        assert!(mini2.pos_in_kmer < 9);
    }

    #[test]
    fn test_minimizer_info_creation() {
        let mini = MinimizerInfo::new(42, 100, 5);
        assert_eq!(mini.value, 42);
        assert_eq!(mini.position, 100);
        assert_eq!(mini.pos_in_kmer, 5);
    }

    #[test]
    fn test_minimizer_consistency_fresh_vs_sequential() {
        // This tests that a fresh iterator gives the same minimizer as
        // processing k-mers sequentially through a string.
        // Uses the failing k-mer from test_small.fa seq1 position 47.
        
        let seq = "ATTTTCAGGATGTTTTCAGGTTCATCATCTCCCTTCTTTGCAGGATAGTAGATAAGATCGCTCATCAACGGATGTTGTGT";
        let k = 31usize;
        let m = 19usize;
        let seed = 1u64;
        let num_kmers = seq.len() - k + 1;
        
        // Process ALL k-mers sequentially, like the build does
        let mut seq_iter = MinimizerIterator::with_seed(k, m, seed);
        seq_iter.set_position(0);
        
        for kmer_pos in 0..num_kmers {
            let kmer_str = &seq[kmer_pos..kmer_pos+k];
            let kmer = Kmer::<31>::from_str(kmer_str).unwrap();
            
            let seq_mini = seq_iter.next(kmer);
            
            // Now compute with a FRESH iterator (like query does)
            let mut fresh_iter = MinimizerIterator::with_seed(k, m, seed);
            let fresh_mini = fresh_iter.next(kmer);
            
            assert_eq!(
                seq_mini.value, fresh_mini.value,
                "Minimizer VALUE mismatch at kmer_pos={}: seq={} fresh={} kmer={}",
                kmer_pos, seq_mini.value, fresh_mini.value, kmer_str
            );
            assert_eq!(
                seq_mini.pos_in_kmer, fresh_mini.pos_in_kmer,
                "Minimizer POS mismatch at kmer_pos={}: seq={} fresh={} kmer={}",
                kmer_pos, seq_mini.pos_in_kmer, fresh_mini.pos_in_kmer, kmer_str
            );
        }
    }
}
