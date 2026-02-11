//! Minimizer tuple computation for the build pipeline
//!
//! This module extracts (minimizer, position, count) tuples from the
//! encoded strings, which are then used to build the sparse/skew index.
//!
//! **Important**: MinimizerTuple is used ONLY during index construction and is
//! NOT stored in the final dictionary. It's an intermediate representation that
//! gets converted into the sparse/skew index structure.
//!
//! ## Parallelism
//!
//! Minimizer tuple extraction is parallelized across SPSS strings using rayon.
//! Each string is processed independently (super-k-mers never cross string
//! boundaries), producing a local `Vec<MinimizerTuple>`. Results are collected
//! and sorted in parallel via `par_sort_unstable()`.
//!
//! ## External Sorting
//!
//! For large datasets that exceed RAM limits, external sorting is available via
//! [`compute_minimizer_tuples_external`]. This follows the C++ implementation:
//! - Each thread has a RAM-bounded buffer
//! - When buffer fills → parallel sort + flush to temp binary file
//! - After all tuples → k-way merge of temp files using loser tree

use std::cmp::Ordering;
use std::sync::Arc;

use rayon::prelude::*;
use tracing::info;

use crate::{
    kmer::{Kmer, KmerBits},
    minimizer::MinimizerIterator,
    spectrum_preserving_string_set::SpectrumPreservingStringSet,
};
use super::config::BuildConfiguration;
use super::external_sort::{ExternalSorter, MinimizerTupleExternal, GIB};

/// A minimizer tuple representing a k-mer's minimizer and its position
///
/// This is the fundamental unit for indexing. Each tuple represents a
/// minimizer occurrence in the input and tracks where it appears.
///
/// # Canonical Mode Orientation Handling
///
/// In canonical mode, when a k-mer's reverse complement minimizer is smaller
/// than the forward minimizer, we use the RC minimizer BUT adjust `pos_in_kmer`
/// to reflect the position in the **forward** k-mer:
///
/// ```text
/// If RC minimizer is chosen:
///   pos_in_kmer = (k - m) - original_rc_pos
/// ```
///
/// This elegantly encodes orientation implicitly without needing a separate flag.
/// The minimizer position is always relative to the forward k-mer representation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinimizerTuple {
    /// The minimizer value (m-mer)
    pub minimizer: u64,
    
    /// Position of the minimizer in the sequence (in bases from string start)
    pub pos_in_seq: u64,
    
    /// Position of the minimizer within the k-mer (0 to k-m)
    /// 
    /// In canonical mode, this is ALWAYS relative to the forward k-mer,
    /// even when the RC minimizer is chosen (adjusted accordingly).
    pub pos_in_kmer: u8,
    
    /// Number of consecutive k-mers that share this minimizer (super-k-mer size)
    pub num_kmers_in_super_kmer: u8,
}

impl MinimizerTuple {
    /// Create a new minimizer tuple
    pub fn new(minimizer: u64, pos_in_seq: u64, pos_in_kmer: u8) -> Self {
        Self {
            minimizer,
            pos_in_seq,
            pos_in_kmer,
            num_kmers_in_super_kmer: 1,
        }
    }
}

impl PartialOrd for MinimizerTuple {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MinimizerTuple {
    fn cmp(&self, other: &Self) -> Ordering {
        // Primary sort by minimizer value
        match self.minimizer.cmp(&other.minimizer) {
            Ordering::Equal => {
                // Secondary sort by position in sequence
                self.pos_in_seq.cmp(&other.pos_in_seq)
            }
            other => other,
        }
    }
}

/// Compute minimizer tuples from an encoded SPSS (parallelized with rayon)
///
/// This function extracts all k-mers from the SPSS, computes their minimizers
/// (with canonical mode handling), and produces sorted, coalesced MinimizerTuples.
///
/// # Parallelism
///
/// Strings are processed in parallel via rayon. Each string is independent:
/// super-k-mers never cross string boundaries, so per-string coalescing is
/// correct. The resulting tuples are merged and sorted in parallel using
/// `par_sort_unstable()`.
///
/// The rayon thread pool size is controlled by the caller (typically
/// `DictionaryBuilder` installs a pool sized to `config.num_threads`).
///
/// # Canonical Mode Orientation
///
/// When canonical mode is enabled, for each k-mer we compute both forward and
/// reverse-complement minimizers. If the RC minimizer is smaller, we use it BUT
/// adjust the position: `pos_in_kmer = (k - m) - rc_pos`. This ensures pos_in_kmer
/// always represents the position relative to the forward k-mer.
///
/// # Returns
///
/// A sorted vector of MinimizerTuples, with consecutive k-mers sharing the same
/// minimizer coalesced into super-k-mers (num_kmers_in_super_kmer > 1).
///
/// IMPORTANT: Coalescing happens DURING extraction (while iterating k-mers in sequence),
/// not after sorting. This matches the C++ implementation exactly.
pub fn compute_minimizer_tuples<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    config: &BuildConfiguration,
) -> Vec<MinimizerTuple>
where
    Kmer<K>: KmerBits,
{
    let k = K;
    let m = config.m;
    
    if k <= m {
        panic!("k must be > m (k={}, m={})", k, m);
    }

    let num_strings = spss.num_strings() as usize;

    // Process strings in parallel — each string produces locally coalesced tuples.
    // rayon's flat_map_iter collects per-string Vec results efficiently.
    let mut tuples: Vec<MinimizerTuple> = (0..num_strings)
        .into_par_iter()
        .flat_map_iter(|string_id| {
            extract_tuples_for_string::<K>(spss, string_id as u64, k, m, config)
        })
        .collect();

    // Parallel sort by (minimizer, pos_in_seq) — exactly like C++
    // This is the only sorting step; coalescing has already happened during extraction
    tuples.par_sort_unstable();

    tuples
}

/// Extract minimizer tuples from a single string in the SPSS
///
/// Processes one string sequentially, maintaining MinimizerIterator state
/// for efficient sliding-window minimizer computation. Coalesces consecutive
/// k-mers sharing the same minimizer into super-k-mers inline.
fn extract_tuples_for_string<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    string_id: u64,
    k: usize,
    m: usize,
    config: &BuildConfiguration,
) -> Vec<MinimizerTuple>
where
    Kmer<K>: KmerBits,
{
    let string_len = spss.string_length(string_id);

    if string_len < k {
        return Vec::new();
    }

    let num_kmers = string_len - k + 1;
    let mut tuples = Vec::new();

    // CREATE MINIMIZER ITERATOR ONCE PER STRING (matching C++ behavior!)
    // Maintains state across k-mers for efficient super-k-mer detection
    let mut minimizer_iter = MinimizerIterator::with_seed(k, m, config.seed);

    // Initialize iterator with absolute string start position (matching C++ behavior)
    let string_begin = spss.string_offset(string_id);
    minimizer_iter.set_position(string_begin);

    // COALESCING DURING EXTRACTION - matching C++ behavior exactly
    // Track the previous minimizer info to detect when it changes
    let mut prev_mini_tuple: Option<MinimizerTuple> = None;
    let mut num_kmers_in_super_kmer = 0u8;

    for kmer_pos in 0..num_kmers {
        // Decode k-mer from SPSS
        let kmer = spss.decode_kmer::<K>(string_id, kmer_pos);

        // Compute forward minimizer using PERSISTENT iterator state
        let forward_minimizer = minimizer_iter.next(kmer);

        let (final_minimizer, final_pos_in_kmer) = if config.canonical {
            // Compute RC minimizer using a FRESH iterator (not sequential).
            // RC k-mers slide in the opposite direction (new base at front,
            // not at end), so the sequential sliding optimization gives wrong
            // results. A fresh iterator always does a full rescan, matching
            // the query-time behavior exactly.
            let kmer_rc = kmer.reverse_complement();
            let mut fresh_rc_iter = MinimizerIterator::with_seed(k, m, config.seed);
            let rc_minimizer = fresh_rc_iter.next(kmer_rc);

            // Choose smaller minimizer (compare by value, which is deterministic)
            if rc_minimizer.value < forward_minimizer.value {
                // RC minimizer wins - adjust position to forward k-mer frame
                let adjusted_pos = (k - m) as u8 - rc_minimizer.pos_in_kmer as u8;
                (rc_minimizer.value, adjusted_pos)
            } else {
                (forward_minimizer.value, forward_minimizer.pos_in_kmer as u8)
            }
        } else {
            // Non-canonical mode: just use forward minimizer
            (forward_minimizer.value, forward_minimizer.pos_in_kmer as u8)
        };

        // Store absolute position of minimizer in concatenated SPSS
        // (matching C++ decoded_offsets approach)
        // absolute_pos = string_begin + kmer_pos + pos_in_kmer
        // During lookup, we recover: kmer_start = absolute_pos - pos_in_kmer
        let absolute_pos_in_seq = string_begin + kmer_pos as u64 + final_pos_in_kmer as u64;

        let current_mini = MinimizerTuple {
            minimizer: final_minimizer,
            pos_in_seq: absolute_pos_in_seq,
            pos_in_kmer: final_pos_in_kmer,
            num_kmers_in_super_kmer: 1, // Will be updated during coalescing
        };

        // Check if this minimizer matches the previous one
        // Per C++ code: check only minimizer and pos_in_seq (NOT pos_in_kmer)
        // See compute_minimizer_tuples.cpp line 99-107
        if let Some(prev) = prev_mini_tuple.take() {
            if current_mini.minimizer == prev.minimizer && current_mini.pos_in_seq == prev.pos_in_seq {
                // Same minimizer occurrence - part of same super-k-mer
                num_kmers_in_super_kmer += 1;
                prev_mini_tuple = Some(prev);
            } else {
                // Different minimizer - save previous and start a new one
                let mut saved = prev;
                saved.num_kmers_in_super_kmer = num_kmers_in_super_kmer;
                tuples.push(saved);

                // Start tracking the new minimizer
                prev_mini_tuple = Some(current_mini);
                num_kmers_in_super_kmer = 1;
            }
        } else {
            // First k-mer in this string
            prev_mini_tuple = Some(current_mini);
            num_kmers_in_super_kmer = 1;
        }
    }

    // Don't forget to save the last accumulated tuple for this string
    if let Some(mut last) = prev_mini_tuple.take() {
        last.num_kmers_in_super_kmer = num_kmers_in_super_kmer;
        tuples.push(last);
    }

    tuples
}

/// Estimate memory needed for in-memory tuple processing
///
/// Returns estimated bytes needed to hold all tuples in memory.
/// Uses heuristic: ~1 tuple per k-mer on average (conservative).
pub fn estimate_memory_bytes(total_kmers: u64) -> u64 {
    // Each tuple is roughly 24 bytes in Rust (with alignment)
    // Plus sorting needs temporary space (~2x)
    total_kmers * 24 * 2
}

/// Check if external sorting is needed based on estimated data size and RAM limit
/// 
/// Returns false if ram_limit_gib is 0 (meaning unlimited/in-memory)
pub fn needs_external_sorting(total_kmers: u64, ram_limit_gib: usize) -> bool {
    if ram_limit_gib == 0 {
        return false; // 0 means unlimited, use in-memory sorting
    }
    let estimated_bytes = estimate_memory_bytes(total_kmers);
    let ram_bytes = (ram_limit_gib as u64) * (GIB as u64);
    estimated_bytes > ram_bytes
}

/// Compute minimizer tuples using external sorting (RAM-bounded)
///
/// This follows the C++ implementation precisely:
/// - Each thread partition has its own buffer (sized by RAM limit)
/// - When buffer fills → parallel sort + flush to temp binary file
/// - After all strings processed → k-way merge of temp files
///
/// Use this for large datasets that exceed available RAM.
///
/// # Arguments
/// * `spss` - The spectrum-preserving string set
/// * `config` - Build configuration (includes RAM limit, threads, tmp dir)
///
/// # Returns
/// Sorted vector of MinimizerTuples (loaded from merged temp file)
pub fn compute_minimizer_tuples_external<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    config: &BuildConfiguration,
) -> Result<Vec<MinimizerTuple>, std::io::Error>
where
    Kmer<K>: KmerBits,
{
    let k = K;
    let m = config.m;
    
    if k <= m {
        panic!("k must be > m (k={}, m={})", k, m);
    }

    let num_threads = if config.num_threads == 0 {
        rayon::current_num_threads()
    } else {
        config.num_threads
    };
    
    let sorter = Arc::new(ExternalSorter::new(
        &config.tmp_dirname,
        config.ram_limit_gib,
        num_threads,
        config.verbose,
    )?);

    let buffer_size = sorter.buffer_size_per_thread();
    info!(
        "External sorting: {} threads, {} GiB RAM limit, {} tuples/buffer",
        num_threads, config.ram_limit_gib, buffer_size
    );

    let num_strings = spss.num_strings() as usize;
    let num_strings_per_thread = num_strings.div_ceil(num_threads);

    // Process strings in thread partitions, each with its own buffer
    // This matches C++ behavior exactly
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(num_threads);
        
        for t in 0..num_threads {
            let begin = t * num_strings_per_thread;
            if begin >= num_strings {
                break;
            }
            let end = (begin + num_strings_per_thread).min(num_strings);
            
            let sorter = Arc::clone(&sorter);
            let handle = scope.spawn(move || {
                process_string_range_external::<K>(
                    spss,
                    begin..end,
                    k,
                    m,
                    config,
                    &sorter,
                    buffer_size,
                )
            });
            handles.push(handle);
        }
        
        // Wait for all threads
        for handle in handles {
            handle.join().expect("Thread panicked")?;
        }
        
        Ok::<(), std::io::Error>(())
    })?;

    // Merge all temp files
    let _merge_result = sorter.merge()?;
    
    // Read merged tuples into memory
    let tuples = sorter.read_merged_tuples()?;
    
    // Cleanup is handled by Drop
    
    Ok(tuples)
}

/// Process a range of strings, flushing to disk when buffer fills
fn process_string_range_external<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    string_range: std::ops::Range<usize>,
    k: usize,
    m: usize,
    config: &BuildConfiguration,
    sorter: &ExternalSorter,
    buffer_size: usize,
) -> Result<(), std::io::Error>
where
    Kmer<K>: KmerBits,
{
    let mut buffer: Vec<MinimizerTupleExternal> = Vec::with_capacity(buffer_size);

    for string_id in string_range {
        let string_len = spss.string_length(string_id as u64);

        if string_len < k {
            continue;
        }

        let num_kmers = string_len - k + 1;
        
        // Create minimizer iterator for this string
        let mut minimizer_iter = MinimizerIterator::with_seed(k, m, config.seed);
        let string_begin = spss.string_offset(string_id as u64);
        minimizer_iter.set_position(string_begin);
        
        // Track previous minimizer for coalescing
        let mut prev_mini: Option<MinimizerTupleExternal> = None;
        let mut num_kmers_in_super_kmer: u8 = 0;

        for kmer_pos in 0..num_kmers {
            let kmer = spss.decode_kmer::<K>(string_id as u64, kmer_pos);
            let forward_minimizer = minimizer_iter.next(kmer);

            let (final_minimizer, final_pos_in_kmer) = if config.canonical {
                let kmer_rc = kmer.reverse_complement();
                let mut fresh_rc_iter = MinimizerIterator::with_seed(k, m, config.seed);
                let rc_minimizer = fresh_rc_iter.next(kmer_rc);

                if rc_minimizer.value < forward_minimizer.value {
                    let adjusted_pos = (k - m) as u8 - rc_minimizer.pos_in_kmer as u8;
                    (rc_minimizer.value, adjusted_pos)
                } else {
                    (forward_minimizer.value, forward_minimizer.pos_in_kmer as u8)
                }
            } else {
                (forward_minimizer.value, forward_minimizer.pos_in_kmer as u8)
            };

            let absolute_pos_in_seq = string_begin + kmer_pos as u64 + final_pos_in_kmer as u64;

            let current = MinimizerTupleExternal {
                minimizer: final_minimizer,
                pos_in_seq: absolute_pos_in_seq,
                pos_in_kmer: final_pos_in_kmer,
                num_kmers_in_super_kmer: 1,
            };

            // Coalescing logic (matching C++: check minimizer and pos_in_seq only)
            // pos_in_kmer changes with each consecutive k-mer in a super-k-mer,
            // so we must NOT compare it here — we keep the pos_in_kmer from the
            // first k-mer in the run (stored in `prev`).
            if let Some(mut prev) = prev_mini.take() {
                if current.minimizer == prev.minimizer 
                    && current.pos_in_seq == prev.pos_in_seq 
                {
                    // Same super-k-mer
                    num_kmers_in_super_kmer = num_kmers_in_super_kmer.saturating_add(1);
                    prev_mini = Some(prev);
                } else {
                    // Different - flush previous
                    prev.num_kmers_in_super_kmer = num_kmers_in_super_kmer;
                    
                    // Check if buffer is full
                    if buffer.len() >= buffer_size {
                        let mut buf = std::mem::take(&mut buffer);
                        sorter.sort_and_flush(&mut buf)?;
                        buffer = buf; // Reuse the allocation
                    }
                    buffer.push(prev);
                    
                    prev_mini = Some(current);
                    num_kmers_in_super_kmer = 1;
                }
            } else {
                prev_mini = Some(current);
                num_kmers_in_super_kmer = 1;
            }
        }

        // Flush last tuple for this string
        if let Some(mut last) = prev_mini.take() {
            last.num_kmers_in_super_kmer = num_kmers_in_super_kmer;
            
            if buffer.len() >= buffer_size {
                let mut buf = std::mem::take(&mut buffer);
                sorter.sort_and_flush(&mut buf)?;
                buffer = buf;
            }
            buffer.push(last);
        }
    }

    // Flush remaining buffer
    if !buffer.is_empty() {
        let mut buf = buffer;
        sorter.sort_and_flush(&mut buf)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_minimizer_tuple_creation() {
        let tuple = MinimizerTuple::new(12345, 100, 5);
        assert_eq!(tuple.minimizer, 12345);
        assert_eq!(tuple.pos_in_seq, 100);
        assert_eq!(tuple.pos_in_kmer, 5);
        assert_eq!(tuple.num_kmers_in_super_kmer, 1);
    }
    
    #[test]
    fn test_minimizer_tuple_ordering() {
        let t1 = MinimizerTuple::new(100, 50, 0);
        let t2 = MinimizerTuple::new(200, 50, 0);
        let t3 = MinimizerTuple::new(100, 100, 0);
        
        assert!(t1 < t2);  // Different minimizer
        assert!(t1 < t3);  // Same minimizer, different position
        assert!(t2 > t1);
    }
    
    #[test]
    fn test_compute_minimizer_tuples_basic() {
        use crate::builder::config::BuildConfiguration;
        use crate::builder::encode::Encoder;
        
        // Build a simple SPSS with a single short sequence
        let config = BuildConfiguration::new(31, 13).unwrap();
        let mut encoder = Encoder::<31>::new();
        
        // Add a sequence that's exactly 31 bases (1 k-mer)
        let sequence = "ACGTACGTACGTACGTACGTACGTACGTACG";
        encoder.add_sequence(sequence.as_bytes()).unwrap();
        
        let spss = encoder.build(config.m);
        
        // Compute minimizer tuples
        let tuples = compute_minimizer_tuples::<31>(&spss, &config);
        
        // Should have exactly 1 tuple (1 k-mer)
        assert_eq!(tuples.len(), 1);
        assert!(tuples[0].pos_in_seq < sequence.len() as u64);
        assert_eq!(tuples[0].num_kmers_in_super_kmer, 1);
    }
    
    #[test]
    fn test_compute_minimizer_tuples_canonical() {
        use crate::builder::config::BuildConfiguration;
        use crate::builder::encode::Encoder;
        
        // Build SPSS in canonical mode
        let mut config = BuildConfiguration::new(31, 13).unwrap();
        config.canonical = true;
        let mut encoder = Encoder::<31>::new();
        
        // Add a sequence
        let sequence = "ACGTACGTACGTACGTACGTACGTACGTACG";
        encoder.add_sequence(sequence.as_bytes()).unwrap();
        
        let spss = encoder.build(config.m);
        
        // Compute minimizer tuples with canonical mode
        let tuples = compute_minimizer_tuples::<31>(&spss, &config);
        
        // Should have tuples with canonical minimizers
        assert!(!tuples.is_empty());
        
        // In canonical mode, pos_in_kmer should be adjusted when RC minimizer wins
        // This is a structural test - the actual values depend on the hashing
        for tuple in &tuples {
            assert!(tuple.pos_in_kmer <= (31 - 13) as u8);
        }
    }
}
