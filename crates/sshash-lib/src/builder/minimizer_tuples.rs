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
/// The minimizer value is the canonical m-mer at the anchored locus (the
/// unified scheme), and `pos_in_kmer` is the locus in the forward frame of
/// the SPSS string; a lookup probes both `pos` and its mirror `k - m - pos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub struct MinimizerTuple {
    /// The minimizer value (m-mer)
    pub minimizer: u64,

    /// `pos_in_seq << 12 | pos_in_kmer << 6 | num_kmers_in_super_kmer`, so
    /// that a tuple is 16 bytes and ordering by `(minimizer, packed)` is
    /// ordering by `(minimizer, pos_in_seq)` (a minimizer occurs once per
    /// position). `pos_in_kmer <= k - m` and `num_kmers_in_super_kmer <=
    /// k - m + 1` fit 6 bits for k <= 63, and positions fit 52 bits.
    packed: u64,
}

impl MinimizerTuple {
    /// Largest supported `pos_in_seq` (exclusive): 2^52 bases.
    pub const MAX_POS_IN_SEQ: u64 = 1 << 52;

    /// Create a new minimizer tuple (one k-mer in its super-k-mer)
    #[inline]
    pub fn new(minimizer: u64, pos_in_seq: u64, pos_in_kmer: u8) -> Self {
        Self::with_kmers(minimizer, pos_in_seq, pos_in_kmer, 1)
    }

    /// Create a tuple for a super-k-mer of `num_kmers_in_super_kmer` k-mers
    #[inline]
    pub fn with_kmers(minimizer: u64, pos_in_seq: u64, pos_in_kmer: u8, num_kmers_in_super_kmer: u8) -> Self {
        debug_assert!(pos_in_seq < Self::MAX_POS_IN_SEQ);
        debug_assert!(pos_in_kmer < 64 && num_kmers_in_super_kmer < 64);
        Self {
            minimizer,
            packed: pos_in_seq << 12 | (pos_in_kmer as u64) << 6 | num_kmers_in_super_kmer as u64,
        }
    }

    /// Position of the minimizer in the sequence (in bases from string start)
    #[inline(always)]
    pub fn pos_in_seq(&self) -> u64 {
        self.packed >> 12
    }

    /// The selected locus within the k-mer window (0 to k-m), in the
    /// forward frame of the stored string.
    #[inline(always)]
    pub fn pos_in_kmer(&self) -> u8 {
        (self.packed >> 6 & 63) as u8
    }

    /// Number of consecutive k-mers that share this minimizer (super-k-mer size)
    #[inline(always)]
    pub fn num_kmers_in_super_kmer(&self) -> u8 {
        (self.packed & 63) as u8
    }

    /// Set the super-k-mer size
    #[inline(always)]
    pub fn set_num_kmers_in_super_kmer(&mut self, n: u8) {
        debug_assert!(n < 64);
        self.packed = self.packed & !63 | n as u64;
    }
}

const _: () = assert!(std::mem::size_of::<MinimizerTuple>() == 16);

/// Compute minimizer tuples from an encoded SPSS (parallelized with rayon)
///
/// This function extracts all k-mers from the SPSS, computes their minimizers
/// (the unified canonical scheme), and produces sorted, coalesced MinimizerTuples.
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
    let total_bases = spss.total_bases();
    assert!(
        total_bases < MinimizerTuple::MAX_POS_IN_SEQ,
        "{total_bases} bases exceed the 2^52 supported by minimizer tuples"
    );

    // Strings are processed in parallel batches, a wave of batches at a time,
    // and each wave is copied (in parallel) into one vector. Tuples are never
    // all held twice, as when collecting per-thread vectors and concatenating
    // them. The vector reserves the expected tuple count plus 25% up front:
    // capacity that is never written costs no memory, while growing it would
    // copy every tuple so far. The order of tuples before sorting does not
    // matter: the sort below is a total order.
    let num_batches = (num_strings / 64).clamp(1, 4096);
    let batch_len = num_strings.div_ceil(num_batches);
    let batches: Vec<std::ops::Range<usize>> = (0..num_strings)
        .step_by(batch_len.max(1))
        .map(|lo| lo..(lo + batch_len).min(num_strings))
        .collect();
    let wave_len = batches.len().div_ceil(16).max(1);
    let total_kmers = total_bases.saturating_sub(num_strings as u64 * (k as u64 - 1));
    let expected = estimate_num_tuples(total_kmers, k, m) as f64 / 1.25; // 2/(w+1) per k-mer
    let tuples_per_base = expected / total_bases.max(1) as f64;

    let mut tuples: Vec<MinimizerTuple> = Vec::with_capacity((expected * 1.25) as usize);
    let mut bases_done = 0u64;
    for wave in batches.chunks(wave_len) {
        let outs: Vec<(Vec<MinimizerTuple>, u64)> = wave
            .par_iter()
            .map(|range| {
                let bases: u64 = range.clone().map(|s| spss.string_length(s as u64) as u64).sum();
                let mut out = Vec::with_capacity((bases as f64 * tuples_per_base * 1.25) as usize + 16);
                for string_id in range.clone() {
                    extract_tuples_for_string::<K>(spss, string_id as u64, k, m, config, &mut out);
                }
                (out, bases)
            })
            .collect();
        let wave_tuples: usize = outs.iter().map(|(o, _)| o.len()).sum();
        bases_done += outs.iter().map(|(_, b)| b).sum::<u64>();
        if tuples.len() + wave_tuples > tuples.capacity() {
            // More tuples than expected: grow to the projection from the
            // density so far (rare; costs one copy of the tuples so far)
            let projected = ((tuples.len() + wave_tuples) as f64 * total_bases as f64
                / bases_done.max(1) as f64
                * 1.1) as usize;
            tuples.reserve_exact(projected.max(tuples.len() + wave_tuples) - tuples.len());
        }
        // Copy the wave into the spare capacity, each batch at its offset
        let mut dst = &mut tuples.spare_capacity_mut()[..wave_tuples];
        let mut parts = Vec::with_capacity(outs.len());
        for (out, _) in &outs {
            let (head, rest) = dst.split_at_mut(out.len());
            parts.push((head, out.as_slice()));
            dst = rest;
        }
        parts.into_par_iter().for_each(|(dst, src)| {
            for (d, s) in dst.iter_mut().zip(src) {
                d.write(*s);
            }
        });
        // SAFETY: the `wave_tuples` slots after `len` were all written above
        unsafe { tuples.set_len(tuples.len() + wave_tuples) };
    }

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
    tuples: &mut Vec<MinimizerTuple>,
) where
    Kmer<K>: KmerBits,
{
    let string_len = spss.string_length(string_id);

    if string_len < k {
        return;
    }

    let num_kmers = string_len - k + 1;

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

        // The unified canonical scheme: one iterator, one order, both strands.
        // `position` is the absolute (tie-broken) minimizer position in the
        // concatenated SPSS; lookup recovers kmer_start = pos_in_seq - pos_in_kmer.
        let kmer_rc = kmer.reverse_complement();
        let mini = minimizer_iter.next(&kmer, &kmer_rc);
        debug_assert_eq!(
            mini.position,
            string_begin + kmer_pos as u64 + mini.pos_in_kmer as u64
        );

        // num_kmers_in_super_kmer is set during coalescing
        let current_mini = MinimizerTuple::new(mini.value, mini.position, mini.pos_in_kmer as u8);

        // Check if this minimizer matches the previous one
        // Per C++ code: check only minimizer and pos_in_seq (NOT pos_in_kmer)
        // See compute_minimizer_tuples.cpp line 99-107
        if let Some(prev) = prev_mini_tuple.take() {
            if current_mini.minimizer == prev.minimizer && current_mini.pos_in_seq() == prev.pos_in_seq() {
                // Same minimizer occurrence - part of same super-k-mer
                num_kmers_in_super_kmer += 1;
                prev_mini_tuple = Some(prev);
            } else {
                // Different minimizer - save previous and start a new one
                let mut saved = prev;
                saved.set_num_kmers_in_super_kmer(num_kmers_in_super_kmer);
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
        last.set_num_kmers_in_super_kmer(num_kmers_in_super_kmer);
        tuples.push(last);
    }
}

/// Expected number of minimizer tuples (super-k-mers), with 25% headroom.
///
/// With random minimizers a super-k-mer starts at a fraction 2/(w + 1) of the
/// k-mers, for w = k - m + 1 windows (human, k = 31, m = 19: 0.154 measured
/// against 0.143 expected).
pub fn estimate_num_tuples(total_kmers: u64, k: usize, m: usize) -> u64 {
    let w = k.saturating_sub(m) as f64 + 1.0;
    (total_kmers as f64 * 2.0 / (w + 1.0) * 1.25) as u64
}

/// Estimated peak memory of the in-memory build: 16-byte tuples, plus about
/// as much again for bucket starts, the MPHF order and the offsets.
pub fn estimate_memory_bytes(total_kmers: u64, k: usize, m: usize) -> u64 {
    estimate_num_tuples(total_kmers, k, m) * 32
}

/// Whether to build through external sorting: when the in-memory build is
/// estimated to exceed `ram_limit_gib` (0 = unlimited), or could hold more
/// tuples than the in-memory path indexes (2^32).
pub fn needs_external_sorting(total_kmers: u64, k: usize, m: usize, ram_limit_gib: usize) -> bool {
    if estimate_num_tuples(total_kmers, k, m) >= u32::MAX as u64 {
        return true;
    }
    if ram_limit_gib == 0 {
        return false; // 0 means unlimited, use in-memory sorting
    }
    let ram_bytes = (ram_limit_gib as u64) * (GIB as u64);
    estimate_memory_bytes(total_kmers, k, m) > ram_bytes
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
/// Sorted vector of MinimizerTuples (loaded from merged temp file).
/// For large datasets, prefer [`compute_minimizer_tuples_external_file`],
/// which returns a file handle instead of materializing all tuples.
pub fn compute_minimizer_tuples_external<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    config: &BuildConfiguration,
) -> Result<Vec<MinimizerTuple>, std::io::Error>
where
    Kmer<K>: KmerBits,
{
    let sorter = run_external_sort::<K>(spss, config)?;
    sorter.read_merged_tuples()
}

/// Compute minimizer tuples using external sorting, returning a file handle.
///
/// Like [`compute_minimizer_tuples_external`], but returns a
/// [`FileTuples`](super::external_sort::FileTuples)
/// instead of materializing all tuples in memory. This is the production path
/// for large datasets — tuples are accessed sequentially via buffered I/O.
pub fn compute_minimizer_tuples_external_file<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    config: &BuildConfiguration,
) -> Result<super::external_sort::FileTuples, std::io::Error>
where
    Kmer<K>: KmerBits,
{
    let sorter = run_external_sort::<K>(spss, config)?;
    sorter.open_merged_file()
}

/// Run the external sort pipeline (shared by both materialized and mmap paths).
///
/// Returns the ExternalSorter after merge is complete.
fn run_external_sort<const K: usize>(
    spss: &SpectrumPreservingStringSet,
    config: &BuildConfiguration,
) -> Result<ExternalSorter, std::io::Error>
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

    // Unwrap the Arc — we're the only holder at this point
    Arc::try_unwrap(sorter).map_err(|_| {
        std::io::Error::other("ExternalSorter still has multiple owners")
    })
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
            // The unified canonical scheme (see extract_tuples_for_string).
            let kmer_rc = kmer.reverse_complement();
            let mini = minimizer_iter.next(&kmer, &kmer_rc);
            debug_assert_eq!(
                mini.position,
                string_begin + kmer_pos as u64 + mini.pos_in_kmer as u64
            );

            let current = MinimizerTupleExternal {
                minimizer: mini.value,
                pos_in_seq: mini.position,
                pos_in_kmer: mini.pos_in_kmer as u8,
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
        assert_eq!(tuple.pos_in_seq(), 100);
        assert_eq!(tuple.pos_in_kmer(), 5);
        assert_eq!(tuple.num_kmers_in_super_kmer(), 1);
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
        assert!(tuples[0].pos_in_seq() < sequence.len() as u64);
        assert_eq!(tuples[0].num_kmers_in_super_kmer(), 1);
    }
    
    #[test]
    fn test_compute_minimizer_tuples_canonical() {
        use crate::builder::config::BuildConfiguration;
        use crate::builder::encode::Encoder;
        
        let config = BuildConfiguration::new(31, 13).unwrap();
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
            assert!(tuple.pos_in_kmer() <= (31 - 13) as u8);
        }
    }
}
