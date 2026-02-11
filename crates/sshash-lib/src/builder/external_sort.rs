//! External sorting for minimizer tuples
//!
//! This module implements RAM-bounded external sorting for large datasets,
//! following the C++ sshash implementation precisely:
//!
//! 1. Each thread has a RAM-bounded buffer
//! 2. When buffer fills → parallel sort + flush to temp binary file
//! 3. After all tuples processed → k-way merge of temp files
//! 4. Merge uses loser tree for >16 files, linear scan for ≤16 files
//!
//! ## Memory Management
//!
//! Buffer size per thread = `(ram_limit_gib * GiB) / (2 * sizeof(tuple) * num_threads)`
//!
//! The factor of 2 accounts for temporary memory during parallel sort.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

use memmap2::Mmap;
use rayon::prelude::*;
use tracing::{debug, info};

use super::minimizer_tuples::MinimizerTuple;

/// Size of `MinimizerTupleExternal` in bytes (packed, no padding)
pub const TUPLE_SIZE_BYTES: usize = 18;

/// Bytes per GiB
pub const GIB: usize = 1024 * 1024 * 1024;

/// Packed minimizer tuple for disk I/O (matches C++ `#pragma pack(push, 2)`)
///
/// Layout: minimizer (8) + pos_in_seq (8) + pos_in_kmer (1) + num_kmers_in_super_kmer (1) = 18 bytes
#[repr(C, packed(2))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MinimizerTupleExternal {
    /// The minimizer hash value
    pub minimizer: u64,
    /// Position in the sequence
    pub pos_in_seq: u64,
    /// Position of the minimizer within the k-mer
    pub pos_in_kmer: u8,
    /// Number of k-mers in the super-k-mer
    pub num_kmers_in_super_kmer: u8,
}

impl MinimizerTupleExternal {
    /// Convert from internal MinimizerTuple
    pub fn from_internal(t: &MinimizerTuple) -> Self {
        Self {
            minimizer: t.minimizer,
            pos_in_seq: t.pos_in_seq,
            pos_in_kmer: t.pos_in_kmer,
            num_kmers_in_super_kmer: t.num_kmers_in_super_kmer,
        }
    }

    /// Convert to internal MinimizerTuple
    pub fn to_internal(&self) -> MinimizerTuple {
        MinimizerTuple {
            minimizer: self.minimizer,
            pos_in_seq: self.pos_in_seq,
            pos_in_kmer: self.pos_in_kmer,
            num_kmers_in_super_kmer: self.num_kmers_in_super_kmer,
        }
    }

    /// Read from bytes (unsafe, assumes correct alignment)
    ///
    /// # Safety
    /// Caller must ensure `bytes` points to a valid `MinimizerTupleExternal`
    #[inline]
    pub unsafe fn from_bytes(bytes: *const u8) -> Self {
        // SAFETY: read_unaligned handles packed/unaligned access
        unsafe { std::ptr::read_unaligned(bytes as *const Self) }
    }

    /// Write to bytes
    pub fn to_bytes(&self) -> [u8; TUPLE_SIZE_BYTES] {
        let mut buf = [0u8; TUPLE_SIZE_BYTES];
        unsafe {
            std::ptr::copy_nonoverlapping(
                self as *const Self as *const u8,
                buf.as_mut_ptr(),
                TUPLE_SIZE_BYTES,
            );
        }
        buf
    }
}

impl PartialOrd for MinimizerTupleExternal {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MinimizerTupleExternal {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // Copy fields to avoid taking references to packed struct fields
        let self_min = self.minimizer;
        let other_min = other.minimizer;
        let self_pos = self.pos_in_seq;
        let other_pos = other.pos_in_seq;
        
        match self_min.cmp(&other_min) {
            std::cmp::Ordering::Equal => self_pos.cmp(&other_pos),
            ord => ord,
        }
    }
}

/// External sorter for minimizer tuples
///
/// Manages RAM-bounded sorting with temp file spillover and k-way merge.
pub struct ExternalSorter {
    /// Temp directory for intermediate files
    tmp_dir: PathBuf,
    /// Run identifier (timestamp-based for uniqueness)
    run_id: u64,
    /// Atomic counter for temp file IDs
    num_files: AtomicU64,
    /// RAM limit in GiB
    ram_limit_gib: usize,
    /// Number of threads
    num_threads: usize,
    /// Verbose logging
    verbose: bool,
}

impl ExternalSorter {
    /// Create a new external sorter
    pub fn new(tmp_dir: impl AsRef<Path>, ram_limit_gib: usize, num_threads: usize, verbose: bool) -> std::io::Result<Self> {
        let tmp_dir = tmp_dir.as_ref().to_path_buf();
        
        // Create temp directory if it doesn't exist
        fs::create_dir_all(&tmp_dir)?;
        
        // Generate unique run ID from timestamp
        let run_id = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos() as u64;
        
        Ok(Self {
            tmp_dir,
            run_id,
            num_files: AtomicU64::new(0),
            ram_limit_gib,
            num_threads,
            verbose,
        })
    }

    /// Calculate buffer size per thread in number of tuples
    ///
    /// Formula: `(ram_limit * GiB) / (2 * sizeof(tuple) * num_threads)`
    /// The factor of 2 accounts for temporary memory during parallel sort.
    pub fn buffer_size_per_thread(&self) -> usize {
        let total_bytes = self.ram_limit_gib * GIB;
        let bytes_per_thread = total_bytes / (2 * self.num_threads.max(1));
        bytes_per_thread / TUPLE_SIZE_BYTES
    }

    /// Get path for a temp file by ID
    fn temp_file_path(&self, id: u64) -> PathBuf {
        self.tmp_dir.join(format!(
            "sshash.tmp.run_{}.minimizers.{}.bin",
            self.run_id, id
        ))
    }

    /// Get path for the final merged file
    fn merged_file_path(&self) -> PathBuf {
        self.tmp_dir.join(format!(
            "sshash.tmp.run_{}.minimizers.bin",
            self.run_id
        ))
    }

    /// Sort a buffer and flush to a temp file
    ///
    /// Returns the file ID. Thread-safe via atomic counter.
    pub fn sort_and_flush(&self, buffer: &mut Vec<MinimizerTupleExternal>) -> std::io::Result<u64> {
        // Parallel sort
        buffer.par_sort_unstable();
        
        // Get unique file ID
        let file_id = self.num_files.fetch_add(1, Ordering::SeqCst);
        let path = self.temp_file_path(file_id);
        
        if self.verbose {
            debug!("Flushing {} tuples to {:?}", buffer.len(), path);
        }
        
        // Write to binary file
        let file = File::create(&path)?;
        let mut writer = BufWriter::with_capacity(1024 * 1024, file);
        
        for tuple in buffer.iter() {
            writer.write_all(&tuple.to_bytes())?;
        }
        
        writer.flush()?;
        buffer.clear();
        
        Ok(file_id)
    }

    /// Number of temp files created
    pub fn num_files(&self) -> u64 {
        self.num_files.load(Ordering::SeqCst)
    }

    /// Merge all temp files into final sorted output
    ///
    /// Returns statistics: (num_minimizers, num_positions, num_super_kmers)
    pub fn merge(&self) -> std::io::Result<MergeResult> {
        let num_files = self.num_files();
        
        if num_files == 0 {
            return Ok(MergeResult::default());
        }
        
        if num_files == 1 {
            // Just rename the single file
            let src = self.temp_file_path(0);
            let dst = self.merged_file_path();
            fs::rename(&src, &dst)?;
            return self.scan_merged_file();
        }
        
        // Multiple files: k-way merge
        info!("Merging {} temp files...", num_files);
        
        let mut merger = FileMergingIterator::new(
            (0..num_files).map(|id| self.temp_file_path(id)).collect(),
        )?;
        
        let merged_path = self.merged_file_path();
        let file = File::create(&merged_path)?;
        let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
        
        let mut result = MergeResult::default();
        let mut prev_minimizer = u64::MAX;
        let mut prev_pos_in_seq = u64::MAX;
        
        while merger.has_next() {
            let tuple = merger.current();
            
            // Track statistics
            if tuple.minimizer != prev_minimizer {
                prev_minimizer = tuple.minimizer;
                result.num_minimizers += 1;
                result.num_positions += 1;
            } else if tuple.pos_in_seq != prev_pos_in_seq {
                result.num_positions += 1;
            }
            prev_pos_in_seq = tuple.pos_in_seq;
            result.num_super_kmers += 1;
            
            writer.write_all(&tuple.to_bytes())?;
            merger.next();
            
            if self.verbose && result.num_super_kmers % 100_000_000 == 0 {
                info!("Merged {} tuples...", result.num_super_kmers);
            }
        }
        
        writer.flush()?;
        drop(merger);
        
        // Remove temp files
        for id in 0..num_files {
            let _ = fs::remove_file(self.temp_file_path(id));
        }
        
        info!(
            "Merge complete: {} minimizers, {} positions, {} super-kmers",
            result.num_minimizers, result.num_positions, result.num_super_kmers
        );
        
        Ok(result)
    }

    /// Scan merged file to compute statistics (for single-file case)
    fn scan_merged_file(&self) -> std::io::Result<MergeResult> {
        let path = self.merged_file_path();
        let file = File::open(&path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        
        let mut result = MergeResult::default();
        let mut prev_minimizer = u64::MAX;
        let mut prev_pos_in_seq = u64::MAX;
        
        let num_tuples = mmap.len() / TUPLE_SIZE_BYTES;
        for i in 0..num_tuples {
            let offset = i * TUPLE_SIZE_BYTES;
            let tuple = unsafe { MinimizerTupleExternal::from_bytes(mmap.as_ptr().add(offset)) };
            
            if tuple.minimizer != prev_minimizer {
                prev_minimizer = tuple.minimizer;
                result.num_minimizers += 1;
                result.num_positions += 1;
            } else if tuple.pos_in_seq != prev_pos_in_seq {
                result.num_positions += 1;
            }
            prev_pos_in_seq = tuple.pos_in_seq;
            result.num_super_kmers += 1;
        }
        
        Ok(result)
    }

    /// Read merged tuples into memory as internal MinimizerTuples
    ///
    /// Call this after `merge()` to get the final sorted tuples.
    pub fn read_merged_tuples(&self) -> std::io::Result<Vec<MinimizerTuple>> {
        let path = self.merged_file_path();
        let file = File::open(&path)?;
        let mmap = unsafe { Mmap::map(&file)? };
        
        let num_tuples = mmap.len() / TUPLE_SIZE_BYTES;
        let mut tuples = Vec::with_capacity(num_tuples);
        
        for i in 0..num_tuples {
            let offset = i * TUPLE_SIZE_BYTES;
            let ext = unsafe { MinimizerTupleExternal::from_bytes(mmap.as_ptr().add(offset)) };
            tuples.push(ext.to_internal());
        }
        
        Ok(tuples)
    }

    /// Remove the merged file (cleanup)
    pub fn remove_merged_file(&self) -> std::io::Result<()> {
        let path = self.merged_file_path();
        if path.exists() {
            fs::remove_file(path)?;
        }
        Ok(())
    }
}

impl Drop for ExternalSorter {
    fn drop(&mut self) {
        // Clean up any remaining temp files
        for id in 0..self.num_files() {
            let _ = fs::remove_file(self.temp_file_path(id));
        }
        let _ = fs::remove_file(self.merged_file_path());
    }
}

/// Result of merge operation
#[derive(Debug, Default, Clone, Copy)]
pub struct MergeResult {
    /// Number of distinct minimizers
    pub num_minimizers: u64,
    /// Total number of positions
    pub num_positions: u64,
    /// Number of super-k-mers
    pub num_super_kmers: u64,
}

/// K-way merge iterator using loser tree
///
/// For ≤16 files: uses linear scan (simpler, cache-friendly for small N)
/// For >16 files: uses loser tree (O(log N) per element)
///
/// Files are memory-mapped for efficient access.
struct FileMergingIterator {
    /// Memory-mapped files (kept alive for iterator pointers)
    #[allow(dead_code)]
    mmaps: Vec<Mmap>,
    /// (current_ptr, end_ptr) for each file
    iterators: Vec<(*const u8, *const u8)>,
    /// Loser tree (only used for >16 files)
    tree: Vec<u32>,
    /// Tree parameters
    tree_begin: usize,
    tree_size: usize,
    /// Current minimum index
    min_idx: usize,
    /// Number of active files
    num_active: usize,
    /// Threshold for switching between linear scan and loser tree
    scan_threshold: usize,
}

impl FileMergingIterator {
    const SCAN_THRESHOLD: usize = 16;

    fn new(paths: Vec<PathBuf>) -> std::io::Result<Self> {
        let num_files = paths.len();
        if num_files == 0 {
            return Ok(Self {
                mmaps: Vec::new(),
                iterators: Vec::new(),
                tree: Vec::new(),
                tree_begin: 0,
                tree_size: 0,
                min_idx: 0,
                num_active: 0,
                scan_threshold: Self::SCAN_THRESHOLD,
            });
        }

        let mut mmaps = Vec::with_capacity(num_files);
        let mut iterators = Vec::with_capacity(num_files);

        for path in &paths {
            let file = File::open(path)?;
            let mmap = unsafe { Mmap::map(&file)? };
            let begin = mmap.as_ptr();
            let end = unsafe { begin.add(mmap.len()) };
            iterators.push((begin, end));
            mmaps.push(mmap);
        }

        let mut merger = Self {
            mmaps,
            iterators,
            tree: Vec::new(),
            tree_begin: 0,
            tree_size: 0,
            min_idx: 0,
            num_active: num_files,
            scan_threshold: Self::SCAN_THRESHOLD,
        };

        if num_files <= Self::SCAN_THRESHOLD {
            merger.compute_min_linear();
        } else {
            // Build loser tree
            let n = num_files;
            let m = 2 * n - 1;
            merger.tree_size = n;
            merger.tree.resize(m, 0);
            merger.tree_begin = (1u64 << (n as f64).log2().ceil() as u64) as usize - 1;
            
            // Initialize leaves
            let mut i = 0;
            while merger.tree_begin + i < m {
                merger.tree[merger.tree_begin + i] = i as u32;
                i += 1;
            }
            let mut j = 0;
            while i < n {
                merger.tree[n - 1 + j] = i as u32;
                i += 1;
                j += 1;
            }
            
            // Build tree bottom-up
            merger.build_tree(0);
            merger.min_idx = merger.tree[0] as usize;
        }

        Ok(merger)
    }

    fn has_next(&self) -> bool {
        self.num_active > 0
    }

    fn current(&self) -> MinimizerTupleExternal {
        debug_assert!(self.num_active > 0);
        let (ptr, _) = self.iterators[self.min_idx];
        unsafe { MinimizerTupleExternal::from_bytes(ptr) }
    }

    fn next(&mut self) {
        if self.num_active == 0 {
            return;
        }

        if self.iterators.len() <= self.scan_threshold {
            self.update_linear();
        } else {
            self.update_loser_tree();
        }
    }

    fn update_linear(&mut self) {
        let (ptr, end) = &mut self.iterators[self.min_idx];
        *ptr = unsafe { ptr.add(TUPLE_SIZE_BYTES) };
        
        if *ptr == *end {
            // This file is exhausted
            self.iterators.swap_remove(self.min_idx);
            self.num_active -= 1;
            if self.num_active == 0 {
                return;
            }
        }
        
        self.compute_min_linear();
    }

    fn compute_min_linear(&mut self) {
        if self.iterators.is_empty() {
            return;
        }
        
        self.min_idx = 0;
        let mut min_tuple = unsafe { MinimizerTupleExternal::from_bytes(self.iterators[0].0) };
        
        for (i, &(ptr, _)) in self.iterators.iter().enumerate().skip(1) {
            let tuple = unsafe { MinimizerTupleExternal::from_bytes(ptr) };
            if tuple < min_tuple {
                min_tuple = tuple;
                self.min_idx = i;
            }
        }
    }

    fn update_loser_tree(&mut self) {
        self.min_idx = self.tree[0] as usize;
        
        let (ptr, end) = &mut self.iterators[self.min_idx];
        *ptr = unsafe { ptr.add(TUPLE_SIZE_BYTES) };
        
        // Calculate leaf position
        let mut p = self.tree_begin + self.min_idx;
        if p >= self.tree.len() {
            p -= self.tree_size;
        }
        
        if *ptr == *end {
            // Mark this file as exhausted
            self.tree[p] = u32::MAX;
            self.num_active -= 1;
        }
        
        // Propagate up the tree
        while p > 0 {
            let is_right_child = (p & 1) == 0;
            let sibling = if is_right_child { p - 1 } else { p + 1 };
            
            let l_idx = if is_right_child { self.tree[sibling] } else { self.tree[p] };
            let r_idx = if is_right_child { self.tree[p] } else { self.tree[sibling] };
            
            let winner = if l_idx == u32::MAX {
                r_idx
            } else if r_idx == u32::MAX {
                l_idx
            } else {
                let l_tuple = unsafe { MinimizerTupleExternal::from_bytes(self.iterators[l_idx as usize].0) };
                let r_tuple = unsafe { MinimizerTupleExternal::from_bytes(self.iterators[r_idx as usize].0) };
                if l_tuple < r_tuple { l_idx } else { r_idx }
            };
            
            let parent = (p - 1) / 2;
            self.tree[parent] = winner;
            p = parent;
        }
        
        self.min_idx = self.tree[0] as usize;
    }

    fn build_tree(&mut self, p: usize) -> u32 {
        if p >= self.tree.len() {
            return u32::MAX;
        }
        if p >= self.tree_size - 1 {
            // Leaf
            return self.tree[p];
        }
        
        let l = self.build_tree(2 * p + 1);
        let r = self.build_tree(2 * p + 2);
        
        let winner = if l == u32::MAX {
            r
        } else if r == u32::MAX {
            l
        } else {
            let l_tuple = unsafe { MinimizerTupleExternal::from_bytes(self.iterators[l as usize].0) };
            let r_tuple = unsafe { MinimizerTupleExternal::from_bytes(self.iterators[r as usize].0) };
            if l_tuple < r_tuple { l } else { r }
        };
        
        self.tree[p] = winner;
        winner
    }
}

// Safety: The raw pointers in iterators point to memory-mapped data that
// lives as long as the mmaps vec. We only access via the iterator methods.
unsafe impl Send for FileMergingIterator {}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_tuple_external_size() {
        // Verify packed size
        assert_eq!(std::mem::size_of::<MinimizerTupleExternal>(), TUPLE_SIZE_BYTES);
    }

    #[test]
    fn test_tuple_roundtrip() {
        let tuple = MinimizerTupleExternal {
            minimizer: 12345,
            pos_in_seq: 67890,
            pos_in_kmer: 5,
            num_kmers_in_super_kmer: 3,
        };
        
        let bytes = tuple.to_bytes();
        let recovered = unsafe { MinimizerTupleExternal::from_bytes(bytes.as_ptr()) };
        
        assert_eq!(tuple, recovered);
    }

    #[test]
    fn test_external_sorter_basic() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();
        
        // Buffer size should be reasonable
        let buf_size = sorter.buffer_size_per_thread();
        assert!(buf_size > 0);
        
        // Test sort_and_flush
        let mut buffer: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal {
                minimizer: 100,
                pos_in_seq: 10,
                pos_in_kmer: 1,
                num_kmers_in_super_kmer: 2,
            },
            MinimizerTupleExternal {
                minimizer: 50,
                pos_in_seq: 20,
                pos_in_kmer: 3,
                num_kmers_in_super_kmer: 1,
            },
        ];
        
        sorter.sort_and_flush(&mut buffer).unwrap();
        assert!(buffer.is_empty());
        assert_eq!(sorter.num_files(), 1);
    }

    #[test]
    fn test_external_sorter_merge() {
        let tmp_dir = TempDir::new().unwrap();
        let sorter = ExternalSorter::new(tmp_dir.path(), 1, 2, false).unwrap();
        
        // Create multiple temp files
        let mut buffer1: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal { minimizer: 10, pos_in_seq: 1, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 },
            MinimizerTupleExternal { minimizer: 30, pos_in_seq: 3, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 },
        ];
        sorter.sort_and_flush(&mut buffer1).unwrap();
        
        let mut buffer2: Vec<MinimizerTupleExternal> = vec![
            MinimizerTupleExternal { minimizer: 20, pos_in_seq: 2, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 },
            MinimizerTupleExternal { minimizer: 40, pos_in_seq: 4, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 },
        ];
        sorter.sort_and_flush(&mut buffer2).unwrap();
        
        assert_eq!(sorter.num_files(), 2);
        
        // Merge
        let result = sorter.merge().unwrap();
        assert_eq!(result.num_super_kmers, 4);
        assert_eq!(result.num_minimizers, 4);
        
        // Read merged tuples
        let tuples = sorter.read_merged_tuples().unwrap();
        assert_eq!(tuples.len(), 4);
        
        // Verify sorted order
        assert_eq!(tuples[0].minimizer, 10);
        assert_eq!(tuples[1].minimizer, 20);
        assert_eq!(tuples[2].minimizer, 30);
        assert_eq!(tuples[3].minimizer, 40);
    }

    #[test]
    fn test_tuple_ordering() {
        let t1 = MinimizerTupleExternal { minimizer: 100, pos_in_seq: 50, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 };
        let t2 = MinimizerTupleExternal { minimizer: 100, pos_in_seq: 60, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 };
        let t3 = MinimizerTupleExternal { minimizer: 200, pos_in_seq: 10, pos_in_kmer: 0, num_kmers_in_super_kmer: 1 };
        
        assert!(t1 < t2);  // Same minimizer, different pos
        assert!(t1 < t3);  // Different minimizer
        assert!(t2 < t3);  // Different minimizer
    }
}
