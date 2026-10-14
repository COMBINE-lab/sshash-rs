//! Partitioned Minimal Perfect Hash Function
//!
//! Wraps multiple PHastR MPHFs, one per partition, with hash-based partition
//! selection matching the C++ `partitioned_phf` design. Keys are assigned to
//! partitions via Lemire fast-range reduction on an independently seeded
//! rapidhash, and each partition's MPHF is built independently (in parallel via
//! rayon).
//!
//! For datasets with fewer than `AVG_PARTITION_SIZE` keys (including all our
//! test indices), only a single partition is created and the `get()` fast path
//! skips the partition hash entirely.

use crate::mphf_config::{
    KeyBytes, Mphf, MphfKey, build_mphf, build_mphf_streaming, build_mphfs, hash_key_bytes, mphf_get,
    mphf_serialized_bytes,
    read_mphf, write_mphf,
};
use rayon::prelude::*;
use std::borrow::Borrow;
use std::io::{self, Read, Write};

/// Average number of keys per partition, matching C++ PTHash `partitioned_phf`.
const AVG_PARTITION_SIZE: usize = 3_000_000;

/// Seed used for partition selection, chosen to differ from the hash seeds
/// PHastR derives from its own seed (0) to ensure statistical independence.
const PARTITION_HASH_SEED: u64 = 0xC6A4_A793_5BD1_E995;

/// A partitioned minimal perfect hash function.
///
/// Splits keys into partitions by hash range (Lemire fast-range reduction),
/// builds one PHastR MPHF per partition, and routes queries transparently.
/// Global indices are `offsets[partition] + inner_mphf.get(key)`.
pub struct PartitionedMphf {
    /// One PHastR MPHF per partition.
    inners: Vec<Mphf>,
    /// Cumulative key counts: `offsets[i]` = total keys in partitions 0..i.
    /// Length = `num_partitions + 1`.
    offsets: Vec<usize>,
    /// Number of partitions.
    num_partitions: u32,
    /// Total number of keys across all partitions.
    num_keys: usize,
}

impl PartitionedMphf {
    /// Build a partitioned MPHF from an owned Vec of keys.
    ///
    /// If `partitioned` is false (or there are fewer than `AVG_PARTITION_SIZE`
    /// keys), a single partition is used — equivalent to a monolithic MPHF with
    /// zero query overhead.
    pub fn build_from_vec<K: MphfKey>(keys: Vec<K>, partitioned: bool) -> Self {
        Self::build_from_slice(&keys, partitioned)
    }

    /// Build a single-partition MPHF by streaming distinct keys, read once
    /// (see [`build_mphf_streaming`]). `num_keys_hint` selects offline
    /// construction for large key sets; `tmp_dir` holds its temporary files.
    pub fn build_streaming<I>(keys: I, num_keys_hint: usize, tmp_dir: &std::path::Path) -> Self
    where
        I: Iterator<Item: MphfKey> + Send,
    {
        let mphf = build_mphf_streaming(keys, num_keys_hint, tmp_dir);
        let num_keys = mphf.len();
        if num_keys == 0 {
            return Self { inners: Vec::new(), offsets: vec![0], num_partitions: 0, num_keys: 0 };
        }
        Self { inners: vec![mphf], offsets: vec![0, num_keys], num_partitions: 1, num_keys }
    }

    /// Build a partitioned MPHF from a slice of keys.
    pub fn build_from_slice<K: MphfKey>(keys: &[K], partitioned: bool) -> Self {
        let num_keys = keys.len();
        if num_keys == 0 {
            return Self {
                inners: Vec::new(),
                offsets: vec![0],
                num_partitions: 0,
                num_keys: 0,
            };
        }

        let num_partitions = if partitioned {
            num_keys.div_ceil(AVG_PARTITION_SIZE).max(1)
        } else {
            1
        };

        if num_partitions == 1 {
            // Single partition: build directly, no partitioning overhead
            let bytes: Vec<K::Bytes> = keys.iter().map(|k| k.key_bytes()).collect();
            let mphf = build_mphf(&bytes);
            return Self {
                inners: vec![mphf],
                offsets: vec![0, num_keys],
                num_partitions: 1,
                num_keys,
            };
        }

        // Multi-partition: hash-and-partition.
        // Count partition sizes first (in parallel) so each partition is
        // allocated exactly, instead of carrying up to 2x growth slack.
        let np = num_partitions as u32;
        let sizes = keys
            .par_chunks(1 << 16)
            .fold(
                || vec![0usize; num_partitions],
                |mut counts, chunk| {
                    for &key in chunk {
                        counts[partition_of(key.key_bytes().borrow(), np)] += 1;
                    }
                    counts
                },
            )
            .reduce(
                || vec![0usize; num_partitions],
                |mut a, b| {
                    a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                    a
                },
            );

        // Assign keys to partitions (in input order)
        let mut partition_keys: Vec<Vec<K::Bytes>> =
            sizes.iter().map(|&s| Vec::with_capacity(s)).collect();
        for &key in keys {
            let bytes = key.key_bytes();
            partition_keys[partition_of(bytes.borrow(), np)].push(bytes);
        }

        // Compute cumulative offsets
        let mut offsets = Vec::with_capacity(num_partitions + 1);
        offsets.push(0);
        for pk in &partition_keys {
            let prev = *offsets.last().unwrap();
            offsets.push(prev + pk.len());
        }

        // Build inner MPHFs, current_num_threads() at a time (each
        // single-threaded internally), freeing each partition's keys as it
        // is built
        let inners = build_mphfs(partition_keys);

        Self {
            inners,
            offsets,
            num_partitions: num_partitions as u32,
            num_keys,
        }
    }

    /// Look up a key and return its global index in [0, num_keys).
    ///
    /// For keys NOT in the build set, the result is arbitrary: PHastR maps every
    /// key into range, so callers must verify the key at the returned index.
    /// An empty MPHF returns `num_keys` (0), i.e., out of range.
    #[inline]
    pub fn get<K: MphfKey>(&self, key: &K) -> usize {
        if self.num_partitions == 1 {
            // Fast path: skip partition hash entirely.
            return mphf_get(&self.inners[0], *key);
        }
        if self.num_partitions == 0 {
            return self.num_keys;
        }
        let bytes = key.key_bytes();
        let p = partition_of(bytes.borrow(), self.num_partitions);
        self.offsets[p] + self.inners[p].get(bytes.borrow())
    }

    /// Total number of keys.
    pub fn num_keys(&self) -> usize {
        self.num_keys
    }

    /// Number of partitions.
    pub fn num_partitions(&self) -> u32 {
        self.num_partitions
    }

    /// Serialize to a writer.
    pub fn write_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        // Header
        writer.write_all(&self.num_partitions.to_le_bytes())?;
        writer.write_all(&(self.num_keys as u64).to_le_bytes())?;

        // Offsets (num_partitions + 1 entries)
        for &off in &self.offsets {
            writer.write_all(&(off as u64).to_le_bytes())?;
        }

        // Inner MPHFs
        for mphf in &self.inners {
            write_mphf(mphf, writer)?;
        }

        Ok(())
    }

    /// Deserialize from a reader.
    pub fn read_from(reader: &mut dyn Read) -> io::Result<Self> {
        let mut buf4 = [0u8; 4];
        let mut buf8 = [0u8; 8];

        reader.read_exact(&mut buf4)?;
        let num_partitions = u32::from_le_bytes(buf4);

        reader.read_exact(&mut buf8)?;
        let num_keys = u64::from_le_bytes(buf8) as usize;

        let num_offsets = num_partitions as usize + 1;
        let mut offsets = Vec::with_capacity(num_offsets);
        for _ in 0..num_offsets {
            reader.read_exact(&mut buf8)?;
            offsets.push(u64::from_le_bytes(buf8) as usize);
        }

        let mut inners = Vec::with_capacity(num_partitions as usize);
        for _ in 0..num_partitions {
            inners.push(read_mphf(reader)?);
        }

        Ok(Self {
            inners,
            offsets,
            num_partitions,
            num_keys,
        })
    }

    /// Exact serialized byte size (as written by [`Self::write_to`]).
    pub fn write_bytes(&self) -> usize {
        let header = 4 + 8; // num_partitions + num_keys
        let offsets = (self.offsets.len()) * 8;
        let mphfs: usize = self.inners.iter().map(mphf_serialized_bytes).sum();
        header + offsets + mphfs
    }
}

/// Compute which partition a key belongs to (Lemire fast-range reduction on
/// a rapidhash seeded independently of the per-partition MPHFs).
#[inline]
fn partition_of(key: &KeyBytes, num_partitions: u32) -> usize {
    let hash = hash_key_bytes(key.as_bytes(), PARTITION_HASH_SEED);
    ((hash as u128 * num_partitions as u128) >> 64) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partition_count_math() {
        // < AVG_PARTITION_SIZE → 1 partition
        assert_eq!(1_000_000usize.div_ceil(AVG_PARTITION_SIZE).max(1), 1);
        // Exactly AVG_PARTITION_SIZE → 1 partition
        assert_eq!(AVG_PARTITION_SIZE.div_ceil(AVG_PARTITION_SIZE).max(1), 1);
        // Just over → 2 partitions
        assert_eq!((AVG_PARTITION_SIZE + 1).div_ceil(AVG_PARTITION_SIZE).max(1), 2);
        // 10M → 4 partitions
        assert_eq!(10_000_000usize.div_ceil(AVG_PARTITION_SIZE).max(1), 4);
    }

    #[test]
    fn test_single_partition_roundtrip() {
        let keys: Vec<u64> = (0..1000).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), true);

        assert_eq!(pmphf.num_partitions(), 1);
        assert_eq!(pmphf.num_keys(), 1000);

        // All keys should get unique indices in [0, 1000)
        let mut indices: Vec<usize> = keys.iter().map(|k| pmphf.get(k)).collect();
        indices.sort();
        indices.dedup();
        assert_eq!(indices.len(), 1000);
        assert!(indices.iter().all(|&i| i < 1000));
    }

    #[test]
    fn test_monolithic_flag() {
        let keys: Vec<u64> = (0..100).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), false);

        assert_eq!(pmphf.num_partitions(), 1);
        assert_eq!(pmphf.num_keys(), 100);

        let mut indices: Vec<usize> = keys.iter().map(|k| pmphf.get(k)).collect();
        indices.sort();
        indices.dedup();
        assert_eq!(indices.len(), 100);
    }

    #[test]
    fn test_serialization_roundtrip() {
        let keys: Vec<u64> = (0..500).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), true);

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();

        let pmphf2 = PartitionedMphf::read_from(&mut buf.as_slice()).unwrap();

        assert_eq!(pmphf.num_partitions(), pmphf2.num_partitions());
        assert_eq!(pmphf.num_keys(), pmphf2.num_keys());

        // Verify same results
        for key in &keys {
            assert_eq!(pmphf.get(key), pmphf2.get(key));
        }
    }

    #[test]
    fn test_empty() {
        let keys: Vec<u64> = Vec::new();
        let pmphf = PartitionedMphf::build_from_vec(keys, true);
        assert_eq!(pmphf.num_partitions(), 0);
        assert_eq!(pmphf.num_keys(), 0);
    }

    #[test]
    fn test_multi_partition_bijection_and_roundtrip() {
        // Just over two partitions' worth of scattered keys
        let n = 2 * AVG_PARTITION_SIZE + 12_345;
        let keys: Vec<u64> = (0..n as u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), true);
        assert_eq!(pmphf.num_partitions(), 3);
        assert_eq!(pmphf.num_keys(), n);

        let mut seen = vec![false; n];
        for k in &keys {
            let i = pmphf.get(k);
            assert!(!seen[i], "index {i} assigned twice");
            seen[i] = true;
        }

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();
        assert_eq!(buf.len(), pmphf.write_bytes());
        let pmphf2 = PartitionedMphf::read_from(&mut buf.as_slice()).unwrap();
        for k in keys.iter().step_by(97) {
            assert_eq!(pmphf.get(k), pmphf2.get(k));
        }
        // Absent keys still land in range
        for k in 0..1000u64 {
            assert!(pmphf.get(&(k | (1 << 63))) < n);
        }
    }

    #[test]
    fn test_u128_keys() {
        let keys: Vec<u128> = (0..5000u128).map(|i| (i << 70) | i).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys.clone(), false);
        let mut indices: Vec<usize> = keys.iter().map(|k| pmphf.get(k)).collect();
        indices.sort_unstable();
        indices.dedup();
        assert_eq!(indices.len(), keys.len());
    }

    #[test]
    fn test_build_independent_of_thread_count() {
        // Large enough that PHastR would split its sweep across threads
        let keys: Vec<u64> = (0..2_000_000u64).map(|i| i.wrapping_mul(0xD6E8_FEB8_6659_FD93)).collect();
        let bytes_with = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let pmphf = pool.install(|| PartitionedMphf::build_from_slice(&keys, false));
            let mut buf = Vec::new();
            pmphf.write_to(&mut buf).unwrap();
            buf
        };
        assert_eq!(bytes_with(1), bytes_with(8));
    }

    #[test]
    fn test_partitioned_build_independent_of_thread_count() {
        let n = 2 * AVG_PARTITION_SIZE + 999;
        let keys: Vec<u64> = (0..n as u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        let bytes_with = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let pmphf = pool.install(|| PartitionedMphf::build_from_slice(&keys, true));
            let mut buf = Vec::new();
            pmphf.write_to(&mut buf).unwrap();
            buf
        };
        assert_eq!(bytes_with(1), bytes_with(8));
    }

    fn streamed_bytes(keys: &[u64], ambient_threads: usize, hint: usize) -> Vec<u8> {
        let tmp = tempfile::tempdir().unwrap();
        let pool = rayon::ThreadPoolBuilder::new().num_threads(ambient_threads).build().unwrap();
        let pmphf = pool.install(|| PartitionedMphf::build_streaming(keys.iter().copied(), hint, tmp.path()));
        assert_eq!(pmphf.num_keys(), keys.len());
        let mut seen = vec![false; keys.len()];
        for k in keys {
            let i = pmphf.get(k);
            assert!(!seen[i]);
            seen[i] = true;
        }
        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();
        buf
    }

    #[test]
    fn test_streaming_build_independent_of_thread_count() {
        // Large enough that PHastR splits its sweep across its pool's threads
        let keys: Vec<u64> = (0..2_000_000u64).map(|i| i.wrapping_mul(0xD6E8_FEB8_6659_FD93)).collect();
        assert_eq!(streamed_bytes(&keys, 1, keys.len()), streamed_bytes(&keys, 16, keys.len()));
    }

    #[test]
    #[ignore = "builds over OFFLINE_MIN_KEYS keys (256 MiB buffers); run explicitly"]
    fn test_streaming_offline_matches_in_memory() {
        let n = crate::mphf_config::OFFLINE_MIN_KEYS;
        let keys: Vec<u64> = (0..n as u64).map(|i| i.wrapping_mul(0x9E37_79B9_7F4A_7C15)).collect();
        // Same keys, same pool size: offline mode must not change the function
        assert_eq!(streamed_bytes(&keys, 4, n), streamed_bytes(&keys, 4, 0));
    }

    #[test]
    fn test_write_bytes_sanity() {
        let keys: Vec<u64> = (0..100).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys, true);

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();

        assert_eq!(buf.len(), pmphf.write_bytes());
    }
}
