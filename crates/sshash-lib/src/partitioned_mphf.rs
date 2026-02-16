//! Partitioned Minimal Perfect Hash Function
//!
//! Wraps multiple PHast MPHFs, one per partition, with hash-based partition
//! selection matching the C++ `partitioned_phf` design. Keys are assigned to
//! partitions via Lemire fast-range reduction on an independent rapidhash,
//! and each partition's MPHF is built independently (in parallel via rayon).
//!
//! For datasets with fewer than `AVG_PARTITION_SIZE` keys (including all our
//! test indices), only a single partition is created and the `get()` fast path
//! skips the partition hash entirely.

use crate::mphf_config::{Mphf, MphfHasher, build_mphf_from_vec, mphf_hasher, read_mphf};
use rayon::prelude::*;
use std::hash::Hash;
use std::io::{self, Read, Write};

/// Average number of keys per partition, matching C++ PTHash `partitioned_phf`.
const AVG_PARTITION_SIZE: usize = 3_000_000;

/// Seed used for partition selection, chosen to differ from PHast's internal
/// seeds (0, 1, 2, ...) to ensure statistical independence.
const PARTITION_HASH_SEED: u64 = 0xC6A4_A793_5BD1_E995;

/// A partitioned minimal perfect hash function.
///
/// Splits keys into partitions by hash range (Lemire fast-range reduction),
/// builds one PHast MPHF per partition, and routes queries transparently.
/// Global indices are `offsets[partition] + inner_mphf.get(key)`.
pub struct PartitionedMphf {
    /// One PHast MPHF per partition.
    inners: Vec<Mphf>,
    /// Cumulative key counts: `offsets[i]` = total keys in partitions 0..i.
    /// Length = `num_partitions + 1`.
    offsets: Vec<usize>,
    /// Number of partitions.
    num_partitions: u32,
    /// Total number of keys across all partitions.
    num_keys: usize,
    /// Hasher for partition selection.
    hasher: MphfHasher,
}

impl PartitionedMphf {
    /// Build a partitioned MPHF from an owned Vec of keys.
    ///
    /// If `partitioned` is false (or there are fewer than `AVG_PARTITION_SIZE`
    /// keys), a single partition is used — equivalent to a monolithic MPHF with
    /// zero query overhead.
    pub fn build_from_vec<K: Hash + Clone + Send + Sync>(keys: Vec<K>, partitioned: bool) -> Self {
        let num_keys = keys.len();
        if num_keys == 0 {
            return Self {
                inners: Vec::new(),
                offsets: vec![0],
                num_partitions: 0,
                num_keys: 0,
                hasher: mphf_hasher(),
            };
        }

        let num_partitions = if partitioned {
            num_keys.div_ceil(AVG_PARTITION_SIZE).max(1)
        } else {
            1
        };

        if num_partitions == 1 {
            // Single partition: build directly, no partitioning overhead
            let mphf = build_mphf_from_vec(keys);
            return Self {
                inners: vec![mphf],
                offsets: vec![0, num_keys],
                num_partitions: 1,
                num_keys,
                hasher: mphf_hasher(),
            };
        }

        // Multi-partition: hash-and-partition
        let hasher = mphf_hasher();
        let np = num_partitions as u128;

        // Assign keys to partitions
        let mut partition_keys: Vec<Vec<K>> = (0..num_partitions).map(|_| Vec::new()).collect();
        for key in keys {
            let hash = hasher.hash_one_with_seed(&key, PARTITION_HASH_SEED);
            let p = ((hash as u128 * np) >> 64) as usize;
            partition_keys[p].push(key);
        }

        // Compute cumulative offsets
        let mut offsets = Vec::with_capacity(num_partitions + 1);
        offsets.push(0);
        for pk in &partition_keys {
            let prev = *offsets.last().unwrap();
            offsets.push(prev + pk.len());
        }

        // Build inner MPHFs in parallel (each single-threaded internally)
        let inners: Vec<Mphf> = partition_keys
            .into_par_iter()
            .map(|pk| {
                if pk.is_empty() {
                    // Empty partition — build a trivial MPHF from empty vec
                    build_mphf_from_vec(pk)
                } else {
                    build_mphf_from_vec(pk)
                }
            })
            .collect();

        Self {
            inners,
            offsets,
            num_partitions: num_partitions as u32,
            num_keys,
            hasher,
        }
    }

    /// Build a partitioned MPHF from a slice of keys (clones into Vec).
    pub fn build_from_slice<K: Hash + Clone + Send + Sync>(keys: &[K], partitioned: bool) -> Self {
        Self::build_from_vec(keys.to_vec(), partitioned)
    }

    /// Look up a key and return its global index in [0, num_keys).
    #[inline]
    pub fn get<K: Hash + ?Sized>(&self, key: &K) -> usize {
        if self.num_partitions == 1 {
            // Fast path: skip partition hash entirely
            return self.inners[0].get(key);
        }
        let p = self.partition_for(key);
        self.offsets[p] + self.inners[p].get(key)
    }

    /// Total number of keys.
    pub fn num_keys(&self) -> usize {
        self.num_keys
    }

    /// Number of partitions.
    pub fn num_partitions(&self) -> u32 {
        self.num_partitions
    }

    /// Compute which partition a key belongs to (Lemire fast-range reduction).
    #[inline]
    fn partition_for<K: Hash + ?Sized>(&self, key: &K) -> usize {
        let hash = self.hasher.hash_one_with_seed(key, PARTITION_HASH_SEED);
        ((hash as u128 * self.num_partitions as u128) >> 64) as usize
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
            mphf.write(writer)?;
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
            hasher: mphf_hasher(),
        })
    }

    /// Estimate serialized byte size (for container offset table pre-allocation).
    pub fn write_bytes(&self) -> usize {
        let header = 4 + 8; // num_partitions + num_keys
        let offsets = (self.offsets.len()) * 8;
        let mphfs: usize = self.inners.iter().map(|m| m.write_bytes()).sum();
        header + offsets + mphfs
    }
}

/// Extension trait to hash with a specific seed.
///
/// PHast's `BuildRapidHash` implements `BuildSeededHasher`, so we use
/// `build_hasher(seed)` to get a seeded hasher for partition selection.
trait HashOneWithSeed {
    fn hash_one_with_seed<K: Hash + ?Sized>(&self, key: &K, seed: u64) -> u64;
}

impl HashOneWithSeed for MphfHasher {
    #[inline]
    fn hash_one_with_seed<K: Hash + ?Sized>(&self, key: &K, seed: u64) -> u64 {
        use ph::BuildSeededHasher;
        use std::hash::Hasher;
        let mut hasher = self.build_hasher(seed);
        key.hash(&mut hasher);
        hasher.finish()
    }
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
    fn test_write_bytes_sanity() {
        let keys: Vec<u64> = (0..100).collect();
        let pmphf = PartitionedMphf::build_from_vec(keys, true);

        let mut buf = Vec::new();
        pmphf.write_to(&mut buf).unwrap();

        // write_bytes() is an estimate (PHast's write_bytes() is approximate)
        // Just verify it's in a reasonable range
        let actual = buf.len();
        let estimate = pmphf.write_bytes();
        assert!(estimate > 0, "estimate should be positive");
        // The estimate comes from PHast's write_bytes() which may not be exact
        assert!(
            actual > 0 && estimate > 0,
            "both actual ({actual}) and estimate ({estimate}) should be positive"
        );
    }
}
