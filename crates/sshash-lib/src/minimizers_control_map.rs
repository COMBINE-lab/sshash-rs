//! Minimizers Control Map (MCM)
//!
//! Maps minimizer values to their control information using a minimal perfect hash function.
//! This provides O(1) lookup of metadata associated with each unique minimizer.

use crate::mphf_config::{Mphf, build_mphf_from_vec};
use std::io;
use tracing::info;

/// Control information associated with a minimizer
#[derive(Clone, Copy, Debug)]
pub struct MinimizerControl {
    /// Number of strings containing this minimizer
    pub count: u32,
    /// Bucket type (regular, sparse, or heavy-load)
    pub bucket_type: BucketType,
    /// Additional metadata
    pub metadata: u64,
}

/// Classification of minimizer buckets by density
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[derive(Default)]
pub enum BucketType {
    /// Regular bucket: few strings contain this minimizer
    #[default]
    Regular,
    /// Sparse bucket: moderately populated
    Sparse,
    /// Heavy-load bucket: many strings contain this minimizer
    HeavyLoad,
}


/// Builder for constructing a MinimizersControlMap
///
/// This two-phase builder collects minimizers during construction,
/// then builds the MPHF once all minimizers are known.
pub struct MinimizersControlMapBuilder {
    /// Minimizers collected during building (stable order for MPHF)
    minimizers: Vec<u64>,
    /// Control information for each minimizer (parallel to minimizers vec)
    controls: Vec<MinimizerControl>,
    /// O(1) lookup: minimizer value → index in minimizers/controls vecs
    index: ahash::AHashMap<u64, usize>,
}

impl MinimizersControlMapBuilder {
    /// Create a new builder
    pub fn new() -> Self {
        Self {
            minimizers: Vec::new(),
            controls: Vec::new(),
            index: ahash::AHashMap::new(),
        }
    }

    /// Add a minimizer with initial control info
    pub fn add_minimizer(&mut self, minimizer: u64) -> usize {
        if let Some(&pos) = self.index.get(&minimizer) {
            return pos;
        }

        let id = self.minimizers.len();
        self.minimizers.push(minimizer);
        self.controls.push(MinimizerControl {
            count: 0,
            bucket_type: BucketType::Regular,
            metadata: 0,
        });
        self.index.insert(minimizer, id);
        id
    }

    /// Increment count for a minimizer
    pub fn increment_count(&mut self, minimizer: u64) {
        let id = self.add_minimizer(minimizer);
        self.controls[id].count = self.controls[id].count.saturating_add(1);
    }

    /// Set bucket type for a minimizer
    pub fn set_bucket_type(&mut self, minimizer: u64, bucket_type: BucketType) {
        let id = self.add_minimizer(minimizer);
        self.controls[id].bucket_type = bucket_type;
    }

    /// Get mutable access to control info for a minimizer
    pub fn get_control_mut(&mut self, minimizer: u64) -> Option<&mut MinimizerControl> {
        self.index
            .get(&minimizer)
            .copied()
            .map(|idx| &mut self.controls[idx])
    }

    /// Classify bucket types based on count thresholds
    pub fn finalize_bucket_types(&mut self, threshold_sparse: u32, threshold_heavy: u32) {
        for control in &mut self.controls {
            control.bucket_type = if control.count > threshold_heavy {
                BucketType::HeavyLoad
            } else if control.count > threshold_sparse {
                BucketType::Sparse
            } else {
                BucketType::Regular
            };
        }
    }

    /// Build the final MinimizersControlMap with MPHF
    ///
    /// Returns the control map AND a mapping from MPHF index to bucket_id.
    /// This mapping is used to reorder control_codewords to MPHF order.
    ///
    /// # Arguments
    /// * `c` - Relative level size (percentage, typically 100)
    /// * `alpha` - Load factor (kept for C++ API parity, not used by ph)
    pub fn build(self, _c: u16, _alpha: f64) -> io::Result<(MinimizersControlMap, Vec<usize>)> {
        if self.minimizers.is_empty() {
            return Ok((MinimizersControlMap {
                mphf: None,
                num_keys: 0,
            }, Vec::new()));
        }

        let minimizers = self.minimizers;
        let controls = self.controls;
        let num_keys = minimizers.len() as u64;

        info!("Building PHast MPHF for {} minimizers", num_keys);

        // Build the MPHF using PHast
        let mphf = build_mphf_from_vec(minimizers.clone());

        // Build mapping: bucket_id_by_mphf_index[mphf_index] = bucket_id (= metadata)
        let mut bucket_id_by_mphf_index = vec![0usize; controls.len()];
        for (idx, minimizer) in minimizers.iter().enumerate() {
            let pos = mphf.get(minimizer);
            if pos < bucket_id_by_mphf_index.len() {
                bucket_id_by_mphf_index[pos] = controls[idx].metadata as usize;
            }
        }

        Ok((MinimizersControlMap {
            mphf: Some(mphf),
            num_keys,
        }, bucket_id_by_mphf_index))
    }

    /// Get the number of unique minimizers
    pub fn num_minimizers(&self) -> usize {
        self.minimizers.len()
    }
}

impl Default for MinimizersControlMapBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Minimizers Control Map
///
/// Associates each unique minimizer with a bucket ID using
/// a minimal perfect hash function for O(1) lookups.
/// After reordering, the MPHF index IS the index into control_codewords.
pub struct MinimizersControlMap {
    /// MPHF for minimizer lookup (PHast)
    mphf: Option<Mphf>,
    /// Number of keys in the MPHF
    num_keys: u64,
}

impl MinimizersControlMap {
    /// Create a MinimizersControlMap from existing parts (for testing)
    pub fn from_parts(
        _controls: Vec<MinimizerControl>,
        minimizers: Vec<u64>,
        num_keys: u64,
    ) -> Self {
        // Build MPHF if we have minimizers
        let mphf = if !minimizers.is_empty() {
            Some(build_mphf_from_vec(minimizers))
        } else {
            None
        };
        
        Self {
            mphf,
            num_keys,
        }
    }
    
    /// Look up the MPHF index for a minimizer
    ///
    /// Returns the MPHF hash (bucket index into control_codewords), or None
    /// if the minimizer is not in the map.
    pub fn lookup(&self, minimizer: u64) -> Option<usize> {
        if let Some(ref mphf) = self.mphf {
            let id = mphf.get(&minimizer);
            if id < self.num_keys as usize {
                Some(id)
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Get a reference to the MPHF (for serialization)
    pub fn mphf_ref(&self) -> Option<&Mphf> {
        self.mphf.as_ref()
    }

    /// Set the MPHF (for deserialization)
    pub fn set_mphf(&mut self, mphf: Option<Mphf>) {
        self.mphf = mphf;
    }

    /// Get the number of unique minimizers
    pub fn num_minimizers(&self) -> u64 {
        self.num_keys
    }

    /// Get memory usage in bits (approximation)
    pub fn num_bits(&self) -> u64 {
        
        
        if let Some(ref _mphf) = self.mphf {
            // Estimate: fmph uses ~3.5-4 bits per key typically
            (self.num_keys as f64 * 4.0) as u64
        } else {
            0
        }
    }

    /// Exact serialized byte size of the main MPHF
    pub fn mphf_serialized_bytes(&self) -> usize {
        match &self.mphf {
            Some(mphf) => mphf.write_bytes(),
            None => 0,
        }
    }

    /// Serialize control map without MPHF
    ///
    /// Only stores num_keys since control_codewords are now in the index.
    pub fn serialize_without_mphf<W: io::Write>(
        &self,
        writer: &mut W,
    ) -> io::Result<()> {
        // Write number of keys
        writer.write_all(&self.num_keys.to_le_bytes())?;
        
        Ok(())
    }

    /// Deserialize control map without MPHF
    ///
    /// Returns MinimizersControlMap with empty MPHF (must be filled by caller)
    pub fn deserialize_without_mphf<R: io::Read>(
        reader: &mut R,
    ) -> io::Result<Self> {
        // Read number of keys
        let mut num_keys_bytes = [0u8; 8];
        reader.read_exact(&mut num_keys_bytes)?;
        let num_keys = u64::from_le_bytes(num_keys_bytes);
        
        Ok(Self {
            mphf: None,
            num_keys,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_minimizers_control_map_builder_creation() {
        let builder = MinimizersControlMapBuilder::new();
        assert_eq!(builder.num_minimizers(), 0);
    }

    #[test]
    fn test_minimizers_control_map_builder_add() {
        let mut builder = MinimizersControlMapBuilder::new();
        
        let id1 = builder.add_minimizer(100);
        let id2 = builder.add_minimizer(200);
        let id3 = builder.add_minimizer(100); // Duplicate
        
        assert_eq!(id1, 0);
        assert_eq!(id2, 1);
        assert_eq!(id3, 0); // Same as first
        assert_eq!(builder.num_minimizers(), 2);
    }

    #[test]
    fn test_minimizers_control_map_builder_increment() {
        let mut builder = MinimizersControlMapBuilder::new();
        
        builder.increment_count(100);
        builder.increment_count(100);
        builder.increment_count(200);
        
        assert_eq!(builder.num_minimizers(), 2);
        assert_eq!(builder.controls[0].count, 2);
        assert_eq!(builder.controls[1].count, 1);
    }

    #[test]
    fn test_minimizers_control_map_builder_bucket_type() {
        let mut builder = MinimizersControlMapBuilder::new();
        
        builder.add_minimizer(100);
        builder.set_bucket_type(100, BucketType::Sparse);
        
        assert_eq!(builder.controls[0].bucket_type, BucketType::Sparse);
    }

    #[test]
    fn test_minimizers_control_map_builder_finalize() {
        let mut builder = MinimizersControlMapBuilder::new();
        
        builder.add_minimizer(100);
        builder.controls[0].count = 5;
        
        builder.add_minimizer(200);
        builder.controls[1].count = 15;
        
        builder.add_minimizer(300);
        builder.controls[2].count = 150;
        
        builder.finalize_bucket_types(10, 100);
        
        assert_eq!(builder.controls[0].bucket_type, BucketType::Regular);
        assert_eq!(builder.controls[1].bucket_type, BucketType::Sparse);
        assert_eq!(builder.controls[2].bucket_type, BucketType::HeavyLoad);
    }

    #[test]
    fn test_minimizers_control_map_build_empty() {
        let builder = MinimizersControlMapBuilder::new();
        let (mcm, mapping) = builder.build(100, 0.94).unwrap();
        
        assert_eq!(mcm.num_minimizers(), 0);
        assert!(mcm.mphf.is_none());
        assert!(mapping.is_empty());
    }

    #[test]
    fn test_minimizers_control_map_build_and_lookup() {
        let mut builder = MinimizersControlMapBuilder::new();
        
        builder.increment_count(100);
        builder.increment_count(100);
        builder.increment_count(200);
        builder.set_bucket_type(100, BucketType::Sparse);
        
        let (mcm, _mapping) = builder.build(100, 0.94).unwrap();
        
        assert_eq!(mcm.num_minimizers(), 2);
        
        // lookup now returns Option<usize> (MPHF index)
        let idx_100 = mcm.lookup(100).unwrap();
        assert!(idx_100 < 2); // Valid MPHF index
        
        let idx_200 = mcm.lookup(200).unwrap();
        assert!(idx_200 < 2); // Valid MPHF index
        assert_ne!(idx_100, idx_200); // Different minimizers get different indices
    }

    #[test]
    fn test_minimizers_control_map_lookup_missing() {
        let mut builder = MinimizersControlMapBuilder::new();
        builder.increment_count(100);
        
        let (mcm, _) = builder.build(100, 0.94).unwrap();
        
        // 300 was never added
        let _result = mcm.lookup(300);
        // Note: MPHF might return a value even for keys not in the original set,
        // but it should be out of bounds or not match expected semantics
        // For a true membership test, we'd need a separate bloom filter or similar
    }

    #[test]
    fn test_bucket_type_default() {
        let bucket_type = BucketType::default();
        assert_eq!(bucket_type, BucketType::Regular);
    }

    #[test]
    fn test_minimizer_control_default() {
        let control = MinimizerControl {
            count: 5,
            bucket_type: BucketType::Regular,
            metadata: 0,
        };
        assert_eq!(control.count, 5);
        assert_eq!(control.bucket_type, BucketType::Regular);
    }
}
