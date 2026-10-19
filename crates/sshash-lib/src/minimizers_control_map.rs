//! Minimizers Control Map (MCM)
//!
//! Maps minimizer values to their control information using a minimal perfect hash function.
//! This provides O(1) lookup of metadata associated with each unique minimizer.

use crate::partitioned_mphf::PartitionedMphf;
use std::io;

/// Minimizers Control Map
///
/// Associates each unique minimizer with a bucket ID using
/// a minimal perfect hash function for O(1) lookups.
/// After reordering, the MPHF index IS the index into control_codewords.
pub struct MinimizersControlMap {
    /// MPHF for minimizer lookup (partitioned PHastR)
    mphf: Option<PartitionedMphf>,
    /// Number of keys in the MPHF
    num_keys: u64,
}

impl MinimizersControlMap {
    /// An empty map (no minimizers).
    pub fn empty() -> Self {
        Self { mphf: None, num_keys: 0 }
    }

    /// Create a MinimizersControlMap from a pre-built MPHF.
    ///
    /// The builder builds the MPHF of the bucket minimizers (see
    /// [`crate::MphfBuild`]) and wraps it here.
    pub fn from_mphf(mphf: PartitionedMphf, num_keys: u64) -> Self {
        Self {
            mphf: Some(mphf),
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
    pub fn mphf_ref(&self) -> Option<&PartitionedMphf> {
        self.mphf.as_ref()
    }

    /// Set the MPHF (for deserialization)
    pub fn set_mphf(&mut self, mphf: Option<PartitionedMphf>) {
        self.mphf = mphf;
    }

    /// Get the number of unique minimizers
    pub fn num_minimizers(&self) -> u64 {
        self.num_keys
    }

    /// Size of the MPHF in bits (its serialized size)
    pub fn num_bits(&self) -> u64 {
        self.mphf_serialized_bytes() as u64 * 8
    }

    /// Exact serialized byte size of the main MPHF
    pub fn mphf_serialized_bytes(&self) -> usize {
        match &self.mphf {
            Some(pmphf) => pmphf.write_bytes(),
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
    fn test_empty() {
        let mcm = MinimizersControlMap::empty();
        assert_eq!(mcm.num_minimizers(), 0);
        assert!(mcm.mphf_ref().is_none());
        assert_eq!(mcm.lookup(100), None);
        assert_eq!(mcm.num_bits(), 0);
    }

    #[test]
    fn test_from_mphf_lookup() {
        let minimizers = [100u64, 200, 300];
        let mcm = MinimizersControlMap::from_mphf(PartitionedMphf::build_single(&minimizers), 3);
        assert_eq!(mcm.num_minimizers(), 3);
        let mut ids: Vec<usize> = minimizers.iter().map(|&m| mcm.lookup(m).unwrap()).collect();
        ids.sort_unstable();
        assert_eq!(ids, [0, 1, 2]);
        // An absent minimizer still lands in range: callers verify the key
        assert!(mcm.lookup(400).unwrap() < 3);
    }

    #[test]
    fn test_serialize_without_mphf_roundtrip() {
        let mcm = MinimizersControlMap::from_mphf(PartitionedMphf::build_single(&[7u64, 8]), 2);
        let mut buf = Vec::new();
        mcm.serialize_without_mphf(&mut buf).unwrap();
        let back = MinimizersControlMap::deserialize_without_mphf(&mut buf.as_slice()).unwrap();
        assert_eq!(back.num_minimizers(), 2);
        assert!(back.mphf_ref().is_none());
    }
}
