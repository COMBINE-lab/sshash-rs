//! Serialization and deserialization support for Dictionary
//!
//! This module provides efficient zero-copy serialization using epserde for sux-rs types
//! (BitFieldVec, etc.) combined with native pthash serialization for MPHF functions.
//!
//! # File Format
//!
//! The serialization uses a two-file approach:
//!
//! **Main Index File** (`index.ssi`):
//! - DictionarySerializationHeader (magic, version, k, m, canonical, num_mphf_partitions)
//! - SpectrumPreservingStringSet (epserde format)
//! - SparseAndSkewIndex (epserde format, excluding MPHF)
//!
//! **MPHF Container File** (`index.ssi.mphf`):
//! ```text
//! MphfContainerHeader
//!   ├─ magic: "SSHIMH01"
//!   ├─ version_major: u32
//!   ├─ version_minor: u32
//!   └─ num_partitions: u32
//! Offset Table ([num_partitions] entries):
//!   ├─ MphfPartitionEntry 0
//!   │  ├─ partition_id: u32
//!   │  ├─ byte_offset: u64
//!   │  └─ byte_size: u64
//!   ├─ MphfPartitionEntry 1
//!   └─ ...
//! Data Section (variable length):
//!   ├─ MPHF partition 0 (raw fmph::GOFunction serialization)
//!   ├─ MPHF partition 1
//!   └─ ...
//! ```
//!
//! # Benefits of Single MPHF Container
//!
//! - **Scalability**: Works with 1 or 1000 partitions equally well (single file)
//! - **Random access**: Offset table enables seeking to any partition
//! - **Memory mappable**: Entire container can be mmap'd
//! - **Efficient**: No per-file overhead, compact layout
//! - **Clean separation**: MPHF container is independent binary format
//!
//! # Zero-Copy Deserialization
//!
//! When deserializing, sux-rs types are handled by epserde:
//! - `BitFieldVec<Vec<usize>>` deserializes as `BitFieldVec<&[usize]>` (ε-copy)
//! - The deserialized Dictionary can be memory-mapped for instant loading

use std::io::{self, Read, Write, Seek, SeekFrom};
use std::path::{Path, PathBuf};

/// Magic bytes for the SSHash index format
const MAGIC: &[u8; 8] = b"SSHIDX01";

/// Magic bytes for the SSHash MPHF container format
const MPHF_MAGIC: &[u8; 8] = b"SSHIMH01";

/// File format version: (major, minor)
/// Increment major on breaking changes, minor on compatible changes
const FORMAT_VERSION: (u32, u32) = (2, 0);
const MPHF_FORMAT_VERSION: (u32, u32) = (1, 0);

/// Header for the serialized Dictionary
#[derive(Clone, Debug)]
pub struct DictionarySerializationHeader {
    /// Magic number for format identification ("SSHIDX01")
    pub magic: [u8; 8],
    /// Format version (major, minor)
    pub version_major: u32,
    /// Format version minor number
    pub version_minor: u32,
    /// K-mer size
    pub k: usize,
    /// Minimizer size
    pub m: usize,
    /// Whether canonical mode is enabled
    pub canonical: bool,
    /// Number of MPHF partitions (for heavy buckets)
    pub num_mphf_partitions: u32,
}

impl DictionarySerializationHeader {
    /// Create a new header
    pub fn new(k: usize, m: usize, canonical: bool, num_mphf_partitions: u32) -> Self {
        Self {
            magic: *MAGIC,
            version_major: FORMAT_VERSION.0,
            version_minor: FORMAT_VERSION.1,
            k,
            m,
            canonical,
            num_mphf_partitions,
        }
    }

    /// Write header to a writer
    pub fn write(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_all(&self.magic)?;
        writer.write_all(&self.version_major.to_le_bytes())?;
        writer.write_all(&self.version_minor.to_le_bytes())?;
        writer.write_all(&(self.k as u64).to_le_bytes())?;
        writer.write_all(&(self.m as u64).to_le_bytes())?;
        writer.write_all(&[self.canonical as u8])?;
        writer.write_all(&self.num_mphf_partitions.to_le_bytes())?;
        Ok(())
    }

    /// Read header from a reader
    pub fn read(reader: &mut dyn Read) -> io::Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;

        if &magic != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid magic number for SSHash index file",
            ));
        }

        let mut version_major_bytes = [0u8; 4];
        let mut version_minor_bytes = [0u8; 4];
        let mut k_bytes = [0u8; 8];
        let mut m_bytes = [0u8; 8];
        let mut canonical_bytes = [0u8; 1];
        let mut num_partitions_bytes = [0u8; 4];

        reader.read_exact(&mut version_major_bytes)?;
        reader.read_exact(&mut version_minor_bytes)?;
        reader.read_exact(&mut k_bytes)?;
        reader.read_exact(&mut m_bytes)?;
        reader.read_exact(&mut canonical_bytes)?;
        reader.read_exact(&mut num_partitions_bytes)?;

        let version_major = u32::from_le_bytes(version_major_bytes);
        let version_minor = u32::from_le_bytes(version_minor_bytes);

        if version_major != FORMAT_VERSION.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Incompatible format version: {}.{}, expected {}.{}",
                    version_major, version_minor, FORMAT_VERSION.0, FORMAT_VERSION.1
                ),
            ));
        }

        Ok(Self {
            magic,
            version_major,
            version_minor,
            k: u64::from_le_bytes(k_bytes) as usize,
            m: u64::from_le_bytes(m_bytes) as usize,
            canonical: canonical_bytes[0] != 0,
            num_mphf_partitions: u32::from_le_bytes(num_partitions_bytes),
        })
    }
}

/// Entry in the MPHF container offset table
#[derive(Clone, Copy, Debug)]
pub struct MphfPartitionEntry {
    /// Partition ID
    pub partition_id: u32,
    /// Byte offset in the container file where this MPHF starts
    pub byte_offset: u64,
    /// Size in bytes of the serialized MPHF
    pub byte_size: u64,
}

impl MphfPartitionEntry {
    /// Write entry to a writer
    fn write(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_all(&self.partition_id.to_le_bytes())?;
        writer.write_all(&self.byte_offset.to_le_bytes())?;
        writer.write_all(&self.byte_size.to_le_bytes())?;
        Ok(())
    }

    /// Read entry from a reader
    fn read(reader: &mut dyn Read) -> io::Result<Self> {
        let mut id_bytes = [0u8; 4];
        let mut offset_bytes = [0u8; 8];
        let mut size_bytes = [0u8; 8];

        reader.read_exact(&mut id_bytes)?;
        reader.read_exact(&mut offset_bytes)?;
        reader.read_exact(&mut size_bytes)?;

        Ok(Self {
            partition_id: u32::from_le_bytes(id_bytes),
            byte_offset: u64::from_le_bytes(offset_bytes),
            byte_size: u64::from_le_bytes(size_bytes),
        })
    }
}

/// Header for the MPHF container file
///
/// The container format is:
/// ```text
/// MphfContainerHeader
/// offset_table: [MphfPartitionEntry; num_partitions]
/// data_section: (serialized MPHF data concatenated)
/// ```
#[derive(Clone, Debug)]
pub struct MphfContainerHeader {
    /// Magic number for format identification ("SSHIMH01")
    pub magic: [u8; 8],
    /// Format version (major, minor)
    pub version_major: u32,
    /// Format version minor number
    pub version_minor: u32,
    /// Number of MPHF partitions in this container
    pub num_partitions: u32,
}

impl MphfContainerHeader {
    /// Create a new MPHF container header
    pub fn new(num_partitions: u32) -> Self {
        Self {
            magic: *MPHF_MAGIC,
            version_major: MPHF_FORMAT_VERSION.0,
            version_minor: MPHF_FORMAT_VERSION.1,
            num_partitions,
        }
    }

    /// Write header to a writer
    pub fn write(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_all(&self.magic)?;
        writer.write_all(&self.version_major.to_le_bytes())?;
        writer.write_all(&self.version_minor.to_le_bytes())?;
        writer.write_all(&self.num_partitions.to_le_bytes())?;
        Ok(())
    }

    /// Read header from a reader
    pub fn read(reader: &mut dyn Read) -> io::Result<Self> {
        let mut magic = [0u8; 8];
        reader.read_exact(&mut magic)?;

        if &magic != MPHF_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid magic number for SSHash MPHF container file",
            ));
        }

        let mut version_major_bytes = [0u8; 4];
        let mut version_minor_bytes = [0u8; 4];
        let mut num_partitions_bytes = [0u8; 4];

        reader.read_exact(&mut version_major_bytes)?;
        reader.read_exact(&mut version_minor_bytes)?;
        reader.read_exact(&mut num_partitions_bytes)?;

        let version_major = u32::from_le_bytes(version_major_bytes);
        let version_minor = u32::from_le_bytes(version_minor_bytes);

        if version_major != MPHF_FORMAT_VERSION.0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Incompatible MPHF format version: {}.{}, expected {}.{}",
                    version_major, version_minor, MPHF_FORMAT_VERSION.0, MPHF_FORMAT_VERSION.1
                ),
            ));
        }

        Ok(Self {
            magic,
            version_major,
            version_minor,
            num_partitions: u32::from_le_bytes(num_partitions_bytes),
        })
    }
}

/// Build the main index file path from a base path
pub fn index_file_path<P: AsRef<Path>>(base: P) -> PathBuf {
    let mut path = base.as_ref().to_path_buf();
    let ext = path.extension().map(|e| e.to_string_lossy().to_string()).unwrap_or_default();
    if ext == "ssi" {
        // Already has .ssi extension
        path
    } else if ext.is_empty() {
        path.set_extension("ssi");
        path
    } else {
        path.set_extension(format!("{ext}.ssi"));
        path
    }
}

/// Build the MPHF container file path from a base path
pub fn mphf_container_path<P: AsRef<Path>>(base: P) -> PathBuf {
    let base_path = index_file_path(base);
    let mut container_path = base_path.clone();
    let filename = format!("{}.mphf", base_path.file_name().unwrap().to_string_lossy());
    container_path.pop();
    container_path.push(filename);
    container_path
}

/// Serialization errors
#[derive(Debug)]
pub enum SerializationError {
    /// I/O error during serialization
    Io(io::Error),
    /// Other serialization error
    Other(String),
}

impl From<io::Error> for SerializationError {
    fn from(err: io::Error) -> Self {
        SerializationError::Io(err)
    }
}

impl std::fmt::Display for SerializationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SerializationError::Io(e) => write!(f, "IO error: {}", e),
            SerializationError::Other(s) => write!(f, "{}", s),
        }
    }
}

impl std::error::Error for SerializationError {}

/// Result type for serialization operations
pub type SerializationResult<T> = Result<T, SerializationError>;

/// Helper functions for MPHF container operations
///
/// Write MPHFs to a container format
///
/// Creates a container file with:
/// - Header with num_partitions
/// - Offset table (partition_id, byte_offset, byte_size) for each partition
/// - Serialized MPHF data concatenated
///
/// Returns the offset table for reference
pub fn write_mphf_container<W: Write + Seek>(
    writer: &mut W,
    mphfs: &[Option<&crate::mphf_config::Mphf>],
) -> io::Result<Vec<MphfPartitionEntry>> {
    let num_partitions = mphfs.len() as u32;

    // Write header
    let header = MphfContainerHeader::new(num_partitions);
    header.write(writer)?;

    // Calculate and write offset table (with placeholders for now)
    let mut offset_table = Vec::new();
    let offset_table_start = writer.stream_position()?;

    // Write placeholder offset table
    for i in 0..num_partitions {
        let entry = MphfPartitionEntry {
            partition_id: i,
            byte_offset: 0, // Will be updated
            byte_size: 0,   // Will be updated
        };
        entry.write(writer)?;
    }

    let _data_start = writer.stream_position()?;

    // Serialize MPHFs and track their positions
    for (partition_id, mphf_opt) in mphfs.iter().enumerate() {
        let byte_offset = writer.stream_position()?;

        if let Some(mphf) = mphf_opt {
            // Serialize the MPHF to a temporary buffer to get the size
            let mut mphf_buffer = Vec::new();
            mphf.write(&mut mphf_buffer)?;
            let byte_size = mphf_buffer.len() as u64;

            // Write the serialized MPHF
            writer.write_all(&mphf_buffer)?;

            // Record the entry
            offset_table.push(MphfPartitionEntry {
                partition_id: partition_id as u32,
                byte_offset,
                byte_size,
            });
        } else {
            // Empty partition
            offset_table.push(MphfPartitionEntry {
                partition_id: partition_id as u32,
                byte_offset,
                byte_size: 0,
            });
        }
    }

    // Go back and write the actual offset table
    writer.seek(SeekFrom::Start(offset_table_start))?;
    for entry in &offset_table {
        entry.write(writer)?;
    }

    // Seek to end for any further writes
    writer.seek(SeekFrom::End(0))?;

    Ok(offset_table)
}

/// Read MPHFs from a container format
///
/// Returns a vector of Option<Mphf> indexed by partition ID
pub fn read_mphf_container<R: Read + Seek>(
    reader: &mut R,
) -> io::Result<Vec<Option<crate::mphf_config::Mphf>>> {
    // Read header
    let header = MphfContainerHeader::read(reader)?;

    // Read offset table
    let mut offset_table = Vec::with_capacity(header.num_partitions as usize);
    for _ in 0..header.num_partitions {
        offset_table.push(MphfPartitionEntry::read(reader)?);
    }

    // Read MPHFs
    let mut mphfs: Vec<Option<crate::mphf_config::Mphf>> = (0..header.num_partitions).map(|_| None).collect();

    for entry in offset_table {
        if entry.byte_size > 0 {
            reader.seek(SeekFrom::Start(entry.byte_offset))?;
            let mphf = crate::mphf_config::read_mphf(reader)?;
            mphfs[entry.partition_id as usize] = Some(mphf);
        }
    }

    Ok(mphfs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_header_roundtrip() {
        let header = DictionarySerializationHeader::new(31, 13, true, 2);

        let mut buffer = Vec::new();
        header.write(&mut buffer).unwrap();

        let header2 = DictionarySerializationHeader::read(&mut buffer.as_slice()).unwrap();

        assert_eq!(header.k, header2.k);
        assert_eq!(header.m, header2.m);
        assert_eq!(header.canonical, header2.canonical);
        assert_eq!(header.num_mphf_partitions, header2.num_mphf_partitions);
    }

    #[test]
    fn test_mphf_container_header_roundtrip() {
        let header = MphfContainerHeader::new(5);
        let mut buffer = Vec::new();
        header.write(&mut buffer).unwrap();

        let header2 = MphfContainerHeader::read(&mut buffer.as_slice()).unwrap();
        assert_eq!(header.num_partitions, header2.num_partitions);
    }

    #[test]
    fn test_mphf_partition_entry_roundtrip() {
        let entry = MphfPartitionEntry {
            partition_id: 3,
            byte_offset: 1024,
            byte_size: 512,
        };

        let mut buffer = Vec::new();
        entry.write(&mut buffer).unwrap();

        let entry2 = MphfPartitionEntry::read(&mut buffer.as_slice()).unwrap();
        assert_eq!(entry.partition_id, entry2.partition_id);
        assert_eq!(entry.byte_offset, entry2.byte_offset);
        assert_eq!(entry.byte_size, entry2.byte_size);
    }

    #[test]
    fn test_file_path_construction() {
        let base = Path::new("/tmp/my_index");
        let index = index_file_path(base);
        assert!(index.to_string_lossy().ends_with("my_index.ssi"));

        let mphf = mphf_container_path(base);
        assert!(mphf.to_string_lossy().contains("my_index.ssi.mphf"));
        assert!(!mphf.to_string_lossy().contains(".mphf.0")); // Single file, no partition ID
    }
}
