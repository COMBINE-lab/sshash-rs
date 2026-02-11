//! Spectrum-Preserving String Set (SPSS)
//!
//! This is the core data structure that stores the strings in a compressed form.
//! Each string is stored as a sequence of 2-bit encoded DNA bases.

use crate::offsets::{EliasFanoOffsets, OffsetsVector};
use crate::kmer::{Kmer, KmerBits};
use std::io::{self, Read, Write};

/// The spectrum-preserving string set
///
/// Stores all strings in a bit-packed format with offsets to each string.
/// This allows for both memory-efficient storage and efficient access patterns.
///
/// Offsets are stored using Elias-Fano encoding (via cseq) for compact
/// representation and O(1) `locate()` via successor queries with Cursor.
pub struct SpectrumPreservingStringSet {
    /// Bit-packed string data (2 bits per base)
    strings: Vec<u8>,
    /// Offsets into the strings vector (in terms of bases, not bytes),
    /// encoded using Elias-Fano for ~2 bits/element vs 64 bits/element.
    offsets: EliasFanoOffsets,
    /// K-mer size of the strings
    k: usize,
    /// M (minimizer) size
    m: usize,
}

impl SpectrumPreservingStringSet {
    /// Create a new empty SPSS
    pub fn new(k: usize, m: usize) -> Self {
        Self {
            strings: Vec::new(),
            offsets: EliasFanoOffsets::from_vec(&[0]),
            k,
            m,
        }
    }
    
    /// Create a new SPSS from existing strings and offsets
    ///
    /// Converts the `OffsetsVector` to Elias-Fano encoding for compact storage.
    ///
    /// # Arguments
    /// * `strings` - Encoded string data (2-bit packed)
    /// * `offsets` - Offset vector for string boundaries (will be converted to EF)
    /// * `k` - K-mer size
    /// * `m` - Minimizer size
    pub fn from_parts(strings: Vec<u8>, offsets: OffsetsVector, k: usize, m: usize) -> Self {
        Self {
            strings,
            offsets: EliasFanoOffsets::from_offsets_vector(offsets),
            k,
            m,
        }
    }

    /// Get the string offsets (begin, end) for a string ID
    pub fn string_offsets(&self, string_id: u32) -> (u64, u64) {
        let id = string_id as usize;
        let begin = self.offsets.access(id);
        let end = self.offsets.access(id + 1);
        (begin, end)
    }

    /// Get the number of strings stored
    pub fn num_strings(&self) -> u64 {
        if !self.offsets.is_empty() {
            (self.offsets.len() - 1) as u64
        } else {
            0
        }
    }
    
    /// Get the starting offset (in bases) of a string
    pub fn string_offset(&self, string_id: u64) -> u64 {
        self.offsets.access(string_id as usize)
    }

    /// Get the k-mer size
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get the minimizer size
    pub fn m(&self) -> usize {
        self.m
    }

    /// Get the total number of bases stored
    pub fn total_bases(&self) -> u64 {
        if !self.offsets.is_empty() {
            self.offsets.access(self.offsets.len() - 1)
        } else {
            0
        }
    }

    /// Locate which string contains a given absolute position.
    /// Returns `(string_id, string_begin)` or None if out of bounds.
    #[inline]
    pub fn locate(&self, absolute_pos: u64) -> Option<(u64, u64)> {
        self.offsets.locate(absolute_pos)
    }

    /// Locate which string contains a given absolute position, returning
    /// `(string_id, string_begin, string_end)` in a single EF traversal.
    /// This is more efficient than calling `locate()` + `string_offsets()`.
    #[inline]
    pub fn locate_with_end(&self, absolute_pos: u64) -> Option<(u64, u64, u64)> {
        self.offsets.locate_with_end(absolute_pos)
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        (self.strings.len() as u64) * 8 + self.offsets.num_bits()
    }

    /// Get the byte size of the packed strings data
    pub fn strings_bytes(&self) -> usize {
        self.strings.len()
    }

    /// Get the byte size of the offsets vector
    pub fn offsets_bytes(&self) -> usize {
        self.offsets.num_bytes() as usize
    }
    
    /// Get the length of a specific string in bases
    pub fn string_length(&self, string_id: u64) -> usize {
        let (begin, end) = self.string_offsets(string_id as u32);
        (end - begin) as usize
    }
    
    /// Decode a k-mer from a specific position in a string
    ///
    /// Uses word-level loads from the packed buffer for efficiency.
    #[inline]
    pub fn decode_kmer<const K: usize>(&self, string_id: u64, kmer_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        let (begin, _end) = self.string_offsets(string_id as u32);
        let start_base = (begin as usize) + kmer_pos;
        
        let byte_offset = start_base / 4;
        let bit_shift = (start_base % 4) * 2;
        let needed_bits = K * 2;
        
        // Read enough bytes to cover K bases + potential misalignment
        // K bases need ceil((K*2 + bit_shift) / 8) bytes
        let needed_bytes = (needed_bits + bit_shift).div_ceil(8);
        
        if needed_bytes <= 8 {
            // Single u64 read suffices (K <= ~28-31 bases, common case)
            let mut buf = [0u8; 8];
            let avail = self.strings.len().saturating_sub(byte_offset).min(8);
            buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
            let raw = u64::from_le_bytes(buf);
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 64 { u64::MAX } else { (1u64 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
        } else {
            // Need u128 for K > 28 or large bit_shift
            let mut buf = [0u8; 16];
            let avail = self.strings.len().saturating_sub(byte_offset).min(16);
            buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
            let raw = u128::from_le_bytes(buf);
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
        }
    }

    /// Decode a k-mer at an absolute base position in the concatenated strings.
    /// Avoids the need for string_id (no binary search needed).
    /// This matches the C++ `util::read_kmer_at` approach with decoded_offsets.
    #[inline]
    pub fn decode_kmer_at<const K: usize>(&self, absolute_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        let byte_offset = absolute_pos / 4;
        let bit_shift = (absolute_pos % 4) * 2;
        let needed_bits = K * 2;
        let needed_bytes = (needed_bits + bit_shift).div_ceil(8);
        
        if needed_bytes <= 8 {
            // Fast path: unaligned u64 load when we have enough bytes.
            // Avoids bounds-checked copy_from_slice + zero-fill overhead.
            // Note: ptr::read_unaligned produces native-endian u64, which on
            // little-endian matches u64::from_le_bytes (the 2-bit encoding
            // assumes LE byte order).
            let raw = if byte_offset + 8 <= self.strings.len() {
                unsafe {
                    std::ptr::read_unaligned(
                        self.strings.as_ptr().add(byte_offset) as *const u64
                    )
                }
            } else {
                // Near end of strings — rare fallback
                let mut buf = [0u8; 8];
                let avail = self.strings.len() - byte_offset;
                buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                u64::from_le_bytes(buf)
            };
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 64 { u64::MAX } else { (1u64 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
        } else {
            let raw = if byte_offset + 16 <= self.strings.len() {
                unsafe {
                    std::ptr::read_unaligned(
                        self.strings.as_ptr().add(byte_offset) as *const u128
                    )
                }
            } else {
                let mut buf = [0u8; 16];
                let avail = self.strings.len() - byte_offset;
                buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                u128::from_le_bytes(buf)
            };
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
        }
    }

    /// Serialize the SPSS to a writer using a custom binary format.
    ///
    /// Format:
    /// - k: u64 (LE)
    /// - m: u64 (LE)
    /// - strings_len: u64 (LE)
    /// - strings: [u8; strings_len]
    /// - offsets: cseq Elias-Fano binary format
    pub fn serialize_to(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_all(&(self.k as u64).to_le_bytes())?;
        writer.write_all(&(self.m as u64).to_le_bytes())?;
        writer.write_all(&(self.strings.len() as u64).to_le_bytes())?;
        writer.write_all(&self.strings)?;
        self.offsets.write_to(writer)?;
        Ok(())
    }

    /// Deserialize an SPSS from a reader.
    pub fn deserialize_from(reader: &mut dyn Read) -> io::Result<Self> {
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let k = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let m = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let strings_len = u64::from_le_bytes(buf8) as usize;
        let mut strings = vec![0u8; strings_len];
        reader.read_exact(&mut strings)?;
        let offsets = EliasFanoOffsets::read_from(reader)?;
        Ok(Self { strings, offsets, k, m })
    }
}

impl std::fmt::Debug for SpectrumPreservingStringSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SpectrumPreservingStringSet")
            .field("k", &self.k)
            .field("m", &self.m)
            .field("num_strings", &self.num_strings())
            .field("total_bases", &self.total_bases())
            .field("num_bits", &self.num_bits())
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::offsets::OffsetsVector;

    /// Helper: build a test SPSS from DNA strings using the same Encoder
    /// logic (2-bit packing with offsets).
    fn build_test_spss(k: usize, m: usize, strings: &[&str]) -> SpectrumPreservingStringSet {
        let mut packed = Vec::new();
        let mut offsets = OffsetsVector::new();
        let mut total_bases: u64 = 0;

        for s in strings {
            for &b in s.as_bytes() {
                let bits = match b {
                    b'A' | b'a' => 0u8,
                    b'C' | b'c' => 1u8,
                    b'G' | b'g' => 3u8,
                    b'T' | b't' => 2u8,
                    _ => panic!("invalid base"),
                };
                let byte_idx = (total_bases as usize) / 4;
                let bit_off = ((total_bases as usize) % 4) * 2;
                if byte_idx >= packed.len() {
                    packed.push(0u8);
                }
                packed[byte_idx] |= bits << bit_off;
                total_bases += 1;
            }
            offsets.push(total_bases);
        }

        SpectrumPreservingStringSet::from_parts(packed, offsets, k, m)
    }

    #[test]
    fn test_spss_creation() {
        let spss = SpectrumPreservingStringSet::new(31, 13);
        assert_eq!(spss.k(), 31);
        assert_eq!(spss.m(), 13);
        assert_eq!(spss.num_strings(), 0);
    }

    #[test]
    fn test_spss_two_strings() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGCA", // 32 chars
        ]);
        assert_eq!(spss.num_strings(), 2);
    }

    #[test]
    fn test_spss_string_offsets() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGC", // 31 chars
        ]);

        let (begin1, end1) = spss.string_offsets(0);
        let (begin2, end2) = spss.string_offsets(1);

        assert_eq!(begin1, 0);
        assert_eq!(end1 - begin1, 31);
        assert_eq!(begin2, 31);
        assert_eq!(end2 - begin2, 31);
    }

    #[test]
    fn test_spss_total_bases() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",  // 31 chars
            "TGCATGCATGCATGCATGCATGCATGCATGC", // 31 chars
        ]);
        assert_eq!(spss.total_bases(), 62);
    }

    #[test]
    fn test_spss_debug() {
        let spss = build_test_spss(31, 13, &[
            "ACGTACGTACGTACGTACGTACGTACGTACG",
        ]);
        let debug_str = format!("{:?}", spss);
        assert!(debug_str.contains("SpectrumPreservingStringSet"));
        assert!(debug_str.contains("k: 31"));
    }
}
