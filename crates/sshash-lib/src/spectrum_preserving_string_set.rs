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
/// Offsets are stored using Elias-Fano encoding (via sux-rs) for compact
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
    pub fn from_parts(mut strings: Vec<u8>, offsets: OffsetsVector, k: usize, m: usize) -> Self {
        // Pad strings with 8 zero bytes to allow safe unaligned u64 reads
        // near the end without out-of-bounds access.
        strings.resize(strings.len() + 8, 0);
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

    /// Get raw access to the packed string data (for buffered reading)
    #[inline]
    pub fn strings_data(&self) -> &[u8] {
        &self.strings
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

    /// Get memory usage in bits (excludes padding bytes)
    pub fn num_bits(&self) -> u64 {
        ((self.strings.len() - 8) as u64) * 8 + self.offsets.num_bits()
    }

    /// Get the byte size of the packed strings data (excludes padding)
    pub fn strings_bytes(&self) -> usize {
        self.strings.len() - 8
    }

    /// Access the raw 2-bit-packed strings buffer. Intended for crates that
    /// need to clone the SPSS contents into their own representation.
    #[inline]
    pub fn strings_raw(&self) -> &[u8] {
        &self.strings
    }

    /// Number of offset entries (`num_strings + 1`).
    #[inline]
    pub fn offsets_len(&self) -> usize {
        self.offsets.len()
    }

    /// Get the offset at index `i` (0..=num_strings). Offsets are in bases.
    #[inline]
    pub fn offset_at(&self, i: usize) -> u64 {
        self.offsets.access(i)
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
            // Need u128 for K > 28 or large bit_shift.
            // For K=63 (needed_bits=126), when bit_shift > 2 we need more
            // than 128 bits from the byte stream. Load 16 bytes as u128
            // plus an extra byte to cover the overflow.
            let mut buf = [0u8; 17];
            let avail = self.strings.len().saturating_sub(byte_offset).min(17);
            buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
            let raw = u128::from_le_bytes(buf[..16].try_into().unwrap());
            let shifted = if bit_shift > 0 {
                let extra = buf[16] as u128;
                (raw >> bit_shift) | (extra << (128 - bit_shift))
            } else {
                raw
            };
            let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
        }
    }

    /// Decode the m-mer (m <= 31) at an absolute base position.
    ///
    /// Used by the lookup path to verify that a bucket's first offset really
    /// holds the queried minimizer (see `Dictionary::lookup_with_minimizer`).
    /// Same branchless unaligned read as the K<=31 arm of
    /// [`Self::decode_kmer_at`]; the 8-byte tail padding makes it safe.
    #[inline(always)]
    pub fn decode_mmer_at(&self, absolute_pos: usize, m: usize) -> u64 {
        debug_assert!(m >= 1 && 2 * m < 64);
        let byte_offset = absolute_pos / 4;
        let bit_shift = (absolute_pos % 4) * 2;
        let (raw, extra) = unsafe {
            let base = self.strings.as_ptr().add(byte_offset);
            // SAFETY: as in decode_kmer_at — `strings` carries 8 bytes of
            // tail padding and a decodable m-mer starts before the logical end.
            (std::ptr::read_unaligned(base as *const u64), *base.add(8) as u64)
        };
        let result = (raw >> bit_shift) | ((extra << 1) << (63 - bit_shift));
        result & ((1u64 << (2 * m)) - 1)
    }

    /// Decode a k-mer at an absolute base position in the concatenated strings.
    #[inline(always)]
    pub fn decode_kmer_at<const K: usize>(&self, absolute_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        let byte_offset = absolute_pos / 4;
        let bit_shift = (absolute_pos % 4) * 2;
        let needed_bits = K * 2;

        if K <= 31 {
            // K<=31: need at most 62+6=68 bits, so a u64 at `byte_offset` plus the
            // one byte the top bits can spill into.
            //
            // The spill is folded in unconditionally. `bit_shift` is `(pos % 4) * 2`,
            // so it is one of {0, 2, 4, 6} and the spill only actually matters for
            // {4, 6} — i.e. a branch that is taken for half of all positions, on a
            // value derived from the read position, so it does not predict. Making
            // it branchless is worth more than the byte load it saves.
            //
            // `(extra << 1) << (63 - bit_shift)` equals the wanted
            // `extra << (64 - bit_shift)` for bit_shift in {2, 4, 6}, and is
            // identically zero for bit_shift == 0 (the `<< 1` clears the only bit
            // that would survive `<< 63`). That dodges the shift-by-64 UB that
            // makes the direct expression unusable without the branch.
            let (raw, extra) = unsafe {
                let base = self.strings.as_ptr().add(byte_offset);
                // SAFETY: `strings` carries 8 bytes of tail padding (see `new` and
                // the deserializer), and a decodable k-mer starts before the logical
                // end, so `byte_offset + 8` is always within the allocation. This is
                // the same maximum address the previous branching form could touch.
                (std::ptr::read_unaligned(base as *const u64), *base.add(8) as u64)
            };
            let result = (raw >> bit_shift) | ((extra << 1) << (63 - bit_shift));
            let mask = (1u64 << needed_bits) - 1;
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(result & mask))
        } else {
            let needed_bytes = (needed_bits + bit_shift).div_ceil(8);
            if needed_bytes <= 8 {
                let raw = if byte_offset + 8 <= self.strings.len() {
                    unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u64
                        )
                    }
                } else {
                    let mut buf = [0u8; 8];
                    let avail = self.strings.len() - byte_offset;
                    buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                    u64::from_le_bytes(buf)
                };
                let shifted = raw >> bit_shift;
                let mask = if needed_bits >= 64 { u64::MAX } else { (1u64 << needed_bits) - 1 };
                Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
            } else {
                // For K=63 (needed_bits=126), when bit_shift > 2 we need more
                // than 128 bits from the byte stream. Load u128 + extra byte.
                let (raw, extra_byte) = if byte_offset + 17 <= self.strings.len() {
                    let r = unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u128
                        )
                    };
                    (r, self.strings[byte_offset + 16])
                } else if byte_offset + 16 <= self.strings.len() {
                    let r = unsafe {
                        std::ptr::read_unaligned(
                            self.strings.as_ptr().add(byte_offset) as *const u128
                        )
                    };
                    let extra = if byte_offset + 16 < self.strings.len() {
                        self.strings[byte_offset + 16]
                    } else {
                        0u8
                    };
                    (r, extra)
                } else {
                    let mut buf = [0u8; 17];
                    let avail = self.strings.len() - byte_offset;
                    buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                    (u128::from_le_bytes(buf[..16].try_into().unwrap()), buf[16])
                };
                let shifted = if bit_shift > 0 {
                    (raw >> bit_shift) | ((extra_byte as u128) << (128 - bit_shift))
                } else {
                    raw
                };
                let mask = if needed_bits >= 128 { u128::MAX } else { (1u128 << needed_bits) - 1 };
                Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u128(shifted & mask))
            }
        }
    }

    /// Serialize the SPSS to a writer using a custom binary format.
    ///
    /// Format:
    /// - k: u64 (LE)
    /// - m: u64 (LE)
    /// - strings_len: u64 (LE)
    /// - strings: [u8; strings_len]
    /// - offsets: epserde Elias-Fano binary format
    pub fn serialize_to<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&(self.k as u64).to_le_bytes())?;
        writer.write_all(&(self.m as u64).to_le_bytes())?;
        // Write logical length (exclude 8-byte padding)
        let logical_len = self.strings.len() - 8;
        writer.write_all(&(logical_len as u64).to_le_bytes())?;
        writer.write_all(&self.strings[..logical_len])?;
        self.offsets.write_to(writer)?;
        Ok(())
    }

    /// Deserialize an SPSS from a reader.
    pub fn deserialize_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf8 = [0u8; 8];
        reader.read_exact(&mut buf8)?;
        let k = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let m = u64::from_le_bytes(buf8) as usize;
        reader.read_exact(&mut buf8)?;
        let strings_len = u64::from_le_bytes(buf8) as usize;
        let mut strings = vec![0u8; strings_len + 8]; // +8 padding for safe unaligned reads
        reader.read_exact(&mut strings[..strings_len])?;
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

    /// The `K <= 31` arm of `decode_kmer_at` folds the spill byte in without a
    /// branch. Pin that it computes exactly what the branching form did, over
    /// every `(bit_shift, K)` combination that can actually occur — `bit_shift`
    /// is `(pos % 4) * 2`, and K is odd in `3..=31` on this arm.
    #[test]
    fn branchless_spill_matches_branching_form() {
        fn branching(raw: u64, extra: u8, bit_shift: usize, needed_bits: usize) -> u64 {
            let mut result = raw >> bit_shift;
            if bit_shift > 2 {
                result |= (extra as u64) << (64 - bit_shift);
            }
            result & ((1u64 << needed_bits) - 1)
        }
        fn branchless(raw: u64, extra: u8, bit_shift: usize, needed_bits: usize) -> u64 {
            let result = (raw >> bit_shift) | (((extra as u64) << 1) << (63 - bit_shift));
            result & ((1u64 << needed_bits) - 1)
        }

        let mut s: u64 = 0x243F_6A88_85A3_08D3;
        for _ in 0..20_000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let extra = (s >> 11) as u8;
            for bit_shift in [0usize, 2, 4, 6] {
                for k in (3..=31).step_by(2) {
                    assert_eq!(
                        branching(s, extra, bit_shift, k * 2),
                        branchless(s, extra, bit_shift, k * 2),
                        "raw={s:#x} extra={extra:#x} bit_shift={bit_shift} k={k}"
                    );
                }
            }
        }
    }

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
