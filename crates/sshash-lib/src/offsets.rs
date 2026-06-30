//! Compact encoding of offsets into a bit-packed string set
//!
//! This module provides efficient storage of offsets using variable-length encoding,
//! reducing the memory footprint of the index.
//!
//! Two representations are available:
//! - `OffsetsVector`: Plain `Vec<u64>`, used during construction
//! - `EliasFanoOffsets`: Elias-Fano with DArray Select, used after
//!   construction and during queries. Provides O(1) random access and fast
//!   `locate()` via successor queries (matches C++ endpoints_sequence).

/// A decoded offset with both absolute and relative information
#[derive(Clone, Copy, Debug)]
pub struct DecodedOffset {
    /// Absolute byte offset into the string data
    pub absolute_offset: u64,
    /// Relative offset (for retrieval purposes)
    pub relative_offset: u64,
}

impl DecodedOffset {
    /// Create a new decoded offset
    pub fn new(absolute_offset: u64, relative_offset: u64) -> Self {
        Self {
            absolute_offset,
            relative_offset,
        }
    }
}

/// Compact vector of offsets
///
/// Stores offsets using a compact representation.
#[derive(Clone, Debug)]
pub struct OffsetsVector {
    /// Raw offset values
    offsets: Vec<u64>,
}

impl OffsetsVector {
    /// Create a new empty offsets vector
    pub fn new() -> Self {
        Self {
            offsets: vec![0], // Start with 0 for the first offset
        }
    }

    /// Add an offset to the vector
    #[inline]
    pub fn push(&mut self, offset: u64) {
        self.offsets.push(offset);
    }

    /// Get the offset at index `i`
    #[inline]
    pub fn access(&self, i: usize) -> u64 {
        assert!(i < self.offsets.len(), "Offset index {} out of bounds", i);
        self.offsets[i]
    }

    /// Decode an offset (currently identity since we don't compress yet)
    #[inline]
    pub fn decode(&self, absolute_offset: u64) -> DecodedOffset {
        DecodedOffset::new(absolute_offset, absolute_offset)
    }

    /// Get the number of offsets
    #[inline]
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Check if the vector is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Get the number of bytes used (approximation for MVP)
    #[inline]
    pub fn num_bytes(&self) -> u64 {
        (self.offsets.len() * 8) as u64
    }

    /// Get the number of bits used (approximation for MVP)
    #[inline]
    pub fn num_bits(&self) -> u64 {
        (self.offsets.len() * 64) as u64
    }

    /// Binary search to find which string contains a given absolute position.
    /// Returns `(string_id, string_begin)` where `offsets[string_id] <= pos < offsets[string_id + 1]`.
    /// This matches the C++ `decoded_offsets::offset_to_id` / Elias-Fano `locate` approach.
    #[inline]
    pub fn locate(&self, pos: u64) -> Option<(u64, u64)> {
        let n = self.offsets.len();
        if n < 2 {
            return None;
        }

        // Use Rust's optimised binary search: partition_point returns the first
        // index where the predicate is false, i.e. the first offset > pos.
        let idx = self.offsets.partition_point(|&x| x <= pos);

        // idx == 0 means pos < offsets[0]: out of bounds.
        if idx == 0 {
            return None;
        }
        let string_id = idx - 1;

        // Validate: pos must be within [offsets[string_id], offsets[string_id + 1])
        if string_id + 1 < n {
            Some((string_id as u64, self.offsets[string_id]))
        } else {
            None
        }
    }

    /// Branchless binary search variant for benchmarking comparison.
    ///
    /// Uses conditional moves instead of branches to avoid branch mispredictions.
    /// The inner loop has no data-dependent branches - only a CMOV.
    #[inline]
    pub fn locate_branchless(&self, pos: u64) -> Option<(u64, u64)> {
        let n = self.offsets.len();
        if n < 2 {
            return None;
        }

        let data = self.offsets.as_slice();

        // Branchless binary search: find rightmost index where data[idx] <= pos
        let mut lo: usize = 0;
        let mut size = n;
        while size > 1 {
            let half = size / 2;
            let mid = lo + half;
            // Branchless: compiler should emit CMOV
            // SAFETY: mid is always < n because lo + half < lo + size <= n
            lo = if unsafe { *data.get_unchecked(mid) } <= pos { mid } else { lo };
            size -= half;
        }

        // lo is now the rightmost index where data[lo] <= pos, or 0 if pos < data[0]
        if unsafe { *data.get_unchecked(lo) } > pos {
            return None;
        }
        if lo + 1 < n {
            Some((lo as u64, unsafe { *data.get_unchecked(lo) }))
        } else {
            None
        }
    }
}

impl Default for OffsetsVector {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_offsets_vector_creation() {
        let offsets = OffsetsVector::new();
        assert_eq!(offsets.len(), 1);
        assert_eq!(offsets.access(0), 0);
    }

    #[test]
    fn test_offsets_vector_push() {
        let mut offsets = OffsetsVector::new();
        offsets.push(100);
        offsets.push(200);
        offsets.push(300);

        assert_eq!(offsets.len(), 4);
        assert_eq!(offsets.access(0), 0);
        assert_eq!(offsets.access(1), 100);
        assert_eq!(offsets.access(2), 200);
        assert_eq!(offsets.access(3), 300);
    }

    #[test]
    fn test_offsets_vector_decode() {
        let offsets = OffsetsVector::new();
        let decoded = offsets.decode(50);
        assert_eq!(decoded.absolute_offset, 50);
        assert_eq!(decoded.relative_offset, 50);
    }

    #[test]
    fn test_decoded_offset_creation() {
        let decoded = DecodedOffset::new(1000, 500);
        assert_eq!(decoded.absolute_offset, 1000);
        assert_eq!(decoded.relative_offset, 500);
    }
}

// ---------------------------------------------------------------------------
// Elias-Fano encoded offsets with DArray Select support
// ---------------------------------------------------------------------------

use crate::darray::DArray;

/// Elias-Fano encoded offsets with DArray-accelerated access and locate.
///
/// Custom Elias-Fano representation with DArray Select1 (for access) and
/// DArray Select0 (for locate/successor queries). Serialized directly
/// without sux-rs dependency.
///
/// Space: ~2 + log(U/N) bits per element (EF low bits) + ~1.5625 bits per
/// element (DArray overhead). For gencode (1.08M strings, 95M universe):
/// ~1.3 MB total vs ~9.7 MB for pre-decoded Vec.
pub struct EliasFanoOffsets {
    num_values: usize,
    low_bit_width: u32,
    low_bits: Vec<u64>,
    high_bits: Vec<u64>,
    high_bits_len: usize,
    num_high_zeros: usize,
    darray1: DArray,
    darray0: DArray,
}

impl EliasFanoOffsets {
    /// Build from a sorted vector of offsets (must start with 0).
    pub fn from_vec(offsets: &[u64]) -> Self {
        let (num_values, low_bit_width, low_bits, high_bits, high_bits_len, num_high_zeros, darray1, darray0) =
            Self::build_darray_ef(offsets);
        Self { num_values, low_bit_width, low_bits, high_bits, high_bits_len, num_high_zeros, darray1, darray0 }
    }

    /// Build from an `OffsetsVector` (consumes it).
    pub fn from_offsets_vector(ov: OffsetsVector) -> Self {
        Self::from_vec(&ov.offsets)
    }

    fn build_darray_ef(values: &[u64]) -> (usize, u32, Vec<u64>, Vec<u64>, usize, usize, DArray, DArray) {
        let n = values.len();
        if n == 0 {
            let empty_bv = vec![0u64; 2];
            let da1 = DArray::build(&empty_bv, 0, false);
            let da0 = DArray::build(&empty_bv, 0, true);
            return (0, 0, Vec::new(), empty_bv, 0, 0, da1, da0);
        }

        let max_val = values[n - 1];

        let low_bit_width = if n <= 1 || max_val == 0 {
            0
        } else {
            let ratio = (max_val + 1) / n as u64;
            if ratio == 0 { 0 } else { 63 - ratio.leading_zeros() }
        };

        let low_mask = if low_bit_width > 0 { (1u64 << low_bit_width) - 1 } else { 0 };

        // Pack low bits
        let total_low_bits = n as u64 * low_bit_width as u64;
        let low_words = ((total_low_bits + 63) / 64) as usize;
        let mut low_bits = vec![0u64; low_words];
        if low_bit_width > 0 {
            for (i, &v) in values.iter().enumerate() {
                set_packed_bits(&mut low_bits, i, low_bit_width, v & low_mask);
            }
        }

        // Build high bitvector
        let max_high = max_val >> low_bit_width;
        let high_bits_len = n + max_high as usize;
        let high_words = high_bits_len / 64 + 2; // +2 for scan padding
        let mut high_bits = vec![0u64; high_words];

        for (i, &v) in values.iter().enumerate() {
            let high = (v >> low_bit_width) as usize;
            let bit_pos = high + i;
            high_bits[bit_pos / 64] |= 1u64 << (bit_pos % 64);
        }

        let num_high_zeros = max_high as usize;
        let darray1 = DArray::build(&high_bits, high_bits_len, false);
        let darray0 = DArray::build(&high_bits, high_bits_len, true);

        (n, low_bit_width, low_bits, high_bits, high_bits_len, num_high_zeros, darray1, darray0)
    }

    #[inline(always)]
    fn get_low_bits(&self, i: usize) -> u64 {
        if self.low_bit_width == 0 {
            return 0;
        }
        let bit_pos = i as u64 * self.low_bit_width as u64;
        let word_idx = (bit_pos >> 6) as usize;
        let bit_offset = (bit_pos & 63) as u32;
        let mask = (1u64 << self.low_bit_width) - 1;
        let mut value = (self.low_bits[word_idx] >> bit_offset) & mask;
        if bit_offset + self.low_bit_width > 64 {
            let overflow_bits = bit_offset + self.low_bit_width - 64;
            value |= (self.low_bits[word_idx + 1] & ((1u64 << overflow_bits) - 1))
                << (self.low_bit_width - overflow_bits);
        }
        value
    }

    /// Get the offset at index `i`.
    #[inline]
    pub fn access(&self, i: usize) -> u64 {
        debug_assert!(i < self.num_values);
        let high = self.darray1.select(&self.high_bits, i as u64, false) - i as u64;
        let low = self.get_low_bits(i);
        (high << self.low_bit_width) | low
    }

    /// Locate which string contains a given absolute position, returning
    /// `(string_id, string_begin, string_end)`.
    #[inline]
    pub fn locate_with_end(&self, pos: u64) -> Option<(u64, u64, u64)> {
        let n = self.num_values;
        if n < 2 {
            return None;
        }

        let h = pos >> self.low_bit_width;

        let start = if h == 0 {
            0
        } else if (h - 1) < self.num_high_zeros as u64 {
            let s0 = self.darray0.select(&self.high_bits, h - 1, true);
            (s0 + 1 - h) as usize
        } else {
            return None;
        };

        // Linear scan from `start`, caching the two boundary offsets as we go.
        // Each `access` is a full Elias-Fano select; the original recomputed
        // `access(string_id)` and `access(string_id + 1)` in the return even
        // though the scan had just evaluated them. We keep the last `<= pos`
        // value (`prev_access`, which becomes `access(string_id)`) and the
        // first `> pos` value that breaks the loop (`break_access`, which is
        // `access(string_id + 1)`), so a successful locate does ~half as many
        // selects. Results are identical (verified by the brute-force stress
        // tests below).
        let mut idx = start;
        let mut prev_access = 0u64;
        let mut break_access = 0u64;
        while idx < n {
            let a = self.access(idx);
            if a > pos {
                break_access = a;
                break;
            }
            prev_access = a;
            idx += 1;
        }

        if idx == 0 {
            return None;
        }
        let string_id = idx - 1;
        if string_id + 1 < n {
            // Reaching here means `idx == string_id + 1 < n`, so the scan broke
            // via `a > pos` and `break_access == access(string_id + 1)`.
            // `access(string_id)` is `prev_access` when the loop advanced at
            // least once (`idx > start`); on an immediate break the index
            // `string_id == start - 1` was never visited, so read it directly.
            let begin = if idx > start { prev_access } else { self.access(string_id) };
            Some((string_id as u64, begin, break_access))
        } else {
            None
        }
    }

    /// Locate which string contains a given absolute position.
    #[inline]
    pub fn locate(&self, pos: u64) -> Option<(u64, u64)> {
        let n = self.num_values;
        if n < 2 {
            return None;
        }

        let h = pos >> self.low_bit_width;

        let start = if h == 0 {
            0
        } else if (h - 1) < self.num_high_zeros as u64 {
            let s0 = self.darray0.select(&self.high_bits, h - 1, true);
            (s0 + 1 - h) as usize
        } else {
            return None;
        };

        let mut idx = start;
        while idx < n && self.access(idx) <= pos {
            idx += 1;
        }

        if idx == 0 {
            return None;
        }
        let string_id = idx - 1;
        if string_id + 1 < n {
            Some((string_id as u64, self.access(string_id)))
        } else {
            None
        }
    }

    /// Get the number of offsets.
    #[inline]
    pub fn len(&self) -> usize {
        self.num_values
    }

    /// Check if the vector is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.num_values == 0
    }

    /// Get the in-memory size in bytes (DArray representation).
    #[inline]
    pub fn num_bytes(&self) -> u64 {
        let low = self.low_bits.len() * 8;
        let high = self.high_bits.len() * 8;
        let da1 = self.darray1.heap_bytes();
        let da0 = self.darray0.heap_bytes();
        (low + high + da1 + da0) as u64
    }

    /// Get the in-memory size in bits.
    #[inline]
    pub fn num_bits(&self) -> u64 {
        self.num_bytes() * 8
    }

    /// Serialize to a writer. Writes the DArray components directly.
    pub fn write_to<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        use crate::darray::{write_u64, write_u32, write_vec_u64};
        write_u64(writer, self.num_values as u64)?;
        write_u32(writer, self.low_bit_width)?;
        write_vec_u64(writer, &self.low_bits)?;
        write_u64(writer, self.high_bits_len as u64)?;
        write_vec_u64(writer, &self.high_bits)?;
        write_u64(writer, self.num_high_zeros as u64)?;
        self.darray1.write_to(writer)?;
        self.darray0.write_to(writer)?;
        Ok(())
    }

    /// Deserialize from a reader. Reads pre-built DArray components directly.
    pub fn read_from<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        use crate::darray::{read_u64, read_u32, read_vec_u64};
        let num_values = read_u64(reader)? as usize;
        let low_bit_width = read_u32(reader)?;
        let low_bits = read_vec_u64(reader)?;
        let high_bits_len = read_u64(reader)? as usize;
        let high_bits = read_vec_u64(reader)?;
        let num_high_zeros = read_u64(reader)? as usize;
        let darray1 = DArray::read_from(reader)?;
        let darray0 = DArray::read_from(reader)?;
        Ok(Self { num_values, low_bit_width, low_bits, high_bits, high_bits_len, num_high_zeros, darray1, darray0 })
    }
}

impl std::fmt::Debug for EliasFanoOffsets {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EliasFanoOffsets")
            .field("len", &self.num_values)
            .field("low_bit_width", &self.low_bit_width)
            .field("darray_bytes", &self.num_bytes())
            .finish()
    }
}

#[inline]
fn set_packed_bits(words: &mut [u64], index: usize, width: u32, value: u64) {
    let bit_pos = index as u64 * width as u64;
    let word_idx = (bit_pos >> 6) as usize;
    let bit_offset = (bit_pos & 63) as u32;
    words[word_idx] |= value << bit_offset;
    if bit_offset + width > 64 {
        let overflow_bits = bit_offset + width - 64;
        words[word_idx + 1] |= value >> (width - overflow_bits);
    }
}

#[cfg(test)]
mod ef_tests {
    use super::*;

    #[test]
    fn test_ef_from_vec() {
        let offsets = vec![0, 100, 200, 300, 400];
        let ef = EliasFanoOffsets::from_vec(&offsets);
        assert_eq!(ef.len(), 5);
        assert_eq!(ef.access(0), 0);
        assert_eq!(ef.access(1), 100);
        assert_eq!(ef.access(2), 200);
        assert_eq!(ef.access(3), 300);
        assert_eq!(ef.access(4), 400);
    }

    #[test]
    fn test_ef_locate() {
        let offsets = vec![0, 100, 200, 300, 400];
        let ef = EliasFanoOffsets::from_vec(&offsets);

        // pos=50 → string 0 (begins at 0)
        assert_eq!(ef.locate(50), Some((0, 0)));
        // pos=100 → exact boundary, string 1 (begins at 100)
        assert_eq!(ef.locate(100), Some((1, 100)));
        // pos=199 → string 1 (begins at 100)
        assert_eq!(ef.locate(199), Some((1, 100)));
        // pos=300 → exact boundary, string 3 (begins at 300)
        assert_eq!(ef.locate(300), Some((3, 300)));
        // pos=399 → string 3 (begins at 300)
        assert_eq!(ef.locate(399), Some((3, 300)));
        // pos=400 → out of range (last element is universe bound)
        assert_eq!(ef.locate(400), None);
    }

    #[test]
    fn test_ef_locate_with_end() {
        let offsets = vec![0, 100, 200, 300, 400];
        let ef = EliasFanoOffsets::from_vec(&offsets);

        // pos=50 → string 0 (begins at 0, ends at 100)
        assert_eq!(ef.locate_with_end(50), Some((0, 0, 100)));
        // pos=100 → exact boundary, string 1 (begins at 100, ends at 200)
        assert_eq!(ef.locate_with_end(100), Some((1, 100, 200)));
        // pos=199 → string 1 (begins at 100, ends at 200)
        assert_eq!(ef.locate_with_end(199), Some((1, 100, 200)));
        // pos=300 → exact boundary, string 3 (begins at 300, ends at 400)
        assert_eq!(ef.locate_with_end(300), Some((3, 300, 400)));
        // pos=399 → string 3 (begins at 300, ends at 400)
        assert_eq!(ef.locate_with_end(399), Some((3, 300, 400)));
        // pos=400 → out of range
        assert_eq!(ef.locate_with_end(400), None);
    }

    #[test]
    fn test_ef_serialization_roundtrip() {
        let offsets = vec![0, 100, 200, 300, 400];
        let ef = EliasFanoOffsets::from_vec(&offsets);

        let mut buf = Vec::new();
        ef.write_to(&mut buf).unwrap();

        let ef2 = EliasFanoOffsets::read_from(&mut &buf[..]).unwrap();
        assert_eq!(ef2.len(), 5);
        for i in 0..5 {
            assert_eq!(ef.access(i), ef2.access(i));
        }
        assert_eq!(ef2.locate(150), Some((1, 100)));
        assert_eq!(ef2.locate_with_end(150), Some((1, 100, 200)));
    }

    /// Stress test: build EF with varying gaps and verify locate_with_end
    /// against a brute-force reference for EVERY position in the range.
    #[test]
    fn test_ef_locate_with_end_stress() {
        // Create offsets with varying gap sizes to exercise different EF bit patterns
        let mut offsets = vec![0u64];
        let gaps = [3, 7, 1, 15, 2, 100, 5, 31, 8, 63, 4, 127, 6, 255, 10, 50,
                    1, 1, 1, 33, 200, 9, 17, 3, 11, 500, 2, 7, 13, 41];
        for &g in gaps.iter() {
            offsets.push(offsets.last().unwrap() + g);
        }
        let ef = EliasFanoOffsets::from_vec(&offsets);

        // Verify access
        for (i, &v) in offsets.iter().enumerate() {
            assert_eq!(ef.access(i), v, "access({i}) mismatch");
        }

        // Reference locate_with_end via brute force
        let universe = *offsets.last().unwrap();
        for pos in 0..=universe {
            let expected = {
                // Find string containing pos: offsets[id] <= pos < offsets[id+1]
                let mut found = None;
                for i in 0..offsets.len() - 1 {
                    if offsets[i] <= pos && pos < offsets[i + 1] {
                        found = Some((i as u64, offsets[i], offsets[i + 1]));
                        break;
                    }
                }
                found
            };
            let got = ef.locate_with_end(pos);
            assert_eq!(got, expected, "locate_with_end({pos}) mismatch");
        }

        // Also test past-the-end
        assert_eq!(ef.locate_with_end(universe), None);
        assert_eq!(ef.locate_with_end(universe + 1), None);
    }

    /// Stress test with large gaps to exercise high-bit patterns
    #[test]
    fn test_ef_locate_large_universe() {
        let offsets: Vec<u64> = vec![0, 1000, 5000, 5001, 5002, 10000, 100000, 100001, 500000];
        let ef = EliasFanoOffsets::from_vec(&offsets);

        // Test all boundary positions and a few interior ones
        let test_positions: Vec<u64> = {
            let mut v = Vec::new();
            for &off in &offsets {
                if off > 0 { v.push(off - 1); }
                v.push(off);
                v.push(off + 1);
            }
            // Some random interior positions
            v.extend_from_slice(&[500, 3000, 5000, 7500, 50000, 200000, 400000]);
            v.sort();
            v.dedup();
            v
        };

        for pos in test_positions {
            let expected = {
                let mut found = None;
                for i in 0..offsets.len() - 1 {
                    if offsets[i] <= pos && pos < offsets[i + 1] {
                        found = Some((i as u64, offsets[i], offsets[i + 1]));
                        break;
                    }
                }
                found
            };
            let got = ef.locate_with_end(pos);
            assert_eq!(got, expected, "locate_with_end({pos}) mismatch");
        }
    }
}
