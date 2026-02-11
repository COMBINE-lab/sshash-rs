//! Compact encoding of offsets into a bit-packed string set
//!
//! This module provides efficient storage of offsets using variable-length encoding,
//! reducing the memory footprint of the index.
//!
//! Two representations are available:
//! - `OffsetsVector`: Plain `Vec<u64>`, used during construction
//! - `EliasFanoOffsets`: Elias-Fano encoding via cseq `Sequence`, used after
//!   construction and during queries. Provides O(1) random access and fast
//!   `locate()` via successor queries with stateful Cursor for cheap adjacent
//!   element access (matches C++ endpoints_sequence).

use epserde::prelude::*;
use dyn_size_of::GetSize;

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
#[derive(Clone, Debug, Epserde)]
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
// Elias-Fano encoded offsets (cseq)
// ---------------------------------------------------------------------------

/// The cseq Elias-Fano sequence type.
pub type CseqSequence = cseq::elias_fano::Sequence;

/// Elias-Fano encoded offsets for compact, fast string boundary lookups.
///
/// Uses cseq's `Sequence` with a stateful `Cursor` interface:
/// - O(1) random access via `get_unchecked(i)` (uses Select1)
/// - O(1) successor query via `geq_cursor(x)` (uses Select0)
/// - Cheap adjacent element access via `Cursor::advance()` / `advance_back()`
///   (just a bit scan + lo fragment read, ~5-10ns)
///
/// This closely matches the C++ `endpoints_sequence` data structure used by SSHash,
/// where `locate()` returns a positioned iterator that reads adjacent elements cheaply.
///
/// Space usage is approximately `2 + log(U/N)` bits per element (Elias-Fano bound),
/// compared to 64 bits per element for `Vec<u64>`.
pub struct EliasFanoOffsets {
    /// Elias-Fano sequence for compact access and locate
    ef: CseqSequence,
}

impl EliasFanoOffsets {
    /// Build from a sorted vector of offsets (must start with 0).
    pub fn from_vec(offsets: &[u64]) -> Self {
        let ef = CseqSequence::with_items_from_slice(offsets);
        Self { ef }
    }

    /// Build from an `OffsetsVector` (consumes it).
    pub fn from_offsets_vector(ov: OffsetsVector) -> Self {
        Self::from_vec(&ov.offsets)
    }

    /// Get the offset at index `i`.
    #[inline]
    pub fn access(&self, i: usize) -> u64 {
        // SAFETY: caller must ensure i < self.len().
        unsafe { self.ef.get_unchecked(i) }
    }

    /// Locate which string contains a given absolute position, returning
    /// `(string_id, string_begin, string_end)` in a single EF traversal.
    ///
    /// Uses cseq's `geq_cursor()` to find the successor, then `advance_back()`
    /// or `advance()` to read the adjacent element cheaply via the positioned
    /// Cursor (just a bit scan, no Select1). This matches the C++
    /// `endpoints_sequence::locate` pattern exactly.
    #[inline]
    pub fn locate_with_end(&self, pos: u64) -> Option<(u64, u64, u64)> {
        let n = self.ef.len();
        if n < 2 {
            return None;
        }

        // geq_cursor returns a Cursor pointing to the first element >= pos.
        // If pos >= all elements, the cursor is past-the-end (is_end()).
        let mut cursor = self.ef.geq_cursor(pos);

        if cursor.is_end() {
            return None;
        }

        // SAFETY: cursor is valid (not past-the-end)
        let val = unsafe { cursor.value_unchecked() };
        let idx = cursor.index();

        if val == pos {
            // Exact hit: pos is at a string boundary → string_id = idx.
            // Need the NEXT element for string_end.
            // Note: advance() returns true even when moving to past-the-end,
            // so we guard with an index check instead.
            if idx + 1 < n {
                cursor.advance();
                let end = unsafe { cursor.value_unchecked() };
                Some((idx as u64, val, end))
            } else {
                // idx is the last element — no string after this boundary
                None
            }
        } else {
            // val > pos → string containing pos starts at idx-1.
            // val IS the end of this string (offsets[idx]).
            // Need the PREVIOUS element for string_begin — cheap via advance_back.
            debug_assert!(idx > 0);
            let end = val;
            cursor.advance_back(); // goes back to idx-1 — just a bit scan
            let string_begin = unsafe { cursor.value_unchecked() };
            Some(((idx - 1) as u64, string_begin, end))
        }
    }

    /// Locate which string contains a given absolute position.
    /// Returns `(string_id, string_begin)` where
    /// `offsets[string_id] <= pos < offsets[string_id + 1]`.
    ///
    /// Uses the Elias-Fano successor query for O(1) lookup,
    /// matching C++ `endpoints_sequence::locate`.
    #[inline]
    pub fn locate(&self, pos: u64) -> Option<(u64, u64)> {
        let n = self.ef.len();
        if n < 2 {
            return None;
        }

        let mut cursor = self.ef.geq_cursor(pos);

        if cursor.is_end() {
            return None;
        }

        let val = unsafe { cursor.value_unchecked() };
        let idx = cursor.index();

        if val == pos {
            // Exact hit: string_id = idx, but only if there's a next element
            if idx + 1 < n {
                Some((idx as u64, val))
            } else {
                None
            }
        } else {
            // val > pos: string containing pos starts at idx - 1
            debug_assert!(idx > 0);
            cursor.advance_back();
            let string_begin = unsafe { cursor.value_unchecked() };
            Some(((idx - 1) as u64, string_begin))
        }
    }

    /// Number of offsets stored.
    #[inline]
    pub fn len(&self) -> usize {
        self.ef.len()
    }

    /// Whether there are no offsets.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.ef.is_empty()
    }

    /// Actual number of bytes used by the Elias-Fano index (including selection structures).
    #[inline]
    pub fn num_bytes(&self) -> u64 {
        self.ef.size_bytes() as u64
    }

    /// Actual number of bits used by the Elias-Fano index.
    #[inline]
    pub fn num_bits(&self) -> u64 {
        self.num_bytes() * 8
    }

    /// Serialize the Elias-Fano offsets to a writer using cseq's binary format.
    pub fn write_to(&self, writer: &mut dyn std::io::Write) -> std::io::Result<()> {
        self.ef.write(writer)
    }

    /// Deserialize Elias-Fano offsets from a reader using cseq's binary format.
    pub fn read_from(reader: &mut dyn std::io::Read) -> std::io::Result<Self> {
        let ef = CseqSequence::read(reader)?;
        Ok(Self { ef })
    }
}

impl std::fmt::Debug for EliasFanoOffsets {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EliasFanoOffsets")
            .field("len", &self.ef.len())
            .finish()
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
