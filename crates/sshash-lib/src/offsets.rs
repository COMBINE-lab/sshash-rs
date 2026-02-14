//! Compact encoding of offsets into a bit-packed string set
//!
//! This module provides efficient storage of offsets using variable-length encoding,
//! reducing the memory footprint of the index.
//!
//! Two representations are available:
//! - `OffsetsVector`: Plain `Vec<u64>`, used during construction
//! - `EliasFanoOffsets`: Elias-Fano encoding via sux-rs `EfSeqDict`, used after
//!   construction and during queries. Provides O(1) random access and fast
//!   `locate()` via successor queries (matches C++ endpoints_sequence).

use epserde::prelude::*;

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
// Elias-Fano encoded offsets (sux-rs)
// ---------------------------------------------------------------------------

use sux::dict::elias_fano::{EliasFanoBuilder, EfSeqDict};
use sux::traits::{IndexedSeq, Succ};
use sux::traits::iter::BidiIterator;
use mem_dbg::{MemSize, SizeFlags};

/// Elias-Fano encoded offsets for compact, fast string boundary lookups.
///
/// Uses sux-rs `EfSeqDict` which provides:
/// - O(1) random access via `get_unchecked(i)` (uses Select1)
/// - O(1) successor query via `succ(x)` (uses Select0)
/// - Bidirectional iterator from successor via `iter_bidi_from_succ(x)`:
///   `next()`/`prev()` use `select_in_word` (single CPU instruction) for
///   cheap adjacent-element access without full Select1 inventory lookup.
///
/// This closely matches the C++ `endpoints_sequence` data structure used by SSHash.
///
/// Space usage is approximately `2 + log(U/N)` bits per element (Elias-Fano bound),
/// compared to 64 bits per element for `Vec<u64>`.
pub struct EliasFanoOffsets {
    /// Elias-Fano sequence for compact access and locate
    ef: EfSeqDict,
}

impl EliasFanoOffsets {
    /// Build from a sorted vector of offsets (must start with 0).
    pub fn from_vec(offsets: &[u64]) -> Self {
        let n = offsets.len();
        let u = if n > 0 { offsets[n - 1] as usize + 1 } else { 1 };
        let mut builder = EliasFanoBuilder::new(n, u);
        for &v in offsets {
            builder.push(v as usize);
        }
        Self { ef: builder.build_with_seq_and_dict() }
    }

    /// Build from an `OffsetsVector` (consumes it).
    pub fn from_offsets_vector(ov: OffsetsVector) -> Self {
        Self::from_vec(&ov.offsets)
    }

    /// Get the offset at index `i`.
    #[inline]
    pub fn access(&self, i: usize) -> u64 {
        // SAFETY: caller must ensure i < self.len().
        unsafe { self.ef.get_unchecked(i) as u64 }
    }

    /// Locate which string contains a given absolute position, returning
    /// `(string_id, string_begin, string_end)`.
    ///
    /// Uses sux-rs `iter_bidi_from_succ()` to find the successor, then
    /// cheap `next()`/`prev()` calls (single-instruction bit scans) to
    /// read adjacent elements without full Select1 inventory lookups.
    #[inline]
    pub fn locate_with_end(&self, pos: u64) -> Option<(u64, u64, u64)> {
        let n = self.ef.len();
        if n < 2 {
            return None;
        }

        // iter_bidi_from_succ returns (index, positioned_iterator) for the
        // first element >= pos. The first next() yields the successor value.
        let (idx, mut iter) = self.ef.iter_bidi_from_succ(pos as usize)?;

        // Get the successor value (cheap: reads from cached word).
        let val = iter.next()?;

        if val == pos as usize {
            // Exact hit: pos is at a string boundary → string_id = idx.
            // Need the NEXT element for string_end (cheap bit scan forward).
            if idx + 1 < n {
                let end = iter.next()? as u64;
                Some((idx as u64, pos, end))
            } else {
                None
            }
        } else {
            // val > pos → string containing pos starts at idx-1.
            // val IS the end of this string. Need begin = offsets[idx-1].
            // prev() undoes the next(), then prev() again gets offsets[idx-1].
            // Both use select_in_word (single CPU instruction per call).
            debug_assert!(idx > 0);
            iter.prev(); // back to idx (returns val, discarded)
            let begin = iter.prev()? as u64; // offsets[idx-1]
            Some(((idx - 1) as u64, begin, val as u64))
        }
    }

    /// Locate which string contains a given absolute position.
    /// Returns `(string_id, string_begin)` where
    /// `offsets[string_id] <= pos < offsets[string_id + 1]`.
    ///
    /// Uses `iter_bidi_from_succ()` + cheap `prev()` bit scans
    /// instead of full Select1 for adjacent element access.
    #[inline]
    pub fn locate(&self, pos: u64) -> Option<(u64, u64)> {
        let n = self.ef.len();
        if n < 2 {
            return None;
        }

        let (idx, mut iter) = self.ef.iter_bidi_from_succ(pos as usize)?;
        let val = iter.next()?;

        if val == pos as usize {
            // Exact hit: string_id = idx, but only if there's a next element
            if idx + 1 < n {
                Some((idx as u64, pos))
            } else {
                None
            }
        } else {
            // val > pos: string containing pos starts at idx - 1
            debug_assert!(idx > 0);
            iter.prev(); // back to idx
            let string_begin = iter.prev()? as u64; // offsets[idx-1]
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
        self.ef.len() == 0
    }

    /// Actual number of bytes used by the Elias-Fano index (including selection structures).
    #[inline]
    pub fn num_bytes(&self) -> u64 {
        self.ef.mem_size(SizeFlags::default()) as u64
    }

    /// Actual number of bits used by the Elias-Fano index.
    #[inline]
    pub fn num_bits(&self) -> u64 {
        self.num_bytes() * 8
    }

    /// Serialize the Elias-Fano offsets to a writer using epserde's binary format.
    pub fn write_to<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        unsafe {
            self.ef.serialize(writer)
                .map_err(std::io::Error::other)?;
        }
        Ok(())
    }

    /// Deserialize Elias-Fano offsets from a reader using epserde's binary format.
    pub fn read_from<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let ef = unsafe {
            EfSeqDict::deserialize_full(reader)
                .map_err(std::io::Error::other)?
        };
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
