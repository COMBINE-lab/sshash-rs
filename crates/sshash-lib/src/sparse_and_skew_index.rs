//! Sparse and Skew Index for k-mer lookup
//!
//! Buckets (one per minimizer, in MPHF order) use one of two layouts:
//! - `Flat` (format 6.0, and the fallback): every super-k-mer offset in one
//!   array, bucket after bucket, delimited by an Elias-Fano of the bucket
//!   boundaries, so each lookup pays an Elias-Fano select.
//! - `Inline` (format 7.0): each bucket's first offset in a direct array.
//!   A bitmap with sampled ranks marks the multi-position buckets, whose
//!   remaining offsets sit in a second array delimited by an Elias-Fano. One
//!   flag per block of [`HEAVY_BLOCK`] buckets marks blocks holding a heavy
//!   bucket.
//!
//! A build keeps whichever layout is smaller in total. A loaded 6.0 index stays
//! flat; [`SparseAndSkewIndex::upgrade_layout`] converts it explicitly.
//!
//! Inline lookup path:
//!   bucket_id = mphf(minimizer), first = first[bucket_id]
//!   singleton: `first` is the only position
//!   multi-position, block without heavy buckets: try `first`; only if that
//!     fails, decode the rest (rank, Elias-Fano pair) and scan it
//!   multi-position, heavy-flagged block: decode the bounds at once, so a
//!     heavy bucket goes straight to the skew index
//! - `SkewIndex`: partitioned MPHFs for heavy buckets (size > MIN_BUCKET_SIZE)

use crate::builder::buckets::{Bucket, ClassifiedBuckets, MIN_BUCKET_SIZE};
use crate::builder::dictionary_builder::BucketMetadata;
use crate::builder::external_sort::FileTuples;
use crate::builder::minimizer_tuples::MinimizerTuple;
use crate::constants::{ceil_log2, MIN_L, MAX_L, INVALID_UINT64};
use crate::kmer::{Kmer, KmerBits};
use crate::partitioned_mphf::PartitionedMphf;
use crate::offsets::EliasFanoOffsets;
use crate::mphf_config::MphfKey;
use sux::bits::bit_field_vec::BitFieldVec;
use tracing::info;
use value_traits::slices::{SliceByValue, SliceByValueMut};

/// Buckets per heavy-block flag of the inline layout.
pub const HEAVY_BLOCK: usize = 64;
/// Bitmap words per sampled rank (one `u64` per 512 buckets).
const RANK_WORDS: usize = 8;
/// Serialized layout tags (format 7.0).
const LAYOUT_FLAT: u64 = 0;
const LAYOUT_INLINE: u64 = 1;

/// Sparse and skew index: where each bucket's super-k-mers are, plus the skew
/// index of the heavy buckets.
pub struct SparseAndSkewIndex {
    layout: Layout,

    /// Skew index for heavy buckets (size > MIN_BUCKET_SIZE)
    pub skew_index: SkewIndex,
}

enum Layout {
    /// `extras[i]` = sum over buckets j < i of (size_j - 1); bucket i spans
    /// `offsets[extras[i] + i .. extras[i + 1] + i + 1]`.
    Flat { extras: EliasFanoOffsets, offsets: BitFieldVec },
    Inline(InlineFirst),
}

struct InlineFirst {
    /// Each bucket's first super-k-mer offset.
    first: BitFieldVec,
    /// Bit i: bucket i has more than one position.
    multi: Vec<u64>,
    /// Set bits of `multi` before each block of [`RANK_WORDS`] words.
    /// Rebuilt on load, not serialized.
    ranks: Vec<u64>,
    /// Bit j: a bucket in `[j * HEAVY_BLOCK, (j + 1) * HEAVY_BLOCK)` is heavy.
    heavy: Vec<u64>,
    /// Prefix sums of (size - 1) over the multi-position buckets.
    rest_bounds: EliasFanoOffsets,
    /// The other offsets of the multi-position buckets, bucket after bucket.
    rest: BitFieldVec,
}

#[inline(always)]
fn is_heavy(size: usize) -> bool {
    size > MIN_BUCKET_SIZE
}

fn sample_ranks(words: &[u64]) -> Vec<u64> {
    let mut acc = 0u64;
    words
        .chunks(RANK_WORDS)
        .map(|block| {
            let r = acc;
            acc += block.iter().map(|w| w.count_ones() as u64).sum::<u64>();
            r
        })
        .collect()
}

/// Bucket sizes of a flat layout, from its boundaries.
fn flat_sizes(extras: &EliasFanoOffsets) -> impl Iterator<Item = usize> + '_ {
    let mut it = extras.iter();
    let mut prev = it.next().unwrap_or(0);
    it.map(move |x| {
        let size = (x - prev) as usize + 1;
        prev = x;
        size
    })
}

impl InlineFirst {
    /// The directory of buckets of the given sizes (MPHF order); `first` and
    /// `rest` are left empty for the caller to fill.
    fn directory(num_buckets: usize, sizes: impl Iterator<Item = usize>) -> Self {
        let mut multi = vec![0u64; num_buckets.div_ceil(64)];
        let mut heavy = vec![0u64; num_buckets.div_ceil(64 * HEAVY_BLOCK)];
        let mut bounds = vec![0u64];
        for (i, size) in sizes.enumerate() {
            if size > 1 {
                multi[i / 64] |= 1 << (i % 64);
                bounds.push(bounds.last().unwrap() + size as u64 - 1);
                if is_heavy(size) {
                    let block = i / HEAVY_BLOCK;
                    heavy[block / 64] |= 1 << (block % 64);
                }
            }
        }
        Self {
            first: BitFieldVec::new(1, 0),
            ranks: sample_ranks(&multi),
            multi,
            heavy,
            rest_bounds: EliasFanoOffsets::from_vec(&bounds),
            rest: BitFieldVec::new(1, 0),
        }
    }

    /// Bytes beyond the offsets themselves.
    fn directory_bytes(&self) -> u64 {
        ((self.multi.len() + self.ranks.len() + self.heavy.len()) * 8) as u64
            + self.rest_bounds.num_bytes()
    }

    #[inline(always)]
    fn is_multi(&self, i: usize) -> bool {
        self.multi[i / 64] >> (i % 64) & 1 != 0
    }

    #[inline(always)]
    fn in_heavy_block(&self, i: usize) -> bool {
        let block = i / HEAVY_BLOCK;
        self.heavy[block / 64] >> (block % 64) & 1 != 0
    }

    /// Multi-position buckets before bucket `i`.
    #[inline(always)]
    fn rank(&self, i: usize) -> usize {
        let w = i / 64;
        let block = w / RANK_WORDS;
        let mut r = self.ranks[block];
        for &x in &self.multi[block * RANK_WORDS..w] {
            r += x.count_ones() as u64;
        }
        (r + (self.multi[w] & ((1u64 << (i % 64)) - 1)).count_ones() as u64) as usize
    }

    /// Span of multi-position bucket `i`, whose first offset is `first`.
    #[inline(always)]
    fn multi_span(&self, i: usize, first: u64) -> BucketSpan {
        let (b, e) = self.rest_bounds.access_pair(self.rank(i));
        BucketSpan::Inline(first, b as usize, e as usize)
    }
}

/// The inline directory for a flat layout whose boundaries are `extras`, when
/// the inline layout is smaller in total (directory, ranks, flags and
/// word-padded offsets); `None` keeps the flat layout.
fn inline_if_smaller(
    extras: &EliasFanoOffsets,
    sizes: impl Iterator<Item = usize>,
    total_offsets: usize,
    width: usize,
    force: bool,
) -> Option<InlineFirst> {
    let num_buckets = extras.len().saturating_sub(1);
    let dir = InlineFirst::directory(num_buckets, sizes);
    let words = |n: usize| (n * width).div_ceil(64) as u64 * 8;
    let inline = dir.directory_bytes() + words(num_buckets) + words(total_offsets - num_buckets);
    let flat = extras.num_bytes() + words(total_offsets);
    info!("  Bucket layout: inline {} bytes, flat {} bytes", inline, flat);
    (force || inline < flat).then_some(dir)
}

/// Where a bucket's super-k-mer offsets are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BucketSpan {
    /// A singleton bucket of the inline layout, with its only offset.
    Single(u64),
    /// Flat layout: `offsets[begin..end]`.
    Flat(usize, usize),
    /// Inline layout: the first offset, then `rest[begin..end]`.
    Inline(u64, usize, usize),
}

impl BucketSpan {
    /// Number of super-k-mers in the bucket.
    #[allow(clippy::len_without_is_empty)] // a bucket is never empty
    #[inline(always)]
    pub fn len(self) -> usize {
        match self {
            BucketSpan::Single(_) => 1,
            BucketSpan::Flat(b, e) => e - b,
            BucketSpan::Inline(_, b, e) => e - b + 1,
        }
    }
}

/// What a lookup learns from one bucket entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Probe {
    /// The whole span.
    Span(BucketSpan),
    /// A light multi-position bucket of the inline layout: its first offset
    /// and id, the rest not decoded ([`SparseAndSkewIndex::rest_span`]).
    First(u64, usize),
}

impl SparseAndSkewIndex {
    /// Create a new empty index
    pub fn new() -> Self {
        Self {
            layout: Layout::Flat {
                extras: EliasFanoOffsets::from_vec(&[0]),
                offsets: BitFieldVec::new(1, 0),
            },
            skew_index: SkewIndex::new(),
        }
    }

    /// Whether the buckets use the inline-first layout (else the flat one).
    pub fn is_inline(&self) -> bool {
        matches!(self.layout, Layout::Inline(_))
    }

    /// Locate a bucket's super-k-mer offsets.
    #[inline(always)]
    pub fn bucket(&self, bucket_id: usize) -> BucketSpan {
        match self.probe(bucket_id) {
            Probe::Span(span) => span,
            Probe::First(first, id) => self.rest_span(id, first),
        }
    }

    /// Read a bucket's entry, leaving the rest of a light multi-position
    /// inline bucket undecoded.
    #[inline(always)]
    pub(crate) fn probe(&self, bucket_id: usize) -> Probe {
        match &self.layout {
            Layout::Flat { extras, .. } => {
                let (b, e) = extras.access_pair(bucket_id);
                Probe::Span(BucketSpan::Flat(b as usize + bucket_id, e as usize + bucket_id + 1))
            }
            Layout::Inline(a) => {
                let first = a.first.index_value(bucket_id) as u64;
                if !a.is_multi(bucket_id) {
                    Probe::Span(BucketSpan::Single(first))
                } else if a.in_heavy_block(bucket_id) {
                    Probe::Span(a.multi_span(bucket_id, first))
                } else {
                    Probe::First(first, bucket_id)
                }
            }
        }
    }

    /// The whole span of a [`Probe::First`] bucket.
    #[inline]
    pub(crate) fn rest_span(&self, bucket_id: usize, first: u64) -> BucketSpan {
        match &self.layout {
            Layout::Inline(a) => a.multi_span(bucket_id, first),
            Layout::Flat { .. } => unreachable!("flat buckets are always decoded whole"),
        }
    }

    /// The `i`-th super-k-mer offset of a bucket (`i < span.len()`).
    #[inline(always)]
    pub fn bucket_offset(&self, span: BucketSpan, i: usize) -> u64 {
        match (span, &self.layout) {
            (BucketSpan::Single(pos), _) => {
                debug_assert_eq!(i, 0);
                pos
            }
            (BucketSpan::Inline(first, _, _), _) if i == 0 => first,
            (BucketSpan::Inline(_, b, _), Layout::Inline(a)) => a.rest.index_value(b + i - 1) as u64,
            (BucketSpan::Flat(b, _), Layout::Flat { offsets, .. }) => offsets.index_value(b + i) as u64,
            _ => unreachable!("span from another layout"),
        }
    }

    /// Convert a flat index (a loaded format 6.0 index, or a build that fell
    /// back to flat) to the inline layout when that is smaller in total.
    /// Returns whether the layout changed. Loading never converts: the
    /// conversion costs a few hundred milliseconds on a transcriptome, so
    /// upgrade once and save the result.
    pub fn upgrade_layout(&mut self) -> bool {
        self.make_inline(false)
    }

    /// Test helper: the same buckets in the inline layout, whatever its size.
    #[cfg(test)]
    pub(crate) fn force_inline(&mut self) {
        self.make_inline(true);
    }

    fn make_inline(&mut self, force: bool) -> bool {
        let Layout::Flat { extras, offsets } = &self.layout else { return false };
        let width = offsets.bit_width();
        let Some(mut a) = inline_if_smaller(extras, flat_sizes(extras), offsets.len(), width, force) else {
            return false;
        };
        let num_buckets = extras.len().saturating_sub(1);
        a.first = BitFieldVec::new(width, num_buckets);
        a.rest = BitFieldVec::new(width, offsets.len() - num_buckets);
        let mut values = offsets.iter();
        let mut dst = 0usize;
        for (i, size) in flat_sizes(extras).enumerate() {
            a.first.set_value(i, values.next().expect("offsets shorter than buckets"));
            for v in values.by_ref().take(size - 1) {
                a.rest.set_value(dst, v);
                dst += 1;
            }
        }
        assert!(values.next().is_none() && dst == a.rest.len(), "bucket sizes do not cover the offsets");
        self.layout = Layout::Inline(a);
        true
    }

    /// Test helper: the same buckets in the flat layout.
    #[cfg(test)]
    pub(crate) fn downgrade_layout(&mut self) {
        let n = self.num_buckets();
        let mut extras = vec![0u64];
        let mut offsets = Vec::new();
        for id in 0..n {
            let span = self.bucket(id);
            extras.push(extras.last().unwrap() + span.len() as u64 - 1);
            offsets.extend((0..span.len()).map(|i| self.bucket_offset(span, i)));
        }
        let width = (64 - offsets.iter().max().copied().unwrap_or(0).leading_zeros() as usize).max(1);
        let mut packed = BitFieldVec::new(width, offsets.len());
        for (i, &v) in offsets.iter().enumerate() {
            packed.set_value(i, v as usize);
        }
        self.layout = Layout::Flat { extras: EliasFanoOffsets::from_vec(&extras), offsets: packed };
    }

    /// Build the sparse and skew index from owned `Bucket` objects (legacy API).
    ///
    /// Prefer [`Self::build_from_classified`] for production use — it avoids
    /// duplicating the tuple payload.
    pub fn build<const K: usize>(
        buckets: &[Bucket],
        mphf_order: Option<&[usize]>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let num_buckets = buckets.len();
        let t0 = std::time::Instant::now();

        let BucketOffsetStarts {
            offset_start_by_orig,
            num_super_kmers_before_bucket,
            total_super_kmers,
            max_bucket_size,
        } = bucket_offset_starts(num_buckets, mphf_order, |b| buckets[b].size());
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size);
        let t1 = std::time::Instant::now();

        let mut offsets = BitFieldVec::new(num_bits_per_offset.max(1), total_super_kmers);
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = Vec::new();
        for (orig_idx, bucket) in buckets.iter().enumerate() {
            let size = bucket.size();
            let idx0 = offset_start_by_orig.index_value(orig_idx);
            // Forward scheme: one tuple per minimizer position (enforced at
            // classification/merge time), so every tuple is one write.
            debug_assert!(bucket.tuples.windows(2).all(|p| p[0].pos_in_seq() != p[1].pos_in_seq()));
            for (off, tuple) in bucket.tuples.iter().enumerate() {
                offsets.set_value(idx0 + off, tuple.pos_in_seq() as usize);
            }
            if size > MIN_BUCKET_SIZE { heavy_buckets_for_skew.push((orig_idx, size)); }
        }
        drop(offset_start_by_orig);
        info!("  Pass 2 (fill offsets): {:.2}s — {} heavy buckets collected",
            t1.elapsed().as_secs_f64(), heavy_buckets_for_skew.len());
        let t2 = std::time::Instant::now();

        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);
        let heavy_refs: Vec<(usize, &Bucket)> = heavy_buckets_for_skew.iter()
            .map(|&(orig_idx, _size)| (orig_idx, &buckets[orig_idx]))
            .collect();
        let skew_index = SkewIndex::build::<K>(&heavy_refs, max_bucket_size, spss);
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Self::flat_upgraded(num_super_kmers_before_bucket, offsets, skew_index)
    }

    /// Build the sparse and skew index from in-place classified buckets.
    ///
    /// This is the production path — avoids tuple duplication.
    /// `mphf_index_by_bucket[b]` is the minimizer MPHF's index of bucket `b`;
    /// the index stores buckets in MPHF order. It is dropped as soon as the
    /// offsets are filled, before the skew index is built.
    pub fn build_from_classified<const K: usize>(
        classified: &ClassifiedBuckets,
        mphf_index_by_bucket: Vec<u32>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        use rayon::prelude::*;
        use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering::Relaxed};
        use sux::bits::bit_field_vec::AtomicBitFieldVec;
        use sux::traits::AtomicBitFieldSlice;

        let num_buckets = classified.num_buckets();
        assert_eq!(mphf_index_by_bucket.len(), num_buckets);
        let t0 = std::time::Instant::now();

        // --- Pass 1: bucket sizes in MPHF order, prefix-summed in place into
        // each bucket's first offset (`starts[num_buckets]` = total), and the
        // EF of super-k-mers before each bucket beyond its first. ---
        let starts: Vec<AtomicU32> = (0..=num_buckets).into_par_iter().map(|_| AtomicU32::new(0)).collect();
        (0..num_buckets).into_par_iter().for_each(|b| {
            starts[mphf_index_by_bucket[b] as usize].store(classified.bucket_size(b) as u32, Relaxed);
        });
        let mut starts: Vec<u32> = starts.into_iter().map(AtomicU32::into_inner).collect();
        let mut total_super_kmers: u64 = 0;
        let mut max_bucket_size: usize = 0;
        for start in starts.iter_mut() {
            let size = *start;
            max_bucket_size = max_bucket_size.max(size as usize);
            *start = total_super_kmers as u32;
            total_super_kmers += size as u64;
        }
        assert!(total_super_kmers < u32::MAX as u64);
        let total_super_kmers = total_super_kmers as usize;
        // extras[i] = sum over buckets j < i of (size_j - 1) = starts[i] - i
        let num_super_kmers_before_bucket = EliasFanoOffsets::from_iter(
            num_buckets + 1,
            (total_super_kmers - num_buckets) as u64,
            starts.iter().enumerate().map(|(i, &s)| s as u64 - i as u64),
        );
        let width = num_bits_per_offset.max(1);
        let inline = inline_if_smaller(
            &num_super_kmers_before_bucket,
            starts.windows(2).map(|w| (w[1] - w[0]) as usize),
            total_super_kmers,
            width,
            false,
        );
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}, {} layout",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size,
            if inline.is_some() { "inline" } else { "flat" });
        let t1 = std::time::Instant::now();

        // --- Pass 2: fill the offsets, buckets in parallel (each writes its
        // own ranges; neighbors may share words, hence atomics). `heads` holds
        // every offset (flat) or each bucket's first one (inline), `rest` the
        // others; bucket i's rest begins at starts[i] - i. ---
        let is_inline = inline.is_some();
        let (heads_len, rest_len) =
            if is_inline { (num_buckets, total_super_kmers - num_buckets) } else { (total_super_kmers, 0) };
        let heads: AtomicBitFieldVec<Vec<AtomicUsize>> = AtomicBitFieldVec::new(width, heads_len);
        let rest: AtomicBitFieldVec<Vec<AtomicUsize>> = AtomicBitFieldVec::new(width, rest_len);
        (0..num_buckets).into_par_iter().with_min_len(1 << 12).for_each(|b| {
            let i = mphf_index_by_bucket[b] as usize;
            let start = starts[i] as usize;
            // Forward scheme: one tuple per minimizer position (enforced at
            // classification time), so every tuple is one write.
            for (off, tuple) in classified.bucket_tuples(b).iter().enumerate() {
                let pos = tuple.pos_in_seq() as usize;
                if !is_inline {
                    heads.set_atomic(start + off, pos, Relaxed);
                } else if off == 0 {
                    heads.set_atomic(i, pos, Relaxed);
                } else {
                    rest.set_atomic(start - i + off - 1, pos, Relaxed);
                }
            }
        });
        drop(starts);
        drop(mphf_index_by_bucket);
        let heads: BitFieldVec = heads.try_into().expect("AtomicUsize and usize have the same alignment");
        let rest: BitFieldVec = rest.try_into().expect("AtomicUsize and usize have the same alignment");
        let layout = match inline {
            Some(mut a) => {
                drop(num_super_kmers_before_bucket);
                a.first = heads;
                a.rest = rest;
                Layout::Inline(a)
            }
            None => Layout::Flat { extras: num_super_kmers_before_bucket, offsets: heads },
        };
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = (0..num_buckets)
            .into_par_iter()
            .filter_map(|b| {
                let size = classified.bucket_size(b);
                (size > MIN_BUCKET_SIZE).then_some((b, size))
            })
            .collect();

        info!("  Pass 2 (fill offsets): {:.2}s — {} heavy buckets collected",
            t1.elapsed().as_secs_f64(), heavy_buckets_for_skew.len());
        let t2 = std::time::Instant::now();

        // Sort heavy buckets by size for partition processing
        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);

        // --- Pass 3: Build skew index for heavy buckets ---
        let skew_index = SkewIndex::build_from_classified::<K>(
            classified,
            &heavy_buckets_for_skew,
            max_bucket_size,
            spss,
        );
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Self { layout, skew_index }
    }

    /// Build the sparse and skew index from file tuples (streaming path).
    ///
    /// This avoids materializing all tuples in memory. Instead it:
    /// - Pass 1: Uses `cached_sizes` + `mphf_order` to build EF and compute
    ///   `offset_start_by_orig`. Then drops `mphf_order` (saves ~3.2 GB).
    /// - Pass 2: Reads file sequentially via `SequentialTupleReader` (no mmap).
    ///   Fills offsets and collects heavy bucket tuples during the forward scan.
    /// - Pass 3: Builds skew index from materialized heavy buckets.
    pub fn build_from_file<const K: usize>(
        file_tuples: &FileTuples,
        bucket_meta: &BucketMetadata,
        mphf_order: Option<Vec<usize>>,
        num_bits_per_offset: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Result<Self, std::io::Error>
    where
        Kmer<K>: KmerBits,
    {
        let num_buckets = bucket_meta.num_buckets();
        let t0 = std::time::Instant::now();

        // --- Pass 1: Build EF + compute per-bucket offset positions ---
        // Then drop mphf_order to free ~3.2 GB.
        let BucketOffsetStarts {
            offset_start_by_orig,
            num_super_kmers_before_bucket,
            total_super_kmers,
            max_bucket_size,
        } = bucket_offset_starts(num_buckets, mphf_order.as_deref(), |b| {
            bucket_meta.cached_sizes[b] as usize
        });
        // Free mphf_order — no longer needed after offset_start_by_orig is computed
        drop(mphf_order);
        info!("  Pass 1 (EF build): {:.2}s — {} buckets, {} total super-kmers, max_bucket_size={}",
            t0.elapsed().as_secs_f64(), num_buckets, total_super_kmers, max_bucket_size);
        let t1 = std::time::Instant::now();

        // --- Pass 2: Sequential file scan → fill offsets + collect heavy tuples ---
        // Uses a SequentialTupleReader (BufReader) — no mmap pages kept resident.
        let mut offsets = BitFieldVec::new(num_bits_per_offset.max(1), total_super_kmers);
        let mut heavy_buckets_for_skew: Vec<(usize, usize)> = Vec::new();
        let mut heavy_tuples: Vec<Vec<MinimizerTuple>> = Vec::new();
        let mut bucket_idx = 0usize;

        let mut reader = file_tuples.sequential_reader()?;
        let mut next_tuple = reader.read_next()?;

        while let Some(first) = next_tuple {
            let minimizer = first.minimizer;
            let size = bucket_meta.cached_sizes[bucket_idx] as usize;
            let is_heavy = size > MIN_BUCKET_SIZE;

            // Forward scheme: one tuple per minimizer position (enforced at
            // merge time), so every tuple is one offset write.
            let mut write_idx = offset_start_by_orig.index_value(bucket_idx);
            let mut prev_pos_in_seq;

            if is_heavy {
                // Collect heavy bucket tuples during the forward scan
                let mut bucket_tuples = Vec::new();

                offsets.set_value(write_idx, first.pos_in_seq as usize);
                write_idx += 1;
                prev_pos_in_seq = first.pos_in_seq;
                bucket_tuples.push(first.to_internal());

                // Read remaining tuples in this bucket
                next_tuple = reader.read_next()?;
                while let Some(t) = next_tuple {
                    if t.minimizer != minimizer { break; }
                    let t_pos = t.pos_in_seq;
                    debug_assert_ne!(t_pos, prev_pos_in_seq);
                    offsets.set_value(write_idx, t_pos as usize);
                    write_idx += 1;
                    prev_pos_in_seq = t_pos;
                    bucket_tuples.push(t.to_internal());
                    next_tuple = reader.read_next()?;
                }

                heavy_buckets_for_skew.push((heavy_tuples.len(), size));
                heavy_tuples.push(bucket_tuples);
            } else {
                // Stream non-heavy bucket tuples without materialization
                offsets.set_value(write_idx, first.pos_in_seq as usize);
                write_idx += 1;
                prev_pos_in_seq = first.pos_in_seq;

                next_tuple = reader.read_next()?;
                while let Some(ext) = next_tuple {
                    if ext.minimizer != minimizer { break; }
                    let ext_pos = ext.pos_in_seq;
                    debug_assert_ne!(ext_pos, prev_pos_in_seq);
                    offsets.set_value(write_idx, ext_pos as usize);
                    write_idx += 1;
                    prev_pos_in_seq = ext_pos;
                    next_tuple = reader.read_next()?;
                }
            }
            let _ = prev_pos_in_seq;

            bucket_idx += 1;
        }
        drop(offset_start_by_orig);

        info!("  Pass 2 (fill offsets from file): {:.2}s — {} heavy buckets materialized",
            t1.elapsed().as_secs_f64(), heavy_tuples.len());
        let t2 = std::time::Instant::now();

        // Sort heavy buckets by size for partition processing
        heavy_buckets_for_skew.sort_by_key(|&(_idx, size)| size);

        // --- Pass 3: Build skew index from materialized heavy bucket tuples ---
        let tuples_iter = heavy_buckets_for_skew.iter()
            .map(|&(heavy_idx, _size)| {
                let tuples = &heavy_tuples[heavy_idx];
                // Forward scheme: one tuple per position, so size == count.
                (tuples.len(), tuples.as_slice())
            });
        let skew_index = SkewIndex::build_inner::<K>(tuples_iter, max_bucket_size, spss);
        info!("  Pass 3 (skew index): {:.2}s", t2.elapsed().as_secs_f64());

        Ok(Self::flat_upgraded(num_super_kmers_before_bucket, offsets, skew_index))
    }

    /// A flat index, upgraded to the inline layout when that is smaller.
    // shortcut: the legacy and external builds lay out flat, then convert, so
    // they briefly hold both copies of the offsets; build inline directly
    // there too if their peak memory matters.
    fn flat_upgraded(extras: EliasFanoOffsets, offsets: BitFieldVec, skew_index: SkewIndex) -> Self {
        let mut index = Self { layout: Layout::Flat { extras, offsets }, skew_index };
        index.upgrade_layout();
        index
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        let offsets_bits = match &self.layout {
            Layout::Flat { offsets, .. } => offsets.len() * offsets.bit_width(),
            Layout::Inline(a) => a.first.len() * a.first.bit_width() + a.rest.len() * a.rest.bit_width(),
        } as u64;
        self.ef_bytes() as u64 * 8 + offsets_bits + self.skew_index.num_bits()
    }

    /// Byte size of the bucket directory: the Elias-Fano bounds, plus, for the
    /// inline layout, the multi-position bitmap, its ranks and the heavy flags
    pub fn ef_bytes(&self) -> usize {
        match &self.layout {
            Layout::Flat { extras, .. } => extras.num_bytes() as usize,
            Layout::Inline(a) => a.directory_bytes() as usize,
        }
    }

    /// Byte size of the super-k-mer offsets (BitFieldVec backing data)
    pub fn offsets_bytes(&self) -> usize {
        match &self.layout {
            Layout::Flat { offsets, .. } => std::mem::size_of_val(offsets.as_slice()),
            Layout::Inline(a) => {
                std::mem::size_of_val(a.first.as_slice()) + std::mem::size_of_val(a.rest.as_slice())
            }
        }
    }

    /// Byte size of skew index data (positions only, excluding MPHFs)
    pub fn skew_index_bytes(&self) -> usize {
        self.skew_index.data_bytes()
    }

    /// Byte size of skew index MPHFs (serialized)
    pub fn skew_mphf_bytes(&self) -> usize {
        self.skew_index.mphf_bytes()
    }

    /// Number of buckets stored
    pub fn num_buckets(&self) -> usize {
        match &self.layout {
            Layout::Flat { extras, .. } => extras.len().saturating_sub(1),
            Layout::Inline(a) => a.first.len(),
        }
    }

    /// Total number of super-kmers (offsets entries)
    pub fn num_offsets(&self) -> usize {
        match &self.layout {
            Layout::Flat { offsets, .. } => offsets.len(),
            Layout::Inline(a) => a.first.len() + a.rest.len(),
        }
    }
}

impl SparseAndSkewIndex {
    /// Serialize without MPHF (for splitting MPHF across separate container)
    pub fn serialize_without_mphf<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> std::io::Result<()> {
        use crate::darray::{write_u64, write_vec_u64};
        use epserde::ser::Serialize;

        // SAFETY (each `serialize`): BitFieldVec serialization preserves data layout
        match &self.layout {
            Layout::Flat { extras, offsets } => {
                write_u64(writer, LAYOUT_FLAT)?;
                extras.write_to(writer)?;
                unsafe { offsets.serialize(writer).map_err(std::io::Error::other)? };
            }
            Layout::Inline(a) => {
                write_u64(writer, LAYOUT_INLINE)?;
                unsafe { a.first.serialize(writer).map_err(std::io::Error::other)? };
                write_vec_u64(writer, &a.multi)?;
                write_vec_u64(writer, &a.heavy)?;
                a.rest_bounds.write_to(writer)?;
                unsafe { a.rest.serialize(writer).map_err(std::io::Error::other)? };
            }
        }

        // Serialize skew index (excluding MPHF)
        self.skew_index.serialize_without_mphf(writer)?;

        Ok(())
    }

    /// Deserialize without MPHF (MPHFs loaded separately from container).
    /// `format_major` is the index header's: 6 has no layout tag and is flat.
    /// Either layout is kept as stored; see [`Self::upgrade_layout`].
    pub fn deserialize_without_mphf<R: std::io::Read>(
        reader: &mut R,
        format_major: u32,
    ) -> std::io::Result<Self> {
        use crate::darray::{read_u64, read_vec_u64};
        use epserde::deser::Deserialize;

        let invalid = |what: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("sparse index: {what}"));
        let tag = if format_major == 6 { LAYOUT_FLAT } else { read_u64(reader)? };
        // SAFETY (each `deserialize_full`): reading data serialized by our own serialize method
        let layout = match tag {
            LAYOUT_FLAT => {
                let extras = EliasFanoOffsets::read_from(reader)?;
                let offsets = unsafe { BitFieldVec::deserialize_full(reader) }.map_err(std::io::Error::other)?;
                let n = extras.len().saturating_sub(1);
                if extras.is_empty() || extras.access(n) as usize + n != offsets.len() {
                    return Err(invalid("bucket bounds do not cover the offsets"));
                }
                Layout::Flat { extras, offsets }
            }
            LAYOUT_INLINE => {
                let first = unsafe { BitFieldVec::deserialize_full(reader) }.map_err(std::io::Error::other)?;
                let multi = read_vec_u64(reader)?;
                let heavy = read_vec_u64(reader)?;
                let rest_bounds = EliasFanoOffsets::read_from(reader)?;
                let rest = unsafe { BitFieldVec::deserialize_full(reader) }.map_err(std::io::Error::other)?;
                let n = first.len();
                let num_multi: u64 = multi.iter().map(|w| w.count_ones() as u64).sum();
                if multi.len() != n.div_ceil(64)
                    || heavy.len() != n.div_ceil(64 * HEAVY_BLOCK)
                    || rest_bounds.len() as u64 != num_multi + 1
                    || rest_bounds.access(num_multi as usize) as usize != rest.len()
                {
                    return Err(invalid("inline directory does not match the offsets"));
                }
                Layout::Inline(InlineFirst { first, ranks: sample_ranks(&multi), multi, heavy, rest_bounds, rest })
            }
            _ => return Err(invalid("unknown bucket layout")),
        };
        let skew_index = SkewIndex::deserialize_without_mphf(reader)?;
        Ok(Self { layout, skew_index })
    }
}

impl Default for SparseAndSkewIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Skew index for heavy buckets (size > MIN_BUCKET_SIZE)
///
/// Heavy buckets are partitioned by size into logarithmic partitions:
/// - Partition 0: buckets of size (64, 128]
/// - Partition 1: buckets of size (128, 256]
/// - Partition 2: buckets of size (256, 512]
/// - etc.
///
/// For each partition, we build a separate MPHF over all k-mers in that partition.
/// `positions[partition][mphf_idx]` stores the within-bucket super-kmer index.
/// The caller resolves to an absolute position via `offsets[begin + within_pos]`.
pub struct SkewIndex {
    /// Vector of MPHFs, one per partition
    /// Each MPHF maps k-mers to positions within the partition
    pub mphfs: Vec<Option<PartitionedMphf>>,

    /// Vector of position arrays, one per partition.
    /// `positions[i]` stores the within-bucket super-kmer index for each k-mer.
    pub positions: Vec<BitFieldVec>,
}

impl SkewIndex {
    /// Serialize without MPHF (for splitting MPHF across separate container)
    pub fn serialize_without_mphf<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> std::io::Result<()> {
        use epserde::ser::Serialize;

        // Serialize number of partitions
        let num_partitions = self.mphfs.len() as u32;
        writer.write_all(&num_partitions.to_le_bytes())?;

        // Serialize positions (one BitFieldVec per partition)
        writer.write_all(&(self.positions.len() as u32).to_le_bytes())?;
        for pos_vec in &self.positions {
            // SAFETY: BitFieldVec serialization preserves data layout
            unsafe { pos_vec
                .serialize(writer)
                .map_err(std::io::Error::other)? };
        }

        Ok(())
    }

    /// Deserialize without MPHF (MPHFs loaded separately)
    pub fn deserialize_without_mphf<R: std::io::Read>(
        reader: &mut R,
    ) -> std::io::Result<Self> {
        use epserde::deser::Deserialize;

        // Deserialize number of partitions
        let mut num_parts_bytes = [0u8; 4];
        reader.read_exact(&mut num_parts_bytes)?;
        let num_partitions = u32::from_le_bytes(num_parts_bytes) as usize;

        // Initialize with empty MPHFs (will be filled later)
        let mphfs = (0..num_partitions).map(|_| None).collect();

        // Deserialize positions
        let mut pos_bytes = [0u8; 4];
        reader.read_exact(&mut pos_bytes)?;
        let num_pos_vecs = u32::from_le_bytes(pos_bytes) as usize;

        let mut positions = Vec::with_capacity(num_pos_vecs);
        for _ in 0..num_pos_vecs {
            // SAFETY: Reading data serialized by our own serialize method
            let pos_vec = unsafe { BitFieldVec::deserialize_full(reader) }
                .map_err(std::io::Error::other)?;
            positions.push(pos_vec);
        }

        Ok(Self {
            mphfs,
            positions,
        })
    }

    /// Build the skew index from heavy buckets (legacy API for tests).
    ///
    /// See [`Self::build_from_classified`] for the production path.
    pub fn build<const K: usize>(
        heavy_buckets: &[(usize, &Bucket)],
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        // Adapt: extract tuple slices from the owned Bucket refs
        let tuples_iter = heavy_buckets.iter()
            .map(|&(_orig_idx, bucket)| (bucket.size(), bucket.tuples.as_slice()));
        Self::build_inner::<K>(tuples_iter, max_bucket_size, spss)
    }

    /// Build the skew index from [`ClassifiedBuckets`] heavy bucket refs.
    pub fn build_from_classified<const K: usize>(
        classified: &ClassifiedBuckets,
        heavy_buckets_for_skew: &[(usize, usize)], // (orig_idx, size), sorted by size
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let tuples_iter = heavy_buckets_for_skew.iter()
            .map(|&(orig_idx, _size)| {
                (classified.bucket_size(orig_idx), classified.bucket_tuples(orig_idx))
            });
        Self::build_inner::<K>(tuples_iter, max_bucket_size, spss)
    }

    /// Core skew index builder.
    ///
    /// `heavy_buckets` yields `(bucket_size, tuple_slice)` pairs sorted by size.
    fn build_inner<'a, const K: usize>(
        heavy_buckets: impl Iterator<Item = (usize, &'a [MinimizerTuple])>,
        max_bucket_size: usize,
        spss: &crate::spectrum_preserving_string_set::SpectrumPreservingStringSet,
    ) -> Self
    where
        Kmer<K>: KmerBits,
    {
        let min_size = MIN_BUCKET_SIZE;

        // Calculate number of partitions
        let log2_max_bucket_size = ceil_log2(max_bucket_size as u64);
        let num_partitions = if max_bucket_size <= min_size {
            0
        } else if max_bucket_size < (1 << MAX_L) {
            log2_max_bucket_size - MIN_L
        } else {
            MAX_L - MIN_L + 1
        };

        if num_partitions == 0 {
            return Self::new();
        }

        info!("    SkewIndex: max_bucket_size={}, num_partitions={}",
            max_bucket_size, num_partitions);

        let mut mphfs = Vec::with_capacity(num_partitions);
        let mut positions = Vec::with_capacity(num_partitions);

        let mut partition_id = 0;
        let mut lower = min_size;
        let mut upper = 2 * lower;
        let mut num_bits_per_pos = MIN_L + 1;

        let mut partition_kmers: Vec<<Kmer<K> as KmerBits>::Storage> = Vec::new();
        // Positions within a bucket (max bucket size < 2^32)
        let mut partition_positions: Vec<u32> = Vec::new();

        for (bucket_size, bucket_tuples) in heavy_buckets {
            // Finalize current partition(s) if this bucket is too large
            while bucket_size > upper {
                if !partition_kmers.is_empty() {
                    let mphf_opt = build_partition_mphf(&partition_kmers);
                    let reordered_positions = if let Some(ref mphf) = mphf_opt {
                        let mut rp = BitFieldVec::new(num_bits_per_pos, partition_positions.len());
                        for (i, &kmer_value) in partition_kmers.iter().enumerate() {
                            let pos = mphf.get(&kmer_value);
                            rp.set_value(pos, partition_positions[i] as usize);
                        }
                        rp
                    } else {
                        BitFieldVec::new(num_bits_per_pos, 0)
                    };
                    mphfs.push(mphf_opt);
                    positions.push(reordered_positions);
                    partition_kmers.clear();
                    partition_positions.clear();
                } else {
                    // Lookup addresses partitions by size band, including empty
                    // bands. Keep their slots so later bands cannot shift down.
                    mphfs.push(None);
                    positions.push(BitFieldVec::new(num_bits_per_pos, 0));
                }

                lower = upper;
                upper = 2 * lower;
                num_bits_per_pos += 1;
                partition_id += 1;

                if partition_id == num_partitions - 1 {
                    upper = max_bucket_size;
                    num_bits_per_pos = log2_max_bucket_size;
                }
            }

            // Add k-mers from this bucket to current partition. The scheme is
            // forward: one tuple per super-kmer, so pos_in_bucket increments
            // per tuple with no dedupe.
            let mut pos_in_bucket = 0u32;

            for tuple in bucket_tuples {
                pos_in_bucket += 1;

                debug_assert!(tuple.pos_in_seq() >= tuple.pos_in_kmer() as u64);
                let starting_kmer_pos =
                    (tuple.pos_in_seq() - tuple.pos_in_kmer() as u64) as usize;

                for kmer_offset in 0..tuple.num_kmers_in_super_kmer() {
                    let kmer: Kmer<K> = spss.decode_kmer_at(starting_kmer_pos + kmer_offset as usize);
                    // The skew MPHF keys on the canonical k-mer, matching the
                    // canonical probe in the lookup path.
                    let rc_value = kmer.reverse_complement().bits();
                    let kmer_value = kmer.bits().min(rc_value);

                    partition_kmers.push(kmer_value);
                    partition_positions.push(pos_in_bucket - 1);
                }
            }
        }

        // Finalize last partition
        if !partition_kmers.is_empty() {
            let mphf_opt = build_partition_mphf(&partition_kmers);
            let reordered_positions = if let Some(ref mphf) = mphf_opt {
                let mut rp = BitFieldVec::new(num_bits_per_pos, partition_positions.len());
                for (i, &kmer_value) in partition_kmers.iter().enumerate() {
                    let pos = mphf.get(&kmer_value);
                    rp.set_value(pos, partition_positions[i] as usize);
                }
                rp
            } else {
                BitFieldVec::new(num_bits_per_pos, 0)
            };
            mphfs.push(mphf_opt);
            positions.push(reordered_positions);
        }

        Self { mphfs, positions }
    }

    /// Lookup a k-mer in the skew index.
    ///
    /// Returns the within-bucket super-kmer index. The caller translates this
    /// to an absolute position via `offsets[begin + result]`.
    ///
    /// # Arguments
    /// * `kmer_value` - The k-mer value (same type used to build the MPHF)
    /// * `log2_bucket_size` - ceil_log2 of the bucket size (determines partition)
    ///
    /// # Returns
    /// The within-bucket super-kmer index, or INVALID_UINT64 if not found
    pub fn lookup<T: MphfKey>(&self, kmer_value: &T, log2_bucket_size: usize) -> u64 {
        if log2_bucket_size <= MIN_L {
            return INVALID_UINT64;
        }

        let partition_id = (log2_bucket_size - (MIN_L + 1))
            .min(self.mphfs.len().saturating_sub(1));

        if partition_id >= self.mphfs.len() {
            return INVALID_UINT64;
        }

        let mphf = match &self.mphfs[partition_id] {
            Some(m) => m,
            None => return INVALID_UINT64,
        };
        let pos = mphf.get(kmer_value);

        if pos >= self.positions[partition_id].len() {
            return INVALID_UINT64;
        }

        self.positions[partition_id].index_value(pos) as u64
    }

    /// Get memory usage in bits
    pub fn num_bits(&self) -> u64 {
        // MPHF memory is estimated at ~4 bits per key
        let total_keys: usize = self.positions.iter().map(|p| p.len()).sum();
        let mphf_bits = (total_keys as f64 * 4.0) as u64;

        let pos_bits: u64 = self.positions.iter()
            .map(|p| (p.len() * p.bit_width()) as u64)
            .sum();

        mphf_bits + pos_bits
    }

    /// Byte size of data (positions only), excluding MPHFs
    pub fn data_bytes(&self) -> usize {
        self.positions.iter()
            .map(|p| std::mem::size_of_val(p.as_slice()))
            .sum()
    }

    /// Byte size of serialized MPHFs
    pub fn mphf_bytes(&self) -> usize {
        self.mphfs.iter()
            .filter_map(|opt| opt.as_ref())
            .map(|pmphf| pmphf.write_bytes())
            .sum()
    }

    /// Create a new empty skew index
    pub fn new() -> Self {
        Self {
            mphfs: Vec::new(),
            positions: Vec::new(),
        }
    }

    /// Get the number of MPHF partitions
    pub fn num_partitions(&self) -> usize {
        self.mphfs.len()
    }

    /// Get a reference to the MPHFs vector
    pub fn mphfs_ref(&self) -> &[Option<PartitionedMphf>] {
        &self.mphfs
    }

    /// Set the MPHFs vector (used during deserialization)
    pub fn set_mphfs(&mut self, mphfs: Vec<Option<PartitionedMphf>>) {
        self.mphfs = mphfs;
    }
}

/// Build a PartitionedMphf for a skew index partition's k-mers.
///
/// Generic over the key type (`u64` for K <= 31, `u128` for K > 31).
/// Skew index partitions are typically small enough that a single (non-partitioned)
/// MPHF is used — partitioning is only relevant for the main minimizer MPHF.
fn build_partition_mphf<T: MphfKey>(kmers: &[T]) -> Option<PartitionedMphf> {
    if kmers.is_empty() {
        return None;
    }

    // Skew index MPHFs are always monolithic (single partition) since they cover
    // a single size-class bucket, not the full minimizer set.
    Some(PartitionedMphf::build_single(kmers))
}

impl Default for SkewIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Bucket layout in the flat offsets array, buckets laid out in MPHF order.
struct BucketOffsetStarts {
    /// Each bucket's first offset, indexed by original bucket, packed at
    /// the width of the largest start (including starts above u32::MAX).
    offset_start_by_orig: BitFieldVec,
    /// Super-k-mers before each bucket, beyond each bucket's first
    num_super_kmers_before_bucket: EliasFanoOffsets,
    total_super_kmers: usize,
    max_bucket_size: usize,
}

/// Lay out buckets `mphf_order[0], mphf_order[1], ...` (identity if `None`)
/// consecutively; `size_of(b)` is original bucket `b`'s size.
fn bucket_offset_starts(
    num_buckets: usize,
    mphf_order: Option<&[usize]>,
    size_of: impl Fn(usize) -> usize,
) -> BucketOffsetStarts {
    let mut extras = Vec::with_capacity(num_buckets + 1);
    extras.push(0u64);
    let mut total_super_kmers: usize = 0;
    let mut max_bucket_size: usize = 0;
    for i in 0..num_buckets {
        let orig_idx = mphf_order.map_or(i, |order| order[i]);
        let size = size_of(orig_idx);
        max_bucket_size = max_bucket_size.max(size);
        total_super_kmers += size;
        let prev = *extras.last().unwrap();
        extras.push(prev + (size as u64).saturating_sub(1));
    }
    // Every nonempty bucket starts below total_super_kmers. Reuse the
    // uncompressed EF prefix sums to scatter starts, without storing a
    // temporary full-width starts vector or re-reading bucket sizes.
    let max_start = if num_buckets == 0 {
        0
    } else {
        extras[num_buckets - 1] as usize + num_buckets - 1
    };
    let width = (usize::BITS - max_start.leading_zeros()).max(1) as usize;
    let mut offset_start_by_orig = BitFieldVec::new(width, num_buckets);
    for (i, &extra) in extras[..num_buckets].iter().enumerate() {
        let orig_idx = mphf_order.map_or(i, |order| order[i]);
        offset_start_by_orig.set_value(orig_idx, extra as usize + i);
    }
    BucketOffsetStarts {
        offset_start_by_orig,
        num_super_kmers_before_bucket: EliasFanoOffsets::from_vec(&extras),
        total_super_kmers,
        max_bucket_size,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::minimizer_tuples::MinimizerTuple;

    #[test]
    fn test_packed_start_widths() {
        for sizes in [vec![], vec![1], vec![8, 1], vec![7, 1], vec![1, 1000]] {
            let n = sizes.len();
            let order: Vec<usize> = (0..n).rev().collect();
            let starts = bucket_offset_starts(n, Some(&order), |b| sizes[b]);
            let mut total = 0;
            let mut max_start = 0usize;
            for &b in &order {
                assert_eq!(starts.offset_start_by_orig.index_value(b), total);
                max_start = total;
                total += sizes[b];
            }
            assert_eq!(starts.total_super_kmers, total);
            assert_eq!(
                starts.offset_start_by_orig.bit_width(),
                (usize::BITS - max_start.leading_zeros()).max(1) as usize
            );
        }
    }

    #[test]
    fn test_bucket_offset_starts_past_u32() {
        // Sizes summing past u32::MAX, permuted so later buckets would overwrite
        // earlier slots if starts wrapped
        let sizes = [3usize << 31, 3 << 31, 1, 2];
        let order = [2usize, 0, 3, 1];
        let starts = bucket_offset_starts(sizes.len(), Some(&order), |b| sizes[b]);
        assert_eq!(
            starts.offset_start_by_orig.iter().collect::<Vec<_>>(),
            vec![1, (3 << 31) + 3, 0, (3 << 31) + 1]
        );
        assert_eq!(starts.total_super_kmers, (3 << 32) + 3);
        assert_eq!(starts.max_bucket_size, 3 << 31);
        for (i, &b) in order.iter().enumerate() {
            let begin = starts.num_super_kmers_before_bucket.access(i) as usize + i;
            assert_eq!(begin, starts.offset_start_by_orig.index_value(b));
        }
    }

    #[test]
    fn test_sparse_index_creation() {
        let index = SparseAndSkewIndex::new();
        assert_eq!(index.num_buckets(), 0);
        assert_eq!(index.num_offsets(), 0);
    }

    /// Every bucket's offsets, through the public API.
    fn contents(index: &SparseAndSkewIndex) -> Vec<Vec<u64>> {
        (0..index.num_buckets())
            .map(|b| {
                let span = index.bucket(b);
                (0..span.len()).map(|i| index.bucket_offset(span, i)).collect()
            })
            .collect()
    }

    /// The same contents in both layouts, and through a 7.0 round trip.
    fn check_layouts(mut index: SparseAndSkewIndex, expected: &[Vec<u64>]) {
        index.force_inline();
        assert!(index.is_inline());
        assert_eq!(contents(&index), expected);
        let mut bytes = Vec::new();
        index.serialize_without_mphf(&mut bytes).unwrap();
        let back = SparseAndSkewIndex::deserialize_without_mphf(&mut bytes.as_slice(), 7).unwrap();
        assert!(back.is_inline());
        assert_eq!(contents(&back), expected);
        index.downgrade_layout();
        assert!(!index.is_inline());
        assert_eq!(contents(&index), expected);
        assert_eq!(index.num_offsets(), expected.iter().map(Vec::len).sum::<usize>());
    }

    fn tiny_spss() -> crate::spectrum_preserving_string_set::SpectrumPreservingStringSet {
        let mut encoder = crate::builder::encode::Encoder::<31>::new();
        encoder.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACG").unwrap();
        encoder.build(13)
    }

    #[test]
    fn test_sparse_index_build_singleton() {
        let buckets = vec![Bucket::new(100, vec![MinimizerTuple::new(100, 50, 0)])];
        let index = SparseAndSkewIndex::build::<31>(&buckets, None, 32, &tiny_spss());
        assert_eq!(index.num_buckets(), 1);
        assert_eq!(contents(&index), [vec![50]]);
        check_layouts(index, &[vec![50]]);
    }

    #[test]
    fn test_sparse_index_build_light() {
        let bucket = Bucket::new(200, vec![
            MinimizerTuple::new(200, 100, 0),
            MinimizerTuple::new(200, 101, 0),
            MinimizerTuple::new(200, 102, 0),
        ]);
        let index = SparseAndSkewIndex::build::<31>(&[bucket], None, 32, &tiny_spss());
        assert_eq!(index.num_buckets(), 1);
        assert_eq!(index.num_offsets(), 3);
        assert_eq!(contents(&index), [vec![100, 101, 102]]);
        check_layouts(index, &[vec![100, 101, 102]]);
    }

    #[test]
    fn test_sparse_index_build_mixed() {
        let buckets = vec![
            Bucket::new(100, vec![MinimizerTuple::new(100, 50, 0)]), // Singleton
            Bucket::new(200, vec![
                MinimizerTuple::new(200, 100, 0),
                MinimizerTuple::new(200, 101, 0),
            ]), // Light (size=2)
            Bucket::new(300, vec![MinimizerTuple::new(300, 200, 0)]), // Singleton
        ];
        let index = SparseAndSkewIndex::build::<31>(&buckets, None, 32, &tiny_spss());
        let expected = [vec![50], vec![100, 101], vec![200]];
        assert_eq!(index.num_buckets(), 3);
        assert_eq!(contents(&index), expected);
        check_layouts(index, &expected);
    }

    /// Mostly singletons pick the inline layout; all-double buckets, where
    /// the bitmap and the second Elias-Fano cost more than they save, fall
    /// back to flat.
    #[test]
    fn test_layout_follows_storage_cost() {
        let spss = tiny_spss();
        let make = |sizes: &[usize]| -> Vec<Bucket> {
            let mut pos = 0u64;
            sizes
                .iter()
                .enumerate()
                .map(|(b, &size)| {
                    let tuples = (0..size)
                        .map(|_| {
                            pos += 3;
                            MinimizerTuple::new(b as u64, pos, 0)
                        })
                        .collect();
                    Bucket::new(b as u64, tuples)
                })
                .collect()
        };
        let mostly_single: Vec<usize> = (0..20_000).map(|i| if i % 20 == 0 { 3 } else { 1 }).collect();
        let index = SparseAndSkewIndex::build::<31>(&make(&mostly_single), None, 20, &spss);
        assert!(index.is_inline());
        let all_double = vec![2usize; 20_000];
        let index = SparseAndSkewIndex::build::<31>(&make(&all_double), None, 20, &spss);
        assert!(!index.is_inline());
        let expected = contents(&index);
        check_layouts(index, &expected);
    }

    /// The probe decodes a multi-position bucket's bounds at once only in a
    /// heavy-flagged block, including the light buckets that share it; ranks
    /// stay right across sample boundaries.
    #[test]
    fn test_inline_probe_and_heavy_blocks() {
        let n = 3 * 512 + 7;
        // Bucket 5 is heavy (block 0); bucket 9 is light in that block;
        // bucket 600 is light in a block without heavy buckets.
        let size = |i: usize| match i {
            5 => MIN_BUCKET_SIZE + 6,
            9 | 600 => 3,
            i if i % 37 == 0 => 2,
            _ => 1,
        };
        let mut a = InlineFirst::directory(n, (0..n).map(size));
        let total: usize = (0..n).map(size).sum();
        a.first = BitFieldVec::new(24, n);
        a.rest = BitFieldVec::new(24, total - n);
        let mut dst = 0;
        for i in 0..n {
            a.first.set_value(i, 1000 * i);
            for j in 1..size(i) {
                a.rest.set_value(dst, 1000 * i + j);
                dst += 1;
            }
        }
        let index = SparseAndSkewIndex { layout: Layout::Inline(a), skew_index: SkewIndex::new() };
        assert!(matches!(index.probe(5), Probe::Span(BucketSpan::Inline(5000, _, _))));
        assert!(matches!(index.probe(9), Probe::Span(BucketSpan::Inline(9000, _, _))));
        assert_eq!(index.probe(600), Probe::First(600_000, 600));
        assert_eq!(index.probe(1), Probe::Span(BucketSpan::Single(1000)));
        for i in 0..n {
            let span = index.bucket(i);
            assert_eq!(span.len(), size(i));
            for j in 0..size(i) {
                assert_eq!(index.bucket_offset(span, j), (1000 * i + j) as u64);
            }
        }
    }

    /// A format 6.0 index (flat layout) loads as stored, without conversion;
    /// `upgrade_layout` converts it explicitly when inline is smaller.
    #[test]
    fn test_reads_format_6_flat_layout() {
        use epserde::ser::Serialize;
        let sizes: Vec<usize> = (0..5000).map(|b| [1, 2, 1, 3, 1, 1, 1, 1, 1, 1][b % 10]).collect();
        let mut extras = vec![0u64];
        let mut offsets_vec = Vec::new();
        let mut expected = Vec::new();
        for (b, &size) in sizes.iter().enumerate() {
            extras.push(extras.last().unwrap() + size as u64 - 1);
            let bucket: Vec<u64> = (0..size).map(|i| (100 * b + 7 * i) as u64).collect();
            offsets_vec.extend(bucket.iter().copied());
            expected.push(bucket);
        }
        let mut offsets = BitFieldVec::new(20, offsets_vec.len());
        for (i, &v) in offsets_vec.iter().enumerate() {
            offsets.set_value(i, v as usize);
        }
        let mut bytes = Vec::new();
        EliasFanoOffsets::from_vec(&extras).write_to(&mut bytes).unwrap();
        unsafe { offsets.serialize(&mut bytes).unwrap() };
        SkewIndex::new().serialize_without_mphf(&mut bytes).unwrap();

        let mut index = SparseAndSkewIndex::deserialize_without_mphf(&mut bytes.as_slice(), 6).unwrap();
        assert!(!index.is_inline());
        assert_eq!(index.num_buckets(), sizes.len());
        assert_eq!(index.num_offsets(), offsets_vec.len());
        assert_eq!(contents(&index), expected);

        // A 7.0 round trip keeps it flat ...
        let mut bytes7 = Vec::new();
        index.serialize_without_mphf(&mut bytes7).unwrap();
        let back = SparseAndSkewIndex::deserialize_without_mphf(&mut bytes7.as_slice(), 7).unwrap();
        assert!(!back.is_inline());
        assert_eq!(contents(&back), expected);

        // ... until an explicit upgrade, which is then stored.
        let flat_bytes = index.ef_bytes() + index.offsets_bytes();
        assert!(index.upgrade_layout());
        assert!(index.is_inline());
        assert!(index.ef_bytes() + index.offsets_bytes() < flat_bytes);
        assert!(!index.upgrade_layout());
        assert_eq!(contents(&index), expected);
        let mut bytes7 = Vec::new();
        index.serialize_without_mphf(&mut bytes7).unwrap();
        let back = SparseAndSkewIndex::deserialize_without_mphf(&mut bytes7.as_slice(), 7).unwrap();
        assert!(back.is_inline());
        assert_eq!(contents(&back), expected);

        // A corrupted layout tag is refused.
        bytes7[0] = 9;
        assert!(SparseAndSkewIndex::deserialize_without_mphf(&mut bytes7.as_slice(), 7).is_err());
    }
}
