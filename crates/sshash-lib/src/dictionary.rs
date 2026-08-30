//! Dictionary - the main SSHash data structure
//!
//! Provides efficient storage and querying of k-mer sets using:
//! - Spectrum-preserving string encoding
//! - Minimizer-based indexing
//! - Sparse and skew index for O(1) lookup

use crate::{
    kmer::{Kmer, KmerBits},
    minimizer::MinimizerInfo,
    minimizers_control_map::MinimizersControlMap,
    sparse_and_skew_index::SparseAndSkewIndex,
    spectrum_preserving_string_set::SpectrumPreservingStringSet,
    constants::{ceil_log2, INVALID_UINT64, MIN_L},
};
use value_traits::slices::SliceByValue;
use tracing::info;

/// The main dictionary structure
///
/// Note: Serialization support is limited due to underlying types.
/// Use custom save/load methods instead.
pub struct Dictionary {
    /// Encoded string storage
    spss: SpectrumPreservingStringSet,

    /// Maps minimizers to control information
    control_map: MinimizersControlMap,

    /// Sparse and skew index for lookups
    index: SparseAndSkewIndex,

    /// K-mer size
    k: usize,

    /// Minimizer size
    m: usize,

    /// Build seed: the minimizer hasher is derived from it, so it is part of
    /// the index identity (stored in the serialized header).
    seed: u64,

    /// Cached hasher for minimizer extraction (avoids re-creating RandomState per query)
    hasher: crate::hasher::DeterministicHasher,
}

/// Successful k-mer lookup with cached locate info to avoid
/// a redundant second Elias-Fano successor query.
#[derive(Clone, Copy)]
struct LocatedHit {
    kmer_offset: u64,
    orientation: i64,
    string_id: u64,
    string_begin: u64,
    string_end: u64,
}

impl LocatedHit {
    #[inline]
    fn into_lookup_result(self, k: usize) -> crate::streaming_query::LookupResult {
        let kmer_id_in_string = self.kmer_offset - self.string_begin;
        let kmer_id = self.kmer_offset - self.string_id * (k as u64 - 1);
        crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset: self.kmer_offset,
            kmer_orientation: self.orientation,
            string_id: self.string_id,
            string_begin: self.string_begin,
            string_end: self.string_end,
            minimizer_found: true,
        }
    }
}

impl Dictionary {
    /// Create a new dictionary from components.
    ///
    /// `seed` must be the seed the index was *built* with: queries hash
    /// minimizers with a hasher derived from it, so a mismatch silently
    /// returns zero hits.
    pub fn new(
        spss: SpectrumPreservingStringSet,
        control_map: MinimizersControlMap,
        index: SparseAndSkewIndex,
        k: usize,
        m: usize,
        seed: u64,
    ) -> Self {
        Self {
            spss,
            control_map,
            index,
            k,
            m,
            seed,
            hasher: crate::hasher::DeterministicHasher::new(seed),
        }
    }

    /// Look up a k-mer's position in the dictionary.
    ///
    /// The index is canonical: a k-mer and its reverse complement share a
    /// bucket, and one probe answers for both strands.
    #[inline]
    pub fn lookup<const K: usize>(&self, kmer: &Kmer<K>) -> u64
    where
        Kmer<K>: KmerBits,
    {
        self.query(kmer).kmer_offset
    }

    /// Alias of [`Self::lookup`]: in the unified canonical index both strands
    /// are equivalent, so there is no strand-specific lookup to skip.
    #[inline]
    pub fn lookup_forward<const K: usize>(&self, kmer: &Kmer<K>) -> u64
    where
        Kmer<K>: KmerBits,
    {
        self.lookup(kmer)
    }

    /// Alias of [`Self::lookup`]; the flag is ignored (the index is canonical,
    /// matching C++ v6 where the modality no longer exists).
    #[inline]
    pub fn lookup_checked<const K: usize>(&self, kmer: &Kmer<K>, _check_reverse_complement: bool) -> u64
    where
        Kmer<K>: KmerBits,
    {
        self.lookup(kmer)
    }

    /// Query a k-mer and return a full [`LookupResult`](crate::streaming_query::LookupResult),
    /// including the orientation (+1 forward / -1 backward) of the match.
    #[inline]
    pub fn query<const K: usize>(&self, kmer: &Kmer<K>) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let kmer_rc = kmer.reverse_complement();
        let mini = self.extract_minimizer(kmer, &kmer_rc);
        self.lookup_with_minimizer(kmer, &kmer_rc, mini, None)
    }

    /// Alias of [`Self::query`]: both strands are equivalent in the canonical
    /// index (the previous strand-specific behavior existed only for the
    /// removed non-canonical modality).
    #[inline]
    pub fn query_forward<const K: usize>(&self, kmer: &Kmer<K>) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        self.query(kmer)
    }

    /// Alias of [`Self::query`]; the flag is ignored (the index is canonical).
    #[inline]
    pub fn query_checked<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        _check_reverse_complement: bool,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        self.query(kmer)
    }

    /// Query a k-mer given as a DNA string.
    pub fn query_from_str<const K: usize>(&self, kmer_str: &str) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        assert_eq!(kmer_str.len(), self.k, "k-mer string length must equal k={}", self.k);
        match Kmer::<K>::from_str(kmer_str) {
            Ok(kmer) => self.query(&kmer),
            Err(_) => crate::streaming_query::LookupResult::not_found(),
        }
    }

    /// Debug: check if control map MPHF is present
    pub fn debug_has_control_map_mphf(&self) -> bool {
        self.control_map.mphf_ref().is_some()
    }

    /// Debug: check if a minimizer exists in the control map
    pub fn debug_control_map_lookup(&self, minimizer: u64) -> bool {
        self.control_map.lookup(minimizer).is_some()
    }

    /// Debug: extract minimizer info for a k-mer
    pub fn debug_extract_minimizer<const K: usize>(&self, kmer: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        self.extract_minimizer(kmer, &kmer.reverse_complement())
    }

    /// Debug: get bucket and locate info for a k-mer
    pub fn debug_bucket_info<const K: usize>(&self, kmer: &Kmer<K>) -> Option<(u64, u64, u64)>
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer, &kmer.reverse_complement());
        let bucket_id = self.control_map.lookup(minimizer_info.value)?;
        if bucket_id >= self.index.num_buckets() {
            return None;
        }
        let (begin, end) = self.index.locate_bucket(bucket_id);
        Some((minimizer_info.value, bucket_id as u64, (end - begin) as u64))
    }

    /// Look up a k-mer and return position + orientation
    /// (+1 forward / -1 backward; `(INVALID_UINT64, 1)` when absent).
    #[inline]
    pub fn lookup_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let res = self.query(kmer);
        (res.kmer_offset, res.kmer_orientation as i8)
    }

    /// Alias of [`Self::lookup_with_orientation`]: both strands are equivalent
    /// in the unified canonical index.
    #[inline]
    pub fn lookup_forward_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        self.lookup_with_orientation(kmer)
    }

    // -----------------------------------------------------------------------
    // The one lookup path (unified canonical scheme)
    // -----------------------------------------------------------------------

    /// True if the m-mer stored at `minimizer_pos` is the minimizer or its
    /// reverse complement. The minimizer's value is the *canonical* m-mer at
    /// the anchored locus, so the text holds either orientation of it.
    #[inline]
    fn mmer_matches(&self, minimizer_pos: u64, value: u64) -> bool {
        let read = self.spss.decode_mmer_at(minimizer_pos as usize, self.m);
        read == value || crate::minimizer::reverse_complement_mmer(read, self.m) == value
    }

    /// Core lookup with a pre-computed minimizer.
    ///
    /// When `cache` is provided, bucket bounds are memoized across k-mers
    /// sharing a minimizer, and for non-heavy buckets of at most
    /// [`BucketCache::MAX_CACHED_POSITIONS`] entries the decoded locate set is
    /// stored in it (at no extra cost — it is decoded here anyway), letting a
    /// streaming query verify follow-up k-mers with the same minimizer via
    /// [`Self::lookup_at_positions`] with no MPHF/EF/offset work.
    ///
    /// `minimizer_found` semantics: `false` means the minimizer is provably
    /// absent from the index (any other k-mer carrying it is absent too);
    /// on the heavy/skew path the probed position can belong to a different
    /// minimizer without implying absence, so `minimizer_found` stays `true`
    /// there even on mismatch (port of the C++ HEAVYLOAD comment).
    pub(crate) fn lookup_with_minimizer<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        mini: MinimizerInfo,
        mut cache: Option<&mut BucketCache>,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let bounds = match cache.as_deref_mut() {
            Some(c) => {
                c.size = 0;
                self.resolve_bucket(mini.value, c)
            }
            None => self
                .control_map
                .lookup(mini.value)
                .filter(|&id| id < self.index.num_buckets())
                .map(|id| self.index.locate_bucket(id)),
        };
        let Some((begin, end)) = bounds else {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        };
        let n = end - begin;
        if n == 0 {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        }

        let log2_size = ceil_log2(n as u64);
        if log2_size > MIN_L {
            // Heavy bucket: the skew index keys on the canonical k-mer.
            let key = std::cmp::min(kmer.bits(), kmer_rc.bits());
            let within_pos = self.index.skew_index.lookup(&key, log2_size);
            if within_pos == INVALID_UINT64 || within_pos as usize >= n {
                return crate::streaming_query::LookupResult::not_found();
            }
            let minimizer_pos = self.index.offsets.index_value(begin + within_pos as usize) as u64;
            if !self.mmer_matches(minimizer_pos, mini.value) {
                // A foreign k-mer probes an arbitrary in-bucket position, so a
                // mismatch here says nothing about the minimizer's presence.
                return crate::streaming_query::LookupResult::not_found();
            }
            match self.lookup_at_position(kmer, kmer_rc, minimizer_pos, mini.pos_in_kmer) {
                Some(hit) => hit.into_lookup_result(self.k),
                None => crate::streaming_query::LookupResult::not_found(),
            }
        } else {
            // Singleton or light bucket: decode the locate set, then verify.
            let mut buf = [0u64; 1 << MIN_L];
            for (i, slot) in buf[..n].iter_mut().enumerate() {
                *slot = self.index.offsets.index_value(begin + i) as u64;
            }
            // The bucket's positions all anchor the same minimizer, so one
            // presence probe decides for the whole set: a mismatch means an
            // out-of-set minimizer collided into this bucket through the MPHF.
            if !self.mmer_matches(buf[0], mini.value) {
                let mut res = crate::streaming_query::LookupResult::not_found();
                res.minimizer_found = false;
                return res;
            }
            if let Some(c) = cache {
                if n <= BucketCache::MAX_CACHED_POSITIONS {
                    c.positions[..n].copy_from_slice(&buf[..n]);
                    c.size = n;
                }
            }
            self.lookup_at_positions(&buf[..n], kmer, kmer_rc, mini)
        }
    }

    /// The verification half of [`Self::lookup_with_minimizer`]: try the
    /// candidate loci of each position of the minimizer's locate set. The
    /// positions must be the decoded locate set of `mini.value`, whose
    /// presence at `positions[0]` has already been checked.
    pub(crate) fn lookup_at_positions<const K: usize>(
        &self,
        positions: &[u64],
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        mini: MinimizerInfo,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        for &minimizer_pos in positions {
            if let Some(hit) = self.lookup_at_position(kmer, kmer_rc, minimizer_pos, mini.pos_in_kmer) {
                return hit.into_lookup_result(self.k);
            }
        }
        crate::streaming_query::LookupResult::not_found()
    }

    /// The minimizer is anchored at locus `pos_in_kmer` of the k-mer and, by
    /// mirror-equivariance, at locus `k - m - pos_in_kmer` of its reverse
    /// complement. So a minimizer occurrence at offset j anchors a k-mer
    /// starting either at `j - pos_in_kmer` or at `j - (k - m - pos_in_kmer)`:
    /// always exactly two candidates, with no case analysis.
    #[inline]
    fn lookup_at_position<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        pos_in_kmer: usize,
    ) -> Option<LocatedHit>
    where
        Kmer<K>: KmerBits,
    {
        if let Some(hit) = self.try_lookup_at(kmer, kmer_rc, minimizer_pos, pos_in_kmer as u64) {
            return Some(hit);
        }
        self.try_lookup_at(kmer, kmer_rc, minimizer_pos, (K - self.m - pos_in_kmer) as u64)
    }

    /// Verify one candidate k-mer start against the SPSS text.
    #[inline]
    fn try_lookup_at<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        pos_in_kmer: u64,
    ) -> Option<LocatedHit>
    where
        Kmer<K>: KmerBits,
    {
        let kmer_pos = minimizer_pos.checked_sub(pos_in_kmer)?;

        let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);

        let orientation: i64 = if stored_kmer.bits() == query_kmer.bits() {
            1
        } else if stored_kmer.bits() == kmer_rc.bits() {
            -1
        } else {
            return None;
        };

        let (string_id, string_begin, string_end) = self.spss.locate_with_end(kmer_pos)?;
        if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
            Some(LocatedHit { kmer_offset: kmer_pos, orientation, string_id, string_begin, string_end })
        } else {
            None
        }
    }

    /// Check if a k-mer exists in the dictionary
    pub fn access<const K: usize>(&self, kmer: &Kmer<K>) -> bool
    where
        Kmer<K>: KmerBits,
    {
        self.lookup(kmer) != INVALID_UINT64
    }

    /// Create a streaming query engine for this dictionary
    pub fn create_streaming_query<const K: usize>(&self) -> crate::streaming_query::StreamingQueryEngine<'_, K>
    where
        Kmer<K>: KmerBits,
    {
        crate::streaming_query::StreamingQueryEngine::new(self)
    }

    /// Get the k-mer size
    pub fn k(&self) -> usize {
        self.k
    }

    /// Get the minimizer size
    pub fn m(&self) -> usize {
        self.m
    }

    /// Whether the index is canonical. Always `true` since the unified
    /// v6-style scheme: a k-mer and its reverse complement share a bucket and
    /// every lookup reports orientation. Kept so downstream code compiles.
    #[deprecated(since = "0.7.0", note = "the index is always canonical; this returns true")]
    pub fn canonical(&self) -> bool {
        true
    }

    /// The build seed the index was constructed (and is queried) with.
    pub fn seed(&self) -> u64 {
        self.seed
    }

    /// Get a reference to the underlying SPSS
    pub fn spss(&self) -> &SpectrumPreservingStringSet {
        &self.spss
    }

    /// Get a reference to the control map
    pub fn control_map_ref(&self) -> &MinimizersControlMap {
        &self.control_map
    }

    /// Get a reference to the sparse and skew index
    pub fn index_ref(&self) -> &SparseAndSkewIndex {
        &self.index
    }

    /// Debug: get reference to SPSS for testing
    #[cfg(test)]
    pub fn debug_spss(&self) -> &SpectrumPreservingStringSet {
        &self.spss
    }

    /// Get the number of strings in the SPSS
    pub fn num_strings(&self) -> u64 {
        self.spss.num_strings()
    }

    /// Get the length of a specific string in bases
    pub fn string_length(&self, string_id: u64) -> usize {
        self.spss.string_length(string_id)
    }

    /// Locate which string contains a given absolute position.
    #[inline]
    pub fn locate_string(&self, absolute_pos: u64) -> Option<(u64, u64)> {
        self.spss.locate(absolute_pos)
    }

    /// Access a k-mer at a given position within a string
    pub fn access_kmer<const K: usize>(&self, string_id: u64, pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        self.spss.decode_kmer::<K>(string_id, pos)
    }

    /// Decode the k-mer at an absolute base position in the SPSS.
    #[inline]
    pub fn kmer_at_pos<const K: usize>(&self, absolute_base_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        self.spss.decode_kmer_at(absolute_base_pos)
    }

    /// Get the number of unique minimizers
    pub fn num_minimizers(&self) -> u64 {
        self.control_map.num_minimizers()
    }

    /// Get total memory usage in bits
    pub fn num_bits(&self) -> u64 {
        self.spss.num_bits()
            + self.control_map.num_bits()
            + self.index.num_bits()
    }

    /// Print a detailed space breakdown of the index
    pub fn print_space_breakdown(&self) {
        let num_kmers = self.spss.total_bases().saturating_sub(
            (self.k as u64 - 1) * self.spss.num_strings()
        ) as f64;

        let strings_bytes = self.spss.strings_bytes() as f64;
        let offsets_bytes = self.spss.offsets_bytes() as f64;
        let ef_bytes = self.index.ef_bytes() as f64;
        let index_offsets_bytes = self.index.offsets_bytes() as f64;
        let skew_bytes = self.index.skew_index_bytes() as f64;
        let mphf_bytes = self.control_map.mphf_serialized_bytes() as f64;
        let skew_mphf_bytes = self.index.skew_mphf_bytes() as f64;

        let total = strings_bytes + offsets_bytes + ef_bytes
            + index_offsets_bytes + skew_bytes
            + mphf_bytes + skew_mphf_bytes;

        let perc = |x: f64| -> f64 { x * 100.0 / total };

        info!("total index size: {} [B] -- {:.5} [MB] ({:.5} [bits/kmer])",
            total as u64, total / 1_000_000.0, total * 8.0 / num_kmers);
        info!("SPACE BREAKDOWN:");
        info!("  mphf: {:.6} [bits/kmer] ({:.5} [bits/key]) -- {:.4}%",
            mphf_bytes * 8.0 / num_kmers,
            mphf_bytes * 8.0 / self.num_minimizers() as f64,
            perc(mphf_bytes));
        info!("  strings_offsets: {:.6} [bits/kmer] -- {:.5}%",
            offsets_bytes * 8.0 / num_kmers, perc(offsets_bytes));
        info!("  num_super_kmers_before_bucket (EF): {:.5} [bits/kmer] -- {:.4}%",
            ef_bytes * 8.0 / num_kmers, perc(ef_bytes));
        info!("  offsets: {:.6} [bits/kmer] -- {:.5}%",
            index_offsets_bytes * 8.0 / num_kmers, perc(index_offsets_bytes));
        info!("  strings: {:.5} [bits/kmer] -- {:.4}%",
            strings_bytes * 8.0 / num_kmers, perc(strings_bytes));
        info!("  skew_index: {:.6} [bits/kmer] -- {:.5}%",
            skew_bytes * 8.0 / num_kmers, perc(skew_bytes));
        info!("  skew_mphfs: {:.6} [bits/kmer] -- {:.5}%",
            skew_mphf_bytes * 8.0 / num_kmers, perc(skew_mphf_bytes));
        info!("  --------------");
        info!("  total: {:.5} [bits/kmer]", total * 8.0 / num_kmers);
    }

    /// Serialize the dictionary to files
    pub fn save<P: AsRef<std::path::Path>>(&self, path: P) -> crate::serialization::SerializationResult<()> {
        use crate::serialization::*;
        use std::io::BufWriter;

        let base_path = path.as_ref();

        // Create main index file
        let index_path = index_file_path(base_path);
        let index_file = std::fs::File::create(&index_path)?;
        let mut index_writer = BufWriter::new(index_file);

        // Write header
        let header = DictionarySerializationHeader::new(
            self.k,
            self.m,
            self.seed,
            (self.index.skew_index.num_partitions() + 1) as u32,
        );
        header.write(&mut index_writer)?;

        // Write SPSS
        self.spss.serialize_to(&mut index_writer)?;

        // Write control map without MPHF
        self.control_map.serialize_without_mphf(&mut index_writer)?;

        // Write index metadata excluding MPHF
        self.index.serialize_without_mphf(&mut index_writer)?;

        // Create MPHF container file
        let mphf_path = mphf_container_path(base_path);
        let mphf_file = std::fs::File::create(&mphf_path)?;
        let mut mphf_writer = BufWriter::new(mphf_file);

        // Write MPHF container (control_map first, then skew_index partitions)
        let mut mphfs: Vec<Option<&crate::partitioned_mphf::PartitionedMphf>> = Vec::with_capacity(
            self.index.skew_index.num_partitions() + 1,
        );
        mphfs.push(self.control_map.mphf_ref());
        mphfs.extend(self.index.skew_index.mphfs_ref().iter().map(|o| o.as_ref()));
        write_mphf_container(&mut mphf_writer, &mphfs)?;

        Ok(())
    }

    /// Deserialize a dictionary from files
    pub fn load<P: AsRef<std::path::Path>>(path: P) -> crate::serialization::SerializationResult<Self> {
        use crate::serialization::*;
        use std::io::BufReader;

        let base_path = path.as_ref();

        // Load main index file
        let index_path = index_file_path(base_path);
        let index_file = std::fs::File::open(&index_path)?;
        let mut index_reader = BufReader::new(index_file);

        // Read header
        let header = DictionarySerializationHeader::read(&mut index_reader)?;

        // Read SPSS
        let spss = SpectrumPreservingStringSet::deserialize_from(&mut index_reader)?;

        // Read control map without MPHF
        let mut control_map = MinimizersControlMap::deserialize_without_mphf(&mut index_reader)?;

        // Read index metadata without MPHF
        let mut index = SparseAndSkewIndex::deserialize_without_mphf(&mut index_reader)?;

        // Load MPHF container
        let mphf_path = mphf_container_path(base_path);
        let mphf_file = std::fs::File::open(&mphf_path)?;
        let mut mphf_reader = BufReader::new(mphf_file);
        let mut mphfs = read_mphf_container(&mut mphf_reader)?;

        // First MPHF is the control_map's; remaining are skew_index partitions
        let control_mphf = if !mphfs.is_empty() { mphfs.remove(0) } else { None };
        control_map.set_mphf(control_mphf);
        index.skew_index.set_mphfs(mphfs);

        Ok(Dictionary {
            spss,
            control_map,
            index,
            k: header.k,
            m: header.m,
            seed: header.seed,
            hasher: crate::hasher::DeterministicHasher::new(header.seed),
        })
    }

    /// Extract the minimizer from a k-mer using the cached hasher. Delegates
    /// to the shared [`crate::minimizer::compute_minimizer`] so builder, point
    /// lookup, and streaming all compute the identical scheme.
    #[inline]
    pub(crate) fn extract_minimizer<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
    ) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        crate::minimizer::compute_minimizer(kmer, kmer_rc, self.m, &self.hasher)
    }

    // --- Streaming query helpers ---

    /// Resolve a minimizer to its bucket bounds, reusing `cache` when the
    /// minimizer is unchanged from the previous k-mer.
    ///
    /// `None` means the minimizer is not in the index at all (the caller must
    /// report `minimizer_found = false`); `Some((begin, end))` are the bounds.
    #[inline(always)]
    fn resolve_bucket(&self, mini_value: u64, cache: &mut BucketCache) -> Option<(usize, usize)> {
        if cache.valid && cache.mini_value == mini_value {
            return cache.bounds;
        }
        let bounds = self
            .control_map
            .lookup(mini_value)
            .filter(|&id| id < self.index.num_buckets())
            .map(|id| self.index.locate_bucket(id));
        cache.valid = true;
        cache.mini_value = mini_value;
        cache.bounds = bounds;
        bounds
    }


}

/// Memoized bucket-bounds resolution for a streaming query.
///
/// Locating a bucket costs an MPHF evaluation (`PartitionedMphf::get`) plus two
/// Elias-Fano accesses (`DArray::select`), and depends *only* on the minimizer
/// value. Consecutive k-mers share a minimizer for runs averaging `(k-m+2)/2`,
/// so most seeds re-derive bounds they just computed: measured at 71.3% of
/// 46.3M seeds on gencode v49 (k=31, m=19), all of them already downstream of
/// the unitig-extension path, which handles a disjoint set of k-mers.
///
/// The mapping minimizer -> bounds is a pure function of an immutable
/// dictionary, so an entry never goes stale and needs no invalidation, not even
/// across reads. Only the *bounds* are cached; the per-k-mer `MinimizerInfo`
/// (which carries the minimizer's position within the k-mer) still flows into
/// the bucket scan untouched, so results are unchanged.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct BucketCache {
    /// Whether `mini_value`/`bounds` hold anything. Needed because any `u64`,
    /// including `u64::MAX`, is a legal minimizer value.
    valid: bool,
    mini_value: u64,
    /// `None` records a minimizer that is absent from the index -- worth
    /// caching too, since that is the other repeated outcome.
    bounds: Option<(usize, usize)>,
    /// Number of decoded positions held in `positions` (0 = nothing cached).
    /// Filled by `lookup_with_minimizer` for non-heavy buckets of at most
    /// `MAX_CACHED_POSITIONS` entries, at no extra decode cost; consumed by
    /// the streaming same-minimizer memo via `lookup_at_positions`.
    pub(crate) size: usize,
    /// The decoded locate set of `mini_value` (valid for `..size`).
    pub(crate) positions: [u64; BucketCache::MAX_CACHED_POSITIONS],
}

impl BucketCache {
    /// Largest bucket whose decoded locate set is cached. Heavy buckets are
    /// never cached; larger light buckets are too rare to be worth more space
    /// (mirrors C++ `bucket_cache::max_cached_size`).
    pub(crate) const MAX_CACHED_POSITIONS: usize = 8;

    /// Invalidate everything (used by `StreamingQuery::reset`).
    pub(crate) fn reset(&mut self) {
        self.valid = false;
        self.size = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::builder::{BuildConfiguration, DictionaryBuilder};

    #[test]
    fn test_dictionary_creation() {
        let spss = SpectrumPreservingStringSet::new(31, 13);
        let control_map = MinimizersControlMap::from_parts(
            Vec::new(),
            Vec::new(),
            0,
        );
        let index = SparseAndSkewIndex::new();

        let dict = Dictionary::new(spss, control_map, index, 31, 13, 1);

        assert_eq!(dict.k(), 31);
        assert_eq!(dict.m(), 13);
        assert_eq!(dict.seed(), 1);
    }

    #[test]
    fn test_dictionary_build_and_lookup() {
        let mut config = BuildConfiguration::new(31, 21).unwrap();
        config.verbose = true;
        let builder = DictionaryBuilder::new(config).unwrap();

        let sequences = vec![
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_string(),
        ];

        eprintln!("\n=== Building dictionary ===");
        let dict = builder.build_from_sequences(sequences.clone()).unwrap();

        assert_eq!(dict.k(), 31);
        assert_eq!(dict.m(), 21);

        eprintln!("\n=== Dictionary info ===");
        eprintln!("Num minimizers: {}", dict.num_minimizers());
        eprintln!("Num buckets: {}", dict.index.num_buckets());
        eprintln!("Num offsets: {}", dict.index.num_offsets());
        eprintln!("SPSS num strings: {}", dict.spss.num_strings());
        eprintln!("SPSS total bases: {}", dict.spss.total_bases());

        // Create a k-mer from the first sequence
        let test_kmer_str = &sequences[0][0..31];
        eprintln!("\n=== Testing lookup for k-mer: {} ===", test_kmer_str);
        let kmer = crate::kmer::Kmer::<31>::from_str(test_kmer_str).unwrap();

        // Extract minimizer
        let mut mini_iter = crate::minimizer::MinimizerIterator::with_seed(31, 21, 1);
        let mini_info = mini_iter.next(&kmer, &kmer.reverse_complement());
        eprintln!("K-mer minimizer: value={}, pos_in_kmer={}", mini_info.value, mini_info.pos_in_kmer);

        // Look up in control map
        if let Some(bucket_id) = dict.control_map.lookup(mini_info.value) {
            eprintln!("Control map lookup: bucket_id={}", bucket_id);

            if bucket_id < dict.index.num_buckets() {
                let (begin, end) = dict.index.locate_bucket(bucket_id);
                eprintln!("Bucket range: [{}, {}), size={}", begin, end, end - begin);
            }
        } else {
            eprintln!("Minimizer NOT found in control map!");
        }

        // Lookup the k-mer
        let result = dict.lookup(&kmer);

        eprintln!("\nLookup result: {} (INVALID={})", result, crate::constants::INVALID_UINT64);
    }

    #[test]
    fn test_unified_index_answers_both_strands() {
        // The unified scheme: one probe answers both strands, with orientation.
        let config = BuildConfiguration::new(7, 5).unwrap();
        let builder = DictionaryBuilder::new(config).unwrap();
        let dict = builder
            .build_from_sequences(vec!["ACGTTGC".to_string()])
            .unwrap();

        let fwd = Kmer::<7>::from_str("ACGTTGC").unwrap();
        let rc = fwd.reverse_complement();
        assert_ne!(fwd.bits(), rc.bits(), "k-mer must not be its own RC");

        // Both orientations found through every query entry point.
        assert!(dict.query::<7>(&fwd).is_found());
        assert!(dict.query::<7>(&rc).is_found());
        assert!(dict.query_forward::<7>(&rc).is_found());
        assert!(dict.query_checked::<7>(&rc, false).is_found());

        // Orientation is reported: forward => +1, reverse => -1, and the
        // matched text position is the same for both strands.
        assert_eq!(dict.lookup_with_orientation::<7>(&fwd).1, 1);
        assert_eq!(dict.lookup_with_orientation::<7>(&rc).1, -1);
        assert_eq!(
            dict.query::<7>(&fwd).kmer_offset,
            dict.query::<7>(&rc).kmer_offset
        );
        assert_ne!(dict.lookup_forward::<7>(&rc), INVALID_UINT64);

        // An absent k-mer (and its RC) is not found, and the minimizer
        // presence flag is meaningful.
        let absent = Kmer::<7>::from_str("AAAAAAA").unwrap();
        assert!(!dict.query::<7>(&absent).is_found());
        assert!(!dict.query::<7>(&absent.reverse_complement()).is_found());
    }
}
