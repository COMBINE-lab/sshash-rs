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

    /// Whether to use canonical mode (k-mer or RC, whichever is smaller)
    canonical: bool,

    /// Cached hasher for minimizer extraction (avoids re-creating RandomState per query)
    hasher: crate::hasher::DeterministicHasher,
}

impl Dictionary {
    /// Create a new dictionary from components
    pub fn new(
        spss: SpectrumPreservingStringSet,
        control_map: MinimizersControlMap,
        index: SparseAndSkewIndex,
        k: usize,
        m: usize,
        canonical: bool,
    ) -> Self {
        Self {
            spss,
            control_map,
            index,
            k,
            m,
            canonical,
            hasher: crate::hasher::DeterministicHasher::new(1),
        }
    }

    /// Look up a k-mer's position in the dictionary.
    #[inline]
    pub fn lookup<const K: usize>(&self, kmer: &Kmer<K>) -> u64
    where
        Kmer<K>: KmerBits,
    {
        let (pos, _) = self.lookup_with_orientation(kmer);
        pos
    }

    /// Query a k-mer and return a full [`LookupResult`](crate::streaming_query::LookupResult).
    #[inline]
    pub fn query<const K: usize>(&self, kmer: &Kmer<K>) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        if self.canonical {
            return self.query_canonical(kmer);
        }

        // Non-canonical mode: try forward k-mer first, then RC
        if let Some(res) = self.query_regular(kmer, 1) {
            return res;
        }
        let kmer_rc = kmer.reverse_complement();
        if let Some(res) = self.query_regular(&kmer_rc, -1) {
            return res;
        }
        crate::streaming_query::LookupResult::not_found()
    }

    /// Internal: regular (non-canonical) lookup returning full LookupResult.
    #[inline]
    fn query_regular<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        orientation: i8,
    ) -> Option<crate::streaming_query::LookupResult>
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let bucket_id = self.control_map.lookup(minimizer_info.value)?;
        if bucket_id >= self.index.num_buckets() {
            return None;
        }

        let kmer_pos = self.lookup_in_bucket::<K>(kmer, minimizer_info, bucket_id);

        if kmer_pos == INVALID_UINT64 {
            return None;
        }

        let (string_id, string_begin, string_end) = self.spss.locate_with_end(kmer_pos)?;
        let kmer_id_in_string = kmer_pos - string_begin;
        let kmer_id = kmer_pos - string_id * (self.k as u64 - 1);

        Some(crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset: kmer_pos,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        })
    }

    /// Internal: canonical lookup returning full LookupResult.
    #[inline]
    fn query_canonical<const K: usize>(
        &self,
        kmer: &Kmer<K>,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let (kmer_offset, orientation) = self.lookup_canonical_with_orientation(kmer);
        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        let (string_id, string_begin, string_end) = match self.spss.locate_with_end(kmer_offset) {
            Some(triple) => triple,
            None => return crate::streaming_query::LookupResult::not_found(),
        };
        if kmer_offset > string_end - self.k as u64 {
            return crate::streaming_query::LookupResult::not_found();
        }
        let kmer_id_in_string = kmer_offset - string_begin;
        let kmer_id = kmer_offset - string_id * (self.k as u64 - 1);

        crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        }
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
        self.extract_minimizer(kmer)
    }

    /// Debug: get bucket and locate info for a k-mer
    pub fn debug_bucket_info<const K: usize>(&self, kmer: &Kmer<K>) -> Option<(u64, u64, u64)>
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let bucket_id = self.control_map.lookup(minimizer_info.value)?;
        if bucket_id >= self.index.num_buckets() {
            return None;
        }
        let (begin, end) = self.index.locate_bucket(bucket_id);
        Some((minimizer_info.value, bucket_id as u64, (end - begin) as u64))
    }

    /// Look up a k-mer and return position + orientation
    #[inline]
    pub fn lookup_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        if self.canonical {
            self.lookup_canonical_with_orientation(kmer)
        } else {
            let (pos, ori) = self.lookup_regular_with_orientation(kmer);
            if pos != INVALID_UINT64 {
                return (pos, ori);
            }
            let kmer_rc = kmer.reverse_complement();
            let (pos_rc, _) = self.lookup_regular_with_orientation(&kmer_rc);
            if pos_rc != INVALID_UINT64 {
                return (pos_rc, -1);
            }
            (INVALID_UINT64, 1)
        }
    }

    #[inline]
    fn lookup_regular_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let minimizer_info = self.extract_minimizer(kmer);
        let bucket_id = match self.control_map.lookup(minimizer_info.value) {
            Some(id) => id,
            None => return (INVALID_UINT64, 1),
        };
        if bucket_id >= self.index.num_buckets() {
            return (INVALID_UINT64, 1);
        }

        let pos = self.lookup_in_bucket::<K>(kmer, minimizer_info, bucket_id);
        (pos, 1)
    }

    #[inline]
    fn lookup_canonical_with_orientation<const K: usize>(&self, kmer: &Kmer<K>) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let kmer_rc = kmer.reverse_complement();
        let mini_fwd = self.extract_minimizer(kmer);
        let mini_rc = self.extract_minimizer(&kmer_rc);

        if mini_fwd.value < mini_rc.value {
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_fwd)
        } else if mini_rc.value < mini_fwd.value {
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_rc)
        } else {
            let res = self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_fwd);
            if res.0 != INVALID_UINT64 {
                return res;
            }
            self.lookup_canonical_with_minimizer::<K>(kmer, &kmer_rc, mini_rc)
        }
    }

    #[inline]
    fn lookup_canonical_with_minimizer<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_info: MinimizerInfo,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(minimizer_info.value) {
            Some(id) => id,
            None => return (INVALID_UINT64, 1),
        };
        if bucket_id >= self.index.num_buckets() {
            return (INVALID_UINT64, 1);
        }

        self.lookup_in_bucket_canonical::<K>(kmer, kmer_rc, minimizer_info, bucket_id)
    }

    // -----------------------------------------------------------------------
    // Core bucket lookup helpers (replaces old 3-way LSB dispatch)
    // -----------------------------------------------------------------------

    /// Look up a k-mer in a bucket using the locate_bucket + size-based dispatch.
    ///
    /// For small buckets (ceil_log2(size) <= MIN_L): linear scan of offsets.
    /// For large buckets: skew index MPHF lookup.
    #[inline]
    fn lookup_in_bucket<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        bucket_id: usize,
    ) -> u64
    where
        Kmer<K>: KmerBits,
    {
        let (begin, end) = self.index.locate_bucket(bucket_id);
        let n = end - begin;
        if n == 0 {
            return INVALID_UINT64;
        }

        let log2_size = ceil_log2(n as u64);
        if log2_size > MIN_L {
            // Heavy bucket: use skew index
            let within_pos = self.index.skew_index.lookup(&kmer.bits(), log2_size);
            if within_pos == INVALID_UINT64 || within_pos as usize >= n {
                return INVALID_UINT64;
            }
            let minimizer_pos = self.index.offsets.index_value(begin + within_pos as usize) as u64;
            self.lookup_from_minimizer_pos::<K>(kmer, minimizer_pos, minimizer_info)
        } else {
            // Singleton or light bucket: linear scan
            self.lookup_bucket_linear::<K>(kmer, minimizer_info, begin, end)
        }
    }

    /// Canonical lookup in a bucket using locate_bucket + size-based dispatch.
    #[inline]
    fn lookup_in_bucket_canonical<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        bucket_id: usize,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        let (begin, end) = self.index.locate_bucket(bucket_id);
        let n = end - begin;
        if n == 0 {
            return (INVALID_UINT64, 1);
        }

        let log2_size = ceil_log2(n as u64);
        if log2_size > MIN_L {
            // Heavy bucket: use skew index with canonical k-mer
            let kmer_canon_value = std::cmp::min(kmer.bits(), kmer_rc.bits());
            let within_pos = self.index.skew_index.lookup(&kmer_canon_value, log2_size);
            if within_pos == INVALID_UINT64 || within_pos as usize >= n {
                return (INVALID_UINT64, 1);
            }
            let minimizer_pos = self.index.offsets.index_value(begin + within_pos as usize) as u64;
            self.lookup_from_minimizer_pos_canonical::<K>(kmer, kmer_rc, minimizer_pos, minimizer_info)
        } else {
            // Singleton or light bucket: linear scan
            self.lookup_bucket_linear_canonical::<K>(kmer, kmer_rc, minimizer_info, begin, end)
        }
    }

    /// Linear scan through offsets[begin..end] for regular lookup.
    #[inline]
    fn lookup_bucket_linear<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        begin: usize,
        end: usize,
    ) -> u64
    where
        Kmer<K>: KmerBits,
    {
        for i in begin..end {
            let minimizer_pos = self.index.offsets.index_value(i) as u64;
            let pos = self.lookup_from_minimizer_pos::<K>(query_kmer, minimizer_pos, minimizer_info);
            if pos != INVALID_UINT64 {
                return pos;
            }
        }
        INVALID_UINT64
    }

    /// Linear scan for canonical lookup (tries both fwd and RC pos_in_kmer).
    #[inline]
    fn lookup_bucket_linear_canonical<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_info: MinimizerInfo,
        begin: usize,
        end: usize,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        for i in begin..end {
            let minimizer_pos = self.index.offsets.index_value(i) as u64;
            let (pos, ori) = self.lookup_from_minimizer_pos_canonical::<K>(
                query_kmer, kmer_rc, minimizer_pos, minimizer_info,
            );
            if pos != INVALID_UINT64 {
                return (pos, ori);
            }
        }
        (INVALID_UINT64, 1)
    }

    // -----------------------------------------------------------------------
    // Low-level position verification helpers (unchanged)
    // -----------------------------------------------------------------------

    #[inline]
    fn lookup_from_minimizer_pos<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        minimizer_pos: u64,
        minimizer_info: MinimizerInfo,
    ) -> u64
    where
        Kmer<K>: KmerBits,
    {
        let kmer_pos = match minimizer_pos.checked_sub(minimizer_info.pos_in_kmer as u64) {
            Some(pos) => pos,
            None => return INVALID_UINT64,
        };

        let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);

        if stored_kmer.bits() == query_kmer.bits() {
            if let Some((_sid, string_begin, string_end)) = self.spss.locate_with_end(kmer_pos) {
                if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
                    return kmer_pos;
                }
            }
        }

        INVALID_UINT64
    }

    /// Canonical lookup from minimizer position
    #[inline]
    fn lookup_from_minimizer_pos_canonical<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        minimizer_info: MinimizerInfo,
    ) -> (u64, i8)
    where
        Kmer<K>: KmerBits,
    {
        // Try forward position first
        let pos_in_kmer_fwd = minimizer_info.pos_in_kmer as u64;
        if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
            query_kmer, kmer_rc, minimizer_pos, pos_in_kmer_fwd
        ) {
            return (pos, orientation);
        }

        // Try RC position (k - m - pos_in_kmer)
        let pos_in_kmer_rc = K as u64 - self.m() as u64 - minimizer_info.pos_in_kmer as u64;
        if let Some((pos, orientation)) = self.try_canonical_lookup_at_pos::<K>(
            query_kmer, kmer_rc, minimizer_pos, pos_in_kmer_rc
        ) {
            return (pos, orientation);
        }

        (INVALID_UINT64, 1)
    }

    /// Try canonical lookup at a specific position with a given pos_in_kmer
    #[inline]
    fn try_canonical_lookup_at_pos<const K: usize>(
        &self,
        query_kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        minimizer_pos: u64,
        pos_in_kmer: u64,
    ) -> Option<(u64, i8)>
    where
        Kmer<K>: KmerBits,
    {
        let kmer_pos = minimizer_pos.checked_sub(pos_in_kmer)?;

        let stored_kmer = self.spss.decode_kmer_at::<K>(kmer_pos as usize);

        let orientation = if stored_kmer.bits() == query_kmer.bits() {
            1i8
        } else if stored_kmer.bits() == kmer_rc.bits() {
            -1i8
        } else {
            return None;
        };

        if let Some((_sid, string_begin, string_end)) = self.spss.locate_with_end(kmer_pos) {
            if kmer_pos >= string_begin && kmer_pos < string_end - self.k as u64 + 1 {
                return Some((kmer_pos, orientation));
            }
        }
        None
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

    /// Check if canonical mode is enabled
    pub fn canonical(&self) -> bool {
        self.canonical
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
            self.canonical,
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
        let mut mphfs: Vec<Option<&crate::mphf_config::Mphf>> = Vec::with_capacity(
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
            canonical: header.canonical,
            hasher: crate::hasher::DeterministicHasher::new(1),
        })
    }

    /// Extract the minimizer from a k-mer using the cached hasher.
    #[inline]
    pub(crate) fn extract_minimizer<const K: usize>(&self, kmer: &Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        self.extract_minimizer_inline::<K>(*kmer)
    }

    /// Inline minimizer extraction without creating a MinimizerIterator.
    #[inline]
    fn extract_minimizer_inline<const K: usize>(&self, kmer: Kmer<K>) -> MinimizerInfo
    where
        Kmer<K>: KmerBits,
    {
        let num_windows = self.k - self.m + 1;
        let mask = (1u64 << (self.m * 2)) - 1;
        let bits = <Kmer<K> as KmerBits>::to_u128(kmer.bits());

        let mut min_hash = u64::MAX;
        let mut min_value = 0u64;
        let mut min_pos = 0usize;

        for i in 0..num_windows {
            let mmer_value = ((bits >> (i * 2)) as u64) & mask;
            let hash = self.hasher.hash_u64(mmer_value);
            if hash < min_hash {
                min_hash = hash;
                min_value = mmer_value;
                min_pos = i;
            }
        }

        MinimizerInfo::new(min_value, min_pos as u64, min_pos)
    }

    // --- Streaming query helpers ---

    /// Regular (non-canonical) lookup with pre-computed minimizer, returning LookupResult.
    pub(crate) fn lookup_regular_streaming<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        mini_info: MinimizerInfo,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(mini_info.value) {
            Some(id) => id,
            None => {
                let mut res = crate::streaming_query::LookupResult::not_found();
                res.minimizer_found = false;
                return res;
            }
        };

        if bucket_id >= self.index.num_buckets() {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        }

        let kmer_offset = self.lookup_in_bucket::<K>(kmer, mini_info, bucket_id);

        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        self.build_lookup_result(kmer_offset, 1)
    }

    /// Canonical lookup with pre-computed minimizer, returning LookupResult.
    pub(crate) fn lookup_canonical_streaming<const K: usize>(
        &self,
        kmer: &Kmer<K>,
        kmer_rc: &Kmer<K>,
        mini_info: MinimizerInfo,
    ) -> crate::streaming_query::LookupResult
    where
        Kmer<K>: KmerBits,
    {
        let bucket_id = match self.control_map.lookup(mini_info.value) {
            Some(id) => id,
            None => {
                let mut res = crate::streaming_query::LookupResult::not_found();
                res.minimizer_found = false;
                return res;
            }
        };

        if bucket_id >= self.index.num_buckets() {
            let mut res = crate::streaming_query::LookupResult::not_found();
            res.minimizer_found = false;
            return res;
        }

        let (kmer_offset, orientation) = self.lookup_in_bucket_canonical::<K>(
            kmer, kmer_rc, mini_info, bucket_id,
        );

        if kmer_offset == INVALID_UINT64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        self.build_lookup_result(kmer_offset, orientation)
    }

    /// Build a full LookupResult from a kmer_offset and orientation.
    #[inline]
    fn build_lookup_result(
        &self,
        kmer_offset: u64,
        orientation: i8,
    ) -> crate::streaming_query::LookupResult {
        let (string_id, string_begin, string_end) = match self.spss.locate_with_end(kmer_offset) {
            Some(triple) => triple,
            None => return crate::streaming_query::LookupResult::not_found(),
        };

        if kmer_offset > string_end - self.k as u64 {
            return crate::streaming_query::LookupResult::not_found();
        }

        let kmer_id_in_string = kmer_offset - string_begin;
        let kmer_id = kmer_offset - string_id * (self.k as u64 - 1);

        crate::streaming_query::LookupResult {
            kmer_id,
            kmer_id_in_string,
            kmer_offset,
            kmer_orientation: orientation,
            string_id,
            string_begin,
            string_end,
            minimizer_found: true,
        }
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

        let dict = Dictionary::new(spss, control_map, index, 31, 13, false);

        assert_eq!(dict.k(), 31);
        assert_eq!(dict.m(), 13);
        assert!(!dict.canonical());
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
        let mini_info = mini_iter.next(kmer);
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
}
