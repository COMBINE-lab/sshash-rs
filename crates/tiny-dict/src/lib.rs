//! TinyDictionary — a fast-path k-mer dictionary for small references.
//!
//! Drop-in alternative to [`sshash_lib::Dictionary`] for references where
//! RAM is plentiful but speed matters. Trades sshash's 5 bits/k-mer
//! minimizer + partitioned-MPHF + rank/select pipeline for a direct
//! hashbrown hash-table lookup — roughly 20 bytes per canonical k-mer,
//! but ~2–3× faster on the hot path because each lookup is a single
//! hash + probe.
//!
//! # Scope
//!
//! - Canonical-mode only (piscem's default).
//! - K ≤ 31 only for step 2 (canonical k-mer fits in a u64). Larger K
//!   needs a u128-keyed variant.
//! - Hashbrown lookup only — no streaming fast path yet. Consecutive-k-mer
//!   extension along a unitig will come in a later step.
//!
//! # Invariants
//!
//! The SPSS guarantees every canonical k-mer appears exactly once, so the
//! hash table needs a single entry per k-mer with no overflow. See the plan
//! notes for the reconstruction of [`LookupResult`] fields from the stored
//! `(string_id, kmer_id_in_string, spss_is_canonical)`.

use hashbrown::HashMap;
use rapidhash::fast::SeedableState;
use sshash_lib::spectrum_preserving_string_set::SpectrumPreservingStringSet;
use sshash_lib::{
    Dictionary, Kmer, KmerBits, KmerDictionary, KmerStreamingQuery, LookupResult,
};

pub mod io;
mod prefilter;

use prefilter::BlockedBloom;

/// BuildHasher for the TinyDictionary hashmap.
///
/// We use `rapidhash::fast::SeedableState::fixed()` rather than the default
/// `RandomState` so the hash of any given key is byte-for-byte identical
/// across runs of the same binary. Save/load relies on this: precomputed
/// hashes stored in `.tdct` files must re-resolve to the same slots on
/// reload. See `io::RAPIDHASH_VERSION` and the load-time probe check for
/// the safeguards that detect cross-version / cross-platform drift.
pub type TinyBuildHasher = SeedableState<'static>;

#[inline]
pub(crate) fn tiny_build_hasher() -> TinyBuildHasher {
    SeedableState::fixed()
}

/// Owned lightweight SPSS used inside [`TinyDictionary`].
///
/// Distinct from sshash's [`SpectrumPreservingStringSet`] in that we store
/// per-string offsets as a plain `Vec<u64>` rather than Elias-Fano. For
/// small references the few MB saved by EF is not worth the access cost.
pub struct TinySpss {
    /// 2-bit-packed DNA, identical byte layout to sshash's SPSS.
    strings: Vec<u8>,
    /// `offsets[i]` = base-index start of string `i`;
    /// `offsets[i+1]` = base-index end (exclusive). Length = `num_strings + 1`.
    offsets: Vec<u64>,
}

impl TinySpss {
    /// Create by copying out of an sshash SPSS.
    pub fn from_sshash(spss: &SpectrumPreservingStringSet) -> Self {
        let strings = spss.strings_raw().to_vec();
        let n = spss.offsets_len();
        let mut offsets = Vec::with_capacity(n);
        for i in 0..n {
            offsets.push(spss.offset_at(i));
        }
        Self { strings, offsets }
    }

    #[inline]
    pub fn num_strings(&self) -> u64 {
        if self.offsets.is_empty() {
            0
        } else {
            (self.offsets.len() - 1) as u64
        }
    }

    #[inline]
    pub fn string_offsets(&self, string_id: u64) -> (u64, u64) {
        let i = string_id as usize;
        (self.offsets[i], self.offsets[i + 1])
    }

    /// Decode the k-mer at `absolute_pos` bases into the concatenated SPSS.
    ///
    /// Port of [`SpectrumPreservingStringSet::decode_kmer_at`] — word-level
    /// unaligned load + shift + mask. Only the `K ≤ 31` fast path is needed
    /// here because TinyDictionary constrains K to ≤ 31.
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
            let raw = if byte_offset + 8 <= self.strings.len() {
                // SAFETY: bounds-checked above; strings are little-endian 2-bit
                // encoding so native-endian u64 read on LE hosts is correct.
                unsafe {
                    std::ptr::read_unaligned(self.strings.as_ptr().add(byte_offset) as *const u64)
                }
            } else {
                let mut buf = [0u8; 8];
                let avail = self.strings.len() - byte_offset;
                buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                u64::from_le_bytes(buf)
            };
            let shifted = raw >> bit_shift;
            let mask = if needed_bits >= 64 {
                u64::MAX
            } else {
                (1u64 << needed_bits) - 1
            };
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(shifted & mask))
        } else {
            // K=31 (and K=30) at bit_shift ∈ {4, 6} need 9 bytes: 62+6 = 68
            // bits exceeds a u64. Load u128 and mask; value still fits in u64
            // because K ≤ 31 caps needed_bits at 62.
            let raw: u128 = if byte_offset + 16 <= self.strings.len() {
                unsafe {
                    std::ptr::read_unaligned(self.strings.as_ptr().add(byte_offset) as *const u128)
                }
            } else {
                let mut buf = [0u8; 16];
                let avail = self.strings.len() - byte_offset;
                buf[..avail].copy_from_slice(&self.strings[byte_offset..byte_offset + avail]);
                u128::from_le_bytes(buf)
            };
            let shifted = raw >> bit_shift;
            let mask = (1u128 << needed_bits) - 1;
            Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64((shifted & mask) as u64))
        }
    }
}

/// Packed hashbrown value: `string_id(31) | spss_is_canonical(1) | kmer_id_in_string(32)`.
///
/// - `string_id`: ≤ 2^31 unitigs (plenty; typical Flex probe set is ~10^5).
/// - `spss_is_canonical`: bit indicating whether the forward k-mer at
///   `(string_id, kmer_id_in_string)` equals the canonical k-mer. Used to
///   reconstruct `kmer_orientation` at query time without decoding.
/// - `kmer_id_in_string`: ≤ 2^32 bases per string (more than enough).
pub(crate) type PackedValue = u64;

const SPSS_CANONICAL_BIT: u64 = 1u64 << 32;
const KMER_ID_MASK: u64 = (1u64 << 32) - 1;

#[inline]
fn pack(string_id: u32, kmer_id_in_string: u32, spss_is_canonical: bool) -> PackedValue {
    let mut v = (string_id as u64) << 33;
    if spss_is_canonical {
        v |= SPSS_CANONICAL_BIT;
    }
    v |= kmer_id_in_string as u64;
    v
}

#[inline]
fn unpack_string_id(v: PackedValue) -> u32 {
    (v >> 33) as u32
}

#[inline]
fn unpack_spss_is_canonical(v: PackedValue) -> bool {
    (v & SPSS_CANONICAL_BIT) != 0
}

#[inline]
fn unpack_kmer_id(v: PackedValue) -> u32 {
    (v & KMER_ID_MASK) as u32
}

/// Fast hashbrown-backed k-mer dictionary.
pub struct TinyDictionary {
    spss: TinySpss,
    k: usize,
    m: usize,
    canonical: bool,
    /// canonical k-mer bits (u64) → packed (string_id, kmer_id_in_string, spss_is_canonical)
    index: HashMap<u64, PackedValue, TinyBuildHasher>,
    /// Negative prefilter over `index`'s key set. Answers the common
    /// "definitely absent" case without touching `index` at all. Derived, not
    /// serialized — see [`prefilter`].
    bloom: BlockedBloom,
}

impl TinyDictionary {
    /// Build a TinyDictionary from an existing sshash [`Dictionary`].
    ///
    /// Iterates every k-mer in the SPSS, canonicalizes it, and inserts into
    /// the hashbrown map. Panics if the source is non-canonical (not yet
    /// supported) or if `K > 31`.
    pub fn from_sshash<const K: usize>(dict: &Dictionary) -> Self
    where
        Kmer<K>: KmerBits,
    {
        assert!(K <= 31, "TinyDictionary currently supports K ≤ 31 only");
        assert!(
            dict.canonical(),
            "TinyDictionary currently supports canonical-mode indexes only"
        );
        assert_eq!(dict.k(), K, "dict.k() must match const generic K");

        let spss = TinySpss::from_sshash(dict.spss());
        let num_strings = spss.num_strings();

        // Pre-size the hashtable: number of k-mers ≈ total_bases - num_strings*(k-1).
        // Cheap upper bound is total_bases; use the exact computation.
        let total_bases = if let Some(&last) = spss.offsets.last() {
            last
        } else {
            0
        };
        let total_kmers = total_bases.saturating_sub(num_strings * (K as u64 - 1));
        let mut index: HashMap<u64, PackedValue, TinyBuildHasher> =
            HashMap::with_capacity_and_hasher(total_kmers as usize, tiny_build_hasher());

        for string_id in 0..num_strings {
            let (begin, end) = spss.string_offsets(string_id);
            let slen = (end - begin) as usize;
            if slen < K {
                continue;
            }
            // Iterate every k-mer in this string.
            for kmer_pos in 0..=(slen - K) {
                let abs_pos = begin as usize + kmer_pos;
                let fw = spss.decode_kmer_at::<K>(abs_pos);
                let canonical = fw.canonical();
                let fw_bits_u64 = <Kmer<K> as KmerBits>::to_u64(fw.bits());
                let canon_bits_u64 = <Kmer<K> as KmerBits>::to_u64(canonical.bits());
                let spss_is_canonical = fw_bits_u64 == canon_bits_u64;
                let packed = pack(string_id as u32, kmer_pos as u32, spss_is_canonical);
                // SPSS invariant guarantees no duplicate canonical k-mers, so
                // any overwrite would indicate a corrupt input.
                let prev = index.insert(canon_bits_u64, packed);
                debug_assert!(
                    prev.is_none(),
                    "SPSS invariant violated: duplicate canonical k-mer"
                );
            }
        }

        let bloom = BlockedBloom::build(index.keys().copied(), index.len());

        Self {
            spss,
            k: K,
            m: dict.m(),
            canonical: dict.canonical(),
            index,
            bloom,
        }
    }

    /// Number of k-mer entries in the hash table.
    #[inline]
    pub fn len(&self) -> usize {
        self.index.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Internal lookup: returns `(string_id, kmer_id_in_string, orientation)`
    /// if found, else `None`. `orientation` is `+1` when the forward query
    /// k-mer matches the SPSS-stored forward k-mer, `-1` when it matches the
    /// reverse complement.
    #[inline]
    fn lookup_core<const K: usize>(&self, kmer_bytes: &[u8]) -> Option<(u32, u32, i64)>
    where
        Kmer<K>: KmerBits,
    {
        let fw = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
        let canonical = fw.canonical();
        let fw_bits = <Kmer<K> as KmerBits>::to_u64(fw.bits());
        let canon_bits = <Kmer<K> as KmerBits>::to_u64(canonical.bits());
        let query_fw_is_canonical = fw_bits == canon_bits;
        self.lookup_core_bits(canon_bits, query_fw_is_canonical)
    }

    /// Bit-level core lookup: caller has already parsed and canonicalized.
    #[inline]
    fn lookup_core_bits(&self, canon_bits: u64, query_fw_is_canonical: bool) -> Option<(u32, u32, i64)> {
        // Negative prefilter: no false negatives, so a `false` here is proof of
        // absence and skips the hashbrown probe entirely. On probe-panel
        // references ~79% of queried k-mers are absent, and this is where that
        // is turned into ~12 instructions instead of a full group probe.
        if !self.bloom.might_contain(canon_bits) {
            return None;
        }
        let packed = *self.index.get(&canon_bits)?;
        let string_id = unpack_string_id(packed);
        let kmer_id = unpack_kmer_id(packed);
        let spss_is_canonical = unpack_spss_is_canonical(packed);
        let orientation: i64 = if spss_is_canonical == query_fw_is_canonical {
            1
        } else {
            -1
        };
        Some((string_id, kmer_id, orientation))
    }

    /// One-shot lookup producing a full [`LookupResult`] identical to what
    /// sshash's [`Dictionary::query`] would return (modulo the minimizer-found
    /// flag, which TinyDictionary always reports as `true`).
    #[inline]
    pub fn lookup<const K: usize>(&self, kmer_bytes: &[u8]) -> LookupResult
    where
        Kmer<K>: KmerBits,
    {
        self.finish_lookup::<K>(self.lookup_core::<K>(kmer_bytes))
    }

    /// Same as [`Self::lookup`] but caller provides the pre-parsed canonical
    /// k-mer bits. Skips the ASCII→2-bit parse entirely.
    #[inline]
    pub fn lookup_bits<const K: usize>(
        &self,
        canonical_bits: u64,
        fw_is_canonical: bool,
    ) -> LookupResult
    where
        Kmer<K>: KmerBits,
    {
        self.finish_lookup::<K>(self.lookup_core_bits(canonical_bits, fw_is_canonical))
    }

    #[inline]
    fn finish_lookup<const K: usize>(&self, core: Option<(u32, u32, i64)>) -> LookupResult
    where
        Kmer<K>: KmerBits,
    {
        match core {
            Some((string_id, kmer_id_in_string, orientation)) => {
                let (string_begin, string_end) = self.spss.string_offsets(string_id as u64);
                let kmer_offset = string_begin + kmer_id_in_string as u64;
                let kmer_id = kmer_offset - (string_id as u64) * (K as u64 - 1);
                LookupResult {
                    kmer_id,
                    kmer_id_in_string: kmer_id_in_string as u64,
                    kmer_offset,
                    kmer_orientation: orientation,
                    string_id: string_id as u64,
                    string_begin,
                    string_end,
                    minimizer_found: true,
                }
            }
            None => LookupResult::not_found(),
        }
    }
}

/// Anchor state carried between consecutive-k-mer lookups.
///
/// When the previous lookup resolved to `(string_id, kmer_id_in_string)`
/// with orientation `orientation`, the next consecutive read k-mer is
/// *expected* to live at the adjacent SPSS position (forward when
/// `orientation = +1`, backward when `-1`). We can verify that in-place
/// by decoding the stored SPSS k-mer at the adjacent position and
/// comparing its bits against the read k-mer — no hashing, no hashmap
/// probe, just two u64 loads and a compare.
#[derive(Clone, Copy)]
struct Anchor {
    string_id: u32,
    kmer_pos_in_string: u32,
    string_begin: u64,
    string_end: u64,
    /// +1 when SPSS-fw at `(string_id, kmer_pos_in_string)` equals query-fw,
    /// -1 when SPSS-fw equals query-rc.
    orientation: i64,
}

/// Streaming query engine for [`TinyDictionary`] with single-unitig
/// extension fast path.
pub struct TinyStreamingQuery<'a, const K: usize>
where
    Kmer<K>: KmerBits,
{
    dict: &'a TinyDictionary,
    anchor: Option<Anchor>,
    /// Total lookups attempted (sum of fast- and slow-path).
    num_lookups: u64,
    /// Full hashbrown probes (slow path).
    num_searches: u64,
    /// Successful anchor-extension hits (fast path — no hashing).
    num_extensions: u64,
}

impl<'a, const K: usize> TinyStreamingQuery<'a, K>
where
    Kmer<K>: KmerBits,
{
    #[inline]
    pub fn new(dict: &'a TinyDictionary) -> Self {
        Self {
            dict,
            anchor: None,
            num_lookups: 0,
            num_searches: 0,
            num_extensions: 0,
        }
    }

    /// Try to satisfy a lookup by extending the current anchor along its
    /// SPSS unitig. The caller provides the pre-parsed forward read k-mer.
    /// Returns `Some(result)` on success; `None` means the caller must fall
    /// back to the hashmap.
    #[inline]
    fn try_extend_fw(&mut self, read_fw: Kmer<K>) -> Option<LookupResult> {
        let anchor = self.anchor.as_mut()?;
        let k_u64 = K as u64;
        let string_len = anchor.string_end - anchor.string_begin;

        // Adjacent SPSS position in the direction implied by orientation.
        let new_pos: u64 = if anchor.orientation > 0 {
            let np = anchor.kmer_pos_in_string as u64 + 1;
            if np + k_u64 > string_len {
                return None;
            }
            np
        } else {
            if anchor.kmer_pos_in_string == 0 {
                return None;
            }
            (anchor.kmer_pos_in_string - 1) as u64
        };

        let abs_pos = anchor.string_begin + new_pos;
        let spss_fw = self.dict.spss.decode_kmer_at::<K>(abs_pos as usize);

        // Expected read-fw bits given SPSS-fw at new_pos and the anchor's orientation.
        let expected_fw = if anchor.orientation > 0 {
            spss_fw
        } else {
            spss_fw.reverse_complement()
        };
        if <Kmer<K> as KmerBits>::to_u64(read_fw.bits())
            != <Kmer<K> as KmerBits>::to_u64(expected_fw.bits())
        {
            return None;
        }

        anchor.kmer_pos_in_string = new_pos as u32;
        let kmer_offset = anchor.string_begin + new_pos;
        let kmer_id = kmer_offset - (anchor.string_id as u64) * (k_u64 - 1);
        Some(LookupResult {
            kmer_id,
            kmer_id_in_string: new_pos,
            kmer_offset,
            kmer_orientation: anchor.orientation,
            string_id: anchor.string_id as u64,
            string_begin: anchor.string_begin,
            string_end: anchor.string_end,
            minimizer_found: true,
        })
    }

    #[inline]
    fn set_anchor_from_result(&mut self, result: &LookupResult) {
        if result.is_found() {
            self.anchor = Some(Anchor {
                string_id: result.string_id as u32,
                kmer_pos_in_string: result.kmer_id_in_string as u32,
                string_begin: result.string_begin,
                string_end: result.string_end,
                orientation: result.kmer_orientation,
            });
        } else {
            self.anchor = None;
        }
    }
}

impl<'a, const K: usize> KmerStreamingQuery for TinyStreamingQuery<'a, K>
where
    Kmer<K>: KmerBits,
{
    const PREFERS_BITS: bool = true;

    #[inline]
    fn reset(&mut self) {
        self.anchor = None;
    }

    #[inline]
    fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        // Parse once, reuse for fast path and slow path.
        let read_fw = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
        self.num_lookups += 1;
        if let Some(result) = self.try_extend_fw(read_fw) {
            self.num_extensions += 1;
            return result;
        }
        self.num_searches += 1;
        let canonical = read_fw.canonical();
        let fw_bits = <Kmer<K> as KmerBits>::to_u64(read_fw.bits());
        let canon_bits = <Kmer<K> as KmerBits>::to_u64(canonical.bits());
        let fw_is_canonical = fw_bits == canon_bits;
        let result = self.dict.lookup_bits::<K>(canon_bits, fw_is_canonical);
        self.set_anchor_from_result(&result);
        result
    }

    #[inline]
    fn lookup_bits(
        &mut self,
        canonical_bits: u64,
        fw_is_canonical: bool,
        _fw_bytes: &[u8],
    ) -> LookupResult {
        // Caller has already parsed+canonicalized the k-mer. If we have an
        // anchor and the fast-path fits, reconstruct read_fw from bits and try;
        // otherwise, skip straight to the hashbrown probe using the bits.
        self.num_lookups += 1;
        if self.anchor.is_some() {
            // Reconstruct FW k-mer to call try_extend_fw.
            let canonical = Kmer::<K>::new(<Kmer<K> as KmerBits>::from_u64(canonical_bits));
            let read_fw = if fw_is_canonical { canonical } else { canonical.reverse_complement() };
            if let Some(result) = self.try_extend_fw(read_fw) {
                self.num_extensions += 1;
                return result;
            }
        }
        self.num_searches += 1;
        let result = self.dict.lookup_bits::<K>(canonical_bits, fw_is_canonical);
        self.set_anchor_from_result(&result);
        result
    }

    #[inline]
    fn num_searches(&self) -> u64 {
        self.num_searches
    }

    #[inline]
    fn num_extensions(&self) -> u64 {
        self.num_extensions
    }

    #[inline]
    fn skip_anchor_along_string(&mut self, read_offset: i32) -> bool {
        // Shift the live anchor by `read_offset` read-positions. Sequence
        // agreement must have been verified by the caller (piscem's
        // check_direct_match). The SPSS offset direction follows the anchor's
        // orientation.
        let anchor = match &mut self.anchor {
            Some(a) => a,
            None => return false,
        };
        let k_u64 = K as u64;
        let string_len = anchor.string_end - anchor.string_begin;
        let signed_spss_off = if anchor.orientation > 0 {
            read_offset as i64
        } else {
            -(read_offset as i64)
        };
        let new_pos = anchor.kmer_pos_in_string as i64 + signed_spss_off;
        if new_pos < 0 {
            self.anchor = None;
            return false;
        }
        let new_pos_u = new_pos as u64;
        if new_pos_u + k_u64 > string_len {
            self.anchor = None;
            return false;
        }
        anchor.kmer_pos_in_string = new_pos_u as u32;
        true
    }
}

impl KmerDictionary for TinyDictionary {
    type Query<'a, const K: usize>
        = TinyStreamingQuery<'a, K>
    where
        Self: 'a,
        Kmer<K>: KmerBits;

    #[inline]
    fn k(&self) -> usize {
        self.k
    }

    #[inline]
    fn m(&self) -> usize {
        self.m
    }

    #[inline]
    fn num_strings(&self) -> u64 {
        self.spss.num_strings()
    }

    #[inline]
    fn canonical(&self) -> bool {
        self.canonical
    }

    #[inline]
    fn kmer_at_pos<const K: usize>(&self, absolute_base_pos: usize) -> Kmer<K>
    where
        Kmer<K>: KmerBits,
    {
        self.spss.decode_kmer_at::<K>(absolute_base_pos)
    }

    #[inline]
    fn create_streaming_query<const K: usize>(&self) -> Self::Query<'_, K>
    where
        Kmer<K>: KmerBits,
    {
        TinyStreamingQuery::new(self)
    }
}
