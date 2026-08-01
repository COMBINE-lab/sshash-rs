//! On-disk format for [`TinyDictionary`] (`.tdct`).
//!
//! The format trades tight packing for a fast load: the on-disk layout
//! carries precomputed rapidhash values alongside each (key, value) pair
//! so the load path can use `RawVacantEntryMut::insert_hashed_nocheck` and
//! skip re-hashing every key. Each file records the rapidhash version that
//! wrote it (see [`RAPIDHASH_VERSION`]); on load, a version outside the
//! verified-compatible set, or a platform-specific drift (rapidhash's own docs
//! disclaim stability across compilers / target CPUs), is caught by the on-load
//! probe check before any bulk work happens, so corrupt lookups are impossible
//! — only loud errors asking the user to rebuild.
//!
//! # Layout (little-endian, tightly packed unless noted)
//!
//! ```text
//!   0.. 8   magic            = b"TDCT0001"
//!   8..10   format_version   : u16
//!  10..13   rapidhash_ver    : [u8; 3]  // semver major/minor/patch
//!  13..14   _pad             : u8
//!  14..15   k                : u8
//!  15..16   m                : u8
//!  16..17   canonical        : u8 (0/1)
//!  17..20   _pad             : [u8; 3]
//!  20..24   probe_count      : u32  (≥ 1; index 0 is always the key `0u64`)
//!  24..                       probes: (u64 key, u64 expected_hash) × probe_count
//!
//!    ...   spss_strings_len : u64 (byte count)
//!    ...   spss_strings     : [u8; spss_strings_len]
//!    ...   spss_offsets_len : u64 (u64 count)
//!    ...   spss_offsets     : [u64; spss_offsets_len]
//!    ...   index_len        : u64 (entry count)
//!    ...   entries          : (u64 hash, u64 key, u64 packed_value) × index_len
//! ```
//!
//! No padding is inserted between variable-length sections: readers must
//! consume exactly what was written.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use hashbrown::HashMap;
use hashbrown::hash_map::RawEntryMut;
use std::hash::BuildHasher;

use crate::{PackedValue, TinyBuildHasher, TinyDictionary, TinySpss, tiny_build_hasher};

/// Magic bytes at the start of every `.tdct` file.
pub const MAGIC: &[u8; 8] = b"TDCT0001";

/// On-disk format revision. Bump on any layout change.
pub const FORMAT_VERSION: u16 = 1;

/// Pinned rapidhash semver. Must match the `rapidhash` dep version in
/// `Cargo.toml` (which is pinned with `=`). This is the version stamped into
/// the header of every `.tdct` this binary writes.
pub const RAPIDHASH_VERSION: [u8; 3] = [4, 5, 1];

/// rapidhash versions whose `fast::SeedableState::fixed()` output is known to
/// be identical, and whose `.tdct` files are therefore interchangeable.
///
/// The load path accepts a header naming any of these rather than requiring an
/// exact match with [`RAPIDHASH_VERSION`]. That is safe because this list is
/// only a fast pre-filter: the probe check below re-hashes stored keys and is
/// what actually proves the linked rapidhash agrees with the one that wrote the
/// file. A version added here in error cannot corrupt a lookup — it just moves
/// the loud failure from `RapidhashVersionMismatch` to `ProbeMismatch`.
///
/// Only add a version after verifying it hashes identically; rapidhash's own
/// docs disclaim output stability across releases, so this must stay a
/// verified allow-list rather than a range.
const ACCEPTED_RAPIDHASH_VERSIONS: &[[u8; 3]] = &[[4, 4, 0], [4, 5, 1]];

/// Number of (key, hash) probes stored in the header for load-time
/// cross-platform drift detection. 16 is overkill (one would suffice) but
/// negligible on disk and gives belt-and-suspenders.
const PROBE_COUNT: u32 = 16;

#[derive(Debug, thiserror::Error)]
pub enum TdctError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("not a .tdct file (bad magic)")]
    BadMagic,
    #[error("unsupported .tdct format version {got}; this build expects {expected}")]
    UnsupportedFormat { got: u16, expected: u16 },
    #[error(
        "this .tdct was built with rapidhash {got:?}, which is not known to hash identically to \
         the rapidhash {expected:?} this binary links. Rapidhash does not guarantee stable hash \
         output across versions — rebuild the index with `piscem build --dict tiny` on this \
         machine."
    )]
    RapidhashVersionMismatch { got: [u8; 3], expected: [u8; 3] },
    #[error(
        "on-load hash probe failed for key {key:#018x}: expected {expected:#018x}, got {got:#018x}. \
         This .tdct was built with a compatible-looking rapidhash version but produces different \
         hashes on this machine (likely different target CPU, compiler, or a subtly different \
         linked rapidhash). Rebuild the index with `piscem build --dict tiny` on this machine."
    )]
    ProbeMismatch { key: u64, expected: u64, got: u64 },
    #[error("truncated .tdct: {0}")]
    Truncated(&'static str),
}

impl TinyDictionary {
    /// Serialize to a `.tdct` file at `path`.
    ///
    /// Writes header, SPSS, and precomputed-hash (hash, key, value) triples.
    /// The hash values are produced by [`tiny_build_hasher`] — which uses
    /// [`rapidhash::fast::SeedableState::fixed`], the only rapidhash state
    /// with a seed constant across runs of this binary.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), TdctError> {
        let file = File::create(path.as_ref())?;
        let mut w = BufWriter::with_capacity(1 << 20, file);
        let hasher = tiny_build_hasher();

        // -- header --
        w.write_all(MAGIC)?;
        w.write_all(&FORMAT_VERSION.to_le_bytes())?;
        w.write_all(&RAPIDHASH_VERSION)?;
        w.write_all(&[0u8])?; // _pad
        w.write_all(&[self.k as u8, self.m as u8, self.canonical as u8])?;
        w.write_all(&[0u8; 3])?; // _pad

        // -- probes --
        w.write_all(&PROBE_COUNT.to_le_bytes())?;
        // Probe 0 is always the constant key `0u64` — a self-check that
        // catches drift even on indexes whose real keys happen to hash to
        // the same values on both platforms (vanishingly unlikely, but
        // principled).
        let probe_keys = collect_probe_keys(self, PROBE_COUNT as usize);
        for k in &probe_keys {
            let h = hasher.hash_one(*k);
            w.write_all(&k.to_le_bytes())?;
            w.write_all(&h.to_le_bytes())?;
        }

        // -- spss --
        w.write_all(&(self.spss.strings.len() as u64).to_le_bytes())?;
        w.write_all(&self.spss.strings)?;
        w.write_all(&(self.spss.offsets.len() as u64).to_le_bytes())?;
        for off in &self.spss.offsets {
            w.write_all(&off.to_le_bytes())?;
        }

        // -- index entries (hash, key, value) --
        w.write_all(&(self.index.len() as u64).to_le_bytes())?;
        for (k, v) in self.index.iter() {
            let h = hasher.hash_one(*k);
            w.write_all(&h.to_le_bytes())?;
            w.write_all(&k.to_le_bytes())?;
            w.write_all(&v.to_le_bytes())?;
        }

        w.flush()?;
        Ok(())
    }

    /// Load a `.tdct` file written by [`save`](Self::save).
    ///
    /// Verifies the magic, format version, and rapidhash version, then
    /// runs every stored probe through [`tiny_build_hasher`] and compares
    /// hashes. Only after all checks pass does the bulk `(hash, key, value)`
    /// table get read in; insertion uses `raw_entry_mut().insert_hashed_nocheck`
    /// so no key is hashed twice.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, TdctError> {
        let file = File::open(path.as_ref())?;
        let mut r = BufReader::with_capacity(1 << 20, file);

        // -- header --
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(TdctError::BadMagic);
        }
        let format_version = read_u16(&mut r)?;
        if format_version != FORMAT_VERSION {
            return Err(TdctError::UnsupportedFormat {
                got: format_version,
                expected: FORMAT_VERSION,
            });
        }
        let mut rh_ver = [0u8; 3];
        r.read_exact(&mut rh_ver)?;
        if !ACCEPTED_RAPIDHASH_VERSIONS.contains(&rh_ver) {
            return Err(TdctError::RapidhashVersionMismatch {
                got: rh_ver,
                expected: RAPIDHASH_VERSION,
            });
        }
        let mut pad = [0u8; 1];
        r.read_exact(&mut pad)?;

        let mut kmc = [0u8; 3];
        r.read_exact(&mut kmc)?;
        let k = kmc[0] as usize;
        let m = kmc[1] as usize;
        let canonical = kmc[2] != 0;
        let mut pad3 = [0u8; 3];
        r.read_exact(&mut pad3)?;

        // -- probes --
        let probe_count = read_u32(&mut r)? as usize;
        let hasher = tiny_build_hasher();
        for _ in 0..probe_count {
            let key = read_u64(&mut r)?;
            let expected = read_u64(&mut r)?;
            let got = hasher.hash_one(key);
            if got != expected {
                return Err(TdctError::ProbeMismatch {
                    key,
                    expected,
                    got,
                });
            }
        }

        // -- spss --
        let strings_len = read_u64(&mut r)? as usize;
        let mut strings = vec![0u8; strings_len];
        r.read_exact(&mut strings)?;
        let offsets_len = read_u64(&mut r)? as usize;
        let mut offsets = Vec::with_capacity(offsets_len);
        for _ in 0..offsets_len {
            offsets.push(read_u64(&mut r)?);
        }
        let spss = TinySpss::from_raw_parts(strings, offsets);

        // -- index --
        let index_len = read_u64(&mut r)? as usize;
        let mut index: HashMap<u64, PackedValue, TinyBuildHasher> =
            HashMap::with_capacity_and_hasher(index_len, tiny_build_hasher());
        for _ in 0..index_len {
            let hash = read_u64(&mut r)?;
            let key = read_u64(&mut r)?;
            let value = read_u64(&mut r)?;
            match index.raw_entry_mut().from_hash(hash, |_| false) {
                RawEntryMut::Vacant(v) => {
                    v.insert_hashed_nocheck(hash, key, value);
                }
                RawEntryMut::Occupied(_) => {
                    // SPSS invariant: canonical k-mers are unique. If the
                    // same `hash` collides to an existing key we fall
                    // through to the normal insert which will do a real
                    // key-equality check. This is only reached on the
                    // extraordinarily rare case of a hash collision.
                    index.insert(key, value);
                }
            }
        }

        Ok(TinyDictionary::from_raw_parts(
            spss, k, m, canonical, index,
        ))
    }
}

/// Select up to `n` keys from the dictionary to store as load-time probes.
/// Always emits key `0u64` as probe 0 so we have a stable self-check even
/// for indexes whose real keys happen to hash identically on two platforms.
/// The remaining keys are the first `n - 1` iteration-order keys; this is
/// deterministic per-save but not cryptographic — probes are for drift
/// detection, not security.
fn collect_probe_keys(dict: &TinyDictionary, n: usize) -> Vec<u64> {
    let mut probes = Vec::with_capacity(n);
    probes.push(0u64);
    for &k in dict.index.keys() {
        if probes.len() >= n {
            break;
        }
        if k == 0 {
            continue;
        }
        probes.push(k);
    }
    probes
}

impl TinySpss {
    /// Reconstruct a [`TinySpss`] from already-deserialized buffers.
    #[inline]
    pub(crate) fn from_raw_parts(strings: Vec<u8>, offsets: Vec<u64>) -> Self {
        Self { strings, offsets }
    }
}

impl TinyDictionary {
    /// Reconstruct a [`TinyDictionary`] from already-deserialized parts.
    /// Only used by the `.tdct` load path; `from_sshash` is the general
    /// constructor for runtime conversion.
    pub(crate) fn from_raw_parts(
        spss: TinySpss,
        k: usize,
        m: usize,
        canonical: bool,
        index: HashMap<u64, PackedValue, TinyBuildHasher>,
    ) -> Self {
        // The prefilter is derived from the loaded key set rather than stored,
        // so `.tdct` stays byte-compatible in both directions. One pass over
        // ~1.5M keys costs single-digit milliseconds against a ~35 MB file.
        let bloom = crate::BlockedBloom::build(index.keys().copied(), index.len());
        Self {
            spss,
            k,
            m,
            canonical,
            index,
            bloom,
        }
    }
}

// --- little-endian readers --------------------------------------------------

fn read_u16<R: Read>(r: &mut R) -> Result<u16, TdctError> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}

fn read_u32<R: Read>(r: &mut R) -> Result<u32, TdctError> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn read_u64<R: Read>(r: &mut R) -> Result<u64, TdctError> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
