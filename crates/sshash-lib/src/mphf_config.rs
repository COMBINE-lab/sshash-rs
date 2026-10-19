//! MPHF (Minimal Perfect Hash Function) type configuration
//!
//! Central module for MPHF type aliases and helpers used across the crate.
//! We use sux's PHastR (PHast with rotated layouts), which matches PHast in
//! space (about 1.94 bits/key with the defaults) while being faster to build
//! and to query.
//!
//! Keys are hashed with rapidhash v3 (through our own `ToSig` implementation
//! rather than sux's default xxh3, which measured 5-15% slower per query).
//! The `v3` module is rapidhash's stable-output API and is CPU-feature
//! independent (unlike ahash, which switches algorithm when AES-NI is
//! available via target-cpu=native, breaking serialized indices).
//!
//! The function is key-type-erased: it is built over the native-endian bytes
//! of each key (see [`KeyBytes`]), so `u64` minimizers and `u64`/`u128` k-mers
//! share one MPHF type; queries inline the hash at a constant key length.

use dsi_progress_logger::no_logging;
use lender::{FallibleLend, FallibleLender, FallibleLending};
use std::borrow::Borrow;
use std::io;
use std::path::Path;
use std::sync::Mutex;
use sux::func::{PHastR, PHastRBuilder};
use sux::utils::ToSig;

/// Our MPHF type — PHastR with the default parameters (8-bit seeds stored in
/// a `Box<[u8]>`, four layouts, slices of length 1024, bucket size 4.25),
/// keyed by the native-endian bytes of the key (see [`MphfKey`]).
pub type Mphf = PHastR<KeyBytes>;

/// The bytes of a key, as hashed by the MPHF.
#[repr(transparent)]
pub struct KeyBytes([u8]);

impl KeyBytes {
    /// Views a byte slice as key bytes.
    #[inline(always)]
    pub fn new(bytes: &[u8]) -> &KeyBytes {
        // SAFETY: `KeyBytes` is a `repr(transparent)` wrapper around `[u8]`
        unsafe { &*(bytes as *const [u8] as *const KeyBytes) }
    }

    /// Returns the underlying bytes.
    #[inline(always)]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }
}

impl ToSig<[u64; 1]> for KeyBytes {
    #[inline(always)]
    fn to_sig(key: impl Borrow<Self>, seed: u64) -> [u64; 1] {
        [hash_key_bytes(&key.borrow().0, seed)]
    }
}

// The key type is part of the ε-serde type hash of `Mphf`: naming the hash
// scheme here makes a serialized MPHF fail to load (rather than silently
// misroute queries) if the hashing of keys ever changes.
impl epserde::traits::TypeHash for KeyBytes {
    fn type_hash(hasher: &mut impl core::hash::Hasher) {
        use core::hash::Hash;
        "sshash::KeyBytes(rapidhash_v3_inline<true,false,false>, seed_cpp)".hash(hasher);
    }
}

/// Seeded rapidhash v3 of the bytes of a key (also used for partitioning).
#[inline(always)]
pub(crate) fn hash_key_bytes(bytes: &[u8], seed: u64) -> u64 {
    use rapidhash::v3::{RapidSecrets, rapidhash_v3_inline};
    rapidhash_v3_inline::<true, false, false>(bytes, &RapidSecrets::seed_cpp(seed))
}

/// Owned key bytes of a fixed width, borrowable as [`KeyBytes`].
#[derive(Clone, Copy)]
pub struct KeyArray<const N: usize>([u8; N]);

impl<const N: usize> Borrow<KeyBytes> for KeyArray<N> {
    #[inline(always)]
    fn borrow(&self) -> &KeyBytes {
        KeyBytes::new(&self.0)
    }
}

/// A key type that can be stored in an [`Mphf`].
///
/// Implemented for the integer types used as keys (`u64` minimizers and
/// `u64`/`u128` k-mers). The MPHF hashes `key_bytes()`, which for these types
/// are its native-endian bytes.
pub trait MphfKey: Copy + Send + Sync {
    /// The byte representation the MPHF is built over.
    type Bytes: Borrow<KeyBytes> + Copy + Send + Sync;

    /// Returns the bytes hashed by the MPHF for this key.
    fn key_bytes(self) -> Self::Bytes;
}

impl MphfKey for u64 {
    type Bytes = KeyArray<8>;

    #[inline(always)]
    fn key_bytes(self) -> KeyArray<8> {
        KeyArray(self.to_ne_bytes())
    }
}

impl MphfKey for u128 {
    type Bytes = KeyArray<16>;

    #[inline(always)]
    fn key_bytes(self) -> KeyArray<16> {
        KeyArray(self.to_ne_bytes())
    }
}

/// Look up a key in an MPHF.
///
/// For keys in the build set the result is a distinct index in
/// [0, `mphf.len()`); for any other key it is arbitrary but still in range
/// (zero for an empty function).
#[inline(always)]
pub fn mphf_get<K: MphfKey>(mphf: &Mphf, key: K) -> usize {
    mphf.get(key.key_bytes().borrow())
}

/// Build an MPHF from the bytes of a set of distinct keys.
///
/// PHastR builds on the current rayon pool, and the function it builds
/// depends on the number of threads of that pool. We build on a dedicated
/// single-thread pool so that indices are byte-identical whatever
/// `num_threads` they are built with. To build many functions in parallel,
/// use [`build_mphfs`], not this function inside a parallel iterator (see
/// there why).
///
/// # Panics
///
/// If the keys contain duplicates (construction cannot succeed then).
pub fn build_mphf<B: Borrow<KeyBytes> + Sync>(keys: &[B]) -> Mphf {
    single_thread_pool().install(|| build_mphf_on_current_pool(keys))
}

/// Build one MPHF per key set, with at most `rayon::current_num_threads()`
/// constructions in flight, each on a single-thread pool (so the result is
/// the same as [`build_mphf`] on each set). Each key set is freed as soon as
/// its function is built.
///
/// Calling [`build_mphf`] from a rayon parallel iterator would not bound
/// concurrency: a worker that installs a job on another pool keeps stealing
/// work from its own pool while waiting, so every key set would start
/// construction at once, each holding ~10 bytes/key of working memory.
/// Here each of `current_num_threads()` plain threads owns one single-thread
/// pool and pulls key sets from a shared queue.
pub fn build_mphfs<B: Borrow<KeyBytes> + Send + Sync>(key_sets: Vec<Vec<B>>) -> Vec<Mphf> {
    let n = key_sets.len();
    let workers = rayon::current_num_threads().clamp(1, n.max(1));
    let queue = Mutex::new(key_sets.into_iter().enumerate());
    let results: Vec<Mutex<Option<Mphf>>> = (0..n).map(|_| Mutex::new(None)).collect();
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| {
                single_thread_pool().install(|| loop {
                    let next = queue.lock().unwrap().next();
                    let Some((i, keys)) = next else { break };
                    let mphf = build_mphf_on_current_pool(&keys);
                    drop(keys);
                    *results[i].lock().unwrap() = Some(mphf);
                })
            });
        }
    });
    results.into_iter().map(|m| m.into_inner().unwrap().expect("every key set is built")).collect()
}

/// Maximum number of chunks PHastR sweeps a level of a streamed MPHF in.
/// PHastR's output depends on the number of sweep chunks, which by default
/// is the size of the current rayon pool: fixing it keeps indices
/// byte-identical across `num_threads`, while construction still runs on
/// (and is bounded by) the caller's pool. Sweeps use at most this many
/// threads (178M keys: 3.1 s with 8 threads, 8.9 s with 1).
pub const STREAMING_SWEEP_CHUNKS: usize = 8;

/// Streamed construction keeps key signatures in temporary files (offline
/// mode) from this many keys on: offline mode uses 256 MiB of fixed buffers,
/// which is more than in-memory construction below ~16M keys.
pub const OFFLINE_MIN_KEYS: usize = 1 << 24;

/// How the minimizer MPHF of an index is built.
///
/// Every strategy produces a valid index in the same file format; they differ
/// in build time, build memory, temporary disk and query speed (and hence in
/// the exact bytes of the index). Each one is deterministic: the index does
/// not depend on `num_threads`.
///
/// Human genome scale (178M minimizers), MPHF construction alone:
///
/// | strategy | build, 1 / 8 / 32 threads | extra memory | temp disk | lookup |
/// |---|---|---|---|---|
/// | `Streaming(Offline::Auto)` | 8.7 / 3.0 / 3.0 s | 0.27 GB | 16 B/key | ~20 ns |
/// | `Streaming(Offline::Never)` | 8.0 / 2.6 / 2.6 s | ~27 B/key | none | ~20 ns |
/// | `InMemory` | 7.2 / 1.0 / 0.9 s | ~8 B/key + ~20 B/key working | none | ~20 ns |
/// | `Partitioned` | 7.7 / 1.4 / 0.7 s | ~8 B/key + 3M keys per thread | none | ~25 ns |
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MphfBuild {
    /// One function, built by streaming the (sorted, distinct) keys once
    /// without materializing them. Construction records are kept in memory or
    /// in temporary files in `BuildConfiguration::tmp_dirname` according to
    /// [`Offline`]. The default.
    Streaming(Offline),
    /// One function, built from an in-memory copy of the keys on the build
    /// pool. Fastest single-function build, no temporary disk.
    InMemory,
    /// Keys split into partitions of ~3M keys, built `num_threads` at a time.
    /// No temporary disk and the fastest build at high thread counts, but
    /// ~25% slower lookups and a slightly larger function.
    Partitioned,
}

impl Default for MphfBuild {
    fn default() -> Self {
        MphfBuild::Streaming(Offline::default())
    }
}

/// Where [`MphfBuild::Streaming`] keeps its construction records (16 bytes
/// per key).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Offline {
    /// In temporary files from `min_keys` keys on, in memory below. Offline
    /// construction has ~256 MiB of fixed buffers, so it only saves memory
    /// above ~16M keys. The default is `min_keys = OFFLINE_MIN_KEYS` (2^24).
    Auto {
        /// Smallest key count built offline.
        min_keys: usize,
    },
    /// Always in temporary files (least memory).
    Always,
    /// Always in memory (no temporary disk).
    Never,
}

impl Default for Offline {
    fn default() -> Self {
        Offline::Auto { min_keys: OFFLINE_MIN_KEYS }
    }
}

impl Offline {
    /// Whether a set of `num_keys` keys is built offline.
    pub fn is_offline(self, num_keys: usize) -> bool {
        match self {
            Offline::Auto { min_keys } => num_keys >= min_keys,
            Offline::Always => true,
            Offline::Never => false,
        }
    }
}

/// Lends the bytes of each key of an iterator, for PHastR's lender-based
/// constructors.
struct KeyLender<I: Iterator<Item: MphfKey>> {
    iter: I,
    cur: Option<<I::Item as MphfKey>::Bytes>,
}

impl<'lend, I: Iterator<Item: MphfKey>> FallibleLending<'lend> for KeyLender<I> {
    type Lend = &'lend <I::Item as MphfKey>::Bytes;
}

impl<I: Iterator<Item: MphfKey>> FallibleLender for KeyLender<I> {
    lender::check_covariance_fallible!();
    type Error = std::convert::Infallible;

    #[inline]
    fn next(&mut self) -> Result<Option<FallibleLend<'_, Self>>, Self::Error> {
        self.cur = self.iter.next().map(MphfKey::key_bytes);
        Ok(self.cur.as_ref())
    }
}

/// Build an MPHF by streaming distinct keys, which are read once and never
/// materialized.
///
/// With `offline`, the 16-byte per-key construction records go to unlinked
/// temporary files in `tmp_dir` (about 1.5 bytes/key of memory, 16 bytes/key
/// of disk); otherwise they are kept in memory. Construction runs on the
/// current rayon pool with [`STREAMING_SWEEP_CHUNKS`] sweep chunks, so the
/// result does not depend on the caller's thread count.
///
/// # Panics
///
/// If the keys contain duplicates, or offline construction fails (e.g., the
/// temporary directory is not writable).
pub fn build_mphf_streaming<I>(keys: I, offline: bool, tmp_dir: &Path) -> Mphf
where
    I: Iterator<Item: MphfKey> + Send,
{
    let mut builder = PHastRBuilder::default()
        .offline(offline)
        .sweep_chunks(STREAMING_SWEEP_CHUNKS);
    if offline {
        if let Err(e) = std::fs::create_dir_all(tmp_dir) {
            tracing::warn!("cannot create MPHF temporary directory {}: {e}", tmp_dir.display());
        }
        builder = builder.temp_dir(tmp_dir);
    }
    Mphf::try_new_with_builder(KeyLender { iter: keys, cur: None }, builder, no_logging![])
        .unwrap_or_else(|e| panic!("streamed MPHF construction failed: {e:#}"))
}

/// Build an MPHF from the bytes of a set of distinct keys on the current
/// rayon pool, with [`STREAMING_SWEEP_CHUNKS`] sweep chunks so that the
/// result does not depend on the pool's size (unlike [`build_mphf`], which
/// builds on a single thread for the same reason).
///
/// # Panics
///
/// If the keys contain duplicates.
pub fn build_mphf_parallel<B: Borrow<KeyBytes> + Sync>(keys: &[B]) -> Mphf {
    let builder = PHastRBuilder::default().sweep_chunks(STREAMING_SWEEP_CHUNKS);
    Mphf::try_par_new_with_builder(keys, builder, no_logging![])
        .unwrap_or_else(|e| panic!("MPHF construction failed on {} keys: {e:#}", keys.len()))
}

fn single_thread_pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("failed to create the MPHF construction thread pool")
}

/// Builds on the current pool, which must have a single thread.
fn build_mphf_on_current_pool<B: Borrow<KeyBytes> + Sync>(keys: &[B]) -> Mphf {
    debug_assert_eq!(rayon::current_num_threads(), 1);
    Mphf::try_par_new_with_builder(keys, PHastRBuilder::default(), no_logging![])
        .unwrap_or_else(|e| panic!("MPHF construction failed on {} keys: {e:#}", keys.len()))
}

/// Serialize an MPHF (ε-serde format). Returns the number of bytes written.
pub fn write_mphf(mphf: &Mphf, writer: &mut impl io::Write) -> io::Result<usize> {
    use epserde::ser::Serialize;
    // SAFETY: the function was built by its constructor
    unsafe { mphf.serialize(writer) }.map_err(io::Error::other)
}

/// Exact serialized size of an MPHF in bytes.
pub fn mphf_serialized_bytes(mphf: &Mphf) -> usize {
    write_mphf(mphf, &mut io::sink()).expect("serializing to a sink cannot fail")
}

/// Read (deserialize) an MPHF written by [`write_mphf`].
pub fn read_mphf(mut reader: &mut dyn io::Read) -> io::Result<Mphf> {
    use epserde::deser::Deserialize;
    // SAFETY: the stream contains an `Mphf` serialized by `write_mphf`; the
    // container header versions guard against foreign data
    unsafe { Mphf::deserialize_full(&mut reader) }.map_err(io::Error::other)
}
