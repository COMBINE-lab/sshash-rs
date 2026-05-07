//! Python bindings for the sshash k-mer dictionary.
//!
//! Exposes three Python classes:
//! - `Dictionary` — load an sshash index, query individual k-mers
//! - `StreamingQuery` — sliding-window queries over DNA sequences
//! - `Hit` — result of a successful k-mer lookup

use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use sshash_lib::builder::{parse_cf_seg, BuildConfiguration, DictionaryBuilder};
use sshash_lib::builder::parse::parse_sequences;
use sshash_lib::{dispatch_on_k, Dictionary, Kmer, KmerBits, LookupResult, StreamingQuery};
use std::sync::Arc;

// ─── Hit ─────────────────────────────────────────────────────────────────────

/// Result of a successful k-mer lookup.
#[pyclass(get_all, frozen, skip_from_py_object)]
#[derive(Clone)]
struct Hit {
    /// Absolute k-mer ID across the entire index
    kmer_id: u64,
    /// 0-based position of this k-mer within its containing string (unitig)
    kmer_id_in_string: u64,
    /// ID of the string (unitig) that contains this k-mer
    string_id: u64,
    /// Orientation: +1 for forward strand, -1 for reverse complement
    orientation: i64,
    /// Start position (in bases) of the containing string
    string_begin: u64,
    /// End position (in bases) of the containing string
    string_end: u64,
}

#[pymethods]
impl Hit {
    fn __repr__(&self) -> String {
        format!(
            "Hit(kmer_id={}, kmer_id_in_string={}, string_id={}, orientation={}, \
             string_begin={}, string_end={})",
            self.kmer_id,
            self.kmer_id_in_string,
            self.string_id,
            self.orientation,
            self.string_begin,
            self.string_end,
        )
    }
}

fn result_to_hit(r: LookupResult) -> Option<Hit> {
    if r.is_found() {
        Some(Hit {
            kmer_id: r.kmer_id,
            kmer_id_in_string: r.kmer_id_in_string,
            string_id: r.string_id,
            orientation: r.kmer_orientation,
            string_begin: r.string_begin,
            string_end: r.string_end,
        })
    } else {
        None
    }
}

// ─── DynEngine ───────────────────────────────────────────────────────────────
//
// Streaming engine trait erased over K. ConcreteEngine<K> holds a StreamingQuery<K>
// (no borrow from Dictionary) and an Arc<Dictionary> (passed by reference at each
// lookup call). No unsafe required.

trait DynEngine: Send {
    fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult;
    fn reset(&mut self);
    fn num_searches(&self) -> u64;
    fn num_extensions(&self) -> u64;
    fn k(&self) -> usize;
    /// Create a fresh engine of the same K backed by the same Arc<Dictionary>.
    fn make_sibling(&self) -> Box<dyn DynEngine + Send>;
}

struct ConcreteEngine<const K: usize>
where
    Kmer<K>: KmerBits,
{
    query: StreamingQuery<K>,
    dict: Arc<Dictionary>,
}

impl<const K: usize> ConcreteEngine<K>
where
    Kmer<K>: KmerBits,
{
    fn new(dict: Arc<Dictionary>) -> Self {
        let query = StreamingQuery::new(dict.k(), dict.m(), dict.canonical());
        ConcreteEngine { query, dict }
    }
}

impl<const K: usize> DynEngine for ConcreteEngine<K>
where
    Kmer<K>: KmerBits,
{
    fn lookup(&mut self, kmer_bytes: &[u8]) -> LookupResult {
        self.query.lookup_with_dict(kmer_bytes, self.dict.as_ref())
    }

    fn reset(&mut self) {
        self.query.reset();
    }

    fn num_searches(&self) -> u64 {
        self.query.num_searches()
    }

    fn num_extensions(&self) -> u64 {
        self.query.num_extensions()
    }

    fn k(&self) -> usize {
        K
    }

    fn make_sibling(&self) -> Box<dyn DynEngine + Send> {
        Box::new(ConcreteEngine::<K>::new(Arc::clone(&self.dict)))
    }
}

// ─── DynDict ─────────────────────────────────────────────────────────────────
//
// Dictionary trait erased over K. Dispatched exactly once at load time via
// dispatch_on_k!. All subsequent queries are plain vtable calls.

trait DynDict: Send + Sync {
    fn lookup_str(&self, kmer: &str) -> PyResult<Option<u64>>;
    fn query_str(&self, kmer: &str) -> PyResult<Option<Hit>>;
    fn contains_str(&self, kmer: &str) -> PyResult<bool>;
    fn make_engine(&self) -> Box<dyn DynEngine + Send>;
    fn k(&self) -> usize;
    fn m(&self) -> usize;
    fn canonical(&self) -> bool;
    fn num_strings(&self) -> u64;
    fn num_bits(&self) -> u64;
    fn save(&self, path: &str) -> PyResult<()>;
}

struct ConcreteDynDict<const K: usize>
where
    Kmer<K>: KmerBits,
{
    inner: Arc<Dictionary>,
}

impl<const K: usize> DynDict for ConcreteDynDict<K>
where
    Kmer<K>: KmerBits,
{
    fn lookup_str(&self, kmer: &str) -> PyResult<Option<u64>> {
        let k = Kmer::<K>::from_string(kmer)
            .map_err(|e| PyValueError::new_err(format!("Invalid k-mer '{}': {}", kmer, e)))?;
        let id = self.inner.lookup::<K>(&k);
        Ok(if id == u64::MAX { None } else { Some(id) })
    }

    fn query_str(&self, kmer: &str) -> PyResult<Option<Hit>> {
        let k = Kmer::<K>::from_string(kmer)
            .map_err(|e| PyValueError::new_err(format!("Invalid k-mer '{}': {}", kmer, e)))?;
        Ok(result_to_hit(self.inner.query::<K>(&k)))
    }

    fn contains_str(&self, kmer: &str) -> PyResult<bool> {
        Ok(self.lookup_str(kmer)?.is_some())
    }

    fn make_engine(&self) -> Box<dyn DynEngine + Send> {
        Box::new(ConcreteEngine::<K>::new(Arc::clone(&self.inner)))
    }

    fn k(&self) -> usize { self.inner.k() }
    fn m(&self) -> usize { self.inner.m() }
    fn canonical(&self) -> bool { self.inner.canonical() }
    fn num_strings(&self) -> u64 { self.inner.num_strings() }
    fn num_bits(&self) -> u64 { self.inner.num_bits() }

    fn save(&self, path: &str) -> PyResult<()> {
        self.inner
            .save(path)
            .map_err(|e| PyIOError::new_err(format!("Failed to save index: {}", e)))
    }
}

fn make_dyn_dict(dict: Dictionary) -> Box<dyn DynDict + Send + Sync> {
    let k = dict.k();
    let arc = Arc::new(dict);
    // Single dispatch on k — all subsequent queries bypass this match.
    dispatch_on_k!(k, K => Box::new(ConcreteDynDict::<K> { inner: arc }))
}

// ─── PyDictionary ─────────────────────────────────────────────────────────────

/// Compressed k-mer dictionary loaded from an sshash index.
///
/// Load with ``Dictionary.load(prefix)`` where ``prefix`` is the path used
/// when building the index (without the ``.ssi`` / ``.ssi.mphf`` extensions).
#[pyclass(name = "Dictionary")]
struct PyDictionary {
    // K is baked into the vtable at load time — no per-query dispatch needed.
    inner: Box<dyn DynDict + Send + Sync>,
}

#[pymethods]
impl PyDictionary {
    /// Load an sshash index from disk.
    ///
    /// :param path: Index prefix (same value passed to ``sshash build --output``).
    /// :raises IOError: If the index files cannot be read.
    #[staticmethod]
    fn load(path: &str) -> PyResult<Self> {
        let dict = Dictionary::load(path)
            .map_err(|e| PyIOError::new_err(format!("Failed to load index '{}': {}", path, e)))?;
        Ok(PyDictionary { inner: make_dyn_dict(dict) })
    }

    /// Save the index to disk.
    ///
    /// :param path: Output prefix.
    fn save(&self, path: &str) -> PyResult<()> {
        self.inner.save(path)
    }

    /// K-mer size used by this index.
    #[getter]
    fn k(&self) -> usize { self.inner.k() }

    /// Minimizer size used by this index.
    #[getter]
    fn m(&self) -> usize { self.inner.m() }

    /// Whether this index was built in canonical mode.
    #[getter]
    fn canonical(&self) -> bool { self.inner.canonical() }

    /// Number of strings (unitigs) in the index.
    #[getter]
    fn num_strings(&self) -> u64 { self.inner.num_strings() }

    /// Total index size in bits.
    #[getter]
    fn num_bits(&self) -> u64 { self.inner.num_bits() }

    /// Return the global k-mer ID, or ``None`` if the k-mer is not in the index.
    ///
    /// :param kmer: DNA string of length ``k``.
    /// :raises ValueError: If ``kmer`` contains invalid characters or has the wrong length.
    fn lookup(&self, kmer: &str) -> PyResult<Option<u64>> {
        self.inner.lookup_str(kmer)
    }

    /// Return a :class:`Hit` with full location information, or ``None`` if not found.
    ///
    /// :param kmer: DNA string of length ``k``.
    /// :raises ValueError: If ``kmer`` contains invalid characters or has the wrong length.
    fn query(&self, kmer: &str) -> PyResult<Option<Hit>> {
        self.inner.query_str(kmer)
    }

    /// Return ``True`` if the k-mer is present in the index.
    ///
    /// :param kmer: DNA string of length ``k``.
    /// :raises ValueError: If ``kmer`` contains invalid characters or has the wrong length.
    fn contains(&self, kmer: &str) -> PyResult<bool> {
        self.inner.contains_str(kmer)
    }

    /// Create a :class:`StreamingQuery` engine for efficient sliding-window queries.
    ///
    /// The engine maintains minimizer state across consecutive k-mers, avoiding
    /// redundant MPHF lookups when querying adjacent positions in a sequence.
    fn streaming_query(&self) -> PyStreamingQuery {
        // No dispatch_on_k! here — ConcreteDynDict::<K>::make_engine() does it directly.
        PyStreamingQuery { engine: self.inner.make_engine() }
    }

    fn __repr__(&self) -> String {
        format!(
            "Dictionary(k={}, m={}, canonical={}, num_strings={})",
            self.inner.k(),
            self.inner.m(),
            self.inner.canonical(),
            self.inner.num_strings(),
        )
    }
}

// ─── PyStreamingQuery ─────────────────────────────────────────────────────────

/// Streaming query engine for efficient consecutive k-mer lookups.
///
/// Maintains minimizer state across calls so that adjacent k-mers in a sequence
/// can often be resolved by extension rather than a full MPHF lookup.
///
/// Obtain via :meth:`Dictionary.streaming_query`.
// `unsendable`: engine holds mutable state; Python's GIL ensures single-thread access.
#[pyclass(name = "StreamingQuery", unsendable)]
struct PyStreamingQuery {
    engine: Box<dyn DynEngine + Send>,
}

#[pymethods]
impl PyStreamingQuery {
    /// Reset the engine state. Call this before processing a new, non-consecutive sequence.
    fn reset(&mut self) {
        self.engine.reset();
    }

    /// Query all k-mers in a sequence, returning a list.
    ///
    /// Resets the engine state before processing. Each entry in the returned list
    /// corresponds to one k-mer position (sliding by 1 base). Positions where the
    /// k-mer is absent return ``None``.
    ///
    /// :param seq: DNA sequence as ``str`` or ``bytes``.
    /// :returns: ``list[Hit | None]`` of length ``max(0, len(seq) - k + 1)``.
    fn query_sequence(&mut self, seq: &[u8]) -> Vec<Option<Hit>> {
        let k = self.engine.k();
        if seq.len() < k {
            return Vec::new();
        }
        self.engine.reset();
        let n = seq.len() - k + 1;
        let mut results = Vec::with_capacity(n);
        for i in 0..n {
            results.push(result_to_hit(self.engine.lookup(&seq[i..i + k])));
        }
        results
    }

    /// Query all k-mers in a sequence, returning a lazy iterator.
    ///
    /// Unlike :meth:`query_sequence`, this does not materialise all results at once,
    /// making it more memory-efficient for very long sequences. Each iteration step
    /// yields ``Hit | None`` for one k-mer position.
    ///
    /// A fresh engine state is used for each call (independent of ``self``'s state).
    ///
    /// :param seq: DNA sequence as ``str`` or ``bytes``.
    /// :returns: Iterator yielding ``Hit | None``.
    fn iter_sequence(&self, py: Python<'_>, seq: Vec<u8>) -> PyResult<Py<SequenceIterator>> {
        let k = self.engine.k();
        // Fresh engine so this iterator's state is independent of self.
        let engine = self.engine.make_sibling();
        Py::new(py, SequenceIterator { engine, seq, pos: 0, k })
    }

    /// Number of full MPHF lookups performed since the last reset.
    #[getter]
    fn num_searches(&self) -> u64 {
        self.engine.num_searches()
    }

    /// Number of k-mers resolved by extension (no MPHF lookup needed).
    #[getter]
    fn num_extensions(&self) -> u64 {
        self.engine.num_extensions()
    }
}

// ─── SequenceIterator ─────────────────────────────────────────────────────────

/// Lazy iterator over k-mer hits in a sequence. Obtained from
/// :meth:`StreamingQuery.iter_sequence`.
// `unsendable`: holds mutable engine state; GIL-bound single-thread access.
#[pyclass(unsendable)]
struct SequenceIterator {
    engine: Box<dyn DynEngine + Send>,
    seq: Vec<u8>,
    pos: usize,
    k: usize,
}

#[pymethods]
impl SequenceIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<Option<Hit>> {
        if slf.pos + slf.k > slf.seq.len() {
            return None; // StopIteration
        }
        // Materialise the k-mer bytes before the mutable borrow of slf.engine.
        // The block ends the immutable borrow of slf.seq before lookup() borrows
        // slf.engine mutably — the two borrows go through the same PyRefMut.
        let kmer_bytes: Vec<u8> = {
            let s: &SequenceIterator = &*slf;
            s.seq[s.pos..s.pos + s.k].to_owned()
        };
        let result = slf.engine.lookup(&kmer_bytes);
        slf.pos += 1;
        Some(result_to_hit(result))
    }
}

// ─── PyBuildConfig ────────────────────────────────────────────────────────────

/// Configuration for building an sshash index.
///
/// Create with ``BuildConfig(k, m)``, adjust optional fields, then call one of
/// the build methods to get a :class:`Dictionary`.
///
/// Example::
///
///     config = sshash.BuildConfig(k=31, m=19)
///     config.canonical = True
///     config.threads = 8
///     dict = config.build(["ACGTACGT...", "TTGCAACCG..."])
#[pyclass(name = "BuildConfig")]
struct PyBuildConfig {
    config: BuildConfiguration,
}

#[pymethods]
impl PyBuildConfig {
    /// Create a new build configuration.
    ///
    /// :param k: K-mer length (must be odd, 3 ≤ k ≤ 63).
    /// :param m: Minimizer length (must be odd, m < k).
    /// :raises ValueError: If k or m are out of range.
    #[new]
    fn new(k: usize, m: usize) -> PyResult<Self> {
        let config = BuildConfiguration::new(k, m)
            .map_err(PyValueError::new_err)?;
        Ok(PyBuildConfig { config })
    }

    /// K-mer length (set at construction time).
    #[getter]
    fn k(&self) -> usize { self.config.k }

    /// Minimizer length (set at construction time).
    #[getter]
    fn m(&self) -> usize { self.config.m }

    /// Whether to build in canonical mode (k-mer and its reverse complement map to the same entry).
    #[getter]
    fn canonical(&self) -> bool { self.config.canonical }

    #[setter(canonical)]
    fn set_canonical(&mut self, val: bool) { self.config.canonical = val; }

    /// Number of threads (0 = all available cores).
    #[getter]
    fn threads(&self) -> usize { self.config.num_threads }

    #[setter(threads)]
    fn set_threads(&mut self, val: usize) { self.config.num_threads = val; }

    /// RAM limit in GiB for external sorting of minimizer tuples.
    #[getter]
    fn ram_limit_gib(&self) -> usize { self.config.ram_limit_gib }

    #[setter(ram_limit_gib)]
    fn set_ram_limit_gib(&mut self, val: usize) { self.config.ram_limit_gib = val; }

    /// Seed for internal hash functions.
    #[getter]
    fn seed(&self) -> u64 { self.config.seed }

    #[setter(seed)]
    fn set_seed(&mut self, val: u64) { self.config.seed = val; }

    /// Whether to print verbose progress messages during building.
    #[getter]
    fn verbose(&self) -> bool { self.config.verbose }

    #[setter(verbose)]
    fn set_verbose(&mut self, val: bool) { self.config.verbose = val; }

    /// Temporary directory used during external sorting.
    #[getter]
    fn tmp_dir(&self) -> String {
        self.config.tmp_dirname.to_string_lossy().into_owned()
    }

    #[setter(tmp_dir)]
    fn set_tmp_dir(&mut self, val: &str) {
        self.config.tmp_dirname = val.into();
    }

    /// Build an index from a list of in-memory DNA sequences.
    ///
    /// The GIL is released for the duration of the build.
    ///
    /// :param sequences: List of DNA strings (only A/C/G/T characters).
    /// :returns: A ready-to-query :class:`Dictionary`.
    /// :raises RuntimeError: If the build fails.
    fn build(&self, py: Python<'_>, sequences: Vec<String>) -> PyResult<PyDictionary> {
        let config = self.config.clone();
        let dict = py.detach(|| {
            DictionaryBuilder::new(config)
                .and_then(|b| b.build_from_sequences(sequences))
        }).map_err(PyRuntimeError::new_err)?;
        Ok(PyDictionary { inner: make_dyn_dict(dict) })
    }

    /// Build an index from a FASTA or FASTQ file (plain or gzip-compressed).
    ///
    /// The GIL is released for the duration of parsing and building.
    ///
    /// :param path: Path to the input file.
    /// :returns: A ready-to-query :class:`Dictionary`.
    /// :raises RuntimeError: If the file cannot be read or the build fails.
    fn build_from_file(&self, py: Python<'_>, path: &str) -> PyResult<PyDictionary> {
        let config = self.config.clone();
        let path = path.to_string();
        let dict = py.detach(|| -> Result<Dictionary, String> {
            let mut sequences: Vec<String> = Vec::new();
            parse_sequences(&path, |_name, seq| {
                sequences.push(String::from_utf8_lossy(seq).into_owned());
                Ok(())
            }).map_err(|e| e.to_string())?;
            DictionaryBuilder::new(config)?.build_from_sequences(sequences)
        }).map_err(PyRuntimeError::new_err)?;
        Ok(PyDictionary { inner: make_dyn_dict(dict) })
    }

    /// Build an index from a Cuttlefish ``.cf_seg`` segment file.
    ///
    /// Returns both the :class:`Dictionary` and a list of segment IDs. The
    /// ``i``-th element of the segment ID list is the original Cuttlefish node
    /// ID for sshash string ID ``i``.
    ///
    /// The GIL is released for the duration of parsing and building.
    ///
    /// :param path: Path to the ``.cf_seg`` file.
    /// :returns: ``(Dictionary, list[int])`` — index and segment ID mapping.
    /// :raises RuntimeError: If the file cannot be read or the build fails.
    fn build_from_cf_seg(
        &self,
        py: Python<'_>,
        path: &str,
    ) -> PyResult<(Py<PyDictionary>, Vec<u64>)> {
        let config = self.config.clone();
        let path = path.to_string();
        let (dict, segment_ids) = py.detach(|| -> Result<(Dictionary, Vec<u64>), String> {
            let cf_data = parse_cf_seg(&path).map_err(|e| e.to_string())?;
            let segment_ids = cf_data.segment_ids;
            let sequences = cf_data.sequences;
            let dict = DictionaryBuilder::new(config)?.build_from_sequences(sequences)?;
            Ok((dict, segment_ids))
        }).map_err(PyRuntimeError::new_err)?;
        let py_dict = Py::new(py, PyDictionary { inner: make_dyn_dict(dict) })?;
        Ok((py_dict, segment_ids))
    }

    fn __repr__(&self) -> String {
        format!(
            "BuildConfig(k={}, m={}, canonical={}, threads={}, ram_limit_gib={})",
            self.config.k,
            self.config.m,
            self.config.canonical,
            self.config.num_threads,
            self.config.ram_limit_gib,
        )
    }
}

// ─── Module ──────────────────────────────────────────────────────────────────

/// Python bindings for the sshash compressed k-mer dictionary.
#[pymodule]
fn sshash(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyDictionary>()?;
    m.add_class::<PyStreamingQuery>()?;
    m.add_class::<Hit>()?;
    m.add_class::<SequenceIterator>()?;
    m.add_class::<PyBuildConfig>()?;
    Ok(())
}
