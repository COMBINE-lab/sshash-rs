# Feature Tracking

## Planned / Under Exploration

### Python Bindings via PyO3

**Status:** Implemented — `crates/sshash-py/` (`cargo check -p sshash-py` passes)

**Goal:** Expose sshash-rs as a Python library so that users can load an sshash index, query individual k-mers, and perform streaming queries over sequences — all from Python — while getting back hit information as convenient Python objects.

**Non-negotiable constraint:** Zero impact on `sshash-lib`. No PyO3 dependency, no `#[cfg(feature = "...")]` guards, no additional compilation overhead for pure-Rust users of the library.

#### Proposed Structure

Add a third crate to the workspace: `crates/sshash-py`. This crate:
- Depends on `sshash-lib` (path dep)
- Depends on `pyo3` (with the `extension-module` feature)
- Is built with `maturin` into a Python wheel (`sshash` package on PyPI)
- Is completely absent from any path that a Rust-only user would take

`sshash-lib` and `sshash-cli` remain untouched.

#### Workspace layout

```
crates/
  sshash-lib/     # unchanged — no new deps, no feature flags
  sshash-cli/     # unchanged
  sshash-py/      # new: PyO3 extension module
    Cargo.toml    # [lib] crate-type = ["cdylib"]
    src/
      lib.rs      # #[pymodule] + #[pyclass] wrappers
    pyproject.toml  # maturin build backend
```

Because `sshash-py` is an optional workspace member, it can either be listed under `[workspace] members` (harmless for Rust-only builds since cdylib crates are never linked into other crates) or kept outside the workspace entirely and pointed at `sshash-lib` via a relative path dep.

#### Python API surface (proposed)

```python
import sshash

# Load a pre-built index (prefix, not full filename)
d = sshash.Dictionary.load("path/to/index")

# Single k-mer lookup — returns int (kmer_id) or None if not found
kmer_id = d.lookup("ACGTACGTACGTACGTACGTACGTACGTACG")

# Full query — returns a Hit namedtuple/dataclass or None
hit = d.query("ACGTACGTACGTACGTACGTACGTACGTACG")
# hit.kmer_id, hit.string_id, hit.orientation, hit.position_in_string

# Streaming query over a sequence (sliding window, faster for sequences)
engine = d.streaming_query()          # creates a StreamingQuery engine
hits = engine.query_sequence(seq_bytes_or_str)
# returns list[Hit | None] (one entry per k-mer position)

# or iterator variant for large sequences:
for hit in engine.iter_sequence(seq):
    ...
```

#### Key implementation decisions to make

1. **How to expose `dispatch_on_k!`**: The Python `Dictionary` object needs to store the runtime `k` and dispatch internally. One option: wrap `Dictionary` in an enum at the Python boundary that holds the already-dispatched `StreamingQuery<K>` per K, created lazily on first `streaming_query()` call.

2. **`LookupResult` / `LocatedHit` as Python types**: Expose as a `#[pyclass]` with `__repr__` and named attributes (not a raw tuple) so they're ergonomic in Python. Consider making them a simple dataclass-like object.

3. **Orientation encoding**: The C++ convention uses a boolean or enum for fw/rc. Expose as a Python `str` `"fw"` / `"rc"` or a Python `bool` — decide based on ergonomics.

4. **Bytes vs str input for sequences**: Accept both `str` and `bytes` for k-mer/sequence arguments (PyO3's `#[pyo3(from_py_with = ...)]` or accept `&PyAny` and downcast).

5. **Error handling**: Map `anyhow::Error` → `PyErr` (likely `PyIOError` for load failures, `PyValueError` for bad k-mer strings).

6. **Maturin vs setuptools-rust**: Maturin is the modern standard and handles `pyproject.toml` + sdist/wheel building cleanly. Prefer maturin.

7. **Python packaging**: Decide on package name (`sshash`, `sshash-rs`, `pysshash`?) and whether to publish to PyPI.

#### Implementation notes

- `lookup_with_dict` on `StreamingQuery<K>` was made `pub` (was `pub(crate)`) — this lets `ConcreteEngine` in `sshash-py` call it without unsafe, by holding an `Arc<Dictionary>` and passing `dict.as_ref()` at each call.
- `#[pyclass(unsendable)]` on `PyStreamingQuery` and `SequenceIterator` — these hold mutable state; `unsendable` opts out of the `Sync` requirement while PyO3's GIL ensures safe single-thread access.
- `__next__` in `SequenceIterator` materialises the k-mer slice into an owned `Vec<u8>` before calling `engine.lookup()`, avoiding a split-borrow through `PyRefMut`.
- `dispatch_on_k!` runs exactly once at `Dictionary::load()` — all subsequent calls are vtable dispatches through `Box<dyn DynDict>`.

#### Open questions / future work

- NumPy array output for `kmer_id` bulk queries (v2)
- PyPI package name decision (`sshash` vs `sshash-rs` vs `pysshash`)
