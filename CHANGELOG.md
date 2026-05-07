# Changelog

## [0.6.0] - 2025-05-07

### Breaking Changes

- **Index format v4.0**: Minimizer hash function changed from rapidhash to a
  multiply-XOR mixer. Indices built with v0.5.0 are incompatible and must be
  rebuilt.
- **sux-rs 0.14**: Updated from sux 0.13.1 to 0.14.0. The public
  `SparseAndSkewIndex.offsets` field type changed from `BitFieldVec<usize>` to
  `BitFieldVec` (Backend-based generics). Downstream code that names this type
  explicitly will need updating. Serialized BitFieldVec data is binary
  compatible — no .ssi format change from this upgrade alone.
- **Orientation widened to i64**: `kmer_orientation` changed from `i8` to `i64`
  in the streaming query API, matching C++ sshash's `int64_t`.

### Performance

- **Custom Elias-Fano with DArray**: Replaced sux-rs `SelectAdaptConst` /
  `EfSeqDict` with a custom Elias-Fano implementation backed by DArray
  (Okanohara & Sadakane 2007) for O(1) Select1/Select0. Memory reduced 7.5x
  (~9.7 MB to ~1.3 MB for gencode v44). Index load time reduced 6x (205 ms to
  33 ms) by serializing DArray directly instead of rebuilding from epserde
  EfSeqDict on load.
- **Streaming query hot path**: Merged validation and k-mer update into a single
  branch; replaced `self.k` with const-generic `K`; restructured `try_extend`
  control flow. Streaming throughput improved from 25.3 to 20.8 ns/kmer on
  gencode v44 (3.14B k-mers).
- **Single-base extension verification**: Consecutive k-mers share K-1 bases, so
  only one new base needs checking against the SPSS text. Buffer caches ~31
  bases per 64-bit read (97% hit rate vs 50% previously). 3.4% throughput
  improvement on generic builds.
- **Faster minimizer hashing**: Multiply-XOR mixer replaces rapidhash for
  minimizer computation, avoiding 128-bit multiply overhead.
- **Base encoding**: Replaced 256-byte lookup table with `(base >> 1) & 3` bit
  trick.

### Dependencies

- `sux`: 0.13.1 → 0.14.0
- `epserde`: 0.12 → 0.12.5
- `mem_dbg`: 0.4.0 → 0.4.1

### Other

- Added `stream-bench` and `point-bench` CLI commands for k-mer query
  benchmarking.
- Enabled thin LTO for release builds.

## [0.5.0] - 2025-04-15

Initial public release.
