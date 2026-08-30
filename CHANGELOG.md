# Changelog

## [0.7.0] - unreleased

### Breaking Changes

- **One indexing modality** (port of C++ SSHash v6.0.0): the minimizer of a
  k-mer x is now the locus minimizing `h(kappa(i))`, where `kappa(i)` is the
  canonical m-mer at locus i, with the centre-closest tie-break (Cologni &
  Pibiri, Proposition 24). The scheme is mirror-equivariant (a k-mer and its
  reverse complement share a bucket: one probe, orientation always reported)
  at plain forward density — the old canonical mode paid a 4/3 super-kmer
  density factor for the same answers. Regular (non-canonical) mode is
  removed: `BuildConfiguration.canonical` and `Dictionary::canonical()` are
  deprecated no-ops, CLI `--canonical` warns.
- **Index format v5.0**: indices built with 0.6.x must be rebuilt. The header
  now records the build seed and a hasher-magic guard (a rapidhash behavior
  change now fails loudly at load instead of silently corrupting minimizer
  selection).
- `Dictionary::new` takes the build seed instead of a canonical flag;
  `SparseAndSkewIndex`/`SkewIndex` builders lose their `canonical` parameter.

### Performance

- Builder: canonical extraction no longer runs a fresh full-rescan RC
  minimizer iterator per k-mer; one iterator serves both strands.
- Point lookup: one minimizer scan and one bucket probe (was two scans and up
  to two probes), with exactly two candidate text positions per occurrence.
- Streaming: same-minimizer memos (C++ c22c897) — minimizer-absent negatives,
  singleton same-occurrence skips, and a cached decoded locate set verified
  without MPHF or offset work. New counters in `StreamingQueryStats`.
- Expected from the C++ measurements: ~22% fewer super-kmers and ~0.5-0.7
  bits/kmer smaller than the old canonical mode, faster builds, and 9-11%
  faster streaming on error-containing reads.

### Fixes

- Wide-K (K >= 33) streaming extension used u64 operations on u128 k-mer
  storage (`try_extend` top-base extract, `append_base` mask): debug panic /
  wrong results in release. Analogue of C++ 87cab5f.
- Queries now hash with the seed the index was built with; previously a
  non-default `BuildConfiguration.seed` produced an index that returned zero
  hits (queries hard-coded seed 1).
- The build now aborts with "the minimizer scheme is not forward" if a
  duplicate (minimizer, position) pair would mis-size the index.

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
