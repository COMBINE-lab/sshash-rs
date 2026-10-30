# Changelog

## [0.8.1] - 2026-10-30

### Performance

Query-side only: the index format and every query result are unchanged.

- **In-word select without PDEP** (`DArray`): on aarch64, and on x86-64
  without BMI2, the k-th set bit of a word is now found from per-byte prefix
  popcounts (one multiplication) and a 2 KiB table, instead of six dependent
  popcounts. In isolation, ~43% lower latency on Apple M5 and ~39% on x86-64
  without BMI2; x86-64 with BMI2 keeps PDEP.
- **Elias–Fano adjacent pairs**: `EliasFanoOffsets::access_pair(i)` reads
  elements `i` and `i + 1` with one select (the second is the next one bit),
  used by `SparseAndSkewIndex` bucket bounds; `locate_with_end` walks
  consecutive ones of the high bits instead of doing a select per element.
- **Branch-free minimizer rescan**: `MinimizerIterator`'s re-scan keeps its
  running minimum in locals chosen with `select_unpredictable`, removing a
  frequently mispredicted branch.

In salmon `--sketch` (SRR13638690 / GENCODE v50) these cut mapping time by
~2–3% on Apple M5 and up to ~3% on x86-64, with byte-identical output.

### Dependencies

- sux is used without its default features (`epserde` and `rayon` only):
  sshash never used sux's gzip/zstd line readers, and their `flate2`/`zstd`
  0.13 dependencies are dropped from dependents' trees (zstd 0.13 carries the
  RustSec advisories of 2026-10-16).

## [0.8.0] - 2026-10-19

### Breaking Changes

- **MPHF is now sux's `PHastR`** (sux 0.15.2) instead of `ph`'s PHast
  (`ph-temp`), which is no longer a dependency. Keys are hashed with
  rapidhash v3 through a key-type-erased `KeyBytes` wrapper.
- **Index format v6.0, MPHF container v3.0**: indices built with 0.7.x or
  earlier must be rebuilt. sux 0.15 also moves to epserde 0.13, whose
  serialized format is incompatible with 0.12.
- Dependencies: sux 0.14 → 0.15.2, epserde 0.12 → 0.13, value-traits
  0.2 → 0.3.
- **One setting for how the minimizer MPHF is built**:
  `BuildConfiguration::mphf_build: MphfBuild` (`Streaming(Offline)`, the
  default, `InMemory` or `Partitioned`; `Offline` is `Auto { min_keys }`,
  `Always` or `Never`) replaces the `partitioned_mphf` boolean. Every
  strategy writes the same format and is deterministic in `num_threads`.
  `MphfBuild` and `Offline` are `#[non_exhaustive]`.
- `BuildConfiguration` is `#[non_exhaustive]`: build it with `new(k, m)` or
  `default()` and set fields, not with a struct literal. Its `canonical`
  field (deprecated and ignored since 0.7.0) is removed; `lambda` is
  deprecated and ignored (PHastR has no such parameter).
- `PartitionedMphf`: `build(keys_fn, num_keys, strategy, tmp_dir)` and
  `build_single(&keys)` replace `build_from_vec` and `build_from_slice`.
  Its key types (and those of `SkewIndex::lookup` and `KmerBits::Storage`)
  are bound by `MphfKey` instead of `Hash`.
- `mphf_config`: `Mphf` is now `PHastR<KeyBytes>`. Removed `MphfHasher`,
  `mphf_hasher()`, `mphf_params()` and `build_mphf_from_slice` /
  `build_mphf_from_vec` / `build_mphf_from_slice_mt`. Added `MphfBuild`,
  `Offline`, the `MphfKey` trait, `KeyBytes`, `KeyArray`, `mphf_get`,
  `build_mphf`, `build_mphf_parallel`, `build_mphfs`,
  `build_mphf_streaming(keys, offline, tmp_dir)`, `write_mphf`,
  `mphf_serialized_bytes`, and the `STREAMING_SWEEP_CHUNKS` /
  `OFFLINE_MIN_KEYS` constants.
- Removed `MinimizersControlMapBuilder`, `MinimizerControl`, the control
  map's `BucketType` and `MinimizersControlMap::from_parts` (use `from_mphf`
  or the new `empty`): the minimizer → bucket hash map they built is gone
  (see Performance).
- Builder internals:
  - `MinimizerTuple` is 16 bytes: fields other than `minimizer` are accessor
    methods (`pos_in_seq()`, `pos_in_kmer()`, `num_kmers_in_super_kmer()`),
    with `with_kmers`, `set_num_kmers_in_super_kmer` and
    `MAX_POS_IN_SEQ` (2^52) added; it now derives `Ord` (by minimizer, then
    position, as before) and `Default`.
  - `BucketRef` and `ClassifiedBuckets::bucket_refs` are removed;
    `ClassifiedBuckets` exposes `bucket_size`, `bucket_minimizer`,
    `bucket_type` and `minimizers`.
  - `SparseAndSkewIndex::build_from_classified` takes the per-bucket MPHF
    index (`Vec<u32>`) instead of `Option<&[usize]>` (bucket id by MPHF
    index).
  - `needs_external_sorting` and `estimate_memory_bytes` take `k` and `m`;
    new `estimate_num_tuples`.
  - New `builder::encode::encode_sequences_parallel` and
    `EliasFanoOffsets::from_iter`.
- Absent keys: PHastR maps every key, present or not, to an index in range,
  where `ph`'s PHast could return an out-of-range sentinel. sshash's query
  paths verify the key at the returned position, so lookups are unchanged,
  but code that calls the MPHF directly must not treat an in-range result as
  membership.

### Added

- `sshash build` options `--mphf streaming|in-memory|partitioned`,
  `--mphf-offline auto|always|never` and `--tmp-dir`.

### Changed

- With the default `MphfBuild`, builds of 2^24 minimizers or more (in-memory
  path included) write the MPHF's construction records (16 bytes per
  minimizer, in unlinked temporary files) to `BuildConfiguration::tmp_dirname`,
  which defaults to the relative directory `sshash_tmp`. The directory is
  created if absent and not removed afterwards. Use
  `MphfBuild::Streaming(Offline::Never)`, `InMemory` or `Partitioned` to avoid
  temporary disk.

### Performance

Human (178M minimizers, 32 threads), sshash construction only: 40.9 s →
14.6 s, peak memory 19.1 GB → 7.2 GB; through piscem: 149.6 s / 26.7 GB with
0.7.1 → 78.2 s / 7.4 GB; `.ssi.mphf` 59.5 MB → 54.4 MB. Mapping speed is
unchanged.

- The minimizer MPHF is built as one function by streaming the minimizers
  into PHastR (the default), offline in `tmp_dirname` from 2^24 minimizers on
  (16 bytes/key of disk, about 0.26 GB of memory for 178M minimizers).
- No strategy builds a minimizer → bucket hash map any more: bucket ids are
  their index in the sorted minimizers. With `partitioned_mphf` in 0.7 that
  map and a copy of the minimizers cost ~63 bytes/minimizer (11 GB and 24 s
  for human); `MphfBuild::Partitioned` now takes 1.2 s and peaks at 8.8 GB
  overall.
- Minimizer tuples are 16 bytes instead of 24, collected into one exactly
  reserved vector a wave at a time (no per-thread vectors concatenated).
- Buckets are stored as `u32` starts into the tuples instead of a 40-byte
  descriptor each (7.1 GB → 0.7 GB for human).
- The sparse index is built in MPHF order from a per-bucket `u32` MPHF
  index computed in parallel; the bucket-id-by-MPHF-index array and the
  Elias-Fano input vector are no longer materialized, and offsets are filled
  in parallel.
- Sequences are 2-bit encoded in parallel into an exactly sized buffer
  (13.4 s → 1.2 s for human).
- Partitioned MPHF construction builds at most `num_threads` partitions at a
  time, each in exactly sized memory.

### Fixes

- Indices are byte-identical whatever `num_threads` they are built with:
  the streamed and in-memory minimizer MPHFs use 8 PHastR sweep chunks, and
  the other MPHFs single-thread pools.
- The RAM estimate that chooses external sorting assumed one 48-byte tuple
  per k-mer (67 GB for human, so builds with the default 8 GiB limit always
  went external); it now uses the expected super-k-mer count. Builds that
  could exceed 2^32 tuples always use external sorting.
- `MinimizersControlMap::num_bits` reports the MPHF's actual size instead of
  a 4 bits/key estimate.

### Dependencies

- Picks up the semver-compatible `Cargo.lock` refresh made on `main` after
  0.7.1 (58 crates); no Rust release had included it.

## [0.7.1] - 2026-08-30

### Fixes

- `Cargo.lock` no longer pins the yanked `chacha20 0.10.0` (transitive via
  `rand 0.10` / `sux`); updated to 0.10.2. A full sweep of the lockfile
  found no other yanked versions. Library code is unchanged.

## [0.7.0] - 2026-08-30

### Breaking Changes

- **One indexing modality** (port of C++ SSHash v6.0.0): the minimizer of a
  k-mer x is now the locus minimizing `h(kappa(i))`, where `kappa(i)` is the
  canonical m-mer at locus i, with the centre-closest tie-break (Cologni &
  Pibiri, Proposition 24). The scheme is mirror-equivariant (a k-mer and its
  reverse complement share a bucket: one probe, orientation always reported)
  at plain forward density — the old canonical mode paid a 4/3 super-kmer
  density factor for the same answers. Regular (non-canonical) mode is
  removed: `BuildConfiguration.canonical` and `Dictionary::canonical()` are
  deprecated no-ops, CLI `--canonical` warns. The forward-only query family
  (`query_forward`, `lookup_forward`, `*_checked(…, false)`, CLI
  `--forward-only`, python `forward_only=`) keeps C++ v6 semantics: the
  ordinary single-probe lookup runs, and a match found in backward
  orientation reports not-found — i.e. answers are restricted to k-mers
  occurring in forward orientation in the indexed strings. (On the old
  canonical mode the flag was documented as ignored; it is now honored.)
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
