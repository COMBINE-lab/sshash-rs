# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

`sshash-rs` is a Rust implementation of **SSHash** — a compressed dictionary for DNA k-mers based on minimal perfect hashing (sux PHastR) and succinct data structures (Elias-Fano, BitFieldVec). It is used as the core k-mer index for `piscem-rs`.

## Commands

```bash
cargo build --release          # Build the sshash CLI binary at target/release/sshash
cargo test --lib               # Run all unit tests (110+)
cargo test                     # Run unit + integration tests
cargo check                    # Fast compile check, no warnings
cargo clippy                   # Lint
```

Run a single test:
```bash
cargo test --lib test_name
cargo test --test build_pipeline_integration  # integration test
```

## Architecture

The workspace has two crates: `sshash-lib` (core library) and `sshash-cli` (binary).

**Build pipeline** (`builder/dictionary_builder.rs`):
1. Parse FASTA/FASTQ → encode to 2-bit SPSS (Spectrum-Preserving String Set) with Elias-Fano offsets
2. Extract minimizer tuples (rolling hash, optional external sort if RAM limit exceeded)
3. Classify buckets: Singleton / Light (linear scan) / Heavy (per-bucket MPHF)
4. Build MinimizersControlMap (MCM): MPHF mapping minimizer → bucket ID
5. Build SparseAndSkewIndex: Elias-Fano bucket boundaries + BitFieldVec per-k-mer offsets + SkewIndex for heavy buckets
6. Serialize to 2-file format: `.ssi` (main index, epserde) + `.ssi.mphf` (MPHF container, custom binary)

**Query path** (`dictionary.rs`, `streaming_query.rs`):
- Extract minimizer → MCM lookup (MPHF) → bucket range (Elias-Fano Succ) → size-based routing → retrieve position from BitFieldVec
- `StreamingQuery<K>`: sliding-window optimization for consecutive k-mers in a sequence
- `LookupResult`: returns `(kmer_id, string_id, orientation, position_in_string)`

**K-mer representation** (`kmer.rs`): `Kmer<K>` is const-generic; K ≤ 31 uses `u64`, K ∈ [33, 63] uses `u128`. Only odd K values in [3, 63] are supported.

## Key Patterns

**Runtime-to-const-generic dispatch** — always use this macro when K is only known at runtime:
```rust
use sshash_lib::dispatch_on_k;
dispatch_on_k!(k, K => {
    let kmer = Kmer::<K>::from_string(s)?;
    dict.lookup::<K>(&kmer)
})
```

**Trait imports required for sux-rs types:**
```rust
use value_traits::slices::{SliceByValue, SliceByValueMut};  // BitFieldVec index_value/set_value
use epserde::ser::Serialize;
use epserde::deser::Deserialize;
use sux::traits::{IndexedSeq, Succ};  // Elias-Fano
use mem_dbg::{MemSize, SizeFlags};  // ef.mem_size(SizeFlags::default())
```

## Critical Dependency Notes

- **MPHF** is sux's `PHastR` over `KeyBytes` (see `mphf_config.rs`), hashed with rapidhash v3, not sux's default xxh3 (5-15% slower queries). PHastR's output depends on its number of sweep chunks, which defaults to the rayon pool size: always build with a fixed `sweep_chunks` (or on a 1-thread pool) so indices don't depend on `num_threads`. Offline construction must get its directory via `PHastRBuilder::temp_dir`, never `TMPDIR`.
- **MPHF strategy** is `BuildConfiguration::mphf_build` (`MphfBuild::{Streaming(Offline), InMemory, Partitioned}`), dispatched once through `PartitionedMphf::build(keys_fn, n, strategy, tmp_dir)` on both build paths. Every strategy must stay deterministic in `num_threads` (`tests/build_paths_and_determinism.rs`); `Partitioned` builds the same partitioned function the removed hash-map path built (verified byte-identical at human scale during 0.8 development).
- **Build memory** (human, 178M minimizers: 7.2 GB peak): `MinimizerTuple` is 16 bytes (pos/pos_in_kmer/num_kmers packed in one `u64`); `ClassifiedBuckets` stores only `u32` bucket starts; step 5 works in MPHF order from a `u32` per-bucket MPHF index. The in-memory path therefore supports < 2^32 tuples and falls back to external sort beyond. Don't reintroduce per-bucket structs, `Vec` growth by `push` for large arrays, or `collect()` of per-thread `Vec`s (2x transients).
- **`rapidhash`** (not `ahash`) is used for all index hashing. `ahash` switches algorithms with AES-NI availability, breaking serialized indices across machines. `rapidhash` is CPU-feature independent — indices are portable.
- **`epserde` version must match `sux`'s version** (currently `0.13` for sux `0.15`). Mismatched versions cause deserialization failures; an epserde major bump breaks every serialized index.
- **`sux = { version = "0.15.2", features = ["epserde"] }`** — the `epserde` feature is required.

## Serialization Format

Two-file format with magic bytes for versioning:
- `<prefix>.ssi` — main index (SPSS + SparseAndSkewIndex) via epserde
- `<prefix>.ssi.mphf` — MPHF partitions (offset table + raw data) via custom binary format

The header magic (`SSHIDX01`, `SSHIMH02`) is checked on load. Version bumps require updating `serialization.rs`.
