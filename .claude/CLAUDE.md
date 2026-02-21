# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

`sshash-rs` is a Rust implementation of **SSHash** — a compressed dictionary for DNA k-mers based on minimal perfect hashing (PHast) and succinct data structures (Elias-Fano, BitFieldVec). It is used as the core k-mer index for `piscem-rs`.

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

- **`ph-temp`** is a COMBINE-lab fork of `bsuccinct-rs` that returns `usize::MAX` for out-of-set keys instead of panicking. The commented-out `ph` git dep in `Cargo.toml` is the upstream; use `ph-temp` until merged.
- **`rapidhash`** (not `ahash`) is used for all index hashing. `ahash` switches algorithms with AES-NI availability, breaking serialized indices across machines. `rapidhash` is CPU-feature independent — indices are portable.
- **`epserde` version must match `sux`'s version** (both currently `0.12`). Mismatched versions cause deserialization failures.
- **`sux = { version = "0.12", features = ["epserde"] }`** — the `epserde` feature is required.

## Serialization Format

Two-file format with magic bytes for versioning:
- `<prefix>.ssi` — main index (SPSS + SparseAndSkewIndex) via epserde
- `<prefix>.ssi.mphf` — MPHF partitions (offset table + raw data) via custom binary format

The header magic (`SSHIDX01`, `SSHIMH02`) is checked on load. Version bumps require updating `serialization.rs`.
