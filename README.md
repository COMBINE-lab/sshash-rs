# SSHash-rs

Rust implementation of [**SSHash**](https://github.com/jermp/sshash): **S**parse and **S**kew **Hash**ing of k-mers.

A compressed dictionary for k-mers (strings of length k over the DNA alphabet {A, C, G, T}), based on minimal perfect hashing and succinct data structures.

## Important note

This project grew out of my (Rob's) perpetual desire for a Rust implementation of the excellent SSHash data structure. We use SSHash heavily 
for tools across the lab, including `piscem`, which we really want to move completely into Rust.  I've wanted this translation to exist for _years_,
but never had the time for it, and couldn't justify it as a full-time PhD student project.

This initial port was almost entirely created by Claude Opus 4.5 and 4.6, over the span of about 2 days.  It was created in a tight interactive loop 
in which I (Rob) was involved, and the guidance toward a reasonable implementation (i.e. actually using appropriate succinct data structures, the 
selection of libraries, and the direction about tackling critical optimizations) was provided by me.  Nonetheless, the code itself was not written 
by me, nor by any *person*.  The same is true of the remainder of the `README` below (though I will be updating that as I continue to work on this 
project). I feel that this aspect is a critical disclosure about the initial release.  I plan to continue to iterate on this 
library, and this will involve both further agent-driven edits, as well as, likely, considerable "manual" edits to the code to improve the 
idioms, ergonomics and structure.  

During the entire process of development, the C++ implementation (which was _not_ written in an AI-assisted manner, and was almost entirely written 
by Giulio) was treated as the "guide star".  While what Claude was able to accomplish in only a couple of days is, in my opinion, **quite impressive**, 
it certainly would not have been possible without such a clear and correct refrence implementation against which to check (and whose internals could
be freely inspected).

Anyway, I hope this library is useful to others in addition to myself and our lab, and I expect it to be developed and improved going forward.

## Intentionally unimplemented features

This library does not currently support the *weighted* dictionary functionality of SSHash. This was not needed for our purpose, and I wanted to remove non-essential features from the initial version.

## Eventually supported but currently incomplete features

Right now, *almost* all of the index works for all odd values of k from 3 to 63. However, there is a code path in the MPHF for the heavy buckets that doesn't currently support k > 31. Hence, _as of right now_, the index can only be built for k <= 31.

## Building

Requires **Rust 1.85+** (edition 2024).

```bash
cargo build --release
```

The binary is `target/release/sshash`.

## Usage

### Build an index

```bash
# From FASTA/FASTQ (optionally gzipped)
sshash build --input sequences.fa.gz --k 31 --m 20 --output index

# From SPSS format (one unitig per line)
sshash build --input unitigs.spss --k 31 --m 20 --output index

# Canonical mode (reverse-complement normalization)
sshash build --input sequences.fa --k 31 --m 20 --output index --canonical
```

### Query k-mers

```bash
# Single k-mer queries from a file (one k-mer per line, or FASTA)
sshash query --index index --query queries.fa --k 31

# Streaming queries (processes whole sequences, faster for high-hit-rate inputs)
sshash query --index index --query queries.fa --k 31 --streaming
```

### Check correctness

```bash
# Verify that every k-mer in the input can be found in the index
sshash check --index index --input sequences.fa.gz --k 31
```

### Benchmark

```bash
sshash bench --index index
```

## Library usage

Add sshash-lib to your `Cargo.toml`:

```toml
[dependencies]
sshash-lib = { path = "crates/sshash-lib" }
```

```rust,no_run
use sshash_lib::{Dictionary, Kmer, KmerBits};

type Kmer31 = Kmer<31>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dict = Dictionary::load("index")?;

    let kmer = Kmer31::from_string("ACGTACGTACGTACGTACGTACGTACGTACG")?;
    let id = dict.lookup::<31>(&kmer);

    // Streaming queries over a sequence
    let mut engine = dict.create_streaming_query::<31>();
    // ...

    Ok(())
}
```

## Project structure

```
Cargo.toml              # Workspace root
crates/
  sshash-lib/           # Core library
    src/
      dictionary.rs     # Dictionary: load, save, lookup, query
      builder/          # Index construction pipeline
      streaming_query.rs
      kmer.rs           # Kmer<K> with const-generic sizing
      minimizer.rs      # Minimizer extraction
      minimizers_control_map.rs  # MPHF-based minimizer→bucket map
      spectrum_preserving_string_set.rs  # SPSS unitig storage
      sparse_and_skew_index.rs   # Bucket dispatch (singleton/light/heavy)
      offsets.rs         # Elias-Fano string boundary offsets
  sshash-cli/           # CLI binary
    src/main.rs
```

## Key dependencies

| Crate | Purpose |
|-------|---------|
| [`ph`](https://github.com/beling/bsuccinct-rs) | PHast minimal perfect hash functions |
| [`sux`](https://crates.io/crates/sux) | `BitFieldVec` compact bitvectors |
| [`cseq`](https://crates.io/crates/cseq) | Elias-Fano monotone sequences (string offsets) |
| [`needletail`](https://crates.io/crates/needletail) | FASTA/FASTQ parsing |
| [`rayon`](https://crates.io/crates/rayon) | Parallel sorting during build |

## Basic Correctness Checks

Verified against the C++ reference implementation on:

- *Salmonella enterica* (4.8M k-mers): **100%** correctness (regular and canonical)
- Human chromosome 1 (204M k-mers): **100%** correctness

## License

MIT

## References

Giulio Ermanno Pibiri. "[Sparse and Skew Hashing of K-Mers](https://doi.org/10.1093/bioinformatics/btac245)." *Bioinformatics*, 2022.
Giulio Ermanno Pibiri and Rob Patro. "[Optimizing sparse and skew hashing: faster k-mer dictionaries](https://www.biorxiv.org/content/10.64898/2026.01.21.700884v1)." BioRxiv, 2026.
