# sshash-lib

Core library for **SSHash-rs**: a compressed k-mer dictionary based on sparse and skew hashing.

## Quick Start

```rust,no_run
use sshash_lib::{Dictionary, Kmer, KmerBits};

type Kmer31 = Kmer<31>;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load a previously built index
    let dict = Dictionary::load("index")?;

    // Single k-mer lookup (returns position or INVALID_UINT64)
    let kmer = Kmer31::from_string("ACGTACGTACGTACGTACGTACGTACGTACG")?;
    let pos = dict.lookup::<31>(&kmer);

    // Streaming queries over a sequence
    let mut engine = dict.create_streaming_query::<31>();

    Ok(())
}
```

## Modules

| Module | Purpose |
|--------|---------|
| `dictionary` | Load, save, lookup, and query k-mers |
| `builder` | Index construction pipeline |
| `streaming_query` | Efficient sequential k-mer processing |
| `kmer` | `Kmer<K>` with const-generic sizing and bit-parallel ops |
| `minimizer` | Minimizer extraction and iteration |
| `minimizers_control_map` | MPHF-based minimizer→bucket mapping |
| `spectrum_preserving_string_set` | SPSS: unitig storage and position lookup |
| `sparse_and_skew_index` | Bucket dispatch (singleton / light / heavy) |
| `offsets` | Elias-Fano encoded string boundary offsets |

## License

MIT
