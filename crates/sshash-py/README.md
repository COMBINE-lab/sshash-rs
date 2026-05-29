# sshash

Python bindings for [**sshash-rs**](https://github.com/COMBINE-lab/sshash-rs) — a compressed dictionary for DNA k-mers based on **S**parse and **S**kew **Hash**ing.

sshash stores a set of k-mers (strings of length *k* over {A, C, G, T}) compactly using minimal perfect hashing and succinct data structures (Elias-Fano, BitFieldVec), and supports fast individual and streaming lookups. It is the k-mer index underlying the [piscem](https://github.com/COMBINE-lab/piscem-rs) read mapper.

## Installation

```bash
pip install sshash
```

## Building an index

### From a FASTA/FASTQ file

```python
import sshash

config = sshash.BuildConfig(k=31, m=19)
config.canonical = True   # k-mer and its reverse complement map to the same entry
config.threads = 8        # parallel build (0 = all cores)
config.verbose = False

dict = config.build_from_file("sequences.fa.gz")
dict.save("my_index")
```

### From a list of sequences in memory

```python
config = sshash.BuildConfig(k=31, m=19)
sequences = ["ACGTACGTACGTACGTACGTACGTACGTACG",
             "TTGCAACCGTTAGCAACGTACGTACGTACGT"]
dict = config.build(sequences)
```

### From a Cuttlefish `.cf_seg` file

When sequences come from [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish), `build_from_cf_seg` also returns a mapping from sshash string IDs back to the original Cuttlefish node IDs:

```python
config = sshash.BuildConfig(k=31, m=19)
dict, segment_ids = config.build_from_cf_seg("unitigs.cf_seg")
# segment_ids[i] is the Cuttlefish node ID for sshash string_id i
```

## Loading and saving

```python
# Save to disk (writes <prefix>.ssi and <prefix>.ssi.mphf)
dict.save("my_index")

# Load from disk
dict = sshash.Dictionary.load("my_index")
```

## Querying

### Single k-mer lookup

```python
# Returns a Hit object, or None if not found
hit = dict.query("ACGTACGTACGTACGTACGTACGTACGTACG")
if hit is not None:
    print(hit.kmer_id)           # global k-mer ID
    print(hit.string_id)         # unitig containing this k-mer
    print(hit.kmer_id_in_string) # position within that unitig
    print(hit.orientation)       # +1 forward, -1 reverse complement

# Just the k-mer ID (faster if location info isn't needed)
kmer_id = dict.lookup("ACGTACGTACGTACGTACGTACGTACGTACG")  # None if absent

# Membership test
present = dict.contains("ACGTACGTACGTACGTACGTACGTACGTACG")
```

By default, a query on a non-canonical (strand-specific) index falls back to the
reverse complement, so the reverse complement of an indexed k-mer is still found
(reported with `orientation == -1`). Pass `forward_only=True` for a strand-specific
query that does **not** perform this fallback:

```python
# Non-canonical index: only match the k-mer in its given orientation
hit = dict.query("ACGTACGTACGTACGTACGTACGTACGTACG", forward_only=True)
kmer_id = dict.lookup("ACGTACGTACGTACGTACGTACGTACGTACG", forward_only=True)
present = dict.contains("ACGTACGTACGTACGTACGTACGTACGTACG", forward_only=True)
```

`forward_only` is ignored for canonical indexes, where both strands are inherently
equivalent.

### Streaming queries over a sequence

The streaming engine maintains minimizer state across consecutive k-mers, avoiding redundant MPHF lookups for adjacent positions. This is significantly faster than calling `query` in a loop when processing full reads or contigs.

```python
engine = dict.streaming_query()

# Query all k-mers in a sequence at once (returns list)
hits = engine.query_sequence("ACGTACGTACGTACGTACGTACGTACGTACGTACGT")
for hit in hits:
    if hit is not None:
        print(hit.kmer_id, hit.string_id)

# Lazy iterator (memory-efficient for long sequences)
for hit in engine.iter_sequence(b"ACGTACGT..."):
    if hit is not None:
        print(hit.kmer_id)

# Efficiency statistics
print(engine.num_searches)    # full MPHF lookups performed
print(engine.num_extensions)  # k-mers resolved by sliding-window extension
```

## Index properties

```python
print(dict.k)           # k-mer length
print(dict.m)           # minimizer length
print(dict.canonical)   # canonical mode flag
print(dict.num_strings) # number of unitigs
print(dict.num_bits)    # total index size in bits
```

## BuildConfig options

| Property | Default | Description |
|---|---|---|
| `canonical` | `False` | Map each k-mer and its reverse complement to the same entry |
| `threads` | `0` | Worker threads during build (`0` = all available cores) |
| `ram_limit_gib` | `8` | RAM budget (GiB) before switching to external sort |
| `seed` | internal | Seed for internal hash functions |
| `verbose` | `True` | Print progress during building |
| `tmp_dir` | `"sshash_tmp"` | Directory for temporary files during external sort |

`k` and `m` are set at construction time and cannot be changed afterwards.

## References

- Giulio Ermanno Pibiri. "[Sparse and Skew Hashing of K-Mers](https://doi.org/10.1093/bioinformatics/btac245)." *Bioinformatics*, 2022.
- Giulio Ermanno Pibiri and Rob Patro. "[Optimizing sparse and skew hashing: faster k-mer dictionaries](https://www.biorxiv.org/content/10.64898/2026.01.21.700884v1)." *bioRxiv*, 2026.

## License

BSD 3-Clause
