# Experimental batching and minimizer fingerprints

The `experimental-query` library feature adds two independent query experiments:
fixed-size batches of independent seeds and an optional minimizer-absence filter.
Neither changes index construction, IDs, orientation, or saved index bytes. The
feature is off by default; enabling the feature alone does not allocate a filter.
The existing streaming engine retains its rolling minimizer, extension, and
same-minimizer caches. Piscem currently experiments with the filter only.

These changes were motivated by *Batch, Prefetch, Skip, Orient, Squeeze: Faster,
Smaller k-mer dictionaries via the Spectral Burrows-Wheeler transform* (Diseth,
Alanko, Mackenzie, and Puglisi, supplied preprint 2026.11.09.765318v1), and an audit
of [sbwt-squeeze-rs at ae450089](https://github.com/anadis504/sbwt-squeeze-rs/tree/ae4500895a5ddc3e70281632332738190aeb174c).
The comparison here is against **this Rust SSHash implementation**, not upstream
C++ SSHash. Reproduction programs and recorded measurements are in
[benchmarks/query-campaign](../benchmarks/query-campaign/README.md).

## API and lifecycle

```toml
[dependencies]
sshash-lib = { path = "../sshash-rs/crates/sshash-lib", features = ["experimental-query"] }
```

- `dict.enable_bucket_fingerprints()` builds the table and returns its payload
  size in bytes. Call it before sharing/borrowing the dictionary for queries.
  Point, batch, and streaming lookups then consult the same filter.
- `dict.disable_bucket_fingerprints()` drops it; disabling twice is harmless.
- `dict.query_batch::<K, B, PREFETCH>(&keys)` returns a `[LookupResult; B]` in
  input order. Every slot is real, including padding in a partial final batch;
  callers discard padded outputs or handle the remainder scalarly. Duplicates
  and both strands are supported. An empty batch is allowed but still requires
  `K == dict.k()`. Intermediates use O(B) stack space and no per-call heap.

Public rustdoc includes compiling load/query examples. Odd k in 3..=63 is
supported as in the rest of the library. This API does not schedule reads or
replace the existing streaming API. Batch sizes 16–64 are starting points for
measurement; large arrays can exhaust the stack.

Loading or building a dictionary leaves its filter disabled. Enabling again
rebuilds it sequentially, retaining the old allocation until the replacement is
ready. Retained payload is one byte per bucket, included in `num_bits()`; allocator
and struct overhead are not included. The filter is derived state and is never
written by `save`. No on-disk version bump or rebuild of existing indices is needed.

## Filter representation and correctness

The filter is an `Option<Box<[u8]>>`, addressed by the **existing minimizer MPHF**.
It is not a Bloom, xor, or cuckoo filter and has no separate addressing hash table.
Each nonempty light bucket stores seven bits of a SplitMix64-finalized canonical
minimizer. The builder obtains that minimizer from the bucket's first anchor.
Tags 0..=127 are ordinary fingerprints; 128 means unconditional bypass. Empty and
heavy buckets (more than 64 anchors) receive the bypass marker.

On a query, the existing MPHF produces a candidate bucket ID. An out-of-range ID
is rejected as before. For an in-range light bucket, a different tag proves that
the query's selected minimizer is not this bucket's key. Because a present key
must map to its own MPHF slot and tag, rejecting a mismatch cannot discard a hit.
A matching tag still follows the original bounds/anchor/text verification. Thus
fingerprint collisions affect speed, never membership or result fields.

Heavy buckets deliberately bypass the filter: their existing failure path can
conservatively return `minimizer_found = true`. Rejecting such a query earlier
would change that field even though membership remained negative. Empty buckets
already produce the original exact negative. The shared `lookup_in_bucket`
retains candidate order, orientation selection, and all negative hints.

With well-mixed, independent seven-bit tags, an absent minimizer landing in a
light bucket passes with approximate probability 1/128. That is a modeling
assumption, not an adversarial guarantee of the deterministic mixer. It is also
**not the false-pass rate of all absent k-mers**: an absent k-mer whose selected
minimizer is genuinely indexed always passes, as does a heavy-bucket lookup.
Importantly, a minimizer absent from the selected-minimizer map may still occur
as an unselected substring elsewhere in the reference.

## Batched lookup and prefetch coverage

The batch computes reverse complements and minima, resolves MPHF IDs together,
filters them, obtains bounds, and finally verifies candidates through the scalar
lookup's shared tail. Separating these stages exposes independent memory loads.
`PartitionedMphf::get_batch` uses sux 0.15.2's experimental, doc-hidden PHastR API
for a single partition; multiple partitions and empty maps retain scalar routing.
It does not regroup queries across partitions or change their MPHF addresses.

With `PREFETCH = true`, the additional stages hint:

1. Elias–Fano select block/subblock inventories and the first low-bit word.
2. The first packed anchor-offset word and the following word.
3. The first light-bucket anchor's sequence region and possible preceding starts.

These hints do not fully pipeline the dependent high-bit select loads, later
anchors, or heavy-bucket lookup. `PREFETCH = false` disables the downstream hints;
PHastR's own batch prefetch remains active. sux uses architecture-dependent hints
(x86 T0 in the measured build); unsupported configurations may treat them as
no-ops. No correctness assumption depends on the hinted data being in cache.

## Measurements and limits

The campaign used an EPYC 9555 (64 cores/128 threads), Rust 1.97.1, thin LTO, and
`-C target-cpu=native -C llvm-args=-x86-branches-within-32B-boundaries`. Queries
were pinned to CPU 24. Before each campaign the host passed three 10-second
quiet samples (host busy <4%, benchmark CPU and SMT sibling <5%). Per-trial
telemetry excluded samples with at least two CPU-equivalents of activity outside
the benchmark CPU or >1% sibling activity; reported cells required five accepted
samples. See the raw records for accepted counts and variation.

The human-index confirmation below reverses filter order within the same
experimental binary. Values are median ns/query; point queries materialize full
results from pre-encoded independent seeds, whereas streaming consumes ASCII
windows. Compare columns within each workload, not point and stream as equal work.

| Human workload | Scalar, no filter | Batch 32, no filter | Stream, no filter | Stream, filter |
|---|---:|---:|---:|---:|
| Reference-derived | 833.65 | 524.69 | 51.01 | 51.51 |
| 1% substitutions | 839.64 | 509.12 | 83.31 | 77.59 |
| 5% substitutions | 761.07 | 384.56 | 128.56 | 103.69 |
| Random | 424.59 | 82.55 | 108.96 | 48.97 |
| Real RNA reads | 172.20 | 92.74 | 21.33 | 19.38 |

Filtered batch-32 random seeds measured 45.69 ns. These results do not establish
that batching every successive window is profitable: streaming already reuses
substantial work. Positive streaming saw about 1% overhead in this confirmation
and approximately 10% in an earlier filter-order run; positive-heavy workloads
need their own measurements. The human table uses 177,989,834 bytes (169.7 MiB)
and took roughly 11.5 seconds to construct. Gencode used 13,228,125 bytes (12.6 MiB)
and roughly 0.21 seconds. Construction time is excluded from query timings.

Reference-derived reads sample SPSS strings uniformly, not genomic positions
uniformly. The human real workload uses the first 100k reads of SRR21186103 mate
1. Older Gencode measurements used a five-partition index; the human and freshly
rebuilt Gencode indices use a single MPHF. Do not attribute those index/build
differences to query changes.

A matched-harness, same-build-path ABBA control compared original and candidate
default builds. It did not reproduce an initial apparent 16% default-path
regression caused by comparing different harness revisions. Small mixed shifts
and drift remain; this is not evidence of zero performance cost from refactoring.
The feature-off path preserves results, but performance still needs target-specific
review. See `results/default-comparison.json` for all four rounds.

### Direct SBWT comparison

The controlled synthetic reference has 20M bases, k=31/m=19, 300-k-mer SPSS
segments, alternating stored orientations, and 10k 200-base reads per workload.
One-orientation SBWT uses its own orientation optimizer. Nine repetitions rotate
the five implementations, with an untimed warmup before each measurement.
All 180 timed trials passed the quiet criteria. This harness consumes hit counts
and IDs, unlike the full-result materialization above, so its absolute ns values
must not be combined with the human table.

| Workload | SSHash stream | SSHash filtered | SBWT two orientations | SBWT one orientation | SBWT FBCS |
|---|---:|---:|---:|---:|---:|
| Reference-derived | 8.550 | 8.679 | 8.191 | 9.244 | 23.346 |
| 1% substitutions | 13.452 | 11.936 | 13.067 | 20.667 | 40.542 |
| 5% substitutions | 18.685 | 15.256 | 17.576 | 21.057 | 38.465 |
| Random | 18.187 | 11.646 | 4.247 | 6.678 | 12.845 |

Units are median ns/input window. The fastest SBWT remains about 2.74x faster
than filtered SSHash on random reads here. Filtered SSHash wins on the two
substitution workloads. This is a synthetic comparison, not a human-index SBWT
result or an end-to-end mapping comparison. SBWT colex ranks and SSHash IDs are
different namespaces: membership is compared between structures, while each
implementation's exact positional outputs are checked against its own oracle.

### Downstream mapping

On Gencode v44, nine alternating runs of 100k paired reads, one mapping worker,
serial FASTQ decoding, and explicit `--dict sshash` gave:

| Workload | Original mapping s | Filter mapping s | Original wall s | Filter wall s |
|---|---:|---:|---:|---:|
| Simulated | 0.493 | 0.483 | 0.541 | 0.748 |
| Real RNA | 0.497 | 0.452 | 0.545 | 0.713 |
| Random | 1.327 | 0.617 | 1.374 | 0.883 |

The roughly 9% real-read mapping-time reduction becomes a slower complete short
job because table construction is paid on load. Random mapping improved about
2.15x excluding initialization and 1.56x including it. All 54 RAD outputs were
byte-identical within their workload, with equal semantic mapping metadata.
Runtime/timing metadata naturally differs. No multi-worker speedup is claimed.

## Why SBWT still wins on absent queries

The main transferable lesson is to do less work per overlapping window, not
merely reduce fingerprint collisions. Section 5 of the paper and
[`restart_search.rs`](https://github.com/anadis504/sbwt-squeeze-rs/blob/ae4500895a5ddc3e70281632332738190aeb174c/api/src/restart_search.rs)
use a failed substring search as an absence certificate. With the index's dummy
nodes, a substring occurring within an indexed k-mer has a nonempty search
interval. If a length-l prefix at window start w fails, every window starting in
`[w + l - k, w]` contains that absent substring (clipped to valid read starts).
For k=31 and l=12, that is 20 absent windows from one failure. This example is
illustrative, not a measured average failure length. The eliminated windows
extend **to the left**; the algorithm looks ahead, then fills remaining gaps
right-to-left. It is not permission to skip the next k-1 windows after any miss.

Our filter certifies absence only from the selected-minimizer map. It cannot
safely certify all windows containing those bases; current same-minimizer
memoization already reuses its narrower negative result. SSHash still advances
rolling state per window and pays the MPHF plus tag access at new minima.

SBWT also batches independent reads, separates continuation and restart queues,
prefetches rank work, and uses cache-local rank layouts (see
[`line_rank.rs`](https://github.com/anadis504/sbwt-squeeze-rs/blob/ae4500895a5ddc3e70281632332738190aeb174c/api/src/line_rank.rs)).
Its callbacks report positives sparsely. Our point-batch API is not yet a
cross-read streaming scheduler or a Piscem seed queue. These mechanisms explain
plausible sources of the gap; the campaign did not isolate their individual
contributions with SBWT ablations.

## Tuning and next experiments

The current byte layout is a measured baseline, **not a proven optimum**. Runtime
on/off control is already available in the library. Before adding arbitrary
fingerprint-width knobs, compare a small set of explicit layouts:

| Candidate layout | Bytes/bucket | Ideal conditional collision rate | Tradeoff |
|---|---:|---:|---|
| Current seven-bit tag + bypass | 1 | 1/128 | Simple mask and byte access |
| 255 tag values + one sentinel | 1 | 1/255 | Same memory; uniform reduction costs/bias need measurement |
| Packed nibble, 15 tags + sentinel | 0.5 | 1/15 | Smaller working set; more full verifications |
| 16-bit entry, 65,535 tags + sentinel | 2 | 1/65,535 | Fewer collisions; twice the tag footprint |

Only the first row is implemented. Rates assume truly absent selected minima
landing in light buckets and a uniform independent tag model. Wider tags do
nothing for genuine indexed minima or heavy-bucket bypasses. Seven bits already
reject about 99.22% of eligible misses under that model, so width alone cannot
explain or erase the remaining SBWT gap. A Bloom or xor filter would require
additional addressing/memory work; neither is automatically better when an MPHF
address is already available.

Priorities for further measured work:

1. Instrument a separate diagnostic build to distinguish true minimizer hits,
   tag collisions, heavy bypasses, rejections, and reused negative windows. Avoid
   adding counters to the timed hot path.
2. Compare the layouts above across positive, negative, real, and larger-than-cache
   workloads. Expose only validated choices, keeping the feature default off.
3. Batch independent read seeds while retaining each read's extension/memoization
   state; measure complete mapping time and latency, not just point-query ns.
4. Evaluate bounded-parallel filter construction or a validated derived sidecar
   to amortize startup. A sidecar must verify index identity before use.
5. Explore a separate substring-absence structure only with a sound certificate
   covering all indexed substrings, a full memory budget, and exact skip tests.

For Piscem, a runtime off/on policy and memory budget would be more useful than
the current compile-time research switch. Automatic enabling needs an evidence-based
amortization policy; short jobs in this campaign would make the wrong choice.

## Exactness coverage and known limitation

Feature tests compare complete `LookupResult` values for scalar/batch and
filter-disabled/enabled streaming paths. They cover k=21,31,33,63, even m,
palindromic minima, both strands, heavy buckets, duplicate keys, partial and empty
batches, empty dictionaries, lowercase, N, and short/empty reads. MPHF tests cover
single/multiple partitions, absent keys, and an empty map. Saving before and after
enabling gives identical `.ssi` and `.ssi.mphf` files; a loaded copy starts without
the filter. CI runs these contracts and the public documentation examples.

The campaign also compared explicit little-endian encodings of every result
field, avoiding unspecified Rust struct padding: 1,708,114,192 bytes for the human
workload, 372,101,968 bytes for the older Gencode index, and 1,092,589,648 bytes for
the fresh Gencode original/default/experimental comparison. Hashes and equality
results are retained in `results/validation.json`. These are tested cases, not an
exhaustive proof over every valid index.

The audit found a pre-existing streaming-extension discrepancy: `kmer_offset`
can remain stale while the other hit coordinates advance. This experiment
preserves that behavior; it does not silently repair it. Streaming parity uses
the original streaming implementation as its full-field oracle, and point parity
uses the original point implementation. The SBWT harness reports the discrepancy
explicitly. Any correction belongs in a separate behavior-changing patch.
