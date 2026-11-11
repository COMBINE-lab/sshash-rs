# Fingerprint layouts and the SBWT negative-query gap

This follow-up adds explicit fingerprint layout choices behind the existing
`experimental-query` feature and extends the campaign with filter attribution
and a scalar/queued × no-skip/look-ahead-skip SBWT comparison. It is based on the
original experiment, [PR #10](https://github.com/COMBINE-lab/sshash-rs/pull/10);
the independent streaming-offset fix is
[PR #11](https://github.com/COMBINE-lab/sshash-rs/pull/11).
The baseline streaming behavior is preserved in these optimization comparisons.

## Configuration

```rust,ignore
use sshash_lib::BucketFingerprintLayout;

let bytes = dict.enable_bucket_fingerprints_with_layout(
    BucketFingerprintLayout::Packed4,
);
// Point, batch, and streaming queries now use the selected layout.
dict.disable_bucket_fingerprints();
```

The existing no-argument `enable_bucket_fingerprints()` continues to select
`Bits7`. Merely compiling the feature does not enable or allocate a filter.
This is an explicit runtime policy on the library dictionary. It does not add
Piscem CLI flags or a workload-dependent automatic policy.

| Layout | Payload for n buckets | Normal tag values | Bypass value | Conditional collision model |
|---|---:|---:|---:|---:|
| `Packed4` | ceil(n/2) bytes | 0–14 | 15 | 1/15 |
| `Bits7` (default) | n bytes | 0–127 | 128 | 1/128 |
| `Byte255` | n bytes | 0–254 | 255 | 1/255 |
| `Wide16` | 2n bytes | 0–65,534 | 65,535 | 1/65,535 |

All layouts mix the canonical minimizer with the same SplitMix64 finalizer.
Bits7 retains the original low-seven-bit mask. Other layouts use modulo reduction
to reserve exactly one sentinel. They are deterministic, not cryptographic; the
collision figures assume well-mixed absent selected minima landing in light
buckets, not arbitrary absent k-mers. The query path dispatches on the selected
layout. Reducing collisions can lose time through arithmetic, dispatch, and cache
footprint, which is why the campaign includes an original-versus-candidate control.

All layouts use the original boxed byte slice; Wide16 encodes tags as
little-endian byte pairs. A separate layout field selects the interpretation.
Packed4's final high nibble is unused when n is odd. Heavy and empty
buckets retain bypass semantics, preserving every existing lookup field. A tag
match always receives exact verification. Fingerprint collisions cannot produce
false hits or false negatives.

Building a layout visits the index once, allocates only the requested payload,
and keeps any old table alive until replacement. The benchmark disables the old
table first to cap peak memory; applications replacing a live table should budget
for both allocations. `num_bits()` includes payload bytes, excluding layout
metadata and allocator overhead. Tables are not serialized; loaded dictionaries start unfiltered.

The data structure is a static fingerprint array addressed by the existing
minimizer MPHF. It is not a separate Bloom, xor, or cuckoo filter. Reusing the
address avoids another membership structure and extra independent random
accesses, but a rejected query still pays for minimizer extraction, MPHF lookup,
and the tag read. A stronger fingerprint only saves downstream bucket decoding
and verification on the small fraction that would collide.

This is not a proof of optimality. A separate approximate set could reject before
the MPHF, or tags could be placed closer to MPHF metadata, but either changes the
access and construction tradeoffs. Neither alternative was benchmarked here.
The four layouts test a narrower question: how much tag width helps once the
MPHF address is already known.

For the measured indices, payload bytes are:

| Index | Buckets | Packed4 | Bits7 / Byte255 | Wide16 |
|---|---:|---:|---:|---:|
| Gencode | 13,228,125 | 6,614,063 | 13,228,125 | 26,456,250 |
| Human | 177,989,834 | 88,994,917 | 177,989,834 | 355,979,668 |

## Attribution and measurement contracts

The `filter_layouts` harness uses existing public index accessors to decode the
candidate bucket's exact minimizer before timing. Each valid window is classified
as an invalid slot, empty bucket, heavy bucket with indexed/absent minimizer, or
light bucket with indexed/absent minimizer. For absent light-bucket minima it
computes whether each layout would pass. It also separates positive/negative
k-mers and consecutive-minimizer run heads. This distinguishes genuine indexed
minima and bypasses from fingerprint collisions without adding hot-path counters.

Run heads are potential bucket resolutions, **not measured filter calls**:
positive streaming can extend across a change of minimizer without touching the
MPHF. The distinction matters when relating these counts to timings. A genuine
indexed minimizer and a heavy bypass pass regardless of fingerprint width.

The harness validates all fields against unfiltered point and streaming vectors
before measuring each layout. Timed streaming consumes found IDs; point and
batch methods materialize full result vectors. Thus point/stream absolute timings
represent different work. Methods/workloads rotate within nine repetitions;
layout order runs forward and reverse. Construction is outside query timing and
reported separately. An identical-harness original-layout control detects costs
introduced by configurable storage.

The SBWT harness holds the two-orientation index, prefix table, rank layout,
reads, and positive-result consumer fixed while comparing four APIs:

| | No look-ahead skipping | Look-ahead skipping |
|---|---|---|
| Scalar | `restart_kmers` | `restart_kmers_fskip` |
| Queued | `restart_kmers_queued` | `restart_kmers_queued_fskip` |

The queued APIs combine scheduling, prefetching, and their rank-operation code
paths; this comparison cannot attribute their gain to scheduling alone. Scalar
`COUNT=true` runs measure restart/continuation work outside timing, while timed
functions use `COUNT=false`. Work counts describe the scalar algorithms; queued
work counts are not inferred to be exactly identical. Every additional API's
callback positions and ranks are compared against independent scalar searches.

## SBWT ablation results

On the original campaign's 20M-base synthetic reference (k=31, m=19 for SSHash),
10k reads of 200 bases per workload, and pinned EPYC 9555 core, median ns/input
window were:

| Workload | Scalar, no skip | Scalar, look-ahead skip | Queued, no skip | Queued, look-ahead skip |
|---|---:|---:|---:|---:|
| Reference-derived | 17.413 | 17.415 | 7.655 | 8.407 |
| 1% substitutions | 36.198 | 22.325 | 24.091 | 13.157 |
| 5% substitutions | 61.164 | 24.782 | 47.478 | 17.659 |
| Random | 50.093 | 7.685 | 37.941 | 4.277 |

Each cell has 6–8 accepted trials from nine rotating repetitions. The full eight-
variant run recorded 288 trials; 228 passed the host/sibling telemetry criteria.
All additional APIs passed complete positional colex-rank comparisons. These are
synthetic results, not a human-index SBWT comparison.

On random reads, skipping improves the scalar implementation 6.52x and the queued
implementation 8.87x. Queuing's combined scheduling/prefetch/code-path effect is
1.32x without skipping and 1.80x with it. Thus the negative-query advantage depends
strongly on the stronger absence certificate and the work it eliminates.

The untimed scalar work counters give a concrete explanation:

| Random-read work | No skip | Look-ahead skip |
|---|---:|---:|
| Input windows | 1,700,000 | 1,700,000 |
| Searches from scratch | 1,700,000 | 139,364 |
| Restart character extensions | 9,426,594 | 770,822 |

There are no successful hits or continuations in this random workload. Searches
fall by 12.20x and character extensions by 12.23x. With the eight-base prefix table,
the skipping version's average failed search length is about 13.53 bases. The
gain comes from doing far fewer searches over overlapping windows; the average
work within each search changes little. The API implementations also differ in
traversal order and sparse callback handling, so the measured ratios describe
those implementations together, not an isolated machine instruction.

## Filter outcome attribution

The following deterministic diagnostic counts use 10k sampled random reads from
the existing Gencode and human indices. Read lengths follow sampled SPSS string
lengths, capped at 200; these are different workloads from the fixed-length SBWT
experiment. Counts are **negative minimizer-run heads**, not timed actual probes.

| Outcome | Gencode | Human |
|---|---:|---:|
| Total negative run heads | 68,397 | 42,389 |
| Truly absent minimum, light bucket | 68,369 | 42,083 |
| Genuinely indexed minimum, light bucket | 27 | 295 |
| Absent minimum, heavy bypass | 1 | 11 |
| Seven-bit collisions among absent light minima | 504 | 320 |
| Byte255 collisions among absent light minima | 254 | 139 |
| Wide16 collisions among absent light minima | 3 | 0 |
| Packed4 collisions among absent light minima | 4,582 | 2,808 |

On the human workload, Bits7 passes 626 of these 42,389 potential probes: 320
collisions, 295 genuine indexed minima, and 11 bypasses. Eliminating every collision
would remove only 0.755% of the original potential probes. More tag bits cannot
eliminate the MPHF/tag access already paid for the other queries. Packed4's measured
conditional collision fraction is 6.67%; Bits7's is 0.760%, close to their models.

Real-read negatives have a different composition. On the human real-read sample,
8,481 of 49,536 negative run heads have genuinely indexed minima (including heavy
buckets), versus only 309 seven-bit collisions among absent light minima. Width
tuning has particularly limited scope there. These diagnostic statements concern
the selected-minimizer map; they do not imply substring absence in the reference.

## Adapting stronger absence certificates

An SBWT failed substring search can certify many overlapping windows absent.
For a k-mer search starting at w that fails after l bases, starts from
`w + l - k` through `w` contain the failed substring and are absent (clipped to
the read). For k=31 and l=12, that is up to 20 windows. These windows lie to the
left of the searched window: the implementation looks ahead, then searches any
uncovered gaps backwards. It is not an unconditional forward jump after a miss.
See the pinned [SBWT implementation](https://github.com/anadis504/sbwt-squeeze-rs/blob/ae4500895a5ddc3e70281632332738190aeb174c/api/src/restart_search.rs#L2091).

A selected-minimizer fingerprint mismatch only certifies absence from the
selected-minimizer map; those bases can still occur unselected in indexed strings.
It therefore cannot justify SBWT's broader skip.

A separate q-mer absence structure would have to cover **all** q-mer substrings
of indexed k-mers on both orientations. Its construction time, memory, and query
cost belong in any comparison. Smaller q makes absence rarer; larger q increases
the universe or the size of an approximate set. Skipping also interacts with the
API: the existing per-window result includes a conservative `minimizer_found`
hint. Proving k-mer absence alone does not reproduce that hint, particularly
when the current lookup routes an absent minimum to a heavy bucket. More than a
field-equality issue is involved: `StreamingQuery::seed_with_dict` uses a false
hint to reject subsequent k-mers carrying the same minimizer. Setting it false
from a substring-only absence certificate could poison that memo and suppress a
later genuine hit with an indexed minimizer. The certificate must not be treated
as proof that the selected minimizer is missing.

The current Piscem mapper does not directly branch on this hint; SSHash uses it
inside its streaming state. A new bulk interface reporting positive positions
and exact hit metadata could avoid returning negative hints, with separate state
rules for skipping. Substituting it into Piscem still requires an explicit
ordering/lookup contract and downstream output validation.

The existing rolling stream cannot simply be replaced with independent batches:
retain successful extension and same-minimizer reuse, then batch the independent
reads that need new seeds. This is a separate scheduler experiment from changing
the fingerprint layout.

## Final layout and control measurements

The palindrome original/prototype/final/final/prototype/original control has
32–36 accepted trials per pooled cell. Relative to PR #10, the final seven-bit
streaming path was 12.720 vs 12.781 ns/window on random reads and 9.784 vs 9.800
on real reads, but 20.387 vs 19.598 on 5%-substitution reads (4.0% slower).
Configurability is therefore not established as cost-free. These observations
support keeping the simpler original layout as the production baseline while
retaining this API as an opt-in experiment.

Final Gencode streaming medians, ns/input window (18 accepted trials per cell):

| Layout | Random | Real | Payload bytes |
|---|---:|---:|---:|
| Off | 25.916 | 11.775 | 0 |
| Bits7 | 12.922 | 9.848 | 13,228,125 |
| Packed4 | 13.640 | 10.175 | 6,614,063 |
| Byte255 | 12.939 | 9.946 | 13,228,125 |
| Wide16 | 13.176 | 9.962 | 26,456,250 |

Packed4 halves tag memory at a 5.6% random-stream and 3.3% real-stream cost on
this index. Wider tags did not improve either workload. Mutation workloads do
not rank layouts uniformly, so these are workload-specific observations.

The human sweep completed all full-field checks and has at least five accepted
trials per pooled cell. However, individual timings vary by up to roughly 3x,
and the reverse Bits7 round retains only one random-stream trial. Pooled human
medians therefore do **not** establish a reliable width ranking. The records
retain each round separately; quiet CPU telemetry is necessary but insufficient
to establish stable memory-system timing. Subsequent experiments aggregate at
least five million query windows per trial to reduce short-interval noise.

## Evidence and validation

[Trial records and attribution](../benchmarks/query-campaign/results/fingerprint-tuning.json)
and [source/binary provenance](../benchmarks/query-campaign/results/fingerprint-tuning-provenance.json)
retain the completed SBWT experiment, matched controls, final layout sweeps,
and both indices' diagnostics. Full-field equality passed in every sweep.
The SBWT harness's SSHash comparator was the first configurable prototype;
its SSHash rows do not establish the final candidate's performance. The four
SBWT APIs are unchanged by that library revision.

All feature-enabled library tests (152 passed, 3 ignored) and Clippy with
warnings denied pass; the default workspace suite passed (154, 3 ignored)
before the feature-only storage revision. New Rust files pass formatting;
the workspace-wide check reports existing differences outside this change.
No new Piscem mapping measurements are included. This remains a research draft
because the default-path overhead and human timing stability are unresolved.
