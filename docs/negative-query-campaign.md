# Negative-query production investigation

This campaign follows the requested order: stateful seed batching, q-mer absence
opportunity, then alternative filter structures. Its corrected baseline is
sshash-rs main `644fbdb` (including the merged streaming-offset fix) plus PR #10's
experimental query code (`b126c74` on this branch). Piscem starts at `9d3bd92`.
The tunable-layout prototype in PR #12 is not the baseline: its matched controls
show a default-path regression on one workload.

## Stage 1: deferred streaming seeds

`DeferredQuery` preserves rolling minimizers, successful extensions, negative
minimizer memos, singleton memos, and decoded anchor caches. `prepare()` returns
a complete result when local work suffices, or suspends on an uncached dictionary
resolution. Repeated heavy buckets retain their already-resolved bounds locally.
`complete_seed_batch()` stages MPHF, tag, EF, anchor and text accesses across
pending lanes and updates each read's original cache and extension state.

A stream borrows one immutable dictionary throughout its lifetime. Reset cancels
pending work; advancing a pending stream panics. Batches reject mixed dictionaries,
ignore nonpending lanes, and support empty/partial batches. The implementation is
feature-gated and mirrors preparation/memo/extension code in a separate module
so the existing streaming hot path remains unchanged. This duplication is an
experimental isolation choice, not a proposed long-term maintenance design.

The scheduler drains each read's local work until it needs a seed, then batches
those seeds and refills finished lanes. No speculative lookups or per-read result
reordering occur. Scalar completion in the same scheduler isolates scheduling
cost; a no-downstream-prefetch variant isolates the extra memory hints. PHastR
batching/prefetch is still present in that ablation.

The harness includes reference-derived reads, 1%/5% substitutions, random reads,
real FASTQ reads, and actual Piscem request traces, captured by [Piscem PR #21](https://github.com/COMBINE-lab/piscem-rs/pull/21).
Captured traces retain resets,
positions, ASCII windows, and every result field. The mapper's skip and revisit
behavior is therefore represented; this is query-component replay, not yet an
end-to-end batched mapper. Trace generation is untimed and single-worker. On
10,000 pairs each of simulated/random/real data, capture produced byte-identical
RAD files against the offset-corrected baseline.

All variants compare every result field and ordering before timing. Timed paths
consume identical hit-ID/hint digests. Trials aggregate at least five million
queries, rotate method/workload order, and use off/on/on/off filter order. The
quiet runner pauses between trials on contention and records per-second CPU/SMT
activity; summaries reject any overlapping busy interval. Build, trace capture,
and correctness checks are outside query timing. Report query-component results
separately from mapper throughput and account for scheduler/lane-state memory.

## Measurement criteria

For q-mer skipping, measure additional seed work avoided after existing streaming
reuse, not just the fraction of absent k-mers. Cover all q-mer substrings of indexed
k-mers canonically, and distinguish a zero-cost absence oracle from a real filter.
Preserve the legacy API: substring absence cannot reproduce `minimizer_found` and
must never be used to poison its negative-minimizer memo. Any sparse positive-hit
interface needs exact positions/metadata and downstream output validation.

Compare complete ZOR (never pure/false-negative mode), binary fuse, register-blocked
Bloom and an exact bitmap where the q-mer universe makes it competitive. Measure
lookup cost, positive overhead, retained payload, peak construction memory, build
time and false passes on the same key/query sets. Production recommendations must
be based on net saved work and mapping behavior, not standalone filter throughput.

## Stage 1 result: do not promote this scheduler

The completed campaign retained 3,731 of 3,744 trials; every summary cell has at
least five accepted samples. Full-field and ordering checks passed for all
methods, widths, filters and workloads. The feature-enabled library suite passed
155 tests (3 ignored); Clippy passed with warnings denied.

Median ns/request with seven-bit fingerprints enabled:

| Workload | Original scalar stream | Best tested batch | Lanes |
|---|---:|---:|---:|
| Gencode real Piscem request trace | 16.65 | 19.77 | 8 |
| Gencode simulated Piscem request trace | 75.31 | 71.59 | 8 |
| Gencode random Piscem request trace | 14.82 | 16.31 | 8 |
| Gencode real exhaustive stream | 10.68 | 13.95 | 8 |
| Human real exhaustive stream | 12.04 | 15.99 | 8 |
| Human random exhaustive stream | 14.88 | 15.48 | 8 |

Without fingerprints, the Gencode random mapper trace improves from 38.74 to
20.05 ns/request with 64 lanes. However, scalar streaming with fingerprints is
already faster at 14.82. Scheduler-only interleaving generally loses, and even
one-lane deferred state has measurable overhead. The simulated mapper trace's
modest gain does not outweigh the real-read regressions. These are component
replay/stream measurements, not end-to-end batched mapping. The production gate
fails, so a broad asynchronous mapper refactor is not justified by this candidate.
The existing independent-seed batch API remains useful for genuinely independent
queries; this result concerns integrating it with stateful read streams.

## Stage 2: exact substring-absence opportunity

The query-specific oracle scans the complete SPSS and records exact canonical
suffix membership for q=12/14/16/18/20. It retains only requested keys, so it is
**not a deployable filter**. A build-side Bloom accelerator cannot affect answers:
all candidates receive an exact key comparison. Cardinality estimates use HLL
(16,384 registers, approximately 0.8% nominal relative standard error).

After local streaming reuse, an absent suffix proves the next K-q+1 consecutive
window starts absent. Skipping discards streaming state, and the study counts
extra seeds caused by that reset. It never installs a negative-minimizer memo.
The experimental callback exposes exact positive results and positions; it does
not reproduce negative `minimizer_found` hints. Legacy APIs remain unchanged.
All positive fields, membership and request order passed differential checks.
1,357 of 1,365 timing trials were accepted; every cell has at least five samples.

| Workload | q | Scalar ns/request | Oracle ns/request | Additional filter budget ns/probe |
|---|---:|---:|---:|---:|
| Gencode random exhaustive | 16 | 13.18 | 4.86 | 106.62 |
| Gencode real Piscem trace | 16 | 16.46 | 12.42 | 37.33 |
| Gencode real Piscem trace | 18 | 15.99 | 11.22 | 43.99 |
| Gencode real exhaustive | 18 | 10.30 | 11.37 | none |
| Human random exhaustive | 18 | 14.54 | 6.05 | 81.84 |
| Human real exhaustive | 18 | 11.62 | 15.73 | none |

For the real Gencode mapper trace, q16 reduces dictionary seeds from 185,915 to
51,343; q18 to 36,969. The budget is `(scalar-oracle)*requests/probes`, an estimate
of additional cost replacing the cheap answer-matrix lookup. It is not a strict
bound: real false positives and cache effects change the execution path. These
oracle timings are optimistic for this implementation, not a universal lower
bound on every possible skip algorithm. The after-miss policy also runs in the
published data; it protects initial positive seeds but forfeits opportunities.

Estimated distinct canonical q16/q18 counts are 82.6M/90.8M for Gencode and
603.6M/1,118.5M for human. Direct-address bitmap payloads are 512 MiB and 8 GiB.
Therefore a substring filter is much larger than a selected-minimizer filter.
Stage 3 must measure the full reference sets and actual filter overhead.

## Stage 3: full-reference filters

The implementation compares upstream complete ZOR with 4/5 hash functions,
binary fuse with 8-bit fingerprints, a four-bit-mask register-blocked Bloom
filter at two requested sizes, and an exact direct-address bitmap. The ZOR
remainder is mandatory: pure/main-only construction would allow false negatives.
Every inserted key is checked after construction. Construction failure aborts
the run; no incomplete filter is accepted.

Pinned sources are [ZOR `2967cc4`](https://github.com/Malfoy/ZOR/tree/2967cc4907de2af753a7208ef2c566fe045c8d41)
and [xorf `60106e9`](https://github.com/ayazhafiz/xorf/tree/60106e9fe010505564468eb159bd4e0ed8af89a6).
The comparison concerns these implementations/configurations, not a proof that
one filter family is universally optimal. No dependency enters the production
library workspace. Binary-fuse construction uses deterministic unused-cell fill;
ZOR uses seed 69 and the default tie-breaking/cycle policy.

The full Gencode key sets contain 83,634,080 canonical q16 keys and 90,881,868
canonical q18 keys. The following q18 results are from the initial full comparison;
scalar controls are paired per backend, because code layout and cache state can
shift the baseline. Payload excludes metadata and the existing dictionary.

| Filter | Payload MiB | Build seconds | Process build HWM GiB | Real Piscem trace scalar → filtered ns/request |
|---|---:|---:|---:|---:|
| Binary fuse 8 | 97.625 | 3.57 | 3.19 | 16.04 → 12.56 |
| Blocked Bloom, requested 8 bits/key | 128 | 1.59 | 1.53 | 16.03 → 12.56 |
| Blocked Bloom, requested 16 bits/key | 256 | 1.68 | 1.66 | 16.20 → 12.84 |
| Complete ZOR 4 | 91.953 | 31.98 | 8.43 | 17.16 → 17.91 |
| Complete ZOR 5 | 89.055 | 34.67 | 9.40 | 16.74 → 17.70 |
| Exact bitmap | 8192 | 0.35 | 9.40 | 16.53 → 12.31 |

Build time excludes reference enumeration, key loading, workload preparation,
fallback-tag setup and correctness checks. HWM is the whole process, including
its dictionary, keys and workloads; pre-build RSS is approximately 1.40 GiB in
these q18 runs. Actual retained filter payload is reported separately. Bitmap
construction time is especially misleading without its allocation size and
reference-key enumeration: it only sets bits from an already sorted key vector.

The q18 ZOR-4 main layer abandoned 2,516,670 keys; ZOR-5 abandoned 1,135,528. Their
complete filters check the main layer OR the remainder, so most absent queries
also pay for remainder lookup. ZOR saved only 6–9% payload versus binary fuse
here, with substantially higher lookup and construction cost. This does not
justify pursuing ZOR as the speed-oriented production candidate.

The filter's key set matters more than its nominal false-positive rate. Q18
substring filters cover roughly 91M Gencode keys, whereas selected-minimizer
filters cover 13,228,125. On the latter set, standalone random-mapper seed probes
cost 10.22 ns for the existing MPHF-plus-seven-bit decision, 2.61 ns for 16 MiB
blocked Bloom, 3.68 ns for 14.22 MiB binary fuse, and 14.51/16.29 ns for ZOR-4/5.
These are repeated-vector throughputs, not dependent latency or net stream
speedups. A pre-MPHF filter still adds a probe before real positive dictionary
seeds; net streaming tests and larger query sets are required below.

### Contracts and cache sensitivity

Substring skipping exposes exact positive results, positions and membership,
not field-identical negative hints. A suffix absence certificate must never be
installed as a negative-minimizer memo. In contrast, an independent filter over
**all selected minimizers** can certify that memo semantically. The experimental
`complete_with_minimizer_filter` method preserves rolling state and clears old
bucket bounds/anchors on rejection. Its predicate must have no false negatives.
It can strengthen a legacy heavy-bucket `minimizer_found=true` miss to false;
therefore this method is explicitly separate from full-field scalar/batch
completion. No arbitrary positive result or SPSS offset can be supplied.

There is no universal win from the q18 prototype. Binary fuse improved synthetic
random streams from 13.23 to 4.83 ns/request, but the larger random mapper trace
only improved 13.94 to 12.50. Its pending-suffix probe throughput was about
23.61 ns on that trace, versus a few ns on smaller working sets. The real mapper
trace improved approximately 22%, but exhaustive real reads slowed 10.31→11.29
and simulated mapper requests slowed 73.74→95.09. After-miss probing reduced the
simulated regression (76.60) while also forfeiting some real-trace gains (13.07).
Reference-derived positive streams also regressed. These are reasons to protect
the existing positive path and test cache footprint, not to enable this prototype
unconditionally or claim end-to-end mapping acceleration.

## Scale checks and replication

A second Gencode run with a separately recorded harness binary confirmed the
q18 real-trace effect: binary fuse 15.94→12.37 ns/request and blocked Bloom
15.91→12.47. The q16 bitmap confirmed at 16.43→13.09. These repetitions establish
the direction on this request trace; they do not turn a replay into mapper timing.

The human reference contains 610,123,646 distinct canonical q16 keys and
1,118,467,309 q18 keys. At q16, exact bitmap payload (512 MiB) is smaller than
binary fuse (654.75 MiB) and blocked Bloom (1024 MiB). At q18 the corresponding
payloads are 8192, 1200 and 2048 MiB. Binary-fuse construction took 28.72/53.67
seconds for q16/q18, excluding enumeration and key loading; q18 process HWM was
32.01 GiB, including about 10.04 GiB of pre-build dictionary/keys/workloads.

| Human workload, 10k reads | q | Scalar → binary fuse ns/request | Scalar → bitmap ns/request |
|---|---:|---:|---:|
| Random exhaustive | 16 | 14.67 → 8.54 | 14.48 → 9.03 |
| Real exhaustive | 16 | 11.55 → 16.88 | 11.56 → 17.07 |
| Random exhaustive | 18 | 14.62 → 5.82 | 14.63 → 6.56 |
| Real exhaustive | 18 | 11.52 → 15.18 | 11.59 → 15.70 |

No actual human-index mapper trace was captured; these are exhaustive streams
against that reference. Construction observations are single runs, not repeated
build-time distributions. The full initial filter comparison accepted 2,678 of
2,688 trials; three ZOR/sub1 cells fell below five accepted samples and are not
used for conclusions. The scale/replication round accepted 2,716 of 2,730; four
cells fell below threshold (Gencode Bloom-q18/sub1 and human Bloom16 standalone
real probes), also excluded. Exact checks passed independently of timing quality.

## Net filtering and larger request sets

The deferred pre-MPHF prototype does not realize the standalone probe advantage.
On the Gencode 10k-pair mapper traces, blocked Bloom changes real requests
16.81→19.57 ns, random requests 13.71→14.34, and simulated requests 74.41→88.12.
Binary fuse also regresses. Positive/exhaustive streams expose substantial wrapper
cost. An outlined-seed control, including an always-pass predicate, is therefore
needed before attributing these losses to the filter data structure itself.

The 100k-read Gencode q-mer runs also undermine extrapolation from small random
workloads: q16 binary fuse is approximately flat on random windows (13.49→13.50),
while q18 binary fuse regresses (13.46→14.28). Q18 real exhaustive requests change
11.33→16.32. This expansion changes both cache footprint and query mix. The q18
filter covers 569,899/1.2M real requests in the 10k prefix (47.5%), but only
2,326,967/12M in the 100k prefix (19.4%). Its standalone pending-suffix throughput
also rises from a few ns to 11.47 ns/probe. These measurements must not be
presented as a pure cache experiment or a workload-independent q-mer speedup.

At 100k reads on human, the deferred prototypes retain a negative-workload win:
blocked Bloom over selected minimizers changes random streams 42.89→33.48 ns;
binary fuse over q18 substrings changes 41.31→19.72. Their payloads are 256 MiB
and 1200 MiB respectively. Both regress on real exhaustive requests (14.35→21.23
and 14.41→19.89). All 1,848 trials in this net/scale round were accepted. Thus
both primary-index size and query mix matter; Gencode and human results must not
be collapsed into one claimed speedup.

The final control keeps the filter in an outlined seed routine and uses a direct
full-result query loop. An always-pass predicate isolates residual experimental
engine/code-generation overhead. On Gencode 10k reads, this control still costs
about 1–2 ns/request versus the original stream (real trace 16.08→18.01; positive
stream 11.97→13.32). It is therefore wrong to infer net benefits from standalone
filter lookup times or to describe this as a zero-overhead production hook.

With actual blocked Bloom, the direct path changes Gencode random exhaustive
queries 13.05→11.62 ns and random mapper requests 13.05→11.81, but the real mapper
trace is essentially flat (16.44→16.38) and simulated requests regress
73.51→85.46. This fixes much of the deferred wrapper loss without establishing a
production win on real mapping. The initial q18 substring filters remain the
stronger measured candidate for that negative-heavy request trace.

The final always-pass controls also record query mix: the first 10k real reads
have 591,730/1.2M positive windows (49.3%), compared with 9,439,781/12M (78.7%) in
the 100k prefix. The actual 10k-pair mapper trace has only 33,640/1,059,809 positive
requests (3.17%), because the mapper's request/skip policy differs from exhaustive
window enumeration. This is why exhaustive-stream regressions cannot alone rule
out a targeted mapper filter, and trace wins cannot establish whole-mapper gains.

On the 100k-read human workload, direct selected-minimizer Bloom improves random
queries 41.77→29.38 ns (29.7%) with 256 MiB payload; binary fuse improves
42.20→32.44 (23.1%) with a smaller payload. Real exhaustive queries still regress
14.27→14.88 (Bloom) and 14.43→15.06 (fuse). The always-pass human control itself
regresses real requests 14.16→15.25, so preserving the original hot-path code
shape remains a necessary integration task. All 1,890 final-control trials were
accepted, and every cell has seven samples.

## Production decision

Keep the existing fixed seven-bit MPHF-indexed tags as the validated baseline.
None of these new prototypes justifies an unconditional default change.

1. **Prioritize q18 binary-fuse absence certificates for a targeted Piscem
   unresolved-seed/positive-hit path.** This is the strongest measured real-trace
   opportunity: approximately 22% lower component query time for about 97.6 MiB
   additional Gencode payload. It offers about 2.1× throughput on large human
   random streams, at 1200 MiB. It must protect positive streaming reuse and
   preserve hit positions, orientation, IDs and offsets. Build/persist the filter
   outside timed mapping, attach it to the correct index identity, and validate
   end-to-end RAD bytes and total wall time before promotion. No such mapper
   throughput gain or new-filter RAD validation is claimed by this PR.
2. **Keep direct pre-MPHF blocked Bloom as the lower-memory alternative.** It
   needs 16 MiB on Gencode and 256 MiB on human in these configurations. Its
   approximately 30% large-human random-query gain is useful, but Gencode real
   mapper requests are flat and positive-heavy regressions remain. The always-pass
   control identifies avoidable integration overhead; standalone 2–4 ns probes
   alone are not production evidence. Binary fuse offers a smaller selected-key
   payload but was slower on the large human negative workload tested here.
3. **Use exact bitmap only where universe density makes it competitive.** Human
   q16 is such a case (512 MiB versus 654.75 MiB binary fuse). Gencode q18 is not:
   8 GiB buys only a small trace-time advantage over a roughly 98 MiB filter.
4. **Do not pursue the tested cross-read scheduler or complete ZOR as the next
   speed-oriented production change.** The scheduler loses on filtered real
   streams; ZOR's small space saving is outweighed by query/build cost. Independent
   point batching remains available for callers that actually have independent
   seeds. These conclusions do not rule out every possible scheduler, ZOR
   specialization, or multi-threaded construction strategy.

The SBWT lesson is chiefly to avoid searching windows that an absent substring
already excludes. A selected-minimizer tag rejects a dictionary route; it cannot
supply that substring certificate. Adding the certificate to SSHash costs a much
larger key set and can discard streaming reuse, which explains why copying a
filter alone does not reproduce SBWT's negative-query advantage. This campaign
compares against our Rust streaming implementation with its existing memos and
seven-bit filter, not upstream C++ SSHash.

## Validation and limits

Across all six measurement rounds, 14,220 of 14,265 trials passed the recorded
quiet-host criteria. Seven exploratory summary cells were below the five-sample
reporting threshold and were excluded from conclusions; all remaining cells have
sufficient samples. Every source/binary snapshot and input/key-file hash is
recorded. Unpacking `results/negative/raw.json.gz` and rerunning the summarizer
reproduced every timing and acceptance decision exactly.

The final feature-enabled library/integration/doc suite passed 156 tests (3
ignored), including both strands, K21/K31/K63, heavy buckets, resets, ambiguous
and lowercase inputs, partial/cancelled batches, selected-minimizer rejection,
and full-field always-pass direct controls. The release filter-harness test also
passed. Standard-target Clippy passed with `-D warnings`. A native-target Clippy
attempt exposed an existing `needless_return` in `darray.rs:263`; that unrelated
code was not changed. Native-target correctness tests and release builds passed.

Deferred scalar/batch completion preserves every result field. The q-mer and
independent-minimizer research paths check membership and every positive field,
plus order. They deliberately do not claim legacy negative-hint equality. Every
inserted key in every complete filter was tested for false negatives. Three
untimed Piscem trace captures produced byte-identical RAD outputs against the
corrected baseline; this validates capture, not an integrated new-filter mapper.

Measurements are single-core on this EPYC host with repeated warmed query sets.
They do not establish multi-threaded throughput, cold-start latency, a universal
q/width choice, or a production workload-selection policy. Reference-derived
reads are sampled by SPSS string, which biases lengths; real mapper traces and
larger real prefixes are reported separately to expose that limitation. The
100k-prefix comparison changes both working set and query mix. No combined
minimizer-plus-q-mer cascade, compressed canonical bitmap, or adaptive policy was
measured. The APIs/code duplication remain experimental and feature-gated.
