# Query campaign reproduction

These two standalone Cargo packages preserve the measurement loops from the
2026-11-11 campaign. Dependency paths were made portable, comments expanded, and
Rust formatting applied after measurement. The library's sequence-prefetch loop
was also rewritten as an equivalent iterator to satisfy Clippy; final contract
tests were rerun, but the campaign timings precede that stylistic cleanup.
The harnesses are outside the production
workspace and each has its own lockfile. Do not time debug builds.

Read the [design and interpretation notes](../../docs/experimental-query.md)
before comparing results. The follow-up filter harness also requires the pinned
external clones documented under Stateful negative-query campaign below. Input datasets and large binary result vectors are not
committed. The SBWT harness generates its own deterministic data and can reproduce
the direct comparison without downloading a biological reference.

## Builds

Run commands below from the **sshash-rs repository root**. Use identical flags,
compiler, harness source, allocator, and target build path for baseline/candidate
comparisons. The recorded run used Rust 1.97.1 and these release flags:

```bash
export RUSTFLAGS='-C target-cpu=native -C llvm-args=-x86-branches-within-32B-boundaries'
cargo build --manifest-path benchmarks/query-campaign/sshash/Cargo.toml \
  --release --locked --features experiments
```

For the direct SBWT experiment, place the audited repository beside sshash-rs:

```bash
git clone https://github.com/anadis504/sbwt-squeeze-rs.git ../sbwt-squeeze-rs
git -C ../sbwt-squeeze-rs checkout ae4500895a5ddc3e70281632332738190aeb174c
cargo build --manifest-path benchmarks/query-campaign/sbwt/Cargo.toml --release --locked
```

The SBWT package patches `simple-sds-sbwt` to the audited repository's bundled
source and pins sux/wide to the measured versions. The locks are part of the
experiment; do not silently update them for a comparison.

## Existing-index SSHash experiment

The program accepts `INDEX_PREFIX OUTPUT_FIELDS [READ_COUNT [READS_FASTX]]`.
The dictionary must have k=31. `OUTPUT_FIELDS` is overwritten; use a new path for
each run. Defaults are 20k sampled reads and nine repetitions. FASTA/FASTQ, including
formats supported by needletail, may supply a final real-read workload.

```bash
# Select an idle CPU and check its SMT sibling before timing. CPU 24 was used
# on the recorded EPYC host; it is not necessarily suitable on another machine.
REPS=9 QUICK=1 FILTER_REVERSE=1 taskset -c 24 \
  benchmarks/query-campaign/sshash/target/release/sshash-query-campaign \
  /data/human/index /tmp/human-filter.fields 100000 /data/reads_1.fastq \
  > /tmp/human-filter.log
```

The generator samples SPSS strings uniformly, then positions within each chosen
string, taking at most 200 bases and alternating read orientation. It creates
reference-derived, 1% substitution, 5% substitution, and independent random-base
reads. Ten independent seed offsets per valid read are pre-encoded outside timing.
Each timed call materializes complete lookup results into a reused vector.
Streaming instead processes every ASCII window and resets state between reads.

`QUICK` selects scalar, batch-32, and streaming; unset it to also measure batches
4/8/16/64 and batch-16 without downstream prefetch. `FILTER_REVERSE` runs filter-on
before filter-off; otherwise the order is off/on. Method order rotates within
each repetition. Exact full-field comparisons run before timings and warm the
paths. The filter is rebuilt once per workload when enabled; its construction
time is reported separately. Avoid `READ_COUNT=0`, zero repetitions, or FASTX
inputs with no valid windows: these are benchmark preconditions, not CLI validation.

The output field file is independent of filter state and method. For each
workload it contains point results followed by streaming results, each prefixed
by a little-endian u64 result count. Every result is eight little-endian u64s:
`kmer_id`, `kmer_id_in_string`, `kmer_offset`, `kmer_orientation`, `string_id`,
`string_begin`, `string_end`, `minimizer_found`. This avoids Rust struct padding.
Compare entire files with `cmp`, not just hit counts. Do not compare streaming
and point vectors directly: see the documented pre-existing offset discrepancy.

Build without `--features experiments` for the feature-off candidate. To compare
the original implementation, use the same harness against commit
`80131cd3806e79d4d74fa7cbceff9e5b8da16a8d`, without the feature. Preserve each binary
before rebuilding, use separate result filenames, and run original/candidate in
ABBA order. The campaign's matched-build control uses exactly this principle;
earlier runs with different harness revisions are not a clean default-path control.

Log records:

- `INDEX`, `READS`, `FILTER`: input dimensions and allocation/construction cost.
- `CHECK_BASE`, `CHECK_EXACT`: result counts and full-field parity assertions.
- `TIME`: workload, filter 0/1, method, repetition, queries, ns/query, wall start,
  wall end. Wall timestamps allow joining with independently recorded CPU telemetry.

## Direct SBWT experiment

```bash
RAYON_NUM_THREADS=4 taskset -c 24 \
  benchmarks/query-campaign/sbwt/target/release/sshash-sbwt-campaign \
  20000000 10000 reorient > /tmp/sshash-sbwt.log
```

Arguments are reference bases (default 200k), reads per workload (default 2k),
and optional `reorient` to enable SBWT's orientation optimizer. Both dictionary
builders use four threads; construction and correctness checks are untimed.
The seeded generator makes 300-k-mer segments with k-1 overlap and alternating
stored strands. Reads are 200 bases, balanced across strands, with 0/1/5/100%
substitution workloads (`sub100` is random in the results table).

The program checks SBWT callback positions and ranks against scalar SBWT searches,
checks SSHash membership against the two-orientation SBWT, and checks filtered
SSHash against unfiltered complete streaming results. Small, empty, ambiguous,
and lowercase reads augment the correctness set. SSHash and SBWT use distinct
ID namespaces; their raw IDs must not be compared across structures.

Timing runs nine repetitions in rotating implementation order, with a separate
untimed warmup before each trial. Two-orientation and one-orientation rank queues
use 64 continuation / 32 restart slots; FBCS uses 32/32. Reverse complements for
the one-orientation queries are prepared outside timing. Timed callbacks consume
counts and rank sums; these only prevent elimination, while exact vector checks
establish correctness separately. SSHash consumes found IDs without materializing
full result vectors. `TIME` reports one trial (its best and median columns are
therefore equal); aggregate across `TRIAL` repetitions for the final median.
`SIZE` reports component sizes with explicit exclusions, not comparable total RSS.

## Quiet-machine protocol and records

The harnesses do **not** gate on system activity themselves. Stop other owned
builds/tests, choose an idle physical core, and monitor `/proc/stat` with the SMT
sibling identified from sysfs. Do not terminate unrelated users' workloads.
For the recorded campaign we required three consecutive 10-second samples with
host busy <4% and CPU 24/88 busy <5%, pinned queries to CPU 24, then sampled once
per second. A trial was rejected if an overlapping sample showed >=2 busy
CPU-equivalents outside CPU 24 or >1% activity on CPU 88. Every reported summary
cell required at least five accepted trials. Avoid running builds alongside tests
or timing; an idle whole-machine average alone can hide a busy SMT sibling.

`results/` preserves these records from the original campaign:

- `provenance.json`: baseline commits, compiler, flags, paper checksum, and measured
  binary checksums. Later documentation/formatting builds have different hashes.
- `samples.csv`: all library trials, including exploratory runs; `quiet_sample`
  records acceptance where reviewed. Blank acceptance is not evidence of quiet.
- `summary.json`: accepted library summaries; the human table uses
  `human-final-reverse`, not a mixture of exploratory runs. `n`, `minimum`, and
  `maximum` expose accepted count and spread.
- `default-comparison.json`: the matched-harness four-round default-path control.
- `sbwt-trials.json`: 180 direct-comparison trials with timestamps and quiet flags.
- `validation.json`: byte-equality checks, sizes, and hashes of full field files.
  Its mapping details refer to `notes/experimental-bucket-fingerprints/results.json`
  in the companion Piscem PR; those records are kept with the downstream change.

The original raw logs, telemetry, datasets, and result vectors are retained in
the local campaign archive, not distributed in this repository. Published records
allow inspection of the reported samples but do not independently reconstruct
the original host's activity. New runs must collect their own telemetry.

## Stateful negative-query campaign (follow-up)

The [negative-query report](../../docs/negative-query-campaign.md) uses corrected
main `644fbdb` plus PR #10, independently of the historical campaign above. In
particular, successful streaming extensions now have corrected absolute offsets.
The legacy full-result contract remains the baseline for deferred batching.

External filter sources are benchmark-only dependencies. Before building this
harness, put these exact revisions under `../audit/` relative to the repo root:

```bash
git clone https://github.com/Malfoy/ZOR ../audit/ZOR
git -C ../audit/ZOR checkout 2967cc4907de2af753a7208ef2c566fe045c8d41
git clone https://github.com/ayazhafiz/xorf ../audit/xorf
git -C ../audit/xorf checkout 60106e9fe010505564468eb159bd4e0ed8af89a6
cargo build --manifest-path benchmarks/query-campaign/sshash/Cargo.toml \
  --release --locked --features experiments
```

The binaries accept the following arguments (k=31, nonempty workloads):

```text
negative_pipeline INDEX READ_COUNT REAL_FASTX [PISCEM_TRACE...]
qmer_opportunity INDEX READ_COUNT REAL_FASTX [PISCEM_TRACE...]
filter_compare INDEX READ_COUNT REAL_FASTX Q KEY_FILE KIND [PISCEM_TRACE...]
```

`REPS` defaults to 9 for the scheduler and 7 for the q-mer/filter studies. Each
trial aggregates at least 5M requests/probes. `QS=12,14,16,18,20` controls oracle
q values. For real filters, first run `KIND=extract` to write all sorted unique
reference keys, then use `fuse8`, `zor4`, `zor5`, `bloom8`, `bloom16`, or `bitmap`.
Q=0 means selected minimizers; the final harness measures both standalone
pending-seed probes and streaming completion through an independent filter;
`mphf7` reproduces the existing MPHF-plus-tag decision for that comparison.
A minimizer filter cannot safely certify absent overlapping substrings. Bitmap
is only for q>0. Q16/Q18 direct-address allocations are 512 MiB/8 GiB; Q20 would
allocate 128 GiB. Extraction also writes 8 bytes per distinct key to disk. Key
files have magic `QKEYS001`, little-endian u64 q/count, then sorted u64 keys.
They do not embed index identity: supplying the matching key file is a harness
precondition. A production sidecar would need index/encoding identity checks.

The existing seven-bit bucket filter stays enabled for scalar/fallback seeds.
`always` tests a suffix only after local streaming work fails; `after-miss` also
requires the previous request to miss. A certificate covers [position,position+
K-q]. Both discard streaming state across skips. Their positive-only callbacks
are an experimental contract, not a replacement for legacy negative hints.

`negative/` contains the actual campaign controllers and untimed capture driver.
They assume the archived workspace layout (`campaign/`, sibling checkouts and
index paths); adapt these paths on another host. Keep the binary per stage: the
initial Gencode full-filter run used the exact source snapshot
`sshash/snapshots/filter_compare-qmer-v1.rs` (copy it to `src/bin/filter_compare.rs`
to rebuild); the subsequent standalone minimizer/human/finalist comparisons use
`sshash/snapshots/filter_compare-v2.rs`. The initial net-minimizer and larger-query workload campaign uses
`sshash/snapshots/filter_compare-v3.rs` and `deferred-v2.rs`. The final outlined
seed-path control uses the current source. Provenance records distinguish those binaries. `sshash/snapshots/deferred-v1.rs` preserves the original deferred implementation
for the first four stages (its hash matches stage 1 provenance). The final binary
adds the documented experimental minimizer-filter completion
method; it does not change ordinary streaming or scalar completion.

Run each controller through the cooperative quiet gate, from the workspace root:

```bash
python3 sshash-negative-campaign/benchmarks/query-campaign/run_filter_quiet.py \
  --logs campaign/logs filter-negative-stage3 \
  python3 campaign/negative-query/run_stage3.py
python3 sshash-negative-campaign/benchmarks/query-campaign/summarize_filters.py \
  campaign/logs campaign/negative-query/all-results.json
```

Use fresh log directories/run labels; the runner refuses existing labels and
the archived controllers assume fresh child-log names. Stop builds/tests before
this command. The runner neither signals nor changes
other users' processes. It gates startup, pauses between trials, and samples CPU
activity every second. The summarizer requires complete telemetry coverage and
rejects any interval overlapping contention. Retain original logs and telemetry
alongside summaries; do not infer quietness from a completed run alone.

Additional CSV records (all untimed except `LAYOUT`):

- `LAYOUT`: round/q, layout, workload, method, repetition, requests, ns/request,
  wall start/end. `probe` measures repeated-vector throughput, not dependent
  latency; its vectors are suffixes/minimizers of baseline pending seeds.
- `OPPORTUNITY`: q, workload, mode, requests, remaining dictionary seeds, filter
  probes, absence certificates, covered requests, locally completed requests.
- `BASELINE_SEEDS`: workload, requests, uncached dictionary seeds.
- `ORACLE_SET`: q, queried distinct keys, present queried keys, HLL estimated
  reference distinct keys, construction-accelerator bytes.
- `ORACLE_BUILD`: scanned SPSS bases, seconds.
- `FILTER_BUILD`: q, kind, distinct keys, actual retained payload bytes, build
  seconds, pre-build RSS, post-build RSS, process HWM (last three in KiB).
  HWM includes dictionary, key vector and workloads; it is not filter-only peak.
- `FILTER_PROBES`: q, kind, workload, probes, truly absent probes, false passes.
  Repeated keys count repeatedly: this is workload-weighted false-positive rate.
- `ZOR_REMAINDER`: q, kind, abandoned main-layer keys, remainder slots. Complete
  ZOR always includes that remainder. Pure mode is never used.
- `EXTRACT`: q, exact distinct keys, extraction seconds, deduplication payload.
- `EXACT`: the passed validation scope. Full-filter runs check every inserted
  key for false negatives, plus every workload's membership/positive fields/order.

Reported build time excludes reference enumeration, key-file loading, input
workload creation, seven-bit fallback setup, and exhaustive correctness checks.
Construction is single-core in these runs. ZOR seed is 69 with default tie/cycle
policy; binary fuse uses xorf's deterministic non-random-fill configuration.
Bloom uses four bits in one 64-bit register with SplitMix hashing; powers of two
round the requested 8/16 bits/key upward, so compare actual bytes, not names.

The follow-up publishes `results/negative/stage*-results.json` with every timing,
quiet flag, summary spread, exact check and untimed diagnostic. `coverage` lists
any missing/insufficient summary cells, including cells with zero accepted
trials. Unlike the historical campaign, `results/negative/raw.json.gz` also
contains raw logs and original gate/CPU telemetry. Recompute acceptance with:

```bash
python3 benchmarks/query-campaign/negative/unpack.py \
  benchmarks/query-campaign/results/negative/raw.json.gz /tmp/negative-raw
python3 benchmarks/query-campaign/summarize_filters.py \
  /tmp/negative-raw /tmp/negative-recomputed.json
```

`stage*-provenance.json` identifies source/binary hashes. Large input files and
key sets remain in the local archive. Stage 5 includes 100,000-read exhaustive
workloads to expose query-working-set sensitivity, plus the original 10,000-pair
mapper traces for net selected-minimizer tests. It compares a pre-MPHF filter
while retaining the existing seven-bit table on fallback; it does not measure
removal of that table, combined minimizer+substring filtering, or a dynamic
workload-selection policy. The default scalar control is paired in every run.

Stage 6 adds `direct` for Q=0 and a `pass` (always-true) ablation. The direct
variant follows the original full-result streaming loop with an outlined
filter-aware seed routine, avoiding a prepare/complete Option boundary. The
pass control compares **all** result fields to the original stream; real filters
compare membership and every positive field. `OPPORTUNITY` rows still describe
the deferred always/after-miss variants, not direct-mode instrumentation.

The final control also prints `HITS`: q, kind, workload, requests, positive
requests, positive-ID digest. This records query mix independently of filtering.
Increasing the prefix from 10k to 100k real reads changes the mapping fraction
as well as working-set size; it is not a controlled cache-only perturbation.
