# Query campaign reproduction

For the subsequent configurable-filter and SBWT skipping/queuing experiments,
see [the follow-up design](../../docs/fingerprint-tuning.md) and the additional
commands at the end of this file. Original measurements below remain historical.

These two standalone Cargo packages preserve the measurement loops from the
2026-11-11 campaign. Dependency paths were made portable, comments expanded, and
Rust formatting applied after measurement. The library's sequence-prefetch loop
was also rewritten as an equivalent iterator to satisfy Clippy; final contract
tests were rerun, but the campaign timings precede that stylistic cleanup.
The harnesses are outside the production
workspace and each has its own lockfile. Do not time debug builds.

Read the [design and interpretation notes](../../docs/experimental-query.md)
before comparing results. Input datasets and large binary result vectors are not
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

The original harnesses do **not** gate on system activity themselves. The filter
follow-up adds an optional cooperative runner, described below. Stop other owned
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

## Fingerprint layout and SBWT ablations

These commands use the same compiler/flags/quiet-host policy as above. Run the
programs one at a time, pin them to an idle physical core, and record CPU telemetry.
The example CPU 24 has SMT sibling 88 on the measured host.

```bash
cargo build --manifest-path benchmarks/query-campaign/sshash/Cargo.toml \
  --bin filter_layouts --features layouts --release --locked
REPS=9 python3 benchmarks/query-campaign/run_filter_quiet.py \
  --logs /tmp/query-logs --cpu 24 filter-layouts \
  benchmarks/query-campaign/sshash/target/release/filter_layouts \
  /data/index 10000 /data/reads_1.fq

cargo build --manifest-path benchmarks/query-campaign/sbwt/Cargo.toml \
  --bin ablation --release --locked
RAYON_NUM_THREADS=4 taskset -c 24 \
  benchmarks/query-campaign/sbwt/target/release/ablation \
  20000000 10000 reorient > /tmp/sbwt-ablation.log
```

`filter_layouts` takes an existing k=31 index, optional read count (default 10k),
and optional real FASTX file. It uses the original campaign's deterministic
SPSS-string sampling and 0/1/5/100% mutation workloads. `DIAG` rows count all
valid-window and minimizer-run-head outcomes before timing; no counters are added
to the library. Each row gives workload, run-head flag, found flag, bucket class,
window count, and pass counts for Packed4/Bits7/Byte255/Wide16. These counts describe
potential filter probes, not actual MPHF calls during successful extension.

`EXACT` confirms complete result-vector equality for scalar, batch-32, and
streaming methods before timing each layout. Point results and streaming results
use their respective unfiltered oracle. Timed streaming consumes IDs; timed point
methods allocate/materialize full result vectors, identically in the control and
candidate builds. Each method receives an untimed warmup. Method and workload
order rotate; layouts run off/7/4/255/16 then 16/255/4/7/off. Construction is untimed
and separately logged as `BUILD`; each transition disables the old table first.
`LAYOUT` records round, layout, workload, method, repetition, query count, ns/query,
and wall start/end for joining with telemetry.

For an original-layout control, build this **same harness source** with
`--features experiments` instead of `layouts`. This selects only off/7/7/off and
compiles without references to the new API, so it also works against PR #10's
library. Point the manifest dependency to that checkout for the original build,
restore it for the candidate, and preserve each binary. Use a balanced order with
`NO_DIAG=1`; the recorded follow-up uses original/prototype/final/final/prototype/
original to distinguish the first configurable prototype from the revised default
fast path. Keep the build path, compiler, flags, and workload arguments identical.
Do not change the production dependency manifest or commit machine-specific paths.

The `ablation` binary retains the original synthetic construction and checks and
adds scalar-no-skip, scalar-look-ahead-skip, and queued-no-skip variants, using
the original queued-look-ahead-skip as the fourth comparison. All eight total
variants rotate over nine repetitions. Additional variants are validated by
complete callback position/rank vectors. Untimed `WORK` rows report workload,
skip flag, characters, continuations, successful continuations, restarts, prefix
table hits, and restart extensions, using the scalar APIs' built-in `COUNT=true`.
Timed paths use `COUNT=false`; counts are not presented as instrumentation of the
queued implementation. The rank representation and prefix table stay fixed in
this comparison; queued APIs also change prefetch and rank-operation code paths.

`summarize_filters.py LOG_DIRECTORY OUTPUT_JSON` retains individual decisions
and aggregates only quiet timings (minimum five per reported cell). It expects
`filter-*.log` files with matching `*-cpu.jsonl` telemetry, using the original
runner schema (`time`, `outside_busy_cpus`, and `sibling.busy`). Without telemetry,
no trial is accepted. Diagnostic and correctness rows are always retained.

The Linux `run_filter_quiet.py` runner identifies SMT siblings from sysfs, waits
for three quiet ten-second samples, pins the child to the requested CPU, and logs
per-second host/sibling activity. It passes `QUERY_QUIET_FILE` to `filter_layouts`;
that harness waits outside validation, warmup, and timing until three consecutive
monitor samples satisfy the timed acceptance criteria. If contention begins
during a trial, the summarizer still rejects it. The runner neither signals nor
changes other processes. Standalone runs without this environment variable do
not pause; the SBWT harness uses the initial gate and post-hoc rejection only.

Use labels beginning `filter-` so the summarizer discovers the logs. A controller
launching several binaries under one runner must associate its common telemetry
with every child log (the recorded controller uses hard links). Build all binaries
before starting the runner; do not compile or test alongside accepted timings.

The follow-up's `results/fingerprint-tuning.json` contains per-trial acceptance,
summary cells, exactness markers, scalar work counters, and filter attribution.
`fingerprint-tuning-provenance.json` distinguishes prototype and final binaries
and records the source revision, hashes, and any measurements still pending.
See [the interpretation](../../docs/fingerprint-tuning.md) before comparing
rows from different harness versions.
