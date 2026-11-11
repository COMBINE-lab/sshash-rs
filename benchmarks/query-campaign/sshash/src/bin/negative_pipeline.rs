//! Stateful cross-read seed batching on full streams and actual Piscem request
//! traces. Correctness materializes full fields outside timing; timed paths all
//! consume the same hit-ID/hint digest. Args: INDEX N REAL_FASTX [TRACE...].
use sshash_lib::{DeferredQuery, Dictionary, LookupResult};
use std::{
    hint::black_box,
    time::{Instant, SystemTime, UNIX_EPOCH},
};
const K: usize = 31;
struct Read {
    bases: Vec<u8>,
    positions: Vec<usize>,
    expected: Vec<LookupResult>,
}
fn next(s: &mut u64) -> usize {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*s >> 32) as usize
}
fn whole(bases: Vec<u8>) -> Read {
    let positions = (0..bases.len().saturating_sub(K - 1)).collect();
    Read {
        bases,
        positions,
        expected: Vec::new(),
    }
}
fn load_trace(path: &str) -> Vec<Read> {
    let mut reads = Vec::new();
    let mut r = whole(Vec::new());
    for line in std::fs::read_to_string(path).unwrap().lines() {
        if line == "R" {
            if !r.positions.is_empty() {
                reads.push(r);
            }
            r = whole(Vec::new());
            continue;
        }
        let p: Vec<_> = line.split(',').collect();
        assert_eq!(p[0], "Q");
        let pos: usize = p[1].parse().unwrap();
        let bytes = p[2].as_bytes();
        assert_eq!(bytes.len(), K);
        if r.bases.len() < pos + K {
            r.bases.resize(pos + K, b'N');
        }
        for j in 0..K {
            let old = r.bases[pos + j];
            assert!(old == b'N' || old == bytes[j]);
            r.bases[pos + j] = bytes[j];
        }
        r.positions.push(pos);
        r.expected.push(LookupResult {
            kmer_id: p[3].parse().unwrap(),
            kmer_id_in_string: p[4].parse().unwrap(),
            kmer_offset: p[5].parse().unwrap(),
            kmer_orientation: p[6].parse().unwrap(),
            string_id: p[7].parse().unwrap(),
            string_begin: p[8].parse().unwrap(),
            string_end: p[9].parse().unwrap(),
            minimizer_found: p[10] == "1",
        });
    }
    if !r.positions.is_empty() {
        reads.push(r);
    }
    reads
}
fn scalar(d: &Dictionary, reads: &[Read], mut emit: impl FnMut(usize, usize, LookupResult)) {
    let mut q = d.create_streaming_query::<K>();
    for (rid, r) in reads.iter().enumerate() {
        q.reset();
        for (j, &pos) in r.positions.iter().enumerate() {
            if j != 0 && pos != r.positions[j - 1] + 1 {
                q.reset();
            }
            emit(rid, j, q.lookup(&r.bases[pos..pos + K]));
        }
    }
}
/// Drain local extension/memo work, suspend on actual unresolved dictionary
/// seeds, then resolve all pending lanes. Refill completed lanes with new reads.
/// No result-dependent speculative queries or within-read reordering occur.
fn pipeline<const B: usize, const BATCH: bool, const PREFETCH: bool>(
    d: &Dictionary,
    reads: &[Read],
    mut emit: impl FnMut(usize, usize, LookupResult),
) -> (usize, usize) {
    let mut lanes: [DeferredQuery<K>; B] = std::array::from_fn(|_| DeferredQuery::new(d));
    let mut ids = [usize::MAX; B];
    let mut cursors = [0; B];
    let mut next_read = 0;
    let (mut seeds, mut batches) = (0, 0);
    loop {
        for i in 0..B {
            loop {
                if ids[i] == usize::MAX {
                    if next_read == reads.len() {
                        break;
                    }
                    ids[i] = next_read;
                    next_read += 1;
                    cursors[i] = 0;
                    lanes[i].reset();
                }
                let rid = ids[i];
                let j = cursors[i];
                let r = &reads[rid];
                if j == r.positions.len() {
                    ids[i] = usize::MAX;
                    continue;
                }
                let pos = r.positions[j];
                if j != 0 && pos != r.positions[j - 1] + 1 {
                    lanes[i].reset();
                }
                if let Some(result) = lanes[i].prepare(&r.bases[pos..pos + K]) {
                    emit(rid, j, result);
                    cursors[i] += 1;
                } else {
                    break;
                }
            }
        }
        let pending = std::array::from_fn::<_, B, _>(|i| lanes[i].is_pending());
        let count = pending.iter().filter(|&&p| p).count();
        if count == 0 {
            break;
        }
        seeds += count;
        batches += 1;
        if BATCH {
            d.complete_seed_batch::<K, B, PREFETCH>(&mut lanes);
        } else {
            for q in lanes.iter_mut().filter(|q| q.is_pending()) {
                q.complete_scalar();
            }
        }
        for i in 0..B {
            if pending[i] {
                emit(ids[i], cursors[i], lanes[i].result());
                cursors[i] += 1;
            }
        }
    }
    (seeds, batches)
}
fn execute(
    d: &Dictionary,
    reads: &[Read],
    method: &str,
    emit: impl FnMut(usize, usize, LookupResult),
) -> (usize, usize) {
    match method {
        "scalar" => {
            scalar(d, reads, emit);
            (0, 0)
        }
        "local1" => pipeline::<1, false, false>(d, reads, emit),
        "interleave32" => pipeline::<32, false, false>(d, reads, emit),
        "batch8" => pipeline::<8, true, true>(d, reads, emit),
        "batch16" => pipeline::<16, true, true>(d, reads, emit),
        "batch32" => pipeline::<32, true, true>(d, reads, emit),
        "batch64" => pipeline::<64, true, true>(d, reads, emit),
        "batch32-no-prefetch" => pipeline::<32, true, false>(d, reads, emit),
        _ => unreachable!(),
    }
}
fn digest(d: &Dictionary, reads: &[Read], method: &str) -> (u64, u64) {
    let (mut found, mut sum) = (0u64, 0u64);
    execute(d, reads, method, |_, _, r| {
        if r.is_found() {
            found += 1;
            sum = sum.wrapping_add(r.kmer_id);
        }
        sum = sum.wrapping_add(r.minimizer_found as u64);
    });
    black_box((found, sum))
}
fn quiet() {
    if let Some(p) = std::env::var_os("QUERY_QUIET_FILE") {
        while std::fs::read_to_string(&p).ok().as_deref() != Some("ready") {
            std::thread::sleep(std::time::Duration::from_millis(100));
        }
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut d = Dictionary::load(&args[1]).unwrap();
    assert_eq!(d.k(), K);
    let n: usize = args[2].parse().unwrap();
    let reps: usize = std::env::var("REPS").unwrap_or("9".into()).parse().unwrap();
    let mut rng = 0x543971;
    let mut positive = Vec::new();
    for i in 0..n {
        let sid = next(&mut rng) % d.num_strings() as usize;
        let begin = d.spss().offset_at(sid) as usize;
        let end = d.spss().offset_at(sid + 1) as usize;
        let len = (end - begin).min(200);
        let start = begin + next(&mut rng) % (end - begin - len + 1);
        let text = d.spss().strings_data();
        let mut r: Vec<_> = (start..start + len)
            .map(|p| b"ACTG"[((text[p / 4] >> (2 * (p % 4))) & 3) as usize])
            .collect();
        if i % 2 == 1 {
            r.reverse();
            for b in &mut r {
                *b = match *b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                };
            }
        }
        positive.push(r);
    }
    let mut workloads: Vec<(String, Vec<Read>)> = vec![(
        "positive".into(),
        positive.iter().cloned().map(whole).collect(),
    )];
    for rate in [1, 5, 100] {
        let reads = positive
            .iter()
            .map(|r| {
                whole(
                    r.iter()
                        .map(|&b| {
                            if next(&mut rng) % 100 >= rate {
                                return b;
                            }
                            if rate == 100 {
                                b"ACGT"[next(&mut rng) % 4]
                            } else {
                                let i = b"ACGT".iter().position(|&x| x == b).unwrap();
                                b"ACGT"[(i + 1 + next(&mut rng) % 3) % 4]
                            }
                        })
                        .collect(),
                )
            })
            .collect();
        workloads.push((format!("sub{rate}"), reads));
    }
    let mut reader = needletail::parse_fastx_file(&args[3]).unwrap();
    let mut real = Vec::new();
    while real.len() < n {
        let Some(r) = reader.next() else {
            break;
        };
        real.push(whole(r.unwrap().seq().to_ascii_uppercase()));
    }
    workloads.push(("real".into(), real));
    for path in &args[4..] {
        workloads.push((
            std::path::Path::new(path)
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .into(),
            load_trace(path),
        ));
    }
    // Baseline vectors are materialized once. Traces also carry the mapper's
    // independently captured full results, checked before any filter is enabled.
    for (name, reads) in &mut workloads {
        let mut want: Vec<Vec<_>> = reads
            .iter()
            .map(|r| vec![LookupResult::not_found(); r.positions.len()])
            .collect();
        scalar(&d, reads, |i, j, r| want[i][j] = r);
        for (i, r) in reads.iter_mut().enumerate() {
            if !r.expected.is_empty() {
                assert_eq!(r.expected, want[i], "trace {name},read {i}");
            }
            r.expected = std::mem::take(&mut want[i]);
        }
        println!(
            "WORKLOAD,{name},{},{}",
            reads.len(),
            reads.iter().map(|r| r.positions.len()).sum::<usize>()
        );
    }
    let methods = [
        "scalar",
        "local1",
        "interleave32",
        "batch8",
        "batch16",
        "batch32",
        "batch64",
        "batch32-no-prefetch",
    ];
    for (round, filtered) in [false, true, true, false].into_iter().enumerate() {
        quiet();
        if filtered {
            d.enable_bucket_fingerprints();
        } else {
            d.disable_bucket_fingerprints();
        }
        let layout = if filtered { "bits7" } else { "off" };
        for (name, reads) in &workloads {
            for method in methods {
                let mut counts = vec![0; reads.len()];
                let work = execute(&d, reads, method, |i, j, r| {
                    assert_eq!(j, counts[i]);
                    assert_eq!(r, reads[i].expected[j], "{name}/{method}/{i}/{j}");
                    counts[i] += 1;
                });
                assert!(
                    counts
                        .iter()
                        .zip(reads)
                        .all(|(n, r)| *n == r.positions.len())
                );
                println!(
                    "EXACT,{round},{layout},{name},{method},seeds={},batches={}",
                    work.0, work.1
                );
            }
        }
        for rep in 0..reps {
            for i in 0..workloads.len() {
                let (name, reads) = &workloads[(i + rep) % workloads.len()];
                let expected = digest(&d, reads, "scalar");
                let count: usize = reads.iter().map(|r| r.positions.len()).sum();
                if count == 0 {
                    continue;
                }
                for j in 0..methods.len() {
                    quiet();
                    let method = methods[(j + rep) % methods.len()];
                    assert_eq!(digest(&d, reads, method), expected);
                    let wall = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    let now = Instant::now();
                    // Aggregate a fixed number of query windows per trial to
                    // reduce sensitivity to very short timing intervals.
                    let inner = 5_000_000_usize.div_ceil(count).max(1);
                    for _ in 0..inner {
                        assert_eq!(digest(&d, black_box(reads), method), expected);
                    }
                    let ns = now.elapsed().as_nanos() as f64 / (count * inner) as f64;
                    let end = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    println!(
                        "LAYOUT,{round},{layout},{name},{method},{rep},{count},{ns},{wall},{end}"
                    );
                }
            }
        }
    }
}
