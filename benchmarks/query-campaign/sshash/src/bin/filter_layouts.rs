//! Fingerprint attribution and layout ablation on an existing k=31 index.
//! Diagnostic anchor reads happen before timing and do not instrument hot paths.
//! Arguments: INDEX [READ_COUNT=10000] [REAL_FASTX]. REPS defaults to 9.
//! With `layouts`, tests all layouts in forward/reverse order. Without it this
//! same source compiles against PR #10 for a matched original-layout control.
use sshash_lib::{Dictionary, Kmer, LookupResult};
use std::{collections::BTreeMap, hint::black_box, time::Instant};
use value_traits::slices::SliceByValue;

const K: usize = 31;

/// Fixed RNG keeps sampled strings, substitutions, and query offsets identical.
fn next(s: &mut u64) -> usize {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*s >> 32) as usize
}

fn hash(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

/// Count potential filter outcomes for all valid windows and minimizer-run heads.
/// Run heads are NOT claimed to be actual filter calls: positive streaming can
/// extend without consulting the MPHF. Random-read run heads approximate probes.
fn diagnose(d: &Dictionary, name: &str, reads: &[Vec<u8>]) {
    let mut counts = BTreeMap::<(bool, bool, &str), [u64; 5]>::new();
    for read in reads {
        let mut previous = None;
        for window in read.windows(K) {
            let Ok(key) = Kmer::<K>::from_str(std::str::from_utf8(window).unwrap()) else {
                previous = None;
                continue;
            };
            let mini = d.debug_extract_minimizer(&key).value;
            let run_head = previous != Some(mini);
            previous = Some(mini);
            let found = d.query(&key).is_found();
            let mut passes = [0; 4];
            let class = match d.control_map_ref().lookup(mini) {
                None => "no_slot",
                Some(id) if id >= d.index_ref().num_buckets() => "no_slot",
                Some(id) => {
                    let (begin, end) = d.index_ref().locate_bucket(id);
                    if begin == end {
                        passes = [1; 4];
                        "empty"
                    } else {
                        let pos = d.index_ref().offsets.index_value(begin);
                        let stored = d.spss().decode_mmer_at(pos, d.m());
                        let stored = stored.min(sshash_lib::reverse_complement_mmer(stored, d.m()));
                        if end - begin > 64 {
                            passes = [1; 4];
                            if mini == stored {
                                "heavy_indexed"
                            } else {
                                "heavy_absent"
                            }
                        } else if mini == stored {
                            passes = [1; 4];
                            "light_indexed"
                        } else {
                            let (a, b) = (hash(mini), hash(stored));
                            passes = [
                                (a % 15 == b % 15) as u64,
                                (a & 127 == b & 127) as u64,
                                (a % 255 == b % 255) as u64,
                                (a % 65535 == b % 65535) as u64,
                            ];
                            "light_absent"
                        }
                    }
                }
            };
            let row = counts.entry((run_head, found, class)).or_default();
            row[0] += 1;
            for i in 0..4 {
                row[i + 1] += passes[i];
            }
        }
    }
    for ((head, found, class), row) in counts {
        println!(
            "DIAG,{name},{head},{found},{class},{},{},{},{},{}",
            row[0], row[1], row[2], row[3], row[4]
        );
    }
}

/// Keep complete-result correctness vectors outside timed ID-consumption loops.
fn stream_results(d: &Dictionary, reads: &[Vec<u8>]) -> Vec<LookupResult> {
    let mut q = d.create_streaming_query::<K>();
    let mut out = Vec::new();
    for r in reads {
        q.reset();
        for w in r.windows(K) {
            out.push(q.lookup(w));
        }
    }
    out
}

fn stream(d: &Dictionary, reads: &[Vec<u8>]) -> (usize, u64) {
    let mut q = d.create_streaming_query::<K>();
    let (mut found, mut sum) = (0, 0_u64);
    for r in reads {
        q.reset();
        for w in r.windows(K) {
            let x = q.lookup(w);
            if x.is_found() {
                found += 1;
                sum = sum.wrapping_add(x.kmer_id);
            }
        }
    }
    black_box((found, sum))
}

/// Always materialize full point results; streaming timing instead consumes IDs.
fn points(d: &Dictionary, keys: &[Kmer<K>], batch: bool) -> Vec<LookupResult> {
    let mut out = Vec::with_capacity(keys.len());
    if batch {
        let mut chunks = keys.chunks_exact(32);
        for c in &mut chunks {
            out.extend(d.query_batch::<K, 32, true>(c.try_into().unwrap()));
        }
        out.extend(chunks.remainder().iter().map(|k| d.query(k)));
    } else {
        out.extend(keys.iter().map(|k| d.query(k)));
    }
    out
}

fn enable(d: &mut Dictionary, layout: usize) -> usize {
    d.disable_bucket_fingerprints(); // cap peak allocation during sweeps
    if layout == 0 {
        return 0;
    }
    #[cfg(feature = "layouts")]
    {
        use sshash_lib::BucketFingerprintLayout::*;
        d.enable_bucket_fingerprints_with_layout([Bits7, Bits7, Packed4, Byte255, Wide16][layout])
    }
    #[cfg(not(feature = "layouts"))]
    {
        assert_eq!(layout, 1);
        d.enable_bucket_fingerprints()
    }
}

/// Optional external quiet monitor: wait outside validation, warmup, and timing.
/// The runner atomically writes `ready` after consecutive quiet samples. Missing
/// or busy status waits; with no path, standalone use retains the original mode.
fn wait_for_quiet(path: Option<&std::path::Path>) {
    if let Some(path) = path {
        while std::fs::read_to_string(path).ok().as_deref() != Some("ready") {
            std::thread::sleep(std::time::Duration::from_millis(100));
        }
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let quiet_file = std::env::var_os("QUERY_QUIET_FILE").map(std::path::PathBuf::from);
    let mut d = Dictionary::load(&args[1]).unwrap();
    assert_eq!(d.k(), K);
    let n: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(10000);
    let reps: usize = std::env::var("REPS")
        .ok()
        .map(|x| x.parse().unwrap())
        .unwrap_or(9);
    assert!(n > 0 && reps > 0);
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
    let mut reads = vec![("positive", positive.clone())];
    for rate in [1, 5, 100] {
        let rows = positive
            .iter()
            .map(|r| {
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
                    .collect()
            })
            .collect();
        reads.push((
            match rate {
                1 => "sub1",
                5 => "sub5",
                _ => "random",
            },
            rows,
        ));
    }
    if let Some(path) = args.get(3) {
        let mut reader = needletail::parse_fastx_file(path).unwrap();
        let mut real = Vec::new();
        while real.len() < n {
            let Some(record) = reader.next() else {
                break;
            };
            real.push(record.unwrap().seq().to_ascii_uppercase());
        }
        reads.push(("real", real));
    }
    let mut workloads = Vec::new();
    for (name, reads) in reads {
        if std::env::var_os("NO_DIAG").is_none() {
            diagnose(&d, name, &reads);
        }
        let mut keys = Vec::new();
        for _ in 0..10 {
            for r in &reads {
                if r.len() >= K {
                    let pos = next(&mut rng) % (r.len() - K + 1);
                    if let Ok(key) = Kmer::from_str(std::str::from_utf8(&r[pos..pos + K]).unwrap())
                    {
                        keys.push(key);
                    }
                }
            }
        }
        let point = points(&d, &keys, false);
        let streamed = stream_results(&d, &reads);
        let count = stream(&d, &reads);
        let windows: usize = reads.iter().map(|r| r.len().saturating_sub(K - 1)).sum();
        println!("WORKLOAD,{name},{},{windows},{}", keys.len(), count.0);
        workloads.push((name, reads, keys, point, streamed, count, windows));
    }
    let layouts = if cfg!(feature = "layouts") {
        vec![0, 1, 2, 3, 4, 4, 3, 2, 1, 0]
    } else {
        vec![0, 1, 1, 0]
    };
    for (round, layout) in layouts.into_iter().enumerate() {
        wait_for_quiet(quiet_file.as_deref());
        let t = Instant::now();
        let bytes = enable(&mut d, layout);
        let label = ["off", "bits7", "packed4", "byte255", "wide16"][layout];
        println!(
            "BUILD,{round},{label},{bytes},{}",
            t.elapsed().as_secs_f64()
        );
        for (name, reads, keys, point, streamed, _, _) in &workloads {
            assert_eq!(points(&d, keys, false), *point);
            assert_eq!(points(&d, keys, true), *point);
            assert_eq!(stream_results(&d, reads), *streamed);
            println!("EXACT,{round},{label},{name},all_fields");
        }
        for rep in 0..reps {
            for i in 0..workloads.len() {
                let (name, reads, keys, _, _, want, windows) =
                    &workloads[(i + rep) % workloads.len()];
                for j in 0..3 {
                    wait_for_quiet(quiet_file.as_deref());
                    let method = (j + rep) % 3;
                    // Warm the selected path separately; allocation is included
                    // for point methods in both control and candidate builds.
                    if method == 0 {
                        black_box(stream(&d, reads));
                    } else {
                        black_box(points(&d, keys, method == 2));
                    }
                    let wall = std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    let start = Instant::now();
                    if method == 0 {
                        assert_eq!(stream(&d, black_box(reads)), *want);
                    } else {
                        black_box(points(&d, black_box(keys), method == 2));
                    }
                    let count = if method == 0 { *windows } else { keys.len() };
                    let ns = start.elapsed().as_nanos() as f64 / count as f64;
                    let end = std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    println!(
                        "LAYOUT,{round},{label},{name},{},{rep},{count},{ns},{wall},{end}",
                        ["stream", "scalar", "batch32"][method]
                    );
                }
            }
        }
    }
}
