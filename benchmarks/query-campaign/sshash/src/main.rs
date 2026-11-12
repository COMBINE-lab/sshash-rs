//! Reproduce independent-seed and streaming ablations on an existing k=31 index.
//! See ../README.md (one level above this Cargo package) for arguments, workload
//! sampling, timing units, environment controls, and the field-file schema.
//! Full-field equality checks warm each path before rotating timed trials.
//! This is a research harness: malformed arguments/inputs fail fast with panics.
use sshash_lib::streaming_query::LookupResult;
use sshash_lib::{Dictionary, Kmer};
use std::{hint::black_box, io::Write, time::Instant};
const K: usize = 31;
/// Deterministic LCG for repeatable workload/offset selection, not cryptography.
fn next(s: &mut u64) -> usize {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*s >> 32) as usize
}
/// Reverse complement a generated ACGT read, preserving any other byte.
fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => *x,
        })
        .collect()
}
/// Explicit field order avoids unspecified padding in raw LookupResult memory.
fn encode(x: LookupResult) -> [u64; 8] {
    [
        x.kmer_id,
        x.kmer_id_in_string,
        x.kmer_offset,
        x.kmer_orientation as u64,
        x.string_id,
        x.string_begin,
        x.string_end,
        x.minimizer_found as u64,
    ]
}
/// Reuse output capacity; reset rolling state between reads, not windows.
fn stream(d: &Dictionary, reads: &[Vec<u8>], out: &mut Vec<LookupResult>) {
    out.clear();
    let mut q = d.create_streaming_query::<K>();
    for r in reads {
        q.reset();
        for w in r.windows(K) {
            out.push(q.lookup(w));
        }
    }
}
/// Evaluate pre-encoded independent seeds and materialize every result field.
fn points(d: &Dictionary, kmers: &[Kmer<K>], out: &mut Vec<LookupResult>, method: usize) {
    out.clear();
    #[cfg(feature = "experiments")]
    {
        // Exact full batches plus a scalar remainder avoid counting padded work.
        fn batch<const B: usize, const P: bool>(
            d: &Dictionary,
            k: &[Kmer<K>],
            out: &mut Vec<LookupResult>,
        ) {
            let mut chunks = k.chunks_exact(B);
            for c in &mut chunks {
                out.extend_from_slice(&d.query_batch::<K, B, P>(c.try_into().unwrap()));
            }
            out.extend(chunks.remainder().iter().map(|k| d.query(k)));
        }
        match method {
            1 => return batch::<4, true>(d, kmers, out),
            2 => return batch::<8, true>(d, kmers, out),
            3 => return batch::<16, true>(d, kmers, out),
            4 => return batch::<32, true>(d, kmers, out),
            5 => return batch::<64, true>(d, kmers, out),
            6 => return batch::<16, false>(d, kmers, out),
            _ => (),
        }
    }
    let _ = method;
    out.extend(kmers.iter().map(|k| d.query(k)));
}
// Construction and workload generation are excluded from the timed query loops.
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut d = Dictionary::load(&args[1]).expect("load existing index");
    assert_eq!(d.k(), K);
    let nr: usize = args.get(3).map(|x| x.parse().unwrap()).unwrap_or(20000);
    let reps: usize = std::env::var("REPS")
        .ok()
        .map(|x| x.parse().unwrap())
        .unwrap_or(9);
    println!(
        "INDEX,{},{},{},{},{}",
        args[1],
        d.num_strings(),
        d.num_minimizers(),
        d.num_bits(),
        nr
    );
    let mut rng = 0x543971u64;
    let mut positive = Vec::new();
    for j in 0..nr {
        let id = next(&mut rng) % d.num_strings() as usize;
        let (begin, end) = (d.spss().offset_at(id), d.spss().offset_at(id + 1));
        let len = (end - begin).min(200) as usize;
        assert!(len >= K);
        let start = begin as usize + next(&mut rng) % (end as usize - begin as usize - len + 1);
        let text = d.spss().strings_data();
        let r: Vec<_> = (start..start + len)
            .map(|p| b"ACTG"[((text[p / 4] >> (2 * (p % 4))) & 3) as usize])
            .collect();
        positive.push(if j % 2 == 0 { r } else { rc(&r) });
    }
    let mut workloads = vec![("positive", positive.clone())];
    for rate in [1, 5, 100] {
        let reads = positive
            .iter()
            .map(|r| {
                r.iter()
                    .map(|&b| {
                        if next(&mut rng) % 100 < rate {
                            if rate == 100 {
                                b"ACGT"[next(&mut rng) % 4]
                            } else {
                                let i = b"ACGT".iter().position(|&x| x == b).unwrap();
                                b"ACGT"[(i + 1 + next(&mut rng) % 3) % 4]
                            }
                        } else {
                            b
                        }
                    })
                    .collect()
            })
            .collect();
        workloads.push((
            match rate {
                1 => "sub1",
                5 => "sub5",
                _ => "random",
            },
            reads,
        ));
    }
    if let Some(path) = args.get(4) {
        let mut reader = needletail::parse_fastx_file(path).unwrap();
        let mut reads = Vec::new();
        while reads.len() < nr {
            let Some(r) = reader.next() else { break };
            let r = r.unwrap();
            reads.push(r.seq().to_ascii_uppercase());
        }
        println!("READS,{path},{}", reads.len());
        workloads.push(("real", reads));
    }
    let mut baseline = std::io::BufWriter::new(std::fs::File::create(&args[2]).unwrap());
    let methods = if cfg!(feature = "experiments") { 7 } else { 1 };
    let names = [
        "scalar",
        "batch4",
        "batch8",
        "batch16",
        "batch32",
        "batch64",
        "batch16-no-prefetch",
    ];
    for (label, reads) in workloads {
        // Independent seeds: one valid window per read, repeated with shuffled
        // offsets to obtain enough timing work without adjacent-window locality.
        let mut seeds = Vec::new();
        for _ in 0..10 {
            for r in &reads {
                if r.len() >= K {
                    let p = next(&mut rng) % (r.len() - K + 1);
                    if let Ok(k) = Kmer::from_str(std::str::from_utf8(&r[p..p + K]).unwrap()) {
                        seeds.push(k);
                    }
                }
            }
        }
        let mut want = Vec::new();
        points(&d, &seeds, &mut want, 0);
        if label == "positive" {
            assert!(
                want.iter().all(|x| x.is_found()),
                "loaded index must find all sampled reference seeds"
            );
        }
        let mut want_stream = Vec::new();
        stream(&d, &reads, &mut want_stream);
        for xs in [&want, &want_stream] {
            baseline
                .write_all(&(xs.len() as u64).to_le_bytes())
                .unwrap();
            for &x in xs {
                for n in encode(x) {
                    baseline.write_all(&n.to_le_bytes()).unwrap();
                }
            }
        }
        println!(
            "CHECK_BASE,{label},seeds={},windows={},found_seeds={},found_windows={}",
            want.len(),
            want_stream.len(),
            want.iter().filter(|x| x.is_found()).count(),
            want_stream.iter().filter(|x| x.is_found()).count()
        );
        let mut out = Vec::with_capacity(want.len().max(want_stream.len()));
        let filters = if cfg!(feature = "experiments") {
            if std::env::var_os("FILTER_REVERSE").is_some() {
                vec![1, 0]
            } else {
                vec![0, 1]
            }
        } else {
            vec![0]
        };
        for filter in filters {
            #[cfg(feature = "experiments")]
            if filter == 1 {
                let t = Instant::now();
                let bytes = d.enable_bucket_fingerprints();
                println!(
                    "FILTER,{label},bytes={bytes},build_s={}",
                    t.elapsed().as_secs_f64()
                );
            }
            for method in 0..methods {
                points(&d, &seeds, &mut out, method);
                assert_eq!(out, want, "{label}: method {method} filter {filter}");
            }
            stream(&d, &reads, &mut out);
            assert_eq!(out, want_stream, "{label}: stream filter {filter}");
            println!("CHECK_EXACT,{label},filter={filter},all_fields=true");
            let order: Vec<usize> = if methods > 1 && std::env::var_os("QUICK").is_some() {
                vec![0, 4, methods]
            } else {
                (0..=methods).collect()
            };
            for rep in 0..reps {
                for j in 0..order.len() {
                    let method = order[(j + rep) % order.len()];
                    let wall = std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    let t = Instant::now();
                    if method == methods {
                        stream(&d, black_box(&reads), &mut out);
                    } else {
                        points(&d, black_box(&seeds), &mut out, method);
                    }
                    black_box(&out);
                    let ns = t.elapsed().as_nanos() as f64 / out.len() as f64;
                    let end = std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    println!(
                        "TIME,{label},{filter},{},{rep},{},{ns:.4},{wall:.6},{end:.6}",
                        if method == methods {
                            "stream"
                        } else {
                            names[method]
                        },
                        out.len()
                    );
                }
            }
            #[cfg(feature = "experiments")]
            d.disable_bucket_fingerprints();
        }
    }
}
