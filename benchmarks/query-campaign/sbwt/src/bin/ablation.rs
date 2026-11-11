//! Factorial ablation: scalar/queued searches, each without/with look-ahead skipping.
//! Same index, rank layout, prefix table, reads, and result consumer for all four.
//! COUNT=true work counters run outside timing; timed APIs use COUNT=false.
//! Independent positional differential checks and a synthetic throughput experiment.
//! Uses the feature-gated SSHash candidate and the pinned, unmodified SBWT source.
//! Timings are NOT the paper's datasets or mapping timings. See the campaign README
//! for argument defaults, representation differences, and timing exclusions.
use sbwt::line_rank::{LineRankInline, LineRankMatrix};
use sbwt::restart_search::{
    RestartStats, restart_kmers, restart_kmers_fskip, restart_kmers_queued,
};
use sbwt::restart_search::{
    SuffixGroupMarks, both_orientations_kmers_queued, restart_kmers_queued_fskip,
};
use sbwt::subset_fixed_correction_sets::SubsetFixedCorrectionSets;
use sbwt::{BitPackedKmerSortingMem, SbwtIndex};
use sbwt::{SubsetMatrix, SubsetSeq};
use sshash_lib::Dictionary;
use sshash_lib::builder::{config::BuildConfiguration, dictionary_builder::DictionaryBuilder};
use std::{hint::black_box, time::Instant};

const K: usize = 31;
type Lines = LineRankMatrix<LineRankInline>;
/// Fixed workload generator; intentionally identical across comparison runs.
fn next(s: &mut u64) -> usize {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*s >> 33) as usize
}
/// Generate an ACGT-only sequence from the deterministic workload RNG.
fn dna(s: &mut u64, n: usize) -> Vec<u8> {
    (0..n).map(|_| b"ACGT"[next(s) % 4]).collect()
}
/// Reverse complement while preserving ambiguous bytes for correctness cases.
fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => *b,
        })
        .collect()
}
/// Convert the same subset matrix to cache-line-local rank storage.
fn lines(index: SbwtIndex<SubsetMatrix>) -> SbwtIndex<Lines> {
    let (m, n, k, c, p) = index.into_parts();
    SbwtIndex::from_parts(Lines::from_matrix_rank_only(m), n, k, c, p)
}
/// Build the compressed FBCS variant with the same k-mers and prefix length 8.
fn fcs(index: &SbwtIndex<SubsetMatrix>) -> SbwtIndex<SubsetFixedCorrectionSets> {
    let rows = (0..4)
        .map(|c| {
            (0..index.n_sets())
                .map(|i| index.sbwt().set_contains(i, c))
                .collect::<bitvec::vec::BitVec<u64>>()
        })
        .collect();
    SbwtIndex::from_subset_seq(
        SubsetFixedCorrectionSets::new_from_bit_vectors(rows),
        index.n_kmers(),
        K,
        8,
    )
}
/// Independent per-window scalar SBWT oracle, using its own colex-rank namespace.
fn expected<S: SubsetSeq>(
    idx: &SbwtIndex<S>,
    reads: &[Vec<u8>],
    both: bool,
) -> Vec<Vec<Option<usize>>> {
    reads
        .iter()
        .map(|r| {
            r.windows(K)
                .map(|w| {
                    idx.search(w)
                        .or_else(|| if both { idx.search(&rc(w)) } else { None })
                        .map(|x| x.start)
                })
                .collect()
        })
        .collect()
}
/// Validate all callback positions/ranks, including omissions and duplicates.
fn check<S: SubsetSeq>(
    name: &str,
    idx: &SbwtIndex<S>,
    marks: &SuffixGroupMarks,
    reads: &[Vec<u8>],
    both: bool,
) {
    let want = expected(idx, reads, both);
    let mut got: Vec<Vec<Option<usize>>> = want.iter().map(|r| vec![None; r.len()]).collect();
    let rcs: Vec<_> = reads.iter().map(|r| rc(r)).collect();
    let mut out = |r: usize, w: usize, x: usize| {
        assert!(got[r][w].is_none(), "{name}: duplicate callback at {r}/{w}");
        got[r][w] = Some(x);
    };
    if name.starts_with("fcs") {
        both_orientations_kmers_queued::<_, 32, 32>(idx, marks, reads, &rcs, &mut out);
    } else if both {
        both_orientations_kmers_queued::<_, 64, 32>(idx, marks, reads, &rcs, &mut out);
    } else {
        restart_kmers_queued_fskip::<_, _, 64, 32>(idx, marks, reads, &mut out);
    }
    for (r, (a, b)) in got.iter().zip(&want).enumerate() {
        assert_eq!(a, b, "{name}: read {r}");
    }
    println!(
        "CHECK {name}: {} reads, {} windows, exact positional colex ranks PASS",
        reads.len(),
        want.iter().map(Vec::len).sum::<usize>()
    );
}
/// Warm once, time once; the caller rotates methods over nine repetitions.
fn measure(name: &str, workload: &str, windows: usize, mut f: impl FnMut() -> (usize, usize)) {
    let want = black_box(f()); // untimed warmup
    let mut ns = Vec::new();
    for _ in 0..1 {
        let start = Instant::now();
        let got = black_box(f());
        ns.push(start.elapsed().as_nanos() as f64 / windows as f64);
        assert_eq!(got, want);
    }
    ns.sort_by(f64::total_cmp);
    println!(
        "TIME,{workload},{name},{windows},{},{:.3},{:.3}",
        want.0, ns[0], ns[0]
    );
}
/// Consume found IDs without full-result materialization, resetting per read.
fn stream(dict: &Dictionary, reads: &[Vec<u8>]) -> (usize, usize) {
    let mut q = dict.create_streaming_query::<K>();
    let (mut found, mut sum) = (0, 0usize);
    for r in reads {
        q.reset();
        for w in r.windows(K) {
            let x = q.lookup(w);
            if x.is_found() {
                found += 1;
                sum = sum.wrapping_add(x.kmer_id as usize);
            }
        }
    }
    (found, sum)
}
/// Exclude reverse-complement preparation; consume sparse callbacks in the timer.
fn sbwt_time<S: SubsetSeq>(
    name: &str,
    workload: &str,
    idx: &SbwtIndex<S>,
    marks: &SuffixGroupMarks,
    reads: &[Vec<u8>],
    both: bool,
) {
    let windows = reads.iter().map(|r| r.len().saturating_sub(K - 1)).sum();
    let rcs: Vec<_> = reads.iter().map(|r| rc(r)).collect(); // RC conversion outside timing, as paper harness
    measure(name, workload, windows, || {
        let (mut found, mut sum) = (0, 0usize);
        let mut out = |_: usize, _: usize, x: usize| {
            found += 1;
            sum = sum.wrapping_add(x);
        };
        if name.starts_with("fcs") {
            both_orientations_kmers_queued::<_, 32, 32>(idx, marks, reads, &rcs, &mut out);
        } else if both {
            both_orientations_kmers_queued::<_, 64, 32>(idx, marks, reads, &rcs, &mut out);
        } else {
            restart_kmers_queued_fskip::<_, _, 64, 32>(idx, marks, reads, &mut out);
        }
        (found, sum)
    });
}
/// Emit sparse positives from the three additional two-orientation variants.
fn ablation_run(
    idx: &SbwtIndex<Lines>,
    marks: &SuffixGroupMarks,
    reads: &[Vec<u8>],
    method: usize,
    mut out: impl FnMut(usize, usize, usize),
) {
    let mut stats = RestartStats::default();
    match method {
        5 => {
            for (r, read) in reads.iter().enumerate() {
                restart_kmers::<_, false>(idx, marks, read, &mut stats, |w, x| {
                    if let Some(x) = x {
                        out(r, w, x);
                    }
                });
            }
        }
        6 => {
            for (r, read) in reads.iter().enumerate() {
                restart_kmers_fskip::<_, false>(idx, marks, read, &mut stats, |w, x| out(r, w, x));
            }
        }
        7 => restart_kmers_queued::<_, _, 64, 32>(idx, marks, reads, |r, w, x| {
            if let Some(x) = x {
                out(r, w, x);
            }
        }),
        _ => unreachable!(),
    }
}

// Build and validate all representations before entering the timing loop.
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let n: usize = args.get(1).map(|x| x.parse().unwrap()).unwrap_or(200_000);
    let nr: usize = args.get(2).map(|x| x.parse().unwrap()).unwrap_or(2000);
    let reorient = args.get(3).is_some_and(|x| x == "reorient");
    assert!(n > 1000 && nr > 0);
    let mut rng = 0x543971u64;
    let genome = dna(&mut rng, n);
    // Partition into 300-kmer strings, overlapping by k-1 bases. Alternate stored
    // orientations to exercise reverse runs and transitions without duplicate kmers.
    let seqs: Vec<Vec<u8>> = (0..n - K + 1)
        .step_by(300)
        .enumerate()
        .map(|(j, i)| {
            let s = genome[i..(i + 300 + K - 1).min(n)].to_vec();
            if j % 2 == 0 { s } else { rc(&s) }
        })
        .collect();
    let mut cfg = BuildConfiguration::new(K, 19).unwrap();
    cfg.num_threads = 4;
    let dict = DictionaryBuilder::new(cfg.clone())
        .unwrap()
        .build_from_sequences(
            seqs.iter()
                .map(|s| String::from_utf8(s.clone()).unwrap())
                .collect(),
        )
        .unwrap();
    let mut filtered = DictionaryBuilder::new(cfg)
        .unwrap()
        .build_from_sequences(
            seqs.iter()
                .map(|s| String::from_utf8(s.clone()).unwrap())
                .collect(),
        )
        .unwrap();
    let filter_bytes = filtered.enable_bucket_fingerprints();
    println!("FILTER_BYTES,{filter_bytes}");
    let (two, _, marks2) = BitPackedKmerSortingMem::new_from_vecs(&seqs, K)
        .n_threads(4)
        .add_rev_comp(true)
        .build_marks(true)
        .run_with_marks();
    let oriented = if reorient {
        let orientations = sbwt::optimize_unitig_orientation(sbwt::VecSeqStream::new(&seqs), K);
        seqs.iter()
            .zip(orientations)
            .map(|(s, o)| match o {
                sbwt::Orientation::Forward => s.clone(),
                sbwt::Orientation::Reverse => rc(s),
            })
            .collect()
    } else {
        seqs.clone()
    };
    let (one, _, marks1) = BitPackedKmerSortingMem::new_from_vecs(&oriented, K)
        .n_threads(4)
        .add_rev_comp(false)
        .build_marks(true)
        .run_with_marks();
    assert_eq!(
        one.n_kmers(),
        n - K + 1,
        "synthetic reference has duplicate oriented kmers"
    );
    assert_eq!(
        two.n_kmers(),
        2 * (n - K + 1),
        "synthetic reference has canonical duplicates"
    );
    let compact = fcs(&one);
    let two = lines(two);
    let one = lines(one);
    let marks1 = marks1.unwrap();
    let marks2 = marks2.unwrap();
    println!(
        "CONFIG bases={n}, reads_per_workload={nr}, k={K}, m=19, seed=0x543971, string_kmers=300, alternating_orientation=true, sbwt_reorient={reorient}"
    );
    println!(
        "SIZE sshash_reported_bits={} two_rank_bytes={} one_rank_bytes={} fcs_bytes={} two_marks_bits={} one_marks_bits={} (SBWT prefix table, 1048576 bytes/index, excluded)",
        dict.num_bits(),
        two.sbwt().line_rank_size_in_bytes(),
        one.sbwt().line_rank_size_in_bytes(),
        compact.sbwt().size_in_bytes(),
        marks2.len(),
        marks1.len()
    );
    let mut checks = Vec::new();
    let mut workloads = Vec::new();
    for rate in [0, 1, 5, 100] {
        let reads: Vec<_> = (0..nr)
            .map(|j| {
                let pos = next(&mut rng) % (n - 200);
                let mut r = if rate == 100 {
                    dna(&mut rng, 200)
                } else {
                    genome[pos..pos + 200].to_vec()
                };
                if rate != 100 {
                    for b in &mut r {
                        if next(&mut rng) % 100 < rate {
                            let x = b"ACGT".iter().position(|x| x == b).unwrap();
                            *b = b"ACGT"[(x + 1 + next(&mut rng) % 3) % 4];
                        }
                    }
                }
                if j % 2 == 0 { r } else { rc(&r) }
            })
            .collect();
        checks.extend(reads.iter().cloned());
        workloads.push((format!("sub{rate}"), reads));
    }
    for len in [0, 1, 30, 31, 32, 63, 64, 65, 199, 200, 257] {
        checks.push(genome[..len].to_vec());
        for at in [0, len / 2, len.saturating_sub(1)] {
            let mut r = genome[..len].to_vec();
            if len > 0 {
                r[at] = b'N';
            }
            checks.push(r);
        }
    }
    checks.push(genome[..200].to_ascii_lowercase());
    // Positional validation, not just a hit count or commutative rank checksum.
    check("two-fskip", &two, &marks2, &checks, false);
    check("one-runs", &one, &marks1, &checks, true);
    check("fcs-runs", &compact, &marks1, &checks, true);
    for r in &checks {
        let mut a = dict.create_streaming_query::<K>();
        let mut b = filtered.create_streaming_query::<K>();
        for w in r.windows(K) {
            assert_eq!(
                a.lookup(w),
                b.lookup(w),
                "filter full-field streaming parity"
            );
        }
    }
    println!("CHECK filter full-field streaming parity PASS");
    let want = expected(&two, &checks, false);
    let mut offset_differences = 0;
    let mut first = None;
    for (r, e) in checks.iter().zip(&want) {
        let mut q = dict.create_streaming_query::<K>();
        for (w, expected) in r.windows(K).zip(e) {
            let a = q.lookup(w);
            let p = dict.query_from_str::<K>(std::str::from_utf8(w).unwrap());
            assert_eq!(a.is_found(), expected.is_some());
            assert_eq!(a.kmer_id, p.kmer_id);
            if a.is_found() {
                assert_eq!(
                    (
                        a.string_id,
                        a.kmer_id_in_string,
                        a.kmer_orientation,
                        a.string_begin,
                        a.string_end
                    ),
                    (
                        p.string_id,
                        p.kmer_id_in_string,
                        p.kmer_orientation,
                        p.string_begin,
                        p.string_end
                    )
                );
                if a.kmer_offset != p.kmer_offset {
                    offset_differences += 1;
                    first.get_or_insert((a.kmer_offset, p.kmer_offset));
                }
            }
        }
    }
    println!(
        "CHECK sshash: positional membership, IDs, string coordinates, orientation PASS; kmer_offset discrepancies={offset_differences}, first={first:?}"
    );
    for method in 5..8 {
        let mut got: Vec<_> = want.iter().map(|x| vec![None; x.len()]).collect();
        ablation_run(&two, &marks2, &checks, method, |r, w, x| {
            assert!(got[r][w].is_none());
            got[r][w] = Some(x);
        });
        assert_eq!(got, want, "positional ablation method {method}");
        println!("CHECK ablation {method}: all positions and ranks PASS");
    }
    for (label, reads) in &workloads {
        for skip in [false, true] {
            let mut stats = RestartStats::default();
            for read in reads {
                if skip {
                    restart_kmers_fskip::<_, true>(&two, &marks2, read, &mut stats, |_, _| {});
                } else {
                    restart_kmers::<_, true>(&two, &marks2, read, &mut stats, |_, _| {});
                }
            }
            println!(
                "WORK,{label},{skip},{},{},{},{},{},{}",
                stats.n_chars,
                stats.n_continuations,
                stats.n_continuations_found,
                stats.n_restarts,
                stats.n_table_hits,
                stats.n_restart_extensions
            );
        }
    }
    // Full exact vectors are checked above, checksums only prevent timing-loop elimination.
    println!("TIME_HEADER,workload,method,windows,found,best_ns_per_window,median_ns_per_window");
    for rep in 0..9 {
        for (label, reads) in &workloads {
            let windows = reads.iter().map(|r| r.len() - K + 1).sum();
            for j in 0..8 {
                let method = (j + rep) % 8;
                println!(
                    "TRIAL,{rep},{label},{method},{}",
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64()
                );
                match method {
                    0 => measure("sshash-stream", label, windows, || stream(&dict, reads)),
                    1 => measure("sshash-filter-stream", label, windows, || {
                        stream(&filtered, reads)
                    }),
                    2 => sbwt_time("two-fskip-64:32", label, &two, &marks2, reads, false),
                    3 => sbwt_time("one-runs-64:32", label, &one, &marks1, reads, true),
                    4 => sbwt_time("fcs-runs-32:32", label, &compact, &marks1, reads, true),
                    5..=7 => {
                        let name = ["two-scalar-noskip", "two-scalar-fskip", "two-queued-noskip"]
                            [method - 5];
                        measure(name, label, windows, || {
                            let (mut n, mut sum) = (0, 0usize);
                            ablation_run(&two, &marks2, reads, method, |_, _, x| {
                                n += 1;
                                sum = sum.wrapping_add(x);
                            });
                            (n, sum)
                        });
                    }
                    _ => unreachable!(),
                }
            }
        }
    }
}
