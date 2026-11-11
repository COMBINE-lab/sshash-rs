//! Shared deterministic workloads and exact scalar oracle for absence studies.
use sshash_lib::{Dictionary, LookupResult};
pub const K: usize = 31;
pub struct Read {
    pub bases: Vec<u8>,
    pub positions: Vec<usize>,
    pub expected: Vec<LookupResult>,
}
fn next(s: &mut u64) -> usize {
    *s = s
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*s >> 32) as usize
}
pub fn whole(bases: Vec<u8>) -> Read {
    let positions = (0..bases.len().saturating_sub(K - 1)).collect();
    Read {
        bases,
        positions,
        expected: Vec::new(),
    }
}
pub fn load_trace(path: &str) -> Vec<Read> {
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
pub fn scalar(d: &Dictionary, reads: &[Read], mut emit: impl FnMut(usize, usize, LookupResult)) {
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

pub fn workloads(d: &Dictionary, args: &[String]) -> Vec<(String, Vec<Read>)> {
    let n: usize = args[2].parse().unwrap();
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
        scalar(d, reads, |i, j, r| want[i][j] = r);
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

    workloads
}
