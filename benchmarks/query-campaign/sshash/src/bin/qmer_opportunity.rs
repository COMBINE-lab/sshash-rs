//! Exact, query-specific q-mer absence oracle and streaming skip opportunity.
//! Scan the complete SPSS, but retain membership only for queried suffix q-mers.
//! This is an oracle, NOT a deployable reference filter. A build-side Bloom
//! accelerator only avoids target-map probes; every actual match is exact.
//! Args match negative_pipeline: INDEX N REAL_FASTX [TRACE...].
#[path = "../common/qmer.rs"]
mod common;
use common::{K, Read};
use sshash_lib::{DeferredQuery, Dictionary, LookupResult};
use std::{
    collections::HashMap,
    hash::{BuildHasherDefault, Hasher},
    hint::black_box,
    time::{Instant, SystemTime, UNIX_EPOCH},
};
#[derive(Default)]
struct Mix(u64);
#[inline]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
impl Hasher for Mix {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = mix(x)
    }
    fn write(&mut self, _: &[u8]) {
        unreachable!("u64 keys only")
    }
}
type Map = HashMap<u64, usize, BuildHasherDefault<Mix>>;
#[inline]
fn mask(h: u64) -> u64 {
    (1 << (h & 63)) | (1 << ((h >> 6) & 63)) | (1 << ((h >> 12) & 63)) | (1 << ((h >> 18) & 63))
}
/// Canonical suffix key in SSHash's ACTG encoding, including lowercase input.
fn suffix(r: &Read, j: usize, q: usize) -> Option<u64> {
    let p = r.positions[j] + K - q;
    let mut key = 0;
    for (i, &b) in r.bases[p..p + q].iter().enumerate() {
        let v = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'T' | b't' => 2,
            b'G' | b'g' => 3,
            _ => return None,
        };
        key |= v << (2 * i);
    }
    Some(key.min(sshash_lib::reverse_complement_mmer(key, q)))
}
struct Oracle {
    q: usize,
    targets: Map,
    present: Vec<bool>,
    bloom: Vec<u64>,
    hll: Vec<u8>,
}
impl Oracle {
    fn new(q: usize, workloads: &[(String, Vec<Read>)]) -> Self {
        let mut targets = Map::default();
        for (_, reads) in workloads {
            for r in reads {
                for j in 0..r.positions.len() {
                    if let Some(key) = suffix(r, j, q) {
                        let n = targets.len();
                        targets.entry(key).or_insert(n);
                    }
                }
            }
        }
        let present = vec![false; targets.len()];
        let mut bloom = vec![0u64; (targets.len() * 16).div_ceil(64).max(1).next_power_of_two()];
        let bmask = bloom.len() - 1;
        for &key in targets.keys() {
            let h = mix(key);
            bloom[((h >> 32) as usize) & bmask] |= mask(h);
        }
        Self {
            q,
            targets,
            present,
            bloom,
            hll: vec![0; 16384],
        }
    }
    #[inline]
    fn observe(&mut self, key: u64) {
        let h = mix(key);
        let bin = (h & 16383) as usize;
        let rank = ((h >> 14).leading_zeros() - 14 + 1) as u8;
        self.hll[bin] = self.hll[bin].max(rank);
        if self.bloom[((h >> 32) as usize) & (self.bloom.len() - 1)] & mask(h) == mask(h) {
            if let Some(&id) = self.targets.get(&key) {
                self.present[id] = true;
            }
        }
    }
    fn answers(&self, reads: &[Read]) -> Vec<Vec<bool>> {
        reads
            .iter()
            .map(|r| {
                (0..r.positions.len())
                    .map(|j| {
                        suffix(r, j, self.q).is_some_and(|key| {
                            self.present
                                [*self.targets.get(&key).expect("unregistered oracle query")]
                        })
                    })
                    .collect()
            })
            .collect()
    }
    fn distinct_estimate(&self) -> f64 {
        let m = self.hll.len() as f64;
        let raw = (0.7213 / (1. + 1.079 / m)) * m * m
            / self
                .hll
                .iter()
                .map(|&r| 2f64.powi(-(r as i32)))
                .sum::<f64>();
        let zero = self.hll.iter().filter(|&&r| r == 0).count();
        if raw < 2.5 * m && zero > 0 {
            m * (m / zero as f64).ln()
        } else {
            raw
        }
    }
}
/// One complete scan proves absence for the registered query set. Canonical
/// enumeration covers both strands, including substrings near SPSS boundaries.
/// HLL cardinalities are estimates (~0.8% nominal relative standard error), used
/// only for planning filter memory; exact membership never depends on HLL.
fn build_oracles(d: &Dictionary, workloads: &[(String, Vec<Read>)]) -> Vec<Oracle> {
    let qs: Vec<usize> = std::env::var("QS")
        .unwrap_or("12,14,16,18,20".into())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    assert!(qs.iter().all(|&q| q > 0 && q <= K));
    let mut oracles: Vec<_> = qs.iter().map(|&q| Oracle::new(q, workloads)).collect();
    let text = d.spss().strings_data();
    let start = Instant::now();
    let mut bases = 0usize;
    for sid in 0..d.num_strings() as usize {
        let begin = d.spss().offset_at(sid) as usize;
        let end = d.spss().offset_at(sid + 1) as usize;
        let (mut fw, mut rc) = (0u64, 0u64);
        for (i, pos) in (begin..end).enumerate() {
            let base = ((text[pos / 4] >> (2 * (pos % 4))) & 3) as u64;
            fw = (fw >> 2) | (base << 62);
            rc = (rc << 2) | (base ^ 2);
            for oracle in &mut oracles {
                let q = oracle.q;
                if i + 1 >= q {
                    let key = (fw >> (64 - 2 * q)).min(rc & ((1u64 << (2 * q)) - 1));
                    oracle.observe(key);
                }
            }
        }
        bases += end - begin;
    }
    println!("ORACLE_BUILD,{bases},{:.6}", start.elapsed().as_secs_f64());
    for o in &oracles {
        println!(
            "ORACLE_SET,{},{},{},{:.0},{}",
            o.q,
            o.targets.len(),
            o.present.iter().filter(|&&v| v).count(),
            o.distinct_estimate(),
            o.bloom.len() * 8
        );
    }
    oracles
}
#[derive(Default, Debug, Clone, Copy)]
struct Work {
    seeds: u64,
    probes: u64,
    certificates: u64,
    covered: u64,
    local: u64,
}
/// New sparse-hit contract: exact positive positions/fields, no negative hints.
/// A rejected suffix q-mer certifies starts [pos,pos+K-q]. Discard stream state
/// when skipping, so no incorrect minimizer memo survives; count the resulting
/// extra seeds as well as saved seeds. Backward requests outside the remembered
/// interval fall back normally. `after_miss` protects initial/positive seeds.
fn skipped(
    d: &Dictionary,
    reads: &[Read],
    answers: &[Vec<bool>],
    q: usize,
    after_miss: bool,
    mut emit: impl FnMut(usize, usize, Option<LookupResult>),
) -> Work {
    let mut work = Work::default();
    let mut query = DeferredQuery::<K>::new(d);
    for (i, r) in reads.iter().enumerate() {
        query.reset();
        let mut last = None;
        let mut cert = None;
        let mut negative = false;
        for (j, &pos) in r.positions.iter().enumerate() {
            if cert.is_some_and(|(lo, hi)| pos >= lo && pos <= hi) {
                if last.is_some() {
                    query.reset();
                    last = None;
                }
                negative = true;
                work.covered += 1;
                emit(i, j, None);
                continue;
            }
            if last != pos.checked_sub(1) {
                query.reset();
            }
            let result = match query.prepare(&r.bases[pos..pos + K]) {
                Some(result) => {
                    work.local += 1;
                    result
                }
                None => {
                    if !after_miss || negative {
                        work.probes += 1;
                        if !answers[i][j] {
                            work.certificates += 1;
                            work.covered += 1;
                            cert = Some((pos, pos + K - q));
                            query.reset();
                            last = None;
                            negative = true;
                            emit(i, j, None);
                            continue;
                        }
                    }
                    work.seeds += 1;
                    query.complete_scalar()
                }
            };
            last = Some(pos);
            negative = !result.is_found();
            emit(i, j, result.is_found().then_some(result));
        }
    }
    work
}
fn baseline_seeds(d: &Dictionary, reads: &[Read]) -> u64 {
    let mut n = 0;
    let mut q = DeferredQuery::<K>::new(d);
    for r in reads {
        q.reset();
        for (j, &p) in r.positions.iter().enumerate() {
            if j != 0 && p != r.positions[j - 1] + 1 {
                q.reset();
            }
            if q.prepare(&r.bases[p..p + K]).is_none() {
                n += 1;
                q.complete_scalar();
            }
        }
    }
    n
}
fn digest(
    d: &Dictionary,
    reads: &[Read],
    answers: &[Vec<bool>],
    q: usize,
    mode: &str,
) -> (u64, u64) {
    let (mut found, mut sum) = (0u64, 0u64);
    let mut add = |r: LookupResult| {
        if r.is_found() {
            found += 1;
            sum = sum.wrapping_add(r.kmer_id);
        }
    };
    if mode == "scalar" {
        common::scalar(d, reads, |_, _, r| add(r));
    } else {
        skipped(d, reads, answers, q, mode == "after-miss", |_, _, r| {
            if let Some(r) = r {
                add(r)
            }
        });
    }
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
    let workloads = common::workloads(&d, &args);
    let oracles = build_oracles(&d, &workloads);
    d.enable_bucket_fingerprints();
    let reps: usize = std::env::var("REPS").unwrap_or("7".into()).parse().unwrap();
    for (name, reads) in &workloads {
        let baseline = baseline_seeds(&d, reads);
        let count: usize = reads.iter().map(|r| r.positions.len()).sum();
        println!("BASELINE_SEEDS,{name},{count},{baseline}");
    }
    for oracle in &oracles {
        let q = oracle.q;
        for (name, reads) in &workloads {
            let answers = oracle.answers(reads);
            let count: usize = reads.iter().map(|r| r.positions.len()).sum();
            if count == 0 {
                continue;
            }
            for mode in ["always", "after-miss"] {
                let mut seen = vec![0; reads.len()];
                let w = skipped(
                    &d,
                    reads,
                    &answers,
                    q,
                    mode == "after-miss",
                    |i, j, result| {
                        assert_eq!(j, seen[i]);
                        seen[i] += 1;
                        let expected = reads[i].expected[j];
                        assert_eq!(result.is_some(), expected.is_found(), "{name}/q{q}/{i}/{j}");
                        if let Some(r) = result {
                            assert_eq!(r, expected);
                        }
                    },
                );
                assert!(seen.iter().zip(reads).all(|(n, r)| *n == r.positions.len()));
                println!(
                    "OPPORTUNITY,{q},{name},{mode},{count},{},{},{},{},{}",
                    w.seeds, w.probes, w.certificates, w.covered, w.local
                );
                println!("EXACT,{q},oracle,{name},{mode}");
            }
            let want = digest(&d, reads, &answers, q, "scalar");
            for rep in 0..reps {
                for j in 0..3 {
                    let mode = ["scalar", "always", "after-miss"][(j + rep) % 3];
                    quiet();
                    assert_eq!(digest(&d, reads, &answers, q, mode), want);
                    let start = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    let now = Instant::now();
                    let inner = 5_000_000_usize.div_ceil(count).max(1);
                    for _ in 0..inner {
                        assert_eq!(digest(&d, black_box(reads), &answers, q, mode), want);
                    }
                    let ns = now.elapsed().as_nanos() as f64 / (count * inner) as f64;
                    let end = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .unwrap()
                        .as_secs_f64();
                    println!(
                        "LAYOUT,{q},oracle-q{q},{name},{mode},{rep},{count},{ns},{start},{end}"
                    );
                }
            }
        }
    }
}
