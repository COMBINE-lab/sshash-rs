//! Complete reference q-mer filters and sparse positive-hit replay.
//! Build/exhaustive no-false-negative checks are outside query timing. Each
//! process retains exactly one filter. CLI: INDEX N FASTX Q KEYS KIND [TRACE...].
//! KIND=extract creates a sorted unique key file; others benchmark full filters.
//! This is a research harness: legacy negative hints are deliberately not exposed.
#[path = "../common/qmer.rs"]
mod common;
use common::{K, Read};
use sshash_lib::{DeferredQuery, Dictionary, LookupResult};
use std::{
    hint::black_box,
    io::{BufReader, BufWriter, Read as IoRead, Write},
    time::{Instant, SystemTime, UNIX_EPOCH},
};
use value_traits::slices::SliceByValue;
use xorf::Filter;

/// SplitMix finalizer: deterministic mixing, not a cryptographic hash.
#[inline(always)]
fn mix(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
#[inline(always)]
fn mask(h: u64) -> u64 {
    (1 << (h & 63)) | (1 << ((h >> 6) & 63)) | (1 << ((h >> 12) & 63)) | (1 << ((h >> 18) & 63))
}
/// Canonical suffix from the already encoded pending k-mer. ACTG uses XOR 2
/// for complementation; canonicalization is identical during construction.
#[inline(always)]
fn suffix_bits(bits: u64, q: usize) -> u64 {
    let key = bits >> (2 * (K - q));
    key.min(sshash_lib::reverse_complement_mmer(key, q))
}
trait Membership {
    fn contains(&self, key: u64) -> bool;
    fn bytes(&self) -> usize;
}
impl Membership for xorf::BinaryFuse8 {
    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        Filter::contains(self, &key)
    }
    fn bytes(&self) -> usize {
        self.fingerprints.len()
    }
}
impl Membership for zor_filter::ZorFilter<u8> {
    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        self.contains(key)
    }
    fn bytes(&self) -> usize {
        self.main_filter().fingerprint_bytes()
            + self.remainder_filter().map_or(0, |r| r.fingerprint_bytes())
    }
}
/// Four hash-derived bits in one 64-bit register; power-of-two register count.
/// Requested bits/key may round upward substantially: report actual payload.
struct Bloom(Vec<u64>);
impl Bloom {
    fn build(keys: &[u64], bpk: usize) -> Self {
        let mut b = Self(vec![
            0;
            (keys.len() * bpk)
                .div_ceil(64)
                .max(1)
                .next_power_of_two()
        ]);
        let m = b.0.len() - 1;
        for &key in keys {
            let h = mix(key);
            b.0[(h >> 32) as usize & m] |= mask(h)
        }
        b
    }
}
impl Membership for Bloom {
    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        let h = mix(key);
        self.0[(h >> 32) as usize & (self.0.len() - 1)] & mask(h) == mask(h)
    }
    fn bytes(&self) -> usize {
        self.0.len() * 8
    }
}
/// Exact direct-address canonical-key bitmap. Unused noncanonical slots remain
/// allocated; no compressed-rank lookup is hidden in the measurement.
struct Bitmap(Vec<u64>);
impl Bitmap {
    fn new(q: usize) -> Self {
        Self(vec![0; (1usize << (2 * q)).div_ceil(64)])
    }
    fn insert(&mut self, key: u64) {
        self.0[key as usize / 64] |= 1 << (key % 64)
    }
    fn build(keys: &[u64], q: usize) -> Self {
        let mut b = Self::new(q);
        for &key in keys {
            b.insert(key)
        }
        b
    }
}
impl Membership for Bitmap {
    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        self.0[key as usize / 64] & (1 << (key % 64)) != 0
    }
    fn bytes(&self) -> usize {
        self.0.len() * 8
    }
}
/// Enumerate every canonical substring within SPSS strings, including ends.
/// A dense temporary bitmap deduplicates without a large hash-table overhead.
fn extract(d: &Dictionary, q: usize, path: &str) {
    if q == 0 {
        let now = Instant::now();
        let index = d.index_ref();
        let mut keys = Vec::with_capacity(index.num_buckets());
        for id in 0..index.num_buckets() {
            let (begin, end) = index.locate_bucket(id);
            if begin < end {
                let pos = index.offsets.index_value(begin);
                let key = d.spss().decode_mmer_at(pos, d.m());
                keys.push(key.min(sshash_lib::reverse_complement_mmer(key, d.m())));
            }
        }
        keys.sort_unstable();
        assert!(keys.windows(2).all(|w| w[0] < w[1]));
        let mut out = BufWriter::new(std::fs::File::create(path).unwrap());
        out.write_all(b"QKEYS001").unwrap();
        out.write_all(&0u64.to_le_bytes()).unwrap();
        out.write_all(&(keys.len() as u64).to_le_bytes()).unwrap();
        for key in &keys {
            out.write_all(&key.to_le_bytes()).unwrap();
        }
        out.flush().unwrap();
        println!(
            "EXTRACT,0,{},{:.6},{}",
            keys.len(),
            now.elapsed().as_secs_f64(),
            keys.len() * 8
        );
        return;
    }
    let start = Instant::now();
    let mut b = Bitmap::new(q);
    let text = d.spss().strings_data();
    for sid in 0..d.num_strings() as usize {
        let begin = d.spss().offset_at(sid) as usize;
        let end = d.spss().offset_at(sid + 1) as usize;
        let (mut fw, mut rc) = (0u64, 0u64);
        for (i, pos) in (begin..end).enumerate() {
            let base = ((text[pos / 4] >> (2 * (pos % 4))) & 3) as u64;
            fw = (fw >> 2) | (base << 62);
            rc = (rc << 2) | (base ^ 2);
            if i + 1 >= q {
                b.insert((fw >> (64 - 2 * q)).min(rc & ((1u64 << (2 * q)) - 1)));
            }
        }
    }
    let n: u64 = b.0.iter().map(|v| v.count_ones() as u64).sum();
    let mut out = BufWriter::new(std::fs::File::create(path).unwrap());
    out.write_all(b"QKEYS001").unwrap();
    out.write_all(&(q as u64).to_le_bytes()).unwrap();
    out.write_all(&n.to_le_bytes()).unwrap();
    for (i, &word) in b.0.iter().enumerate() {
        let mut w = word;
        while w != 0 {
            let key = (i * 64 + w.trailing_zeros() as usize) as u64;
            out.write_all(&key.to_le_bytes()).unwrap();
            w &= w - 1;
        }
    }
    out.flush().unwrap();
    println!(
        "EXTRACT,{q},{n},{:.6},{}",
        start.elapsed().as_secs_f64(),
        b.bytes()
    );
}
fn load_keys(path: &str, q: usize) -> Vec<u64> {
    let mut f = BufReader::new(std::fs::File::open(path).unwrap());
    let mut buf = [0u8; 8];
    f.read_exact(&mut buf).unwrap();
    assert_eq!(&buf, b"QKEYS001");
    f.read_exact(&mut buf).unwrap();
    assert_eq!(u64::from_le_bytes(buf), q as u64);
    f.read_exact(&mut buf).unwrap();
    let n = u64::from_le_bytes(buf) as usize;
    let mut keys = Vec::with_capacity(n);
    for _ in 0..n {
        f.read_exact(&mut buf).unwrap();
        keys.push(u64::from_le_bytes(buf));
    }
    assert!(keys.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(f.read(&mut buf).unwrap(), 0);
    keys
}

/// Existing seven-bit baseline reproduced in the harness: the same MPHF, mixer,
/// and bypass marker for heavy/empty buckets. Payload excludes the already
/// required MPHF. This times the pre-bucket decision, not a full seed lookup.
struct Tags<'a> {
    d: &'a Dictionary,
    tags: Vec<u8>,
}
impl<'a> Tags<'a> {
    fn build(d: &'a Dictionary) -> Self {
        let index = d.index_ref();
        let tags = (0..index.num_buckets())
            .map(|id| {
                let (begin, end) = index.locate_bucket(id);
                if begin == end || end - begin > 64 {
                    return 128;
                }
                let key = d
                    .spss()
                    .decode_mmer_at(index.offsets.index_value(begin), d.m());
                (mix(key.min(sshash_lib::reverse_complement_mmer(key, d.m()))) & 127) as u8
            })
            .collect();
        Self { d, tags }
    }
}
impl Membership for Tags<'_> {
    #[inline(always)]
    fn contains(&self, key: u64) -> bool {
        self.d.control_map_ref().lookup(key).is_some_and(|id| {
            let t = self.tags[id];
            t == 128 || t == (mix(key) & 127) as u8
        })
    }
    fn bytes(&self) -> usize {
        self.tags.len()
    }
}

/// Linux process RSS/HWM in KiB. HWM includes dictionary, keys and workloads;
/// report the pre-build RSS separately instead of calling it filter-only memory.
fn memory(field: &str) -> usize {
    std::fs::read_to_string("/proc/self/status")
        .unwrap()
        .lines()
        .find(|l| l.starts_with(field))
        .unwrap()
        .split_whitespace()
        .nth(1)
        .unwrap()
        .parse()
        .unwrap()
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
    filter: &impl Membership,
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
                        let key = suffix_bits(query.pending_seed().unwrap().0.as_u64(), q);
                        if !filter.contains(key) {
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
fn digest(
    d: &Dictionary,
    reads: &[Read],
    filter: &impl Membership,
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
        skipped(d, reads, filter, q, mode == "after-miss", |_, _, r| {
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

/// Classify and time the actual pending-seed suffixes. These are repeated-vector
/// throughput measurements, not dependent-query latency. Exact membership is
/// determined outside timing using the sorted full reference key set.
fn probe<F: Membership>(
    filter: &F,
    keys: &[u64],
    d: &Dictionary,
    reads: &[Read],
    q: usize,
    name: &str,
    kind: &str,
    reps: usize,
) {
    let mut query = DeferredQuery::<K>::new(d);
    let mut probes = Vec::new();
    for r in reads {
        query.reset();
        for (j, &p) in r.positions.iter().enumerate() {
            if j != 0 && p != r.positions[j - 1] + 1 {
                query.reset()
            }
            if query.prepare(&r.bases[p..p + K]).is_none() {
                let (kmer, minimizer) = query.pending_seed().unwrap();
                probes.push(if q == 0 {
                    minimizer.value
                } else {
                    suffix_bits(kmer.as_u64(), q)
                });
                query.complete_scalar();
            }
        }
    }
    let (mut absent, mut fp) = (0, 0);
    for &key in &probes {
        let yes = keys.binary_search(&key).is_ok();
        let hit = filter.contains(key);
        assert!(!yes || hit);
        if !yes {
            absent += 1;
            fp += usize::from(hit)
        }
    }
    println!(
        "FILTER_PROBES,{q},{kind},{name},{},{absent},{fp}",
        probes.len()
    );
    if probes.is_empty() {
        return;
    }
    let want = probes.iter().filter(|&&key| filter.contains(key)).count();
    let inner = 5_000_000_usize.div_ceil(probes.len());
    for rep in 0..reps {
        quiet();
        let start = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        let now = Instant::now();
        for _ in 0..inner {
            assert_eq!(
                black_box(
                    probes
                        .iter()
                        .filter(|&&key| filter.contains(black_box(key)))
                        .count()
                ),
                want
            );
        }
        let ns = now.elapsed().as_nanos() as f64 / (inner * probes.len()) as f64;
        let end = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        println!(
            "LAYOUT,{q},{kind}-q{q},{name},probe,{rep},{},{ns},{start},{end}",
            probes.len()
        );
    }
}
fn run<F: Membership>(
    filter: F,
    keys: Vec<u64>,
    d: &Dictionary,
    workloads: &[(String, Vec<Read>)],
    q: usize,
    kind: &str,
    build_seconds: f64,
    rss_before: usize,
) {
    println!(
        "FILTER_BUILD,{q},{kind},{},{},{build_seconds},{rss_before},{},{}",
        keys.len(),
        filter.bytes(),
        memory("VmRSS:"),
        memory("VmHWM:")
    );
    for &key in &keys {
        assert!(filter.contains(key), "false negative for {key}");
    }
    println!("EXACT,{q},{kind},all-inserted,{}", keys.len());
    let reps: usize = std::env::var("REPS").unwrap_or("7".into()).parse().unwrap();
    for (name, reads) in workloads {
        probe(&filter, &keys, d, reads, q, name, kind, reps);
        if q == 0 {
            continue;
        }
        let count: usize = reads.iter().map(|r| r.positions.len()).sum();
        if count == 0 {
            continue;
        }
        for mode in ["always", "after-miss"] {
            let mut seen = vec![0; reads.len()];
            let w = skipped(
                d,
                reads,
                &filter,
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
                "OPPORTUNITY,{q},{name},{kind}-{mode},{count},{},{},{},{},{}",
                w.seeds, w.probes, w.certificates, w.covered, w.local
            );
            println!("EXACT,{q},{kind},{name},{mode}");
        }
        let want = digest(d, reads, &filter, q, "scalar");
        for rep in 0..reps {
            for j in 0..3 {
                let mode = ["scalar", "always", "after-miss"][(j + rep) % 3];
                quiet();
                assert_eq!(digest(d, reads, &filter, q, mode), want);
                let start = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs_f64();
                let now = Instant::now();
                let inner = 5_000_000_usize.div_ceil(count).max(1);
                for _ in 0..inner {
                    assert_eq!(digest(d, black_box(reads), &filter, q, mode), want);
                }
                let ns = now.elapsed().as_nanos() as f64 / (count * inner) as f64;
                let end = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs_f64();
                println!("LAYOUT,{q},{kind}-q{q},{name},{mode},{rep},{count},{ns},{start},{end}");
            }
        }
    }
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    assert!(args.len() >= 7, "INDEX N FASTX Q KEYFILE KIND [TRACE...]");
    let mut d = Dictionary::load(&args[1]).unwrap();
    assert_eq!(d.k(), K);
    let q: usize = args[4].parse().unwrap();
    assert!(q <= 20);
    let kind = &args[6];
    assert!(q != 0 || kind != "bitmap", "bitmap requires q > 0");
    if kind == "extract" {
        extract(&d, q, &args[5]);
        return;
    }
    let keys = load_keys(&args[5], q);
    let mut workload_args = args[..4].to_vec();
    workload_args.extend_from_slice(&args[7..]);
    let workloads = common::workloads(&d, &workload_args);
    d.enable_bucket_fingerprints();
    let rss = memory("VmRSS:");
    let start = Instant::now();
    match kind.as_str() {
        "mphf7" => {
            assert_eq!(q, 0);
            let f = Tags::build(&d);
            let elapsed = start.elapsed().as_secs_f64();
            run(f, keys, &d, &workloads, q, kind, elapsed, rss)
        }
        "fuse8" => {
            let f = xorf::BinaryFuse8::try_from(keys.as_slice())
                .expect("binary fuse construction failed");
            let elapsed = start.elapsed().as_secs_f64();
            run(f, keys, &d, &workloads, q, kind, elapsed, rss)
        }
        "zor4" | "zor5" => {
            let config = zor_filter::FilterConfig {
                num_hashes: if kind == "zor4" { 4 } else { 5 },
                ..Default::default()
            };
            let out = zor_filter::ZorFilter::<u8>::build_with_config(&keys, &config)
                .expect("complete ZOR construction failed");
            let elapsed = start.elapsed().as_secs_f64();
            println!(
                "ZOR_REMAINDER,{q},{kind},{},{}",
                out.main_abandoned_keys.len(),
                out.remainder_total_slots.unwrap_or(0)
            );
            drop(out.main_abandoned_keys);
            run(out.filter, keys, &d, &workloads, q, kind, elapsed, rss)
        }
        "bloom8" | "bloom16" => {
            let f = Bloom::build(&keys, if kind == "bloom8" { 8 } else { 16 });
            let elapsed = start.elapsed().as_secs_f64();
            run(f, keys, &d, &workloads, q, kind, elapsed, rss)
        }
        "bitmap" => {
            let f = Bitmap::build(&keys, q);
            let elapsed = start.elapsed().as_secs_f64();
            run(f, keys, &d, &workloads, q, kind, elapsed, rss)
        }
        _ => panic!("unknown filter"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn suffix_encoding_and_complete_filters() {
        let seq = &b"ACTGGTCAACGTACGTGGGTACCGTACCTGAA"[..K];
        let k = sshash_lib::Kmer::<K>::from_ascii_unchecked(seq);
        for q in [12, 14, 16, 18, 20] {
            let mut bits = 0u64;
            for (i, &base) in seq[K - q..].iter().enumerate() {
                bits |= (b"ACTG".iter().position(|&v| v == base).unwrap() as u64) << (2 * i)
            }
            assert_eq!(
                suffix_bits(k.as_u64(), q),
                bits.min(sshash_lib::reverse_complement_mmer(bits, q))
            );
        }
        let keys: Vec<_> = (0..10000u64).map(|i| mix(i)).collect();
        for bpk in [8, 16] {
            let f = Bloom::build(&keys, bpk);
            assert!(keys.iter().all(|&k| f.contains(k)));
        }
        let f = xorf::BinaryFuse8::try_from(keys.as_slice()).unwrap();
        assert!(keys.iter().all(|&k| Membership::contains(&f, k)));
        for n in [4, 5] {
            let out = zor_filter::ZorFilter::<u8>::build_with_config(
                &keys,
                &zor_filter::FilterConfig {
                    num_hashes: n,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(keys.iter().all(|&k| out.filter.contains(k)));
        }
        let keys: Vec<_> = (0..10000).collect();
        let b = Bitmap::build(&keys, 7);
        assert!((0..16384).all(|k| b.contains(k) == (k < 10000)));
    }
}
