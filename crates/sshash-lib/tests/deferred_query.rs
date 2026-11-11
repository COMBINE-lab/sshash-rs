#![cfg(feature = "experimental-query")]
//! Differential state-machine tests: every field and per-read order must match
//! the corrected scalar stream, including resets, invalid bases and heavy hints.
use sshash_lib::{BuildConfiguration, DeferredQuery, DictionaryBuilder, Kmer, KmerBits};
fn dna(n: usize, state: &mut u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[(*state >> 32) as usize % 4]
        })
        .collect()
}
fn check<const K: usize>(m: usize)
where
    Kmer<K>: KmerBits,
{
    let mut rng = 77;
    let texts: Vec<_> = (0..12)
        .map(|_| String::from_utf8(dna(400, &mut rng)).unwrap())
        .collect();
    let mut dict = DictionaryBuilder::new(BuildConfiguration::new(K, m).unwrap())
        .unwrap()
        .build_from_sequences(texts.clone())
        .unwrap();
    let mut reads: Vec<_> = texts.iter().map(|s| s.as_bytes().to_vec()).collect();
    for r in reads.clone() {
        reads.push(
            r.iter()
                .rev()
                .map(|b| match b {
                    b'A' => b'T',
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => b'A',
                })
                .collect(),
        );
    }
    reads.push(dna(497, &mut rng));
    let mut invalid = reads[0].clone();
    invalid[50] = b'N';
    reads.push(invalid);
    reads.push(reads[1].to_ascii_lowercase());
    for filtered in [false, true] {
        if filtered {
            dict.enable_bucket_fingerprints();
        }
        let expected: Vec<Vec<_>> = reads
            .iter()
            .map(|r| {
                let mut q = dict.create_streaming_query::<K>();
                r.windows(K)
                    .enumerate()
                    .map(|(i, w)| {
                        if i == 99 {
                            q.reset();
                        }
                        q.lookup(w)
                    })
                    .collect()
            })
            .collect();
        for (id, read) in reads.iter().enumerate() {
            let mut direct = DeferredQuery::<K>::new(&dict);
            for (pos, window) in read.windows(K).enumerate() {
                if pos == 99 {
                    direct.reset();
                }
                assert_eq!(direct.lookup_filtered(window, |_| true), expected[id][pos]);
            }
        }
        for batched in [false, true] {
            let mut lanes: [DeferredQuery<K>; 4] =
                std::array::from_fn(|_| DeferredQuery::new(&dict));
            let mut ids = [usize::MAX; 4];
            let mut positions = [0; 4];
            let mut next = 0;
            let mut actual = vec![Vec::new(); reads.len()];
            loop {
                for i in 0..4 {
                    loop {
                        if ids[i] == usize::MAX {
                            if next == reads.len() {
                                break;
                            }
                            ids[i] = next;
                            next += 1;
                            positions[i] = 0;
                            lanes[i].reset();
                        }
                        let id = ids[i];
                        let pos = positions[i];
                        if pos + K > reads[id].len() {
                            ids[i] = usize::MAX;
                            continue;
                        }
                        if pos == 99 {
                            lanes[i].reset();
                        }
                        match lanes[i].prepare(&reads[id][pos..pos + K]) {
                            Some(r) => {
                                actual[id].push(r);
                                positions[i] += 1;
                            }
                            None => break,
                        }
                    }
                }
                if !lanes.iter().any(|q| q.is_pending()) {
                    break;
                }
                let pending = std::array::from_fn::<_, 4, _>(|i| lanes[i].is_pending());
                if batched {
                    dict.complete_seed_batch::<K, 4, true>(&mut lanes);
                } else {
                    for q in lanes.iter_mut().filter(|q| q.is_pending()) {
                        q.complete_scalar();
                    }
                }
                for i in 0..4 {
                    if pending[i] {
                        actual[ids[i]].push(lanes[i].result());
                        positions[i] += 1;
                    }
                }
            }
            assert_eq!(
                actual, expected,
                "filtered={filtered},batched={batched},k={K},m={m}"
            );
        }
    }
}
#[test]
fn complete_fields_memos_strands_and_resets() {
    check::<21>(3);
    check::<31>(18);
    check::<63>(31);
}

#[test]
fn empty_dictionary_partial_batches_and_cancelled_seed() {
    let d = sshash_lib::Dictionary::new(
        sshash_lib::spectrum_preserving_string_set::SpectrumPreservingStringSet::new(3, 2),
        sshash_lib::minimizers_control_map::MinimizersControlMap::empty(),
        sshash_lib::sparse_and_skew_index::SparseAndSkewIndex::new(),
        3,
        2,
        1,
    );
    let mut lanes: [DeferredQuery<3>; 4] = std::array::from_fn(|_| DeferredQuery::new(&d));
    assert!(lanes[0].prepare(b"ACG").is_none());
    assert!(lanes[1].prepare(b"TGC").is_none());
    lanes[1].reset(); // cancelled lanes are never completed or emitted
    assert!(!lanes[2].prepare(b"A").unwrap().is_found());
    d.complete_seed_batch::<3, 4, false>(&mut lanes);
    assert_eq!(
        lanes[0].result(),
        d.query(&Kmer::<3>::from_str("ACG").unwrap())
    );
    let mut empty: [DeferredQuery<3>; 0] = [];
    d.complete_seed_batch::<3, 0, true>(&mut empty);
}

/// An exact selected-minimizer set is a safe oracle for the independent filter.
/// It may strengthen negative hints, so compare membership and all positive
/// fields while preserving rolling state through negative runs and later hits.
#[test]
fn independent_minimizer_filter_preserves_positive_contract() {
    use std::collections::HashSet;
    let mut rng = 987;
    let texts: Vec<_> = (0..8).map(|_| dna(701, &mut rng)).collect();
    let strings: Vec<_> = texts
        .iter()
        .map(|s| String::from_utf8(s.clone()).unwrap())
        .collect();
    for m in [3, 18] {
        let d = DictionaryBuilder::new(BuildConfiguration::new(31, m).unwrap())
            .unwrap()
            .build_from_sequences(strings.clone())
            .unwrap();
        let selected: HashSet<_> = texts
            .iter()
            .flat_map(|s| s.windows(31))
            .map(|w| {
                d.debug_extract_minimizer(&Kmer::<31>::from_ascii_unchecked(w))
                    .value
            })
            .collect();
        let mut reads = texts.clone();
        reads.push(dna(3001, &mut rng));
        let mut mixed = dna(400, &mut rng);
        mixed.extend_from_slice(&texts[0]);
        mixed.extend_from_slice(&dna(400, &mut rng));
        reads.push(mixed);
        for read in reads {
            let mut original = d.create_streaming_query::<31>();
            let mut filtered = DeferredQuery::<31>::new(&d);
            let mut direct = DeferredQuery::<31>::new(&d);
            let mut pass = DeferredQuery::<31>::new(&d);
            for (j, w) in read.windows(31).enumerate() {
                if j == 199 {
                    original.reset();
                    filtered.reset();
                    direct.reset();
                    pass.reset();
                }
                let want = original.lookup(w);
                let got = filtered.prepare(w).unwrap_or_else(|| {
                    filtered.complete_with_minimizer_filter(|key| selected.contains(&key))
                });
                let direct_result = direct.lookup_filtered(w, |key| selected.contains(&key));
                assert_eq!(direct_result.is_found(), want.is_found());
                if want.is_found() {
                    assert_eq!(direct_result, want);
                }
                assert_eq!(pass.lookup_filtered(w, |_| true), want);
                assert_eq!(got.is_found(), want.is_found());
                if want.is_found() {
                    assert_eq!(got, want);
                }
            }
        }
    }
}
