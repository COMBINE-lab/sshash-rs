//! Absolute SPSS base offsets must follow successful streaming extensions.
//!
//! Exercise both directions, nonzero string starts, buffer refills, resets,
//! failed extensions, and recovery after ambiguous bases. Expected coordinates
//! come from the input strings as well as the independent point-query path.

use sshash_lib::builder::{config::BuildConfiguration, dictionary_builder::DictionaryBuilder};
use sshash_lib::{Kmer, KmerBits};

fn sequences() -> Vec<String> {
    let mut state = 42_u64;
    [157, 193]
        .into_iter()
        .map(|len| {
            (0..len)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[(state >> 33) as usize % 4] as char
                })
                .collect()
        })
        .collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!("test sequences contain only ACGT"),
        })
        .collect()
}

fn check_offsets<const K: usize>()
where
    Kmer<K>: KmerBits,
{
    let sequences = sequences();
    let mut cfg = BuildConfiguration::new(K, 13).unwrap();
    cfg.num_threads = 2;
    let dict = DictionaryBuilder::new(cfg)
        .unwrap()
        .build_from_sequences(sequences.clone())
        .unwrap();
    let mut engine = dict.create_streaming_query::<K>();
    let mut string_begin = 0_u64;
    for (string_id, seq) in sequences.iter().enumerate() {
        for reverse in [false, true] {
            let read = if reverse {
                reverse_complement(seq.as_bytes())
            } else {
                seq.as_bytes().to_vec()
            };
            engine.reset();
            let extensions_before = engine.num_extensions();
            for (i, window) in read.windows(K).enumerate() {
                let relative = if reverse { seq.len() - K - i } else { i } as u64;
                let expected_offset = string_begin + relative;
                let point = dict.query(&Kmer::<K>::from_ascii_unchecked(window));
                assert!(point.is_found());
                assert_eq!(point.string_id, string_id as u64);
                assert_eq!(point.kmer_offset, expected_offset);
                let streamed = engine.lookup(window);
                assert_eq!(
                    streamed.kmer_offset, expected_offset,
                    "K={K}, string={string_id}, reverse={reverse}, window={i}"
                );
                assert_eq!(streamed, point, "all positive result fields must agree");
                assert_eq!(
                    streamed.kmer_offset,
                    streamed.string_begin + streamed.kmer_id_in_string
                );
            }
            // Every lookup after the first must extend; a test consisting only
            // of dictionary seeds would miss the original stale-offset bug.
            assert_eq!(
                engine.num_extensions() - extensions_before,
                (seq.len() - K) as u64
            );
        }
        string_begin += seq.len() as u64;
    }

    // A substitution forces an extension failure, then positive reseeding.
    // N forces invalid-window resets; concatenation crosses a stored boundary.
    let mut mutated = sequences[1].as_bytes().to_vec();
    mutated[80] = if mutated[80] == b'A' { b'C' } else { b'A' };
    let joined = [sequences[0].as_bytes(), b"N", sequences[1].as_bytes()].concat();
    for read in [mutated, joined] {
        engine.reset();
        let mut misses = 0;
        let mut recovered = false;
        for window in read.windows(K) {
            let point = dict.query_from_str::<K>(std::str::from_utf8(window).unwrap());
            let streamed = engine.lookup(window);
            assert_eq!(streamed.kmer_id, point.kmer_id);
            if point.is_found() {
                assert_eq!(streamed, point);
                recovered |= misses > 0;
            } else {
                assert_eq!(streamed.kmer_offset, u64::MAX);
                misses += 1;
            }
        }
        assert!(
            misses > 0 && recovered,
            "exercise a miss and positive recovery"
        );
    }
}

#[test]
fn streaming_offsets_k31() {
    check_offsets::<31>();
}

#[test]
fn streaming_offsets_k63() {
    check_offsets::<63>();
}
