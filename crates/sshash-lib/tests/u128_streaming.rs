//! Streaming-vs-point equivalence for K > 31 (u128 k-mer storage).
//!
//! Regression test for the wide-K extension bug: `try_extend` used to extract
//! the newest base via `to_u64(bits) >> (2*(K-1))`, which truncates u128
//! storage and overflows the u64 shift for K >= 33 (panic in debug, masked
//! garbage in release). Analogue of C++ sshash commit 87cab5f. Streaming both
//! strands forces forward and backward extensions through the buggy paths.

use sshash_lib::builder::config::BuildConfiguration;
use sshash_lib::builder::dictionary_builder::DictionaryBuilder;

/// Deterministic base sequence; avoids depending on `rand`'s API surface.
fn lcg_sequences(seed: u64, num: usize, len: usize) -> Vec<String> {
    let mut state = seed | 1;
    let mut next = || {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (state >> 33) as usize
    };
    (0..num)
        .map(|_| (0..len).map(|_| b"ACGT"[next() % 4] as char).collect())
        .collect()
}

fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!("non-ACGT base"),
        })
        .collect()
}

fn check_streaming_equals_point<const K: usize>(m: usize)
where
    sshash_lib::Kmer<K>: sshash_lib::KmerBits,
{
    let sequences = lcg_sequences(0x5b_5b + K as u64, 8, 300);
    let config = BuildConfiguration::new(K, m).expect("valid config");
    let dict = DictionaryBuilder::new(config)
        .expect("builder")
        .build_from_sequences(sequences.clone())
        .expect("build");

    // Stream each sequence on both strands; the RC pass anchors matches in
    // backward orientation, exercising backward extension.
    for seq in &sequences {
        for strand in [seq.clone(), reverse_complement(seq)] {
            let mut engine = dict.create_streaming_query::<K>();
            let bytes = strand.as_bytes();
            for window in bytes.windows(K) {
                let streamed = engine.lookup(window);
                let point = dict.query_from_str::<K>(std::str::from_utf8(window).unwrap());
                assert_eq!(streamed.kmer_id, point.kmer_id, "kmer_id, K={K}");
                assert_eq!(
                    streamed.kmer_id_in_string, point.kmer_id_in_string,
                    "kmer_id_in_string, K={K}"
                );
                assert_eq!(
                    streamed.kmer_orientation, point.kmer_orientation,
                    "kmer_orientation, K={K}"
                );
                assert_eq!(streamed.string_id, point.string_id, "string_id, K={K}");
                assert_eq!(streamed.string_begin, point.string_begin, "string_begin, K={K}");
                assert_eq!(streamed.string_end, point.string_end, "string_end, K={K}");
                assert!(point.is_found(), "every input k-mer must be present, K={K}");
            }
            // Extensions must actually have happened, or this test proves nothing.
            assert!(
                engine.num_extensions() > 0,
                "streaming never extended (K={K}); the extension path is untested"
            );
        }
    }
}

#[test]
fn streaming_equals_point_k33() {
    check_streaming_equals_point::<33>(13);
}

#[test]
fn streaming_equals_point_k47() {
    check_streaming_equals_point::<47>(17);
}

#[test]
fn streaming_equals_point_k63() {
    check_streaming_equals_point::<63>(21);
}
