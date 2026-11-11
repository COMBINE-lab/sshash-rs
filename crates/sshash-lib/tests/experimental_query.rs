#![cfg(feature = "experimental-query")]
//! Differential contracts for the opt-in query paths. Compare every result
//! field, including conservative negative hints, and the serialized bytes.
//! Small minimizers exercise heavy buckets; even m exercises palindromic minima;
//! k > 31 exercises the u128 representation. These are correctness tests, not
//! timing assertions. The scalar/streaming baselines are kept separate because
//! their pre-existing kmer_offset behavior can differ on successful extensions.
use sshash_lib::builder::{config::BuildConfiguration, dictionary_builder::DictionaryBuilder};
use sshash_lib::kmer::KmerBits;
use sshash_lib::{BucketFingerprintLayout, Dictionary, Kmer};

fn dna(n: usize, rng: &mut u64) -> String {
    (0..n)
        .map(|_| {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(*rng >> 32) as usize % 4] as char
        })
        .collect()
}

fn check<const K: usize>(m: usize)
where
    Kmer<K>: KmerBits,
{
    let mut rng = 937;
    let sequences: Vec<_> = (0..40).map(|_| dna(503, &mut rng)).collect();
    let mut cfg = BuildConfiguration::new(K, m).unwrap();
    cfg.num_threads = 2;
    let mut d = DictionaryBuilder::new(cfg)
        .unwrap()
        .build_from_sequences(sequences.clone())
        .unwrap();
    let mut queries = Vec::new();
    let mut reads = sequences.clone();
    reads.push(dna(10007, &mut rng));
    for r in &reads {
        for w in r.as_bytes().windows(K).step_by(7) {
            let k = Kmer::<K>::from_ascii_unchecked(w);
            queries.extend([k, k.reverse_complement()]);
        }
    }
    let want: Vec<_> = queries.iter().map(|k| d.query(k)).collect();
    reads.push(sequences[0].to_ascii_lowercase());
    reads.push(format!("{}N{}", &sequences[0][..100], &sequences[1][..100]));
    reads.push(String::new());
    reads.push("AC".into());
    let stream = |d: &Dictionary| {
        reads
            .iter()
            .flat_map(|r| {
                let mut q = d.create_streaming_query::<K>();
                r.as_bytes()
                    .windows(K)
                    .map(|w| q.lookup(w))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let want_stream = stream(&d);
    let tmp = tempfile::tempdir().unwrap();
    let before = tmp.path().join("before");
    let after = tmp.path().join("after");
    d.save(&before).unwrap();
    let base_bits = d.num_bits();
    for layout in [
        None,
        Some(BucketFingerprintLayout::Packed4),
        Some(BucketFingerprintLayout::Bits7),
        Some(BucketFingerprintLayout::Byte255),
        Some(BucketFingerprintLayout::Wide16),
    ] {
        let n = d.num_minimizers() as usize;
        let bytes = match layout {
            None => 0,
            Some(BucketFingerprintLayout::Packed4) => n.div_ceil(2),
            Some(BucketFingerprintLayout::Wide16) => n * 2,
            _ => n,
        };
        if let Some(layout) = layout {
            assert_eq!(d.enable_bucket_fingerprints_with_layout(layout), bytes);
        }
        assert_eq!(d.num_bits(), base_bits + bytes as u64 * 8);
        assert_eq!(
            stream(&d),
            want_stream,
            "stream fields changed, k={K}, m={m}"
        );
        assert_eq!(queries.iter().map(|k| d.query(k)).collect::<Vec<_>>(), want);
        for (i, c) in queries.chunks(16).enumerate() {
            let mut batch = [queries[0]; 16];
            batch[..c.len()].copy_from_slice(c);
            for got in [
                d.query_batch::<K, 16, true>(&batch),
                d.query_batch::<K, 16, false>(&batch),
            ] {
                assert_eq!(&got[..c.len()], &want[i * 16..i * 16 + c.len()]);
            }
        }
        assert!(d.query_batch::<K, 0, true>(&[]).is_empty());
        d.save(&after).unwrap();
        for ext in ["ssi", "ssi.mphf"] {
            assert_eq!(
                std::fs::read(before.with_extension(ext)).unwrap(),
                std::fs::read(after.with_extension(ext)).unwrap()
            );
        }
    }
    d.save(&after).unwrap();
    for ext in ["ssi", "ssi.mphf"] {
        assert_eq!(
            std::fs::read(before.with_extension(ext)).unwrap(),
            std::fs::read(after.with_extension(ext)).unwrap()
        );
    }
    let loaded = Dictionary::load(&after).unwrap();
    assert_eq!(
        queries.iter().map(|k| loaded.query(k)).collect::<Vec<_>>(),
        want
    );
    d.disable_bucket_fingerprints();
    assert_eq!(stream(&d), want_stream);
}

#[test]
fn exact_fields_strands_heavy_buckets_and_serialized_bytes() {
    check::<21>(3); // few minimizers: exercises heavy buckets
    check::<31>(18); // even m includes self reverse-complement minimizers
    check::<33>(19); // cross the u64/u128 boundary
    check::<63>(31); // maximum supported k
}

#[test]
fn empty_dictionary_and_empty_batch() {
    let mut d = Dictionary::new(
        sshash_lib::spectrum_preserving_string_set::SpectrumPreservingStringSet::new(3, 2),
        sshash_lib::minimizers_control_map::MinimizersControlMap::empty(),
        sshash_lib::sparse_and_skew_index::SparseAndSkewIndex::new(),
        3,
        2,
        1,
    );
    let keys = [Kmer::<3>::from_str("ACG").unwrap(); 4];
    let want = keys.map(|k| d.query(&k));
    assert_eq!(d.query_batch::<3, 4, true>(&keys), want);
    assert_eq!(d.enable_bucket_fingerprints(), 0);
    assert_eq!(d.query_batch::<3, 4, true>(&keys), want);
    assert!(d.query_batch::<3, 0, true>(&[]).is_empty());
}

#[test]
fn self_reverse_complement_minimizer() {
    let sequence = "ATATATATATATATATATATATATATATATA";
    let cfg = BuildConfiguration::new(31, 18).unwrap();
    let mut d = DictionaryBuilder::new(cfg)
        .unwrap()
        .build_from_sequences(vec![sequence.into()])
        .unwrap();
    let k = Kmer::<31>::from_str(sequence).unwrap();
    let m = d.debug_extract_minimizer(&k).value;
    assert_eq!(m, sshash_lib::minimizer::reverse_complement_mmer(m, 18));
    let keys = [
        k,
        k.reverse_complement(),
        Kmer::from_str("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA").unwrap(),
        k,
    ];
    let want = keys.map(|k| d.query(&k));
    d.enable_bucket_fingerprints();
    assert_eq!(d.query_batch::<31, 4, true>(&keys), want);
}
