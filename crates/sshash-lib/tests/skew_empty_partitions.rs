//! Empty skew size bands must not shift the partitions used by later buckets.

use sshash_lib::{BuildConfiguration, Dictionary, DictionaryBuilder, Kmer};

const K: usize = 31;

fn sequences() -> Vec<String> {
    // This fixed input has heavy buckets in (64, 128] and above 256, but none
    // in (128, 256]. Omitting that interior skew slot loses 6,218 source queries
    // across the two strands before any streaming or serialization checks.
    let mut state = 813_u64;
    (0..40)
        .map(|_| {
            (0..600)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    b"ACGT"[((state >> 33) & 3) as usize] as char
                })
                .collect()
        })
        .collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!("the fixture contains only ACGT"),
        })
        .collect()
}

fn check_queries(dict: &Dictionary, reads: &[Vec<u8>]) {
    // Source membership is an independent oracle: agreement between two
    // lookup paths alone would not detect a shared false negative.
    let misses = reads
        .iter()
        .flat_map(|read| read.windows(K))
        .filter(|window| {
            !dict
                .query(&Kmer::<K>::from_ascii_unchecked(window))
                .is_found()
        })
        .count();
    assert_eq!(misses, 0, "source k-mers must be found on both strands");

    let mut stream = dict.create_streaming_query::<K>();
    for read in reads {
        stream.reset();
        for window in read.windows(K) {
            let point = dict.query(&Kmer::<K>::from_ascii_unchecked(window));
            assert_eq!(stream.lookup(window), point);
        }
    }
}

fn check_build(force_external_sort: bool) {
    let sequences = sequences();
    let reads: Vec<_> = sequences
        .iter()
        .flat_map(|seq| [seq.as_bytes().to_vec(), reverse_complement(seq.as_bytes())])
        .collect();
    let dir = tempfile::tempdir().unwrap();
    let mut cfg = BuildConfiguration::new(K, 3).unwrap();
    cfg.num_threads = 3;
    cfg.force_external_sort = force_external_sort;
    cfg.tmp_dirname = dir.path().join("sort");
    std::fs::create_dir_all(&cfg.tmp_dirname).unwrap();
    let dict = DictionaryBuilder::new(cfg)
        .unwrap()
        .build_from_sequences(sequences)
        .unwrap();

    let index = dict.index_ref();
    let sizes: Vec<_> = (0..index.num_buckets())
        .map(|id| {
            let (begin, end) = index.locate_bucket(id);
            end - begin
        })
        .collect();
    assert!(sizes.iter().any(|size| (65..=128).contains(size)));
    assert!(sizes.iter().any(|&size| size > 256));
    assert!(!sizes.iter().any(|size| (129..=256).contains(size)));
    check_queries(&dict, &reads);

    let prefix = dir.path().join("index");
    dict.save(&prefix).unwrap();
    let loaded = Dictionary::load(&prefix).unwrap();
    check_queries(&loaded, &reads);
}

#[test]
fn classified_build_preserves_empty_skew_size_bands() {
    check_build(false);
}

#[test]
fn external_build_preserves_empty_skew_size_bands() {
    check_build(true);
}
