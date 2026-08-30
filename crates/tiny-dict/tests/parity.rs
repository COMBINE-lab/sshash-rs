//! Parity test: TinyDictionary vs sshash Dictionary.
//!
//! Builds a small canonical sshash Dictionary from synthetic sequences,
//! converts to TinyDictionary, and asserts that every canonical k-mer in
//! the SPSS returns identical LookupResult fields from both dictionaries.
//! Also samples random absent k-mers.

use rand::{rngs::StdRng, RngExt, SeedableRng};
use sshash_lib::builder::{BuildConfiguration, DictionaryBuilder};
use sshash_lib::{Kmer, KmerDictionary, KmerStreamingQuery, LookupResult};
use tiny_dict::TinyDictionary;

const K: usize = 23;
const M: usize = 11;

fn random_dna(rng: &mut StdRng, len: usize) -> String {
    const BASES: &[u8] = b"ACGT";
    (0..len)
        .map(|_| BASES[rng.random_range(0..4)] as char)
        .collect()
}

fn build_sshash(sequences: Vec<String>) -> sshash_lib::Dictionary {
    let config = BuildConfiguration {
        k: K,
        m: M,
        num_threads: 1,
        verbose: false,
        ..BuildConfiguration::default()
    };
    let builder = DictionaryBuilder::new(config).expect("build config");
    builder
        .build_from_sequences(sequences)
        .expect("sshash build")
}

/// Compare the subset of LookupResult fields that TinyDictionary can
/// reconstruct. `minimizer_found` is skipped (TinyDict always sets true).
fn assert_result_equal(a: &LookupResult, b: &LookupResult, kmer: &str) {
    assert_eq!(a.is_found(), b.is_found(), "is_found mismatch at {kmer}");
    if !a.is_found() {
        return;
    }
    assert_eq!(a.kmer_id, b.kmer_id, "kmer_id mismatch at {kmer}");
    assert_eq!(
        a.kmer_id_in_string, b.kmer_id_in_string,
        "kmer_id_in_string mismatch at {kmer}"
    );
    assert_eq!(
        a.kmer_offset, b.kmer_offset,
        "kmer_offset mismatch at {kmer}"
    );
    assert_eq!(
        a.kmer_orientation, b.kmer_orientation,
        "kmer_orientation mismatch at {kmer}"
    );
    assert_eq!(a.string_id, b.string_id, "string_id mismatch at {kmer}");
    assert_eq!(
        a.string_begin, b.string_begin,
        "string_begin mismatch at {kmer}"
    );
    assert_eq!(a.string_end, b.string_end, "string_end mismatch at {kmer}");
}

#[test]
fn tiny_dict_parity_all_kmers() {
    let mut rng = StdRng::seed_from_u64(0xC0FFEE);
    // A handful of small sequences, each well above K.
    let sequences: Vec<String> = (0..5).map(|_| random_dna(&mut rng, 200)).collect();

    let sshash = build_sshash(sequences.clone());
    let tiny = TinyDictionary::from_sshash::<K>(&sshash);

    assert_eq!(tiny.k(), sshash.k());
    assert_eq!(tiny.m(), sshash.m());
    assert!(tiny.canonical());

    let mut found_count = 0u64;
    for seq in &sequences {
        let bytes = seq.as_bytes();
        if bytes.len() < K {
            continue;
        }
        for start in 0..=(bytes.len() - K) {
            let kmer_bytes = &bytes[start..start + K];
            let kmer_str = std::str::from_utf8(kmer_bytes).unwrap();
            let fw = Kmer::<K>::from_ascii_unchecked(kmer_bytes);

            let ss_res = sshash.query::<K>(&fw);
            let tiny_res = tiny.lookup::<K>(kmer_bytes);

            if ss_res.is_found() {
                found_count += 1;
            }
            assert_result_equal(&ss_res, &tiny_res, kmer_str);
        }
    }
    assert!(found_count > 0, "no k-mers found — bad test");
}

#[test]
fn tiny_dict_parity_random_absent() {
    let mut rng = StdRng::seed_from_u64(0xBEEFCAFE);
    // Build from one sequence, then query many random k-mers. Vast
    // majority will be absent.
    let seq = random_dna(&mut rng, 500);
    let sshash = build_sshash(vec![seq]);
    let tiny = TinyDictionary::from_sshash::<K>(&sshash);

    for _ in 0..500 {
        let kmer = random_dna(&mut rng, K);
        let kmer_bytes = kmer.as_bytes();
        let fw = Kmer::<K>::from_ascii_unchecked(kmer_bytes);
        let ss_res = sshash.query::<K>(&fw);
        let tiny_res = tiny.lookup::<K>(kmer_bytes);
        assert_result_equal(&ss_res, &tiny_res, &kmer);
    }
}

#[test]
fn tiny_dict_tdct_roundtrip() {
    let mut rng = StdRng::seed_from_u64(0xD07D07);
    let sequences: Vec<String> = (0..4).map(|_| random_dna(&mut rng, 300)).collect();
    let sshash = build_sshash(sequences.clone());
    let tiny = TinyDictionary::from_sshash::<K>(&sshash);

    let tmp = tempfile::NamedTempFile::new().expect("tempfile");
    let path = tmp.path().to_path_buf();
    tiny.save(&path).expect("save");
    let loaded = TinyDictionary::load(&path).expect("load");

    assert_eq!(loaded.k(), tiny.k());
    assert_eq!(loaded.m(), tiny.m());
    assert_eq!(loaded.canonical(), tiny.canonical());

    for seq in &sequences {
        let bytes = seq.as_bytes();
        if bytes.len() < K {
            continue;
        }
        for start in 0..=(bytes.len() - K) {
            let kmer_bytes = &bytes[start..start + K];
            let kmer_str = std::str::from_utf8(kmer_bytes).unwrap();
            let a = tiny.lookup::<K>(kmer_bytes);
            let b = loaded.lookup::<K>(kmer_bytes);
            assert_result_equal(&a, &b, kmer_str);
        }
    }
}

#[test]
fn tiny_dict_streaming_query_matches_one_shot() {
    let mut rng = StdRng::seed_from_u64(42);
    let sequences: Vec<String> = (0..3).map(|_| random_dna(&mut rng, 150)).collect();
    let sshash = build_sshash(sequences.clone());
    let tiny = TinyDictionary::from_sshash::<K>(&sshash);

    let mut q = tiny.create_streaming_query::<K>();
    for seq in &sequences {
        let bytes = seq.as_bytes();
        if bytes.len() < K {
            continue;
        }
        q.reset();
        for start in 0..=(bytes.len() - K) {
            let kmer_bytes = &bytes[start..start + K];
            let kmer_str = std::str::from_utf8(kmer_bytes).unwrap();
            let streamed = q.lookup(kmer_bytes);
            let one_shot = tiny.lookup::<K>(kmer_bytes);
            assert_result_equal(&streamed, &one_shot, kmer_str);
        }
    }
}

/// The rapidhash version recorded in a `.tdct` header is checked against an
/// allow-list of versions verified to hash identically, not against the exact
/// version this binary links. Rewriting the header's three version bytes is
/// enough to exercise both sides of that check: the payload is untouched, so a
/// file naming an accepted version must still load and answer identically,
/// while one naming an unaccepted version must be refused up front.
#[test]
fn tdct_rapidhash_version_allow_list() {
    let mut rng = StdRng::seed_from_u64(0x4A91D0);
    let sequences: Vec<String> = (0..4).map(|_| random_dna(&mut rng, 300)).collect();
    let sshash = build_sshash(sequences.clone());
    let tiny = TinyDictionary::from_sshash::<K>(&sshash);

    let tmp = tempfile::NamedTempFile::new().expect("tempfile");
    let path = tmp.path().to_path_buf();
    tiny.save(&path).expect("save");
    let original = std::fs::read(&path).expect("read back");

    // Bytes 10..13 are the rapidhash semver; everything else stays byte-identical.
    let rewrite_version = |ver: [u8; 3]| {
        let mut bytes = original.clone();
        bytes[10..13].copy_from_slice(&ver);
        let f = tempfile::NamedTempFile::new().expect("tempfile");
        std::fs::write(f.path(), &bytes).expect("write");
        f
    };

    // 4.4.0 is on the allow-list: a file stamped with it still loads, and every
    // k-mer answers exactly as it does through the freshly built dictionary.
    let accepted = rewrite_version([4, 4, 0]);
    let loaded = TinyDictionary::load(accepted.path()).expect("accepted version must load");
    for seq in &sequences {
        let bytes = seq.as_bytes();
        if bytes.len() < K {
            continue;
        }
        for start in 0..=(bytes.len() - K) {
            let kmer_bytes = &bytes[start..start + K];
            let kmer_str = std::str::from_utf8(kmer_bytes).unwrap();
            let a = tiny.lookup::<K>(kmer_bytes);
            let b = loaded.lookup::<K>(kmer_bytes);
            assert_result_equal(&a, &b, kmer_str);
        }
    }

    // Anything not on the list is rejected before any bulk work, so widening the
    // check cannot degrade into silently trusting an unknown hash.
    for ver in [[4, 3, 0], [5, 0, 0]] {
        let rejected = rewrite_version(ver);
        match TinyDictionary::load(rejected.path()) {
            Err(tiny_dict::io::TdctError::RapidhashVersionMismatch { got, .. }) => {
                assert_eq!(got, ver, "error should report the version it found");
            }
            Err(other) => panic!("rapidhash {ver:?} must fail the version check, got {other:?}"),
            Ok(_) => panic!("rapidhash {ver:?} must be refused, but the file loaded"),
        }
    }
}
