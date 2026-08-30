//! Dictionary save/load round-trip and format-version gating.

use sshash_lib::builder::config::BuildConfiguration;
use sshash_lib::builder::dictionary_builder::DictionaryBuilder;
use sshash_lib::kmer::Kmer;

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

#[test]
fn save_load_roundtrip_preserves_every_lookup() {
    const K: usize = 21;
    let sequences = lcg_sequences(7, 12, 400);
    let config = BuildConfiguration::new(K, 11).unwrap();
    let dict = DictionaryBuilder::new(config)
        .unwrap()
        .build_from_sequences(sequences.clone())
        .unwrap();

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("rt_index");
    dict.save(&prefix).unwrap();
    let loaded = sshash_lib::Dictionary::load(&prefix).unwrap();

    assert_eq!(loaded.k(), dict.k());
    assert_eq!(loaded.m(), dict.m());
    assert_eq!(loaded.seed(), dict.seed());
    assert_eq!(loaded.num_strings(), dict.num_strings());

    // Every input k-mer (both strands) answers identically pre- and post-save.
    for seq in &sequences {
        for w in seq.as_bytes().windows(K) {
            let kmer = Kmer::<K>::from_ascii_unchecked(w);
            for probe in [kmer, kmer.reverse_complement()] {
                let a = dict.query::<K>(&probe);
                let b = loaded.query::<K>(&probe);
                assert_eq!(a, b, "roundtrip lookup mismatch");
                assert!(a.is_found());
            }
        }
    }

    // Absent k-mers agree too (sampled).
    let absent = lcg_sequences(99, 50, K);
    for s in &absent {
        let kmer = Kmer::<K>::from_ascii_unchecked(s.as_bytes());
        assert_eq!(dict.query::<K>(&kmer), loaded.query::<K>(&kmer));
    }
}

#[test]
fn nondefault_seed_roundtrip() {
    // Regression for the latent seed bug: a non-default build seed must
    // produce a queryable index, before and after save/load.
    const K: usize = 15;
    let sequences = lcg_sequences(11, 4, 200);
    let mut config = BuildConfiguration::new(K, 9).unwrap();
    config.seed = 0xC0FFEE;
    let dict = DictionaryBuilder::new(config)
        .unwrap()
        .build_from_sequences(sequences.clone())
        .unwrap();
    assert_eq!(dict.seed(), 0xC0FFEE);

    let kmer = Kmer::<K>::from_ascii_unchecked(&sequences[0].as_bytes()[..K]);
    assert!(dict.query::<K>(&kmer).is_found(), "non-default seed index must answer");

    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("seed_index");
    dict.save(&prefix).unwrap();
    let loaded = sshash_lib::Dictionary::load(&prefix).unwrap();
    assert_eq!(loaded.seed(), 0xC0FFEE);
    assert!(loaded.query::<K>(&kmer).is_found());
}

#[test]
fn rejects_old_format_major_with_rebuild_message() {
    // Craft a header that mimics a (4,0)-era file: same magic, old major.
    // The loader must reject it with an actionable "rebuilt" message.
    use std::io::Write;
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("old_index");
    let ssi = prefix.with_extension("ssi");
    {
        let mut f = std::fs::File::create(&ssi).unwrap();
        f.write_all(b"SSHIDX01").unwrap();
        f.write_all(&4u32.to_le_bytes()).unwrap(); // old major
        f.write_all(&0u32.to_le_bytes()).unwrap();
        f.write_all(&31u64.to_le_bytes()).unwrap(); // k
        f.write_all(&19u64.to_le_bytes()).unwrap(); // m
        f.write_all(&[1u8]).unwrap(); // canonical
        f.write_all(&1u32.to_le_bytes()).unwrap(); // partitions
    }
    // The mphf file must exist for load to get that far — but the header is
    // read first, so its absence is fine.
    let msg = match sshash_lib::Dictionary::load(&prefix) {
        Ok(_) => panic!("loading a (4,0)-era index must fail"),
        Err(e) => e.to_string(),
    };
    assert!(
        msg.contains("rebuilt"),
        "old-format rejection must tell the user to rebuild, got: {msg}"
    );
}
