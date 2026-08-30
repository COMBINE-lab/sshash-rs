//! The in-memory and external-sort build paths must produce identical
//! dictionaries, and building twice from the same input must be
//! byte-identical on disk (multi-threaded builds included).

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

fn build(
    sequences: Vec<String>,
    external: bool,
    tmp: &std::path::Path,
) -> sshash_lib::Dictionary {
    const K: usize = 21;
    let mut config = BuildConfiguration::new(K, 11).unwrap();
    config.force_external_sort = external;
    config.tmp_dirname = tmp.to_path_buf();
    config.verbose = false;
    DictionaryBuilder::new(config)
        .unwrap()
        .build_from_sequences(sequences)
        .unwrap()
}

fn save_bytes(dict: &sshash_lib::Dictionary, dir: &std::path::Path, name: &str) -> (Vec<u8>, Vec<u8>) {
    let prefix = dir.join(name);
    dict.save(&prefix).unwrap();
    (
        std::fs::read(prefix.with_extension("ssi")).unwrap(),
        std::fs::read(dir.join(format!("{name}.ssi.mphf"))).unwrap(),
    )
}

#[test]
fn external_sort_path_matches_in_memory() {
    const K: usize = 21;
    let sequences = lcg_sequences(3, 16, 500);
    let dir = tempfile::tempdir().unwrap();

    let mem = build(sequences.clone(), false, dir.path());
    let ext = build(sequences.clone(), true, dir.path());

    // Identical answers on every input k-mer, both strands, plus absent probes.
    for seq in &sequences {
        for w in seq.as_bytes().windows(K) {
            let kmer = Kmer::<K>::from_ascii_unchecked(w);
            for probe in [kmer, kmer.reverse_complement()] {
                let a = mem.query::<K>(&probe);
                let b = ext.query::<K>(&probe);
                assert_eq!(a, b, "in-memory vs external-sort lookup mismatch");
                assert!(a.is_found());
            }
        }
    }
    for s in &lcg_sequences(77, 100, K) {
        let kmer = Kmer::<K>::from_ascii_unchecked(s.as_bytes());
        assert_eq!(mem.query::<K>(&kmer), ext.query::<K>(&kmer));
    }

    // And identical serialized bytes.
    let (mem_ssi, mem_mphf) = save_bytes(&mem, dir.path(), "mem");
    let (ext_ssi, ext_mphf) = save_bytes(&ext, dir.path(), "ext");
    assert_eq!(mem_ssi, ext_ssi, ".ssi bytes differ between build paths");
    assert_eq!(mem_mphf, ext_mphf, ".ssi.mphf bytes differ between build paths");
}

#[test]
fn build_is_deterministic() {
    let sequences = lcg_sequences(5, 16, 500);
    let dir = tempfile::tempdir().unwrap();

    let a = build(sequences.clone(), false, dir.path());
    let b = build(sequences, false, dir.path());

    let (a_ssi, a_mphf) = save_bytes(&a, dir.path(), "a");
    let (b_ssi, b_mphf) = save_bytes(&b, dir.path(), "b");
    assert_eq!(a_ssi, b_ssi, "two builds of the same input differ (.ssi)");
    assert_eq!(a_mphf, b_mphf, "two builds of the same input differ (.ssi.mphf)");
}
