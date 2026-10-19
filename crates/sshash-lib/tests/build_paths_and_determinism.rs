//! The in-memory and external-sort build paths must produce identical
//! dictionaries, and building twice from the same input must be
//! byte-identical on disk (multi-threaded builds included).

use sshash_lib::builder::config::BuildConfiguration;
use sshash_lib::builder::dictionary_builder::DictionaryBuilder;
use sshash_lib::kmer::Kmer;
use sshash_lib::{MphfBuild, Offline};

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
    build_with(sequences, external, MphfBuild::default(), 0, tmp)
}

fn build_with(
    sequences: Vec<String>,
    external: bool,
    mphf_build: MphfBuild,
    num_threads: usize,
    tmp: &std::path::Path,
) -> sshash_lib::Dictionary {
    const K: usize = 21;
    let mut config = BuildConfiguration::new(K, 11).unwrap();
    config.force_external_sort = external;
    config.mphf_build = mphf_build;
    config.num_threads = num_threads;
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

const STRATEGIES: [MphfBuild; 5] = [
    MphfBuild::Streaming(Offline::Auto { min_keys: 1 << 24 }),
    MphfBuild::Streaming(Offline::Always),
    MphfBuild::Streaming(Offline::Never),
    MphfBuild::InMemory,
    MphfBuild::Partitioned,
];

/// Every MPHF strategy, on both build paths and at any thread count, builds
/// a dictionary that answers every query the same way; each strategy's bytes
/// do not depend on the build path or `num_threads`, and offline streaming
/// builds the same function as in-memory streaming.
#[test]
fn mphf_strategies_agree_and_are_deterministic() {
    const K: usize = 21;
    let sequences = lcg_sequences(11, 24, 600);
    let probes: Vec<Kmer<K>> = sequences
        .iter()
        .flat_map(|s| s.as_bytes().windows(K).map(Kmer::<K>::from_ascii_unchecked).collect::<Vec<_>>())
        .chain(lcg_sequences(78, 200, K).iter().map(|s| Kmer::<K>::from_ascii_unchecked(s.as_bytes())))
        .collect();
    let dir = tempfile::tempdir().unwrap();

    let reference = build(sequences.clone(), false, dir.path());
    let mut streaming_bytes = None;
    for (si, strategy) in STRATEGIES.into_iter().enumerate() {
        let mut bytes = None;
        for (external, threads) in [(false, 1), (false, 4), (true, 1), (true, 4)] {
            let dict = build_with(sequences.clone(), external, strategy, threads, dir.path());
            for probe in &probes {
                for p in [*probe, probe.reverse_complement()] {
                    assert_eq!(reference.query::<K>(&p), dict.query::<K>(&p), "{strategy:?} answers differ");
                }
            }
            let b = save_bytes(&dict, dir.path(), &format!("s{si}-{external}-{threads}"));
            match &bytes {
                None => bytes = Some(b),
                Some(first) => assert!(
                    *first == b,
                    "{strategy:?}: bytes depend on the build path or num_threads (external={external}, threads={threads})"
                ),
            }
        }
        if matches!(strategy, MphfBuild::Streaming(_)) {
            match &streaming_bytes {
                None => streaming_bytes = bytes,
                Some(first) => assert!(Some(first) == bytes.as_ref(), "{strategy:?}: offline mode changed the index"),
            }
        }
    }
}
