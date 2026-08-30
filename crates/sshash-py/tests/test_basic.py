"""Self-contained runtime tests for the sshash python bindings.

Run with: maturin develop --release && pytest crates/sshash-py/tests
No external data needed — the index is built from sequences generated here.
"""
import random
import warnings

import pytest
import sshash

COMP = str.maketrans("ACGT", "TGCA")
K, M = 21, 11


def rc(s):
    return s.translate(COMP)[::-1]


@pytest.fixture(scope="module")
def seqs():
    rng = random.Random(0x5B)
    return ["".join(rng.choice("ACGT") for _ in range(400)) for _ in range(6)]


@pytest.fixture(scope="module")
def dct(seqs):
    cfg = sshash.BuildConfig(k=K, m=M)
    cfg.threads = 1
    cfg.verbose = False
    return cfg.build(list(seqs))


def test_shape(dct, seqs):
    assert dct.k == K and dct.m == M
    assert dct.canonical is True
    assert dct.num_strings == len(seqs)


def test_point_queries_both_strands(dct, seqs):
    for sid, s in enumerate(seqs):
        for i in range(0, len(s) - K + 1, 5):
            kmer = s[i : i + K]
            h = dct.query(kmer)
            assert h is not None, f"missing k-mer seq {sid} pos {i}"
            assert h.string_id == sid and h.kmer_id_in_string == i
            assert h.orientation == 1
            hr = dct.query(rc(kmer))
            assert hr is not None and hr.kmer_id == h.kmer_id
            assert hr.orientation == -1 and hr.string_id == sid
            assert dct.contains(kmer) and dct.contains(rc(kmer))


def test_absent_kmer(dct):
    assert dct.query("A" * K) is None
    assert not dct.contains("A" * K)


def test_streaming_equals_point(dct, seqs):
    sq = dct.streaming_query()
    for s in seqs:
        sq.reset()
        hits = sq.query_sequence(s.encode())
        assert len(hits) == len(s) - K + 1
        for i, h in enumerate(hits):
            p = dct.query(s[i : i + K])
            assert h is not None and p is not None
            assert (h.kmer_id, h.string_id, h.kmer_id_in_string, h.orientation) == (
                p.kmer_id, p.string_id, p.kmer_id_in_string, p.orientation)
    assert sq.num_extensions > 0


def test_forward_only_restricts_orientation(dct, seqs):
    # C++ v6 semantics (check_reverse_complement=False): the single-probe
    # lookup runs as usual, but a match found in backward orientation
    # reports absent. Forward-orientation k-mers still answer.
    kmer = seqs[0][:K]
    assert dct.query(kmer, forward_only=True) is not None
    assert dct.lookup(kmer, forward_only=True) is not None
    assert dct.contains(kmer, forward_only=True)
    assert dct.query(rc(kmer), forward_only=True) is None
    assert dct.lookup(rc(kmer), forward_only=True) is None
    assert not dct.contains(rc(kmer), forward_only=True)
    # The unrestricted query still answers both strands.
    assert dct.query(rc(kmer)) is not None
    # And forward_only emits no warnings — it is a real, supported option.
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        dct.query(kmer, forward_only=True)
        dct.query(rc(kmer), forward_only=True)
    assert not [x for x in w if issubclass(x.category, DeprecationWarning)]


def test_canonical_setter_deprecation():
    cfg = sshash.BuildConfig(k=K, m=M)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        cfg.canonical = False  # ignored either way
    dep = [x for x in w if issubclass(x.category, DeprecationWarning)]
    assert len(dep) == 1
    assert cfg.canonical is True


def test_save_load_roundtrip(dct, seqs, tmp_path):
    prefix = str(tmp_path / "pydict")
    dct.save(prefix)
    re = sshash.Dictionary.load(prefix)
    for s in seqs[:2]:
        for i in range(0, len(s) - K + 1, 17):
            a, b = dct.query(s[i : i + K]), re.query(s[i : i + K])
            assert (a.kmer_id, a.string_id, a.orientation) == (b.kmer_id, b.string_id, b.orientation)
