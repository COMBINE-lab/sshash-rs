#!/usr/bin/env bash
# Old-vs-new equivalence oracle for the sshash-rs v6 port.
#
# Gates, per input x (k,m):
#  1. old point dump == new point dump, byte-identical (the semantic oracle:
#     SPSS text and string numbering are unchanged by the scheme).
#  2. new stream == new point on the six C++-compared fields (kmer_offset is
#     not maintained by extensions, in C++ either).
#  3. for K <= 31 only: old stream dump == new stream dump byte-identical.
#     (0.6.x streaming at K >= 33 is known-broken — wide-K u64 truncation,
#     fixed by this port — so old stream is not a reference there.)
set -uo pipefail
OLD=/scratch4/tmp/sshash-oracle-target/release/sshash
NEW=/scratch2/rob/salmon-rewrite/sshash-rs/target/release/sshash
DATA=/scratch3/rob/claude-13757/-scratch2-rob-salmon-rewrite/b70468e3-8092-445c-9c16-1431b9a1f53a/scratchpad/sshash-cpp/data
WORK=/scratch4/tmp/sshash-v6-oracle
cd "$WORK"
FAIL=0

six() { cut -f1-4,6-9 "$1"; }

check_queries() {
  local name=$1 k=$2 qfile=$3 tag=$4
  RUST_LOG=warn $OLD dump -i "old_$name" -q "$qfile" -o "old_${name}_${tag}_point.tsv" --point 2>/dev/null
  RUST_LOG=warn $NEW dump -i "new_$name" -q "$qfile" -o "new_${name}_${tag}_point.tsv" --point 2>/dev/null
  RUST_LOG=warn $OLD dump -i "old_$name" -q "$qfile" -o "old_${name}_${tag}_stream.tsv" 2>/dev/null
  RUST_LOG=warn $NEW dump -i "new_$name" -q "$qfile" -o "new_${name}_${tag}_stream.tsv" 2>/dev/null
  local n; n=$(wc -l < "new_${name}_${tag}_point.tsv")
  if cmp -s "old_${name}_${tag}_point.tsv" "new_${name}_${tag}_point.tsv"; then
    echo "  $tag/point old==new: IDENTICAL ($n lines)"
  else
    echo "  $tag/point old==new: *** DIFFERS ***"; FAIL=1
  fi
  if diff <(six "new_${name}_${tag}_point.tsv") <(six "new_${name}_${tag}_stream.tsv") > /dev/null; then
    echo "  $tag new stream==point (6 fields): IDENTICAL"
  else
    echo "  $tag new stream==point: *** DIFFERS ***"; FAIL=1
  fi
  if [ "$k" -le 31 ]; then
    if cmp -s "old_${name}_${tag}_stream.tsv" "new_${name}_${tag}_stream.tsv"; then
      echo "  $tag/stream old==new: IDENTICAL"
    else
      echo "  $tag/stream old==new: *** DIFFERS ***"; FAIL=1
    fi
  else
    if diff <(six "old_${name}_${tag}_point.tsv") <(six "new_${name}_${tag}_stream.tsv") > /dev/null; then
      echo "  $tag new stream == old point (6 fields): IDENTICAL (old stream known-broken at K>31)"
    else
      echo "  $tag new stream vs old point: *** DIFFERS ***"; FAIL=1
    fi
  fi
}

run_case() {
  local name=$1 input=$2 k=$3 m=$4 reads=$5
  echo "=== $name (k=$k m=$m) ==="
  RUST_LOG=warn $OLD build -i "$input" -k "$k" -m "$m" --canonical -o "old_$name" -r 0 2>/dev/null
  RUST_LOG=warn $NEW build -i "$input" -k "$k" -m "$m" -o "new_$name" -r 0 2>/dev/null
  check_queries "$name" "$k" "$input" self
  [ -n "$reads" ] && check_queries "$name" "$k" "$reads" reads
}

run_case sal_m13 "$DATA/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz" 31 13 "$DATA/queries/SRR5833294.10K.fastq.gz"
run_case sal_m16 "$DATA/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz" 31 16 "$DATA/queries/SRR5833294.10K.fastq.gz"
run_case sal_m20 "$DATA/unitigs_stitched/salmonella_enterica_k31_ust.fa.gz" 31 20 "$DATA/queries/salmonella_enterica.fasta.gz"
run_case ecoli_m16 "$DATA/unitigs_stitched/ecoli1_k31_ust.fa.gz" 31 16 "$DATA/queries/ecoli1.fasta.gz"
run_case sek47_m21 "$DATA/unitigs_stitched/se.ust.k47.fa.gz" 47 21 "$WORK/reads_err_k47.fa"
run_case sek63_m21 "$DATA/unitigs_stitched/se.ust.k63.fa.gz" 63 21 "$WORK/reads_err_k63.fa"
if [ "$FAIL" -eq 0 ]; then echo "ORACLE: ALL GATES PASS"; else echo "ORACLE: FAILURES PRESENT"; exit 1; fi
