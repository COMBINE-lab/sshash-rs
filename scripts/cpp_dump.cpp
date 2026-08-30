// Scratchpad-only per-kmer dump driver for cross-implementation diffing.
// Prints the same TSV as sshash-rs `sshash dump --point`.
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "include/dictionary_types.hpp"
#include "src/dictionary.cpp"

using namespace sshash;

static bool valid_window(const char* p, uint64_t k) {
    for (uint64_t i = 0; i != k; ++i) {
        char c = p[i];
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: cpp_dump <index> <query.fa (plain)> <out.tsv>\n";
        return 1;
    }
    dictionary_type dict;
    essentials::load(dict, argv[1]);
    const uint64_t k = dict.k();

    std::ifstream in(argv[2]);
    if (!in.good()) { std::cerr << "bad query file\n"; return 1; }
    std::FILE* out = std::fopen(argv[3], "w");
    std::fprintf(out,
        "seq_idx\tkmer_pos\tkmer_id\tkmer_id_in_string\tkmer_offset\tkmer_orientation\tstring_id\tstring_begin\tstring_end\n");

    std::string line, seq;
    long seq_idx = -1;
    auto flush_seq = [&]() {
        if (seq.empty()) return;
        for (char& c : seq) c = std::toupper(c);
        if (seq.size() >= k) {
            for (uint64_t i = 0; i + k <= seq.size(); ++i) {
                lookup_result res;
                if (valid_window(seq.data() + i, k)) {
                    res = dict.lookup(seq.data() + i);
                } else {
                    res = lookup_result();
                }
                std::fprintf(out, "%ld\t%llu\t%llu\t%llu\t%llu\t%lld\t%llu\t%llu\t%llu\n",
                    seq_idx, (unsigned long long)i,
                    (unsigned long long)res.kmer_id,
                    (unsigned long long)res.kmer_id_in_string,
                    (unsigned long long)res.kmer_offset,
                    (long long)res.kmer_orientation,
                    (unsigned long long)res.string_id,
                    (unsigned long long)res.string_begin,
                    (unsigned long long)res.string_end);
            }
        }
        seq.clear();
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            flush_seq();
            ++seq_idx;
        } else {
            seq += line;
        }
    }
    flush_seq();
    std::fclose(out);
    return 0;
}
