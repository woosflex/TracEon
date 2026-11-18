// Experimental FASTA-only adapter for the mapping pilot; not a mapper index.
#include "Cache.h"
#include "BenchmarkUtils.h"
#ifdef TRACEON_SUBSET_REMOTE
#include "TraceonClient.h"
#endif
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <unordered_set>
#include <vector>

int main(int argc, char** argv) {
    const bool remote = argc == 5 && std::string(argv[1]) == "--remote";
    if (argc != 3 && !remote) {
        std::cerr << "Usage: traceon_subset cache.traceon ids.txt > subset.fa\n"
                     "       traceon_subset --remote host port ids.txt > subset.fa\n"
                     "Remote mode emits sequences only; no FASTQ qualities or header descriptions.\n";
        return 2;
    }
    try {
        std::ifstream input(argv[remote ? 4 : 2]);
        if (!input) throw std::runtime_error("cannot open ID list");
        std::vector<std::string> ids;
        std::unordered_set<std::string> seen;
        for (std::string id; std::getline(input, id);) {
            if (!id.empty() && id.back() == '\r') id.pop_back();
            if (id.empty() || id.find_first_of(" \t\r") != std::string::npos)
                throw std::runtime_error("expected one nonempty identifier per line");
            if (!seen.insert(id).second) throw std::runtime_error("duplicate ID: " + id);
            ids.push_back(std::move(id));
        }
        if (input.bad()) throw std::runtime_error("failed reading ID list");
        if (ids.empty()) throw std::runtime_error("empty ID list");
        if (remote) {
#ifdef TRACEON_SUBSET_REMOTE
            const auto port = benchmarkUnsigned(argv[3], "port");
            if (port == 0 || port > 65535) throw std::runtime_error("port must be 1..65535");
            TracEon::TraceonClient client(argv[2], static_cast<uint16_t>(port));
            std::vector<std::string> sequences;
            for (const auto& id : ids) {
                auto sequence = client.getView(id);
                if (!sequence) throw std::runtime_error("missing ID: " + id);
                sequences.push_back(std::move(*sequence));
            }
            for (size_t i = 0; i < ids.size(); ++i)
                std::cout << '>' << ids[i] << '\n' << sequences[i] << '\n';
#else
            throw std::runtime_error("remote support is disabled in this build");
#endif
        } else {
            TracEon::Cache cache;
            cache.restore(argv[1]);
            const auto format = cache.getDetectedFormat();
            if (format != TracEon::FileFormat::DNA_FASTA &&
                format != TracEon::FileFormat::RNA_FASTA &&
                format != TracEon::FileFormat::PROTEIN_FASTA)
                throw std::runtime_error("pilot accepts FASTA only; FASTQ qualities must not be discarded");
            for (const auto& id : ids)
                if (!cache.hasSequence(id)) throw std::runtime_error("missing ID: " + id);
            // Emit in request order; the caller must choose source order for seqtk parity.
            for (const auto& id : ids) std::cout << '>' << id << '\n' << cache.getView(id) << '\n';
        }
        std::cout.flush();
        if (!std::cout) throw std::runtime_error("failed writing FASTA");
    } catch (const std::exception& error) {
        std::cerr << "traceon_subset: " << error.what() << '\n';
        return 1;
    }
}
