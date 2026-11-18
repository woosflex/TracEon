// Broader control, derived from resident_record_worker.cpp; adds public pyfastx API.
#include <Python.h>
// Controlled local benchmark protocol; no public service interface is added.
#include "TraceonClient.h"
#include <htslib/faidx.h>
#include <openssl/evp.h>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>

using Clock = std::chrono::steady_clock;
static double seconds(Clock::duration d) { return std::chrono::duration<double>(d).count(); }
static double stamp(Clock::time_point t) { return seconds(t.time_since_epoch()); }
static double cpu() {
    rusage r{};
    if (getrusage(RUSAGE_SELF, &r)) throw std::runtime_error("getrusage failed");
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec + (r.ru_utime.tv_usec + r.ru_stime.tv_usec) / 1e6;
}
static void update(EVP_MD_CTX* ctx, const void* p, size_t n) {
    if (EVP_DigestUpdate(ctx, p, n) != 1) throw std::runtime_error("digest update failed");
}
static void field(EVP_MD_CTX* ctx, const char* p, size_t n) {
    unsigned char size[8];
    for (unsigned i = 0; i != 8; ++i) size[i] = uint64_t(n) >> (i * 8);
    update(ctx, size, 8); update(ctx, p, n);
}
static unsigned integer(const std::string& s, unsigned maximum) {
    size_t end = 0;
    if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos) throw std::runtime_error("invalid integer");
    auto n = std::stoul(s, &end);
    if (end != s.size() || n == 0 || n > maximum) throw std::runtime_error("integer out of range");
    return n;
}
struct HtsRecord {
    std::unique_ptr<char, decltype(&free)> sequence{nullptr, free}, quality{nullptr, free};
    size_t length = 0;
    bool found = false;
};


struct PyRelease { void operator()(PyObject* p) const { Py_XDECREF(p); } };
using PyOwned = std::unique_ptr<PyObject, PyRelease>;
static PyOwned checked(PyObject* p) {
    if (!p) { PyErr_Print(); throw std::runtime_error("pyfastx operation failed"); }
    return PyOwned(p);
}
struct PyRecord {
    PyOwned record, sequence, quality;
    const char *seq = "", *qual = "";
    Py_ssize_t n = 0, q = 0;
    bool found = false;
};

int main(int argc, char** argv) try {
    // backend reads fai fastq(0/1) port snapshot
    if (argc != 7) throw std::runtime_error("usage: worker indexed|remote reads fai 0|1 port snapshot");
    std::string backend = argv[1];
    if (backend != "indexed" && backend != "remote" && backend != "pyfastx") throw std::runtime_error("invalid backend");
    if (std::string(argv[4]) != "0" && std::string(argv[4]) != "1") throw std::runtime_error("invalid format");
    const bool fastq = std::string(argv[4]) == "1";
    auto init = Clock::now();
    std::unique_ptr<faidx_t, decltype(&fai_destroy)> index(nullptr, fai_destroy);
    std::unique_ptr<TracEon::TraceonClient> client;
    PyOwned pyindex;
    if (backend == "pyfastx") {
        // Interpreter/module initialization belongs to this backend's setup/PSS.
        Py_Initialize();
        auto module = checked(PyImport_ImportModule("pyfastx"));
        auto type = checked(PyObject_GetAttrString(module.get(), fastq ? "Fastq" : "Fasta"));
        auto args = checked(Py_BuildValue("(s)", argv[2]));
        auto kwargs = checked(fastq ? Py_BuildValue("{s:s}", "index_file", argv[3])
                                   : Py_BuildValue("{s:s,s:O}", "index_file", argv[3], "uppercase", Py_False));
        pyindex = checked(PyObject_Call(type.get(), args.get(), kwargs.get()));
    } else if (backend == "indexed") {
        index.reset(fai_load3_format(argv[2], argv[3], nullptr, 0, fastq ? FAI_FASTQ : FAI_FASTA));
        if (!index) throw std::runtime_error("index load failed");
    } else {
        client = std::make_unique<TracEon::TraceonClient>("127.0.0.1", integer(argv[5], 65535), 30000, 5000, argv[6]);
    }
    std::cout << std::setprecision(17) << "{\"ready_seconds\":" << seconds(Clock::now() - init) << "}\n" << std::flush;
    std::string command;
    while (std::getline(std::cin, command)) {
        std::istringstream in(command);
        std::string word, ids_path, batch_arg, latency_path, extra;
        if (!std::getline(in, word, '\t') || word != "PREP" || !std::getline(in, ids_path, '\t') ||
            !std::getline(in, batch_arg, '\t') || !std::getline(in, latency_path, '\t') || std::getline(in, extra, '\t'))
            throw std::runtime_error("invalid PREP");
        auto batch_size = integer(batch_arg, 256);
        std::ifstream file(ids_path);
        if (!file) throw std::runtime_error("ID input failed");
        std::vector<std::string> ids;
        std::string id;
        while (std::getline(file, id)) {
            if (id.empty() || id.find('\r') != std::string::npos) throw std::runtime_error("invalid ID");
            ids.push_back(id);
        }
        if (!file.eof() || ids.empty()) throw std::runtime_error("invalid ID file");
        std::cout << "PREPARED\n" << std::flush;
        if (!std::getline(std::cin, command) || command != "RUN") throw std::runtime_error("expected RUN");
        std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> digest(EVP_MD_CTX_new(), EVP_MD_CTX_free);
        if (!digest || EVP_DigestInit_ex(digest.get(), EVP_sha256(), nullptr) != 1) throw std::runtime_error("digest init failed");
        std::vector<double> times;
        times.reserve((ids.size() + batch_size - 1) / batch_size);
        uint64_t bytes = 0, found = 0;
        double lookup_seconds = 0, initial_cpu = cpu();
        const auto started = Clock::now();
        for (size_t offset = 0; offset < ids.size(); offset += batch_size) {
            const auto begin = Clock::now();
            const auto count = std::min<size_t>(batch_size, ids.size() - offset);
            std::vector<std::string> keys(ids.begin() + offset, ids.begin() + offset + count);
            std::vector<std::optional<TracEon::DecodedFastqRecord>> remote;
            std::vector<HtsRecord> local;
            std::vector<PyRecord> python;
            if (client) {
                if (batch_size > 1) remote = client->getBatch(keys, fastq);
                else if (fastq) remote.push_back(client->getFastqRecord(keys[0]));
                else {
                    auto seq = client->getView(keys[0]);
                    remote.push_back(seq ? std::optional<TracEon::DecodedFastqRecord>({std::move(*seq), {}}) : std::nullopt);
                }
                if (remote.size() != count) throw std::runtime_error("remote count mismatch");
            } else if (pyindex) {
                python.resize(count);
                for (size_t j = 0; j < count; ++j) {
                    auto& r = python[j];
                    auto key = checked(PyUnicode_FromStringAndSize(keys[j].data(), keys[j].size()));
                    r.record.reset(PyObject_GetItem(pyindex.get(), key.get()));
                    if (!r.record) {
                        if (PyErr_ExceptionMatches(PyExc_KeyError)) { PyErr_Clear(); continue; }
                        PyErr_Print(); throw std::runtime_error("pyfastx lookup failed");
                    }
                    r.sequence = checked(PyObject_GetAttrString(r.record.get(), "seq"));
                    r.seq = PyUnicode_AsUTF8AndSize(r.sequence.get(), &r.n);
                    if (!r.seq) { PyErr_Print(); throw std::runtime_error("sequence conversion failed"); }
                    if (fastq) {
                        r.quality = checked(PyObject_GetAttrString(r.record.get(), "qual"));
                        r.qual = PyUnicode_AsUTF8AndSize(r.quality.get(), &r.q);
                        if (!r.qual || r.n != r.q) throw std::runtime_error("quality conversion failed");
                    }
                    r.found = true;
                }
            } else {
                local.resize(count);
                for (size_t j = 0; j < count; ++j) {
                    auto& r = local[j];
                    hts_pos_t n, actual;
                    r.sequence.reset(faidx_fetch_seq64(index.get(), keys[j].c_str(), 0, HTS_POS_MAX, &n));
                    if (n == -2) continue;
                    if (!r.sequence || n < 0) throw std::runtime_error("sequence retrieval failed");
                    r.found = true; r.length = n;
                    if (fastq) {
                        r.quality.reset(faidx_fetch_qual64(index.get(), keys[j].c_str(), 0, HTS_POS_MAX, &actual));
                        if (!r.quality || actual != n) throw std::runtime_error("quality retrieval failed");
                    }
                }
            }
            const auto verify_start = Clock::now();
            for (size_t j = 0; j < count; ++j) {
                const bool exists = client ? bool(remote[j]) : (pyindex ? python[j].found : local[j].found);
                const unsigned char status = exists ? 1 : 0;
                field(digest.get(), keys[j].data(), keys[j].size());
                update(digest.get(), &status, 1);
                if (!exists) continue;
                ++found;
                const char* seq = client ? remote[j]->sequence.data() : (pyindex ? python[j].seq : local[j].sequence.get());
                const char* qual = client ? remote[j]->quality.data() : (pyindex ? python[j].qual : local[j].quality.get());
                size_t n = client ? remote[j]->sequence.size() : (pyindex ? python[j].n : local[j].length);
                size_t q = fastq ? (client ? remote[j]->quality.size() : (pyindex ? python[j].q : local[j].length)) : 0;
                field(digest.get(), seq ? seq : "", n); field(digest.get(), qual ? qual : "", q);
                bytes += n + q;
            }
            const auto verify_end = Clock::now();
            // Include owned-result destruction; SHA-256 verification is excluded
            // from the lookup timer but included in the primary active window.
            decltype(remote){}.swap(remote);
            decltype(local){}.swap(local);
            decltype(python){}.swap(python);
            decltype(keys){}.swap(keys);
            double elapsed = seconds(Clock::now() - begin) - seconds(verify_end - verify_start);
            lookup_seconds += elapsed;
            times.push_back(elapsed * 1e6);
        }
        const auto ended = Clock::now();
        const double used_cpu = cpu() - initial_cpu;
        unsigned char hash[EVP_MAX_MD_SIZE]; unsigned n = 0;
        if (EVP_DigestFinal_ex(digest.get(), hash, &n) != 1 || n != 32) throw std::runtime_error("digest final failed");
        std::ostringstream hex;
        for (unsigned i = 0; i < n; ++i) hex << std::hex << std::setfill('0') << std::setw(2) << unsigned(hash[i]);
        std::ofstream latency(latency_path);
        for (auto t : times) latency << std::setprecision(17) << t << '\n';
        latency.close();
        if (!latency) throw std::runtime_error("latency output failed");
        std::cout << "{\"start\":" << stamp(started) << ",\"end\":" << stamp(ended)
                  << ",\"lookup_seconds\":" << lookup_seconds << ",\"cpu_seconds\":" << used_cpu
                  << ",\"records\":" << ids.size() << ",\"found\":" << found << ",\"bytes\":" << bytes
                  << ",\"requests\":" << times.size() << ",\"sha256\":\"" << hex.str() << "\"}\n" << std::flush;
    }
    if (!std::cin.eof()) throw std::runtime_error("stdin failed");
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
