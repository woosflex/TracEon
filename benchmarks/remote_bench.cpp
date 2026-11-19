// benchmarks/remote_bench.cpp — latency/throughput benchmark for TracEon
// remote read access (v2.3.0 "Harpe").
//
// Three modes:
//
//   --mode=local  In-process getView() lookups against a locally loaded v4
//                 cache. This is the zero-copy baseline (no network).
//
//   --mode=remote TCP lookups using --ids (no local cache) or legacy --file
//                 enumeration. --expect-snapshot pins the dataset; --fastq
//                 retrieves sequence and quality. Timings include owned-result
//                 destruction and CRC verification, but exclude handshake/warmup.
//   --mode=codec  Response allocation/copy/CRC costs without socket I/O.
//
//   --mode=serve  Load a v4 cache and serve it over TCP (this process embeds
//                 the server; used by the Docker testbed and manual runs).
//
// Usage:
//   remote_bench --mode=local --file=cache.traceon [--keys=N --count=N --warmup=N --seed=S]
//   remote_bench --mode=remote --host=H --port=P --file=cache.traceon [--threads=T --keys=N --count=N --warmup=N --seed=S]
//   remote_bench --mode=serve --file=cache.traceon [--host=0.0.0.0 --port=9876]
//   remote_bench --mode=remote --ids=ids.txt --expect-snapshot=SHA256 --latencies=request-us.txt
//   remote_bench --mode=codec --bytes=10000 --count=10000 --latencies=codec-us.tsv
//
// Output: a single tab-separated summary line (easy to grep) plus a human
// readable table with median/p95/p99 latency in microseconds and aggregate
// ops/s. Latencies are per-request round-trips measured with steady_clock;
// the median is taken over all measured requests; p95/p99 are the tail.
// With --batch-size=2..256, count/warmup/ops/s remain REQUEST units. Each
// latency includes batch-key construction, transport, decoding and destruction.

#define _GNU_SOURCE 1

#include "Cache.h"
#include "BenchmarkUtils.h"
#include <limits>
#include <exception>
#include "TraceonClient.h"
#include "TraceonServer.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    std::string mode = "local";
    std::string host = "127.0.0.1";
    uint16_t port = 9876;
    std::string file;
    std::string ids, snapshot, latencies;
    bool fastq = false;
    uint64_t batch_size = 1;
    uint64_t bytes = 1000;
    uint64_t keys = 10000;
    uint64_t count = 10000;
    uint64_t warmup = 1000;
    int threads = 1;
    uint64_t seed = 42;
};

Options parseArgs(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        // Accept both "--opt=value" and "--opt value" forms (== is the
        // documented style and the one used by the Docker testbed).
        std::string arg = argv[i];
        std::string inline_val;
        const std::size_t eq = arg.find('=');
        const bool has_inline = (eq != std::string::npos);
        if (has_inline) {
            inline_val = arg.substr(eq + 1);
            arg = arg.substr(0, eq);
        }
        auto val = [&]() -> std::string {
            if (has_inline) return inline_val;
            if (i + 1 >= argc) {
                std::fprintf(stderr, "missing value for %s\n", arg.c_str());
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--mode" || arg == "-m") o.mode = val();
        else if (arg == "--host") o.host = val();
        else if (arg == "--port" || arg == "-p") {
            const auto n = benchmarkUnsigned(val(), arg);
            if (n > 65535) throw std::invalid_argument("--port must be <= 65535");
            o.port = static_cast<uint16_t>(n);
        }
        else if (arg == "--file" || arg == "-f") o.file = val();
        else if (arg == "--ids") o.ids = val();
        else if (arg == "--expect-snapshot") {
            o.snapshot = val();
            if (!TracEon::proto::validSnapshotId(o.snapshot))
                throw std::invalid_argument("expected a lowercase snapshot SHA-256");
        }
        else if (arg == "--latencies") o.latencies = val();
        else if (arg == "--fastq" && !has_inline) o.fastq = true;
        else if (arg == "--batch-size") o.batch_size = benchmarkUnsigned(val(), arg);
        else if (arg == "--bytes") o.bytes = benchmarkUnsigned(val(), arg);
        else if (arg == "--keys" || arg == "-k") o.keys = benchmarkUnsigned(val(), arg);
        else if (arg == "--count" || arg == "-c") o.count = benchmarkUnsigned(val(), arg);
        else if (arg == "--warmup") o.warmup = benchmarkUnsigned(val(), arg);
        else if (arg == "--threads" || arg == "-t") {
            const auto n = benchmarkUnsigned(val(), arg);
            if (n == 0 || n > static_cast<uint64_t>(std::numeric_limits<int>::max()))
                throw std::invalid_argument("--threads must be a positive int");
            o.threads = static_cast<int>(n);
        }
        else if (arg == "--seed") o.seed = benchmarkUnsigned(val(), arg);
        else {
            std::fprintf(stderr, "unknown option: %s\n", arg.c_str());
            std::exit(2);
        }
    }
    if (o.mode != "local" && o.mode != "remote" && o.mode != "serve" && o.mode != "codec") {
        std::fprintf(stderr, "remote_bench: --mode must be local|remote|serve|codec\n");
        std::exit(2);
    }
    if (o.mode != "codec" && o.file.empty() && (o.mode != "remote" || o.ids.empty()))
        throw std::invalid_argument("--file is required, or --ids for remote mode");
    if (!o.ids.empty() && (o.mode != "remote" || !o.file.empty()))
        throw std::invalid_argument("--ids is remote-only and cannot be combined with --file");
    if ((o.fastq || !o.snapshot.empty()) && o.mode != "remote")
        throw std::invalid_argument("--fastq and --expect-snapshot are remote-only");
    if (o.mode == "serve" && !o.latencies.empty())
        throw std::invalid_argument("--latencies is not available in serve mode");
    if (o.mode == "codec" && (o.bytes > TracEon::proto::kMaxMessageSize - 10 ||
                            o.threads != 1 || !o.file.empty()))
        throw std::invalid_argument("codec requires one thread, no cache, and a payload within the frame limit");
    if (o.batch_size == 0 || o.batch_size > 256 || (o.batch_size > 1 && o.mode != "remote"))
        throw std::invalid_argument("batch size must be 1–256; batching requires remote mode");
    if (o.keys == 0 || o.count == 0)
        throw std::invalid_argument("--keys and --count must be positive");
    return o;
}

std::vector<std::string> readIds(const std::string& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot open ID file: " + file);
    std::vector<std::string> ids;
    std::string id;
    while (std::getline(input, id)) {
        if (!id.empty() && id.back() == '\r') id.pop_back();
        if (id.empty() || id.find_first_of(" \t\r\n") != id.npos || id.find('\0') != id.npos)
            throw std::invalid_argument("invalid ID file entry");
        ids.push_back(id);
    }
    if (input.bad() || ids.empty()) throw std::runtime_error("empty or unreadable ID file");
    return ids;
}

std::vector<std::string> loadKeys(const std::string& file) {
    // Loading the cache on the client exists only to enumerate the REAL key
    // set (the server serves a cache with exactly these keys). The load time
    // is never part of any measured window.
    TracEon::Cache cache;
    cache.restore(file);
    return cache.getAllKeys();
}

struct LatencyStats {
    double median_us = 0, p95_us = 0, p99_us = 0, ops = 0, total_ms = 0;
};

LatencyStats summarize(std::vector<double>& times_us, double wall_s) {
    std::sort(times_us.begin(), times_us.end());
    LatencyStats s;
    if (times_us.empty()) return s;
    auto at = [&](double p) -> double {
        size_t idx = static_cast<size_t>(p * (times_us.size() - 1));
        return times_us[idx];
    };
    s.median_us = at(0.50);
    s.p95_us = at(0.95);
    s.p99_us = at(0.99);
    s.total_ms = wall_s * 1000.0;
    s.ops = wall_s > 0 ? static_cast<double>(times_us.size()) / wall_s : 0;
    return s;
}

void report(const std::string& mode, const Options& o, std::vector<double>& times_us, double wall_s) {
    if (!o.latencies.empty()) {
        std::ofstream out(o.latencies);
        out << std::setprecision(17);
        for (double value : times_us) out << value << '\n';
        out.close();
        if (!out) throw std::runtime_error("cannot write latency samples: " + o.latencies);
    }
    const LatencyStats s = summarize(times_us, wall_s);
    std::printf(
        "TRACEON_BENCH\tmode=%s\tkeys=%llu\tcount=%llu\twarmup=%llu\tthreads=%d\t"
        "median_us=%.2f\tp95_us=%.2f\tp99_us=%.2f\twall_s=%.9f\tops=%.0f\n",
        mode.c_str(),
        (unsigned long long)o.keys, (unsigned long long)o.count,
        (unsigned long long)o.warmup, o.threads,
        s.median_us, s.p95_us, s.p99_us, wall_s, s.ops);
}

// Pure response path: allocation/copy + CRC, without sockets, lookup, or scheduling.
int runCodec(const Options& o) {
    namespace proto = TracEon::proto;
    const std::vector<char> payload(static_cast<size_t>(o.bytes), 'A');
    std::vector<double> encode, frame, receive_copy, verify, owned;
    auto elapsed = [](auto a, auto b) { return std::chrono::duration<double, std::micro>(b - a).count(); };
    auto one = [&](bool record) {
        const auto t0 = Clock::now();
        auto body = proto::encodeOk(proto::kOk, payload);
        const auto t1 = Clock::now();
        auto wire = proto::encodeFrame(proto::MsgType::Ok, body);
        const auto t2 = Clock::now();
        // readFrame's body copy, after the transport has supplied a complete frame.
        std::vector<char> received(wire.begin() + proto::kFrameHeaderSize + 1, wire.end());
        const auto t3 = Clock::now();
        auto [status, bytes] = proto::decodeOk(received);
        const auto t4 = Clock::now();
        const std::string result(bytes);
        const auto t5 = Clock::now();
        if (status != proto::kOk || result.size() != payload.size() ||
            !std::equal(result.begin(), result.end(), payload.begin()))
            throw std::runtime_error("codec payload mismatch");
        if (record) {
            encode.push_back(elapsed(t0, t1)); frame.push_back(elapsed(t1, t2));
            receive_copy.push_back(elapsed(t2, t3)); verify.push_back(elapsed(t3, t4));
            owned.push_back(elapsed(t4, t5));
        }
    };
    for (auto* values : {&encode, &frame, &receive_copy, &verify, &owned}) values->reserve(o.count);
    for (uint64_t i = 0; i < o.warmup; ++i) one(false);
    for (uint64_t i = 0; i < o.count; ++i) one(true);
    if (!o.latencies.empty()) {
        std::ofstream out(o.latencies);
        out << std::setprecision(17)
            << "encode_crc_us\tframe_copy_us\treceive_copy_us\tverify_crc_us\towned_copy_us\n";
        for (size_t i = 0; i < encode.size(); ++i)
            out << encode[i] << '\t' << frame[i] << '\t' << receive_copy[i] << '\t'
                << verify[i] << '\t' << owned[i] << '\n';
        out.close();
        if (!out) throw std::runtime_error("cannot write codec samples");
    }
    auto median = [](auto& values) { return summarize(values, 0).median_us; };
    std::cout << std::setprecision(9) << "TRACEON_CODEC\tbytes=" << o.bytes << "\tcount=" << o.count
              << "\tencode_crc_us=" << median(encode) << "\tframe_copy_us=" << median(frame)
              << "\treceive_copy_us=" << median(receive_copy) << "\tverify_crc_us=" << median(verify)
              << "\towned_copy_us=" << median(owned) << '\n';
    return 0;
}

// ── local mode ───────────────────────────────────────────────────────────────
int runLocal(const Options& o) {
    TracEon::Cache cache;
    cache.restore(o.file);
    std::vector<std::string> keys = cache.getAllKeys();
    if (keys.empty()) { std::fprintf(stderr, "cache is empty\n"); return 1; }

    std::mt19937_64 rng(o.seed);
    const size_t nkeys = std::min<size_t>(keys.size(), static_cast<size_t>(o.keys));
    std::uniform_int_distribution<size_t> dist(0, keys.size() - 1);
    std::vector<std::string> sample(nkeys);
    for (size_t i = 0; i < nkeys; ++i) sample[i] = keys[dist(rng)];
    if (sample.empty()) { std::fprintf(stderr, "no sampled keys\n"); return 1; }

    auto bench = [&](uint64_t n) {
        std::vector<double> times(n);
        volatile size_t sink = 0;
        for (uint64_t i = 0; i < n; ++i) {
            const std::string& k = sample[static_cast<size_t>(i % sample.size())];
            auto t0 = Clock::now();
            std::string_view v = cache.getView(k);
            auto t1 = Clock::now();
            sink += v.size();
            times[static_cast<size_t>(i)] = std::chrono::duration<double, std::micro>(t1 - t0).count();
        }
        (void)sink;
        return times;
    };

    (void)bench(o.warmup); // warmup
    const auto start = Clock::now();
    std::vector<double> times = bench(o.count);
    report("local", o, times, std::chrono::duration<double>(Clock::now() - start).count());
    return 0;
}

// ── remote mode ──────────────────────────────────────────────────────────────
int runRemote(const Options& o) {
    // Load once here only to enumerate keys; the load itself is untimed.
    std::vector<std::string> keys = o.ids.empty() ? loadKeys(o.file) : readIds(o.ids);
    if (keys.empty()) { std::fprintf(stderr, "cache is empty\n"); return 1; }

    std::mt19937_64 rng(o.seed);
    const size_t nkeys = std::min<size_t>(keys.size(), static_cast<size_t>(o.keys));
    std::uniform_int_distribution<size_t> dist(0, keys.size() - 1);
    std::vector<std::string> sample(nkeys);
    for (size_t i = 0; i < nkeys; ++i) sample[i] = keys[dist(rng)];

    if (sample.empty()) throw std::invalid_argument("no sampled keys");
    std::atomic<size_t> failures{0};
    std::atomic<bool> stop{false};
    std::atomic<int> ready{0};
    std::atomic<bool> go{false};
    std::vector<std::vector<double>> times(o.threads);
    std::vector<Clock::time_point> finished(o.threads);
    std::vector<std::exception_ptr> errors(o.threads);
    std::vector<uint64_t> payload_bytes(o.threads);
    std::vector<std::thread> threads;
    for (int tid = 0; tid < o.threads; ++tid) {
        threads.emplace_back([&, tid] {
            bool announced = false;
            try {
                TracEon::TraceonClient client(o.host, o.port, TracEon::proto::kReadTimeoutMs, 5000, o.snapshot);
                auto query = [&](size_t index) -> uint64_t {
                    if (o.batch_size > 1) {
                        std::vector<std::string> batch;
                        batch.reserve(o.batch_size);
                        for (size_t j = 0; j < o.batch_size; ++j)
                            batch.push_back(sample[(index + j) % sample.size()]);
                        auto values = client.getBatch(batch, o.fastq);
                        uint64_t bytes = 0;
                        for (const auto& value : values) {
                            if (!value) ++failures;
                            else bytes += value->sequence.size() + value->quality.size();
                        }
                        return bytes;
                    }
                    const auto& key = sample[index];
                    if (o.fastq) {
                        auto value = client.getFastqRecord(key);
                        if (!value) { ++failures; return 0; }
                        return value->sequence.size() + value->quality.size();
                    }
                    auto value = client.getView(key);
                    if (!value) { ++failures; return 0; }
                    return value->size();
                };
                auto& local_times = times[tid];
                const uint64_t count = o.count / o.threads + (static_cast<uint64_t>(tid) < o.count % o.threads);
                local_times.reserve(static_cast<size_t>(count));
                for (uint64_t i = 0; i < o.warmup; ++i) {
                    (void)query((i + tid) % sample.size());
                }
                announced = true;
                ready.fetch_add(1, std::memory_order_release);
                while (!go.load(std::memory_order_acquire)) {
                    if (stop.load(std::memory_order_acquire)) return;
                    std::this_thread::yield();
                }
                for (uint64_t i = 0; i < count; ++i) {
                    const auto index = (i * o.threads + tid) % sample.size();
                    const auto t0 = Clock::now();
                    payload_bytes[tid] += query(index);
                    const auto t1 = Clock::now();
                    local_times.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
                }
                finished[tid] = Clock::now();
            } catch (...) {
                errors[tid] = std::current_exception();
                if (!announced) ready.fetch_add(1, std::memory_order_release);
            }
        });
    }
    struct ThreadsGuard {
        std::vector<std::thread>& ths;
        std::atomic<bool>& s;
        ~ThreadsGuard() {
            s.store(true, std::memory_order_release);
            for (auto& th : ths) {
                if (th.joinable()) th.join();
            }
        }
    } guard{threads, stop};
    while (ready.load(std::memory_order_acquire) != o.threads) std::this_thread::yield();
    const auto start = Clock::now();
    go.store(true, std::memory_order_release);
    for (auto& thread : threads) thread.join();
    for (const auto& error : errors) if (error) std::rethrow_exception(error);
    const double wall_s = std::chrono::duration<double>(
        *std::max_element(finished.begin(), finished.end()) - start).count();
    std::vector<double> all_times;
    all_times.reserve(static_cast<size_t>(o.count));
    for (auto& local_times : times)
        all_times.insert(all_times.end(), local_times.begin(), local_times.end());

    if (failures.load() != 0) {
        std::fprintf(stderr, "FAIL: %zu remote lookups returned NOT_FOUND\n", failures.load());
        return 1;
    }
    report("remote", o, all_times, wall_s);
    std::cout << "TRACEON_REQUESTS\tbatch_size=" << o.batch_size << "\tlatency_unit=request\n";
    uint64_t total_bytes = 0;
    for (auto bytes : payload_bytes) total_bytes += bytes;
    std::cout << "TRACEON_PAYLOAD\tbytes=" << total_bytes << "\tfastq=" << o.fastq << '\n';
    return 0;
}

// ── serve mode ───────────────────────────────────────────────────────────────
std::atomic<bool> g_stop{false};

void handleSignal(int) { g_stop.store(true); }

int runServe(const Options& o) {
    TracEon::Cache cache;
    cache.restore(o.file);
    std::fprintf(stderr, "loaded %s: %zu records\n", o.file.c_str(), cache.size());

    TracEon::TraceonServer server(cache, o.port, o.host);
    server.start();
    std::printf("TRACEON_SERVER READY host=%s port=%u records=%zu\n",
                o.host.c_str(), server.port(), cache.size());
    std::fflush(stdout);

    std::signal(SIGINT, handleSignal);
    std::signal(SIGTERM, handleSignal);
    while (!g_stop.load() && server.isRunning()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    server.stop();
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options o = parseArgs(argc, argv);
        if (o.mode == "local") return runLocal(o);
        if (o.mode == "remote") return runRemote(o);
        if (o.mode == "codec") return runCodec(o);
        return runServe(o);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "remote_bench: %s\n", e.what());
        return 1;
    }
}
