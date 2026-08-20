// benchmarks/remote_bench.cpp — latency/throughput benchmark for TracEon
// remote read access (v2.3.0 "Harpe").
//
// Three modes:
//
//   --mode=local  In-process getView() lookups against a locally loaded v4
//                 cache. This is the zero-copy baseline (no network).
//
//   --mode=remote TCP client-server lookups against a TraceonServer. The
//                 cache is loaded locally ONLY to enumerate the real key set
//                 (not timed); every measured lookup travels over the socket
//                 and is CRC32C-verified on the client. One connection per
//                 thread (the recommended deployment pattern).
//
//   --mode=serve  Load a v4 cache and serve it over TCP (this process embeds
//                 the server; used by the Docker testbed and manual runs).
//
// Usage:
//   remote_bench --mode=local --file=cache.traceon [--keys=N --count=N --warmup=N --seed=S]
//   remote_bench --mode=remote --host=H --port=P --file=cache.traceon [--threads=T --keys=N --count=N --warmup=N --seed=S]
//   remote_bench --mode=serve --file=cache.traceon [--host=0.0.0.0 --port=9876]
//
// Output: a single tab-separated summary line (easy to grep) plus a human
// readable table with median/p95/p99 latency in microseconds and aggregate
// ops/s. Latencies are per-request round-trips measured with steady_clock;
// the median is taken over all measured requests; p95/p99 are the tail.

#define _GNU_SOURCE 1

#include "Cache.h"
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
        else if (arg == "--port" || arg == "-p") o.port = static_cast<uint16_t>(std::atoi(val().c_str()));
        else if (arg == "--file" || arg == "-f") o.file = val();
        else if (arg == "--keys" || arg == "-k") o.keys = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--count" || arg == "-c") o.count = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--warmup") o.warmup = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--threads" || arg == "-t") o.threads = std::atoi(val().c_str());
        else if (arg == "--seed") o.seed = std::strtoull(val().c_str(), nullptr, 10);
        else {
            std::fprintf(stderr, "unknown option: %s\n", arg.c_str());
            std::exit(2);
        }
    }
    if (o.file.empty()) {
        std::fprintf(stderr, "remote_bench: --file=<cache.traceon> is required\n");
        std::exit(2);
    }
    if (o.mode != "local" && o.mode != "remote" && o.mode != "serve") {
        std::fprintf(stderr, "remote_bench: --mode must be local|remote|serve\n");
        std::exit(2);
    }
    if (o.mode == "remote" && o.threads < 1) o.threads = 1;
    return o;
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

LatencyStats summarize(std::vector<double>& times_us) {
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
    double sum = 0;
    for (double t : times_us) sum += t;
    s.total_ms = sum / 1000.0;
    s.ops = static_cast<double>(times_us.size()) / (s.total_ms / 1000.0);
    return s;
}

void report(const std::string& mode, const Options& o, const std::vector<double>& times_us) {
    const LatencyStats s = summarize(const_cast<std::vector<double>&>(times_us));
    std::printf(
        "TRACEON_BENCH\tmode=%s\tkeys=%llu\tcount=%llu\twarmup=%llu\tthreads=%d\t"
        "median_us=%.2f\tp95_us=%.2f\tp99_us=%.2f\tops=%.0f\n",
        mode.c_str(),
        (unsigned long long)o.keys, (unsigned long long)o.count,
        (unsigned long long)o.warmup, o.threads,
        s.median_us, s.p95_us, s.p99_us, s.ops);
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
    std::vector<double> times = bench(o.count);
    report("local", o, times);
    return 0;
}

// ── remote mode ──────────────────────────────────────────────────────────────
int runRemote(const Options& o) {
    // Load once here only to enumerate keys; the load itself is untimed.
    std::vector<std::string> keys = loadKeys(o.file);
    if (keys.empty()) { std::fprintf(stderr, "cache is empty\n"); return 1; }

    std::mt19937_64 rng(o.seed);
    const size_t nkeys = std::min<size_t>(keys.size(), static_cast<size_t>(o.keys));
    std::uniform_int_distribution<size_t> dist(0, keys.size() - 1);
    std::vector<std::string> sample(nkeys);
    for (size_t i = 0; i < nkeys; ++i) sample[i] = keys[dist(rng)];

    const uint64_t per_thread_count = (o.count + static_cast<uint64_t>(o.threads) - 1)
                                      / static_cast<uint64_t>(o.threads);
    std::atomic<size_t> failures{0};
    std::vector<double> all_times;
    std::mutex times_mu;
    all_times.reserve(static_cast<size_t>(per_thread_count * static_cast<uint64_t>(o.threads)));

    auto worker = [&](int tid) {
        // One connection per thread (documented deployment pattern).
        TracEon::TraceonClient client(o.host, o.port);
        std::vector<double> local_times;
        local_times.reserve(static_cast<size_t>(per_thread_count));

        // Warmup on this connection.
        for (uint64_t i = 0; i < o.warmup; ++i) {
            const std::string& k = sample[(i + static_cast<uint64_t>(tid)) % sample.size()];
            auto v = client.getView(k);
            if (!v) failures.fetch_add(1, std::memory_order_relaxed);
        }
        for (uint64_t i = 0; i < per_thread_count; ++i) {
            const std::string& k = sample[(i * static_cast<uint64_t>(o.threads) + static_cast<uint64_t>(tid)) % sample.size()];
            auto t0 = Clock::now();
            auto v = client.getView(k);
            auto t1 = Clock::now();
            local_times.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
            if (!v) failures.fetch_add(1, std::memory_order_relaxed);
        }
        std::lock_guard<std::mutex> lock(times_mu);
        all_times.insert(all_times.end(), local_times.begin(), local_times.end());
    };

    std::vector<std::thread> threads;
    for (int t = 0; t < o.threads; ++t) threads.emplace_back(worker, t);
    for (auto& t : threads) t.join();

    if (failures.load() != 0) {
        std::fprintf(stderr, "FAIL: %zu remote lookups returned NOT_FOUND\n", failures.load());
        return 1;
    }
    report("remote", o, all_times);
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
    Options o = parseArgs(argc, argv);
    if (o.mode == "local") return runLocal(o);
    if (o.mode == "remote") return runRemote(o);
    return runServe(o);
}