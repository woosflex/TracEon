// benchmarks/read_scaling.cpp — lock-free getView() read throughput vs thread
// count (v2.3.0 "Harpe", slice 2: multi-core read scaling).
//
// Measures aggregate ops/s and the scaling factor for concurrent ZERO-COPY
// getView() lookups against ONE in-memory immutable cache:
//
//   * The cache is built in-process via SmartStrategy::addEntry() — the exact
//     build path the unit tests use (id/seq/qual stored in manual_store_,
//     non-owning string_view keys, data_ready_ published at the end). The
//     build phase is NOT part of any measured window.
//   * RECORDS distinct records are added; the read set is KEYS (default: all)
//     of them, shuffled per thread with per-thread seeds so threads do NOT
//     walk the index in lockstep (no artificial line-collision phase lock).
//   * Each worker thread loops COUNT times over its shuffled key set. All
//     workers are released simultaneously through a start barrier, then each
//     measures its own wall time. Aggregate ops/s = total lookups / max
//     per-thread elapsed (the coordinator's wall window from barrier release
//     to the last thread's finish).
//
// Honest framing for this host (Intel Core Ultra 5 125H, 14 cores / 18
// threads, single NUMA node, L3 = 18 MiB): with the default 100,000 records
// the read working set (index bucket array ~8 MB + key strings ~3 MB) fits
// in L3, so lookups are L3-resident after warmup. This measures the scaling
// of the lock-free read path (one acquire load of data_ready_ + an ankerl
// hash probe + zero-copy return), NOT DRAM bandwidth. Nothing is pinned:
// thread placement is left to the OS scheduler (P/E-core hybrid), reported
// as-is.
//
// Usage:
//   read_scaling [--records=N] [--keys=N] [--count=N] [--warmup=N]
//                [--threads=1,2,4,8,14] [--reps=R] [--seed=S]
//
// Output: one TRACEON_BENCH line per (thread-count, repetition) — easy to
// grep — plus a readable summary table (threads -> median ops/s -> scaling
// vs 1 thread -> median per-thread latency -> hit rate).

#include "SmartStrategy.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Options {
    uint64_t records = 100000;   // records added to the in-memory cache
    uint64_t keys = 0;           // read-set size (0 = all records)
    uint64_t count = 20;         // passes over the key set per thread
    uint64_t warmup = 2;         // untimed warmup passes per thread
    std::vector<int> threads = {1, 2, 4, 8, 14};
    uint64_t reps = 3;           // repetitions per thread count (min/median)
    uint64_t seed = 42;
};

Options parseArgs(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string inline_val;
        const std::size_t eq = arg.find('=');
        if (eq != std::string::npos) {
            inline_val = arg.substr(eq + 1);
            arg = arg.substr(0, eq);
        }
        auto val = [&]() -> std::string {
            if (!inline_val.empty()) return inline_val;
            if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", arg.c_str()); std::exit(2); }
            return std::string(argv[++i]);
        };
        if (arg == "--records") o.records = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--keys") o.keys = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--count") o.count = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--warmup") o.warmup = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--reps") o.reps = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--seed") o.seed = std::strtoull(val().c_str(), nullptr, 10);
        else if (arg == "--threads") {
            o.threads.clear();
            const std::string s = val();
            size_t pos = 0;
            while (pos < s.size()) {
                size_t comma = s.find(',', pos);
                o.threads.push_back(std::atoi(s.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos).c_str()));
                if (comma == std::string::npos) break;
                pos = comma + 1;
            }
            if (o.threads.empty() || o.threads[0] < 1) { std::fprintf(stderr, "--threads needs >= 1 positive values\n"); std::exit(2); }
        }
        else { std::fprintf(stderr, "unknown option: %s\n", arg.c_str()); std::exit(2); }
    }
    return o;
}

// Deterministic pseudo-genomic data: id "seq_<i>", SEQ_LEN random ACGT bases,
// QUAL_LEN 'I' quality chars (Phred-high, like a real Illumina read).
constexpr size_t SEQ_LEN = 300;
constexpr size_t QUAL_LEN = 300;

std::string makeSeq(std::mt19937_64& rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    std::string s(SEQ_LEN, 'A');
    for (char& c : s) c = bases[rng() & 3u];
    return s;
}

// ── parallel harness ─────────────────────────────────────────────────────────
struct WorkerStats {
    double elapsed_s = 0.0;
    uint64_t lookups = 0;
    uint64_t hits = 0;
    uint64_t sink = 0; // checksum of v.size() — proves real reads, blocks DCE
};

// One pass over the key set for this thread (used for both warmup and timing).
inline WorkerStats runPass(const TracEon::SmartStrategy& cache,
                           const std::vector<std::string>& keys,
                           const std::vector<size_t>& order,
                           uint64_t passes) {
    WorkerStats s;
    s.lookups = static_cast<uint64_t>(order.size()) * passes;
    for (uint64_t p = 0; p < passes; ++p) {
        for (size_t i = 0; i < order.size(); ++i) {
            const std::string_view k(keys[order[i]]);
            const std::string_view v = cache.getView(k);
            s.sink += static_cast<uint64_t>(v.size());
            s.hits += static_cast<uint64_t>(!v.empty());
        }
    }
    return s;
}

std::vector<WorkerStats> runThreads(int nthreads,
                                    const TracEon::SmartStrategy& cache,
                                    const std::vector<std::string>& keys,
                                    uint64_t count, uint64_t warmup,
                                    uint64_t seed) {
    std::atomic<bool> go{false};
    std::vector<WorkerStats> stats(nthreads);
    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(nthreads));

    for (int t = 0; t < nthreads; ++t) {
        pool.emplace_back([&, t]() {
            // Private shuffled traversal order: per-thread RNG so threads do
            // not walk the index in lockstep (seed = base seed + t).
            std::vector<size_t> order(keys.size());
            for (size_t i = 0; i < order.size(); ++i) order[i] = i;
            std::mt19937_64 rng(seed + static_cast<uint64_t>(t));
            std::shuffle(order.begin(), order.end(), rng);

            // Untimed warmup: touch every key, populate TLB/caches.
            for (uint64_t w = 0; w < warmup; ++w) (void)runPass(cache, keys, order, 1);

            // Spin on the start barrier so all threads begin within nanoseconds
            // of each other (no coordinator drift in the measured window).
            while (!go.load(std::memory_order_acquire)) { /* spin */ }

            const auto t0 = Clock::now();
            WorkerStats s = runPass(cache, keys, order, count);
            const auto t1 = Clock::now();
            s.elapsed_s = std::chrono::duration<double>(t1 - t0).count();
            stats[static_cast<size_t>(t)] = s;
        });
    }
    go.store(true, std::memory_order_release);
    for (auto& th : pool) th.join();
    return stats;
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

int runBenchmark(const Options& o) {
    // ── Build phase (not timed) ────────────────────────────────────────────────
    std::fprintf(stderr, "read_scaling: building %llu-record in-memory cache (addEntry path) ...\n",
                 (unsigned long long)o.records);
    TracEon::SmartStrategy cache(TracEon::IndexMode::GENOME);
    {
        std::mt19937_64 rng(o.seed);
        std::string seq_buf(SEQ_LEN, 'A'), qual_buf(QUAL_LEN, 'I');
        for (uint64_t i = 0; i < o.records; ++i) {
            std::string id = "seq_" + std::to_string(i);
            seq_buf = makeSeq(rng);
            cache.addEntry(id, seq_buf, qual_buf);
        }
    }
    const uint64_t nkeys = o.keys == 0 ? o.records : std::min(o.keys, o.records);
    std::vector<std::string> keys;
    keys.reserve(static_cast<size_t>(nkeys));
    for (uint64_t i = 0; i < nkeys; ++i) keys.emplace_back("seq_" + std::to_string(i));

    // Sanity: every read-set key is really present (distinct, non-empty reads).
    {
        uint64_t missing = 0;
        for (const auto& k : keys) if (cache.getView(std::string_view(k)).empty()) ++missing;
        if (missing != 0) {
            std::fprintf(stderr, "FAIL: %llu/%llu read-set keys missing from cache\n",
                         (unsigned long long)missing, (unsigned long long)keys.size());
            return 1;
        }
    }
    std::fprintf(stderr, "read_scaling: cache ready — %llu records, read set = %llu distinct keys, "
                 "working set (index+keys) ~ %.1f MiB (host L3 = 18 MiB)\n",
                 (unsigned long long)o.records, (unsigned long long)nkeys,
                 (nkeys * 64.0 / 0.8 + nkeys * 32.0) / (1024.0 * 1024.0));

    // ── Measurement ───────────────────────────────────────────────────────────
    std::printf("TRACEON_BENCH\tbench=read_scaling\trecords=%llu\tkeys=%llu\tcount=%llu\twarmup=%llu\treps=%llu\tseed=%llu\n",
                (unsigned long long)o.records, (unsigned long long)nkeys,
                (unsigned long long)o.count, (unsigned long long)o.warmup,
                (unsigned long long)o.reps, (unsigned long long)o.seed);

    struct PerThreadCount {
        std::vector<double> ops;        // one per rep
        std::vector<double> med_thread_us;
        std::vector<double> hit_rate;
        std::vector<uint64_t> sink;
    };
    std::vector<PerThreadCount> results(o.threads.size());

    for (size_t ti = 0; ti < o.threads.size(); ++ti) {
        const int nthreads = o.threads[ti];
        for (uint64_t r = 0; r < o.reps; ++r) {
            std::vector<WorkerStats> st = runThreads(nthreads, cache, keys, o.count, o.warmup, o.seed);
            const uint64_t total_lookups = o.count * nkeys * static_cast<uint64_t>(nthreads);
            const uint64_t total_hits = [&]{ uint64_t h = 0; for (auto& s : st) h += s.hits; return h; }();
            const uint64_t total_sink = [&]{ uint64_t s = 0; for (auto& w : st) s += w.sink; return s; }();

            double wall = 0.0;
            std::vector<double> per_thread_us;
            for (auto& s : st) {
                wall = std::max(wall, s.elapsed_s);
                per_thread_us.push_back((s.lookups > 0)
                    ? s.elapsed_s / static_cast<double>(s.lookups) * 1e6 : 0.0);
            }
            const double ops = static_cast<double>(total_lookups) / wall;
            const double med_us = median(per_thread_us);
            const double hit = static_cast<double>(total_hits) / static_cast<double>(total_lookups);
            results[ti].ops.push_back(ops);
            results[ti].med_thread_us.push_back(med_us);
            results[ti].hit_rate.push_back(hit);
            results[ti].sink.push_back(total_sink);

            std::printf("TRACEON_BENCH\tbench=read_scaling\tthreads=%d\trep=%llu\tops=%.0f\treads_per_us=%.3f\tmed_thread_us=%.3f\thit_rate=%.4f\tsink=%llu\n",
                        nthreads, (unsigned long long)(r + 1), ops, ops / 1e6, med_us, hit,
                        (unsigned long long)total_sink);
            std::fflush(stdout);
        }
    }

    // ── Summary table (median across reps; scaling vs a real 1-thread baseline) ──
    // The scaling factor must ALWAYS be anchored to a measured 1-thread run. If
    // the user's thread list didn't include 1, run a single extra 1-thread rep
    // purely to serve as the baseline (not printed as a config row).
    double base = 0.0;
    for (size_t ti = 0; ti < o.threads.size(); ++ti) {
        if (o.threads[ti] == 1) { base = median(results[ti].ops); break; }
    }
    bool base_was_extra = false;
    if (base <= 0.0) {
        std::vector<WorkerStats> st = runThreads(1, cache, keys, o.count, o.warmup, o.seed);
        const uint64_t tl = o.count * nkeys;
        double wall = 0.0;
        for (auto& s : st) wall = std::max(wall, s.elapsed_s);
        base = static_cast<double>(tl) / wall;
        base_was_extra = true;
    }
    std::printf("\nread_scaling summary — thread count -> ops/s (median [min]) -> scaling vs 1 thread -> median per-thread latency (us/op) -> hit%%\n");
    for (size_t ti = 0; ti < o.threads.size(); ++ti) {
        const int n = o.threads[ti];
        const double med_ops = median(results[ti].ops);
        const double min_ops = *std::min_element(results[ti].ops.begin(), results[ti].ops.end());
        const double scale = base > 0.0 ? med_ops / base : 0.0;
        const double med_us = median(results[ti].med_thread_us);
        const double hit = median(results[ti].hit_rate) * 100.0;
        const uint64_t sink0 = results[ti].sink.front();
        const bool sink_stable = std::all_of(results[ti].sink.begin(), results[ti].sink.end(),
                                             [&](uint64_t s) { return s == sink0; });
        std::printf("threads=%2d  ops/s=%12.0f  [min %12.0f]  scale=%.2fx  med_thread=%.3f us/op  hit=%.2f%%  sink_stable=%s\n",
                    n, med_ops, min_ops, scale, med_us, hit, sink_stable ? "yes" : "NO");
    }
    std::printf("read_scaling: 1-thread baseline (median ops/s) = %.0f%s; host threads = %u; nothing pinned (OS scheduler placement)\n",
                base, base_was_extra ? " [extra baseline rep]" : "", std::thread::hardware_concurrency());
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    const Options o = parseArgs(argc, argv);
    return runBenchmark(o);
}
