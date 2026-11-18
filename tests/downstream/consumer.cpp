#include <TracEon.h>
#include <kmerindex_c_api.h>
#include <stdexcept>
#ifdef WITH_REMOTE
#include <TraceonServer.h>
#include <TraceonClient.h>
#endif

void check_consumer() {
    auto* index = kmerindex_create();
    if (!index) throw std::runtime_error("cannot create index");
    const int inserted = kmerindex_insert(index, 2, 42);
    uint64_t stored_key = 0;
    const auto* value = kmerindex_get(index, 2, &stored_key);
    const bool correct = inserted == 1 && value && *value == 42 && stored_key == 2;
    kmerindex_destroy(index);
    if (!correct) throw std::runtime_error("kmer mismatch");
    TracEon::Cache cache;
    cache.set("id", "ACGT");
    cache.save("consumer.traceon");
    cache.clearCache();
    cache.restore("consumer.traceon");
    if (cache.getView("id") != "ACGT") throw std::runtime_error("cache mismatch");
    cache.save("consumer-mapped.traceon", TracEon::CacheStorage::Mapped);
    cache.clearCache();
    cache.restore("consumer-mapped.traceon");
    if (cache.getView("id") != "ACGT") throw std::runtime_error("mapped cache mismatch");
#ifdef WITH_REMOTE
    TracEon::TraceonServer server(cache, 0);
    server.start();
    TracEon::TraceonClient client("127.0.0.1", server.port(), 1000, 1000, server.snapshotInfo().id);
    if (client.getView("id") != "ACGT") throw std::runtime_error("remote mismatch");
#endif
}
