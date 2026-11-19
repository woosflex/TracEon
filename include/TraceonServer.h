#ifndef TRACEON_TRACEON_SERVER_H
#define TRACEON_TRACEON_SERVER_H

/**
 * @file TraceonServer.h
 * @brief TCP server for TracEon remote read access (v2.3.0 "Harpe").
 *
 * Serves key-based lookups over a minimal length-prefixed TCP protocol
 * (see include/TraceonProto.h and docs/architecture/ADR-006-traceon-remote-access.md).
 *
 * The server takes a reference to an already-loaded TracEon::Cache and NEVER
 * mutates it: every request is dispatched through the cache's lock-free,
 * zero-copy read path (getView()/hasSequence()). The immutable-after-load and
 * lock-free-read invariants (ADR-001) are preserved exactly — the same snapshot
 * can serve concurrent clients up to its configured admission limit.
 *
 * Lifecycle:
 *   TraceonServer server(cache, port);   // port 0 = OS-assigned ephemeral
 *   server.start();                      // bind + listen + accept thread
 *   uint16_t p = server.port();          // actual bound port
 *   ...                                  // clients connect and query
 *   server.stop();                       // close listener + clients, join
 *
 * Concurrency model: thread-per-connection. Reads are lock-free on the
 * immutable cache, so concurrent clients do not contend with each other
 * (the only serialization is the accept loop and the per-connection socket).
 *
 * @note POSIX-only (v2.3.0 scope). Trusted-network-only: no auth is
 * performed; do not expose the port to untrusted networks.
 */

#include <cstdint>
#include <cstddef>
#include <string>
#include <thread>
#include <mutex>
#include <list>
#include <memory>
#include <atomic>
#include <set>
#include "TraceonProto.h"

namespace TracEon {

class Cache;

class TraceonServer {
public:
    /**
     * @param cache   Reference to a loaded cache. The caller keeps ownership
     *                and must keep the cache alive (and loaded) for the
     *                lifetime of the server. The server never mutates it.
     * @param port    TCP port to bind. 0 lets the OS pick an ephemeral port
     *                (read it back via port() after start()).
     * @param host    Bind address ("127.0.0.1" for loopback, "0.0.0.0" to accept
     *                connections from other hosts).
     * @param recv_timeout_ms  Per-frame read/write deadline, including idle clients.
     * @param max_connections Maximum retained connection workers (must be > 0).
     *        Excess connections close before HELLO; admitted clients are unaffected.
     */
    explicit TraceonServer(Cache& cache, uint16_t port = 9876,
                           std::string host = "127.0.0.1",
                           int recv_timeout_ms = proto::kReadTimeoutMs,
                           std::size_t max_connections = 64);
    ~TraceonServer();

    TraceonServer(const TraceonServer&) = delete;
    TraceonServer& operator=(const TraceonServer&) = delete;

    /** Bind + listen and start the accept loop. Throws on double-start. */
    void start();

    /**
     * Graceful shutdown: stop accepting, wake/close all client sockets, join
     * the accept loop and every connection thread. Idempotent; safe to call
     * from any thread; the server can be start()ed again afterwards.
     */
    void stop();

    /** Actual bound port (valid after start(); with port 0 the OS-assigned
     *  ephemeral port is returned). 0 if not running. */
    uint16_t port() const;

    /** True while the accept loop is running. */
    bool isRunning() const;

    /** Content identity computed before start() publishes readiness. */
    proto::SnapshotInfo snapshotInfo() const;

private:
    void acceptLoop();
    void handleClient(int fd);
    void unregisterClient(int fd);

    Cache& cache_;
    std::string host_;
    uint16_t port_;
    int recv_timeout_ms_;
    std::size_t max_connections_;
    proto::SnapshotInfo snapshot_;

    int listen_fd_ = -1;
    int sockfd_control_ = -1; // wake pipe read end used to unblock poll() on the accept loop
    int sockfd_wake_ = -1;    // wake pipe write end

    std::mutex lifecycle_mu_; // Serializes complete start/stop operations, including joins.
    mutable std::mutex mu_;
    std::thread accept_thread_;
    struct Worker {
        std::atomic<bool> done{false};
        std::thread thread; // Destroy/join before the completion flag.
        Worker() = default;
        ~Worker() {
            if (thread.joinable()) thread.join();
        }
        Worker(Worker&& o) noexcept : done(o.done.load()), thread(std::move(o.thread)) {}
        Worker& operator=(Worker&& o) noexcept {
            if (this != &o) {
                if (thread.joinable()) thread.join();
                done.store(o.done.load());
                thread = std::move(o.thread);
            }
            return *this;
        }
    };
    std::list<Worker> workers_; // Stable addresses until the accept loop reaps them.
    std::set<int> client_fds_;
    bool running_ = false;
    bool stopping_ = false;
};

} // namespace TracEon

#endif // TRACEON_TRACEON_SERVER_H
