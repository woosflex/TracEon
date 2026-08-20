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
 * can serve an arbitrary number of concurrent client connections.
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
#include <string>
#include <thread>
#include <mutex>
#include <vector>
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
     * @param recv_timeout_ms  Per-read idle timeout for client connections.
     */
    explicit TraceonServer(Cache& cache, uint16_t port = 9876,
                           std::string host = "127.0.0.1",
                           int recv_timeout_ms = proto::kReadTimeoutMs);
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

private:
    void acceptLoop();
    void handleClient(int fd);
    void registerClient(int fd);
    void unregisterClient(int fd);

    Cache& cache_;
    std::string host_;
    uint16_t port_;
    int recv_timeout_ms_;

    int listen_fd_ = -1;
    int sockfd_control_ = -1; // wake pipe read end used to unblock poll() on the accept loop
    int sockfd_wake_ = -1;    // wake pipe write end

    mutable std::mutex mu_;
    std::thread accept_thread_;
    std::vector<std::thread> workers_;
    std::set<int> client_fds_;
    bool running_ = false;
    bool stopping_ = false;
};

} // namespace TracEon

#endif // TRACEON_TRACEON_SERVER_H