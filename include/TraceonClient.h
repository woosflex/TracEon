#ifndef TRACEON_TRACEON_CLIENT_H
#define TRACEON_TRACEON_CLIENT_H

/**
 * @file TraceonClient.h
 * @brief TCP client for TracEon remote read access (v2.3.0 "Harpe").
 *
 * Connects to a TraceonServer and performs key-based lookups over the
 * minimal length-prefixed protocol (include/TraceonProto.h,
 * docs/architecture/ADR-006-traceon-remote-access.md).
 *
 * Every OK response is CRC32C-verified on the client (include/Crc32c.h).
 * A checksum mismatch throws std::runtime_error — corrupt data is never
 * returned to the caller.
 *
 * Thread safety: a single TraceonClient connection is safe for concurrent
 * use from multiple threads (each request/response pair is serialized by an
 * internal mutex). For peak throughput, open one connection per thread and
 * keep each connection single-threaded (the benchmark does this).
 *
 * @note POSIX-only (v2.3.0 scope). The protocol is trusted-network-only:
 * no authentication is performed.
 */

#include <cstdint>
#include <string>
#include <optional>
#include <atomic>
#include <mutex>

#include "TraceonProto.h"
#include "SmartStrategy.h" // FileFormat / IndexMode enums

namespace TracEon {

/** Server statistics as reported by the STATS message (OK_STATS). */
struct TraceonStats {
    uint64_t  entries = 0;
    FileFormat format = FileFormat::UNKNOWN;
    IndexMode  index_mode = IndexMode::GENOME;
    uint64_t  arena_bytes = 0;
};

class TraceonClient {
public:
    /**
     * Connect to `host:port` and perform the HELLO handshake.
     * @throws std::runtime_error on connect failure (including connection
     *         refused), handshake failure, or protocol/version mismatch.
     */
    explicit TraceonClient(std::string host, uint16_t port,
                           int recv_timeout_ms = proto::kReadTimeoutMs,
                           int connect_timeout_ms = 5000);
    ~TraceonClient();

    TraceonClient(const TraceonClient&) = delete;
    TraceonClient& operator=(const TraceonClient&) = delete;

    /**
     * Fetch a record's sequence by key. Returns std::nullopt when the key is
     * not present (NOT_FOUND). Throws std::runtime_error on protocol errors
     * or a CRC32C mismatch on the wire (connection is closed afterwards).
     */
    std::optional<std::string> getView(const std::string& key);

    /** True iff the key exists in the served cache. */
    bool has(const std::string& key);

    /** Fetch server statistics (entries, format, index mode, arena bytes). */
    TraceonStats stats();

    /** Send BYE (best-effort) and close the connection. Idempotent. */
    void close();

    /** True while the connection is open. */
    bool isOpen() const;

private:
    // Serializes each request/response pair so concurrent callers on one
    // connection cannot interleave frames. Held only for the duration of a
    // single round trip — never across two requests.
    mutable std::mutex io_mu_;
    std::atomic<bool> open_{false};
    int fd_ = -1;
    int recv_timeout_ms_;
};

} // namespace TracEon

#endif // TRACEON_TRACEON_CLIENT_H