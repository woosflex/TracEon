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
 * internal mutex). close() publishes a closing state before waiting for the
 * one admitted request; queued/new requests are rejected, and the admitted
 * request is bounded by one absolute operation deadline. Destruction still
 * requires object-lifetime quiescence: no thread may call a member after the
 * destructor begins. For throughput, use one connection per thread.
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
#include "DecodedRecordTypes.h"

namespace TracEon {

class TraceonClient {
public:
    /**
     * Connect to `host:port` and perform the HELLO handshake.
     * @throws std::runtime_error on connect failure (including connection
     *         refused), handshake failure, or protocol/version mismatch.
     */
    // recv_timeout_ms is retained for source compatibility; in 3.0.0-dev it
    // bounds the complete HELLO or GET operation (write plus response read).
    explicit TraceonClient(std::string host, uint16_t port,
                           int recv_timeout_ms = proto::kReadTimeoutMs,
                           int connect_timeout_ms = 5000,
                           std::string expected_snapshot = {},
                           uint8_t protocol_version = proto::kProtocolVersion);
    ~TraceonClient();

    TraceonClient(const TraceonClient&) = delete;
    TraceonClient& operator=(const TraceonClient&) = delete;

    /**
     * Fetch a record's sequence by key. Returns std::nullopt when the key is
     * not present (NOT_FOUND). Throws std::runtime_error on protocol errors
     * or a CRC32C mismatch on the wire (connection is closed afterwards).
     */
    std::optional<std::string> getView(const std::string& key);

    /** v2 FASTQ only. Owned sequence and quality; nullopt means a missing ID.
     * Unsupported capability throws before sending a request. */
    std::optional<DecodedFastqRecord> getFastqRecord(const std::string& key);

    /** Opt-in bounded batch, in request order with one optional result per ID.
     * Requires v2 batch capability; FASTQ additionally requires quality support.
     * Empty/oversized key lists and oversized response frames are errors. */
    std::vector<std::optional<DecodedFastqRecord>> getBatch(
        const std::vector<std::string>& keys, bool fastq = false);

    /** Immutable handshake metadata; absent only in explicit v1 mode. */
    const std::optional<proto::SnapshotInfo>& snapshotInfo() const { return snapshot_; }

    /** Gracefully drain at most one admitted request, send BYE best-effort, and close. Idempotent. */
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
    int operation_timeout_ms_;
    std::optional<proto::SnapshotInfo> snapshot_;

    std::optional<DecodedFastqRecord> getRecord(const std::string& key, bool fastq);
    proto::Frame exchangeLocked(proto::MsgType type, const std::vector<char>& body,
                                proto::IoDeadline deadline);

    // io_mu_ must be held. send_bye is false for failed transports/protocols.
    void closeLocked(bool send_bye) noexcept;
};

} // namespace TracEon

#endif // TRACEON_TRACEON_CLIENT_H
