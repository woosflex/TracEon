#pragma once
/**
 * @file TraceonProto.h
 * @brief Wire protocol for TracEon remote read access (v2.3.0 "Harpe").
 *
 * A minimal length-prefixed TCP protocol with NO external dependencies
 * (POSIX sockets). See docs/architecture/ADR-006-traceon-remote-access.md
 * for the full specification.
 *
 * Framing
 * -------
 * Every message on the wire is:
 *
 *     [ u64 LE payload_len ][ payload ]
 *
 * `payload_len` counts the leading message-type byte plus the body, so it is
 * always in [1, kMaxMessageSize]. The maximum frame size is 64 MiB; a peer
 * that advertises a larger payload is a protocol violation (BAD_REQUEST) and
 * the connection is closed without attempting to allocate it.
 *
 * Payload layout (first byte = message type UUID):
 *
 *   Client -> Server
 *     HELLO 0x01  body: "TRO-RPC" (7 bytes) + protocol version (u8)
 *     GET   0x02  body: record key bytes (raw string)
 *     BYE   0x05  body: (empty)
 *     GET_FASTQ 0x06 (v2 only) body: record key bytes
 *
 *   Server -> Client
 *     OK       0x81  body: status(u8) | plen(u32 LE) | payload[plen] | crc32c(u32 LE)
 *     ERROR    0x82  body: status(u8) | mlen(u32 LE) | message[mlen]
 *
 * v1 HELLO OK payload is empty. v2 HELLO OK payload contains SnapshotInfo
 * (encodeSnapshot below). GET_FASTQ OK payload is sequence_len(u32 LE),
 * quality_len(u32 LE), sequence bytes, quality bytes; lengths must match.
 * See ADR-008 for content identity and explicit v1 compatibility.
 *
 * Integrity
 * ---------
 * Every OK-type response carries a trailing CRC32C (Castagnoli, via
 * include/Crc32c.h — the same checksum family as the v4/v5 cache formats).
 * The CRC covers everything in the message body up to (but excluding) the
 * checksum field itself. The client recomputes and compares; a mismatch is a
 * connection error (the client throws and closes). ERROR frames are not CRC'd
 * (they only ever carry error codes + short human-readable text).
 *
 * Status codes (OK body status byte / ERROR body status byte):
 *   0 OK, 1 NOT_FOUND, 2 BAD_REQUEST
 *
 * @note POSIX-only (v2.3.0 scope): this header uses <sys/socket.h> fd
 * helpers. Windows support is out of scope for the remote slice.
 */

#ifndef _WIN32
// POSIX socket headers MUST be included at global scope BEFORE the
// namespace opens — otherwise their declarations would be emitted inside
// TracEon::proto and the ::qualified calls below would not find them.
#include <sys/socket.h>
#include <poll.h>
#include <unistd.h>
#endif

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <optional>
#include <string_view>
#include <string>
#include <vector>
#include <stdexcept>
#include <chrono>
#include <limits>

#include "Crc32c.h"
#include "DecodedRecordTypes.h"

namespace TracEon::proto {

// ── Constants ────────────────────────────────────────────────────────────────
inline constexpr uint8_t  kProtocolVersion = 2;
inline constexpr size_t   kMaxMessageSize  = 64u * 1024u * 1024u; // 64 MiB cap
inline constexpr size_t   kFrameHeaderSize = sizeof(uint64_t);    // 8-byte prefix
inline constexpr std::string_view kMagic = "TRO-RPC";

// Message types.
enum class MsgType : uint8_t {
    Hello    = 0x01,
    Get      = 0x02,
    Bye      = 0x05,
    GetFastq = 0x06, // v2 only; 0x03/0x04 remain reserved
    GetBatch = 0x07, // v2 + kBatch capability
    Ok       = 0x81,
    Error    = 0x82,
};

// Status / error codes.
enum Status : uint8_t {
    kOk          = 0,
    kNotFound    = 1,
    kBadRequest  = 2,
};

// A decoded frame: the message type plus its raw body bytes (excluding the
// leading msg_type byte, which the decoder already strips).
struct Frame {
    MsgType      type;
    std::vector<char> body;
};

inline constexpr uint32_t kSequence = 1;
inline constexpr uint32_t kFastq = 2;
inline constexpr uint32_t kBatch = 4;
inline constexpr size_t kMaxBatchRecords = 256;
struct SnapshotInfo {
    uint32_t capabilities = kSequence;
    uint64_t records = 0;
    std::string id; // 64 lowercase hex digits: canonical content SHA-256
};

inline bool validSnapshotId(std::string_view id) {
    return id.size() == 64 && id.find_first_not_of("0123456789abcdef") == id.npos;
}

// ── Little-endian serialization helpers ─────────────────────────────────────
inline void append_u8(std::vector<char>& v, uint8_t x) { v.push_back(static_cast<char>(x)); }

inline void append_le32(std::vector<char>& v, uint32_t x) {
    for (int i = 0; i < 4; ++i) v.push_back(static_cast<char>((x >> (8 * i)) & 0xFFu));
}

inline void append_le64(std::vector<char>& v, uint64_t x) {
    for (int i = 0; i < 8; ++i) v.push_back(static_cast<char>((x >> (8 * i)) & 0xFFu));
}

inline uint32_t read_le32(const char* p) {
    uint32_t x = 0;
    for (int i = 0; i < 4; ++i) x |= (static_cast<uint32_t>(static_cast<uint8_t>(p[i])) << (8 * i));
    return x;
}

inline uint64_t read_le64(const char* p) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i) x |= (static_cast<uint64_t>(static_cast<uint8_t>(p[i])) << (8 * i));
    return x;
}

// v2 HELLO OK payload: version(u8), capabilities(u32), count(u64), SHA-256(32 bytes).
inline std::vector<char> encodeSnapshot(const SnapshotInfo& info) {
    if (!validSnapshotId(info.id)) throw std::invalid_argument("invalid snapshot SHA-256");
    std::vector<char> payload;
    append_u8(payload, kProtocolVersion);
    append_le32(payload, info.capabilities);
    append_le64(payload, info.records);
    auto nibble = [](char c) { return c <= '9' ? c - '0' : c - 'a' + 10; };
    for (size_t i = 0; i < info.id.size(); i += 2)
        append_u8(payload, static_cast<uint8_t>((nibble(info.id[i]) << 4) | nibble(info.id[i + 1])));
    return payload;
}

inline SnapshotInfo decodeSnapshot(std::string_view payload) {
    if (payload.size() != 45 || static_cast<uint8_t>(payload[0]) != kProtocolVersion)
        throw std::runtime_error("traceon protocol: malformed v2 snapshot metadata");
    SnapshotInfo info;
    info.capabilities = read_le32(payload.data() + 1);
    info.records = read_le64(payload.data() + 5);
    if (!(info.capabilities & kSequence))
        throw std::runtime_error("traceon protocol: missing sequence capability");
    constexpr char hex[] = "0123456789abcdef";
    for (size_t i = 13; i < payload.size(); ++i) {
        const auto byte = static_cast<uint8_t>(payload[i]);
        info.id += hex[byte >> 4];
        info.id += hex[byte & 15];
    }
    return info;
}

// Batch request: fastq(u8), count(u32), repeated key_len(u32) + key bytes.
inline std::vector<char> encodeBatchRequest(const std::vector<std::string>& keys, bool fastq) {
    if (keys.empty() || keys.size() > kMaxBatchRecords)
        throw std::invalid_argument("batch requires 1 to 256 keys");
    size_t size = 5;
    for (const auto& key : keys) {
        const size_t remaining = kMaxMessageSize - 1 - size;
        if (remaining < 4 || key.size() > remaining - 4)
            throw std::length_error("batch request exceeds frame limit");
        size += 4 + key.size();
    }
    std::vector<char> body;
    body.reserve(size);
    append_u8(body, fastq);
    append_le32(body, static_cast<uint32_t>(keys.size()));
    for (const auto& key : keys) {
        append_le32(body, static_cast<uint32_t>(key.size()));
        body.insert(body.end(), key.begin(), key.end());
    }
    return body;
}

inline std::pair<bool, std::vector<std::string_view>> decodeBatchRequest(const std::vector<char>& body) {
    if (body.size() < 5 || static_cast<uint8_t>(body[0]) > 1)
        throw std::runtime_error("malformed batch request");
    const auto count = read_le32(body.data() + 1);
    if (count == 0 || count > kMaxBatchRecords)
        throw std::runtime_error("batch requires 1 to 256 keys");
    std::vector<std::string_view> keys;
    keys.reserve(count);
    size_t offset = 5;
    for (uint32_t i = 0; i < count; ++i) {
        if (body.size() - offset < 4) throw std::runtime_error("truncated batch key length");
        const size_t size = read_le32(body.data() + offset);
        offset += 4;
        if (size > body.size() - offset) throw std::runtime_error("truncated batch key");
        keys.emplace_back(body.data() + offset, size);
        offset += size;
    }
    if (offset != body.size()) throw std::runtime_error("trailing batch request bytes");
    return {body[0] != 0, std::move(keys)};
}

// Batch OK payload: count(u32), then status(u8), seq_len(u32), qual_len(u32),
// sequence and quality per record. Missing entries carry status=NOT_FOUND and zero lengths.
inline std::vector<std::optional<DecodedFastqRecord>> decodeBatchPayload(
        std::string_view payload, size_t expected_count, bool fastq) {
    if (expected_count == 0 || expected_count > kMaxBatchRecords || payload.size() < 4 ||
        read_le32(payload.data()) != expected_count)
        throw std::runtime_error("batch response count mismatch");
    std::vector<std::optional<DecodedFastqRecord>> records;
    records.reserve(expected_count);
    size_t offset = 4;
    for (size_t i = 0; i < expected_count; ++i) {
        if (payload.size() - offset < 9) throw std::runtime_error("truncated batch record");
        const auto status = static_cast<uint8_t>(payload[offset]);
        const size_t seq_size = read_le32(payload.data() + offset + 1);
        const size_t qual_size = read_le32(payload.data() + offset + 5);
        offset += 9;
        if (status == kNotFound) {
            if (seq_size || qual_size) throw std::runtime_error("missing batch record has payload");
            records.emplace_back(std::nullopt);
            continue;
        }
        if (status != kOk || (fastq ? seq_size != qual_size : qual_size != 0) ||
            seq_size > payload.size() - offset || qual_size > payload.size() - offset - seq_size)
            throw std::runtime_error("invalid batch record status or lengths");
        records.emplace_back(DecodedFastqRecord{std::string(payload.substr(offset, seq_size)),
                                               std::string(payload.substr(offset + seq_size, qual_size))});
        offset += seq_size + qual_size;
    }
    if (offset != payload.size()) throw std::runtime_error("trailing batch response bytes");
    return records;
}

// ── Frame encode / decode (pure, testable without a socket) ─────────────────
/** Encode a frame to the full wire representation (prefix + type + body). */
inline std::vector<char> encodeFrame(MsgType type, const std::vector<char>& body) {
    const size_t payload_len = 1 + body.size();
    if (payload_len > kMaxMessageSize) {
        throw std::length_error("traceon protocol: message body exceeds "
                                "kMaxMessageSize (64 MiB)");
    }
    std::vector<char> wire;
    wire.reserve(kFrameHeaderSize + payload_len);
    append_le64(wire, payload_len);
    append_u8(wire, static_cast<uint8_t>(type));
    wire.insert(wire.end(), body.begin(), body.end());
    return wire;
}

// ── OK / ERROR body builders and validators ─────────────────────────
/**
 * Build an OK body: status(u8) | plen(u32 LE) | payload | crc32c(u32 LE).
 * The CRC covers status + plen + payload (everything before the checksum).
 */
inline std::vector<char> encodeOk(uint8_t status, const std::vector<char>& payload) {
    if (payload.size() > kMaxMessageSize - 10)
        throw std::length_error("OK payload exceeds frame limit");
    std::vector<char> body;
    body.reserve(1 + 4 + payload.size() + 4);
    append_u8(body, status);
    append_le32(body, static_cast<uint32_t>(payload.size()));
    body.insert(body.end(), payload.begin(), payload.end());
    // CRC over everything written so far (excludes the CRC field itself).
    TracEon::Crc32c c;
    c.update(body.data(), body.size());
    append_le32(body, c.finalize());
    return body;
}

/**
 * Decode + verify an OK body. Throws std::runtime_error if the checksum does
 * not match (CRC_MISMATCH condition). Returns the status and payload slice.
 */
inline std::pair<uint8_t, std::string_view> decodeOk(const std::vector<char>& body) {
    if (body.size() < 9) { // status(1) + plen(4) + crc(4) minimum
        throw std::runtime_error("traceon protocol: malformed OK frame (short body)");
    }
    const size_t plen = read_le32(body.data() + 1);
    if (body.size() != 1 + 4 + plen + 4) {
        throw std::runtime_error("traceon protocol: malformed OK frame (length mismatch)");
    }
    const uint8_t status = static_cast<uint8_t>(body[0]);
    const char* payload = body.data() + 1 + 4;
    const uint32_t wire_crc = read_le32(body.data() + 1 + 4 + plen);

    TracEon::Crc32c c;
    c.update(body.data(), 1 + 4 + plen); // status + plen + payload
    const uint32_t computed = c.finalize();
    if (computed != wire_crc) {
        throw std::runtime_error("traceon protocol: CRC32C mismatch in OK frame "
                                 "(corrupt response on the wire)");
    }
    return {status, std::string_view(payload, plen)};
}

/** Build an ERROR body: status(u8) | mlen(u32 LE) | message. */
inline std::vector<char> encodeError(uint8_t status, const std::string& message) {
    std::vector<char> body;
    body.reserve(1 + 4 + message.size());
    append_u8(body, status);
    append_le32(body, static_cast<uint32_t>(message.size()));
    body.insert(body.end(), message.begin(), message.end());
    return body;
}

/** Decode an ERROR body into {status, message}. */
inline std::pair<uint8_t, std::string> decodeError(const std::vector<char>& body) {
    if (body.size() < 5) {
        throw std::runtime_error("traceon protocol: malformed ERROR frame");
    }
    const uint8_t status = static_cast<uint8_t>(body[0]);
    const size_t mlen = read_le32(body.data() + 1);
    if (body.size() != 1 + 4 + mlen) {
        throw std::runtime_error("traceon protocol: malformed ERROR frame (length mismatch)");
    }
    return {status, std::string(body.data() + 1 + 4, mlen)};
}

// ── POSIX socket helpers (Linux/macOS; v2.3.0 scope) ─────────────────────────
#ifndef _WIN32

inline constexpr int kReadTimeoutMs = 30000; // default complete-operation deadline
using IoClock = std::chrono::steady_clock;
using IoDeadline = IoClock::time_point;

inline IoDeadline deadlineAfter(int timeout_ms) {
    return IoClock::now() + std::chrono::milliseconds(timeout_ms);
}

inline int remainingMillis(IoDeadline deadline) {
    const auto remaining = std::chrono::duration_cast<std::chrono::milliseconds>(
        deadline - IoClock::now()).count();
    if (remaining <= 0) return 0;
    return remaining > std::numeric_limits<int>::max()
        ? std::numeric_limits<int>::max() : static_cast<int>(remaining);
}

/** Read exactly n bytes before one caller-owned absolute deadline. */
inline bool readFullUntil(int fd, char* buf, size_t n, IoDeadline deadline) {
    size_t got = 0;
    while (got < n) {
        const int remaining = remainingMillis(deadline);
        if (remaining == 0) return false;
        struct pollfd p{fd, POLLIN, 0};
        const int pr = ::poll(&p, 1, remaining);
        if (pr == 0) return false;
        if (pr < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        const ssize_t r = ::recv(fd, buf + got, n - got, MSG_DONTWAIT);
        if (r == 0) return false;
        if (r < 0) {
            if (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK) continue;
            return false;
        }
        got += static_cast<size_t>(r);
    }
    return true;
}

inline bool readFull(int fd, char* buf, size_t n, int timeout_ms = kReadTimeoutMs) {
    return readFullUntil(fd, buf, n, deadlineAfter(timeout_ms));
}

/** Write all bytes before one caller-owned absolute deadline. */
inline bool writeAllUntil(int fd, const char* buf, size_t n, IoDeadline deadline) {
    size_t sent = 0;
    while (sent < n) {
        const int remaining = remainingMillis(deadline);
        if (remaining == 0) return false;
        struct pollfd p{fd, POLLOUT, 0};
        const int pr = ::poll(&p, 1, remaining);
        if (pr == 0) return false;
        if (pr < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (p.revents & (POLLERR | POLLHUP | POLLNVAL)) return false;
#if defined(MSG_NOSIGNAL)
        const int flags = MSG_NOSIGNAL | MSG_DONTWAIT;
#else
        const int flags = MSG_DONTWAIT;
#endif
        const ssize_t w = ::send(fd, buf + sent, n - sent, flags);
        if (w < 0) {
            if (errno == EINTR || errno == EAGAIN || errno == EWOULDBLOCK) continue;
            return false;
        }
        if (w == 0) return false;
        sent += static_cast<size_t>(w);
    }
    return true;
}

inline bool writeAll(int fd, const char* buf, size_t n,
                     int timeout_ms = kReadTimeoutMs) {
    return writeAllUntil(fd, buf, n, deadlineAfter(timeout_ms));
}

inline std::optional<Frame> readFrameUntil(int fd, IoDeadline deadline) {
    char hdr[kFrameHeaderSize];
    if (!readFullUntil(fd, hdr, sizeof(hdr), deadline)) return std::nullopt;
    const uint64_t payload_len = read_le64(hdr);
    if (payload_len == 0 || payload_len > kMaxMessageSize) {
        throw std::runtime_error("traceon protocol: implausible payload length "
                                 "(0 or > 64 MiB)");
    }
    char type_byte;
    if (!readFullUntil(fd, &type_byte, 1, deadline))
        throw std::runtime_error("traceon protocol: truncated frame payload");
    Frame f;
    f.type = static_cast<MsgType>(static_cast<uint8_t>(type_byte));
    size_t remaining = payload_len - 1;
    // Read body in bounded chunks so untrusted length headers cannot force 64 MiB allocations upfront.
    char chunk[65536];
    while (remaining > 0) {
        const size_t to_read = std::min(sizeof(chunk), remaining);
        if (!readFullUntil(fd, chunk, to_read, deadline))
            throw std::runtime_error("traceon protocol: truncated frame payload");
        if (f.body.size() + to_read > f.body.capacity()) {
            const size_t needed = f.body.size() + to_read;
            const size_t target = std::min(static_cast<size_t>(payload_len - 1),
                                           std::max(needed, f.body.capacity() * 2));
            f.body.reserve(target);
        }
        f.body.insert(f.body.end(), chunk, chunk + to_read);
        remaining -= to_read;
    }
    return f;
}

inline std::optional<Frame> readFrame(int fd, int timeout_ms = kReadTimeoutMs) {
    return readFrameUntil(fd, deadlineAfter(timeout_ms));
}

inline bool writeFrameUntil(int fd, MsgType type, const std::vector<char>& body,
                            IoDeadline deadline) {
    const std::vector<char> wire = encodeFrame(type, body);
    return writeAllUntil(fd, wire.data(), wire.size(), deadline);
}

inline bool writeFrame(int fd, MsgType type, const std::vector<char>& body,
                       int timeout_ms = kReadTimeoutMs) {
    return writeFrameUntil(fd, type, body, deadlineAfter(timeout_ms));
}

#endif // !_WIN32

} // namespace TracEon::proto
