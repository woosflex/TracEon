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
 *     HAS   0x03  body: record key bytes (raw string)
 *     STATS 0x04  body: (empty)
 *     BYE   0x05  body: (empty)
 *
 *   Server -> Client
 *     OK       0x81  body: status(u8) | plen(u32 LE) | payload[plen] | crc32c(u32 LE)
 *     ERROR    0x82  body: status(u8) | mlen(u32 LE) | message[mlen]
 *     OK_STATS 0x83  body: entries(u64 LE) | format(u8) | index_mode(u8)
 *                       | arena_bytes(u64 LE) | crc32c(u32 LE)
 *
 * Integrity
 * ---------
 * Every OK-type response carries a trailing CRC32C (Castagnoli, via
 * include/Crc32c.h — the same checksum family as the v4 binary cache format).
 * The CRC covers everything in the message body up to (but excluding) the
 * checksum field itself. The client recomputes and compares; a mismatch is a
 * connection error (the client throws and closes). ERROR frames are not CRC'd
 * (they only ever carry error codes + short human-readable text).
 *
 * Status codes (OK body status byte / ERROR body status byte):
 *   0 OK, 1 NOT_FOUND, 2 BAD_REQUEST, 3 INTERNAL, 4 CRC_MISMATCH
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

#include "Crc32c.h"

namespace TracEon::proto {

// ── Constants ────────────────────────────────────────────────────────────────
inline constexpr uint8_t  kProtocolVersion = 1;
inline constexpr size_t   kMaxMessageSize  = 64u * 1024u * 1024u; // 64 MiB cap
inline constexpr size_t   kFrameHeaderSize = sizeof(uint64_t);    // 8-byte prefix
inline constexpr std::string_view kMagic = "TRO-RPC";

// Message types.
enum class MsgType : uint8_t {
    Hello    = 0x01,
    Get      = 0x02,
    Has      = 0x03,
    Stats    = 0x04,
    Bye      = 0x05,
    Ok       = 0x81,
    Error    = 0x82,
    OkStats  = 0x83,
};

// Status / error codes.
enum Status : uint8_t {
    kOk          = 0,
    kNotFound    = 1,
    kBadRequest  = 2,
    kInternal    = 3,
    kCrcMismatch = 4,
};

// A decoded frame: the message type plus its raw body bytes (excluding the
// leading msg_type byte, which the decoder already strips).
struct Frame {
    MsgType      type;
    std::vector<char> body;
};

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

/**
 * Decode a frame from a buffer that begins at the 8-byte length prefix.
 * Returns {frame, bytes_consumed}. Throws std::runtime_error when the
 * advertised payload length exceeds kMaxMessageSize.
 */
inline std::pair<Frame, size_t> decodeFrame(const char* buf, size_t avail) {
    if (avail < kFrameHeaderSize) {
        throw std::runtime_error("traceon protocol: truncated frame header");
    }
    const uint64_t payload_len = read_le64(buf);
    if (payload_len == 0 || payload_len > kMaxMessageSize) {
        throw std::runtime_error("traceon protocol: implausible payload length "
                                 "(0 or > 64 MiB)");
    }
    const size_t total = static_cast<size_t>(payload_len) + kFrameHeaderSize;
    if (avail < total) {
        throw std::runtime_error("traceon protocol: truncated frame payload");
    }
    Frame f;
    f.type  = static_cast<MsgType>(static_cast<uint8_t>(buf[kFrameHeaderSize]));
    f.body.assign(buf + kFrameHeaderSize + 1,
                  buf + kFrameHeaderSize + static_cast<size_t>(payload_len));
    return {std::move(f), total};
}

// ── OK / ERROR / OK_STATS body builders and validators ───────────────────────
/**
 * Build an OK body: status(u8) | plen(u32 LE) | payload | crc32c(u32 LE).
 * The CRC covers status + plen + payload (everything before the checksum).
 */
inline std::vector<char> encodeOk(uint8_t status, const std::vector<char>& payload) {
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

/** Build an OK_STATS body with the trailing CRC32C over its fixed fields. */
inline std::vector<char> encodeOkStats(uint64_t entries, uint8_t format,
                                       uint8_t index_mode, uint64_t arena_bytes) {
    std::vector<char> body;
    body.reserve(8 + 1 + 1 + 8 + 4);
    append_le64(body, entries);
    append_u8(body, format);
    append_u8(body, index_mode);
    append_le64(body, arena_bytes);
    TracEon::Crc32c c;
    c.update(body.data(), body.size());
    append_le32(body, c.finalize());
    return body;
}

struct OkStats {
    uint64_t entries = 0;
    uint8_t  format = 0xFF;
    uint8_t  index_mode = 0;
    uint64_t arena_bytes = 0;
};

inline OkStats decodeOkStats(const std::vector<char>& body) {
    if (body.size() != 8 + 1 + 1 + 8 + 4) {
        throw std::runtime_error("traceon protocol: malformed OK_STATS frame");
    }
    const uint32_t wire_crc = read_le32(body.data() + 8 + 1 + 1 + 8);
    TracEon::Crc32c c;
    c.update(body.data(), 8 + 1 + 1 + 8); // entries + format + mode + arena_bytes
    if (c.finalize() != wire_crc) {
        throw std::runtime_error("traceon protocol: CRC32C mismatch in OK_STATS frame");
    }
    OkStats s;
    s.entries      = read_le64(body.data());
    s.format       = static_cast<uint8_t>(body[8]);
    s.index_mode   = static_cast<uint8_t>(body[9]);
    s.arena_bytes  = read_le64(body.data() + 10);
    return s;
}

// ── POSIX socket helpers (Linux/macOS; v2.3.0 scope) ─────────────────────────
#ifndef _WIN32

inline constexpr int kReadTimeoutMs = 30000; // default per-read idle timeout

/**
 * Read exactly n bytes from fd. Blocks until n bytes are available, EOF, or
 * `timeout_ms` elapses with no progress. Returns false on EOF/error/timeout
 * (caller distinguishes via consumed count when partial progress happened).
 */
inline bool readFull(int fd, char* buf, size_t n, int timeout_ms = kReadTimeoutMs) {
    size_t got = 0;
    while (got < n) {
        struct pollfd p{fd, POLLIN, 0};
        int pr = ::poll(&p, 1, timeout_ms);
        if (pr == 0) return false;                      // idle timeout
        if (pr < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        ssize_t r = ::recv(fd, buf + got, n - got, 0);
        if (r == 0) return false;                       // EOF
        if (r < 0) {
            if (errno == EINTR) continue;
            return false;                               // connection error
        }
        got += static_cast<size_t>(r);
    }
    return true;
}

/** Write all bytes; returns false on error. Suppresses SIGPIPE. */
inline bool writeAll(int fd, const char* buf, size_t n) {
    size_t sent = 0;
    while (sent < n) {
#if defined(MSG_NOSIGNAL)
        ssize_t w = ::send(fd, buf + sent, n - sent, MSG_NOSIGNAL);
#else
        ssize_t w = ::send(fd, buf + sent, n - sent, 0);
#endif
        if (w < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        sent += static_cast<size_t>(w);
    }
    return true;
}

/**
 * Read one frame from fd. Returns std::nullopt on clean EOF at a frame
 * boundary (peer closed between frames). Throws std::runtime_error on
 * protocol violations (truncated header/payload, oversize payload) and
 * returns nullopt on idle timeout/connection error (readFull failure).
 */
inline std::optional<Frame> readFrame(int fd, int timeout_ms = kReadTimeoutMs) {
    char hdr[kFrameHeaderSize];
    if (!readFull(fd, hdr, sizeof(hdr), timeout_ms)) return std::nullopt;
    const uint64_t payload_len = read_le64(hdr);
    if (payload_len == 0 || payload_len > kMaxMessageSize) {
        throw std::runtime_error("traceon protocol: implausible payload length "
                                 "(0 or > 64 MiB)");
    }
    std::vector<char> payload(payload_len); // bounded by kMaxMessageSize
    if (!readFull(fd, payload.data(), payload.size(), timeout_ms)) {
        throw std::runtime_error("traceon protocol: truncated frame payload");
    }
    Frame f;
    f.type = static_cast<MsgType>(static_cast<uint8_t>(payload[0]));
    f.body.assign(payload.begin() + 1, payload.end());
    return f;
}

/** Write one frame to fd; returns false on error. */
inline bool writeFrame(int fd, MsgType type, const std::vector<char>& body) {
    const std::vector<char> wire = encodeFrame(type, body);
    return writeAll(fd, wire.data(), wire.size());
}

#endif // !_WIN32

} // namespace TracEon::proto