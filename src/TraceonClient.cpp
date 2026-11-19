#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "TraceonClient.h"
#include "TraceonProto.h"

#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <poll.h>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <csignal>

namespace TracEon {

namespace {

void setNonBlocking(int fd) {
    int flags = ::fcntl(fd, F_GETFL, 0);
    ::fcntl(fd, F_SETFL, flags | O_NONBLOCK);
}

void setBlocking(int fd) {
    int flags = ::fcntl(fd, F_GETFL, 0);
    ::fcntl(fd, F_SETFL, flags & ~O_NONBLOCK);
}

// Mark the connection dead and throw — used by getView() when the
// read side of the transport fails mid-pair.
[[noreturn]] void connectionDead() {
    throw std::runtime_error("TraceonClient: connection lost during request");
}

} // namespace

TraceonClient::TraceonClient(std::string host, uint16_t port,
                             int recv_timeout_ms, int connect_timeout_ms,
                             std::string expected_snapshot, uint8_t protocol_version)
    : operation_timeout_ms_(recv_timeout_ms > 0 ? recv_timeout_ms
                                           : proto::kReadTimeoutMs) {
#ifndef _WIN32
    (void)::signal(SIGPIPE, SIG_IGN);
#endif
    if (protocol_version != 1 && protocol_version != proto::kProtocolVersion)
        throw std::invalid_argument("TraceonClient: supported protocol versions are 1 and 2");
    if (!expected_snapshot.empty() &&
        (protocol_version == 1 || !proto::validSnapshotId(expected_snapshot)))
        throw std::invalid_argument("TraceonClient: snapshot verification requires v2 and a lowercase SHA-256");
    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) {
        throw std::runtime_error("TraceonClient: socket() failed: " +
                                 std::string(std::strerror(errno)));
    }
    int one = 1;
    ::setsockopt(fd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
#if defined(SO_NOSIGPIPE)
    ::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    // Resolve host (accepts both an IPv4 literal like "127.0.0.1" and a
    // hostname like the Docker bridge service name "server") via getaddrinfo.
    // inet_pton alone would reject hostnames, which breaks containerized
    // deployments where the peer is reached by service name.
    struct addrinfo hints{};
    hints.ai_family = AF_INET; // matches the server's AF_INET bind
    hints.ai_socktype = SOCK_STREAM;
    struct addrinfo* res = nullptr;
    if (::getaddrinfo(host.c_str(), nullptr, &hints, &res) != 0 || res == nullptr) {
        ::close(fd);
        throw std::invalid_argument("TraceonClient: invalid host: " + host);
    }
    std::memcpy(&addr, res->ai_addr, sizeof(addr));
    addr.sin_port = htons(port);
    ::freeaddrinfo(res);

    // Non-blocking connect + poll with timeout so a dead host fails fast
    // instead of hanging on the kernel's multi-minute SYN retry.
    setNonBlocking(fd);
    int rc = ::connect(fd, reinterpret_cast<const sockaddr*>(&addr), sizeof(addr));
    if (rc != 0 && errno != EINPROGRESS) {
        std::string msg = std::string("TraceonClient: connect(") + host + ":" +
                          std::to_string(port) + ") failed: " + std::strerror(errno);
        ::close(fd);
        throw std::runtime_error(msg);
    }
    if (rc != 0) {
        struct pollfd p{fd, POLLOUT, 0};
        int pr = ::poll(&p, 1, connect_timeout_ms > 0 ? connect_timeout_ms : 5000);
        int err = 0;
        socklen_t len = sizeof(err);
        if (pr <= 0) {
            err = (pr == 0) ? ETIMEDOUT : errno;
        } else if (::getsockopt(fd, SOL_SOCKET, SO_ERROR, &err, &len) < 0 || err != 0 ||
                   (p.revents & (POLLERR | POLLHUP | POLLNVAL))) {
            if (err == 0) err = errno ? errno : ECONNREFUSED;
        }
        if (err != 0) {
            std::string msg = std::string("TraceonClient: connect(") + host + ":" +
                              std::to_string(port) + ") timed out/refused: " +
                              std::strerror(err);
            ::close(fd);
            throw std::runtime_error(msg);
        }
    }
    setNonBlocking(fd);
    fd_ = fd;

    // ── HELLO handshake ───────────────────────────────────────────────────
    try {
        const auto deadline = proto::deadlineAfter(operation_timeout_ms_);
        std::vector<char> body(proto::kMagic.size() + 1);
        std::memcpy(body.data(), proto::kMagic.data(), proto::kMagic.size());
        body.back() = static_cast<char>(protocol_version);
        if (!proto::writeFrameUntil(fd_, proto::MsgType::Hello, body, deadline)) {
            throw std::runtime_error("TraceonClient: failed to send HELLO");
        }
        std::optional<proto::Frame> rsp = proto::readFrameUntil(fd_, deadline);
        if (!rsp) {
            throw std::runtime_error("TraceonClient: no HELLO response (server closed)");
        }
        if (rsp->type == proto::MsgType::Error) {
            auto [st, msg] = proto::decodeError(rsp->body);
            throw std::runtime_error("TraceonClient: handshake rejected (status " +
                                     std::to_string(st) + "): " + msg);
        }
        if (rsp->type != proto::MsgType::Ok) {
            throw std::runtime_error("TraceonClient: unexpected handshake response");
        }
        auto [st, payload] = proto::decodeOk(rsp->body); // also verifies CRC
        if (st != proto::kOk) {
            throw std::runtime_error("TraceonClient: handshake status " +
                                     std::to_string(st));
        }
        if (protocol_version == 1) {
            if (!payload.empty()) throw std::runtime_error("TraceonClient: unexpected v1 handshake payload");
        } else {
            snapshot_ = proto::decodeSnapshot(payload);
            if (!expected_snapshot.empty() && snapshot_->id != expected_snapshot)
                throw std::runtime_error("TraceonClient: snapshot mismatch: expected " +
                                         expected_snapshot + ", received " + snapshot_->id);
        }
    } catch (...) {
        ::close(fd_);
        fd_ = -1;
        throw;
    }
    open_.store(true, std::memory_order_release);
}

TraceonClient::~TraceonClient() { close(); }

bool TraceonClient::isOpen() const { return open_.load(std::memory_order_acquire); }

std::optional<std::string> TraceonClient::getView(const std::string& key) {
    auto record = getRecord(key, false);
    if (!record) return std::nullopt;
    return std::move(record->sequence);
}

std::optional<DecodedFastqRecord> TraceonClient::getFastqRecord(const std::string& key) {
    if (!snapshot_ || !(snapshot_->capabilities & proto::kFastq))
        throw std::logic_error("TraceonClient: server does not advertise FASTQ retrieval");
    return getRecord(key, true);
}

proto::Frame TraceonClient::exchangeLocked(proto::MsgType type, const std::vector<char>& body,
                                           proto::IoDeadline deadline) {
    if (!proto::writeFrameUntil(fd_, type, body, deadline)) connectionDead();
    auto response = proto::readFrameUntil(fd_, deadline);
    if (!response) connectionDead();
    if (response->type != proto::MsgType::Ok && response->type != proto::MsgType::Error)
        throw std::runtime_error("TraceonClient: unexpected response type");
    return std::move(*response);
}

std::vector<std::optional<DecodedFastqRecord>> TraceonClient::getBatch(
        const std::vector<std::string>& keys, bool fastq) {
    if (!snapshot_ || !(snapshot_->capabilities & proto::kBatch) ||
        (fastq && !(snapshot_->capabilities & proto::kFastq)))
        throw std::logic_error("TraceonClient: server does not advertise requested batch capability");
    if (keys.empty() || keys.size() > proto::kMaxBatchRecords)
        throw std::invalid_argument("batch requires 1 to 256 keys");
    if (!open_.load(std::memory_order_acquire))
        throw std::logic_error("TraceonClient: connection is closed");
    std::lock_guard<std::mutex> lock(io_mu_);
    if (!open_.load(std::memory_order_acquire) || fd_ < 0)
        throw std::logic_error("TraceonClient: connection is closed");
    try {
        const auto deadline = proto::deadlineAfter(operation_timeout_ms_);
        const auto body = proto::encodeBatchRequest(keys, fastq);
        auto response = exchangeLocked(proto::MsgType::GetBatch, body, deadline);
        if (response.type == proto::MsgType::Error) {
            const auto [status, message] = proto::decodeError(response.body);
            throw std::runtime_error("TraceonClient: batch error (status " + std::to_string(status) + "): " + message);
        }
        const auto [status, payload] = proto::decodeOk(response.body);
        if (status != proto::kOk) throw std::runtime_error("TraceonClient: invalid batch status");
        return proto::decodeBatchPayload(payload, keys.size(), fastq);
    } catch (...) {
        closeLocked(false);
        throw;
    }
}

std::optional<DecodedFastqRecord> TraceonClient::getRecord(const std::string& key, bool fastq) {
    // Fast admission check before queueing on io_mu_. close() publishes false
    // before it waits, so calls arriving or queued behind the active request
    // are rejected instead of extending shutdown indefinitely.
    if (!open_.load(std::memory_order_acquire))
        throw std::logic_error("TraceonClient: connection is closed");
    std::lock_guard<std::mutex> lock(io_mu_);
    if (!open_.load(std::memory_order_acquire) || fd_ < 0)
        throw std::logic_error("TraceonClient: connection is closed");
    try {
        const auto deadline = proto::deadlineAfter(operation_timeout_ms_);
        if (key.size() > proto::kMaxMessageSize - 1)
            throw std::length_error("TraceonClient: key exceeds frame limit");
        std::vector<char> body(key.begin(), key.end());
        auto rsp = exchangeLocked(fastq ? proto::MsgType::GetFastq : proto::MsgType::Get, body, deadline);
        if (rsp.type == proto::MsgType::Error) {
            auto [st, msg] = proto::decodeError(rsp.body);
            if (st == proto::kNotFound) return std::nullopt;
            throw std::runtime_error("TraceonClient: server error (status " +
                                     std::to_string(st) + "): " + msg);
        }
        auto [st, payload] = proto::decodeOk(rsp.body); // CRC32C verified here
        if (st != proto::kOk) {
            throw std::runtime_error("TraceonClient: GET status " + std::to_string(st));
        }
        if (!fastq) return DecodedFastqRecord{std::string(payload), {}};
        if (payload.size() < 8)
            throw std::runtime_error("TraceonClient: short FASTQ payload");
        const size_t sequence_size = proto::read_le32(payload.data());
        const size_t quality_size = proto::read_le32(payload.data() + 4);
        if (sequence_size != quality_size || sequence_size > (payload.size() - 8) / 2 ||
            payload.size() - 8 != sequence_size * 2)
            throw std::runtime_error("TraceonClient: FASTQ sequence and quality lengths differ");
        return DecodedFastqRecord{std::string(payload.substr(8, sequence_size)),
                                 std::string(payload.substr(8 + sequence_size))};
    } catch (...) {
        // A request/response pair cannot be resumed after a transport,
        // framing, protocol, or checksum failure. Close while io_mu_ is held
        // so close() cannot race the descriptor or interleave BYE.
        closeLocked(false);
        throw;
    }
}

void TraceonClient::closeLocked(bool send_bye) noexcept {
    open_.store(false, std::memory_order_release);
    const int fd = fd_;
    fd_ = -1;
    if (fd < 0) return;
    if (send_bye) {
        // Best-effort BYE. Encoding allocates, so contain every exception to
        // keep destructor-driven close noexcept in practice.
        try {
            std::vector<char> empty;
            (void)proto::writeFrame(fd, proto::MsgType::Bye, empty, operation_timeout_ms_);
        } catch (...) {
        }
    }
    ::shutdown(fd, SHUT_RDWR);
    ::close(fd);
}

void TraceonClient::close() {
    // Publish closing before waiting for the one admitted request. This is a
    // graceful bounded drain: the active request keeps its absolute deadline;
    // queued/new requests observe false and do not become admitted.
    open_.store(false, std::memory_order_release);
    std::lock_guard<std::mutex> lock(io_mu_);
    if (fd_ < 0) return;
    closeLocked(true);
}

} // namespace TracEon
