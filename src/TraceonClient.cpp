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

// Mark the connection dead and throw — used by getView()/stats() when the
// read side of the transport fails mid-pair.
[[noreturn]] void connectionDead() {
    throw std::runtime_error("TraceonClient: connection lost during request");
}

// Used by has(): a transport failure on a membership check must not be
// reported as a definitive "absent" — but we cannot be sure of the answer.
// The honest behavior is to throw: the caller asked for a definitive
// boolean and the transport cannot provide one.
[[noreturn]] bool miss() {
    throw std::runtime_error("TraceonClient: connection lost during HAS request");
}

} // namespace

TraceonClient::TraceonClient(std::string host, uint16_t port,
                             int recv_timeout_ms, int connect_timeout_ms)
    : recv_timeout_ms_(recv_timeout_ms > 0 ? recv_timeout_ms
                                           : proto::kReadTimeoutMs) {
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
        if (pr <= 0 || (p.revents & (POLLERR | POLLHUP | POLLNVAL))) {
            std::string msg = std::string("TraceonClient: connect(") + host + ":" +
                              std::to_string(port) + ") timed out/refused: " +
                              std::strerror(errno);
            ::close(fd);
            throw std::runtime_error(msg);
        }
    }
    setBlocking(fd);
    fd_ = fd;

    // ── HELLO handshake ───────────────────────────────────────────────────
    try {
        std::vector<char> body(proto::kMagic.begin(), proto::kMagic.end());
        body.push_back(static_cast<char>(proto::kProtocolVersion));
        if (!proto::writeFrame(fd_, proto::MsgType::Hello, body)) {
            throw std::runtime_error("TraceonClient: failed to send HELLO");
        }
        std::optional<proto::Frame> rsp = proto::readFrame(fd_, recv_timeout_ms_);
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
        auto [st, _] = proto::decodeOk(rsp->body); // also verifies CRC
        if (st != proto::kOk) {
            throw std::runtime_error("TraceonClient: handshake status " +
                                     std::to_string(st));
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
    std::lock_guard<std::mutex> lock(io_mu_);
    if (!open_.load(std::memory_order_acquire) || fd_ < 0) {
        throw std::logic_error("TraceonClient: connection is closed");
    }
    std::vector<char> body(key.begin(), key.end());
    if (!proto::writeFrame(fd_, proto::MsgType::Get, body)) connectionDead();

    auto rsp = proto::readFrame(fd_, recv_timeout_ms_);
    if (!rsp) connectionDead();
    if (rsp->type == proto::MsgType::Error) {
        auto [st, msg] = proto::decodeError(rsp->body);
        if (st == proto::kNotFound) return std::nullopt;
        throw std::runtime_error("TraceonClient: server error (status " +
                                 std::to_string(st) + "): " + msg);
    }
    if (rsp->type != proto::MsgType::Ok) {
        throw std::runtime_error("TraceonClient: unexpected response type");
    }
    auto [st, payload] = proto::decodeOk(rsp->body); // CRC32C verified here
    if (st != proto::kOk) {
        throw std::runtime_error("TraceonClient: GET status " + std::to_string(st));
    }
    return std::string(payload);
}

bool TraceonClient::has(const std::string& key) {
    std::lock_guard<std::mutex> lock(io_mu_);
    if (!open_.load(std::memory_order_acquire) || fd_ < 0) {
        throw std::logic_error("TraceonClient: connection is closed");
    }
    std::vector<char> body(key.begin(), key.end());
    if (!proto::writeFrame(fd_, proto::MsgType::Has, body)) miss();
    auto rsp = proto::readFrame(fd_, recv_timeout_ms_);
    if (!rsp) miss();
    if (rsp->type == proto::MsgType::Error) {
        auto [st, msg] = proto::decodeError(rsp->body);
        throw std::runtime_error("TraceonClient: server error (status " +
                                 std::to_string(st) + "): " + msg);
    }
    if (rsp->type != proto::MsgType::Ok) {
        throw std::runtime_error("TraceonClient: unexpected response type");
    }
    auto [st, payload] = proto::decodeOk(rsp->body); // CRC32C verified here
    if (st != proto::kOk) {
        throw std::runtime_error("TraceonClient: HAS status " + std::to_string(st));
    }
    return !payload.empty() && payload[0] != 0;
}

TraceonStats TraceonClient::stats() {
    std::lock_guard<std::mutex> lock(io_mu_);
    if (!open_.load(std::memory_order_acquire) || fd_ < 0) {
        throw std::logic_error("TraceonClient: connection is closed");
    }
    if (!proto::writeFrame(fd_, proto::MsgType::Stats, {})) {
        throw std::runtime_error("TraceonClient: failed to send STATS");
    }
    auto rsp = proto::readFrame(fd_, recv_timeout_ms_);
    if (!rsp) throw std::runtime_error("TraceonClient: no STATS response");
    if (rsp->type == proto::MsgType::Error) {
        auto [st, msg] = proto::decodeError(rsp->body);
        throw std::runtime_error("TraceonClient: server error (status " +
                                 std::to_string(st) + "): " + msg);
    }
    if (rsp->type != proto::MsgType::OkStats) {
        throw std::runtime_error("TraceonClient: unexpected response type for STATS");
    }
    const proto::OkStats s = proto::decodeOkStats(rsp->body); // CRC32C verified
    TraceonStats out;
    out.entries     = s.entries;
    out.format      = static_cast<FileFormat>(s.format);
    out.index_mode  = static_cast<IndexMode>(s.index_mode);
    out.arena_bytes = s.arena_bytes;
    return out;
}

void TraceonClient::close() {
    bool expected = true;
    if (!open_.compare_exchange_strong(expected, false, std::memory_order_acq_rel)) {
        return; // already closed
    }
    int fd = fd_;
    fd_ = -1;
    if (fd >= 0) {
        // Best-effort BYE; the server answers OK and closes. If the peer is
        // already gone, MSG_NOSIGNAL keeps this from killing our process.
        std::vector<char> empty;
        proto::writeFrame(fd, proto::MsgType::Bye, empty);
        ::shutdown(fd, SHUT_RDWR);
        ::close(fd);
    }
}

} // namespace TracEon