#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "TraceonServer.h"
#include "TraceonProto.h"
#include "Cache.h"

#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <fcntl.h>
#include <poll.h>
#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <utility>

namespace TracEon {

namespace {

// Wake the accept loop's poll() by writing a byte to the wake pipe.
void wakeFd(int wake_write_fd) {
    if (wake_write_fd >= 0) {
        char b = 1;
        ssize_t r = ::write(wake_write_fd, &b, 1);
        (void)r; // best-effort
    }
}

// Best-effort: mark a socket to wake blocked recv() (shutdown) then close it.
void closeSocket(int* fd_ptr) {
    int fd = *fd_ptr;
    *fd_ptr = -1;
    if (fd >= 0) {
        ::shutdown(fd, SHUT_RDWR); // wakes a blocked recv on another thread
        ::close(fd);
    }
}

} // namespace

TraceonServer::TraceonServer(Cache& cache, uint16_t port, std::string host,
                             int recv_timeout_ms)
    : cache_(cache), host_(std::move(host)), port_(port),
      recv_timeout_ms_(recv_timeout_ms > 0 ? recv_timeout_ms
                                           : proto::kReadTimeoutMs) {}

TraceonServer::~TraceonServer() { stop(); }

void TraceonServer::start() {
    std::lock_guard<std::mutex> lock(mu_);
    if (running_) {
        throw std::logic_error("TraceonServer::start(): already running");
    }

    // ── Wake pipe: lets stop() unblock the accept loop immediately ────────
    int p[2] = {-1, -1};
    if (::pipe(p) != 0) {
        throw std::runtime_error("TraceonServer: pipe() failed: " +
                                 std::string(std::strerror(errno)));
    }
    sockfd_control_ = p[0];
    sockfd_wake_ = p[1];

    // ── Listen socket ────────────────────────────────────────────────────
    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) {
        ::close(sockfd_control_); ::close(sockfd_wake_);
        sockfd_control_ = sockfd_wake_ = -1;
        throw std::runtime_error("TraceonServer: socket() failed: " +
                                 std::string(std::strerror(errno)));
    }
    int one = 1;
    ::setsockopt(fd, SOL_SOCKET, SO_REUSEADDR, &one, sizeof(one));
#if defined(SO_NOSIGPIPE)
    ::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port_);
    if (inet_pton(AF_INET, host_.c_str(), &addr.sin_addr) != 1) {
        ::close(fd); ::close(sockfd_control_); ::close(sockfd_wake_);
        sockfd_control_ = sockfd_wake_ = -1;
        throw std::invalid_argument("TraceonServer: invalid bind host: " + host_);
    }

    if (::bind(fd, reinterpret_cast<const sockaddr*>(&addr), sizeof(addr)) != 0) {
        std::string msg = std::string("TraceonServer: bind(") + host_ + ":" +
                          std::to_string(port_) + ") failed: " +
                          std::strerror(errno);
        ::close(fd); ::close(sockfd_control_); ::close(sockfd_wake_);
        sockfd_control_ = sockfd_wake_ = -1;
        throw std::runtime_error(msg);
    }
    if (::listen(fd, 128) != 0) {
        std::string msg = std::string("TraceonServer: listen() failed: ") +
                          std::strerror(errno);
        ::close(fd); ::close(sockfd_control_); ::close(sockfd_wake_);
        sockfd_control_ = sockfd_wake_ = -1;
        throw std::runtime_error(msg);
    }

    // Record the actual bound port (port 0 → OS-assigned).
    socklen_t alen = sizeof(addr);
    if (::getsockname(fd, reinterpret_cast<sockaddr*>(&addr), &alen) == 0) {
        port_ = ntohs(addr.sin_port);
    }

    listen_fd_ = fd;
    running_ = true;
    stopping_ = false;

    accept_thread_ = std::thread(&TraceonServer::acceptLoop, this);
}

bool TraceonServer::isRunning() const {
    std::lock_guard<std::mutex> lock(mu_);
    return running_;
}

uint16_t TraceonServer::port() const {
    std::lock_guard<std::mutex> lock(mu_);
    return running_ ? port_ : 0;
}

void TraceonServer::acceptLoop() {
    for (;;) {
        struct pollfd pfds[2];
        pfds[0] = {listen_fd_, POLLIN, 0};
        pfds[1] = {sockfd_control_, POLLIN, 0};
        int pr = ::poll(pfds, 2, -1);
        if (pr < 0) {
            if (errno == EINTR) continue;
            break; // fatal
        }
        {
            std::lock_guard<std::mutex> lock(mu_);
            if (stopping_) break;
        }
        if (pfds[1].revents & POLLIN) break; // stop() woke us

        if (!(pfds[0].revents & POLLIN)) continue;

        sockaddr_in caddr{};
        socklen_t clen = sizeof(caddr);
        int cfd = ::accept(listen_fd_, reinterpret_cast<sockaddr*>(&caddr), &clen);
        if (cfd < 0) {
            if (errno == EINTR) continue;
            {
                std::lock_guard<std::mutex> lock(mu_);
                if (stopping_) break;
            }
            continue; // transient accept error (e.g. EMFILE) — keep serving
        }
        int one = 1;
        ::setsockopt(cfd, IPPROTO_TCP, TCP_NODELAY, &one, sizeof(one));
#if defined(SO_NOSIGPIPE)
        ::setsockopt(cfd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif

        // Thread-per-connection. The worker is tracked so stop() can join it.
        // The stopping_ re-check + spawn happen atomically under mu_ so a
        // connection accepted just as stop() begins is either joined by stop()
        // (spawned before stopping_) or dropped cleanly here — never leaked as
        // an unjoined/untracked thread.
        try {
            std::lock_guard<std::mutex> lock(mu_);
            if (stopping_) {
                ::close(cfd); // shutdown in progress — don't spawn a worker
                return;
            }
            workers_.emplace_back(&TraceonServer::handleClient, this, cfd);
        } catch (const std::system_error&) {
            ::close(cfd); // thread creation failed — drop the connection
        }
    }
}

void TraceonServer::registerClient(int fd) {
    std::lock_guard<std::mutex> lock(mu_);
    client_fds_.insert(fd);
}

void TraceonServer::unregisterClient(int fd) {
    std::lock_guard<std::mutex> lock(mu_);
    client_fds_.erase(fd);
}

void TraceonServer::handleClient(int fd) {
    registerClient(fd);
    bool clean = false;
    try {
        // ── HELLO handshake ───────────────────────────────────────────────
        auto first = proto::readFrame(fd, recv_timeout_ms_);
        if (!first) { clean = true; return; }
        if (first->type != proto::MsgType::Hello) {
            proto::writeFrame(fd, proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest,
                                                 "expected HELLO handshake"));
            return;
        }
        if (first->body.size() != proto::kMagic.size() + 1 ||
            std::memcmp(first->body.data(), proto::kMagic.data(), proto::kMagic.size()) != 0) {
            proto::writeFrame(fd, proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest, "bad HELLO magic"));
            return;
        }
        const uint8_t ver = static_cast<uint8_t>(first->body[proto::kMagic.size()]);
        if (ver != proto::kProtocolVersion) {
            proto::writeFrame(fd, proto::MsgType::Error,
                              proto::encodeError(
                                  proto::kBadRequest,
                                  "unsupported protocol version " +
                                      std::to_string(ver) +
                                      " (server speaks " +
                                      std::to_string(proto::kProtocolVersion) + ")"));
            return;
        }
        if (!proto::writeFrame(fd, proto::MsgType::Ok, proto::encodeOk(proto::kOk, {}))) {
            return; // peer went away during handshake
        }

        // ── Request dispatch loop ─────────────────────────────────────────
        while (auto req = proto::readFrame(fd, recv_timeout_ms_)) {
            bool close_after = false;
            switch (req->type) {
                case proto::MsgType::Get: {
                    const std::string key(req->body.begin(), req->body.end());
                    if (cache_.hasSequence(key)) {
                        const std::string_view seq = cache_.getView(key);
                        std::vector<char> payload(seq.begin(), seq.end());
                        if (!proto::writeFrame(fd, proto::MsgType::Ok,
                                               proto::encodeOk(proto::kOk, payload))) {
                            return;
                        }
                    } else {
                        if (!proto::writeFrame(fd, proto::MsgType::Error,
                                               proto::encodeError(proto::kNotFound,
                                                                  "record key not found"))) {
                            return;
                        }
                    }
                    break;
                }
                case proto::MsgType::Has: {
                    const std::string key(req->body.begin(), req->body.end());
                    std::vector<char> payload(1, cache_.hasSequence(key) ? 1 : 0);
                    if (!proto::writeFrame(fd, proto::MsgType::Ok,
                                           proto::encodeOk(proto::kOk, payload))) {
                        return;
                    }
                    break;
                }
                case proto::MsgType::Stats: {
                    uint64_t entries = static_cast<uint64_t>(cache_.size());
                    uint8_t fmt = static_cast<uint8_t>(cache_.getDetectedFormat());
                    uint8_t mode = static_cast<uint8_t>(cache_.getIndexMode());
                    uint64_t arena = static_cast<uint64_t>(cache_.getArenaBytes());
                    if (!proto::writeFrame(fd, proto::MsgType::OkStats,
                                           proto::encodeOkStats(entries, fmt, mode, arena))) {
                        return;
                    }
                    break;
                }
                case proto::MsgType::Bye:
                    proto::writeFrame(fd, proto::MsgType::Ok, proto::encodeOk(proto::kOk, {}));
                    return;
                default:
                    proto::writeFrame(fd, proto::MsgType::Error,
                                      proto::encodeError(proto::kBadRequest,
                                                         "unknown message type"));
                    close_after = true;
                    break;
            }
            if (close_after) return;
        }
        clean = true;
    } catch (const std::exception& e) {
        // Protocol violations (bad length prefix, corrupt frame) — send
        // an ERROR when the connection is still usable, then close.
        try {
            proto::writeFrame(fd, proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest, e.what()));
        } catch (...) {
        }
    }
    unregisterClient(fd);
    ::shutdown(fd, SHUT_RDWR);
    ::close(fd);
    (void)clean;
}

void TraceonServer::stop() {
    std::vector<std::thread> joinable;
    {
        std::lock_guard<std::mutex> lock(mu_);
        if (!running_) return;
        stopping_ = true;

        // Wake the accept loop by writing to the wake pipe. We deliberately
        // do NOT close the listener / pipe fds here: closing an fd that
        // another thread is blocked in poll() on does not reliably wake poll
        // on Linux (undefined behaviour) and can leave the accept loop stuck
        // forever, deadlocking stop()'s join. The fds are closed only after
        // the accept loop has fully exited, below.
        wakeFd(sockfd_wake_);

        // Wake every blocked worker read (best-effort): shutdown() reliably
        // unblocks a recv() on the same socket. A worker that has already
        // exited may have closed its fd — racing a stale fd number here is
        // benign: the server is shutting down and no new sockets are created.
        for (int fd : client_fds_) ::shutdown(fd, SHUT_RDWR);

        // Collect joinable threads, then release the lock BEFORE joining
        // (workers may need mu_ inside unregisterClient()).
        if (accept_thread_.joinable()) joinable.push_back(std::move(accept_thread_));
        for (auto& t : workers_) {
            if (t.joinable()) joinable.push_back(std::move(t));
        }
        workers_.clear();
        client_fds_.clear();
        running_ = false;
        stopping_ = false;
    }
    // Join first: the accept loop must have returned before we tear down the
    // listener and wake pipe, so we never close an fd it is still poll()ing.
    for (auto& t : joinable) t.join();
    closeSocket(&listen_fd_);
    closeSocket(&sockfd_control_);
    closeSocket(&sockfd_wake_);
}

} // namespace TracEon