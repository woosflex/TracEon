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
#include <algorithm>
#include <csignal>
#include <openssl/evp.h>

namespace TracEon {

namespace {

proto::SnapshotInfo identifySnapshot(Cache& cache) {
    proto::SnapshotInfo info;
    info.capabilities |= proto::kBatch;
    const auto format = cache.getDetectedFormat();
    if (format == FileFormat::DNA_FASTQ || format == FileFormat::RNA_FASTQ ||
        format == FileFormat::PROTEIN_FASTQ) info.capabilities |= proto::kFastq;
    auto keys = cache.getAllKeys();
    std::sort(keys.begin(), keys.end());
    info.records = keys.size();
    std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> digest(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!digest || EVP_DigestInit_ex(digest.get(), EVP_sha256(), nullptr) != 1)
        throw std::runtime_error("cannot initialize snapshot SHA-256");
    auto update = [&](std::string_view bytes) {
        if (!bytes.empty() && EVP_DigestUpdate(digest.get(), bytes.data(), bytes.size()) != 1)
            throw std::runtime_error("cannot compute snapshot SHA-256");
    };
    auto length = [&](uint64_t n) {
        char bytes[8];
        for (int i = 0; i < 8; ++i) bytes[i] = static_cast<char>((n >> (8 * i)) & 255);
        update({bytes, sizeof(bytes)});
    };
    auto field = [&](std::string_view bytes) { length(bytes.size()); update(bytes); };
    update("TracEon snapshot v1\n");
    const char format_byte = static_cast<char>(format);
    update({&format_byte, 1});
    length(info.records);
    // ponytail: startup copies/sorts IDs; add a canonical cache digest on disk if this cost dominates setup.
    for (const auto& key : keys) {
        field(key);
        const auto rec = cache.getRecordView(key);
        field(rec ? rec->sequence : std::string_view{});
        field(rec ? rec->quality : std::string_view{});
    }
    unsigned char bytes[EVP_MAX_MD_SIZE];
    unsigned size = 0;
    if (EVP_DigestFinal_ex(digest.get(), bytes, &size) != 1 || size != 32)
        throw std::runtime_error("cannot finalize snapshot SHA-256");
    constexpr char hex[] = "0123456789abcdef";
    for (unsigned i = 0; i < size; ++i) {
        info.id += hex[bytes[i] >> 4];
        info.id += hex[bytes[i] & 15];
    }
    return info;
}

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
                             int recv_timeout_ms, std::size_t max_connections)
    : cache_(cache), host_(std::move(host)), port_(port),
      recv_timeout_ms_(recv_timeout_ms > 0 ? recv_timeout_ms
                                           : proto::kReadTimeoutMs),
      max_connections_(max_connections) {
    if (max_connections_ == 0)
        throw std::invalid_argument("TraceonServer: max_connections must be positive");
}

TraceonServer::~TraceonServer() { stop(); }

void TraceonServer::start() {
#ifndef _WIN32
    (void)::signal(SIGPIPE, SIG_IGN);
#endif
    std::lock_guard<std::mutex> lifecycle_lock(lifecycle_mu_);
    std::lock_guard<std::mutex> lock(mu_);
    if (running_) {
        throw std::logic_error("TraceonServer::start(): already running");
    }
    snapshot_ = identifySnapshot(cache_);

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

    // A readiness notification can become stale before accept().
    const int flags = ::fcntl(fd, F_GETFL, 0);
    if (flags < 0 || ::fcntl(fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        running_ = false;
        closeSocket(&listen_fd_);
        closeSocket(&sockfd_control_);
        closeSocket(&sockfd_wake_);
        throw std::runtime_error("TraceonServer: cannot make listener nonblocking");
    }
    try {
        accept_thread_ = std::thread(&TraceonServer::acceptLoop, this);
    } catch (...) {
        running_ = false;
        closeSocket(&listen_fd_);
        closeSocket(&sockfd_control_);
        closeSocket(&sockfd_wake_);
        throw;
    }
}

bool TraceonServer::isRunning() const {
    std::lock_guard<std::mutex> lock(mu_);
    return running_;
}

uint16_t TraceonServer::port() const {
    std::lock_guard<std::mutex> lock(mu_);
    return running_ ? port_ : 0;
}

proto::SnapshotInfo TraceonServer::snapshotInfo() const {
    std::lock_guard<std::mutex> lock(mu_);
    if (!running_) throw std::logic_error("TraceonServer: snapshot is not serving");
    return snapshot_;
}

void TraceonServer::acceptLoop() {
    for (;;) {
        std::list<Worker> finished;
        {
            std::lock_guard<std::mutex> lock(mu_);
            for (auto it = workers_.begin(); it != workers_.end();) {
                // done is published after the worker releases every server lock.
                if (it->done.load(std::memory_order_acquire)) {
                    finished.splice(finished.end(), workers_, it++);
                } else {
                    ++it;
                }
            }
            if (stopping_) break;
        }
        finished.clear();
        struct pollfd pfds[2];
        pfds[0] = {listen_fd_, POLLIN, 0};
        pfds[1] = {sockfd_control_, POLLIN, 0};
        int pr = ::poll(pfds, 2, 100);
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
        int cflags = ::fcntl(cfd, F_GETFL, 0);
        if (cflags >= 0) ::fcntl(cfd, F_SETFL, cflags | O_NONBLOCK);
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
            if (workers_.size() >= max_connections_) {
                // Workers may have finished while poll() waited for this client.
                for (auto it = workers_.begin(); it != workers_.end();) {
                    if (it->done.load(std::memory_order_acquire)) {
                        finished.splice(finished.end(), workers_, it++);
                    } else {
                        ++it;
                    }
                }
            }
            if (workers_.size() >= max_connections_) {
                ::close(cfd); // No worker/buffer allocation and no blocking overload response.
                continue;
            }
            // Register before launch so stop() cannot miss an admitted socket.
            workers_.emplace_back();
            auto& worker = workers_.back();
            try {
                client_fds_.insert(cfd);
                worker.thread = std::thread([this, cfd, &worker] {
                    handleClient(cfd);
                    worker.done.store(true, std::memory_order_release);
                });
            } catch (...) {
                client_fds_.erase(cfd);
                workers_.pop_back();
                throw;
            }
        } catch (const std::exception&) {
            ::close(cfd); // Allocation/thread creation failed; keep serving.
        }
    }
}

void TraceonServer::unregisterClient(int fd) {
    std::lock_guard<std::mutex> lock(mu_);
    client_fds_.erase(fd);
}

void TraceonServer::handleClient(int fd) {
    auto send = [&](proto::MsgType type, const std::vector<char>& body) {
        return proto::writeFrame(fd, type, body, recv_timeout_ms_);
    };
    try {
        // All early returns below exit this lambda, never the function, so
        // the unregister/shutdown/close epilogue after the try/catch always
        // runs exactly once per connection (no descriptor or bookkeeping leak).
        [&] {
        // ── HELLO handshake ───────────────────────────────────────────────
        auto first = proto::readFrame(fd, recv_timeout_ms_);
        if (!first) return;
        if (first->type != proto::MsgType::Hello) {
            send(proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest,
                                                 "expected HELLO handshake"));
            return;
        }
        if (first->body.size() != proto::kMagic.size() + 1 ||
            std::memcmp(first->body.data(), proto::kMagic.data(), proto::kMagic.size()) != 0) {
            send(proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest, "bad HELLO magic"));
            return;
        }
        const uint8_t ver = static_cast<uint8_t>(first->body[proto::kMagic.size()]);
        if (ver != 1 && ver != proto::kProtocolVersion) {
            send(proto::MsgType::Error,
                              proto::encodeError(
                                  proto::kBadRequest,
                                  "unsupported protocol version " +
                                      std::to_string(ver) +
                                      " (server speaks " +
                                      std::to_string(proto::kProtocolVersion) + ")"));
            return;
        }
        if (!send(proto::MsgType::Ok, proto::encodeOk(proto::kOk,
                ver == 1 ? std::vector<char>{} : proto::encodeSnapshot(snapshot_)))) {
            return; // peer went away during handshake
        }

        // ── Request dispatch loop ─────────────────────────────────────────
        while (auto req = proto::readFrame(fd, recv_timeout_ms_)) {
            switch (req->type) {
                case proto::MsgType::GetBatch: {
                    if (ver != 2) throw std::runtime_error("batch retrieval requires protocol v2");
                    const auto [fastq, keys] = proto::decodeBatchRequest(req->body);
                    if (fastq && !(snapshot_.capabilities & proto::kFastq))
                        throw std::runtime_error("FASTQ batch requires a FASTQ snapshot");
                    std::vector<char> payload;
                    proto::append_le32(payload, static_cast<uint32_t>(keys.size()));
                    constexpr size_t limit = proto::kMaxMessageSize - 10;
                    auto room = [&](size_t bytes) {
                        if (bytes > limit - payload.size())
                            throw std::length_error("batch response exceeds frame limit; reduce batch size");
                        const size_t required = payload.size() + bytes;
                        if (required > payload.capacity())
                            payload.reserve(std::min(limit, std::max(required, payload.capacity() * 2)));
                    };
                    for (const auto key_view : keys) {
                        const auto rec = cache_.getRecordView(key_view);
                        if (!rec) {
                            room(9);
                            proto::append_u8(payload, proto::kNotFound);
                            proto::append_le32(payload, 0);
                            proto::append_le32(payload, 0);
                            continue;
                        }
                        const auto sequence = rec->sequence;
                        const size_t remaining = limit - payload.size();
                        if (remaining < 9 || sequence.size() > (remaining - 9) / (fastq ? 2 : 1))
                            throw std::length_error("batch response exceeds frame limit; reduce batch size");
                        const auto quality = fastq ? rec->quality : std::string_view{};
                        if (fastq && sequence.size() != quality.size())
                            throw std::runtime_error("FASTQ sequence and quality lengths differ");
                        room(9 + sequence.size() + quality.size());
                        proto::append_u8(payload, proto::kOk);
                        proto::append_le32(payload, static_cast<uint32_t>(sequence.size()));
                        proto::append_le32(payload, static_cast<uint32_t>(quality.size()));
                        payload.insert(payload.end(), sequence.begin(), sequence.end());
                        payload.insert(payload.end(), quality.begin(), quality.end());
                    }
                    if (!send(proto::MsgType::Ok, proto::encodeOk(proto::kOk, payload))) return;
                    break;
                }
                case proto::MsgType::Get:
                case proto::MsgType::GetFastq: {
                    const bool fastq = req->type == proto::MsgType::GetFastq;
                    if (fastq && (ver == 1 || !(snapshot_.capabilities & proto::kFastq)))
                        throw std::runtime_error("FASTQ retrieval requires v2 and a FASTQ snapshot");
                    const std::string_view key(req->body.data(), req->body.size());
                    const auto rec = cache_.getRecordView(key);
                    if (rec) {
                        const std::string_view seq = rec->sequence;
                        // Validate size before copying either sequence or quality.
                        const auto limit = fastq ? (proto::kMaxMessageSize - 18) / 2
                                                 : proto::kMaxMessageSize - 10;
                        if (seq.size() > limit)
                            throw std::length_error("record exceeds remote frame limit");
                        std::vector<char> payload;
                        if (fastq) {
                            const auto quality = rec->quality;
                            if (seq.size() != quality.size())
                                throw std::runtime_error("FASTQ sequence and quality lengths differ");
                            payload.reserve(8 + seq.size() + quality.size());
                            proto::append_le32(payload, static_cast<uint32_t>(seq.size()));
                            proto::append_le32(payload, static_cast<uint32_t>(quality.size()));
                            payload.insert(payload.end(), seq.begin(), seq.end());
                            payload.insert(payload.end(), quality.begin(), quality.end());
                        } else payload.assign(seq.begin(), seq.end());
                        if (!send(proto::MsgType::Ok,
                                               proto::encodeOk(proto::kOk, payload))) {
                            return;
                        }
                    } else {
                        if (!send(proto::MsgType::Error,
                                               proto::encodeError(proto::kNotFound,
                                                                  "record key not found"))) {
                            return;
                        }
                    }
                    break;
                }
                case proto::MsgType::Bye:
                    send(proto::MsgType::Ok, proto::encodeOk(proto::kOk, {}));
                    return;
                default:
                    send(proto::MsgType::Error,
                                      proto::encodeError(proto::kBadRequest,
                                                         "unknown message type"));
                    return;
            }
        }
        }(); // end of scoped request body; epilogue below always runs
    } catch (const std::exception& e) {
        // Protocol violations (bad length prefix, corrupt frame) — send
        // an ERROR when the connection is still usable, then close.
        try {
            send(proto::MsgType::Error,
                              proto::encodeError(proto::kBadRequest, e.what()));
        } catch (...) {
        }
    }
    unregisterClient(fd);
    ::shutdown(fd, SHUT_RDWR);
    ::close(fd);
}

void TraceonServer::stop() {
    std::lock_guard<std::mutex> lifecycle_lock(lifecycle_mu_);
    std::list<Worker> workers;
    std::thread accept_thread;
    {
        std::lock_guard<std::mutex> lock(mu_);
        if (!running_) return;
        stopping_ = true;
        wakeFd(sockfd_wake_);
        for (int fd : client_fds_) ::shutdown(fd, SHUT_RDWR);
        accept_thread = std::move(accept_thread_);
        workers.splice(workers.end(), workers_);
    }
    // Workers need mu_ to unregister; keep only the lifecycle lock while joining.
    if (accept_thread.joinable()) accept_thread.join();
    workers.clear();
    {
        std::lock_guard<std::mutex> lock(mu_);
        closeSocket(&listen_fd_);
        closeSocket(&sockfd_control_);
        closeSocket(&sockfd_wake_);
        running_ = false;
    }
}

} // namespace TracEon
