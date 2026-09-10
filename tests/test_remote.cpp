// tests/test_remote.cpp — remote read access domain (v2.3.0 "Harpe")
// (wire-protocol codecs, server/client lifecycle, client-server round-trip
//  against a real v4 cache, concurrency).
//
// Conventions match the other per-domain tests:
//   - single TU, no CATCH_CONFIG_MAIN (links Catch2::Catch2WithMain)
//   - run from build/ (relative temp-file convention)
//   - POSIX-only (the remote slice never builds on Windows)

#include "TestHelpers.h"

#include "TraceonProto.h"
#include "TraceonServer.h"
#include "TraceonClient.h"

#include <utility>

#ifndef _WIN32
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#endif

using namespace TracEon;

namespace {

// Bind an ephemeral loopback socket, read back the OS-assigned port, close
// it, and return that port. Best-effort "nothing is listening here" probe
// for connect-refused tests (tiny race window; the caller retries).
uint16_t unusedPort() {
    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) throw std::runtime_error("unusedPort: socket() failed");
    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    addr.sin_port = 0;
    if (::bind(fd, reinterpret_cast<const sockaddr*>(&addr), sizeof(addr)) != 0) {
        ::close(fd);
        throw std::runtime_error("unusedPort: bind() failed");
    }
    socklen_t alen = sizeof(addr);
    if (::getsockname(fd, reinterpret_cast<sockaddr*>(&addr), &alen) != 0) {
        ::close(fd);
        throw std::runtime_error("unusedPort: getsockname() failed");
    }
    const uint16_t p = ntohs(addr.sin_port);
    ::close(fd);
    return p;
}

} // namespace

// Raw TCP connect to 127.0.0.1:port for protocol-violation tests. Returns
// -1 (after closing) on any failure instead of throwing.
static int connectRaw(uint16_t port) {
    int fd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (fd < 0) return -1;
    sockaddr_in a{};
    a.sin_family = AF_INET;
    a.sin_port = htons(port);
    if (::inet_pton(AF_INET, "127.0.0.1", &a.sin_addr) != 1 ||
        ::connect(fd, reinterpret_cast<const sockaddr*>(&a), sizeof(a)) != 0) {
        ::close(fd);
        return -1;
    }
    return fd;
}

// ── Protocol unit tests (pure codecs — no socket) ───────────────────────────────
TEST_CASE("remote protocol — framing and message codecs", "[remote][protocol]") {
    SECTION("encodeFrame rejects bodies over the 64 MiB cap") {
        std::vector<char> big(proto::kMaxMessageSize, 'A'); // 1 + 64 MiB > cap
        REQUIRE_THROWS_AS(proto::encodeFrame(proto::MsgType::Get, big),
                          std::length_error);
    }

    SECTION("encodeOk/decodeOk round-trip with CRC") {
        std::vector<char> payload{'G', 'A', 'T', 'T', 'A', 'C', 'A'};
        const std::vector<char> body = proto::encodeOk(proto::kOk, payload);
        auto [status, view] = proto::decodeOk(body);
        REQUIRE(status == proto::kOk);
        REQUIRE(view == std::string_view(payload.data(), payload.size()));
    }

    SECTION("decodeOk throws on a forged CRC byte") {
        std::vector<char> body = proto::encodeOk(proto::kOk, {'A', 'C'});
        body.back() ^= 0x5A; // corrupt the trailing CRC32C byte
        REQUIRE_THROWS_AS(proto::decodeOk(body), std::runtime_error);
    }

    SECTION("decodeOk rejects malformed bodies") {
        REQUIRE_THROWS_AS(proto::decodeOk(std::vector<char>(4, 0)),
                          std::runtime_error); // shorter than the 9-byte minimum
        std::vector<char> body = proto::encodeOk(proto::kOk, {});
        body.push_back(0); // plen says 0 but the body is one byte longer
        REQUIRE_THROWS_AS(proto::decodeOk(body), std::runtime_error);
    }

    SECTION("encodeError/decodeError round-trip") {
        const std::vector<char> body = proto::encodeError(proto::kNotFound,
                                                          "record key not found");
        auto [status, msg] = proto::decodeError(body);
        REQUIRE(status == proto::kNotFound);
        REQUIRE(msg == "record key not found");
        REQUIRE_THROWS_AS(proto::decodeError(std::vector<char>(3, 0)),
                          std::runtime_error);
    }
}

// ── Server + client lifecycle, round-trip, concurrency ──────────────────────
TEST_CASE("remote client-server integration", "[remote][integration]") {
    // Populate a real v4 cache through the standard pipeline: write a small
    // FASTA, parse it (loadFile), serialize to v4 (save), and memory-map it
    // back (restore) — the exact deployment shape the server is designed to
    // serve.
    const std::string fa = "remote_test_input.fasta";
    const std::string bin = "remote_test_cache.traceon";
    const std::vector<std::pair<std::string, std::string>> records = {
        {"seq1", "GATTACAGATTACAGATTACA"},
        {"seq2", "CGCGCGCGCGCGCGCGCGCG"},
        {"seq3", "ACGTACGTACGTACGTACGT"},
    };
    {
        std::ofstream out(fa);
        for (const auto& [id, seq] : records)
            out << '>' << id << '\n' << seq << '\n';
    }
    {
        Cache builder;
        builder.loadFile(fa);
        builder.save(bin);
    }
    Cache cache;
    cache.restore(bin);
    REQUIRE(cache.size() == records.size());
    REQUIRE(cache.getDetectedFormat() == FileFormat::DNA_FASTA);
    REQUIRE(cache.getIndexMode() == IndexMode::GENOME);

    SECTION("ephemeral port, HELLO handshake, idempotent BYE/stop, restart") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();
        REQUIRE(server.isRunning());
        const uint16_t p1 = server.port();
        REQUIRE(p1 != 0); // OS-assigned ephemeral port

        TraceonClient client("127.0.0.1", p1);
        REQUIRE(client.isOpen()); // HELLO handshake completed

        client.close(); // BYE
        REQUIRE(!client.isOpen());
        client.close(); // idempotent — must not throw
        REQUIRE(!client.isOpen());

        server.stop();
        REQUIRE(!server.isRunning());
        server.stop(); // idempotent
        REQUIRE(!server.isRunning());

        // Restart on the same port (SO_REUSEADDR) with a fresh handshake.
        server.start();
        REQUIRE(server.port() == p1);
        TraceonClient again("127.0.0.1", p1);
        REQUIRE(again.isOpen());
        again.close();
        server.stop();
    }

    SECTION("connect-refused and invalid-host fail fast") {
        bool refused = false;
        for (int i = 0; i < 5 && !refused; ++i) {
            const uint16_t p = unusedPort();
            try {
                TraceonClient c("127.0.0.1", p); // nothing is listening
                c.close(); // (something grabbed the port — retry)
            } catch (const std::runtime_error&) {
                refused = true;
            }
        }
        REQUIRE(refused);
        REQUIRE_THROWS_AS(TraceonClient("256.256.256.256", 9876),
                          std::invalid_argument);
    }

    SECTION("server rejects a wrong HELLO protocol version") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();
        const uint16_t p = server.port();

        int fd = connectRaw(p);
        REQUIRE(fd >= 0);

        std::vector<char> body(proto::kMagic.begin(), proto::kMagic.end());
        body.push_back(99); // unsupported protocol version
        REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, body));
        auto rsp = proto::readFrame(fd, 5000);
        REQUIRE(rsp.has_value());
        REQUIRE(rsp->type == proto::MsgType::Error);
        auto [st, msg] = proto::decodeError(rsp->body);
        REQUIRE(st == proto::kBadRequest);
        REQUIRE(msg.find("unsupported protocol version") != std::string::npos);

        ::close(fd);
        server.stop();
    }

    // Regression tests for the handleClient() epilogue: every early-exit
    // path (handshake rejection, protocol violation, clean EOF, BYE) must
    // run the unregister/shutdown/close epilogue exactly once — the server
    // must stop() cleanly with no leaked/jammed descriptors or threads.
    SECTION("raw-socket epilogue paths: handshake rejects, violations, EOF, BYE") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();
        const uint16_t p = server.port();

        // (a) HELLO with a corrupted magic string.
        {
            int fd = connectRaw(p);
            REQUIRE(fd >= 0);
            std::vector<char> body{'X', 'X', 'X', '-', 'R', 'P', 'C', 1};
            REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, body));
            auto rsp = proto::readFrame(fd, 5000);
            REQUIRE(rsp.has_value());
            REQUIRE(rsp->type == proto::MsgType::Error);
            auto [st, msg] = proto::decodeError(rsp->body);
            REQUIRE(st == proto::kBadRequest);
            REQUIRE(msg == "bad HELLO magic");
            ::close(fd);
        }

        // (b) First frame is not HELLO (GET before handshake).
        {
            int fd = connectRaw(p);
            REQUIRE(fd >= 0);
            REQUIRE(proto::writeFrame(fd, proto::MsgType::Get, {'k'}));
            auto rsp = proto::readFrame(fd, 5000);
            REQUIRE(rsp.has_value());
            REQUIRE(rsp->type == proto::MsgType::Error);
            auto [st, msg] = proto::decodeError(rsp->body);
            REQUIRE(st == proto::kBadRequest);
            REQUIRE(msg == "expected HELLO handshake");
            ::close(fd);
        }

        // (c) Connect then immediate EOF (no bytes): clean server-side exit.
        {
            int fd = connectRaw(p);
            REQUIRE(fd >= 0);
            ::close(fd); // EOF at frame boundary — not a protocol violation
        }

        // (d) payload_len == 0: implausible-length protocol violation. This
        // is the only client-triggerable path through the catch → epilogue
        // (readFrame throws; the server replies ERROR with the message).
        {
            int fd = connectRaw(p);
            REQUIRE(fd >= 0);
            const char zeros[8] = {0};
            REQUIRE(proto::writeAll(fd, zeros, sizeof(zeros)));
            auto rsp = proto::readFrame(fd, 5000);
            REQUIRE(rsp.has_value());
            REQUIRE(rsp->type == proto::MsgType::Error);
            auto [st, msg] = proto::decodeError(rsp->body);
            REQUIRE(st == proto::kBadRequest);
            REQUIRE(msg.find("implausible payload length") != std::string::npos);
            ::close(fd);
        }

        // (e) Unknown message type after a valid handshake.
        // (f) Explicit BYE round-trip: OK, then the server closes.
        for (int variant = 0; variant < 2; ++variant) {
            int fd = connectRaw(p);
            REQUIRE(fd >= 0);
            std::vector<char> hello(proto::kMagic.begin(), proto::kMagic.end());
            hello.push_back(static_cast<char>(proto::kProtocolVersion));
            REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, hello));
            auto hs = proto::readFrame(fd, 5000);
            REQUIRE(hs.has_value());
            REQUIRE(hs->type == proto::MsgType::Ok);

            if (variant == 0) {
                REQUIRE(proto::writeFrame(fd, static_cast<proto::MsgType>(0x7F), {}));
                auto rsp = proto::readFrame(fd, 5000);
                REQUIRE(rsp.has_value());
                REQUIRE(rsp->type == proto::MsgType::Error);
                auto [st, msg] = proto::decodeError(rsp->body);
                REQUIRE(st == proto::kBadRequest);
                REQUIRE(msg == "unknown message type");
                // Connection must be closed after the violation (EOF next).
                REQUIRE_FALSE(proto::readFrame(fd, 5000).has_value());
            } else {
                REQUIRE(proto::writeFrame(fd, proto::MsgType::Bye, {}));
                auto rsp = proto::readFrame(fd, 5000);
                REQUIRE(rsp.has_value());
                REQUIRE(rsp->type == proto::MsgType::Ok);
                auto [st, payload] = proto::decodeOk(rsp->body);
                REQUIRE(st == proto::kOk);
                REQUIRE(payload.empty());
            }
            ::close(fd);
        }

        // The point of the epilogue guarantee: after all of the above, stop()
        // must join every worker cleanly and promptly (no hang, no abort).
        server.stop();
        REQUIRE(!server.isRunning());
    }

    SECTION("GET round-trip against the real cache") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();
        TraceonClient client("127.0.0.1", server.port());

        // Exact sequences for known keys (zero-copy getView on the server).
        REQUIRE(client.getView("seq1") == records[0].second);
        REQUIRE(client.getView("seq2") == records[1].second);
        REQUIRE(client.getView("seq3") == records[2].second);
        // A missing key is a definitive NOT_FOUND — never garbage.
        REQUIRE(!client.getView("missing").has_value());
        REQUIRE(!client.getView("").has_value());

        client.close();
        server.stop();
    }

    SECTION("concurrent readers — one shared client (documented thread-safe)") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();
        TraceonClient client("127.0.0.1", server.port());

        constexpr int kThreads = 4;
        constexpr int kIters = 1000;
        std::atomic<bool> ok{true};
        std::atomic<size_t> misses{0};
        std::vector<std::thread> threads;
        for (int t = 0; t < kThreads; ++t) {
            threads.emplace_back([&, t] {
                for (int i = 0; i < kIters; ++i) {
                    const auto& [id, seq] = records[(t + i) % records.size()];
                    auto v = client.getView(id);
                    if (!v) {
                        misses.fetch_add(1, std::memory_order_relaxed);
                        ok.store(false);
                    } else if (*v != seq) {
                        ok.store(false);
                    }
                }
            });
        }
        for (auto& th : threads) th.join();
        REQUIRE(ok.load());
        REQUIRE(misses.load() == 0);

        client.close();
        server.stop();
    }

    SECTION("concurrent readers — one client per thread (deployment pattern)") {
        TraceonServer server(cache, 0, "127.0.0.1");
        server.start();

        constexpr int kThreads = 4;
        constexpr int kIters = 1000;
        std::atomic<bool> ok{true};
        std::vector<std::thread> threads;
        for (int t = 0; t < kThreads; ++t) {
            threads.emplace_back([&, t] {
                TraceonClient client("127.0.0.1", server.port());
                for (int i = 0; i < kIters; ++i) {
                    const auto& [id, seq] = records[(t + i) % records.size()];
                    auto v = client.getView(id);
                    if (!v || *v != seq) ok.store(false);
                }
            });
        }
        for (auto& th : threads) th.join();
        REQUIRE(ok.load());

        server.stop();
    }

    fs::remove(fa);
    fs::remove(bin);
}
