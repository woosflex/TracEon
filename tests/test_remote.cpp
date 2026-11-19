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
#include <fcntl.h>
#include <csignal>
#endif

using namespace TracEon;

namespace {

#ifndef _WIN32
struct RemoteSigpipeIgnore {
    RemoteSigpipeIgnore() {
        (void)::signal(SIGPIPE, SIG_IGN);
    }
} g_remote_sigpipe_ignore;
#endif

static int configureTestSocket(int fd) {
    if (fd >= 0) {
#if defined(SO_NOSIGPIPE)
        int one = 1;
        (void)::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif
        int flags = ::fcntl(fd, F_GETFL, 0);
        if (flags >= 0) (void)::fcntl(fd, F_SETFL, flags | O_NONBLOCK);
    }
    return fd;
}

// Bind an ephemeral loopback socket, read back the OS-assigned port, close
// it, and return that port. Best-effort "nothing is listening here" probe
// for connect-refused tests (tiny race window; the caller retries).
uint16_t unusedPort() {
    int fd = configureTestSocket(::socket(AF_INET, SOCK_STREAM, 0));
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
#if defined(SO_NOSIGPIPE)
    int one = 1;
    (void)::setsockopt(fd, SOL_SOCKET, SO_NOSIGPIPE, &one, sizeof(one));
#endif
    sockaddr_in a{};
    a.sin_family = AF_INET;
    a.sin_port = htons(port);
    if (::inet_pton(AF_INET, "127.0.0.1", &a.sin_addr) != 1 ||
        ::connect(fd, reinterpret_cast<const sockaddr*>(&a), sizeof(a)) != 0) {
        ::close(fd);
        return -1;
    }
    configureTestSocket(fd);
    return fd;
}

// ── Protocol unit tests (pure codecs — no socket) ───────────────────────────────
TEST_CASE("remote protocol — framing and message codecs", "[remote][protocol]") {
    SECTION("encodeFrame rejects bodies over the 64 MiB cap") {
        std::vector<char> big(proto::kMaxMessageSize, 'A'); // 1 + 64 MiB > cap
        REQUIRE_THROWS_AS(proto::encodeFrame(proto::MsgType::Get, big),
                          std::length_error);
    }

    SECTION("writeAll obeys its total deadline under backpressure") {
        int sv[2];
        REQUIRE(::socketpair(AF_UNIX, SOCK_STREAM, 0, sv) == 0);
        configureTestSocket(sv[0]);
        configureTestSocket(sv[1]);
        int small = 4096;
        REQUIRE(::setsockopt(sv[0], SOL_SOCKET, SO_SNDBUF, &small, sizeof(small)) == 0);
        const std::vector<char> body(2u * 1024u * 1024u, 'X');
        const auto start = std::chrono::steady_clock::now();
        REQUIRE_FALSE(proto::writeAll(sv[0], body.data(), body.size(), 50));
        const auto elapsed = std::chrono::steady_clock::now() - start;
        REQUIRE(elapsed < std::chrono::seconds(2));
        ::close(sv[0]);
        ::close(sv[1]);
    }

    SECTION("readFrame uses one absolute deadline across header and payload") {
        int sv[2];
        REQUIRE(::socketpair(AF_UNIX, SOCK_STREAM, 0, sv) == 0);
        configureTestSocket(sv[0]);
        configureTestSocket(sv[1]);
        const auto wire = proto::encodeFrame(proto::MsgType::Get,
                                              std::vector<char>(32, 'X'));
        std::jthread trickle([&](std::stop_token st) {
            (void)::send(sv[1], wire.data(), proto::kFrameHeaderSize, MSG_NOSIGNAL);
            for (size_t i = proto::kFrameHeaderSize; i < wire.size() && !st.stop_requested(); ++i) {
                std::this_thread::sleep_for(std::chrono::milliseconds(15));
                if (st.stop_requested()) break;
                if (::send(sv[1], wire.data() + i, 1, MSG_NOSIGNAL) != 1) break;
            }
            ::close(sv[1]);
        });
        const auto start = std::chrono::steady_clock::now();
        REQUIRE_THROWS_AS(proto::readFrame(sv[0], 70), std::runtime_error);
        REQUIRE(std::chrono::steady_clock::now() - start < std::chrono::seconds(1));
        trickle.request_stop();
        ::close(sv[0]);
        trickle.join();
    }

    SECTION("readFrame handles untrusted huge frame length without allocating upfront") {
        {
            int sv[2];
            REQUIRE(::socketpair(AF_UNIX, SOCK_STREAM, 0, sv) == 0);
            configureTestSocket(sv[0]);
            configureTestSocket(sv[1]);
            std::vector<char> wire;
            proto::append_le64(wire, proto::kMaxMessageSize);
            REQUIRE(::write(sv[1], wire.data(), wire.size()) == static_cast<ssize_t>(wire.size()));
            ::close(sv[1]);
            try {
                (void)proto::readFrame(sv[0], 500);
                FAIL("expected exception");
            } catch (const std::runtime_error& e) {
                REQUIRE(std::string(e.what()).find("truncated frame payload") != std::string::npos);
            }
            ::close(sv[0]);
        }
        {
            int sv[2];
            REQUIRE(::socketpair(AF_UNIX, SOCK_STREAM, 0, sv) == 0);
            configureTestSocket(sv[0]);
            configureTestSocket(sv[1]);
            std::vector<char> wire;
            proto::append_le64(wire, proto::kMaxMessageSize);
            proto::append_u8(wire, static_cast<uint8_t>(proto::MsgType::Get));
            REQUIRE(::write(sv[1], wire.data(), wire.size()) == static_cast<ssize_t>(wire.size()));
            ::close(sv[1]);
            try {
                (void)proto::readFrame(sv[0], 500);
                FAIL("expected exception");
            } catch (const std::runtime_error& e) {
                REQUIRE(std::string(e.what()).find("truncated frame payload") != std::string::npos);
            }
            ::close(sv[0]);
        }
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

// ── Client failure lifecycle ────────────────────────────────────────────────
TEST_CASE("remote client closes after transport, protocol, or CRC failure",
          "[remote][client][failure]") {
    using FaultResponse = std::optional<std::pair<proto::MsgType, std::vector<char>>>;
    auto exercise = [](FaultResponse response, bool fastq = false, bool handshake_fault = false, bool batch = false) {
        int listener = configureTestSocket(::socket(AF_INET, SOCK_STREAM, 0));
        REQUIRE(listener >= 0);
        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
        addr.sin_port = 0;
        REQUIRE(::bind(listener, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == 0);
        REQUIRE(::listen(listener, 1) == 0);
        socklen_t addr_len = sizeof(addr);
        REQUIRE(::getsockname(listener, reinterpret_cast<sockaddr*>(&addr), &addr_len) == 0);
        const uint16_t port = ntohs(addr.sin_port);

        std::jthread peer([listener, response = std::move(response), fastq, handshake_fault, batch] {
            struct pollfd p{listener, POLLIN, 0};
            int pr = ::poll(&p, 1, 5000);
            if (pr <= 0 || !(p.revents & POLLIN)) {
                ::close(listener);
                return;
            }
            const int fd = configureTestSocket(::accept(listener, nullptr, nullptr));
            if (fd >= 0) {
                try {
                    auto hello = proto::readFrame(fd, 5000);
                    if (hello && hello->type == proto::MsgType::Hello) {
                        if (handshake_fault) {
                            if (response) (void)proto::writeFrame(fd, response->first, response->second);
                        } else {
                            (void)proto::writeFrame(fd, proto::MsgType::Ok, proto::encodeOk(proto::kOk,
                                proto::encodeSnapshot({proto::kSequence | proto::kFastq | proto::kBatch, 1, std::string(64, 'a')})));
                            auto get = proto::readFrame(fd, 5000);
                            if (get && get->type == (batch ? proto::MsgType::GetBatch : fastq ? proto::MsgType::GetFastq : proto::MsgType::Get) && response)
                                (void)proto::writeFrame(fd, response->first, response->second);
                        }
                    }
                } catch (...) {
                }
                ::shutdown(fd, SHUT_RDWR);
                ::close(fd);
            }
            ::close(listener);
        });

        if (handshake_fault) {
            REQUIRE_THROWS_AS(TraceonClient("127.0.0.1", port), std::runtime_error);
            peer.join();
            return;
        }
        TraceonClient client("127.0.0.1", port);
        REQUIRE(client.isOpen());
        if (batch) REQUIRE_THROWS_AS(client.getBatch({"seq"}, fastq), std::runtime_error);
        else if (fastq) REQUIRE_THROWS_AS(client.getFastqRecord("seq"), std::runtime_error);
        else REQUIRE_THROWS_AS(client.getView("seq"), std::runtime_error);
        REQUIRE_FALSE(client.isOpen());
        REQUIRE_THROWS_AS(client.getView("seq"), std::logic_error);
        REQUIRE_NOTHROW(client.close());
        peer.join();
    };

    SECTION("CRC mismatch") {
        auto corrupt = proto::encodeOk(proto::kOk, {'A', 'C'});
        corrupt.back() ^= 0x55;
        exercise(std::make_pair(proto::MsgType::Ok, std::move(corrupt)));
    }
    SECTION("unexpected response type") {
        exercise(std::make_pair(proto::MsgType::Hello, std::vector<char>{}));
    }
    SECTION("transport EOF") {
        exercise(std::nullopt);
    }
    SECTION("FASTQ malformed lengths and CRC") {
        for (const auto& payload : std::vector<std::vector<char>>{
                {}, {1,0,0,0,2,0,0,0,'A','I','I'},
                {1,0,0,0,1,0,0,0,'A'}, {0,0,0,0,0,0,0,0,'X'},
                {-1,-1,-1,-1,-1,-1,-1,-1}})
            exercise(std::make_pair(proto::MsgType::Ok, proto::encodeOk(proto::kOk, payload)), true);
        auto corrupt = proto::encodeOk(proto::kOk, {1,0,0,0,1,0,0,0,'A','I'});
        corrupt.back() ^= 1;
        exercise(std::make_pair(proto::MsgType::Ok, corrupt), true);
    }
    SECTION("batch malformed results and CRC close the connection") {
        for (const auto& payload : std::vector<std::vector<char>>{
                {}, {0,0,0,0}, {1,0,0,0,0,1,0,0,0,2,0,0,0,'A','I','I'},
                {1,0,0,0,1,1,0,0,0,0,0,0,0,'A'}, {1,0,0,0,9,0,0,0,0,0,0,0,0},
                {1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1},
                {1,0,0,0,0,0,0,0,0,0,0,0,0,'X'}})
            exercise(std::make_pair(proto::MsgType::Ok, proto::encodeOk(proto::kOk, payload)), true, false, true);
        auto corrupt = proto::encodeOk(proto::kOk, {1,0,0,0,0,0,0,0,0,0,0,0,0});
        corrupt.back() ^= 1;
        exercise(std::make_pair(proto::MsgType::Ok, corrupt), true, false, true);
    }
    SECTION("v2 handshake rejects legacy empty payload, malformed metadata and CRC") {
        exercise(std::make_pair(proto::MsgType::Ok, proto::encodeOk(proto::kOk, {})), false, true);
        auto metadata = proto::encodeSnapshot({proto::kSequence, 1, std::string(64, 'a')});
        metadata[0] = 99;
        exercise(std::make_pair(proto::MsgType::Ok, proto::encodeOk(proto::kOk, metadata)), false, true);
        metadata[0] = 2;
        metadata[1] = 0;
        exercise(std::make_pair(proto::MsgType::Ok, proto::encodeOk(proto::kOk, metadata)), false, true);
        metadata[1] = 1;
        auto corrupt = proto::encodeOk(proto::kOk, metadata);
        corrupt.back() ^= 1;
        exercise(std::make_pair(proto::MsgType::Ok, corrupt), false, true);
    }
}

TEST_CASE("remote v2 identifies content and returns complete FASTQ", "[remote][snapshot][fastq]") {
    const std::string fq = "remote_v2.fastq", bin = "remote_v2.traceon";
    { std::ofstream out(fq); out << "@a description\nACGTN\n+description\n!@+IJ\n@empty\n\n+\n\n"; }
    Cache cache;
    cache.loadFile(fq);
    cache.save(bin);
    TraceonServer server(cache, 0);
    server.start();
    const auto info = server.snapshotInfo();
    REQUIRE(info.records == 2);
    REQUIRE(info.capabilities == (proto::kSequence | proto::kFastq | proto::kBatch));
    REQUIRE(proto::validSnapshotId(info.id));
    REQUIRE_THROWS_AS(TraceonClient("127.0.0.1", server.port(), 1000, 1000, std::string(64, '0')),
                      std::runtime_error);
    TraceonClient client("127.0.0.1", server.port(), 1000, 1000, info.id);
    REQUIRE(client.snapshotInfo()->id == info.id);
    const auto record = client.getFastqRecord("a");
    REQUIRE(record);
    REQUIRE(record->sequence == "ACGTN");
    REQUIRE(record->quality == "!@+IJ");
    REQUIRE(client.getView("a") == record->sequence);
    const auto empty = client.getFastqRecord("empty");
    REQUIRE(empty);
    REQUIRE(empty->sequence.empty());
    REQUIRE(empty->quality.empty());
    REQUIRE_FALSE(client.getFastqRecord("missing"));
    REQUIRE(client.getFastqRecord("a")->quality == record->quality);
    const std::vector<std::string> keys{"a", "empty", "missing", "a"};
    const auto batch = client.getBatch(keys, true);
    REQUIRE(batch.size() == keys.size());
    REQUIRE(batch[0]->sequence == record->sequence);
    REQUIRE(batch[0]->quality == record->quality);
    REQUIRE(batch[1]);
    REQUIRE(batch[1]->sequence.empty());
    REQUIRE(batch[1]->quality.empty());
    REQUIRE_FALSE(batch[2]);
    REQUIRE(batch[3]->quality == batch[0]->quality);
    const auto sequences = client.getBatch(keys);
    REQUIRE(sequences[0]->sequence == record->sequence);
    REQUIRE(sequences[0]->quality.empty());
    REQUIRE_FALSE(sequences[2]);
    REQUIRE_THROWS_AS(client.getBatch({}), std::invalid_argument);
    REQUIRE_THROWS_AS(client.getBatch(std::vector<std::string>(257, "a")), std::invalid_argument);
    REQUIRE(client.getBatch(std::vector<std::string>(256, "a"), true).size() == 256);
    std::atomic<bool> batch_ok{true};
    {
        std::vector<std::jthread> readers;
        for (int i = 0; i < 2; ++i) readers.emplace_back([&] {
            try {
                for (int j = 0; j < 20; ++j) {
                    const auto values = client.getBatch(keys, true);
                    if (!values[0] || values[0]->quality != "!@+IJ" || values[2]) batch_ok = false;
                }
            } catch (...) { batch_ok = false; }
        });
    }
    REQUIRE(batch_ok.load());
    TraceonClient legacy("127.0.0.1", server.port(), 1000, 1000, {}, 1);
    REQUIRE_FALSE(legacy.snapshotInfo());
    REQUIRE(legacy.getView("a") == record->sequence);
    REQUIRE_THROWS_AS(legacy.getFastqRecord("a"), std::logic_error);
    REQUIRE_THROWS_AS(legacy.getBatch({"a"}), std::logic_error);
    REQUIRE(legacy.isOpen());
    // A raw v1 peer cannot bypass capability checks by sending the v2 opcode.
    const int fd = connectRaw(server.port());
    REQUIRE(fd >= 0);
    REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, {'T','R','O','-','R','P','C',1}));
    auto hello = proto::readFrame(fd);
    REQUIRE(hello);
    REQUIRE(proto::decodeOk(hello->body).second.empty());
    REQUIRE(proto::writeFrame(fd, proto::MsgType::GetFastq, {'a'}));
    auto error = proto::readFrame(fd);
    REQUIRE(error);
    REQUIRE(error->type == proto::MsgType::Error);
    ::close(fd);
    server.stop();

    Cache restored;
    restored.restore(bin);
    TraceonServer restored_server(restored, 0);
    restored_server.start();
    REQUIRE(restored_server.snapshotInfo().id == info.id);
    restored_server.stop();
    // Ordering and discarded header metadata do not change logical identity.
    { std::ofstream out(fq); out << "@empty renamed\n\n+\n\n@a\nACGTN\n+\n!@+IJ\n"; }
    cache.loadFile(fq);
    server.start();
    REQUIRE(server.snapshotInfo().id == info.id);
    server.stop();
    // A quality-only change must change the identity, even with identical IDs/sequences.
    { std::ofstream out(fq); out << "@empty\n\n+\n\n@a\nACGTN\n+\n!!!!!\n"; }
    cache.loadFile(fq);
    server.start();
    REQUIRE(server.snapshotInfo().id != info.id);
    server.stop();
    std::remove(fq.c_str());
    std::remove(bin.c_str());
}

TEST_CASE("batch limits reject malformed requests and oversized responses", "[remote][batch]") {
    Cache cache;
    cache.set("large", std::string(16 * 1024 * 1024, 'A'));
    TraceonServer server(cache, 0);
    server.start();
    for (const auto& body : std::vector<std::vector<char>>{
            {}, {2,1,0,0,0,0,0,0,0}, {0,0,0,0,0}, {0,1,1,0,0},
            {0,1,0,0,0}, {0,1,0,0,0,9,0,0,0,'x'}, {0,1,0,0,0,0,0,0,0,'x'}}) {
        const int fd = connectRaw(server.port());
        REQUIRE(fd >= 0);
        REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, {'T','R','O','-','R','P','C',2}));
        REQUIRE(proto::readFrame(fd));
        REQUIRE(proto::writeFrame(fd, proto::MsgType::GetBatch, body));
        const auto error = proto::readFrame(fd);
        REQUIRE(error);
        REQUIRE(error->type == proto::MsgType::Error);
        REQUIRE(proto::decodeError(error->body).first == proto::kBadRequest);
        REQUIRE_FALSE(proto::readFrame(fd));
        ::close(fd);
    }
    TraceonClient client("127.0.0.1", server.port());
    REQUIRE_THROWS_AS(client.getBatch({"large"}, true), std::logic_error);
    REQUIRE_THROWS_AS(client.getBatch({"large", "large", "large", "large"}), std::runtime_error);
    REQUIRE_FALSE(client.isOpen());
    TraceonClient recovery("127.0.0.1", server.port());
    REQUIRE(recovery.getBatch({"missing"})[0] == std::nullopt);
    REQUIRE(recovery.getView("large")->size() == 16 * 1024 * 1024);
}

// ── Server + client lifecycle, round-trip, concurrency ──────────────────────
TEST_CASE("remote client-server integration", "[remote][integration]") {
    // Populate a real v5 cache through the standard pipeline: write a small
    // FASTA, parse it, serialize it, then fully restore and validate it — the
    // deployment shape the server is designed to
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
            } catch (const std::runtime_error& e) {
                refused = true;
                REQUIRE(std::string(e.what()).find("Connection refused") != std::string::npos);
                REQUIRE(std::string(e.what()).find("Operation now in progress") == std::string::npos);
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

    SECTION("close waits for a proven in-flight request, with bounded I/O") {
        int listener = configureTestSocket(::socket(AF_INET, SOCK_STREAM, 0));
        REQUIRE(listener >= 0);
        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
        addr.sin_port = 0;
        REQUIRE(::bind(listener, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) == 0);
        REQUIRE(::listen(listener, 1) == 0);
        socklen_t addr_len = sizeof(addr);
        REQUIRE(::getsockname(listener, reinterpret_cast<sockaddr*>(&addr), &addr_len) == 0);
        const uint16_t port = ntohs(addr.sin_port);
        std::atomic<bool> request_seen{false};
        std::atomic<bool> peer_ok{true};
        std::jthread peer([&] {
            struct pollfd p{listener, POLLIN, 0};
            int pr = ::poll(&p, 1, 5000);
            if (pr <= 0 || !(p.revents & POLLIN)) {
                peer_ok.store(false);
                ::close(listener);
                return;
            }
            const int fd = configureTestSocket(::accept(listener, nullptr, nullptr));
            if (fd < 0) { peer_ok.store(false); ::close(listener); return; }
            try {
                auto hello = proto::readFrame(fd, 1000);
                if (!hello || hello->type != proto::MsgType::Hello ||
                    !proto::writeFrame(fd, proto::MsgType::Ok,
                                      proto::encodeOk(proto::kOk, proto::encodeSnapshot(
                                          {proto::kSequence, 1, std::string(64, 'a')})), 1000)) {
                    peer_ok.store(false);
                } else {
                    auto get = proto::readFrame(fd, 1000);
                    if (!get || get->type != proto::MsgType::Get) peer_ok.store(false);
                    else {
                        request_seen.store(true, std::memory_order_release);
                        std::this_thread::sleep_for(std::chrono::milliseconds(100));
                        if (!proto::writeFrame(fd, proto::MsgType::Ok,
                                               proto::encodeOk(proto::kOk,
                                                   {'G','A','T','T','A'}), 1000))
                            peer_ok.store(false);
                    }
                }
            } catch (...) { peer_ok.store(false); }
            ::shutdown(fd, SHUT_RDWR);
            ::close(fd);
            ::close(listener);
        });

        TraceonClient client("127.0.0.1", port, 1000, 1000);
        std::optional<std::string> result;
        std::atomic<bool> reader_ok{true};
        std::jthread reader([&] {
            try { result = client.getView("seq1"); }
            catch (...) { reader_ok.store(false); }
        });
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
        while (!request_seen.load(std::memory_order_acquire) &&
               std::chrono::steady_clock::now() < deadline)
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        REQUIRE(request_seen.load(std::memory_order_acquire));
        std::atomic<bool> queued_started{false};
        std::atomic<bool> queued_rejected{false};
        std::jthread queued([&] {
            queued_started.store(true, std::memory_order_release);
            try { (void)client.getView("queued"); }
            catch (const std::logic_error&) { queued_rejected.store(true); }
            catch (...) {}
        });
        const auto q_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
        while (!queued_started.load(std::memory_order_acquire) &&
               std::chrono::steady_clock::now() < q_deadline)
            std::this_thread::yield();
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        client.close(); // publishes closing, then drains only the admitted request
        reader.join();
        queued.join();
        peer.join();
        REQUIRE(reader_ok.load());
        REQUIRE(queued_rejected.load());
        REQUIRE(peer_ok.load());
        REQUIRE(result == "GATTA");
        REQUIRE_FALSE(client.isOpen());
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

TEST_CASE("Remote workers drain during connection churn and overlapping lifecycle calls",
          "[audit][remote]") {
    TracEon::Cache cache;
    cache.set("r", "ACGT");
    TracEon::TraceonServer server(cache, 0);
    server.start();
    for (int i = 0; i < 128; ++i) {
        TracEon::TraceonClient client("127.0.0.1", server.port());
        REQUIRE(client.getView("r") == "ACGT");
    }
    std::jthread first_stop([&] { server.stop(); });
    std::jthread second_stop([&] { server.stop(); });
    first_stop.join();
    second_stop.join();
    REQUIRE_FALSE(server.isRunning());
    server.start();
    TracEon::TraceonClient client("127.0.0.1", server.port());
    REQUIRE(client.getView("r") == "ACGT");
    client.close();
    server.stop();
}

TEST_CASE("Remote admission cap preserves existing clients and recovers capacity", "[remote][admission]") {
    Cache cache;
    cache.set("r", "ACGT");
    REQUIRE_THROWS_AS(TraceonServer(cache, 0, "127.0.0.1", 1000, 0), std::invalid_argument);
    TraceonServer server(cache, 0, "127.0.0.1", 1000, 2);
#ifdef __linux__
    auto threads = [] {
        return std::distance(fs::directory_iterator("/proc/self/task"), fs::directory_iterator{});
    };
    const auto before = threads();
#endif
    server.start();
    TraceonClient first("127.0.0.1", server.port());
    TraceonClient second("127.0.0.1", server.port());
    for (int i = 0; i < 32; ++i) {
        REQUIRE_THROWS_AS(TraceonClient("127.0.0.1", server.port(), 200), std::runtime_error);
        REQUIRE(first.getView("r") == "ACGT");
        REQUIRE(second.getView("r") == "ACGT");
    }
#ifdef __linux__
    REQUIRE(threads() <= before + 3); // One accept thread, two admitted workers.
#endif
    second.close();
    bool recovered = false;
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
    while (!recovered && std::chrono::steady_clock::now() < deadline) {
        try {
            TraceonClient replacement("127.0.0.1", server.port(), 200);
            recovered = replacement.getView("r") == "ACGT";
        } catch (const std::runtime_error&) {
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
    }
    REQUIRE(recovered);
    first.close();
    server.stop();
}

TEST_CASE("An idle unhandshaken client cannot hold admission forever", "[remote][admission]") {
    Cache cache;
    cache.set("r", "ACGT");
    TraceonServer server(cache, 0, "127.0.0.1", 80, 1);
    server.start();
    const int idle = connectRaw(server.port());
    REQUIRE(idle >= 0);
    // Wait for its read deadline, then prove that a new client can query.
    std::this_thread::sleep_for(std::chrono::milliseconds(150));
    TraceonClient replacement("127.0.0.1", server.port(), 500);
    REQUIRE(replacement.getView("r") == "ACGT");
    ::close(idle);
    replacement.close();
    server.stop();
}

TEST_CASE("Slow readers release their admission slot at the configured deadline", "[remote][admission]") {
    Cache cache;
    cache.set("large", std::string(4 * 1024 * 1024, 'A'));
    TraceonServer server(cache, 0, "127.0.0.1", 80, 1);
    server.start();
    const int fd = connectRaw(server.port());
    REQUIRE(fd >= 0);
    int small = 1024;
    REQUIRE(::setsockopt(fd, SOL_SOCKET, SO_RCVBUF, &small, sizeof(small)) == 0);
    std::vector<char> hello(proto::kMagic.size() + 1);
    std::memcpy(hello.data(), proto::kMagic.data(), proto::kMagic.size());
    hello.back() = proto::kProtocolVersion;
    REQUIRE(proto::writeFrame(fd, proto::MsgType::Hello, hello));
    REQUIRE(proto::readFrame(fd, 500).has_value());
    REQUIRE(proto::writeFrame(fd, proto::MsgType::Get, {'l', 'a', 'r', 'g', 'e'}));
    bool recovered = false;
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
    while (!recovered && std::chrono::steady_clock::now() < deadline) {
        try {
            TraceonClient replacement("127.0.0.1", server.port(), 200);
            recovered = true;
        } catch (const std::runtime_error&) {
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
    }
    ::close(fd);
    REQUIRE(recovered);
    server.stop();
}
