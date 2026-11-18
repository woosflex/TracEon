// Supported cache-once service interface. The libraries own parsing and transport.
#include "Cache.h"
#include "TraceonClient.h"
#include "TraceonServer.h"

#include <charconv>
#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <thread>

namespace {
volatile std::sig_atomic_t stopped = 0;
void stop(int) { stopped = 1; }

unsigned number(const std::string& text, unsigned maximum, bool zero = false) {
    unsigned value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() ||
        value > maximum || (!zero && value == 0))
        throw std::invalid_argument("invalid numeric argument: " + text);
    return value;
}

void usage() {
    std::cout <<
        "Usage:\n"
        "  traceon build INPUT CACHE [--mapped]\n"
        "  traceon serve CACHE [--host 127.0.0.1] [--port 9876]\n"
        "                [--max-connections 64] [--timeout-ms 5000]\n"
        "  traceon query [--host 127.0.0.1] [--port 9876] [--timeout-ms 5000]\n"
        "                [--expect-snapshot SHA256] [--fastq] [--protocol 2] [--batch-size 1]\n"
        "                (--ids FILE|- | ID...)\n"
        "  traceon info [--host 127.0.0.1] [--port 9876] [--timeout-ms 5000]\n"
        "Query writes FASTA by default, or FASTQ with --fastq, in request order.\n"
        "Descriptions/plus-line metadata are not retained. --protocol 1 explicitly\n"
        "enables legacy sequence-only access without snapshot verification.\n"
        "Exit 3 means a missing ID, 2 invalid arguments,\n"
        "1 an I/O or protocol failure. Check exit status before using output.\n";
}

void checkId(const std::string& id) {
    if (id.empty() || id.find_first_of(" \t\r\n") != std::string::npos ||
        id.find('\0') != std::string::npos)
        throw std::invalid_argument("expected a nonempty ID without whitespace or NUL");
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") { usage(); return 0; }
        if (argc < 2) { usage(); return 2; }
        const std::string command = argv[1];
        if (command == "build") {
            if (argc != 4 && !(argc == 5 && std::string(argv[4]) == "--mapped"))
                throw std::invalid_argument("build requires INPUT CACHE [--mapped]");
            TracEon::Cache cache;
            cache.loadFile(argv[2]);
            if (!cache.size()) throw std::runtime_error("input contains no records");
            cache.save(argv[3], argc == 5 ? TracEon::CacheStorage::Mapped : TracEon::CacheStorage::Lz4);
            std::cerr << "wrote " << argv[3] << ": " << cache.size() << " records\n";
            return 0;
        }
        if (command != "serve" && command != "query" && command != "info")
            throw std::invalid_argument("unknown command: " + command);
        const bool serve = command == "serve";
        const bool query = command == "query";
        if (serve && argc < 3) throw std::invalid_argument("serve requires CACHE");
        std::string host = "127.0.0.1", ids_file, expected_snapshot;
        unsigned port = 9876, timeout = 5000, connections = 64, protocol = 2, batch_size = 1;
        bool ids_option = false, literals = false, fastq = false;
        std::vector<std::string> ids;
        for (int i = serve ? 3 : 2; i < argc; ++i) {
            const std::string option = argv[i];
            auto value = [&]() -> std::string {
                if (++i == argc) throw std::invalid_argument("missing value for " + option);
                return argv[i];
            };
            if (!literals && option == "--") { literals = true; continue; }
            if (!literals && option == "--host") host = value();
            else if (!literals && option == "--port") port = number(value(), 65535, serve);
            else if (!literals && option == "--timeout-ms")
                timeout = number(value(), std::numeric_limits<int>::max());
            else if (!literals && serve && option == "--max-connections")
                connections = number(value(), std::numeric_limits<unsigned>::max());
            else if (!literals && query && option == "--expect-snapshot") {
                expected_snapshot = value();
                if (!TracEon::proto::validSnapshotId(expected_snapshot))
                    throw std::invalid_argument("expected a 64-digit lowercase snapshot SHA-256");
            } else if (!literals && query && option == "--fastq") fastq = true;
            else if (!literals && query && option == "--batch-size")
                batch_size = number(value(), TracEon::proto::kMaxBatchRecords);
            else if (!literals && query && option == "--protocol") protocol = number(value(), 2);
            else if (!literals && query && option == "--ids") {
                if (ids_option) throw std::invalid_argument("--ids may occur only once");
                ids_file = value();
                ids_option = true;
            } else if (query && (literals || option.rfind("--", 0) != 0)) ids.push_back(option);
            else throw std::invalid_argument("unexpected argument: " + option);
        }
        if (host.empty()) throw std::invalid_argument("host must not be empty");
        if (serve) {
            // Install handlers before publishing readiness, so supervisors can stop immediately.
            std::signal(SIGINT, stop);
            std::signal(SIGTERM, stop);
            TracEon::Cache cache;
            cache.restore(argv[2]);
            TracEon::TraceonServer server(cache, static_cast<uint16_t>(port), host,
                                          static_cast<int>(timeout), connections);
            server.start();
            std::cout << "TRACEON_SERVER READY host=" << host << " port=" << server.port()
                      << " records=" << cache.size() << " snapshot=" << server.snapshotInfo().id
                      << " protocol=2" << std::endl;
            if (!std::cout) throw std::runtime_error("cannot publish server readiness");
            while (!stopped && server.isRunning())
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            server.stop();
            return 0;
        }
        if (command == "info") {
            TracEon::TraceonClient client(host, static_cast<uint16_t>(port), static_cast<int>(timeout));
            const auto& info = *client.snapshotInfo();
            std::cout << "protocol=2 snapshot=" << info.id << " records=" << info.records
                      << " sequence=1 fastq=" << bool(info.capabilities & TracEon::proto::kFastq)
                      << " batch=" << bool(info.capabilities & TracEon::proto::kBatch) << '\n';
            std::cout.flush();
            if (!std::cout) throw std::runtime_error("failed writing snapshot metadata");
            return 0;
        }
        if (protocol == 1 && (fastq || !expected_snapshot.empty() || batch_size > 1))
            throw std::invalid_argument("FASTQ, batching and snapshot verification require protocol 2");
        if ((ids_option && !ids.empty()) || (!ids_option && ids.empty()))
            throw std::invalid_argument("query requires either --ids FILE|- or one or more IDs");
        std::ifstream input;
        if (ids_option && ids_file != "-") {
            input.open(ids_file);
            if (!input) throw std::runtime_error("cannot open ID file: " + ids_file);
        }
        for (const auto& id : ids) checkId(id);
        TracEon::TraceonClient client(host, static_cast<uint16_t>(port), static_cast<int>(timeout),
                                      5000, expected_snapshot, static_cast<uint8_t>(protocol));
        if (fastq && !(client.snapshotInfo()->capabilities & TracEon::proto::kFastq))
            throw std::invalid_argument("--fastq requires a FASTQ snapshot");
        if (batch_size > 1 && !(client.snapshotInfo()->capabilities & TracEon::proto::kBatch))
            throw std::invalid_argument("server does not support batching");
        auto writeRecord = [&](const std::string& id, const std::string& sequence, const std::string& quality) {
            std::cout << (fastq ? '@' : '>') << id << '\n' << sequence << '\n';
            if (fastq) std::cout << "+\n" << quality << '\n';
            if (!std::cout) throw std::runtime_error("failed writing query output");
        };
        auto emit = [&](const std::string& id) {
            checkId(id);
            if (fastq) {
                auto record = client.getFastqRecord(id);
                if (!record) { std::cerr << "missing ID: " << id << '\n'; return false; }
                writeRecord(id, record->sequence, record->quality);
            } else {
                auto sequence = client.getView(id);
                if (!sequence) { std::cerr << "missing ID: " << id << '\n'; return false; }
                writeRecord(id, *sequence, {});
            }
            if (!std::cout) throw std::runtime_error("failed writing query output");
            return true;
        };
        std::vector<std::string> pending;
        auto flushBatch = [&]() {
            const auto records = client.getBatch(pending, fastq);
            for (size_t i = 0; i < records.size(); ++i) {
                if (!records[i]) { std::cerr << "missing ID: " << pending[i] << '\n'; return false; }
                writeRecord(pending[i], records[i]->sequence, records[i]->quality);
            }
            pending.clear();
            return true;
        };
        auto queryId = [&](const std::string& id) {
            if (batch_size == 1) return emit(id);
            checkId(id);
            pending.push_back(id);
            return pending.size() < batch_size || flushBatch();
        };
        if (ids_option) {
            auto& source = ids_file == "-" ? std::cin : input;
            std::string id;
            bool any = false;
            while (std::getline(source, id)) {
                if (!id.empty() && id.back() == '\r') id.pop_back();
                any = true;
                if (!queryId(id)) return 3;
            }
            if (source.bad()) throw std::runtime_error("failed reading IDs");
            if (!any) throw std::invalid_argument("empty ID list");
        } else for (const auto& id : ids) if (!queryId(id)) return 3;
        if (!pending.empty() && !flushBatch()) return 3;
        std::cout.flush();
        if (!std::cout) throw std::runtime_error("failed writing query output");
        return 0;
    } catch (const std::invalid_argument& error) {
        std::cerr << "traceon: " << error.what() << '\n';
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "traceon: " << error.what() << '\n';
        return 1;
    }
}
