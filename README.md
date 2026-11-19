# TracEon: prepare once, query genomic records locally

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![Release](https://img.shields.io/badge/release-v2.3.0-blueviolet.svg)](RELEASE_v2.3.0.md)
[![Pre-Release](https://img.shields.io/badge/pre--release-v3.0.0--dev-orange.svg)](RELEASE_v3.0.0-dev.md)

TracEon is a C++20 library for repeated whole-record lookup by FASTA/FASTQ
identifier. **Prepare one immutable cache, open it in independent processes,
and query locally.** Sequence views borrow existing storage; FASTQ access can
return owned sequence and quality strings.

The development candidate supports two persistence choices:

- **Mapped:** uncompressed, portable, CRC32C-validated payloads remain file-backed.
  Independent readers share payload pages through the OS; each builds a private
  index. No daemon or network transfer is required for local lookup.
- **LZ4:** smaller artifacts, restored into a private decoded arena per process.
  This remains the default for existing one-argument save/build calls.

```bash
# Build once; also available in a core-only installation.
./build/traceon_make_cache reads.fastq reads.traceon --mapped
```

Each C++ process calls `cache.restore("reads.traceon")` and then `getView(id)`
or `getFastqRecord(id)`. Restore validates every payload byte and rebuilds the
index: it is not instantaneous. Mapped files must remain immutable while open;
TracEon's atomic writer safely replaces the pathname on POSIX.

The [current mapped-storage validation](docs/mapped-storage-validation.md)
checks the implementation and its measured trade-offs. The
[historical rerun](docs/earliest-claims-rerun.md) and
[mapped-cache study](docs/mapped-cache-study.md) support a local repeated-access
benefit on the tested inputs, including shared mapped payload memory. Indexed
readers remain strong alternatives and can use less memory. Preparation costs,
record sizes, access patterns, and reuse determine whether caching helps.
See the [roadmap](ROADMAP.md) for the next workflow validation and binding work.

An optional TCP service provides owned remote records and bounded batching.
Its [persistent-reader](docs/resident-record-study.md) and
[broader competitor](docs/broader-record-study.md) studies include substantial
losses and memory trade-offs. Those service measurements are separate from
local borrowed-view access. TracEon does not claim a new alignment paradigm or
replacement of mapper lookup engines; the
[native mapper control](docs/shared-minimap-index-pilot.md) did not establish
such an advantage.

The latest stable release is **v2.3.0**. The current pre-release development
candidate is **v3.0.0-dev** ([release notes](RELEASE_v3.0.0-dev.md)). It introduces default cache writes
to portable v5 and tightens parser, persistence, k-mer, and client-lifecycle
contracts. See [CHANGELOG.md](CHANGELOG.md) and
[ADR-007](docs/architecture/ADR-007-traceon-v5-portable-format.md).

## Scope and data model

TracEon is designed for repeated whole-sequence-payload lookup after an explicit
build or restore phase. It is not a streaming transformer, coordinate index, or
out-of-core database.

- FASTA header descriptions are discarded; the first token is the identifier.
- FASTQ plus-line metadata is discarded.
- `Cache::getView()` returns only a borrowed sequence view.
- `Cache::get()` returns an owned sequence string; `getFastqRecord()` returns
  owned sequence and quality payloads.
- `Cache::hasSequence()` distinguishes an empty sequence from a missing key.
- The hash index is private to each reader. Payloads use a decoded arena or
  a read-only mapping; mapped pages are managed by the OS. Restore reads all
  payload bytes for validation, so this is not a selective out-of-core index.

Duplicate identifiers in parsed text use first-record-wins behavior. Restored
binary payloads containing duplicate identifiers are rejected as corrupt rather
than silently changing their declared record set.

## Released v2.3.0 versus 3.0.0-dev

| Area | Released v2.3.0 | Current 3.0.0-dev candidate |
|---|---|---|
| Default cache write | v4 | v5 |
| Cache read | v4 | v4 migration plus v5 |
| NGS key persistence | implementation hash in v4 | full identifiers in v5 |
| Parser integrity | malformed suffix could publish a valid prefix | structural validation before publication |
| Cache EOF | trailing physical bytes accepted | physical EOF required |
| Format metadata | inferred after restore | exact `FileFormat` byte in v5 |
| Remote client timeout | segment-level behavior | one absolute operation deadline |

The malformed-prefix publication defect persisted through development commit
`9b6349eddfc9fb009b1088ed78302140889885a9`; it was not fixed in v2.0.0.
The 3.0.0-dev candidate adds the repair and regression tests. This distinction
is important when interpreting historical release notes.

## Build

Fresh builds default to `WITH_NATIVE_INSTRUCTIONS=OFF`, avoiding host-native
zlib-ng compiler flags. For an existing build directory that previously enabled
them, configure with `-DWITH_NATIVE_INSTRUCTIONS=OFF` explicitly. Local tuning
with `ON` is supported, but those binaries may require the builder's CPU features.
System runtime dependencies still constrain binary portability.

Requirements:

- CMake 3.20+
- C++20 compiler
- POSIX threads
- OpenSSL development headers/library for the optional remote service's SHA-256
  identity (`libssl-dev` on Debian/Ubuntu, `openssl@3` with Homebrew)
- Git/network access on first configure for commit-pinned zlib-ng,
  unordered_dense, and Catch2 dependencies

```bash
git clone https://github.com/woosflex/TracEon.git
cd TracEon
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

Remote components default to ON on POSIX. Set `TRACEON_BUILD_SERVER=OFF` for a
core-only build without OpenSSL. To enable them explicitly:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DTRACEON_BUILD_SERVER=ON
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

Install the libraries, headers, and tools into a local prefix:

```bash
cmake --install build --prefix "$PWD/build/install" --component TracEon
```

An external CMake project can consume the installation without the source tree:

```cmake
find_package(TracEon 3 REQUIRED COMPONENTS core)
target_link_libraries(my_worker PRIVATE TracEon::core)
```

Pass `-DCMAKE_PREFIX_PATH=/path/to/install` when configuring that project.
`TracEon::core`, `TracEon::remote`, and the experimental `TracEon::kmer` export their
required include paths, C++20 requirement, and link dependencies. Static codecs
are bundled privately with their licenses; remote builds also need system OpenSSL
Crypto. The installation is relocatable and supports shared-library consumers.
`-DBUILD_TESTING=OFF` omits tests and the Catch2 download; a core-only build installs
the libraries and `traceon_make_cache` without OpenSSL or the remote CLI.
This remains a development package, not a published 3.0.0 release.

## Local API

```cpp
#include "Cache.h"
#include <iostream>

int main() {
    TracEon::Cache cache;
    cache.loadFile("reads.fastq");

    if (cache.hasSequence("read-1")) {
        std::string_view sequence = cache.getView("read-1");
        auto owned = cache.getFastqRecord("read-1");
        if (owned) {
            std::cout << owned->sequence << '\n';
            std::cout << owned->quality << '\n';
        }
    }

    cache.save("reads.traceon", TracEon::CacheStorage::Mapped); // shared payload pages
    cache.restore("reads.traceon");   // reads v4 or v5
}
```

### Borrowed-view lifecycle

`getView()` returns a non-owning `std::string_view`. It is valid only while the
same snapshot remains installed. The application must stop using every borrowed
view before `clearCache()`, reload, or destruction begins. The writer mutex does
not provide reader reclamation. Use `get()` when ownership must cross uncertain
component or task lifetimes.

Manual `set()`/`addEntry()` calls belong to the build phase. They are rejected
after a load until `clearCache()` starts a new build phase.

## Input integrity

The 3.0.0-dev parser accepts the project's documented structural dialect:

- FASTQ records use four lines: `@` header, sequence, `+` separator, quality;
- sequence and quality lengths must match;
- a zero-length read is valid only when both payload lines are empty;
- FASTA records require a line-start `>` header and a sequence payload;
- malformed/truncated input throws before any prefix is published; and
- GZIP truncation, invalid concatenated members, and trailing data are rejected.

An all-empty valid FASTQ is classified as `PROTEIN_FASTQ`, the approved subtype
fallback when no alphabet evidence exists. When any non-empty sequence exists,
format inference samples non-empty sequence data and preserves FASTQ framing.

## Binary cache

The current candidate writes v5 and reads v4/v5. v1/v2/v3 are rejected.

- v5 uses explicit little-endian integer fields and full identifiers in both
  index modes.
- Codec 0 stores raw mapped payloads; codec 1 stores one LZ4 Frame.
- CRC32C detects accidental corruption; it is not authentication.
- The declared stored payload must end at physical EOF.
- Mapped restore validates the full payload and builds an index over borrowed
  file-backed fields. The mapping remains live until clear/reload/destruction.
- LZ4 restore decompresses into owned memory, validates it, builds the index,
  and releases its temporary file mapping. Neither path publishes partial data.
- Earlier development v5 readers reject codec 0; rebuild readers before using
  mapped artifacts. LZ4 v5 and v4 migration reads remain supported.
- Never modify or truncate an open mapped file in place: the OS may invalidate
  readers or deliver changed bytes after initial validation. Use atomic cache
  replacement; existing POSIX readers keep the old inode, new readers open the
  replacement. Clear/reload/destruction still invalidate borrowed views.
- Record iteration/order is unspecified, so byte-identical saves are not
  guaranteed.
- v4 remains a migration reader on compatible legacy platforms.

See [ADR-005](docs/architecture/ADR-005-traceon-v4-binary-format.md) and
[ADR-007](docs/architecture/ADR-007-traceon-v5-portable-format.md).

Create a cache from the command line:

```bash
./build/traceon_make_cache reads.fastq reads.traceon --mapped
# Omit --mapped for the compatible compressed default.
```

## Experimental remote access

Remote access is optional, POSIX-only, unauthenticated, and unencrypted. Use it
only on loopback or an administratively isolated trusted network.

The supported service commands build once, then serve the resident snapshot:

```bash
./build/traceon build reads.fasta reads.traceon
./build/traceon serve reads.traceon --host 127.0.0.1 --port 9876 --max-connections 64
# In any number of independent worker processes, subject to the connection cap:
./build/traceon query --host 127.0.0.1 --port 9876 read-1 read-2 > selected.fasta
./build/traceon query --port 9876 --ids ids.txt > selected.fasta
./build/traceon info --port 9876
# Pin the expected digest from your trusted server launch or dataset manifest:
./build/traceon query --port 9876 --expect-snapshot "$SNAPSHOT_SHA256" --ids ids.txt > selected.fasta
# When serving a FASTQ cache:
./build/traceon query --port 9876 --fastq --ids ids.txt > selected.fastq
# Optional bounded batches; default is one request per ID.
./build/traceon query --port 9876 --batch-size 64 --ids ids.txt > selected.fasta
```

`serve` restores once, computes the content digest, and emits `TRACEON_SERVER READY`
with that digest only after listening. Hashing is startup work, never per query.
Port `0` chooses an available port, reported in that readiness line. SIGINT or
SIGTERM stops admission, closes connections, and joins workers before the cache
is destroyed. Clients need only endpoint and IDs, not the cache/source files.

`query` keeps one connection, emits FASTA by default in request order, and
accepts `--ids -` for one ID per standard-input line. It returns 0 on success,
3 for a missing ID, 2 for invalid arguments, and 1 for an I/O/protocol failure.
Output streams incrementally: if a later ID fails, earlier records can already
be present. Always check exit status before accepting an output file. Duplicate
requests are emitted in order; empty present sequences remain distinct from
missing IDs. `--fastq` emits equal-length sequence and quality from FASTQ snapshots;
it rejects FASTA snapshots. Header descriptions and original plus-line metadata
are not retained; FASTQ output uses a bare `+` line.

Protocol v2 (the default) includes snapshot SHA-256, record count, and supported
operations in the CRC-checked handshake. `--expect-snapshot` rejects a mismatch
before any query/output. The digest identifies logical content, including quality
and detected format, independently of record order or compression; it is not the
cache file's SHA-256 and provides no peer authentication. Startup temporarily
copies/sorts IDs and streams payloads into the digest. See the exact identity and
wire specification in [ADR-008](docs/architecture/ADR-008-snapshot-fastq-protocol.md).
The server also accepts v1 GET clients. Connecting to an old v1 server requires
explicit `query --protocol 1`, which cannot verify snapshots or retrieve FASTQ;
there is no silent protocol downgrade.

The default admission cap is 64 retained connection workers. Excess connections
close before HELLO without creating another worker or request buffer; admitted
clients continue normally. Capacity returns when workers finish and are reaped.
`--timeout-ms` sets the server's per-frame read/write deadline, including idle
clients and slow readers; on `query` it sets the complete request deadline.

```cpp
#include "TraceonClient.h"

TracEon::TraceonClient client("127.0.0.1", 9090);
auto sequence = client.getView("read-1");
if (sequence.has_value()) {
    // present; payload may itself be empty
}
// For a FASTQ snapshot: client.getFastqRecord("read-1") returns sequence + quality.
// To pin identity, pass the expected digest as the constructor's fifth argument.
client.close();
```

A connection serializes request/response pairs. In 3.0.0-dev, one absolute
deadline covers a complete request, and `close()` publishes closing before it
waits: queued/new calls are rejected while at most one admitted request drains.
Object destruction still requires application-level lifetime quiescence.
Remote responses are owned strings, not end-to-end zero-copy views.

The server has no authentication, TLS, authorization, or rate limit and uses a
bounded thread-per-connection model. Completed workers are reaped during service; start/stop calls serialize through
shutdown so a restart cannot reuse descriptors before the old workers exit.

## Experimental k-mer C API

The opt-in `traceon_kmer` target provides minimap2-compatible minimizer-key
semantics. The first lookup freezes the index. Insert/reserve and freeze are
synchronized so no admitted mutation can cross publication. Returned interior
pointers remain valid only while the frozen handle exists. See
`include/kmerindex_c_api.h` and ADR-001.

## Testing and fuzzing

```bash
ctest --test-dir build --output-on-failure
```

The full Linux build registers 21 CTest entries, including service CLI, atomic
save failures, and a relocated external-consumer build when Python is available. Counts vary with build options and
are not assertion counts or coverage measurements. Optional libFuzzer targets cover the
binary-cache loader, GZIP, FASTA/FASTQ, k-mer encoding/API, and TRKI surfaces.
See [fuzz/README.md](fuzz/README.md).

## Performance evidence

The [installed shared-cache study](docs/indexed-mapping-study.md) compares
snapshot-pinned retrieval with samtools indexed access on three public example
workloads. It validates FASTA/FASTQ and alignment parity, and records setup,
affinity, raw timing, and process memory. Observed retrieval gains came with
memory/setup costs and did not establish a general end-to-end mapping speedup.
See the [workflow](docs/research-workflow.md) and
[sustained resource checks](docs/service-resource-study.md).

Older reports under `benchmarks/` and `PERFORMANCE.md` are retained
for provenance but were not produced by the corrected publication-grade harness
and must not be read as current comparative evidence. `read_scaling` and
`remote_bench` are diagnostic benchmark executables, not CTest entries.

A valid comparison must pin inputs and tool versions, fail on subprocess errors,
verify equal hit counts and payload hashes, separate construction/restore/query
phases, report raw repeated samples, and distinguish borrowed from owned output.

## Architecture records

- [ADR-001: lock-free reads and lifecycle](docs/architecture/ADR-001-lock-free-reads.md)
- [ADR-002: GZIP integration](docs/architecture/ADR-002-gzip-integration.md)
- [ADR-003: SIMD parsing and hash map](docs/architecture/ADR-003-simd-parsing-hash-map.md)
- [ADR-004: parallel GZIP and cache](docs/architecture/ADR-004-parallel-gzip-lz4-cache.md)
- [ADR-005: v4 streaming cache](docs/architecture/ADR-005-traceon-v4-binary-format.md)
- [ADR-006: remote access](docs/architecture/ADR-006-traceon-remote-access.md)
- [ADR-007: portable v5](docs/architecture/ADR-007-traceon-v5-portable-format.md)

## Known limitations

- The complete decoded dataset and index must fit memory.
- Borrowed views require reader quiescence across lifecycle changes.
- Cache replacement is atomic on POSIX; see the persistence guarantees below.
- v4 migration payloads retain legacy native-endian limitations.
- Windows does not build the optional POSIX remote component.
- Install/export packaging is not yet implemented.
- Comparative performance remains unvalidated.

### Persistence guarantees

Saves write a unique temporary file beside the destination, synchronize file
data, and replace the destination only when serialization succeeds. POSIX saves
then synchronize the parent directory. Pre-replacement failures preserve the
previous cache and clean up the temporary file. A killed process may leave an
unreferenced `.traceon-tmp-*` file; it never publishes that partial file as the
destination. Existing file permissions are preserved on POSIX; new caches are
created with private permissions. Symlink/device destinations are rejected.

An error explicitly saying `replacement completed` means the new complete cache
is visible but directory synchronization failed, so crash durability is uncertain.
It does not mean the old cache was retained. Successful durability depends on
the filesystem/storage honoring synchronization. Windows synchronizes file data
before replacement, but directory crash durability is not established there.

## Contributing and license

See [CONTRIBUTING.md](CONTRIBUTING.md). TracEon is distributed under the
[MIT License](LICENSE).
