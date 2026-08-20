# ADR-006: Remote Read Access over TCP (`.traceon` RPC)

**Date:** v2.3.0 (2026) — "Harpe" slice 1
**Status:** Accepted (v2.3.0)
**Deciders:** Adnan Raza (Woosflex)
**Depends on:** ADR-001 (lock-free reads / immutable-after-load),
ADR-005 (`.traceon` v4 binary format + CRC32C)

## Context

TracEon is a high-performance genomic **data cache** — a process-local
library for a single cache snapshot. Consumers that want to serve the same
loaded snapshot to many processes (a pipeline controller, a web backend, a
set of worker jobs) previously had to load the `.traceon` file independently
in each process. For large caches (multi-GB arenas), that multiplies load
time and resident memory, and re-encodes shuffling data around via the
filesystem.

We need a small, dependency-free way to expose a **loaded, immutable snapshot**
to remote clients with the same guarantees the library already gives local
readers: zero-copy views, no per-request allocation on the read path, and
CRC-verified integrity so a corrupted payload is never handed to a caller.

## Decision

Add a minimal **length-prefixed TCP RPC** (POSIX sockets only) built directly
on the existing lock-free read path (`Cache::getView()` / `hasSequence()`).
A `TraceonServer` holds a reference to one already-loaded `Cache` and never
mutates it; every request dispatches through the cache's lock-free,
immutable-after-load read path (ADR-001). Clients are `TraceonClient`
instances that verify every OK response with CRC32C (the same checksum
family as the v4 binary format, ADR-005).

### Framing

Every message on the wire is a single frame:

```
[ u64 LE payload_len ][ payload ]
```

- `payload_len` counts the leading **message-type byte plus the body**, so it
  is always in `[1, kMaxMessageSize]`.
- `kMaxMessageSize = 64 MiB` (a per-record sequence is bounded by the v4
  format's max record size; a peer advertising a larger payload is a protocol
  violation — `BAD_REQUEST` — and the frame is rejected **before any
  allocation**).
- All multi-byte fields are little-endian.

### Message types and layouts

First payload byte = message-type UUID:

| Dir | Type | ID | Body |
|-----|------|----|------|
| C→S | HELLO | `0x01` | `"TRO-RPC"` (7 B) + protocol version (u8) |
| C→S | GET | `0x02` | record key bytes (raw string) |
| C→S | HAS | `0x03` | record key bytes (raw string) |
| C→S | STATS | `0x04` | (empty) |
| C→S | BYE | `0x05` | (empty) |
| S→C | OK | `0x81` | `status(u8) \| plen(u32) \| payload[plen] \| crc32c(u32)` |
| S→C | ERROR | `0x82` | `status(u8) \| mlen(u32) \| message[mlen]` |
| S→C | OK_STATS | `0x83` | `entries(u64) \| format(u8) \| index_mode(u8) \| arena_bytes(u64) \| crc32c(u32)` |

### Integrity model (CRC32C)

- Every **OK-type** response (`OK`, `OK_STATS`) carries a trailing
  CRC32C (Castagnoli, poly `0x1EDC6F41` — the SCTP/iSCSI variant, via
  `include/Crc32c.h`, the same family as the v4 cache format) over all body
  bytes **up to but excluding** the checksum field itself.
- The **client** recomputes and compares on every response. A mismatch throws
  `std::runtime_error` and the connection is closed — **corrupt data is never
  returned** to the caller.
- `ERROR` frames are **not** CRC'd: they only carry a status code plus short
  human-readable text.
- This is accidental-corruption detection, **not authentication** (a deliberate
  attacker can craft collisions).

### Status codes

Shared status byte for OK / ERROR bodies:

| Code | Meaning |
|------|---------|
| `0` | OK |
| `1` | NOT_FOUND |
| `2` | BAD_REQUEST |
| `3` | INTERNAL |
| `4` | CRC_MISMATCH |

### Concurrency model

- **Thread-per-connection**: the accept loop spawns one worker thread per
  client. Reads are lock-free on the immutable cache, so concurrent clients
  do **not** contend with each other (the only serialization is the accept
  loop and the per-connection socket). The server never mutates the cache —
  the same snapshot serves an arbitrary number of connections.
- **One-buffer snapshot**: the server serves exactly the one `Cache` its
  constructor was given. No reload, no per-request lookup cache, no cross-
  request optimistic concurrency to reason about — reads go straight to the
  lock-free `getView()`/`hasSequence()`.
- **Client thread safety**: a single `TraceonClient` serializes each
  request/response pair with an internal mutex, so a connection is safe for
  concurrent use from multiple threads. For peak throughput the deployment
  pattern is **one connection per thread** (the benchmark does this).
- **Shutdown**: `stop()` wakes the accept loop via a self-pipe (write, then
  join, then close the fds — see the lifecycle fix note below), shuts down
  worker sockets to unblock blocked reads, and joins every worker thread.
  `stop()` is idempotent and a server may be `start()`ed again.

### POSIX-only scope

v2.3.0 remote access is **POSIX-only** (`sys/socket.h`, `poll.h`, `unistd.h`).
The `traceon_remote` library is never linked into `traceon_core` (core stays
free of socket code — remote links core, not vice versa) and is gated by
`TRACEON_BUILD_SERVER` (ON by default on non-Windows). Windows support is out
of scope for this slice.

### Security: trusted-network-only

The protocol performs **no authentication or encryption** — it is designed
for a trusted private network (and Docker bridge/testbed). A server binds
`127.0.0.1` by default; bind `0.0.0.0` only when the network is trusted. **Do
not expose the port to the public internet.** If authentication/encryption is
ever required, layer TLS in front of (or in place of) this framing.

## Consequences

### Positive
- Zero dependencies — no gRPC/HTTP/TLS stack to vendor or secure.
- Remote reads preserve the library's invariants: immutable snapshot + lock-
  free dispatch + CRC-verified integrity on every OK response.
- One loaded snapshot serves many clients/processes (no per-process re-load,
  no re-encoding through the filesystem).
- Ephemeral-port (`port 0`) lifecycle makes tests and wrappers trivial.

### Negative
- TCP round-trip latency and per-request copy into a `std::string` on the
  client (remote clients can't get a zero-copy view across a socket).
- No auth/encryption — deploy on trusted networks only.
- POSIX-only for now; Windows clients/servers are out of scope.

## Alternatives Considered

1. **gRPC**: mature but heavy (protobuf toolchain, codegen, HTTP/2, a large
   dependency graph) — overkill for a single immutable look-up endpoint.
   Rejected: violates the zero-dependency, header-light design.
2. **Plain HTTP/REST**: easy to debug, but pulls in an HTTP server/framework
   and parsing overhead, and offers no native CRC integrity enforcement on
   every response. Rejected for the same dependency reason.
3. **Shared-memory / mmap fan-out**: fastest possible fan-out (each process
   mmaps the same `.traceon`), but has no connection/lifecycle model (readers
   must manage snapshot validity themselves), no request/response contract,
   and doesn't serve a cache that only lives in one process's arena today.
   Kept as a future direction, not this slice.

## References
- `include/TraceonProto.h` — framing + codecs (single header).
- `include/TraceonServer.h` / `src/TraceonServer.cpp` — server.
- `include/TraceonClient.h` / `src/TraceonClient.cpp` — CRC-verified client.
- `benchmarks/remote_bench.cpp` — local/remote/serve benchmark.
- `tests/test_remote.cpp` — protocol + integration + concurrency tests.
