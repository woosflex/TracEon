# ADR-008: Snapshot identity and FASTQ access (protocol v2)

Date: 2026-09-21. Status: implemented in the unreleased 3.0.0 candidate.
Supersedes ADR-006 only where specified below; framing, GET semantics, CRC,
deadlines, and trusted-network deployment remain unchanged.

## Snapshot identity

Before listening, each server computes a SHA-256 digest of the logical cache.
Clients can pin this digest when connecting, so a different dataset with the same
IDs cannot silently supply the wrong records. Identity covers the complete served
snapshot, not just the requested subset. A restart or change of compression/index
mode preserves identity when logical content and detected format are identical.

The canonical digest input is, in this exact order:

1. ASCII bytes `TracEon snapshot v1` followed by LF (no NUL).
2. One detected-format byte, using `FileFormat` values (0–5, or 255 for UNKNOWN).
3. Record count as u64 little-endian.
4. Records sorted lexicographically by unsigned ID bytes. For each record, ID,
   sequence, and quality, each encoded as u64 little-endian byte length followed
   by its bytes. Empty fields have a zero length and no bytes.

Descriptors and original FASTQ plus lines are not cached and cannot affect the
digest. A cache assembled with UNKNOWN format is distinct from a parsed FASTA
cache, even if its sequences match. This digest is not the compressed file hash.
SHA-256 is implemented with OpenSSL EVP (Crypto library, remote target only).
Hashing streams one record at a time after copying/sorting the keys; temporary
memory includes all copied IDs and one copied quality. Include this startup cost
in deployment/amortization results. No digest work occurs on the query path.

The caller must keep the cache immutable while the server runs, including during
startup hashing. The metadata is stable for that server run. Neither the digest
nor CRC authenticates an untrusted peer; pin a digest obtained from a trusted
launch or manifest, rather than trusting discovery from an unknown endpoint.

## Protocol versions and capabilities

HELLO still contains `TRO-RPC` and a version byte. Servers accept versions 1 and 2;
all others receive BAD_REQUEST and close. V1 has its original empty OK handshake
payload and sequence-only GET/BYE operations. Removed opcodes 0x03/0x04 remain
reserved. GET and its NOT_FOUND/empty-sequence behavior do not change.

V2 HELLO returns the existing CRC-protected OK envelope with a 45-byte payload:

| Field | Encoding |
|---|---|
| Negotiated version | u8, exactly 2 |
| Capabilities | u32 LE; bit 0 GET, bit 1 GET_FASTQ, bit 2 GET_BATCH |
| Record count | u64 LE |
| Canonical snapshot SHA-256 | 32 raw digest bytes |

GET is always supported; GET_FASTQ is advertised only for detected FASTQ formats.
Clients ignore unknown capability bits but require GET and the exact v2 metadata
shape/version. Public APIs expose the digest as 64 lowercase hexadecimal digits.

New clients default to v2. They reject an empty v1 response, invalid metadata,
checksum failure, or an expected-digest mismatch, and close before returning a
usable connection. There is no automatic downgrade. The optional final C++
constructor argument / CLI `query --protocol 1` explicitly supports old v1
servers, with no metadata or FASTQ. Snapshot pinning is incompatible with v1.

## Complete FASTQ access

GET_FASTQ uses the new opcode 0x06 and the same raw-key request body as GET.
It is legal only on v2 connections to a FASTQ snapshot. Success uses the existing
CRC-protected OK envelope, containing:

```
sequence_len(u32 LE) | quality_len(u32 LE) | sequence | quality
```

Both lengths must match, and the payload must contain exactly those bytes.
The maximum sequence length is floor((64 MiB - 18) / 2), accounting for both
payloads and all framing within the existing payload-length cap. Servers check
sequence size before constructing a response and validate quality length;
clients validate both lengths and CRC before returning an owned record.

NOT_FOUND stays distinct from a present empty record (two zero lengths).
Unsupported operations and malformed frames receive BAD_REQUEST and close.
Local client capability failures send nothing and leave the connection usable.
The CLI emits four-line FASTQ using the queried ID and a bare `+` line. It preserves
sequence/quality bytes, request order, and duplicate requests; original header
and plus-line descriptions are intentionally unavailable.

## Bounded batch access

GET_BATCH (0x07) is an additive v2 capability. Requests contain a FASTQ flag
(u8, 0 or 1), a count (u32 LE, 1–256), then repeated key lengths (u32 LE)
and key bytes. FASTQ batches also require the GET_FASTQ capability.

The CRC-protected OK payload contains the count (u32 LE), followed by each
record's status (u8), sequence length (u32 LE), quality length (u32 LE), sequence,
and quality. Entries preserve request order and duplicates. NOT_FOUND entries
have zero lengths; present empty records use OK. Sequence-only entries have no
quality bytes; FASTQ lengths must match. Trailing bytes are invalid.

Both directions retain the 64 MiB frame limit. Oversized responses fail the
whole batch with BAD_REQUEST and close the connection; callers must reduce the
batch size. The client uses one operation deadline for the complete exchange.
Unsupported capabilities are rejected locally without closing the connection.

`getBatch(keys, fastq)` returns ordered optional owned records. CLI batching is
opt-in through `query --batch-size 2..256`; the default remains scalar access.
The CLI stops output at the first missing ID with exit status 3, though later
IDs in that same batch have already been retrieved. No automatic fallback or
retry occurs. Batch latency is latency per batch, not per-record tail latency.

## Verification

`test_remote` covers identity after restore/reordering, quality-only changes,
wrong expected digests, v1 compatibility/rejection of v2-only operations, complete
and empty FASTQ records, misses, malformed lengths, metadata, CRC, and transport
failure closure. `service_cli` independently recomputes the canonical digest with
Python hashlib and tests pinned independent processes after deleting cache/source
files, CLI metadata, FASTQ output, and rejection before output on mismatch.

OpenSSL API reference: [EVP digest operations](https://docs.openssl.org/3.0/man3/EVP_DigestInit/).
