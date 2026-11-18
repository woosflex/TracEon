# ADR-007: Portable `.traceon` v5 format

**Status:** Accepted for the 3.0.0-dev candidate
**Supersedes for writes:** [ADR-005 v4](ADR-005-traceon-v4-binary-format.md)
**Compatibility:** v4 remains readable; v1/v2/v3 remain unsupported
**Decision date:** not assigned
**Deciders:** [TODO: record maintainers who approved the release decision]

## Context

The v4 payload writes native-endian counts and lengths and, in NGS mode,
persists an implementation hash. C++ guarantees a `std::hash` result only for
the duration of one execution, and a hash-only map cannot retain two distinct
identifiers with the same value. v4 also lacked an exact detected-format byte
and accepted physical bytes after its declared frame. These properties prevent
a portable collision-safe interchange contract.

## Decision

`saveBinary()` writes v5. `loadBinary()` reads v4 and v5. Both readers require
the declared stored payload to end at physical EOF. v5 codec 0 adds explicit
mapped persistence; the one-argument writer retains codec 1 (LZ4). v4 is retained only as a
migration reader; this deliberately amends ADR-005's original permissive EOF
behavior.

The v5 header is 27 bytes:

| Offset | Size | Field |
|---:|---:|---|
| 0 | 4 | bytes `54 52 4f 05` (`TRO\x05`) |
| 4 | 1 | codec (`00` = raw mapped payload, `01` = LZ4 Frame) |
| 5 | 1 | index mode (`00` GENOME, `01` NGS) |
| 6 | 1 | exact `FileFormat` value |
| 7 | 8 | uncompressed payload length, little-endian |
| 15 | 8 | stored payload length, little-endian |
| 23 | 4 | CRC32C, little-endian |
| 27 | variable | raw payload or one LZ4 frame, ending at physical EOF |

`FileFormat` byte values are: `00` DNA_FASTA, `01` RNA_FASTA, `02`
PROTEIN_FASTA, `03` DNA_FASTQ, `04` RNA_FASTQ, `05` PROTEIN_FASTQ, and `ff`
UNKNOWN. Other values are invalid.

For codec 0 the stored and logical lengths must be identical. Other codec
values are rejected. Older development v5 readers reject codec 0 rather than
misinterpreting it; released v2.3.0 readers support neither v5 codec.

The uncompressed payload uses this field encoding in both index modes:

1. record count (`u64` little-endian);
2. for each record, identifier byte length (`u32` little-endian) and bytes;
3. sequence byte length (`u32` little-endian) and bytes; and
4. quality byte length (`u32` little-endian) and bytes.

Identifier, sequence, and quality are length-prefixed opaque byte strings; no
NUL terminator, text normalization, or character encoding conversion is added
by the container. Their interpretation follows the parser/API contract.
Duplicate identifiers make the payload invalid.

Hash-map record iteration order is unspecified. Therefore v5 has canonical
field and integer encoding, but two saves of the same logical snapshot are not
promised byte-identical. Deterministic record ordering is explicitly out of
scope.

No implementation hash is persisted. Both in-memory modes retain complete
identifier views as map keys, so normal hash-table equality resolves collisions.
The modes remain distinct API/format values for compatibility.

## Integrity checksum

CRC-32C uses the Castagnoli polynomial `0x1EDC6F41` (reflected polynomial
`0x82F63B78`), reflected input/output, initial value `0xffffffff`, and final XOR
`0xffffffff`. The running checksum covers the full uncompressed payload first,
then header bytes `[0,23)`; the checksum field is excluded. See
[`include/Crc32c.h`](../../include/Crc32c.h) and the CRC-32C definition used by
iSCSI in [RFC 3720 section 12.1](https://www.rfc-editor.org/rfc/rfc3720#section-12.1).
CRC32C detects accidental corruption; it does not authenticate the file.

## v4 migration behavior

The v4 reader retains legacy native-endian payload interpretation. For legacy
NGS payloads it consumes but does not trust the persisted hash and reconstructs
the current index from each stored identifier. CRC-valid duplicate identifiers
are rejected. A v4 frame must now reach physical EOF. Existing v4 caches can be
read and re-saved as v5 on their original compatible platform.

LZ4 restore temporarily maps the compressed file, allocates the full decoded arena,
decompresses and validates it, constructs the index, releases the mapping, and
only then publishes the snapshot. Streaming bounds additional codec buffers; it
does not bound total decoded memory.

Mapped restore validates the same CRC32C and parses the same checked portable
payload, then retains the read-only mapping as backing for identifier, sequence,
and quality views. The index remains private per process. Validation reads every
payload byte; this is neither instantaneous nor an O(1) restore. `getArenaBytes()`
counts the retained mapped file including its header, not physical/private memory.

The file must remain immutable while mapped. Save uses the existing atomic writer;
on POSIX replacement keeps existing readers on the previous inode. External in-place
writes or truncation violate the contract and can cause changed data or SIGBUS.
Reload, clear, and destruction retain the existing reader-quiescence requirement.

## Alternatives considered

- **Keep v4 writes:** rejected because native payload fields and persisted hashes
  cannot support the approved portability/collision contract.
- **Persist a newly specified numeric hash:** rejected because collision chains
  and a stable hash specification add complexity without removing identifiers
  from the payload.
- **Sort records for byte-deterministic output:** deferred; it requires extra
  memory/work and is unnecessary for logical interoperability.
- **Drop v4 reads:** rejected to retain a migration path.
- **Allow extension trailers:** rejected; future extensions require a new format
  version rather than ambiguous trailing data.

## Consequences

- Default writes are breaking and require the 3.0.0 development identity.
- Newly written caches preserve exact detected format and portable field values.
- v4 remains migration-only and is not retroactively portable.
- Appended physical bytes and duplicate identifiers are rejected.
- No crash-atomic replacement guarantee is added by this ADR.
