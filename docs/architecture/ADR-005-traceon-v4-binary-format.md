# ADR-005: `.traceon` v4 Binary Format & CRC32C Integrity

**Date:** v2.0.0 (2026)
**Status:** Accepted (v2.0.0)
**Deciders:** Adnan Raza (Woosflex)
**Supersedes:** v3 (`TRO\x03`, LZ4 Frame, no checksum), v2 (`TRO\x02`, LZ4 block),
v1 (`TRO\x01`, uncompressed) — **all legacy readers removed in v2.0.0**.

## Context

The pre-v4 `.traceon` formats had no integrity protection: a truncated or
bit-flipped cache either decompressed to garbage or failed with an incidental
error, and the v2 int-truncation bugs (2^32+N) silently bypassed size checks.
The v2.0.0 design review (`outputs/traceon-v2-design-review.md`, Q1) mandated a
whole-payload **CRC32C** with explicit lengths and a required end condition, and
a clean break with **no legacy readers** (Q5).

## Decision

### Wire format (all multi-byte fields little-endian)

```
offset 0  magic          "TRO\x04"           (4 bytes)
offset 4  codec flags    0x01 = LZ4 Frame;   (1 byte)
                         bits 1..7 reserved, must be 0
offset 5  index mode     0 = GENOME, 1 = NGS (1 byte)
offset 6  logical length  uncompressed payload length (u64 LE)
offset 14 frame length    LZ4 Frame length    (u64 LE)
offset 22 CRC32C          Castagnoli, u32 LE  (4 bytes)
offset 26 LZ4 Frame       (frame length bytes)
```

Total header size: **26 bytes**. The frame is a standard LZ4 Frame written by
`LZ4F_compressBegin/Update/End` with a 1 MiB streaming window (bounded memory
regardless of dataset size).

### Checksum definition

CRC-32C (reflected in/out, poly `0x1EDC6F41`, init `0xFFFFFFFF`, final XOR
`0xFFFFFFFF` — the SCTP/iSCSI variant) over:

```
CRC32C( uncompressed_logical_payload || header[0..22) )
```

- The **entire** uncompressed logical payload, then the canonical header fields
  — everything up to but **excluding** the checksum field itself.
- The payload-first ordering is deliberate: `frame length` is only known once
  the frame is complete, so both save and load feed payload-then-header
  incrementally. Save updates the accumulator as serialized chunks pass through
  the LZ4F compressor; load updates it as decompressed chunks are written to
  `text_arena_`. **No second full-payload pass, no extra allocation.**
- This is accidental-corruption detection, **not authentication** (a deliberate
  attacker can construct collisions).

### Implementation & dispatch

`include/Crc32c.h` provides three byte-identical implementations:

| Path | Guard | Selection |
|---|---|---|
| x86 SSE4.2 `crc32` (`_mm_crc32_u*`) | `TRACEON_HAS_AVX2` + GCC/Clang, `target("sse4.2")` attribute | runtime `__builtin_cpu_supports("sse4.2")` |
| AArch64 `crc32cx/crc32cw/crc32ch/crc32cb` | `__ARM_FEATURE_CRC32` (AArch64 always defines `TRACEON_HAS_NEON`) | compile-time on AArch64 |
| portable byte-wise table | always | fallback |

Equivalence between fallback and hardware paths is tested over many lengths.

### Reader validation (hardening semantics preserved from v2/v3)

1. File ≥ 26 bytes, else "too small to contain header".
2. Magic `"TRO"` + version `0x04`; v1/v2/v3 → "unsupported cache version;
   regenerate with v2.0.0"; other magics → "Invalid binary cache magic bytes".
3. Codec flags must be exactly `0x01`; index mode must be 0 or 1.
4. `logical length < 8` (payload always starts with the 8-byte record count),
   empty frame, or `logical length > SIZE_MAX` → implausible, rejected
   **before allocation** (the v2 INT_MAX-truncation class cannot occur in the
   64-bit-clean LZ4F API; these checks plus the OOM guard are its equivalent).
5. `frame length > remaining bytes` → truncated (subtraction-form bounds check:
   `n > end - ptr`, no pointer-arithmetic overflow possible).
6. OOM guard (`logical length > available/2`) + `bad_alloc` catch around the
   arena resize.
7. Streaming decompression must consume exactly `frame length` bytes and end
   exactly at the LZ4 Frame terminator (`hint == 0`), and must produce exactly
   `logical length` bytes — truncation and wrong-length detection.
8. CRC32C over payload + header must equal the stored value — mutation
   detection.
9. Record count bounded by remaining payload bytes (≥ 12 B/record GENOME,
   ≥ 20 B NGS) before `reserve()`.
10. **Failure atomicity**: any failure resets the cache to the pristine empty
    state; `data_ready_` is never published for an invalid load.

### Payload layout (unchanged from v1–v3)

```
count (u64, native endian) | per record:
  GENOME: id_len(u32) id | seq_len(u32) seq | qual_len(u32) qual
  NGS:    hash(u64) | id_len(u32) id | seq_len(u32) seq | qual_len(u32) qual
```

## Consequences

### Positive

- Modified/truncated/crafted caches fail loudly with specific errors; partial
  data is never served.
- Streaming CRC adds no second pass and no peak-memory growth; hardware CRC is
  ~13 GB/s class on the reference CPUs.
- Removing v1/v2/v3 readers deletes the v2 int-truncation and v3 framing edge
  cases entirely.

### Negative

- v1/v2/v3 caches must be **regenerated** after upgrading (clean break; no
  public adopters identified at the time of the decision).
- 4 header bytes + CRC computation cost on every save/restore (negligible vs
  LZ4 itself).

## Alternatives Considered

- XXH64 (8 bytes, no hardware instruction; would couple the format to the
  vendored LZ4 tree's internal xxHash).
- BLAKE3 (cryptographic strength unnecessary for storage-corruption detection;
  larger field and dependency cost).
- Per-record checksums (4 B/record ≈ 26 MB for a 6.5M-record WGS cache;
  localization without recovery did not justify the size/complexity).
- Legacy readers kept for migration (rejected: clean-break decision, Q5).
