# tests/AGENTS.md — TracEon Test Suite

## Purpose

Catch2-based unit test suite for TracEon. All tests must run from the `build/`
directory (relative-path file fixtures). The combined target is `unit_tests`
(128 test cases / 4195 assertions — the AGENTS.md root contract).

## Layout

`tests/SmartStrategyTests.cpp` was split into per-domain translation units so
each test executable is small enough to start under AddressSanitizer on hosts
where the ~42 MB combined `unit_tests` binary fails at ASan init
(`InternalMmapVector ENOMEM` — a binary-size/ASLR interaction). The combined
target links ALL split `.cpp` files together, preserving the 128/4195
contract.

| File | Domain | Tests |
|---|---|---|
| `TestHelpers.h` | Shared includes + inline fixtures (v4 blob crafting, MT fixture builders, gzip byte helpers) | — |
| `test_parser_fastq.cpp` | FASTQ framing, empty quality/seq, MT chunk classifier, bug1 (`@` in quality) / bugB (empty lines) repros, `+`-leading quality | 20 |
| `test_parser_fasta.cpp` | FASTA normalize, header-only, abutting header, terminator, CRLF, MT/prereserve | 24 |
| `test_gzip.cpp` | GZIP truncation, trailing garbage, Z_BUF_ERROR, parallel member decode, zip-bomb guard | 13 |
| `test_binary_cache.cpp` | `.traceon` v4 format, CRC32C, mutation/truncation hardening, count bounds, legacy rejection, smart compression | 26 |
| `test_lifecycle.cpp` | immutable-after-load, addEntry/set, clearCache, reader quiescence, failure atomicity, `[oom]` guard | 10 |
| `test_api_misc.cpp` | initial state, getView/size/keys, duplicate keys, SIMD, NGS mode, format detection | 20 |

The older `CacheTests.cpp` (4), `FastqTests.cpp` (1), and `MapDefsTests.cpp`
(10) remain as separate TUs inside `unit_tests`.

## Executables

- `unit_tests` — combined suite (`tests/main.cpp` provides `CATCH_CONFIG_MAIN`;
  links `Catch2::Catch2WithMain`). CTest name `TracEonTests`.
- `test_parser_fastq`, `test_parser_fasta`, `test_gzip`, `test_binary_cache`,
  `test_lifecycle`, `test_api_misc` — one per domain, each a self-contained
  Catch2 executable (`Catch2::Catch2WithMain`). All registered in CTest.

## Running

```sh
cd build
./unit_tests            # expect: All tests passed (4195 assertions in 128 test cases)
ctest -R "test_parser|test_gzip|test_binary_cache|test_lifecycle|test_api_misc"
```

### Under AddressSanitizer

The `[oom]` test (`setrlimit(RLIMIT_AS)`) is incompatible with ASan
shadow-mapping and aborts — a documented sanitizer limitation, not a defect.
Filter it with the Catch2 tag:

```sh
cd build-asan
./test_parser_fastq "~[oom]"
./test_parser_fasta "~[oom]"
./test_gzip "~[oom]"
./test_binary_cache "~[oom]"
./test_lifecycle "~[oom]"
./test_api_misc "~[oom]"
```

## Contracts

- Every `TEST_CASE` lives in exactly one translation unit (moved exactly once
  during the split).
- Do not add tests to `TestHelpers.h` — it holds helpers only.
- Keep each domain file small enough to start under ASan on this host
  (individual executables ~34–35 MB ASan; combined ~45 MB).
