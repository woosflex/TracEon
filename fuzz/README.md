# TracEon Fuzzing (v2.1.0 workstream)

OSS-Fuzz-style libFuzzer targets closing the last never-run verification
gate: the loader/parser/kmer code paths under sanitizers with arbitrary
bytes.

## Targets (each a self-contained `fuzz/*.cpp`)

| Target            | Code under test                                                                 |
|-------------------|----------------------------------------------------------------------------------|
| `fuzz_v4_loader`  | `.traceon` v4 binary loader — `SmartStrategy::loadBinary` / `Cache::restore`     |
| `fuzz_gzip_loader`| GZIP load path — `loadGzipFile` + `loadFile` extension/magic detection           |
| `fuzz_fastq`      | FASTQ text parsing (`loadFile`, GENOME + NGS modes)                              |
| `fuzz_fasta`      | FASTA text parsing + in-place `normalizeFastaArena` (GENOME + NGS modes)         |
| `fuzz_kmer_encode`| `encode_kmer` (vuln-0001 contract: no OOB read, no bad shift)                   |
| `fuzz_kmer_api`   | kmerindex C API (vuln-0003: no exception/UB across the C boundary)              |
| `fuzz_trki`       | `KmerReferenceIndex::mmap_open` (vuln-0006: crafted TRKI must reject, no SEGV)   |

Every target uses the standard entry point:

```cpp
int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);
```

plus a portable file-fed `main` guarded by
`#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION` (in `fuzz_common.h`):

- **clang `-fsanitize=fuzzer`** defines `FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION`
  and links libFuzzer's own driver (the harness is just the target function).
- **GCC / plain builds** compile the guarded main, which reads each argv file
  and feeds it to the target, exiting 0 on all inputs. This lets the exact
  same sources be built and exercised on clang-free hosts (this repo's dev
  box is GCC-only) with the seed corpus as plain files:

  ```bash
  g++ -std=c++20 -O1 -g -include cstdint -I include -I fuzz \
      fuzz/fuzz_v4_loader.cpp build/libtraceon_core.a build/libtraceon_kmer.a \
      build/third_party/lz4/build/cmake/liblz4.a build/_deps/zlib-ng-build/libz.a \
      -lpthread -o /tmp/fuzz_v4_loader
  /tmp/fuzz_v4_loader fuzz/corpus/v4_loader/*   # must exit 0
  ```

The loaders under test take file paths, so the path-based harnesses write
each input to a per-pid scratch dir (`/tmp/traceon_fuzz_<pid>`) before
calling the API.

## Build (clang, sanitizers)

```bash
cmake -B build-fuzz -DTRACEON_BUILD_FUZZERS=ON \
      -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=RelWithDebInfo \
      -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-fuzz -j 4 --target fuzz_v4_loader fuzz_gzip_loader \
      fuzz_fastq fuzz_fasta fuzz_kmer_encode fuzz_kmer_api fuzz_trki
./build-fuzz/fuzz_v4_loader -max_total_time=120 fuzz/corpus/v4_loader/
```

- Requires clang (`-fsanitize=fuzzer` is clang-only); CMake enforces this.
- The core/kmer libraries are instrumented too (whole-program ASan/UBSan).
- Fuzzers are NOT registered with CTest; CI runs them with explicit budgets.
- The normal build (`TRACEON_BUILD_FUZZERS=OFF`, the default) never compiles
  `fuzz/*.cpp`.

CI: `.github/workflows/fuzz.yml` — push/PR to main + weekly schedule.

## Seed corpus (`fuzz/corpus/<target>/`)

Reproducible: `fuzz/tools/gen_corpus.py` regenerates every non-v4 seed and
`fuzz/tools/gen_v4_corpus.cpp` regenerates the v4 loader seeds (it crafts
LZ4-compressed payloads + CRC32C, so it links the project's lz4:

```bash
ANKERL=$(find build build-fuzz -path '*ankerl_unordered_dense-src/include' 2>/dev/null | head -1)
g++ -std=c++20 -O2 -I include -I third_party/lz4/lib -I "$ANKERL" \
    fuzz/tools/gen_v4_corpus.cpp build/libtraceon_core.a \
    build/third_party/lz4/build/cmake/liblz4.a \
    build/_deps/zlib-ng-build/libz.a -lpthread -o /tmp/gen_v4_corpus \
    && /tmp/gen_v4_corpus fuzz/corpus/v4_loader
```

Contents per target:

- **v4_loader/** — valid GENOME + NGS round-trip caches, plus the crafted
  malformed variants from `tests/test_binary_cache.cpp`: legacy v1/v2/v3
  magic, unknown version, bad codec flags, invalid mode, truncated header /
  frame / tail, wrong logical length (2^64-1, <8), 2^64-1 frame length,
  mutated CRC / header bytes, count bombs (5e8, 3-with-1-record,
  2-with-junk), old `MMAP` format, empty frame.
- **gzip_loader/** — valid small gzip FASTQ + FASTA, truncated member,
  trailing garbage, concatenated valid+truncated, non-gzip text, binary
  garbage, empty file.
- **fastq/** — `test_data/simple.fastq` + crafted framing variants: empty
  quality, empty sequence, '+'-leading quality, '@' in quality, CRLF,
  missing trailing newline, partial records, binary garbage.
- **fasta/** — `test_data/simple.fasta` + multi-line, CRLF, header-only,
  abutting-header, missing trailing newline, garbage.
- **kmer_encode/** — `(int32 k, seq bytes)` inputs incl. k=0/33, short views,
  non-ACGT bases.
- **kmer_api/** — empty / tiny / random bytes driving
  create→reserve→insert→get→iterate→destroy.
- **trki/** — valid small index + every vuln-0006 PoC: `n_slots=2^60`
  56-byte evil file, out-of-range slot offsets, truncated/oversized, bad
  magic/version/k, non-pow2 n_slots, n_keys>n_slots, empty 40-byte index.

The >10MB multithreaded-parser seeds are generated at CI time (not committed
— they'd bloat the repo); see the "Generate large corpus seeds" step in
`fuzz.yml`.
