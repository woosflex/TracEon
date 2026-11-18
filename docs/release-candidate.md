# Local release candidate — 2026-09-22

This is an unreleased working-tree candidate, not a tagged or published release.
The full roadmap is not complete.

A fresh Release build in `/tmp/traceon-release-check-20260922` passed all 21
CTest entries, including remote protocol/CLI, persistence fault injection,
benchmark-output validation, and relocated installed-consumer checks. Only
pinned dependency sources were reused; no compiled dependency artifacts were
reused. Existing benchmark checks are now registered in the default Linux suite.

Candidate artifacts are in `build/release-candidate-2026-09-22/`:

| Artifact | SHA-256 |
|---|---|
| `source.tar.gz` | `265014f32c61c0afb700eca314ec03f768177da3969ae0b0223d829ab0d161a7` |
| `linux-local-install.tar.gz` | `e991a13b8f33e95c93cc0eeafdda58d3fdfa3f77d4ac8478d2444ab29fb2ce78` |

The source archive contains 463 working-tree source/documentation/data files,
including uncommitted changes, with individual hashes in `manifest.json`.
It excludes Git metadata, ignored build outputs, and fetched dependency sources.
CMake pins dependency revisions; a fresh download requires network access.
The archive predates this release note and is identified by its checksum rather
than a claim to represent all subsequent working-tree changes.

Every archived source and installed file was checked against the manifest.
The installed archive was extracted into a new directory; a downstream consumer
compiled and ran against that extracted package, and its CLI help check passed.
The binary package uses the host's system runtime libraries, including OpenSSL;
it is not claimed portable across Linux distributions or CPUs.

[Preserved manifests and verification logs](results/release-checks-2026-09-22.tar.gz)
record the source commit, dirty-tree state, per-file hashes, compiler, platform,
and test output. A separate [Docker check](results/docker-validation-2026-09-22.json)
validated four independent clients and batched exact output across a bridge.

Still required: consolidate the final source revision, execute hosted CI, repeat
relevant evidence against that revision, prepare publishable release artifacts,
and complete the manuscript and author declarations. No public visibility,
publication, or archival DOI is implied by this local packaging step.

Subsequent portability check: the build had forcibly enabled zlib-ng's
host-native instructions. The default is now OFF, and release/Docker builds
explicitly disable them. Reconfiguring the fresh build with this setting produced
227 compile commands with no `-march=native` or `-mcpu=native`, and all 21 tests
passed. The candidate archives above predate this fix; their checksums and
historical evidence remain unchanged.

The corrected current installation archive is
`build/release-portable-2026-09-22/traceon-Linux-x86_64.tar.gz`, SHA-256
`1e6cf6de5830748763171d0216b1e3d11b9cb678606af4baf10b5e4f25b5607e`.
Its checksum and every installed-file hash were verified, and a consumer built
and ran against the extracted archive. [Corrected-build verification logs](results/portable-release-checks-2026-09-22.tar.gz)
include the installation manifest. This removes forced host-native codec flags;
it does not remove system-runtime dependencies or establish cross-platform
compatibility. A final tagged revision and hosted checks remain required.

The corrected core-only configuration also built successfully with
`BUILD_TESTING=OFF`, `TRACEON_BUILD_SERVER=OFF`,
`CMAKE_DISABLE_FIND_PACKAGE_OpenSSL=TRUE`, and `WITH_NATIVE_INSTRUCTIONS=OFF`.
Its installation was relocated and a downstream consumer compiled and ran
successfully. This check covers users who need the local library without the
service or its OpenSSL dependency.

All six existing libFuzzer targets also passed bounded ASan/UBSan checks with
leak detection enabled: cache loader, gzip loader, FASTQ, FASTA, k-mer encoding,
and k-mer API. Each ran with an eight-second fuzz budget on copied corpora;
FASTA/FASTQ inputs included large seeds above the multithreaded parser threshold.
[Logs, executable hashes, commands, and execution counts](results/fuzz-verification-2026-09-22.tar.gz)
are preserved. An initial sandbox run failed in LeakSanitizer's ptrace-incompatible
shutdown; the same empty-input artifact and all six campaigns passed outside
that restriction. This is bounded testing, not proof that no defects remain.

The remote suite additionally passed under ThreadSanitizer: 629 assertions in
nine test cases, exit zero, with no race report. The Debug build used
`-fsanitize=thread` for compilation and linking and disabled host-native codec
instructions. [Build/test logs and executable hash](results/remote-tsan-2026-09-22.tar.gz)
are preserved. This result covers the remote suite, not every library operation.

The complete Docker wrapper subsequently passed with host-native tuning disabled.
It now starts the cache-builder dependency once, waits for server health, and
runs clients with `--no-deps`. Scalar/batch exact-output checks and both benchmark
commands completed; the dedicated test project's containers were removed. The
wrapper log is `build/docker-wrapper-validation-2026-09-22.log`. Its smoke timings
are not controlled comparative performance evidence. This wrapper correction
postdates the initial 475-file review snapshot.
