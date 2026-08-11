# TracEon v2.0.0 Design Review

**Scope:** the current `/home/adnanraza/agent_workspace/TracEon` checkout, the maintainer's stated profile, and primary/official sources checked during this review.

**Release premise:** Per the maintainer's stated project context, no public adopters have been identified, so v2.0.0 may make a clean break. This does **not** make the changes SemVer-minor; it makes the migration cost acceptable. The current branch is `fix/p0-integrity` with uncommitted changes; `main` is tagged v1.3.0, and `dev` is an ancestor of the current fix branch. No code was changed by this review and no benchmark/test result below is presented as newly run unless explicitly marked.

## Decision summary

| Question | v2.0.0 decision |
|---|---|
| Q1 checksum | Mandatory whole-logical-payload **CRC32C**, with explicit lengths and a required end condition; no per-record checksums in the default format. |
| Q2 clear/reload | Keep **documented external reader quiescence** for the existing naked `std::string_view` API; retain write-side locking for writer serialization. Do not ship a version-check-only or mutex-only reclamation fix. |
| Q3 k-mer C API | **Defer to v2.1.0 or a separately versioned experimental target.** |
| Q4 index performance | **Freeze the current measured wins** for v2.0.0; optimize later behind a benchmark and memory/correctness gate. |
| Q5 version | **v2.0.0 is correct.** Fold unreleased v1.4.0 work into it; do not create a public v1.4.0 tag first. |
| Q6 Git/release | Split the stack into 5–6 logical commits, fast-forward the already-ancestral `dev`, cut `release/2.0.0`, merge to `main`, then tag the tested merge commit. |

---

## Q1. Checksum for `.traceon` v4

### Recommendation

**Use CRC32C as the mandatory whole-payload integrity check, calculated over canonical header fields plus the uncompressed logical payload; use length fields and complete-frame/end-of-file validation to detect truncation. Do not add per-record checksums by default.**

### Evidence

- The official xxHash benchmark table reports, on its stated Intel i7-9700K reference system, approximately **19.4 GB/s for XXH64**, **13.0 GB/s for CRC32C with SSE4.2**, and **4.4 GB/s for BLAKE3 with AVX2**. These are directional implementation benchmarks, not TracEon measurements. The table also makes an important distinction: the high SIMD numbers belong to **XXH3**, not classic XXH64. [xxHash documentation](https://xxhash.com/doc/v0.8.2/index.html)
- CRC32C has direct hardware support on the target families: Intel documents CRC32/SSE4.2 intrinsics, and AArch64 documents CRC32C instructions operating on 8/16/32/64-bit operands. This gives CRC32C a small, mature streaming implementation with scalar fallback. [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html), [Arm CRC32C instruction](https://developer.arm.com/documentation/ddi0602/2026-06/Base-Instructions/CRC32CB--CRC32CH--CRC32CW--CRC32CX--CRC32C-checksum-)
- CRC32C is a 32-bit accidental-corruption checksum: its field costs 4 bytes and the simple random-mismatch model gives an approximately `2^-32` undetected probability. XXH64 costs 8 bytes and gives a larger approximately `2^-64` random-mismatch space, but remains a non-cryptographic hash and lacks the same classic CRC error-detection specification. BLAKE3's default output is 32 bytes and is cryptographic, but that strength is unnecessary for ordinary storage corruption and increases format/dependency cost. [IETF RFC 4960](https://www.rfc-editor.org/rfc/rfc4960.html), [xxHash](https://github.com/Cyan4973/xxHash), [BLAKE3](https://github.com/BLAKE3-team/BLAKE3)
- None of the three algorithms knows that a missing tail was supposed to exist. A checksum is useful only if the reader also requires the declared logical length, compressed length/frame termination, and a complete checksum-bearing end condition. A truncated file that also loses its footer cannot be distinguished from a valid prefix by a digest alone. This is a format/framing requirement, not an algorithm feature.

### Trade-offs

- CRC32C gives up the stronger random-collision margin of XXH64 and the adversarial/tamper resistance of BLAKE3. It must not be described as authentication; a deliberate attacker can construct CRC collisions.
- CRC32C's hardware acceleration is not guaranteed on every supported CPU. Provide a portable table/folding fallback and runtime dispatch. Do not publish the xxHash reference benchmark as a TracEon result.
- XXH64 may have a superficially low integration cost because the vendored LZ4 tree contains xxHash sources, but relying on that transitive/internal API couples the `.traceon` format to a compression dependency. Classic XXH64 also should not inherit XXH3's SIMD claims.
- A whole-payload checksum detects that the cache is invalid but does not identify the first bad record. Per-record checksums would improve localization and partial-recovery diagnostics, at the cost of 4 bytes per record, more checksum calls, and more format complexity. For the repository's record-dense WGS profile, the historical v1.4 benchmark table reports about 6.5 million records for a 1 GB synthetic input ([`benchmarks/v1.4.0_benchmark_report.md`, WGS table](file:///home/adnanraza/agent_workspace/TracEon/benchmarks/v1.4.0_benchmark_report.md)); 4 bytes per record would be roughly 26 MB of checksum fields before alignment/other overhead. This is a derived estimate, not a v2 measurement.
- A file-level BLAKE3-256 sidecar remains a reasonable future option for content identity, hostile transport, or reproducible artifact addressing. It should not be the hot-path default merely because it is stronger.

### Integration notes

- Replace the current v3 magic/reader with a v4 header in [`src/SmartStrategy.cpp`](file:///home/adnanraza/agent_workspace/TracEon/src/SmartStrategy.cpp). A concrete layout should include: `TRO\x04`, format/codec flags, index mode, logical uncompressed payload length, compressed frame length or exact frame boundary, and a checksum field with an explicitly specified byte order and CRC32C parameters.
- Update `serializePayload()`/`saveBinary()` so the CRC state is updated as serialized chunks are passed to LZ4F. Update `loadBinary()` while decompressed chunks are written to `text_arena_`. This avoids a second full-payload pass and preserves the current bounded streaming design.
- Check the canonical header fields that affect interpretation plus the logical payload; exclude the checksum field itself. Require exact decompressed length, exact frame completion, and checksum equality before publishing `data_ready_`.
- Do not checksum every FASTA/FASTQ record in v2.0.0. If later requirements need corruption localization or partial restore, add fixed-size chunk checksums, preferably around the existing 1 MiB streaming boundary, as a separately measured format extension.
- Add known-answer tests, single-byte mutations at header/middle/tail, truncation at every frame boundary, wrong logical length, wrong record count, and a portable-fallback-vs-hardware-dispatch equivalence test. The current repository's `loadBinary()` tests cover v1/v2/v3 behavior and must be replaced with v4-only tests under the no-legacy policy.

---

## Q2. `clearCache()`/reload while a reader holds the old arena

### Recommendation

**Choose (a): retain the documented external-quiescence contract for the current naked `std::string_view` API, with a debug-only lifecycle diagnostic; keep the write-side mutex to serialize writers but do not claim it protects lock-free readers.**

### Evidence

- In the current implementation, `getView()` acquire-loads `data_ready_` and then reads `file_cache_` without taking `cache_mutex_` ([`src/SmartStrategy.cpp:1652–1667`](file:///home/adnanraza/agent_workspace/TracEon/src/SmartStrategy.cpp#L1652-L1667)), while `clearInternal()` clears the map and then releases `text_arena_`, `manual_store_`, and the mmap ([`src/SmartStrategy.cpp:89–112`](file:///home/adnanraza/agent_workspace/TracEon/src/SmartStrategy.cpp#L89-L112)). The release/acquire flag publishes initialization; it does not wait for a reader that already observed `true` and is not a reclamation mechanism. This supports the proposed external-quiescence contract; it is not evidence that concurrent clear/reload is currently safe. [`ADR-001`](file:///home/adnanraza/agent_workspace/TracEon/docs/architecture/ADR-001-lock-free-reads.md)
- `std::string_view` is non-owning. Even if a reader-side version check observes the same generation before and after lookup, the caller can retain the returned view and dereference it after the writer destroys the old arena. A version check detects some races during lookup; it cannot extend the lifetime of bytes after return. [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines.html)
- A real epoch/RCU design requires a grace period: removal/publication is separated from reclamation until all pre-existing readers have exited. Linux's RCU documentation states this explicitly. For TracEon, the protected object would need to be the complete immutable snapshot—map plus arena/mmap ownership—not just a version integer. [Linux RCU concepts](https://docs.kernel.org/RCU/rcu.html)
- An atomic `shared_ptr<const Snapshot>` could provide ownership, but C++ specifies atomicity without promising universal lock-freedom; the specialization's `is_always_lock_free` is implementation-defined. It also requires returning an owner-bearing `SnapshotView` or handle rather than discarding ownership and returning only a naked `string_view`. [C++ draft atomic shared_ptr](https://eel.is/c++draft/util.smartptr.atomic.shared)

### Trade-offs

- The library gives up transparent hot reload while a view may still be in use. Applications must stop readers before `clearCache()`, reload, or destruction. This is a real ergonomics cost, but it is honest and preserves the core promise: lock-free reads against an immutable, live snapshot.
- A debug assertion cannot prove that an arbitrary caller will not use a returned `string_view` later. A debug reader counter can catch overlap with the lookup call if readers use a documented scope, but it cannot make the current API safe after the lookup returns. Treat it as a diagnostic, not synchronization.
- Option (b) is the right future architecture only if the API changes to an epoch-bearing `SnapshotView`, callback-scoped read API, or owner-carrying handle. A bare versioned arena is insufficient.
- Option (c), a write-side mutex alone, serializes `clearCache()` against other writers but does not exclude current lock-free `getView()` calls and does not protect a view after the mutex is released.

### Integration notes

- Update [`include/Cache.h`](file:///home/adnanraza/agent_workspace/TracEon/include/Cache.h), [`include/SmartStrategy.h`](file:///home/adnanraza/agent_workspace/TracEon/include/SmartStrategy.h), README, and ADR-001 with the precise rule: reads are safe concurrently only while the same loaded snapshot remains installed; a returned view is invalid once `clearCache()`, reload, or destruction begins.
- Keep `cache_mutex_` around writer operations and build-phase mutation. Do not add a shared lock to `getView()` as a purported solution: it would only protect the lookup call, not later use of the returned view.
- Add a debug-only `TRACEON_DEBUG_LIFECYCLE` counter/scope if useful for catching coordinated lookup overlap, but state its limitation in the code comment. Add sequential clear/reload and use-after-clear documentation tests; do not add a test that treats concurrent clear plus retained `string_view` as supported.
- If hot reload becomes a product requirement, make it a separate API design: immutable `Snapshot`, atomic publication, retired-snapshot reclamation, and `SnapshotView` ownership. Do not retrofit an epoch counter while retaining the current return type.

---

## Q3. Experimental k-mer C API

### Recommendation

**Defer the experimental k-mer C API to v2.1.0 (or a separately versioned experimental target); do not include it in the v2.0.0 stable core release.**

### Evidence

- The files are **ignored/untracked WIP in this checkout**: [`include/kmerindex_c_api.h`](file:///home/adnanraza/agent_workspace/TracEon/include/kmerindex_c_api.h), [`src/kmerindex_c_api.cpp`](file:///home/adnanraza/agent_workspace/TracEon/src/kmerindex_c_api.cpp), [`include/KmerIndex.h`](file:///home/adnanraza/agent_workspace/TracEon/include/KmerIndex.h), and [`include/KmerReferenceIndex.h`](file:///home/adnanraza/agent_workspace/TracEon/include/KmerReferenceIndex.h). They are not in `CORE_SOURCES`/test sources ([`CMakeLists.txt`, lines 35–48](file:///home/adnanraza/agent_workspace/TracEon/CMakeLists.txt)) or public umbrella wiring. Shipping them would be a new ABI/build contract, not simply packaging existing v2 work.
- The C shim has the security-review hazard identified by the maintainer: exported `extern "C"` functions can call `new`, `reserve`, and `try_emplace` without translating exceptions to C status codes. A C caller has no C++ exception contract. The implementation also lacks general null-handle/output-pointer validation. Microsoft documents the need for an explicit error-code boundary when C++ exception-based code interfaces with non-exceptional code. [Microsoft exception-boundary guidance](https://learn.microsoft.com/en-us/cpp/cpp/how-to-interface-between-exceptional-and-non-exceptional-code?view=msvc-170), [C++ language linkage](https://en.cppreference.com/w/cpp/language/language_linkage)
- `kmerindex_iter_begin()`/`iter_next()` mutate iterator state through `const_cast`; the cursor is stored in the shared index object, so iteration is not reentrant and is not safe for concurrent readers. `kmerindex_get()` returns an interior pointer into the map while the API still exposes insertion/reserve operations that may rehash. The header's “no inserts after build” comment is a caller convention, not an enforced phase transition.
- `KmerReferenceIndex` introduces a separate `TRKI` mmap format ([`include/KmerReferenceIndex.h`](file:///home/adnanraza/agent_workspace/TracEon/include/KmerReferenceIndex.h), header/`mmap_open()` implementation). It therefore has its own bounds, overflow, endianness, word-size, mapping-lifetime, and corruption-test obligations. Combining it with proposed `.traceon` v4 would multiply release-surface risk while the core cache is already changing format and lifecycle behavior.

### Trade-offs

- v2.0.0 will not provide the minimap2 integration or persistent k-mer reference index. That delays a potentially valuable bioinformatics use case.
- Deferral costs another release cycle and may require a separate ABI/versioning decision, but it prevents an unbuilt, untested C boundary from becoming part of the stable contract.
- Keeping it out means v2.0.0 can make a clean break in the cache core without also freezing minimap2-specific key semantics, iterator behavior, and `TRKI` persistence.

### Integration notes

- Leave the WIP files out of `TracEon.h`, the installed public target, and v2.0.0 release notes. Keep the work on `feature/kmer-c-api` or a dedicated experimental target rather than silently shipping ignored files.
- Before v2.1.0, define a C error model: every exported function must be `noexcept` at the ABI boundary, catch all C++ exceptions, return status codes, and expose `kmerindex_last_error()` or caller-provided error storage. Define null behavior explicitly.
- Replace shared mutable iteration state with caller-owned iterator state or a callback/visitor API. Enforce a build/freeze transition before returning interior pointers, or return values/copies where pointer stability cannot be guaranteed.
- Add a C-compiled minimap2 integration target, ABI compile tests, malformed-input tests, sanitizer/UBSan/ASan tests, and a benchmark that measures the actual cross-language call boundary.

---

## Q4. Further index optimization before release

### Recommendation

**Freeze performance at the current validated wins for v2.0.0; do not ship a custom probe strategy, map replacement, or merge redesign without a separate benchmark-gated change.**

### Evidence

- The maintainer's current profile attributes approximately 65% of a 1000 MB load to `ankerl::unordered_dense` emplacement/content hashing/open-addressing insertion; GZIP wall time is flat and appears inflate/memory-bound on the hybrid-core Meteor Lake host. This is the review request's current measurement, not an independently rerun result.
- The current branch already includes the high-leverage, low-surface-area changes: `std::string_view` keys, wyhash-style transparent hashing, and pre-reservation in [`MapDefs.h`](file:///home/adnanraza/agent_workspace/TracEon/include/MapDefs.h) and [`SmartStrategy.cpp`](file:///home/adnanraza/agent_workspace/TracEon/src/SmartStrategy.cpp). The maintainer reports approximately +33% insert throughput and -14% to -28% parse-load time for these wins, but those exact deltas were not found in the checked benchmark reports; treat them as input observations, not independently verified release evidence, until a raw benchmark artifact is attached.
- Amdahl's law puts a ceiling on the likely release benefit. If 65% is the fixed insertion phase, a 25% insertion speedup yields about **17.6% faster total load** (`1/(0.35+0.65/1.25)`); a 50% insertion speedup yields about **21.7% faster total load** (`1/(0.35+0.65/1.5)`); eliminating that phase entirely would cap total speedup at about 2.86×. Therefore a risky map rewrite should not be accepted without measuring a material end-to-end gain, not just a faster microbenchmark.
- The existing multithreaded parser creates thread-local maps and then inserts their entries into the destination map. A merge redesign could be promising, but it changes ownership, duplicate semantics, string-view lifetime, reserve behavior, and serialization order. It is not a safe “small optimization” inside an integrity-focused release.

### Trade-offs

- v2.0.0 leaves measurable insertion time on the table. The 65% profile means a future optimization could matter, especially for record-dense FASTQ.
- Freezing avoids changing a third-party hash-table contract, increasing memory unexpectedly, introducing rehash/pointer bugs, or invalidating the zero-copy/immutable-after-load proof immediately before release.
- The decision is not “never optimize.” It is “do not combine speculative optimization with a new binary format, stricter corruption handling, and lifecycle-contract changes.”

### Integration notes

- Add the final Release benchmark as a v2.0.0 gate: parse/load time, per-stage timing if available, restore time, lookup throughput, peak RSS, record count, and checksum overhead on FASTA/FASTQ/GZIP and both index modes. Pin CPU affinity or at least record hybrid-core placement, CPU model, RAM, compiler, flags, OS, input hashes, and repetitions.
- Run the project's regression script from a clean Release build. The repository mandates Release-mode validation in [`AGENTS.md`](file:///home/adnanraza/agent_workspace/TracEon/AGENTS.md); current historical reports are not proof for the final dirty tree.
- Open a post-v2 optimization branch with one hypothesis per change: destination-map merge, alternative map, or custom flat table. Require correctness, duplicate-key semantics, `string_view` lifetime, save/restore, RSS, and end-to-end measurements before accepting it.
- If a future map/layout change alters public headers, serialized representation, or view lifetime, make it a major-format/API decision rather than hiding it under a patch release.

---

## Q5. Is v2.0.0 the right version?

### Recommendation

**Yes. Release this stack as v2.0.0, with `.traceon` v4 as the only readable binary format and an explicit migration note that all old caches must be regenerated. Do not publish an intermediate v1.4.0 tag.**

### Evidence

- SemVer states that the major component **MUST** increase for backward-incompatible changes, while minor and patch releases are for backward-compatible additions and fixes. The requested behavior changes are externally observable even where function signatures remain unchanged. [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html)
- The proposed v2.0.0 changes include: corrupt/truncated GZIP now throws instead of serving partial data; malformed binary caches are rejected; `set()`/`addEntry()` after load throws `std::logic_error`; the key representation changes to non-owning `string_view`; and the binary format changes from v3 to proposed v4 with no legacy reader. These are precisely the sort of behavioral, ABI/source, and data-format boundaries that deserve a major release.
- Git currently has a public `v1.3.0` tag, while the checkout has v1.4.0 metadata in [`include/TracEon.h`](file:///home/adnanraza/agent_workspace/TracEon/include/TracEon.h), README/release documents, and an `[Unreleased]` P0 section in [`CHANGELOG.md`](file:///home/adnanraza/agent_workspace/TracEon/CHANGELOG.md), but no observed v1.4.0 tag. Treating v1.4.0 as released would create a provenance gap.
- The no-adopter premise reduces migration cost but does not change the semantic meaning of the version number. v2.0.0 communicates that old source assumptions and old `.traceon` files cannot be assumed to work.

### Trade-offs

- Existing v1.3.0 users, if any appear later, must recompile and regenerate caches. Old `.traceon` files will not be readable; this is intentional.
- A major version may look heavier than a bug-fix release for the parser hardening, but hiding a format and behavior break under 1.x would create a misleading compatibility promise.
- Removing legacy readers simplifies `loadBinary()` and its attack surface, but eliminates a convenient rolling migration path. The release notes must say this plainly.

### Integration notes

- Set `TracEon::VERSION` and public documentation to `2.0.0`; choose a new codename or explicitly retire “Caliburn.” Remove stale v1.4.0 claims and dates. Update any generated/package metadata and add a CI check that the version header and tag agree.
- Change the magic to `TRO\x04`; make `loadBinary()` reject v1/v2/v3 with a clear “unsupported pre-v2 cache; regenerate” error. Remove legacy parsing branches and legacy-compatibility tests, replacing them with v4 round-trip and rejection tests.
- The CHANGELOG release entry should have a prominent **Breaking changes / migration** section before Added/Changed/Fixed. Keep an `Unreleased` section above it for post-release fixes, following Keep a Changelog's guidance that changelogs are curated for humans and should make breaking changes obvious. [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
- Explicitly document the `getView()` lifetime rule from Q2, the post-load mutation exception, strict GZIP rejection, malformed-cache rejection, cache regeneration command/path, C++20 requirement, and the fact that the experimental k-mer C API is not part of v2.0.0.

### Required breaking-section contents

1. **Binary cache:** `.traceon` v4 only; v1/v2/v3 are rejected; regenerate caches after upgrading.
2. **Input integrity:** truncated/corrupt GZIP and trailing garbage now throw `std::runtime_error`; partial data is never published.
3. **Mutation/lifecycle:** `set()`/`addEntry()` after load throws `std::logic_error`; call `clearCache()` only after all readers stop using views.
4. **View/storage:** `getView()` remains zero-copy and non-owning; its result must not be used after clear/reload/destruction.
5. **Source/ABI:** `std::string_view`-backed index keys and any advanced-header layout changes require recompilation; no ABI stability is promised across 1.x→2.x.
6. **Deferred surface:** k-mer C API and `TRKI` are experimental and not shipped.

---

## Q6. Git flow and release mechanics

### Recommendation

**Use 5–6 logical commits on the existing fix branch, fast-forward the already-ancestral `dev`, cut a release branch, merge the tested release commit to `main`, and place annotated `v2.0.0` on that exact merge commit.**

### Evidence

- The observed topology is favorable: `main`/`origin/main` point at v1.3.0, `dev` points at the unreleased v1.4.0 work, and `dev` is already an ancestor of `fix/p0-integrity`. Therefore a second merge/cherry-pick of dev is unnecessary and risks duplicate or conflicting history.
- The current fix branch is dirty with tracked source/header/test/docs changes and untracked artifacts such as benchmark scripts, `outputs/`, and `strix_runs/`. A release tag must identify committed source and documentation, not an implicit working-tree stack.
- Keep a Changelog explicitly distinguishes a curated release changelog from a noisy commit log. Logical commits should aid review and bisectability; the final release notes should summarize user-visible changes rather than reproduce commit history. [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)

### Trade-offs

- Five or six commits require deliberate staging because the current hardening changes are interrelated. A single mega-commit is easier to create but harder to review, bisect, or revert.
- A merge commit on `main` preserves the release boundary but produces a non-linear history. Fast-forwarding `dev` first keeps the already-linear ancestry; the release merge can then be the one explicit promotion point.
- Cleaning stale v1.4.0 documents during the release adds documentation work, but leaving them inconsistent is worse than the extra commit.

### Integration notes: concrete sequence

1. **Inventory and preserve the dirty worktree.** On `fix/p0-integrity`, review `git status --short`, `git diff --stat`, and `git diff --check`. Keep `outputs/`, `strix_runs/`, generated binaries, benchmark data, and unrelated scripts out of the release unless they are deliberately tracked release artifacts. Do not reset the branch before the stack is safely committed or otherwise backed up.
2. **Split the stack into 5–6 logical commits** (exact names may vary):
   - `fix(parser): harden FASTA/FASTQ framing and gzip integrity`
   - `fix(binary): add v4 framing/checksum and reject malformed payloads`
   - `fix(lifecycle): enforce immutable-after-load and clear contract`
   - `perf(index): zero-copy keys, hashing, and reservation`
   - `test: cover integrity, lifecycle, v4 round-trip, and rejection cases`
   - `docs(release): reconcile v2.0.0 metadata and migration notes`
   Keep tests with their implementation when that improves review; six is a target, not a requirement to manufacture artificial boundaries.
3. **Integrate the current fix branch without replaying dev.** After the commits, verify `git merge-base --is-ancestor dev fix/p0-integrity`. Then fast-forward `dev` to the fix branch (`git switch dev && git merge --ff-only fix/p0-integrity`). If dev has moved in the meantime, use a normal reviewed merge or rebase rather than blindly cherry-picking the same commits twice.
4. **Cut the release branch only after v4 code/tests/docs are present.** `git switch -c release/2.0.0 dev`. On this branch, freeze features. Resolve all version strings, remove v1.4.0 release claims, finalize `CHANGELOG.md`, and ensure the release notes contain the breaking section above.
5. **Run the release gates from a clean Release build:** `cmake -B build-release -S . -DCMAKE_BUILD_TYPE=Release`, `cmake --build build-release -j`, `ctest --test-dir build-release --output-on-failure`, `git diff --check`, sanitizer/fuzz malformed-cache tests, and the tracked benchmark/regression workflow from `build-release`/`benchmarks` as required by the project contracts. The repository's CHANGELOG reports 115 test cases/4076 assertions, but this review did not rerun them; regenerate the count/result on the final v2 tree.
6. **Promote and tag.** Merge `release/2.0.0` to `main` with the repository's normal reviewed PR policy (a non-fast-forward merge is useful for a visible release boundary). Create an annotated tag `v2.0.0` on the exact tested `main` commit, then push the branch and tag. Never move an existing release tag.
7. **Back-merge release-only fixes.** Merge the release branch or any release-only documentation/version corrections back into `dev` after tagging. Record the final commit, tag, test command, benchmark environment, and known residual risks in the release notes.

### Release gate checklist

- [ ] No uncommitted release-relevant changes; generated artifacts excluded intentionally.
- [ ] `TRO\x04` save/restore round trip works; v1/v2/v3 rejection is tested.
- [ ] CRC32C known-answer, mutation, truncation, length, header, and fallback-dispatch tests pass.
- [ ] GZIP truncation/trailing-garbage rejection remains intact.
- [ ] `set()`/`addEntry()` immutable-after-load behavior and clear/reload precondition are documented and tested.
- [ ] Full Release build/tests pass; sanitizer and malformed-input tests run.
- [ ] Final performance profile is regenerated on the actual release commit; no unsupported “verified” claim is copied from the dirty v1.4.0 documents.
- [ ] `include/TracEon.h`, README, CHANGELOG, release notes, and tag all say 2.0.0.
- [ ] K-mer C API is absent from the stable target and explicitly listed as deferred.

---

## Prioritized v2.0.0 engineering sprint

### P0 — integrity and format correctness

1. Freeze the v4 wire specification: header fields, endian rules, codec/frame lengths, CRC32C parameters, checksum coverage, complete-file/truncation rule, and error taxonomy.
2. Implement streaming CRC32C with x86 SSE4.2/AArch64 dispatch plus a tested portable fallback; update save/load without a second full-payload allocation.
3. Remove v1/v2/v3 readers and legacy tests under the explicit no-backward-compatibility decision; add v4 corruption, truncation, wrong-length, wrong-mode, and mutation tests.
4. Re-audit all parser and binary-loader hardening on the final v4 path, including failure atomicity: invalid loads must not publish partial state; define which malformed FASTA/FASTQ cases warn-and-skip versus throw.

### P1 — lifecycle/API contract

5. Make the external reader-quiescence rule prominent in `Cache.h`, `SmartStrategy.h`, README, and ADR-001; keep the writer mutex but do not market it as reader reclamation.
6. Add debug-only lifecycle diagnostics and sequential clear/reload tests; explicitly test/document that retained views become invalid at the lifecycle boundary.
7. Reconcile the `std::string_view` key ownership proof, map-clear-before-arena-release ordering, and all `IndexMode` save/restore paths.

### P1 — release hygiene and verification

8. Split and commit the current stack; exclude WIP k-mer files and unrelated generated artifacts.
9. Update all version/codename/magic comments and release documentation to 2.0.0, with a prominent migration section.
10. Run clean Release build, full Catch2 suite, sanitizer/UBSan/ASan checks, malformed-cache fuzz/property tests, and benchmark/regression gates. Record hardware and exact commands.

### P2 — deliberate follow-up, not release blockers

11. Defer the k-mer C API until its exception boundary, null/error semantics, iterator ownership, pointer stability, C-compiled integration, and `TRKI` format are designed and tested.
12. Benchmark one index optimization hypothesis at a time after the release: destination merge, alternative map, or custom flat table. Accept only changes that improve end-to-end load/lookup while preserving RSS, integrity, serialization, and zero-copy lifetime invariants.

## Sources

### Primary/official external sources

- [Semantic Versioning 2.0.0](https://semver.org/spec/v2.0.0.html)
- [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/)
- [xxHash documentation and benchmark table](https://xxhash.com/doc/v0.8.2/index.html)
- [xxHash official repository](https://github.com/Cyan4973/xxHash)
- [BLAKE3 official repository](https://github.com/BLAKE3-team/BLAKE3)
- [Intel Intrinsics Guide](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
- [Intel CRC32C API documentation](https://www.intel.com/content/www/us/en/docs/ipp/developer-guide-reference/2026-0/crc32-crc32c.html)
- [Arm A64 CRC32C instruction](https://developer.arm.com/documentation/ddi0602/2026-06/Base-Instructions/CRC32CB--CRC32CH--CRC32CW--CRC32CX--CRC32C-checksum-)
- [IETF RFC 4960, SCTP CRC32c](https://www.rfc-editor.org/rfc/rfc4960.html)
- [Linux kernel RCU concepts](https://docs.kernel.org/RCU/rcu.html)
- [C++ draft atomic shared_ptr specialization](https://eel.is/c++draft/util.smartptr.atomic.shared)
- [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines.html)
- [Microsoft: interfacing exceptional and non-exceptional code](https://learn.microsoft.com/en-us/cpp/cpp/how-to-interface-between-exceptional-and-non-exceptional-code?view=msvc-170)
- [cppreference: language linkage](https://en.cppreference.com/w/cpp/language/language_linkage)

### Local primary sources inspected

- [`AGENTS.md`](file:///home/adnanraza/agent_workspace/TracEon/AGENTS.md)
- [`CHANGELOG.md`](file:///home/adnanraza/agent_workspace/TracEon/CHANGELOG.md)
- [`CMakeLists.txt`](file:///home/adnanraza/agent_workspace/TracEon/CMakeLists.txt)
- [`include/Cache.h`](file:///home/adnanraza/agent_workspace/TracEon/include/Cache.h)
- [`include/SmartStrategy.h`](file:///home/adnanraza/agent_workspace/TracEon/include/SmartStrategy.h)
- [`include/MapDefs.h`](file:///home/adnanraza/agent_workspace/TracEon/include/MapDefs.h)
- [`include/TracEon.h`](file:///home/adnanraza/agent_workspace/TracEon/include/TracEon.h)
- [`src/SmartStrategy.cpp`](file:///home/adnanraza/agent_workspace/TracEon/src/SmartStrategy.cpp)
- [`src/kmerindex_c_api.cpp`](file:///home/adnanraza/agent_workspace/TracEon/src/kmerindex_c_api.cpp)
- [`include/kmerindex_c_api.h`](file:///home/adnanraza/agent_workspace/TracEon/include/kmerindex_c_api.h)
- [`include/KmerReferenceIndex.h`](file:///home/adnanraza/agent_workspace/TracEon/include/KmerReferenceIndex.h)
- [`docs/architecture/ADR-001-lock-free-reads.md`](file:///home/adnanraza/agent_workspace/TracEon/docs/architecture/ADR-001-lock-free-reads.md)
- [`benchmarks/real_world_benchmark_report.md`](file:///home/adnanraza/agent_workspace/TracEon/benchmarks/real_world_benchmark_report.md)
- [`benchmarks/v1.4.0_benchmark_report.md`](file:///home/adnanraza/agent_workspace/TracEon/benchmarks/v1.4.0_benchmark_report.md)
- [`tests/SmartStrategyTests.cpp`](file:///home/adnanraza/agent_workspace/TracEon/tests/SmartStrategyTests.cpp)

### Verification status

- **Performed:** repository status/history inspection, source/header/CMake inspection, official-source retrieval, and transparent Amdahl arithmetic for the stated 65% profile.
- **Not performed:** no code edits, no clean Release rebuild, no full test run, no sanitizer/fuzz run, no checksum microbenchmark on Meteor Lake, and no concurrent clear/reload stress test. The report therefore recommends these as release gates rather than claiming they passed.
