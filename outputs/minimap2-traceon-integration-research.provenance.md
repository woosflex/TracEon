# Provenance: TracEon/minimap2 integration research

Canonical artifact: `outputs/minimap2-traceon-integration-research.md`

## Scope

Six-question design review of replacing minimap2's per-bucket `khash_t(idx)` table with TracEon's `traceon_kmer` C API. Sources were queried through current web search/fetch tools and local repository inspection. No production patch or benchmark result was claimed.

## Delegated evidence

- Q1/Q2: official minimap2 GitHub metadata, release v2.31/r1302, v2.26–v2.31 `index.c` refs, license, workflow/fork sources.
- Q3: canonical `index.c` operation inventory and local adapter mismatch audit.
- Q4/Q5: GCC/static-link evidence, C99 header check, mixed C/C++ example, minimap2 alternatives and generic hash benchmarks.
- Q6: official minimap2 Makefile/README/cookbook/CI, public NA12878/CHM13 sources, benchmark protocol.
- Independent reviewer pass identified and corrected: allocating `thread_local std::string` in `noexcept` error paths; incomplete freeze synchronization; overstatement of a two-probe upstream lookup; incomplete Makefile recipe; iterator alignment/lifetime concerns; partial `.mmi` failure handling; and unsupported maintainer-preference wording.

## Local commands and outputs

1. `gcc -std=c99 -pedantic-errors -Wall -Wextra -Werror -Iinclude -fsyntax-only /tmp/traceon_c99_probe.c` → `C99_HEADER_PASS`.
2. `make -C benchmarks/minimap2_src -B -n minimap2` → dry-run showed C compilation, C++ shim compilation, archive, and current GCC final link.
3. Strict local adapter compile at TracEon commit `9566b24d`, GCC 16.1.1:

   `cd benchmarks/minimap2_src && gcc -std=c99 -Wall -Wextra -Werror -DHAVE_KALLOC -I../../include -I../../build/_deps/ankerl_unordered_dense-src/include -c index.c -o /tmp/traceon-index.o`

   → nonzero; primary errors were iterator signature mismatch at `mm_idx_stat`, `mm_idx_cal_max_occ`, and `mm_idx_dump`; additional copied-source warnings were also promoted by `-Werror`.

4. `git ls-remote --tags https://github.com/lh3/minimap2.git` → v2.26–v2.31 SHAs; master points to v2.31 commit.
5. Official GitHub API repository contents query → no root `CONTRIBUTING.md` or AI policy file in the observed tree.

## Not run / unverified

- No corrected minimap2 build.
- No stock-versus-TracEon mapping run.
- No `.mmi` cross-compatibility test.
- No performance number, benchmark table, or human dataset download.
- No ThreadSanitizer run; the frozen-flag race is a source-level finding requiring a fix and test.
- No allocation-failure test for the diagnostic path; the report explicitly marks the current no-throw claim as incomplete.

## Key pinned sources

- https://github.com/lh3/minimap2/releases/tag/v2.31
- https://github.com/lh3/minimap2/commit/3c28777e7e2dcc90f825de1b9f17a89cca7d4452
- https://github.com/lh3/minimap2/blob/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/index.c
- https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/Makefile
- https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/README.md
- https://raw.githubusercontent.com/lh3/minimap2/3c28777e7e2dcc90f825de1b9f17a89cca7d4452/cookbook.md
- https://github.com/woosflex/TracEon/blob/9566b24dfee04e7d582fa082744cf37c8d59b99b/include/kmerindex_c_api.h
- https://github.com/woosflex/TracEon/blob/9566b24dfee04e7d582fa082744cf37c8d59b99b/src/kmerindex_c_api.cpp
