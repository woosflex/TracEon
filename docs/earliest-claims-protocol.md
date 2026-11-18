# Earliest claim rerun — prospective protocol, 2026-09-22

Primary historical claim: 499d7f8 README reports50,000×1,000-base generated
FASTA: scan83.83ms, cache construction174.47ms, cached lookup0.0013ms,
>64,000×lookup speedup and amortization after3requests.

The original generated file/seed, wrapping, hardware details and compiler flags
are not preserved with that claim. Reconstruct deterministic ACGT input with
50,000 unique records of1,000bases; seed20260922, both single-line and80-column
wrapping. This is a reconstruction on the current local host, not an exact
reproduction of historical machine timing. Warm OS page cache; no cache eviction.

Archive exact source499d7f8 and build diagnostics. Use same-day fa95f1d follow-up
for runnable measurements; Cache.cpp and both benchmark drivers must be byte-
identical between commits. The follow-up repairs strategy implementation/linkage.
Record any modern compiler compatibility flags, including forced cstdint.
No current TracEon implementation is substituted. Two builds: O0 and O3/NDEBUG.

Run original unmodified drivers at first/middle/last IDs, six paired repetitions
per layout/build/target, alternating arm order3/3. Preserve all stdout/stderr,
reported timing and process results. The microsecond single-lookup figure is
sensitive to timer/cold-path effects; report ranges and medians, not just maxima.

Add a separate correctness-aware driver reusing the original full-file scanner
and historical Cache. Independently compute length-framed SHA256 over IDs and
all returned sequence bytes. Verify load and save/restore for all50,000 records.
Use O3 single-line input for six paired64-random-ID cohorts: full scans versus
cache versus reusable HTSlib1.24. Match owned-sequence semantics and common digest
consumption; setup separately, FAI construction recorded separately. This is a
stronger comparator control, not the original benchmark. All outputs must match.

Observe amortization directly on single-line input for both builds: six paired
runs each at1/2/3/4/8/16 requests. Each cache arm includes a fresh full parse and
cache construction; scanning repeats the original full-file search. Timing sums
setup plus verified retrieval, excluding process startup, ID loading and digest
initialization. Order balanced3/3 per count. No corpus-wide break-even claim.

Preserve source/binary/input/protocol hashes, compiler/host metadata, scripts,
raw measurements and all outcomes. Readiness smoke: compile both builds; verify
small plain and wrapped inputs against independently computed expected data
before timed full-scale runs. No network, cloud resources or production edits.
