# Contributing to TracEon

Thanks for considering a contribution. TracEon is a high-performance, zero-copy,
lock-free C++20 genomic data cache, and the bar for what lands in `main` is
deliberately high: architectural invariants, performance, and byte-identical
correctness are non-negotiable.

## Code of Conduct

All participants must follow the [Code of Conduct](CODE_OF_CONDUCT.md). Be
professional and assume good intent — feedback is about the code, not about you.

## Developer Certificate of Origin (DCO)

All contributions must certify that you have the right to submit the code under
the project's MIT license. By contributing you agree to the
[Developer Certificate of Origin](https://developercertificate.org/) (v1.1):

```
Developer Certificate of Origin
Version 1.1

By making a contribution to this project, I certify that:

(a) The contribution was created in whole or in part by me and I have the right
    to submit it under the open source license indicated in the file; or

(b) The contribution is based upon previous work, to the best of my knowledge,
    covered under an appropriate open source license and I have the right under
    that license to submit that work with modifications, ...; or

(c) The contribution was provided directly to me by some other person who
    certified (a) or (b) and I have not modified it.

(d) I understand and agree that this project and the contribution are public
    and that a record of the contribution (including all personal information I
    submit with it, including my sign-off) is maintained indefinitely and may be
    redistributed consistent with this project or the open source license(s)
    involved.
```

Sign off every commit with `Signed-off-by: Your Name <you@example.com>` (use
`git commit -s`). Commits without a sign-off may be rejected.

## How to Contribute

1. **Fork** the repository and clone your fork.
2. **Create a topic branch** with a descriptive name:
   - `fix/...` for bug fixes, `feat/...` for features, `perf/...` for performance.
3. Make a **small, focused** change — one PR = one logical change. Small PRs
   review faster.
4. **Run the tests** (see Requirements below) and **add tests for your change**.
5. **Commit** with a clear message (see Commit Conventions) and push to your fork.
6. **Open the PR** early if you want early feedback; mark it `[WIP]` if unfinished.

## Requirements

- **Build:** Release mode is mandatory for any measurement. `cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j$(nproc)`.
- **Tests:** New features must include unit tests (Catch2, in `tests/`).
  Run from `build/`: `ctest` (or `./unit_tests`).
- **Style:** Modern C++20 (RAII, `std::string_view`, no raw `new` in core paths).
- **Performance / invariants:** If you touch core logic, run the regression check:
  ```bash
  python benchmarks/check_regression.py --baseline main
  ```
  and confirm you haven't broken an architectural invariant (see AGENTS.md /
  `docs/architecture/` ADRs): zero-copy reads, immutable-after-load, lock-free
  read path, pre-reserved maps.

## Commit Conventions

- Title: concise, under 50 chars, imperative mood ("Fix X", not "Fixed X").
- Description: wrap at 72 chars; explain the **why** ("So that we...") and
  reference the issue/PR (e.g. `Closes: #15`).
- Group related changes into one commit; keep unrelated changes separate.
- Prefix with the subsystem when it adds context (e.g. `kmer:`, `docs:`, `fuzz:`).

## Pull Request Checklist

- [ ] Tests added/updated and passing (`ctest`)
- [ ] Regression check run if core logic changed
- [ ] Docs/ADRs updated if behavior or architecture changed
- [ ] Commit(s) signed off (`git commit -s`)
- [ ] Branch rebased/tidied; avoid force-pushing after review begins

## Reporting Issues

Please include:
- **System details** (OS, CPU, compiler + version).
- A minimal **reproducible example**.
- Steps to reproduce and expected vs actual behavior.

Security issues: do **not** open a public issue — see [SECURITY.md](SECURITY.md).
