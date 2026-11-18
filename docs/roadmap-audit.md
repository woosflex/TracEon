# Roadmap completion audit — 2026-09-22

The full roadmap is **not complete**. This audit preserves the original four
remaining requirements; successful local checks do not establish release or
submission readiness.

| Requirement | Evidence inspected | Finding |
|---|---|---|
| Measure and remove demonstrated query overhead | Scalar profiling, bounded batch protocol/CLI tests, balanced indexed/scalar/batch studies, batch latency archive | Implemented and measured in development builds. Batching remains optional; long-read latency rises with concurrency. Final-revision reproduction remains open. |
| Controlled multi-process evidence | SIRV, short-read, and long-read manifests, raw measurements, payload/SAM checks, paired variation, setup projections, service resources | Development evidence exists with explicit CPU, page-cache, memory, and sampling limits. It does not establish a general mapping speedup. No final tagged-revision evidence exists yet. |
| Reproducible release | Fresh non-native Release build, 21 passing tests, core-only relocated consumer, corrected installation archive, Docker wrapper, six ASan/UBSan fuzz targets, remote TSan suite | Local verification passed. Working tree remains uncommitted; hosted multi-platform CI and final release artifacts are not verified at an approved revision. |

The artifact audit read eight evidence archives containing 925 file members and
checked successful server exits where recorded in summary objects. This verifies artifact
availability and readable records; it does not independently reproduce every
experiment or certify the manuscripts' scientific adequacy.

## Decisions and external evidence still needed

- The user explicitly requested leaving all changes uncommitted for review.
  Final-revision hosted CI and publication therefore remain deferred; do not
  interpret goal continuation as permission to commit or publish.
- The author confirmed Claude Code/Codex scaffolding and manual review of current
  code and drafts. The disclosure also records the broader assistance performed
  in this session. Exact historical model/version identifiers and competing
  interests still need confirmation. Saveena Solanki's recorded code contributions
  are acknowledged; sole authorship, ORCID, no affiliation, and no funding stand.
- The repository is currently private. Its earlier public-development interval
  remains unverified; commit timestamps alone do not prove public visibility.
- Hosted CI, public release, public research references, and archival publication
  remain future actions. They must not be represented as accomplished by local
  packaging, template changes, or draft prose.

The initial review snapshot is preserved in
`build/review-ready-2026-09-22/`; the subsequent Docker wrapper correction and
its release-note update are explicitly documented. Refresh the review snapshot
before committing the final approved candidate rather than silently treating the
older snapshot as identical to the current worktree.
