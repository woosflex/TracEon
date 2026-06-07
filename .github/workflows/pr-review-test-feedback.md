---
emoji: 🧪
description: Review pull requests, run tests, and post actionable feedback.
on:
  pull_request:
    types: [opened, synchronize, reopened, ready_for_review]
  workflow_dispatch:
permissions:
  contents: read
  issues: read
  pull-requests: read
  actions: read
  checks: read
tools:
  github:
    mode: gh-proxy
    toolsets: [default]
network:
  allowed:
    - defaults
    - threat-detection
safe-outputs:
  add-comment:
    max: 1
    hide-older-comments: true
---

# PR Review and Test Feedback

## Task

Review the triggering pull request, execute repository tests, and post one concise feedback comment.

## Instructions

1. Inspect the PR title, description, changed files, and discussion context.
2. Run tests exactly with:
   - `cmake -B build -S . -DCMAKE_BUILD_TYPE=Release`
   - `cmake --build build --target unit_tests -j 4`
   - `cd build && ctest --output-on-failure`
3. Analyze changed code for correctness risks, regressions, missing tests, and edge cases.
4. Post exactly one `add-comment` output containing:
   - overall assessment
   - test status and key output
   - prioritized actionable findings
   - clear next steps for the PR author
5. If no issues are found, still post a brief success summary and clearly state no blocking issues.
6. If tests fail or context is limited, explain the limitation and provide best-effort review feedback.

## Safe Outputs

- Use configured safe outputs for visible actions.
- Use `noop` only when no visible action is needed.
