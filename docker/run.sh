#!/usr/bin/env bash
# TracEon remote-access testbed driver (v2.3.0 "Harpe").
#
# Builds the server + tools in Release (TRACEON_BUILD_SERVER=ON), generates a
# small current-development v5 cache from the bundled sample FASTA, serves it over a Docker bridge
# network, and runs the remote client benchmark against the server hostname.
# Prints both the remote numbers and an in-process (local) baseline.
#
# Requires Docker + Docker Compose v2.
set -euo pipefail
cd "$(dirname "$0")"

echo "==> [1/4] Building TracEon remote images (Release, TRACEON_BUILD_SERVER=ON)"
docker compose build

echo "==> [2/4] Building the cache once and starting the resident server"
trap 'echo "==> Tearing down"; docker compose down >/dev/null 2>&1 || true' EXIT
docker compose up --wait --wait-timeout 60 server

echo "==> [3/4] Checking independent query processes and exact output"
docker compose run --rm --no-deps client bash -ec '
  pids=""
  for worker in 1 2 3 4; do
    traceon query --host server --ids sample.ids > "/tmp/query-$worker.fasta" &
    pids="$pids $!"
  done
  for pid in $pids; do wait "$pid"; done
  for worker in 1 2 3 4; do cmp sample.fasta "/tmp/query-$worker.fasta"; done
  traceon query --host server --batch-size 64 --ids sample.ids > /tmp/batched.fasta
  cmp sample.fasta /tmp/batched.fasta
'
echo "==> [4/4] Benchmarking"
echo "--- remote (client container -> server over bridge network) ---"
docker compose run --rm --no-deps client
echo "--- local in-process baseline (same file, no network) ---"
docker compose run --rm --no-deps client-local

echo
echo "Done. Median remote latency is the round-trip from a client container to"
echo "the server container; local is the zero-copy getView() baseline."
