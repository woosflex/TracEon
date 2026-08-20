#!/usr/bin/env bash
# TracEon remote-access testbed driver (v2.3.0 "Harpe").
#
# Builds the server + tools in Release (TRACEON_BUILD_SERVER=ON), generates a
# small v4 cache from the bundled sample FASTA, serves it over a Docker bridge
# network, and runs the remote client benchmark against the server hostname.
# Prints both the remote numbers and an in-process (local) baseline.
#
# Requires Docker + Docker Compose v2.
set -euo pipefail
cd "$(dirname "$0")"

echo "==> [1/4] Building TracEon remote images (Release, TRACEON_BUILD_SERVER=ON)"
docker compose build

echo "==> [2/4] Generating v4 cache from docker/data/sample.fasta"
docker compose run --rm make-cache

echo "==> [3/4] Starting server on the bridge network (hostname: server:9876)"
docker compose up -d server
trap 'echo "==> Tearing down"; docker compose down >/dev/null 2>&1 || true' EXIT

# Give the healthcheck a moment so the client's depends_on:service_healthy is clean.
echo "==> [4/4] Benchmarking"
echo "--- remote (client container -> server over bridge network) ---"
docker compose run --rm client
echo "--- local in-process baseline (same file, no network) ---"
docker compose run --rm client-local

echo
echo "Done. Median remote latency is the round-trip from a client container to"
echo "the server container; local is the zero-copy getView() baseline."
