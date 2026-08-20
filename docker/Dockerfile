# syntax=docker/dockerfile:1
# TracEon remote-access testbed image (v2.3.0 "Harpe").
#
# Multi-stage: a builder stage compiles the remote server + tools in Release
# with TRACEON_BUILD_SERVER=ON; a slim runtime stage carries just the two
# binaries + the bundled sample FASTA. Build context must be the repo root
# (docker compose sets it to .. so third_party/lz4 and the sources resolve).
FROM ubuntu:24.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY . .
# Release is MANDATORY for meaningful numbers; the server/client libs are
# gated on TRACEON_BUILD_SERVER (default ON, kept explicit here).
RUN cmake -B build -DCMAKE_BUILD_TYPE=Release -DTRACEON_BUILD_SERVER=ON \
    && cmake --build build -j"$(nproc)" --target remote_bench traceon_make_cache

FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 ca-certificates bash \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/traceon
# Binaries land at the CMake build root (not in subdirectories).
COPY --from=builder /src/build/remote_bench ./remote_bench
COPY --from=builder /src/build/traceon_make_cache ./traceon_make_cache
COPY docker/data/sample.fasta ./sample.fasta

ENV PATH="/opt/traceon:${PATH}"
