#!/usr/bin/env python3
"""Regenerate the non-v4 TracEon fuzz seed corpus.

Reproduces every seed under fuzz/corpus/ EXCEPT the v4 loader seeds (which
need real LZ4-Frame compression + CRC32C — see gen_v4_corpus.cpp).

Usage:  python3 fuzz/tools/gen_corpus.py [corpus_root]
Default corpus_root: fuzz/corpus (run from the repo root).

Deterministic: all random content uses a fixed seed.
"""
import gzip
import os
import random
import shutil
import struct
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CORPUS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "fuzz", "corpus")

rng = random.Random(0x5EED)

def ensure(d):
    os.makedirs(d, exist_ok=True)

def write(path, data):
    with open(path, "wb") as f:
        f.write(data)

def le32(v):
    return struct.pack("<I", v)

def le64(v):
    return struct.pack("<Q", v)

def main():
    # ── kmer_encode/ (int32 LE k + sequence bytes) ───────────────────────────
    d = os.path.join(CORPUS, "kmer_encode")
    ensure(d)
    def kmer_seed(name, k, seq):
        write(os.path.join(d, name), struct.pack("<i", k) + seq)
    kmer_seed("k4_acgt.bin", 4, b"ACGTACGTACGTACGTACGT")
    kmer_seed("k32_acgt.bin", 32, b"ACGT" * 8)
    kmer_seed("k1_a.bin", 1, b"A")
    kmer_seed("k0.bin", 0, b"ACGT")                    # invalid k
    kmer_seed("k33.bin", 33, b"ACGT" * 9)              # invalid k
    kmer_seed("k_neg1.bin", -1, b"ACGT")               # invalid k
    kmer_seed("k_intmax.bin", 0x7FFFFFFF, b"ACGT")     # absurd k
    kmer_seed("k5_short.bin", 5, b"ACG")               # short view
    kmer_seed("k4_nonacgt.bin", 4, b"ACGTNNNN")        # N is invalid
    kmer_seed("k4_empty.bin", 4, b"")                  # empty view
    kmer_seed("k0_empty.bin", 0, b"")
    kmer_seed("k4_garbage.bin", 4, bytes(rng.randrange(256) for _ in range(64)))

    # ── kmer_api/ (arbitrary bytes) ──────────────────────────────────────────
    d = os.path.join(CORPUS, "kmer_api")
    ensure(d)
    write(os.path.join(d, "empty.bin"), b"")
    write(os.path.join(d, "zeros.bin"), b"\x00" * 64)
    write(os.path.join(d, "ff.bin"), b"\xFF" * 64)
    write(os.path.join(d, "small.bin"), bytes(range(16)))
    write(os.path.join(d, "random_4k.bin"),
          bytes(rng.randrange(256) for _ in range(4096)))

    # ── fastq/ ───────────────────────────────────────────────────────────────
    d = os.path.join(CORPUS, "fastq")
    ensure(d)
    shutil.copyfile(os.path.join(ROOT, "test_data", "simple.fastq"),
                    os.path.join(d, "simple.fastq"))
    write(os.path.join(d, "empty_quality.fastq"),
          b"@r1\nACGT\n+\n\n@r2\nTGCA\n+\nIIII\n")
    write(os.path.join(d, "empty_seq.fastq"),
          b"@r1\n\n+\nIIII\n@r2\nACGT\n+\n####\n")
    write(os.path.join(d, "plus_leading_quality.fastq"),
          b"@r1\nACGT\n+\n++II\n@r2\nTGCA\n+\n####\n")
    write(os.path.join(d, "at_in_quality.fastq"),
          b"@r1\nACGT\n+\n@@@D\n@r2\nTGCA\n+\nIIII\n")
    write(os.path.join(d, "crlf.fastq"),
          b"@r1\r\nACGT\r\n+\r\nIIII\r\n@r2\r\nTGCA\r\n+\r\n####\r\n")
    write(os.path.join(d, "no_trailing_newline.fastq"),
          b"@r1\nACGT\n+\nIIII")
    write(os.path.join(d, "partial.fastq"), b"@r1\nACGT\n+")
    write(os.path.join(d, "garbage.fastq"),
          bytes(rng.randrange(32, 127) for _ in range(512)))
    # Gzip-magic-routing seed for loadFile detection.
    with open(os.path.join(ROOT, "test_data", "simple.fastq"), "rb") as f:
        plain = f.read()
    with open(os.path.join(d, "gzip_routed.fastq.gz"), "wb") as f:
        f.write(gzip.compress(plain, mtime=0))

    # ── fasta/ ───────────────────────────────────────────────────────────────
    d = os.path.join(CORPUS, "fasta")
    ensure(d)
    shutil.copyfile(os.path.join(ROOT, "test_data", "simple.fasta"),
                    os.path.join(d, "simple.fasta"))
    write(os.path.join(d, "multiline.fasta"),
          b">chr1\nACGT\nACGT\nACGT\n>chr2\nTGCA\n")
    write(os.path.join(d, "crlf.fasta"), b">chr1\r\nACGT\r\n>chr2\r\nTGCA\r\n")
    write(os.path.join(d, "header_only.fasta"), b">only_header\n>next\nACGT\n")
    write(os.path.join(d, "abutting_header.fasta"), b">a\nACGT>b\nTGCA\n")
    write(os.path.join(d, "no_trailing_newline.fasta"), b">chr1\nACGTACGT")
    write(os.path.join(d, "partial.fasta"), b">chr1\nACG")
    write(os.path.join(d, "garbage.fasta"),
          bytes(rng.randrange(32, 127) for _ in range(512)))

    # ── gzip_loader/ ─────────────────────────────────────────────────────────
    d = os.path.join(CORPUS, "gzip_loader")
    ensure(d)
    def gz(path, payload, mtime=0):
        with open(path, "wb") as f:
            f.write(gzip.compress(payload, mtime=mtime))
    gz(os.path.join(d, "valid_fastq.fastq.gz"), plain)
    with open(os.path.join(ROOT, "test_data", "simple.fasta"), "rb") as f:
        plain_fa = f.read()
    gz(os.path.join(d, "valid_fasta.fasta.gz"), plain_fa)
    valid = gzip.compress(plain, mtime=0)
    write(os.path.join(d, "truncated.fastq.gz"), valid[: len(valid) // 2])
    write(os.path.join(d, "trailing_garbage.fastq.gz"), valid + b"GARBAGE")
    write(os.path.join(d, "not_gzip.gz"),
          b"hello world, this is not gzip data at all")
    write(os.path.join(d, "garbage_binary.gz"),
          bytes(rng.randrange(256) for _ in range(1024)))
    write(os.path.join(d, "empty.gz"), b"")
    # Concatenated: valid member + truncated tail member.
    write(os.path.join(d, "concat_valid_truncated.gz"),
          valid + valid[: len(valid) // 3])
    # Concatenated: two valid members (single-stream path via gzread).
    write(os.path.join(d, "concat_two_valid.gz"), valid + valid)

    print(f"corpus regenerated under {CORPUS}")


if __name__ == "__main__":
    main()
