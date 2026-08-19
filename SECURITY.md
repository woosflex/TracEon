# Security Policy

## Reporting a Vulnerability

Please do **not** open a public issue for security vulnerabilities. Report
privately so the maintainer can address it before disclosure.

**Report via email:** adnanraza3435@gmail.com

Please include:
- The affected version(s) and platform(s).
- A description of the vulnerability and its impact.
- A minimal reproduction (input + steps) where possible.
- If known, a suggested fix.

You will receive an acknowledgement, and the maintainer will work with you
toward a coordinated disclosure. Do not disclose the issue publicly until a fix
is available, unless otherwise agreed.

## Scope

TracEon parses untrusted genomic data (FASTA/FASTQ, GZIP, and the `.traceon`
binary cache format). Anything that accepts untrusted input — parsers,
decompression, and cache deserialization — is the primary security surface.
Crashes, memory-safety issues, unbounded allocation, or OOM in these paths are
treated as security issues.

## Supported Versions

| Version       | Supported          |
| ------------- | ------------------ |
| latest release (v2.2.x) | ✅ |
| older releases | ❌ — upgrade recommended |

Orphaned/pre-release branches are not supported.
