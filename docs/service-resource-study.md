# Sustained service resource checks, 2026-09-21

The installed development service completed three cycles of connection churn,
overload, idle peers, slow readers, maximum-size partial frames, and abrupt
connection termination. This is a bounded local stress experiment, not evidence
for exposing the unauthenticated protocol to an untrusted network.

```bash
python3 benchmarks/service_resources.py \
  --traceon build/research-install/bin/traceon \
  --connections 4 --cycles 3 --phase-seconds 2 \
  --output build/service-resources-final-2026-09-21
```

The fixture contains a four-base record and an 8 MiB record. The admission cap
was four, with a 250 ms frame deadline. Idle and slow peers remain open past
the deadline; a fresh snapshot-pinned CLI query must succeed before the test
closes them. Partial-frame peers advertise the legal 64 MiB maximum and then
stall. Other peers disconnect during headers/HELLO or advertise an oversized
frame that must be rejected. Under deliberate overload, admitted clients keep
returning the expected sequence while excess connections close before HELLO.

| Observation | Result |
|---|---:|
| Resource samples (10 ms target interval) | 4,142 |
| Recorded connection attempts | 205,501 |
| Deliberate overload rejections | 25,072 |
| Successful raw sequence queries | 179,517 |
| Peak threads | 6 (main + accept + 4 workers) |
| Idle threads after every phase | 2 |
| Idle open descriptors after every phase | 6 |
| Peak server PSS | 266.62 MiB |
| Peak server + load-generator PSS | 287.36 MiB |
| Server PSS after first complete cycle | 10.61 MiB |
| Post-warmup retained-PSS growth | 0 KiB |

The server started at 9.96 MiB PSS. After the first cycle its PSS was 10,861 KiB;
the next two cycles also ended at 10,861 KiB. A bounded amount of retained
allocator state is expected; there was no accumulating per-connection retention
in this run. Every phase returned to the baseline thread/descriptor counts, and
the process exited successfully after SIGTERM.

Memory is sampled, so very brief peaks can be missed. The aggregate includes the
Python load generator/controller; it is not a comparison with a mapping worker.
The 64 MiB frame cap applies per connection, so a larger admission cap permits a
larger memory envelope, including temporary framing copies. The observed peak is
not a universal memory limit or a guarantee for all payload sizes/platforms.

Raw samples and phase counts are in
`build/service-resources-final-2026-09-21/resources.jsonl` and `summary.json`.
The tested installed binary SHA-256 is
`dae5ca6ea5f89d5f6ea61c32009eadb82eafe2da9f80db8a21e65720dcb76b4a`.
This was an unreleased working tree. Repeat the experiment at the final release
revision; do not silently attribute these observations to a later binary.

The [evidence bundle](results/development-study-2026-09-21.tar.gz) includes the
final raw samples and the harness sources matching their recorded hashes.
