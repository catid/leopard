# Default-off GF16 split-butterfly cache-block candidate

Bead: `leopard-79h.38.5.4.9`. Base commit: `2aca7f5`.
Date: 2026-09-06. Status: correctness-validated experiment; **performance
unknown, not promoted**. The first diagnostic attempt failed its isolation
gate before timing the candidate. No historical exact-Leopard1 gap is closed.

## Mechanism and cost hypothesis

`LEO2_EXPERIMENT_GF16_SPLIT_CACHE_BLOCK=ON` changes only the AVX-512VL nibble
backend's non-fused GF16 radix-four range callbacks above 8192 bytes per row.
Each independent four-row group completes both existing split layers over
8-KiB segments before advancing its byte offset. The final compact ALTMAP
tile stays intact because every preceding segment is a multiple of 64 bytes.

At the historical 32-KiB execution-row size, the four-row working set goes
from 128 KiB to 32 KiB per segment. The local Threadripper 9980X reports
48 KiB of L1 data cache per core. Reduced L2/L1 traffic is a hypothesis, not a
measured cache-miss result: associativity, other resident data, and compiler
code shape still matter.

This is not the previously rejected smaller **end-to-end** encoding tile,
register-fused radix-four payload, or pre-expanded whole-transform executor.
Public scratch geometry, source-copy policy, transform stages, field products,
and arithmetic load/store counts are unchanged. For a nonzero-skew group,
four split calls become sixteen at 32-KiB rows; extra table preparation and
call overhead are explicit risks. This candidate does not optimize the final
split-XOR accumulation callback or the field driver's distance-one stages.

GCC 13.3 Release disassembly contains the 0x2000 chunk loop, calls to the
existing two-way kernels, and a remainder path in both range callbacks. The
GF16-only AVX-512 object grows from 22,252 to 24,140 text bytes (+1,888);
data/BSS remain 416/8 bytes. No allocation, mutable global state, or new ISA
requirement is added to the hot path. Pure AVX2, scalar, SSSE3, GFNI, and
default-off AVX-512 routing are unchanged.

## Current-host scope correction

The July 0.976778x result is historical, not a current-HEAD measurement.
Current production already has a narrowly gated AUTO GF16/GFNI encoder for
AMD family 1Ah/model 08h, K1000/R200, 64-KiB full-output, native-layout,
single-thread calls. The hookless `leopard2_auto_gf16_gfni_production` test
actually passed on this 9980X; it did not skip. This cache-block candidate
does **not** change that GFNI path. Future comparisons must identify the
selected operation backend, including explicit AVX-512 and other-host AUTO
paths, before attributing the old deficit or a new improvement.

## Deterministic validation

The new production-linked executable has two serial CTests: `_ranges` and
`_encoder`. Both pass in these GF16-only builds:

| Build | Range cases | Full encoder shapes | Peak test cgroup bytes |
| --- | ---: | ---: | ---: |
| GCC Release, experiment ON | 2,048 | 6 | 144,048,128 |
| GCC Release, default OFF | 2,048 | 6 | 144,392,192 |
| GCC ASan+UBSan+LSan, ON, `-O1 -g1` | 2,048 | 6 | 175,370,240 |

Every final scope had all six `memory.events` counters zero and zero swap.
Configure/build jobs used the canonical lock, one build job, 512-MiB limit,
and no swap. Tests used 256 MiB and no swap. The initial full sanitizer build
peaked at 324,304,896 bytes; the default configure/build at 221,233,152 bytes.
Release uses `-Werror -Wall -Wextra`; sanitizer failures are non-recovering.
The sanitizer quarantine is 8 MiB, with 64-KiB thread-local quarantine.

Range coverage includes both directions, both fused hints, all eight
zero-skew combinations, distances 1/4/16/64, and 16 lengths from zero through
65,598 bytes. Unaligned rows have exact allocation ends, not readable SIMD
padding. Small and one-pair cases use independent scalar table arithmetic;
larger ranges compare unchanged per-group AVX-512 calls. Full encoder checks
compare every requested parity byte against an operation table whose range
callbacks invoke the unchanged per-group implementation. Shapes are
K255/R129/B8192, K257/R200/B8256, K1000/R199/B32768,
K1000/R200/B32768, K1000/R200/B65536, and K4096/R512/B4096.
Source hashes must stay unchanged; the reused work slab is poisoned between
the control and candidate calls. This is not a public API/decoder release
matrix or an independently linked Leopard1 oracle.

Initial failures remain in the scratch logs:

- Registration initially fell inside the dual-field test guard, so the
  GF16-only target was missing. It now registers outside that guard.
- The first encoder oracle nulled mandatory non-default Ops callbacks and
  crashed at a null function pointer; GDB identified the harness defect.
  Explicit unchanged callbacks replaced those nulls; candidate code did not
  change in response.
- The first combined sanitizer process hit its 256-MiB cap at the largest
  encoder shape (`run-u283305.scope`, `Result=oom-kill`). No final peak/counter
  record survived, so none is inferred. Slab allocation, one reused transform
  workspace, and separate range/encoder processes reduced memory without
  removing any case or raising limits. All final modes pass.

Local raw logs, build metadata, and disassembly:
`/tmp/leopard-gf16-cache-block.PaIO57`.
Retained source/log/artifact copies and checksum manifest:
`.research/leopard-79h/gf16-split-cache-block.ZmP8hQ`.

Final source SHA-256:

- `Leopard2BackendAVX2.cpp`:
  `911fc1dfb2113118b002c8c58776d520c397ea9a1b183900733989554e10e7aa`
- `test_gf16_split_cache_block.cpp`:
  `694906f0ca9c99bdf51b154dcfa6cf5a591b3c5a9342740961661b17173e6900`

Final test executable SHA-256 (Release ON, Release OFF, sanitizer ON):

```
359287f6bb3c63874cfd7a99a6c19f9097d1c83cf4d69801124abea85c557236
01a515f35924ce019e9583613fd327ed72ac8c30f0cef50b9afbba048ba570c0
0825a1dbdeee58a091b7875c250f375f1e1d1418841c718b57398b8d6918129f
```

Review provenance is Codex self-review plus deterministic tests under the
user's Claude opt-out, not independent-model `CONVERGED`. The unavailable
`deli-auto-research` skill was not replaced by a claimed automatic watchdog.
Existing v19 preregistration, runtime/build handoff, contamination contracts,
and exact-main promotion gates remain unchanged and incomplete. No historical
benchmark executable was run, no timing was analyzed, and no AUTO gate was
widened by this experiment.

## First diagnostic screen: invalid, attempt budget exhausted

Commit `26e3984` preregistered and pushed a same-source ON/OFF screen before
measurement: six fixed cells, three OFF/ON/ON/OFF rounds, 21 samples per
process, one attempt, local CPU4 with sibling68. Only the AVX-512 archive
member differs between the two production libraries. This is a candidate
filter, not an exact-Leopard1 comparison or production qualification.

All 12 untimed public-API preflights matched output/input hashes, scratch,
and selected routes. The first timed invocation, cell0/round0/slot0/OFF,
observed two non-idle sibling jiffies. The runner stopped **before any ON
timing**. The journal has `complete:false`, one invocation, and no analysis
key. No ratios or performance inference are admissible, and the one-attempt
budget is exhausted; this plan must not be retried or moved to another CPU.

The failed scope peaked at 131,006,464 bytes under 256 MiB; all six memory
event counters and swap were zero. The 73-file read-only evidence bundle
includes frozen binaries/libraries/source, all build/check logs, the failed
journal, and subsequent passive host snapshots:
`.research/leopard-79h/gf16-split-screen-work-failed.Y3csBw`.

- Outer `SHA256SUMS`: `50d2b2d7cf69270351ae8a0468a1e70ce2c62580a14512162f57e4cf8b0dd674`
- Failed journal: `688ea3a6d852f8ff0f4a68d2242a5bdd17d4172ca729902253204a08ced75e25`
- Scope log: `32e9784073646fca3db0107df397f6bebab8e557cb4c1d01b650c7a2bebb5954`

An OBS streaming thread was subsequently observed on sibling68 with affinity
0-127. This does not identify the source of the exact two ticks. Both SSH
servers also had active workloads. A passive snapshot window found no idle
sibling half on ripper, but several on foureyes; those historical snapshots
do not reserve a core or predict future quietness. No other process was
stopped or moved. The current candidate stays OFF.
