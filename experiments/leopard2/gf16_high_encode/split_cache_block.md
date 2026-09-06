# Rejected GF16 split-butterfly cache-block candidate

Bead: `leopard-79h.38.5.4.9`. Base commit: `2aca7f5`.
Date: 2026-09-06. Status: **rejected by the preregistered server screen**.
The first workstation attempt was invalid; the separate foureyes attempt
completed with no target reaching 5%. No historical exact-Leopard1 gap is
closed. Experimental production code has been removed; parity tests remain.

## Mechanism and cost hypothesis

The now-removed `LEO2_EXPERIMENT_GF16_SPLIT_CACHE_BLOCK=ON` changed only the AVX-512VL nibble
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

## Server diagnostic screen: rejected

The distinct foureyes plan and its runner were committed and pushed at
`899eb04` before measurement. CPU22/sibling86 was selected using passive
activity evidence, not benchmark results; both topology files confirm the
physical pair. The server is a Threadripper PRO 9985WX, family26/model8,
kernel6.8.0-138-generic. Other jobs remained active and untouched. This does
not establish performance on the historical 9950X3D host.

The new plan retained all six cells, three ABBA rounds, 21 samples per
process, one attempt, and the prior decision rules. It added a fixed
10-second passive sibling check. This observed 190,909 non-idle jiffies
before and after a 10.000051393-second window. All 12 untimed preflights
matched; all 72 timed invocations observed zero sibling activity and matched
the frozen source/archive/executable hashes and workload identities.

Ratios are OFF time / ON time; values above one favor the experiment.
These are diagnostic point estimates, **not confidence intervals**.

| Cell | Route and workload | Role | OFF / ON |
| --- | --- | --- | ---: |
| 0 | AVX-512, K1000/R200, 32 KiB | Target | 1.001299 |
| 1 | AVX-512, K1000/R200, 64 KiB | Target | 0.995270 |
| 2 | AVX-512, K1000/R199, 64 KiB | Neighbor | 0.993820 |
| 3 | AVX-512, K4096/R512, 4 KiB | Unchanged control | 1.001698 |
| 4 | AVX2, K1000/R200, 64 KiB | Unchanged control | 1.013187 |
| 5 | AUTO GFNI, K1000/R200, 64 KiB | Unchanged control | 0.996994 |

All control aggregates lie inside `[1/1.02,1.02]`; the neighbor passes its
floor. Neither target reaches 1.05, so the committed decision is
`reject_for_this_screen`. The 64-KiB target was below one in all three
rounds, but this screen does not establish a statistically significant
slowdown. The result rejects this implementation for further promotion,
not all possible cache-blocking techniques or other processors. There is no
exact-Leopard1 comparison, speedup claim, or widening of AUTO routing here.

The scope exited zero after 42.30 seconds, peaked at 131,469,312 bytes under
256 MiB, and recorded all six memory event counters and swap as zero.
A separately written stdlib replay, without importing the runner, checked
every raw stdout/stderr and all pinned inputs and recomputed all 18 round
ratios and six aggregates. It passed in normal and optimized Python modes.
This is deterministic replay plus Codex self-review, not independent-model
review.

- Preregistered plan SHA: `3a0c14b0a20ca68419646b5517bd090a262ff95d7635be2cf7f0c47319448181`
- Frozen input manifest SHA: `a0b511de5148216b7ec09b32c0630899bb580b05a48b8787db6f54db7b1118f3`
- Raw journal SHA: `e46f9b0f326e5b6fab8827a7b29ac281c5bdf82c6695141f201fbe7e301bb45e`
- Scope log SHA: `09f032a981dc67a143e20f226a12ced199c522e5c3b948548dab8d8f9410540a`

The failed workstation samples contribute nothing to this analysis. Both
attempts and the original candidate source remain recoverable from their
archives and Git history. The production helper and experiment option are
removed. The 2,048 range cases and six full-encoder parity cases are retained
as `tests/leopard2/test_gf16_split_range.cpp`, with CTests
`leopard2_gf16_split_ranges` and `leopard2_gf16_split_encoder`.
`Leopard2BackendAVX2.cpp` now exactly matches pre-experiment commit `2aca7f5`
(SHA `e154659949041d5739b4c85767ac8783d8c65bb5bb95be5cd437cdb9023203ff`).

Fresh cleanup validation passed all three focused CTests in both Release and
ASan+UBSan+LSan: all retained range/encoder cases plus the hookless production
AUTO-GFNI route, which ran rather than skipped. Release test peak was
144,924,672 bytes; sanitizer test peak was 175,357,952 bytes, both under
256 MiB. The largest fresh build peaked at 326,172,672 bytes under 512 MiB
with one compile job. All six event counters and swap were zero. The fresh
Release `libleopard.a` is byte-for-byte identical to the screen's default-OFF
archive (`24141e430ba35048548b122d26b1e77ed2914d1119e0ea74c1a6b7e30a8f2acf`).

The complete 208-file read-only screen/cleanup evidence bundle is retained
locally and on ripper at
`.research/leopard-79h/gf16-split-screen-foureyes.bpRucZ`, outer
`SHA256SUMS` digest
`b65f0563e7afc965f28770a0df9c15ce3319a2647c2d0cb95ee9520c57d8cee8`.
It contains all raw invocations, frozen inputs and plan, host snapshot,
separate replay source/logs, and fresh cleanup sources/build metadata,
libraries, test executables, and validation logs. Original v19 artifacts and
incomplete qualification gates are untouched. Follow-up Bead
`leopard-79h.38.5.4.10` covers a current-route diagnostic against standalone
Leopard1; the overall residual-gap task remains open.
