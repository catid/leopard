# Native / pure-AVX2 Leopard1 / explicit-AVX2 Leopard2 screen

Bead: `leopard-79h.38.5.4.18.1`. Preregistration dated 2026-09-09.
No measurement has been made under this protocol at preregistration.

## Question

The [binary audit](avx2_isa_attribution.md) proves the original native Leopard1
GF16 AVX2-intrinsic path uses EVEX, high YMM registers and ternary logic that
Leopard2's explicit AVX2 contract excludes. Measure the native build-policy
effect separately from the remaining same-ISA implementation difference.
Keep the original native Leopard1 as the product comparison target.

This is not a single-instruction ablation: the existing pure-AVX2 Leopard1
profile changes `-march=native` to
`-march=x86-64 -mtune=generic -mavx2 -mno-avx512f`. All original codec source
bytes remain at `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`; all four objects
and both fields are included. Leopard2 uses the unchanged production archive
from `3a2f064`, with AVX2 explicitly requested throughout. AUTO/GFNI policy
and its completed qualification are not being remeasured.

## Fixed cells and protocol

| Cell | K | R | Bytes | Role |
| --- | ---: | ---: | ---: | --- |
| 0 | 1000 | 200 | 65536 | Remaining explicit-AVX2 deficit |
| 1 | 1000 | 200 | 32768 | Deficit boundary with explicit AVX2 |
| 2 | 1000 | 199 | 65536 | Deficit boundary with explicit AVX2 |
| 3 | 4096 | 512 | 4096 | Neighbor |

Each cell has three rounds of six four-process ABBA comparisons, in the
order fixed in `avx2_isa_screen_plan.json`: native/pure, native/current,
pure/current, same-native, same-pure, same-current. Same-profile slots execute
the identical path/inode; there is no relinked placebo. Total: **288 timed
children and 12 untimed server preflights**.

Each child allocates and initializes outside clocks, performs one checked
encode, four warmups and 21 measured encodes. Each sample covers one public
full-output encode call and the common checked-call wrapper. No plan setup,
source generation, hash, guard check, diagnostic route query, output copy or
file I/O is timed. Leopard1 parity is read from the first R work buffers;
Leopard2 uses separate output buffers. Both use the same original-data seed
`20260906`, shape and parity semantics, with process-local one-thread OpenMP.

A process contributes the median of its 21 samples. Each ABBA round uses
`sqrt((A0/B1)*(A3/B2))`; the aggregate is the geometric mean of three rounds.
For durations A/B, the resulting throughput ratio is B/A. Thus the three
reported ratios are `pure_over_native`, `current_over_native`, and
`current_over_pure`. Do not combine separately measured ratios to claim an
exact causal fraction of the previous deficit.

All **12 same-profile aggregate controls** must remain in `[1/1.02, 1.02]`
before any directional interpretation. The same-current neighbor is included
in that gate. Differences between native/pure builds at the neighbor are
attribution measurements, not regressions of a changed Leopard2 candidate.
There is no candidate in this experiment and no promotion decision, CI,
cross-host claim or v19 qualification. Any subsequent optimization retains
the separate 5% target / 2% controls-and-neighbors promotion requirements.

Run only on local `work` (9980X, family26/model8), CPU26 with sibling90;
controller CPU0. Hold the canonical build/test/timing lock and CPU-pair lease.
Require the recorded Slipgate/OBS shutdown state before and after, ten seconds
of passive zero-sibling activity, and zero sibling non-idle ticks around every
timed invocation. All three executables, archives, source, expected identities
and protocol are frozen and rehashed after every child.

Exactly **one attempt**, at `frozen/../attempt1`. No overwrite, retry,
CPU/host swap, dropped round or pooling with previous attempts. Failure stops
the attempt and retains its records; incomplete data gets no analysis. A
control failure makes the complete screen inconclusive. No unrelated affinity
changes, worker SSH hosts, Claude, subagents or host-wide cache drops.

## Pre-timing validation

The new drivers were compiled against three pinned Release archives and the
existing fully instrumented, both-field Leopard2 ASan/UBSan archive. Only
driver objects were rebuilt. Strict GCC warnings passed for the drivers.
All four variants passed ordinary checks plus link-wrapped `--check` and
26-encode `--exercise` paths that abort on a driver benchmark-clock read.
Each clock guard was deliberately triggered and returned86 before reading
the clock. These guards are not linked into timing executables.

**48 positive native records**, **139,198,464 full parity comparison bytes**
over16 file comparisons, four clock-guard rejections and32 malformed CLI
rejections passed. This includes all four cells for native/pure/current and
sanitized current, source-preservation and outer allocation guards. The
earlier pure-L1 attribution also checked GF8; this screen is GF16 only and
does not claim newly instrumented pure-L1 sanitizer or decode coverage.

Builds used a 512MiB/no-swap scope, serial compilation and compiler GC10/4096,
peaking at145,870,848 bytes. Native preparation used256MiB/no-swap and a
30-second child CPU limit, peaking at160,559,104 bytes. All six memory-event
counters and swap were zero. No benchmark clocks were read in preparation.
Protocol and independent-replay unit tests cover analytical ratios, controls,
identity/type drift, incomplete/reordered data and immutable-input inventory;
they run in normal and optimized Python.

The collector requires an ancestor of the fetched topic branch and verifies
that the frozen protocol/driver/collector sources match that commit. The
commit must actually be pushed before launching the attempt. The independent
replayer imports no collector, executes no codec, reconstructs all72 round
ratios and24 aggregates from raw stdout, checks conditions/resources/pins,
and re-compares every preparation parity byte to the retained native oracle.

## Completed result

Preregistration `a5f55792fcc955a5aa4b997848b3ba58252b611b` was pushed before
the sole attempt. All288 timed children and12 preflights completed. Every
timed sibling delta was zero; passive ticks stayed570324 over10,000,119,040ns.
All12 aggregate same-profile controls passed, ranging from0.990702 to1.008323.

Throughput ratios (larger means the numerator is faster):

| Cell | Pure L1 / native L1 | L2 AVX2 / native L1 | L2 AVX2 / pure L1 |
| --- | ---: | ---: | ---: |
| K1000/R200/64KiB | 0.838397 | 0.988986 | 1.184011 |
| K1000/R200/32KiB | 0.786287 | 0.916724 | 1.166437 |
| K1000/R199/64KiB | 0.827890 | 0.972627 | 1.189710 |
| K4096/R512/4KiB | 0.792902 | 0.997242 | 1.228300 |

Leopard2's explicit AVX2 implementation beats the AVX2-restricted Leopard1
comparator by16.6–22.8%, in the same direction in every round of every cell.
The generic hypothesis that Leopard2 is slower under the same AVX2 ceiling
is not supported by these measurements. Native Leopard1 is substantially
faster than its own restricted build, consistent with the static codegen
finding; do not infer a unique instruction-level cause or exact causal share.

The original native target remains relevant. Native L1 has higher aggregate
throughput than explicit L2 at R200/32KiB (about9.1%) and R199/64KiB
(about2.8%, with variable rounds). The64KiB R200 target has mixed rounds
`0.988224, 1.001007, 0.977862` for L2/native and is treated as near parity,
not a stable win or reproduced3.12% gap. The4KiB neighbor is also near parity.
Do not pool these results with the earlier driver or call a difference
between experiments a production regression: this run is an unchanged-codec
three-build attribution contrast with a new, common driver.

Two individual native-control rounds were outside2%: cell0/round2 at1.024378
and cell2/round1 at0.974963. The preregistered gate was on all12 three-round
aggregate controls, which passed. Retain those rounds; there are no CIs and
no claim that every individual control round was within2%.

The timing scope took121.22s and peaked at129,486,848 bytes under256MiB,
all six memory-event counters and swap zero. Independent replay passed in
normal and optimized Python: all18 frozen inputs,72 round ratios,24
aggregates,302 stdout/stderr pairs, shutdown conditions, resources,48
preparation checks and139,198,464 full parity bytes. Replay peaks were
150,106,112 and16,076,800 bytes under256MiB, all counters/swap zero.

The complete numerical result is
[`results/avx2_isa_screen_20260909.json`](results/avx2_isa_screen_20260909.json).
Raw workspace: `/tmp/leopard-avx2-isa-screen.xsFlEg`; the retained read-only
copy and its outer manifest are recorded in Beads. Key hashes:

- Attempt journal: `9358f12f16525ea1c04c4a336412fc828d6ea940d420a7b3df39cc8d9e1213f3`.
- Timing resource log: `01c79f55f80000202cca904b76f4011f7c1fe1a8c3c966050950d73fbebcf5ef`.
- Frozen pins: `14b646005c3d0ddee8e3418a14db5b61bbd62a288ccba7ba166931f68aabaf10`.

This completes the attribution experiment, not the broader performance goal.
There was no production code change or promotion. The completed AUTO GFNI
improvements at the two boundary cells remain separate. A next AVX2 candidate
should address a concrete kernel cost while preserving the ISA contract;
the inspected two-way inverse loop reloads four spilled table vectors per
64-byte iteration. A reduced-live-range product schedule is a new hypothesis
to test, not a measured improvement. Keep the existing tables, split
transform, byte tiling and public API semantics; do not repeat rejected broad
fusion/cache-block/copy-removal experiments or claim API overhead attribution.
