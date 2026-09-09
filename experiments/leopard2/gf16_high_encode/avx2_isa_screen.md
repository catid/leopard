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
