# AVX2 forward/accumulating runtime control: focused checks pass

Bead `leopard-79h.38.5.4.18.3.1` remains in progress. This milestone is
**untimed and experiment-only**. Full native Leopard1 parity, public frontend
route/count/cost-boundary qualification and the performance gate remain open.
No production source changed; no improvement is claimed from instruction counts.

The qualified compile-time mode3 prototype (`b81f0a5`) is now selectable OFF/ON
within one exact binary. It changes forward and accumulating inverse pair
schedules only. The earlier inverse-only candidate is not included, and its
below-5% result must not be added to predict this candidate's performance.

## Implementation and code generation

The exact-source overlay consumes the earlier qualified static source
`be5a0775...`, preserves its arithmetic and scalar tails, and introduces
two-template-body wrappers. A default-OFF diagnostic Boolean is changed only
when all codec operations are quiescent and remains immutable during concurrent
execution. Selection occurs once per pair call, outside the 64-byte loop.
The existing inverse in-place implementation remains unchanged.

The experiment uses the current `45e2eff` both-field archives as its base, not
an older production core. Each candidate archive replaces the AVX2 member,
retains the other 23 members byte-for-byte, and adds one explicit experimental
control member (25 members total). No production archive or file is modified.

The actual untraced Release object (`0b9ab2b9...`) has these 64-byte loops:

| Loop | Runtime OFF instructions / stack references | Runtime ON instructions / stack references |
| --- | ---: | ---: |
| In-place forward | 40 / 4 | 36 / 0 |
| Prepared-range forward | 40 / 4 | 36 / 0 |
| Accumulating inverse | 45 / 4 | 39 / 0 |
| In-place inverse, unchanged | 40 / 4 | 40 / 4 |

All retain eight byte shuffles. Actual-object replay finds the two distinct
runtime loops for each affected function, one unchanged inverse loop, no
runtime-control reference or observer call inside those loops, and no excluded
EVEX/high-YMM/ZMM/ternary/GFNI instruction anywhere in the Release object.
The saved disassembly also binds ordinary forward ON to the nonzero-control
branch at `0x3c81c` and accumulating OFF to the zero-control branch at `0x3d47f`.

**Runtime OFF is not pristine production.** In particular, the earlier static
production prepared-range loop was 39 instructions / 3 stack references,
whereas runtime OFF is 40 / 4. Both the wrapper overhead and register allocation
can affect cost. A subsequent performance protocol must retain the original
Leopard2 product comparator as well as original native Leopard1; OFF/ON alone
must not be presented as improvement over the shipped product. No timing
protocol has yet been reviewed or registered for this runtime candidate.

## Focused semantic and control checks

All 105 invocations completed as expected: **87 positive checks and 18 CLI
refusals**, across plain Release, traced Release and full ASan/UBSan/LSan.
Both OFF and ON execute the same focused test matrix in every profile:

- 65,535 forward logs and 65,535 accumulating logs against the scalar oracle;
  918 boundary accumulations with input preservation, cancellation and guards.
- 384 forward-range cases plus 66,147 unchanged inverse pair cases and
  256 inverse-range cases.
- Eight public full/subset/bounds shapes, three decode round trips and
  four-thread both-field round trips against the unchanged GFNI reference.
- Default-OFF, setter rejection/preservation, both counter families/states,
  block accounting, thread-local isolation and both overflow paths.

The existing focused harness is not edited. Its private generated include
renames only the final entry point, preserving its nested legacy entry renames
and all checks. The retained original and generated include are pinned and
their exact relation is checked by replay.

Trace counters are absent from the untraced Release path. Focused trace totals
include backend initialization; they are not yet isolated public-workload
counts or timing attribution. Worker thread-local totals are not aggregated.
Traced OFF/ON and Release/sanitizer active-state totals match for every selector.
The independent replayer checks output inventories, workload/guard records,
selected/inactive counter states, counter types, source and archive members.

## Sanitizer finding and retained failures

The first full build's sanitizer control unit reported null member access in
the direct external-TLS overflow path. Release matrices had passed, but that
sanitizer failure stopped qualification before any sanitizer codec test.
The retained assembly showed an external TLS initialization/null-check sequence;
the underlying toolchain cause has not been isolated, so this is not a proven
compiler-bug attribution.

Counter storage now has internal linkage in its defining translation unit and
is accessed through an out-of-line reference function. Full sanitizer and
overflow tests remain enabled—there is no suppression. All corrected tests
pass. The **untraced Release object and archive are byte-identical before and
after this diagnostic-only fix**, so it does not change the candidate codegen.

Earlier preparation failures (copied read-only directory permissions) and a
focused-wrapper nested-main compile failure also remain retained. Those stopped
before codec execution. The private-directory fix never changes retained input
permissions. The 736-file sealed bundle includes all four superseded roots,
their source copies, binaries where built, terminal logs and failures.

## Replay, resources and remaining work

Eight pure overlay/parser/record adversarial tests pass normally and under
`-O`. Collector-free raw and sealed replays each pass in both Python modes;
all four replay results are byte-identical. Replays execute archive inspection
and disassembly tools, never codecs or benchmark clocks.

| Accepted scope | Peak bytes | Limit | Memory events / swap |
| --- | ---: | --- | --- |
| Corrected build | 255,401,984 | 512 MiB | all zero |
| Corrected native matrix | 187,101,184 | 256 MiB | all zero |
| Actual-object audit | 57,163,776 | 256 MiB | all zero |
| Raw replay/tests | 78,954,496 | 256 MiB | all zero |
| Retention and sealed replay | 192,954,368 | 256 MiB | all zero |

Every accepted scope exited zero. Work was local, serial under the canonical
lock, with per-child CPU limits and no swap. No Claude, subagents, remote workers,
host settings or unrelated process affinities were used or changed. Codex
self-review and deterministic/adversarial checks follow the user's Claude
opt-out; this is not independent-model `CONVERGED`.

Raw: `/tmp/leopard-avx2-adjacent-runtime-tls.mjeJwQ`.
Read-only: `.research/leopard-79h/avx2-adjacent-runtime-focused.yMDquC`.
Manifest: `86432e598901a3edd5c532b31b0d1b8b1591f28d725d1365870ffdf87fb8ad83`
(736 files, 210,828,145 bytes). Delivery/retention logs:
`/tmp/leopard-adjacent-runtime-delivery.v7Kupu`.
The [machine-readable result](results/avx2_adjacent_runtime_focused_20260909.json)
contains archive/object identities, actual loops, observations and open gates.

Next, complete the runtime candidate's clock-free public frontend: full original
native Leopard1 parity, original Leopard2 binding, all affected/unaffected
cells, isolated public route/call counts and clock-abort/cost-boundary checks.
Then—and only after a separately reviewed, committed and pushed preregistration—
run a direct performance experiment. Preserve 5% target, 2% control/neighbor,
zero-sibling and single-attempt gates. No retries or pooling of the consumed
inverse-only or GFNI attempts, and no default promotion from this milestone.
