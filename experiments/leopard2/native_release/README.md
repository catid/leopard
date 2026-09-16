# Current native release comparison — preparation

Tracked by Beads `leopard-79h.57.12` (clock-free checks: `.12.1`).
This directory contains **no new performance results**. The timing plan below
must be committed and pushed before its one fixed attempt can be launched.
It does not replace or relax any historical experiment's contract.

The comparator is original Leopard1 commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, compiled using its native
Release policy, not the ISA-restricted atlas build. The initial Leopard2
candidate is production commit `e35b1f0c3a5ef9f241ad3b174b0fbfdb61bbccab`
with the default backends enabled. No R199/32-KiB experiment is retried and
its unqualified AUTO route remains disabled.

## Fixed encode qualification matrix

| ID | K | R | Bytes/shard | Leopard2 request | Selection reason |
|---|---:|---:|---:|---|---|
| copy | 1 | 1 | 4096 | AUTO | Known validation-overhead risk |
| small | 16 | 8 | 64 | AUTO | Small-message dispatch |
| gf8-high | 240 | 16 | 65536 | AUTO | High-rate GF8 |
| gf8-balanced | 128 | 128 | 65536 | AUTO | Balanced GF8 |
| gf16-inflation | 200 | 50 | 65536 | AUTO | Rounded parent crosses field boundary |
| gf16-gfni-region | 1000 | 200 | 65536 | AUTO | Previously qualified GFNI region |
| gf16-explicit-avx2 | 1000 | 200 | 65536 | AVX2 | Known explicit-backend deficit risk |
| gf16-large | 4096 | 512 | 4096 | AUTO | Larger shard count |

These are predefined coverage categories, not assertions that the current
release wins or loses. All eight remain in the eventual report, including
failed or inconclusive cases. Explicit AVX2 is a separate diagnostic, never
substituted for the AUTO product comparison.

`encode_probe.cpp --check N [new_parity_file]` runs two ordinary public full
encodes on identical deterministic input, verifies every input byte is unchanged,
and compares repeated output byte-for-byte. Compile it separately against each
unmodified archive, defining `LEO_NATIVE_RELEASE_CODEC_COMMIT` to that archive's
source revision; define `LEO_NATIVE_RELEASE_BASELINE` only for Leopard1. This
macro is a label, not proof: retain the archive, source and build identities too.
Compare the complete parity files between implementations, not just their hashes.

The default probe has no timing interface. An optional `clock_guard.cpp` link uses
`--wrap=clock_gettime,--wrap=gettimeofday,--wrap=clock` to reject direct clock
references from the executable/static codec. It does not claim to interpose
private clock calls inside shared libraries, or other clocks/instructions
such as `omp_get_wtime`, `time` or `rdtsc`. Context-backend IDs are not
operation-specific ISA attestations; the probe makes no measured-GFNI claim.

Workspace reporting accounts for the APIs' different layouts: Leopard1 parity
occupies the first R work rows; Leopard2 parity is separate from scratch.
`workspace_bytes + separate_output_bytes` is comparable caller-owned encode
storage **excluding input, codec/context allocations, global tables, pointer-array
metadata and probe verification copies**. It is not peak process memory or total
library memory.

Builds use the canonical lock, 512 MiB/no swap and one compiler job. Checks use
the same lock, 256 MiB/no swap and one process at a time. Freeze copies of every
archive/executable and record SHA-256 before and after native checks.

`check_encode.py` consumes a new read-only bundle containing `main`, `current`,
`main-guard`, `current-guard`, `guard-control`, `main.a`, `current.a`, and their
`SHA256SUMS`. It checks 32 probe processes, all full parity files, three positive
clock-guard controls, 20 invalid CLI calls and artifact hashes/metadata around
each child. Use a fresh output directory; retain failed output too. The caller
must hold the canonical lock and enforce the scope's resource limits. The
checker does not establish the archive build/source relationship by itself;
retain configure/build logs, compile commands and source identities separately.

Before any performance run, publish a separate qualified, fixed timing method
and artifact manifest: exact cost boundaries, workload order, sample/group
counts, same-binary controls, reserved CPU/sibling, resource and minimum timer
window gates, retained failure policy, and analysis. Setup, one-shot decoding
and reused decoding are not measured by this encode-only probe. Existing
`--attest-source` benchmark calls still run timers; do not call those clock-free
preflights. No performance claim or production promotion follows from these
correctness checks alone.

## 2026-09-16 qualification result

[Qualification summary](qualification.json): all eight shapes passed, including
32 plain/guarded processes and 41,030,144 parity bytes per variant. Three clock
guard positive controls and 20 invalid CLI calls passed. Native checks peaked
at 201,424,896 bytes under 256 MiB; fresh serial builds peaked at 395,300,864
bytes under 512 MiB. Build/check scope memory-event counters and swap values
were zero. The separate successful evidence-retention job reached the 256-MiB
cap and recorded 491 `memory.max` events, with no OOM or swap; it is not an
all-zero resource run. The retained archive contains the pre-retention report;
the linked public summary records this final retention footer as well.
Pure contract/mutation tests passed 6/6 in both normal and optimized Python;
local read-only static review found no remaining concrete defect.

The observed buffer totals do **not** show a general encode-memory saving:
for example, K240/R16/64-KiB uses 3,152,000 bytes of Leopard2 scratch plus output
versus 2,097,152 bytes of Leopard1 work including output. K1000/R200/64-KiB uses
29,915,712 versus 33,554,432 bytes. Neither comparison includes library-managed
tables or setup allocations. These are current API storage requirements, not
timings or total-memory claims. Failed/inconclusive historical timing attempts
remain unchanged; the release timing and final archive gates remain open.

## Qualified grouped timer and one-attempt plan

Beads `leopard-79h.57.12.2`; [exact plan](timing_plan_v1.json).
`encode_timing.cpp` adds separate `--measure`, `--synthetic` and `--exercise`
interfaces while retaining the default clock-free probe. Each implementation's
single driver object is linked with steady, synthetic and aborting clocks.
Synthetic-only public API wrappers independently count actual calls at every
clock endpoint. The steady executable has no public-call wrapper.

The eight workloads above use respectively 4,194,304; 1,048,576; 256; 256; 256;
16; 16; and 128 encodes per group. Each process performs two untimed parity
checks, four warmup calls and nine timed groups. Every group must last at least
20 ms. Allocations, input generation, scratch poisoning, parity comparisons and
hashing are outside the timed region; the public encode, result check, loop and
timer-adapter overhead are included, with no subtraction. This measures reused
encoding, not setup, one-shot latency or decoding.

The fixed cell-major schedule has three rounds per cell. Each round runs
native/current/current/native, native/native/native/native and
current/current/current/current, for 288 measured processes after 16 clock-free
preflights. Use the median normalized group time per process, the square root of
`slot0*slot3/(slot1*slot2)` per round, and the geometric mean across three rounds.
Every same-binary aggregate must lie within `[1/1.02, 1.02]`; any failure makes
the entire campaign inconclusive. An advantage or deficit additionally needs
all three comparison rounds to agree in direction and an aggregate beyond the
reciprocal-2% boundary. No confidence intervals or production promotion follow.

Controller CPU 0, worker CPU 52 and its SMT sibling 116 are fixed. A single
10-second passive gate and every measured child require zero sibling busy ticks;
each child must also record positive worker busy ticks. Resource limits are
256 MiB/no swap with zero memory-event counters, and each child has a 180-second
CPU/wall limit, 256-MiB address-space limit and 1-MiB output-file limit. The runner
owns both canonical and pair locks: **do not wrap it in another `flock`**.
An existing attempt directory is never reused. Retain failures; do not adjust
groups, choose another CPU, resume, pool runs or analyze a partial attempt.

`timing_plan_v1.json` pins the immutable archives/executables, qualified records,
collector sources, observed dynamic libraries, launcher tools, host and exact
compiled source/object identities. File hashes and metadata are checked around
every child. This is not a hermetic runtime or an operation-specific ISA audit.
The compiled header retained in the qualification archive differs from the
tracked header only in a comment clarifying timer overhead; its actual build
hash is recorded separately. No library sources changed after `e35b1f0`.

### Qualification evidence (2026-09-16)

- 94 native qualification invocations passed: 2 helper unit programs (26 cases
  each, including ASan+UBSan), 48 clock-free checks, 16 synthetic full-group
  runs, 4 abort-clock exercises, 2 positive abort guards, 10 clock faults and
  12 invalid CLI cases. All full synthetic parity files equal the previously
  qualified native files. No real benchmark clock was read.
- The corrected exact-inventory guard audited all 94 saved invocations and
  complete parity again in normal and optimized Python without codec reruns.
  Static review found and fixed replay observation reuse and missing-event
  acceptance; 17 contract tests plus 4 full mocked collector/replay tests pass
  both Python modes, including retained failure and no-retry checks.
- Default non-timing probe syntax checks pass for both implementations. Timer
  build peak: 67,526,656/512 MiB; native qualification: 175,947,776/256 MiB;
  final tests/audit: 42,819,584/256 MiB; default syntax: 47,173,632/512 MiB;
  launch-envelope checks: 132,366,336/256 MiB. All events and swap were zero.
  The final five-test collector/preregistration gate (including exact tracked
  plan-to-source hashes) also passes both Python modes: 34,353,152/256 MiB,
  zero events/swap; its log is retained beside the qualification archive.
- Read-only qualification archive SHA-256:
  `be0dc1bfbe6e4d4ff8b3703b193f12d8bb27d9ee1f125ea407941ddbd19721f2`.
  Retention peak: 194,641,920/256 MiB, zero events/swap. Its exact local path
  and the 12-input campaign manifest are in the plan. Raw artifacts remain
  local, excluded from source releases; method and pins are tracked in Git.

Local read-only Codex review found no remaining blocking defect. Claude was
explicitly opted out; this is not an independent-model `CONVERGED` claim.

After publication, run the fixed attempt with the actual preregistration commit:

```sh
systemd-run --user --scope --expand-environment=no \
  -p MemoryMax=256M -p MemorySwapMax=0 \
  bash experiments/leopard2/gf16_high_encode/tower_public_scope.sh \
  python3 -B experiments/leopard2/native_release/run_timing.py run \
  --plan experiments/leopard2/native_release/timing_plan_v1.json \
  --preregistration-commit COMMITTED_AND_PUSHED_40_HEX_COMMIT
```

`run_timing.py replay --plan PLAN --output RETAINED_ATTEMPT --bundle ARTIFACTS`
checks canonical raw observation labels, qualified preflights, full schedule,
resource/CPU observations and final analysis. Run replay normally and with
`python3 -O` before reporting measured results. The plan is machine-specific;
reproduction elsewhere requires its own declared host and fresh qualification,
not an unrecorded change to this attempt.
