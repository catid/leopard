# AVX2 adjacent scheduling: public frontend qualification

Bead `leopard-79h.38.5.4.18.3.1`: **clock-free qualification passed**.
This is an experiment-only, clock-free qualification of the forward and
accumulating runtime candidate from [the focused milestone](avx2_adjacent_runtime.md).
There is no new performance result or production/default change.

## Comparator and source binding

Five separately linked profiles retain the exact already-qualified archives:
original native Leopard1, original current Leopard2, runtime Release,
traced runtime Release, and full ASan/UBSan/LSan runtime Leopard2.
No backend is rebuilt for this frontend milestone. Original Leopard2 and
both runtime Release profiles share **one compiled driver object**, not
separately compiled copies of similar source. Small link-specific shims
supply archive identity and control/count access. Native and sanitizer driver
objects are compiled separately; every profile reuses its driver object across
plain, synthetic-clock, and abort-clock links.

The original eight `avx2_pair_screen.cpp` cells and seed `20260906` are retained.
Cell8 adds a one-item batch of cell0. In particular, GF8 cell7 explicitly requests
AVX2, not AUTO. Original Leopard2 accepts `PPPP`, native accepts `NNNN`, and
runtime profiles accept `0110`, `1001`, `0000`, or `1111`.

**Runtime OFF is not pristine production**: its prepared-range loop is 40
instructions/4 stack references versus production's 39/3. Any subsequent
performance comparison must retain original Leopard2 and native Leopard1;
OFF/ON alone cannot establish a shipped-product improvement. The earlier
inverse-only scheduling change is not included in this candidate.

## Isolated public-workload counts

The first of four preflight public calls may initialize caches. Counters reset
before each subsequent isolated probe, then reset again before warmup/exercise.
Each of the three isolated probes is checked against the independent structural
traversal and edge-count model, not against counts copied from another runtime.
The one-item batch must execute the same pair work as its ordinary counterpart.

| Cell | K / R / bytes | Requested backend and actual path | Forward / accumulating pairs per public call |
| --- | --- | --- | ---: |
| 0 | 1000 / 200 / 32768 | AVX2 | 665 / 384 |
| 1 | 1000 / 199 / 65536 | AVX2 | 1330 / 768 |
| 2 | 1000 / 200 / 65536 | AVX2 | 1330 / 768 |
| 3 | 4096 / 512 / 4096 | AVX2 | 1793 / 1792 |
| 4 | 1000 / 199 / 32768 | AUTO → AVX2 | 665 / 384 |
| 5 | 1000 / 200 / 32768 | explicit GFNI | 0 / 0 |
| 6 | 1000 / 200 / 32768 | AUTO → GFNI | 0 / 0 |
| 7 | 17 / 7 / 64, GF8 | explicit AVX2 | 0 / 0 |
| 8 | cell0, one-item batch | AVX2 | 665 / 384 |

Counts of 64-byte blocks are also verified, for both families and OFF/ON states:
512 blocks per pair for the 32 KiB tiles, or 64 for cell3. The 64 KiB public
calls already include both tiled passes in their pair counts. Zero counts in
GFNI/GF8 mean the modified GF16 AVX2 functions are bypassed, not zero codec work.
These counts are not CPU-time shares or measured speedups. The R19932 AUTO
GFNI extension remains OFF; its consumed inconclusive attempts are not retried.

## Public-call and clock boundary

The qualified `PairedGroupTiming.h` loop is reused unchanged: start endpoint,
N complete public calls, end endpoint. Public result checks, per-call counters,
and repetition-loop overhead are inside the span. State selection and route
inspection, sample normalization/storage, per-slot bookkeeping, and full output
comparison/guards are outside. Between-group full parity checks touch output
memory and therefore affect cache state; this must be part of any future method.

Four preflight calls precede four warmup passes and 21 exercise passes, each
with four schedule slots. Group1 means 104 public calls; group256, allowed only
for GF8 cell7, means 25,604. There are 84 synthetic spans and 168 endpoints.
Grouped averages are not single-call latency samples. Qualifying these two
group sizes does not choose a timing protocol.

The external public-entry witness verifies identical codec, input/output arrays,
shard pointers and scratch pointers, plus exact API/state order and counts.
It observes the adjacent-scheduling control, not the unrelated GFNI extension.
Synthetic elapsed values are deterministically `257 + 17*i` nanoseconds and
must normalize exactly. Abort-clock links terminate before a benchmark clock;
the same links also complete ordinary clock-free exercises. Plain executables
reject `--clock-guard` rather than accidentally reading their real clock.
No `--measure` invocation is made, including as a negative CLI test.

## Validation and remaining gate

All 449 native records matched expectations: 307 positive checks (305 public
executions and two 37-case grouped-timer units), 42 pre-clock aborts, 20 malformed
synthetic-clock rejections, and 80 CLI refusals. There are 11,760 qualified
synthetic spans. Two canary checks and two expected ASan poisoned-boundary
read rejections additionally exercise this frontend's actual Buffer class.

Full-byte comparison checks 296 parity files / **1,744,826,816 bytes** against
new original-native executions. Those native outputs also reproduce nine
prior-native parity files / 60,981,696 bytes (cell8 maps to prior cell0).
Nine pure adversarial tests pass normally and with Python `-O`, both in the
worktree and sealed copy. Raw and sealed collector-free replays pass in both
modes; all four output hashes are
`4e328eb98f873e2b7abf8255676dc44b34dd27af3e7d16436933032ac816c307`.

Production source and existing evidence remain untouched. The previous focused
matrix is carried forward by exact archive and record identities, not rerun.
The shared Release driver object is
`94128723910cea77359a781767219dcc847fd15b14f40ed91972e172a4fec5da`.

| Scope | Peak bytes | Cap | Resource outcome |
| --- | ---: | --- | --- |
| Frontend build | 203,186,176 | 512 MiB | exit0, all events0, swap0 |
| Native matrix | 190,586,880 | 256 MiB | exit0, all events0, swap0 |
| Guard build | 86,859,776 | 512 MiB | exit0, all events0, swap0 |
| Guard checks | 21,233,664 | 256 MiB | exit0, all events0, swap0 |
| Raw replay/tests | 95,191,040 | 256 MiB | exit0, all events0, swap0 |
| Retention/sealed replay | 268,435,456 | 256 MiB | exit0, **6,738 max events**, no OOM or swap |

Evidence copying/replay reached the memory cap; it is explicitly **not an
all-zero resource run**. The cap was not raised, and all four replay outputs
and the post-replay full-file manifest agree. This does not alter the separate
all-zero native/build results or establish benchmark stability.

Raw: `/tmp/leopard-adjacent-public.84pxnc`.
Read-only: `.research/leopard-79h/avx2-adjacent-public.NFkTeG` (1,401 files,
1,920,371,965 bytes). Manifest:
`f3e58727918ded2e06b62c250ff1b19bde5d7860c107dbf9a8166a67fa0f4dfa`.
Delivery logs: `/tmp/leopard-adjacent-public-delivery.KCD7Pf`.
See the [machine-readable result](results/avx2_adjacent_public_20260909.json).

Work is local-only, serial under `/tmp/leopard-gf8-authoritative.lock`, with
512 MiB builds, 256 MiB checks/replays, no swap, and per-child CPU limits.
No Claude, subagents, SSH workers, host-setting changes, or unrelated affinity
changes are used. Review is Codex self-review plus deterministic/adversarial
checks under the user's Claude opt-out, not independent-model `CONVERGED`.

Next is `leopard-79h.38.5.4.18.3.2`: a separately reviewed, committed **and
pushed** direct timing preregistration before any benchmark clock. Preserve
5% target, 2% control/neighbor, zero-sibling and single-attempt gates, and both
original product comparators. Performance and production integration remain
separate, uncompleted gates; the broader Leopard2 performance goal stays open.
