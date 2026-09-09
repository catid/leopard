# Same-binary AVX2 pair-scheduling screen

Tracking: `leopard-79h.38.5.4.18.2`. This is a new, single-attempt diagnostic
screen of the [validated scheduling prototype](avx2_pair_schedule.md).
Production integration and broad qualification are separate work.

## Actual runtime-controlled candidate

The private source overlay [avx2_pair_runtime.patch](avx2_pair_runtime.patch)
instantiates original and scheduled inverse-pair bodies. A quiescent-only
experimental boolean selects the body outside the 64-byte loop. The existing
tables, transform, tiling, source-copy policy, forward/scalar behavior, other
backends and public wire format are retained. OFF and ON use **one identical
executable path/inode**; only the mode argument differs. OFF includes the new
selector and code layout and is not a byte-identical pristine-production build.

The actual uninstrumented runtime AVX2 object is
`f94d2b6bf1199f9c325fdb719b1ff32394d16d8b2e27fb077a8fb82c6d45b581`.
Disassembly confirms one selector read before either vector loop: OFF has
40 instructions/four stack table loads; ON has38 instructions/no stack
references. Both have eight shuffles. ON has two additional register moves
relative to the earlier compile-time prototype; only the new binary will be
timed. A full-object scan finds no EVEX, high YMM, ZMM, ternary logic or GFNI.
Static instruction counts are not time shares or speed estimates.

Diagnostic Release and ASan/UBSan variants have thread-local branch and
64-byte-block counters. **The timing member has none of this instrumentation.**
The timing driver's trace getter is checked outside clocks and returns zeros.
Diagnostic binaries reject `--measure` before initializing a codec.

## Preparation completed before timing

Both modes passed the full focused Release and ASan/UBSan/leak matrix from
the prototype: all65,535 ordinary multiplier logs plus612 pair boundary cases,
256 scalar-oracle split ranges, eight public subset/bounds shapes, three decode
round trips and four-thread both-field sharing. Neither these tests nor the
ordinary parity checks read a benchmark clock.

The new common driver supports eight cells, single-encode `--check`, and
26-encode clock-free `--exercise` matching the one checked encode, four warmups
and21 samples of the future timed path. All four builds (native L1, normal L2,
traced L2, sanitized L2) passed checks and exercises; normal builds also passed
unwrapped checks. Thread-local route counts multiplied by exactly26 in each
exercise, with no cross-mode or unchanged-neighbor activity.

| Cell | K/R/bytes | L2 request | Prepared calls per encode | 64-byte blocks |
| --- | --- | --- | ---: | ---: |
| 0 | 1000/200/32768 | AVX2 | 3672 | 1880064 |
| 1 | 1000/199/65536 | AVX2 | 7344 | 3760128 |
| 2 | 1000/200/65536 | AVX2 | 7344 | 3760128 |
| 3 | 4096/512/4096 | AVX2 | 12544 | 802816 |
| 4 | 1000/199/32768 | AUTO, remains AVX2 | 3672 | 1880064 |
| 5 | 1000/200/32768 | GFNI | 0 | 0 |
| 6 | 1000/200/32768 | AUTO, qualified GFNI | 0 | 0 |
| 7 | 17/7/64, GF8 | AVX2 | 0 | 0 |

Counts are equal in OFF and ON; the selected counter changes. These are
dynamic operation counts, not measured time attribution.

Preparation produced **184 positive native records**, three deliberate
clock-abort exits86, four diagnostic timing refusals and20 malformed-CLI
rejections. There are **52 full parity comparisons totaling361,368,192 bytes**:
all new profiles/modes match the separately linked original native Leopard1,
and four original-oracle files bind the new native driver to earlier retained
reference bytes. The native archive remains exact source6e5725eb with its
original native compiler target; this is not an ISA-matched comparison.

The candidate archives preserve all23 other production members byte for byte,
replace only the AVX2 member and add one experimental control member. Both
fields and full sanitizer instrumentation remain enabled. Read-only preparation
replay re-compared all raw records, mode/exercise counts, archive members and
parity bytes. Seven analytical/adversarial protocol tests passed in normal
and optimized Python, including all19 controls, every target/neighbor gate,
ordering, sample/type mutations and comparison with the independent replayer.

Resource peaks, all six memory-event counters and swap zero:

- Serial codec/control builds:256,970,752 bytes under512MiB.
- Strict driver builds:137,785,344 bytes under512MiB.
- Native preparation:214,880,256 bytes under256MiB; child CPU limit30s.
- Read-only preparation replay:36,085,760 bytes under256MiB.

Codex self-review and deterministic checks supply the review evidence.
The user's Claude opt-out remains in effect; no independent-model
`CONVERGED` result is claimed. No production source is modified.

## Preregistered attempt

The immutable plan is [avx2_pair_screen_plan.json](avx2_pair_screen_plan.json).
Commit and publish it with all protocol and candidate source files **before
any benchmark clocks**, then fetch and verify its topic-branch ancestry.

Local host `work`, AMD family26/model8, CPU26/sibling90, controller CPU0.
Hold the canonical build/test/timing lock and CPU-pair lease. Require the
recorded Slipgate/OBS services disabled/inactive and containers stopped with
restart=no before/after; require ten seconds of passive zero sibling activity
and zero sibling non-idle ticks around every timed invocation. No unrelated
affinity changes or remote workers.

Three rounds per cell. Each round uses OFF/ON/ON/OFF, same-OFF and same-ON.
Cells0–2 additionally use native/ON/ON/native and same-native. Each process
has21 samples; use its median and each ABBA ratio
`sqrt((A0/B1)*(A3/B2))`, then the three-round geometric mean. Ratios represent
throughput of B over A. Total: **360 timed children,19 untimed preflights,
90 round ratios and30 aggregate ratios**. Initialization, filling, hashing,
route checks, guards, output copying and file I/O are outside clocks; the
common checked public encode call is inside.

All19 same-path aggregate controls and the three unchanged neighbors5–7
must remain within `[1/1.02,1.02]`. A failure makes the screen inconclusive.
Affected neighbors3–4 must not regress by more than2%; improvements there
are allowed because the candidate actually changes their AVX2 paths.
**All three targets0–2 must improve by at least5%** for a positive candidate
screen. Keep every round, including individual control outliers; the gate is
on aggregates. There are no confidence intervals, production promotion or
v19 conclusions from this diagnostic alone.

Exactly one attempt at `frozen/../attempt1`. No overwrite, retry, host/CPU swap,
dropped round, old-plan reuse or cross-experiment pooling. An incomplete
attempt has no performance conclusion. The collector and independent replayer
rehash the frozen binaries, archives, expectations and source files. Frozen
inventory:20 files; pins SHA
`3ed2c034c4348b731d6de3b2a1d3c29d22a822332a3eccb18372e37a50884c03`;
plan SHA `9807c751da5c98a67a6ba3a95908f851c6d926eb07ab02559625b495cb8f0153`.

Raw workspace: `/tmp/leopard-avx2-pair-runtime.W1LY8Q`. Results will be added
after the sole attempt and independent replay. The broader performance goal
remains open regardless of this experiment's result.
