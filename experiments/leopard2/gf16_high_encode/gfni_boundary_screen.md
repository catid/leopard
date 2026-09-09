# GFNI at the two measured AUTO deficit boundaries

Bead: `leopard-79h.38.5.4.17`. Date: 2026-09-09. Status: preregistration.

The completed post-Slipgate diagnostic found current AUTO throughput at
0.971964x Leopard1 for K=1000/R=200/32 KiB and 0.956348x for
K=1000/R=199/64 KiB. Both use AVX2. In contrast, the qualified exact
R=200/64-KiB AUTO/GFNI target is already 1.415585x Leopard1. The original
`8da10dd`/`39a9282` route campaign qualified that exact cell and inactive
neighbors; it did not establish that explicit GFNI loses at these boundaries.

This new experiment tests existing explicit GFNI, without changing the codec,
AUTO predicate, kernel, API, flags, tables, or build configuration. It is not a
retry of either exhausted earlier current-route plan and uses none of their
samples. Prior source-copy/cache-block and sub-5-percent kernel experiments
are not repeated. All work is local; neither SSH worker host is accessed.

## Correctness before clocks

The 158-entry original correctness/build bundle was rehashed in full. All
current codec source inputs still match its pinned source, commit
`36dc0c8f66604b8d974468e687c51e6183ecb61d`. Its unchanged, both-field,
hookless/default Release archive is
`259c2c9aa3f51b1f941eac84270ad37c39c2a5b3cc06b1676a414f8088ad88e6`.
The separately linked Leopard1 executable and native archive remain pinned to
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`; its build is not ISA-matched to
an explicitly restricted Leopard2 backend.

Only new small diagnostic drivers were compiled. New untimed runs compared
every output byte for AUTO, explicit AVX2 and explicit GFNI against Leopard1
at both target shapes: 58,785,792 comparison bytes. All six driver records
also match runs linked against the unchanged both-field ASan/UBSan archive
(`1cac86a52866496279b4b8a08158ebc71438e14b062332f3516ed6beece88b58`),
with leak detection enabled.

Eight additional guarded shapes passed in Release and ASan/UBSan/LSan:
the two targets aligned, the same targets misaligned by one byte, two
even-byte tails misaligned by two bytes, and small GF8/GF16 scalar-tail
cases. Each checks six full/subset masks against every byte of explicit AVX2
parity, preserves unrequested outputs and inputs, checks per-shard/scratch
canaries (and ASan red zones), and rejects short scratch without mutation.
GF16 odd physical lengths are also rejected. ASan's sub-eight-byte shadow
granularity is not a claim of byte-exact poisoned left red zones for unaligned
addresses; canaries still cover those bytes for writes.

Build scope: peak 130,965,504 bytes / 512 MiB, all six events zero, swap zero.
Complete native checks: 15.07 seconds, peak 256,364,544 bytes / 256 MiB,
all six events zero, swap zero. Seven pure collector tests pass normally and
under `python -O`, including invalid protocol/identity/sample/order/isolation
and control-failure tests. Review is Codex self-review and deterministic
checks under the user's Claude opt-out, not independent-model `CONVERGED`.

## Frozen measurement protocol

`gfni_boundary_screen_plan.json` is committed and pushed before clocks.
It pins all executable, archive and expected-record digests. A fresh
lane-owned read-only input directory binds these plus collector, dependencies,
drivers, plan and shutdown checker to the preregistration commit. Inputs are
rehashed before the attempt and after every child. Builds, correctness and
timing are serialized through `/tmp/leopard-gf8-authoritative.lock`; timing
also holds the CPU-pair lease. No other thread's affinity is changed.

- Host: local `work`, Threadripper 9980X, model 08h. CPU 26, sibling 90,
  controller 0. Slipgate's three services remain inactive/disabled and its two
  containers remain stopped with restart `no`, checked before and after.
- Exactly one attempt, two cells in listed order, three rounds per cell.
  Each round contains AUTO/GFNI/GFNI/AUTO, L1/GFNI/GFNI/L1, then the same GFNI
  executable/argument/inode as an A/B/B/A control: 72 timed children total.
- Eight untimed workload checks precede a ten-second zero-sibling passive
  gate. Every timed child's sibling non-idle delta must be zero. Any failure
  stops the attempt; no partial ratios, retry, CPU swap, trimming or pooling.
- Each child uses aligned, unique source shards (seed 20260906), one initial
  encode, four warmups, and 21 samples of one public full-output encode.
  Setup, allocation, hashing and parity dumps are outside clocks. Leopard1
  exposes its first R work buffers without an artificial final copy.
- Ratios derive from median process samples, geometric ABBA round ratios and
  geometric means of three rounds. Both aggregate same-GFNI controls must be
  within `[1/1.02, 1.02]` before any directional inference. A target must gain
  at least 5% over current AUTO with all three rounds above 1 to proceed.
- Native scope: 256 MiB, swap zero; each codec child has a 30-second CPU cap.
  All six memory-event counters must be zero. A positive result is only a
  reason to qualify a future bounded AUTO candidate; no confidence interval,
  unchanged-neighbor qualification, production promotion, broad CPU claim or
  v19 closure is asserted. Explicit AVX2 must remain explicit AVX2.

Raw preparation is at `/tmp/leopard-gfni-boundary.iIwWWx`; the original
reference is `.research/leopard-79h/gf16-current-route-failed.STc10h`, outer
manifest SHA-256
`e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`.
