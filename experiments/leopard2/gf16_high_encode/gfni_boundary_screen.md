# GFNI at the two measured AUTO deficit boundaries

Bead: `leopard-79h.38.5.4.17`. Date: 2026-09-09. Status: positive diagnostic;
bounded AUTO implementation/qualification remains open.

Implementation follow-up: the [default-off AUTO candidate](auto_gfni_boundary_candidate.md)
now passes 60 focused Release/sanitizer checks. Its new AUTO performance gate
is still pending; the explicit-backend results below are unchanged.

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

## Results

Preregistration `d1d3cbfbe6e2fa94460c32bc374f1a0950aaf2e3` was pushed before
the sole launch. All eight untimed checks and 72 timed children passed.
Sibling 90 remained at 570109 non-idle jiffies through the 10.000091881-second
passive window, and every timed invocation had delta zero. Both shutdown
snapshots matched. Native scope exited zero after 39.32 seconds, with peak
127,823,872 bytes / 256 MiB, all six memory events zero, and swap zero.

Ratios are GFNI throughput divided by the named comparator's throughput;
greater than one favors GFNI. These are diagnostic ratios, not confidence
intervals or evidence for another CPU.

| K/R/shard bytes | GFNI / current AUTO | GFNI / Leopard1 | Same-GFNI control |
| --- | ---: | ---: | ---: |
| 1000/200/32768 | 1.544798 | 1.462174 | 0.997978 |
| 1000/199/65536 | 1.471256 | 1.441381 | 1.010200 |

Both targets clear the 5-percent screen gate, with all three target rounds
above one. Both aggregate controls meet the preregistered 2-percent bound.
The R=199 control's individual round 1 is 1.031611; the specified gate is on
the aggregate, not every individual round. No threshold was changed after
measurement, and all rounds are retained.

`replay_gfni_boundary_screen.py` independently derives the 18 round ratios
and six aggregates from raw ordered samples, without importing the collector
or executing a codec. Normal and optimized Python replay passes rehash all
12 frozen inputs, verify raw stdout/stderr, the eight preflights, 58,785,792
full parity-comparison bytes, guarded Release/sanitizer records, shutdown
snapshots and resource counters. Each also rejects 15 semantic/claim
mutations and checks that a failed control suppresses both decisions.

This identifies a substantial opportunity at exactly the observed deficits:
the existing explicit backend is 54.5% and 47.1% faster than current AUTO,
and 46.2% and 44.1% faster than Leopard1 in this screen. It does not yet make
AUTO faster: the production predicate is unchanged. Next qualify a selector
that adds only R=200/32 KiB and R=199/64 KiB to the existing exact target,
while preserving host, thread, profile, layout, flags, full-output and API
gates. Retain explicit backend requests. A fresh same-binary candidate/control
campaign must also check the existing target and inactive neighbors against
the 2-percent regression bound; fallback/qualification failures, concurrency,
partial outputs, decode and batch exclusions need deterministic coverage.

This diagnostic attempt is consumed. Do not rerun it, pool old samples,
relabel it as selector qualification, or claim broader v19 closure.

The next implementation is tracked in `leopard-79h.38.5.4.17.1`; the parent
remains open. The local immutable evidence bundle is
`.research/leopard-79h/gfni-boundary.COQ16A`, including source snapshots,
both-field archive/driver artifacts, all raw checks/timings, build/resource
logs, pure tests and independent replays. No second-host copy is made under
the user's local-only instruction.

The sealed bundle contains 279 manifest entries (150,580,568 bytes including
the manifest). Its outer `SHA256SUMS` hash is
`cff77de2515b10d316dba74a8f029441c9d4cdc3f49db09ad85c17ed380c1e5f`.
The frozen-copy replay passed locally. Retention stayed below 256 MiB with
all six memory events and swap zero. This paragraph and the matching
structured-result manifest metadata are supplemental to the sealed report
snapshot, so no immutable evidence is rewritten.
