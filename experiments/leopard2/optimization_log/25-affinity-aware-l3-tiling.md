# 25 — Affinity-aware GF16 L3 tiling

**Disposition: SHIPPED — hybrid cache policy; LOW side-256 diagnostic retained
but not promoted**

## Question

The original GF16 AVX2 tiling policy was calibrated on a 32 MiB L3 domain.  This
host also exposes a 96 MiB L3 domain, where the same fixed policy can split a
working set that would otherwise fit.  The experiment asked whether a context
can use the smallest L3 visible through its creation-time affinity without
regressing either cache domain.

This is execution geometry only.  It does not change the field, profile,
coordinate map, parity bytes, or erasure mathematics.

## Candidates

- **F**, fixed policy: the pre-experiment 32 MiB calibration.
- **D**, proportional policy: clamp detected L3 to 32--96 MiB, grow the retained
  live-set target to half of L3, and move the generic crossover with cache
  capacity.
- **H**, hybrid policy: clamp detected L3 to 32--96 MiB, retain the calibrated
  16 MiB execution target, and use `max(effective L3, 64 MiB)` as the generic
  crossover.  Existing side-256 and side-512 entries remain cache-gated and
  their requested pass widths scale with effective L3.
- **S**, LOW side-256 candidate: exercise a 16 KiB pass on the measured low-rate
  side-256 cell.  The campaign used exactly one nine-erasure mask.  Production
  does not extrapolate that result to the other nine-erasure schedules.

Only an explicitly requested AVX2 context consumes the affinity-derived GF16
policy.  AUTO, GFNI, scalar, SSSE3, and AVX-512 retain the conservative 32 MiB
calibration until they have their own evidence.

## Authoritative method

All measurements used the fail-closed
`experiments/leopard2/l3_tiling/run_abba.py` runner:

- explicit GF16 and AVX2, one thread;
- three ABBA rounds with nine retained samples per invocation;
- the canonical build/test/timing lock and a physical-pair lease;
- clean, exact-commit worktrees and lane-owned frozen executables;
- build-source attestation and SHA-256 verification before and after timing;
- zero observed non-idle sibling samples; and
- matching input, parity, and recovered-output digests in every completed cell.

CPU 4 with reserved sibling 20 is on the 96 MiB L3 domain.  CPU 14 with
reserved sibling 30 is on the 32 MiB domain.  A reported ratio is left-lane
time divided by right-lane time, so values above one favor the right lane.

## Rejected proportional policy

The F-to-D campaign was valid but rejected.  Its proportional pass grew too
wide and lost the headroom needed by the transform and surrounding state.

| Cell | Operation | F/D | 95% CI | D slowdown |
| --- | --- | ---: | ---: | ---: |
| high live set 96 MiB | encode | 0.7928x | 0.7721--0.8141 | 26.1% |
| high live set 96 MiB | decode | 0.8285x | 0.8070--0.8507 | 20.7% |
| high live set 256 MiB | encode | 0.9444x | 0.9229--0.9664 | 5.9% |

D did improve high side-512 64 KiB encode by 1.1095x, but one favorable region
does not justify those neighboring regressions.

## Promoted hybrid policy

The F-to-H campaign covered eight boundary and calibrated-side cells.  It was
accepted with no credible regression over 2%; all artifact and output digests
matched.

| Cell | Encode F/H (95% CI) | Decode F/H (95% CI) |
| --- | ---: | ---: |
| high live set 96 MiB | 0.9997x (0.9902--1.0092) | 1.0016x (0.9972--1.0059) |
| high live set 256 MiB | 1.0027x (0.9917--1.0139) | 0.9991x (0.9976--1.0007) |
| low live set 128 MiB | 0.9970x (0.9774--1.0170) | 0.9994x (0.9947--1.0041) |
| low live set 256 MiB | 0.9959x (0.9860--1.0060) | 0.9983x (0.9910--1.0056) |
| high side 256, 64 KiB | 1.0285x (1.0235--1.0335) | 1.0245x (1.0202--1.0289) |
| high side 512, 64 KiB | **1.1066x** (1.0936--1.1198) | 1.0222x (1.0149--1.0296) |
| high side 512, 128 KiB | 1.0052x (0.9966--1.0139) | 1.0425x (1.0300--1.0552) |
| high side 512, maximum loss | **1.0946x** (1.0661--1.1239) | 1.0312x (1.0293--1.0332) |

The hybrid therefore removes D's large-cache regressions while retaining the
calibrated side wins.

## LOW side-256 diagnostic

On CPU 14, D-to-S for K=200, R=800, loss count 9, and 64 KiB shards measured:

| Operation | D/S | 95% CI |
| --- | ---: | ---: |
| encode | 1.0089x | 0.9913--1.0268 |
| decode | **1.1855x** | 1.1455--1.2268 |

The encode result is neutral and the measured decode clears the 5% numerical
threshold.  However, the deterministic seed exercised only missing originals
`1,9,22,31,62,112,125,130,148`.  Sampling 1,000 other nine-erasure masks found
39 distinct pruned schedules and a 1,740--1,960 butterfly range.  Pattern shape
is therefore a material unmeasured variable.  Production retains the generic
LOW policy until a counterbalanced multi-pattern campaign qualifies a useful
region; it does not key a permanent dispatch rule to one arbitrary bitmap.

## Implementation and scratch contract

At context construction, an explicit AVX2 request snapshots the creating
thread's affinity and uses the smallest readable L3 among allowed CPUs.  An
unknown, malformed, smaller, or unsupported topology falls back to 32 MiB.
The observation is quantized to MiB and capped at 96 MiB.  The immutable policy
then uses:

- effective L3: 32--96 MiB;
- retained live-set target: 16 MiB; and
- generic tiling threshold: `max(effective L3, 64 MiB)`.

Balanced byte passes divide the aligned prefix as evenly as possible.  Scratch
reserves the widest resulting pass rather than the requested nominal pass, so
no later pass can overrun a slot.  Known plans use their retained row count;
the pattern-independent codec query decides tiling with 2P or 2T base rows
while allocating its conservative maximum rows.

The schema-v4 operation model and production-linked probe independently report
plan and codec pass counts and widest-pass bytes, cover the GF8 side table as
well as GF16 cache policy, and preserve actual batch execution row counts
separately from conservative public scratch rows.

## Correctness and validation

The final integrated tree passed:

- 23/23 schema-v4 unit tests under both ordinary Python and `python -O`;
- 397 operation-model self-test checks under each Python mode;
- the production-linked scratch, selector, and cache-policy cross-check;
- two independent 100,000-case randomized tiling/scratch invariant sweeps;
- 104/104 GNU Release CTests;
- the current focused Clang ASan+UBSan gate (5/5), following a full integrated
  97/97 sanitizer run;
- strict GCC and Clang builds plus GF8-only and GF16-only 7/7 reduced-field
  tests; and
- 40/40 lab containment repetitions and 5/5 fake-counter fixture repetitions
  under full unrelated CPU contention.

The final audit found and resolved three evidence/model-scope issues before
landing: the LOW single-mask extrapolation was withdrawn, impossible GF8
parents are rejected by the scratch helper, and retained executable recipes
must contain exactly one nonempty command.  The build/API and post-fix reviews
found no unresolved defect.  The production archive contains no test-only
hooks, and the private C++ build-tree introspection used by the probe does not
enter the public C ABI or installed header.

## Evidence identity

The compact, tracked
`experiments/leopard2/l3_tiling/results/checkpoint.json` records every campaign
digest, lane commit/tree/artifact identity, decisive interval, and the one
fail-closed runner bootstrap result.  Raw bundles remain generated artifacts
because they include large retained-sample and process-observation streams.

The promoted H lane is commit
`cb2ef79cc180166326acd80a3d798392ce5157fe`, tree
`9d48810b25b227db8a27e3081be6ccd48e3276a6`.  Its frozen benchmark SHA-256 is
`536e06555f60de490a1dfe74295aa7610d27d3c6987cfa3596370758c1748a86`.
The reviewed production/model integration is commit `1a2e7c2`; the positive
LOW diagnostic is intentionally absent from its dispatch policy.

## Reproduction

From clean worktrees at the checkpointed commits:

    python3 experiments/leopard2/l3_tiling/run_abba.py \
      --lane F /tmp/leopard-l3-F \
        ec927d8ae13c0ad2ddf236c98d21165dcb876744 \
        build/release/bench_leopard2 \
      --build-lane H /tmp/leopard-l3-H \
        cb2ef79cc180166326acd80a3d798392ce5157fe \
      --comparison F H --cpu 4 --sibling 20 \
      --output /tmp/leopard-l3-evidence-FH

The D-to-S campaign uses CPU 14/sibling 30 and
`--cells low-side-256-loss9`.  The runner rejects mutable or mismatched
binaries, source identities, sibling contamination, digest disagreement,
truncated subprocess output, and incomplete campaigns rather than emitting a
performance conclusion.
