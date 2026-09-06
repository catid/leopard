# GFNI source-staging timing filter

Bead: `leopard-79h.38.5.4.11`. Date: 2026-09-06.
Status: **rejected by the completed, clean preregistered diagnostic**.

The correctness milestone `0d0ac30` proves removal of the explicit input-copy
pass at AUTO/GFNI K1000/R200/64 KiB, not a speedup. This new filter asks whether
the existing out-of-place GFNI first stage is faster than the current copy-first
policy. It is not a retry of the exhausted current-versus-Leopard1 `.10` plan.

One new executable handles OFF and ON, with identical code layout and archive.
The only codec call intercepted is the existing source-policy function. There
are no libc-copy wrappers. A bounded 64-record trace remains in both modes;
this is a diagnostic, not a zero-instrumentation production benchmark. The
front end delegates the unchanged existing `current_route_screen.cpp` loop:
one initial encode, four additional warmups and 21 timed public encodes.
It checks all pass arguments and counts after the loop, then prints one trace.
The maximum is 52 records; overflow refuses before delegating to the codec.

The public workload, arithmetic, tiling, output binding, archive and production
sources remain unchanged. Only the previously validated narrow source-policy
predicate can alter a pass. Same-current controls use the identical executable
path, inode and arguments. `--exercise` separately runs 26 clock-free check
workloads to test full trace capacity; it does not model the timing loop's
allocation or warmup behavior.

## Untimed checks

- All 24 Release/sanitizer × OFF/ON × six-cell records match the original
  workload. Twelve Release parity files match every byte of exact Leopard1,
  totaling 122,028,032 bytes.
- Four 26-encode target exercises pass, including both sanitizer modes, with
  exactly 52 matching passes and either zero or 52 changed policies.
- Both 16- and 64-record capacities pass overflow-before-delegation, reset and
  neighbor checks in Release and ASan+UBSan+LSan.
- Six rebuilt default copy probes match all old workload/copy/policy records
  and parity files. Sixteen malformed or invalid-timing requests refuse before
  entering a workload. An inherited `OMP_NUM_THREADS=64`/dynamic setting is
  overridden locally by the driver, without changing another process.
- Seven pure collector tests pass with normal Python and `-O`.

The initial build's final binary-equality assertion failed: changing the guard
from `calls == 16` to `calls >= capacity` changes its branch encoding. Retained
disassembly shows exactly that comparison/jump difference in the wrapper.
All compilation commands succeeded; the failed assertion and original recipe
are retained. The corrected recipe removes that invalid equality expectation;
the six functional comparisons above replace it. The default 16-record limit
and all codec/resource checks remain enforced. Build peak was 88,977,408 bytes
under 512 MiB. All successful native
checks stayed under 256 MiB, maximum 141,807,616 bytes, with all six memory-event
counters and swap zero.

## Fixed decision rule

The committed JSON fixes foureyes CPU22/sibling86, a ten-second passive gate,
one attempt, six cells, three OFF/ON ABBA rounds and three same-current ABBA
rounds per cell: 144 timed processes. The target needs at least 5% aggregate
gain and a gain in all three rounds. All six same-current and five mechanically
unchanged OFF/ON controls must lie within the 2% equivalence band. Otherwise
the filter rejects the candidate or is inconclusive, as specified in the plan.
No partial inference, retry, CPU substitution or movement of other workloads.

The filter never claims confidence intervals, a production promotion, an exact
Leopard1 speed comparison or completion of v19 qualification. At preregistration,
the experiment and parent performance-gap tasks remained open. Review is Codex self-review and
deterministic/sanitizer validation under the user's Claude opt-out, not an
independent-model `CONVERGED` result.

## Completed result

Preregistration `0af2f8d` was pushed before the sole foureyes attempt. All twelve
untimed checks passed. The passive sibling counter stayed at 191797 over
10.000066061 seconds, and every one of the 144 timed invocations observed zero
sibling work. All workload and source-policy records matched their fixed
expectations. The result is `reject_for_this_screen`.

Ratios below are current time / staged time; larger than one favors staging.
The control column compares two identical OFF invocations. These are geometric
means of three round contrasts, not confidence intervals.

| Cell / route | Current / staged | Same-current control |
| --- | ---: | ---: |
| K1000/R200/64 KiB, AUTO/GFNI target | 0.998812 | 0.993443 |
| Same shape, explicit AVX2 | 1.006968 | 1.014013 |
| Same shape, explicit AVX-512 | 0.998230 | 1.001817 |
| K1000/R200/32 KiB, AUTO | 1.001989 | 1.001379 |
| K1000/R199/64 KiB, AUTO | 1.006352 | 0.994118 |
| K4096/R512/4 KiB, AUTO | 1.002449 | 0.998415 |

All eleven **aggregate** controls satisfy the predeclared 2% equivalence band;
the plan did not require each individual control round to satisfy that band.
The target's three ratios are 0.996089, 1.000549 and 0.999804. Its aggregate is
essentially flat and nowhere near the required 5% gain. This does not establish
a statistically significant slowdown, or prove that memory traffic never
matters. It rejects this particular copy-first-policy replacement on this host
and shape. The explicit copy count alone was not a useful speedup predictor.

The standalone stdlib `replay_gfni_source_stage_screen.py` imports no collector
and launches no codec. Normal and optimized Python independently rehashed all
eleven frozen inputs, checked all raw workload/trace records and ordering,
recomputed every round/aggregate contrast using log medians, and reproduced
the decision. Its largest scope peak was 12,333,056 bytes under 256 MiB.
The actual server attempt peaked at 133,509,120 bytes under 256 MiB and exited
zero, with all six memory-event counters and swap zero. Before/after hashes
match on both machines. Production codec/header/CMake differences from
`36dc0c8` remain empty. No policy was promoted.

The 528-file, read-only 203-MiB bundle is retained locally and on ripper at
`.research/leopard-79h/gfni-source-stage-screen.Mp5s60`:

- Outer `SHA256SUMS`: `9198d5f445417a342fb76fd7db7ca4f1df7028ad3bb774e4cec7f415d0b3d75d`.
- Attempt journal: `94aee1bcd75dfc2276d7a7721bda974a26237ac88d4d8cccac614d37c0864294`.
- Scope log: `3bd9f105b6a773b15b3f31fbac1ca4efae028acc21f11cc831ca9e467cf67358`.

The bundle includes the original failed binary-equality expectation, corrected
recipe, all native validation, immutable sources/executables/archive, every raw
server output and trace, independent replays and reference-bundle identities.

The `.11` candidate evaluation is closed as a completed negative experiment;
the parent performance gap remains open. Next is `.38.5.4.12`: attribute actual
current GFNI callback work/cost before selecting a materially new optimization.
A bounded `perf stat -e instructions:u -- true` capability probe failed because
foureyes has `perf_event_paranoid=4`; no kernel setting or capability was changed.
Start with clock-free software callback counts, not another source-staging or
cache-block retry. The exhausted `.10` Leopard1 plan and v19 requirements remain
unchanged; this result makes no current-Leopard1 speed or gap-closure claim.
