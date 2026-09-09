# AVX2 adjacent scheduling: completed negative screen

Bead `leopard-79h.38.5.4.18.3.3`. **No production promotion.** The candidate
reduced spills and instructions but delivered only about 0.6% more throughput
than original Leopard2 on the four targets. Native Leopard1 remains faster.
The single attempt is consumed; no retry, pooling, trimming, or relaxed gate.

The [fixed method](avx2_adjacent_screen_method.md), plan, collector and replayer
were committed **and pushed** as
`0b3ff749f5e14d2feb278e702365cba980c14b7b` before benchmark clocks. Those files
remain unchanged. Actual execution completed 36 untimed preflights and all
690 timed processes / 57,960 spans, with all identities, stopped-service
conditions and zero-sibling checks passing. No new codec build was needed.

## Decision and target results

The fixed decision is `reject_unchanged_neighbor`, **not**
`inconclusive_controls`: all 31 cross-process and 124 within-process aggregate
stability controls pass. GF8 cell7 original/ON is 1.0619582165, outside the
symmetric 2% unchanged-path interval. This is a positive throughput shift,
not a slowdown. Its cause has not been isolated; code layout or ASLR is not
established as the cause.

Independently, every target falls short of all 5% performance gates, including
the paired OFF/ON gates, the direct original-Leopard2 gate, and the native-L1
gate. A ratio below 1 for native/ON means native is faster. Ratios below are
cost ratios, equivalently candidate throughput divided by comparator throughput.

| Target (K1000) | Original / ON | Paired 0110 | Paired 1001 | Native / ON | Native throughput lead |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cell0: R200, 32 KiB | 1.006119 | 1.008235 | 1.008667 | 0.928423 | 7.710% |
| Cell1: R199, 64 KiB | 1.006520 | 1.006804 | 1.006565 | 0.973040 | 2.771% |
| Cell2: R200, 64 KiB | 1.005977 | 1.006804 | 1.006877 | 0.977280 | 2.325% |
| Cell8: cell0, one-item batch | 1.005904 | 1.007435 | 1.009194 | 0.928740 | 7.673% |

All corresponding target original/ON and paired round ratios are above 1;
all target native/ON rounds are below 1. Native cell8 uses its ordinary encoder,
not a native batch API. All five neighbors remain in the record:

| Neighbor | Original / ON | Paired 0110 | Paired 1001 | Fixed neighbor gate |
| --- | ---: | ---: | ---: | --- |
| Cell3: AVX2 K4096/R512/4 KiB, affected | 1.007869 | 1.013116 | 1.012403 | Pass |
| Cell4: AUTO AVX2 K1000/R199/32 KiB, affected | 1.007703 | 1.006513 | 1.009578 | Pass |
| Cell5: explicit GFNI, unchanged | 1.003904 | 1.000065 | 0.997333 | Pass |
| Cell6: AUTO GFNI, unchanged | 1.001615 | 1.000746 | 0.998901 | Pass |
| Cell7: explicit AVX2 GF8 K17/R7/64 B, unchanged | 1.061958 | 1.000854 | 0.999720 | Fail: original / ON |

GF8 original/ON rounds are 1.0508154788, 1.0784286943, 1.0568281969.
Original/OFF, diagnostic only, is 1.0408742603 (rounds 1.0444405128,
0.9560463835, 1.1293593060). Runtime OFF is not pristine production, as
already documented in qualification. Neither its overhead nor the prior
inverse-only scheduling gains can be multiplied into this direct result.

The following individual control rounds exceed the symmetric 2% interval and
are retained, even though the preregistered aggregate gates pass. Round indices
are zero-based; no sample or round was excluded.

| Cell | Control | Round | Ratio |
| --- | --- | ---: | ---: |
| 0 | same native | 2 | 1.0250409749 |
| 7 | same OFF | 0 | 0.9730451203 |
| 7 | same OFF | 2 | 1.0646342849 |
| 7 | same original | 2 | 0.9584841442 |

## Validation and retained evidence

Collector-free raw and read-only-copy replays pass normally and with Python
`-O`. All four output hashes equal
`2641966763d332a809e3f95a22e0da09ecf61a77f53a0bd08d8405e253c7784b`.
They reconstruct all 195 aggregate and 585 round ratios from raw outputs,
verify all 40 frozen inputs and the qualified shared Release driver, and check
the exact resource/isolation records. The [machine-readable result](results/avx2_adjacent_screen_20260909.json)
retains all aggregate and round ratios, not just this table's selected columns.

| Scope | Peak bytes | Cap | Outcome |
| --- | ---: | --- | --- |
| Timing | 149,774,336 | 256 MiB | exit0, all six events0, swap0 |
| Raw normal/optimized replay | 36,732,928 | 256 MiB | exit0, all six events0, swap0 |
| Retention and sealed normal/optimized replay | 65,667,072 | 256 MiB | exit0, all six events0, swap0 |

Preparation retained its 36 clock-free checks, 20 pure tests per Python mode,
initial oversized-JSON fixture failure and initial source snapshots. The
separate earlier public-qualification retention hit 6,738 memory `max` events
without OOM or swap; it is not relabeled all-zero or part of this timing scope.
Correctness qualification in `bd4c175` and `a6634c7` remains unchanged, not rerun.

Raw root: `/tmp/leopard-adjacent-screen.04vrZl`.
Read-only local bundle: `.research/leopard-79h/avx2-adjacent-screen-negative.Ap1qSk`,
1,587 files / 14,371,671 bytes, plus its outer manifest. Manifest SHA-256:
`d5e7f9794132973bc0d49d68bcfa38573e98f3cfa3d07a28dff3d7957b7179b4`.
Every retained file was compared to the raw copy where applicable and the full
manifest rechecked after replay. Delivery logs:
`/tmp/leopard-adjacent-screen-delivery.4XPb9Y`.

This is a local `work` result, not an ISA-matched causal ablation, confidence
interval analysis, cross-host conclusion, or authoritative v19 campaign.
Codex self-review and deterministic/adversarial checks replace Claude under
the user's opt-out; no independent-model `CONVERGED` claim. No Claude, subagents,
SSH workers, unrelated process/affinity changes, or host-setting changes.

## Next direction

The spill-scheduling candidate stays experiment-only. The next useful avenue
is a materially different reduction in GF16 multiplication work: first test
whether the actual Cantor basis supports a GF256 tower representation and
three-product multiplication. This is an untimed algebra/code-generation
hypothesis, not an implemented optimization or predicted speedup. Conversion
cost, extra XORs, table setup and whole-codec compatibility must be evaluated.
The subsequent [untimed algebra and isolated-codegen proof](tower_algebra.md)
now passes; whole-transform integration and conversion costs remain open.
Any later timing requires separate qualification, method review and a new
committed-and-pushed preregistration. The broader performance goal stays open.
