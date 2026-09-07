# Combined GFNI inverse-boundary experiment

Bead: `leopard-79h.38.5.4.16`. Date: 2026-09-07.
Status: initial correctness and structural checks pass; **not timed or promoted**.

The separate first-stage and terminal screens observed 1.033262x and
1.042655x, respectively, below the fixed 5% threshold. This new four-state
experiment tests their composition. It does not add or multiply those ratios,
retry an exhausted campaign, or change any qualification threshold.

## Change and scope

`gfni_combined.patch` overlays only the frozen base `36dc0c8` field source.
A default-zero unsigned pass mask is propagated through private encoder
helpers. Bit 1 forwards nonaccumulating distance-one groups to the existing
GFNI range entry; bit 2 selects the unchanged `.14` accumulating kernel only
at the final later-message-block stage. Neither kernel is rewritten here.

The predicate remains exact: qualified GFNI, K1000/R200, all 200 outputs,
transform side 256, a present sparse descriptor with zero blocks, 32-KiB
pass bytes, and 64-KiB source-policy bytes. Source copying, zero filling,
tiling, forward transform, backend selection, public API and Ops ABI remain
unchanged. Low-profile calls retain the default-zero mask. Explicit GFNI at
65,538 bytes selects the two aligned prefix passes only; its tail is unchanged.
AUTO at that size and partial-output cases remain unselected.

The mode hook and bounded trace are **single-thread diagnostic machinery**,
not a concurrent production API. The check executable refuses timing requests.
Its callback and terminal linker wrappers are correctness instrumentation;
they must not be used for performance measurements.

## Measured structural counts, not CPU-time shares

The independent replay starts with the original traversal model, derives
the two substitutions, and compares every typed field and bucket. It does
not normalize boolean counts, merge duplicate buckets, or ignore zero rows.
Bucket ordering alone is irrelevant. The combined mode has no `ifft2` bucket;
that absence is explicitly tested.

| Mode | First stage | Terminal | Ops callbacks | External terminal calls |
| --- | --- | --- | ---: | ---: |
| 0 | Original | Original | 4,132 | 0 |
| 1 | Fused | Original | 2,632 | 0 |
| 2 | Original | Fused | 2,596 | 6 |
| 3 | Fused | Fused | 1,096 | 6 |

These are counts for the K1000/R200/64-KiB GFNI public workload. The first
substitution replaces 2,000 pair callbacks with 500 four-way range callbacks.
The second removes 768 ordinary and 768 accumulating pair callbacks, replacing
them with six external calls covering 384 four-way groups. All five fixed
neighbor/backend cells retain their original callback schedules in every mode.
Fewer calls do not imply a proportional speedup or a measured traffic reduction.

## Validation evidence

- All 48 fixed records (six cells, four modes, Release and sanitizer builds)
  match the original workload records. All 24 full Release parity files match
  independently linked exact Leopard1 `6e5725eb`, totaling 244,056,064 bytes.
  Sanitizer traces match Release byte-for-byte.
- Per build, the unchanged first-stage scalar-pair fixture passes 480 cases
  and the terminal scalar-layer/XOR fixture passes 512. These cover exact
  allocation ends, unaligned pointers, zero skews and tails; the terminal
  fixture also checks readonly inputs and repeated-accumulation cancellation.
- Fourteen directed public cases per build exercise all four modes, varied
  scratch/output initialization, partial outputs, adjacent K/R/byte sizes,
  explicit backends and native odd-length refusal. Target short-scratch and
  source/output-overlap errors leave inputs/outputs unchanged in every mode.
- Seventeen negative predicate neighbors, all four masks, atomic 16-pass
  overflow/reset behavior and nine mode-parser guards pass. Twenty-four
  malformed or timing requests refuse before entering the public workload.
- Ten new pure model tests plus ten original-model/resource tests pass both
  normally and with Python `-O`. Retained-only replay agrees in both modes;
  it executes no codec and imports no timing collector.
- The unchanged whole-archive Release ISA checker passes. Sanitized builds
  pass ASan/UBSan/LSan correctness. Existing project policy separates these
  checks; no additional whole-sanitizer ISA gate is imposed or relabeled.

All 80 native correctness scopes stay below 256 MiB, peaking at 157,872,128
bytes. Ancillary checks peak at 48,529,408 bytes. The serialized build peaks
at 172,507,136 bytes under 512 MiB. All six memory-event counters for those
builds/checks and swap are zero. No monolithic API test, broad controller or
subagent was run.

## Provenance and limitations

Lane-owned copies of the complete `.14` dual-field archives replace only
`LeopardFF16.cpp.o`, using its recorded portable-field compiler commands.
Every other archive member (23), including the GFNI kernel member, and member
order are byte-identical before and after checks. All frozen build inputs
and executables are rehashed afterward. The exact patch applies cleanly to
the separately extracted baseline and reproduces the compiled field diff.
Production root source/header/CMake differences from `36dc0c8` remain empty.

- Field patch: `c94ac12d22b0d83643efbe69604f9b1b74472f5c44d676ad390b0eedb4b26cca`.
- Release archive: `aa3621ca58da46f22bba976f89f258ddc1317f72e74cc2a90602e909856d2c8c`.
- Sanitizer archive: `95736e841c06e90991ed9b0b8c85b31865b815f9438f68adbedd38e30c40a7b6`.

The machine-readable result is `results/gfni_combined_correctness_20260907.json`.
Raw source, scripts, build logs, binaries, parity, traces and resource logs
originate in `/tmp/leopard-gfni-combined.EMX869`. The read-only bundle
`.research/leopard-79h/gfni-combined.7FNBT8` includes the baseline parity and
replay dependencies; its outer manifest and second-host verification are
recorded in Beads. Original baseline and separate `.13`/`.14` bundles remain intact.

Post-copy clarification: the separate artifact-copy scope reached its 256-MiB
cap and recorded 3,332 `memory.max` events, with zero OOM, kill or swap events.
Its bytes and manifest verified, but it is not an all-zero resource result.
This copy log is retained separately; no codec was executed in that scope.
The bundled report predates this clarification. Native correctness/build
resource results above are unchanged.

Review is Codex self-review plus deterministic/adversarial/sanitizer checks,
under the user's Claude opt-out, not independent-model `CONVERGED`.
This is narrow experimental qualification, not full build-metadata/v19 or
production correctness qualification, and not a Leopard1 speed comparison.

The next phase needs a separately qualified timing front end with no
per-callback wrappers, a capacity-64 trace check, and a fresh preregistered
balanced four-state comparison. Retain the 5% target, 2% aggregate control
band and zero-sibling-work gates. Any selected winner still needs hook-free
production and exact-Leopard1 qualification. The task and parent goal stay open.
