# GFNI final inverse accumulation experiment

Bead: `leopard-79h.38.5.4.14`. Date: 2026-09-07.
Status: focused correctness and structural milestone passed; experiment-only,
untimed. The extra sanitizer-archive scan was subsequently resolved as outside
the existing production audit's scope; see `gfni_terminal_sanitizer_policy.md`.
The terminal evaluation and parent performance goal remain incomplete.

The preceding first-stage forwarding screen measured a directional 1.033262x
OFF/ON ratio, below its 5% gate, and was not promoted. This experiment leaves
that first stage unchanged and instead addresses the terminal inverse stage
of later source blocks, which still materializes its first layer in scratch
before accumulating its second layer.

## Distinction from prior experiments

The retained `.5` negative record says the old temporary source patches were
reverted, so an exact source audit of those candidates is unavailable. Its
direct AVX-512 terminal candidate was 1.65% slower and its packed-ZMM version
22.6% slower, with nibble-table rebroadcast/lane reduction costs identified.
Those are historical diagnostic negatives, not current isolated GFNI results.

The shipped GFNI backend already uses compact affine multiplication tables;
this experiment does not repeat that packing optimization. The Ops interface
has four-way GF8 accumulation but no equivalent GF16 entry. The new GFNI-only
entry uses the existing affine helper and changes the terminal dataflow,
without introducing nibble tables, ZMM packing, smaller tiles or source staging.

## Source, activation and build boundary

Production source, headers and CMake remain unchanged from base codec
`36dc0c8f66604b8d974468e687c51e6183ecb61d`. The field overlay
`gfni_terminal.patch` has SHA-256
`a731146b01c5013e084370862923ace33d1b49eb99a8cfc68780a3478d558ea5`.
It propagates a once-per-pass default-off boolean through private encoder
helpers and replaces only `accumulate_split_stage`. Nonaccumulating/first-stage
dispatch, source copies, tiling, forward arithmetic and low-profile calls stay
unchanged. There is no public API or Ops ABI change.

The predicate matches GFNI, K1000/R200, all 200 outputs requested, side256,
32-KiB execution bytes, 64-KiB source-policy bytes, and a present sparse-plan
descriptor with zero blocks. This is an aligned-prefix predicate: explicit
GFNI at 65,538 public bytes also selects two prefix passes and leaves the
two-byte tail alone; AUTO at that neighbor stays AVX2 and does not select it.
Both cases are tested. The global diagnostic trace is bounded to sixteen
passes and used only with one thread; concurrent diagnostic use is not claimed.

`gfni_terminal_backend.cpp` includes the unchanged production GFNI translation
unit and adds `LeoGFNIFinalAccumulate`. Four rows pass through both inverse
layers using the existing `AVX2FF16MultiplyAddPair` and packed GFNI tables,
then XOR into four disjoint accumulator rows. Inputs remain readonly. A compact
even residual below 64 bytes uses four bounded stack buffers, the mature split
helper and XOR stores. Empty byte/distance calls touch no pointers. Odd GF16
payload support and overlapping input/accumulator ranges are not claimed.

Full dual-field Release and ASan+UBSan archives come from the pinned
`gf16-current-route-failed.STc10h` bundle. Only FF16 and GFNI members are rebuilt,
using their recorded original compile commands with source/output/include
paths relocated. Both portable-field SIMD exclusions and GFNI's
`-mavx2 -mgfni -mno-avx512f` remain. All other 22 members and archive order
are byte-identical, rechecked after validation. Sources, binaries and archives
are lane-owned and frozen; their build hashes were reverified after all checks.

The check driver accepts only `--check cell --terminal=0|--terminal=1
[parity_file]`. It delegates the unchanged public workload and records pass
selection, the qualified Ops callback histogram, and separately wrapped
terminal-kernel calls. These counters are diagnostic only, with no timing
path. The workload JSON's commit identifies the base; the overlay, added
kernel and frozen artifact hashes identify the actual candidate.

## Verified result

| Target structure | OFF | ON |
| --- | ---: | ---: |
| First-stage ordinary inverse pairs | 2,000 | 2,000 |
| Terminal ordinary inverse pairs | 768 | 0 |
| Terminal accumulating inverse pairs | 768 | 0 |
| External terminal-kernel calls | 0 | 6 |
| Total Ops callbacks (external calls excluded) | 4,132 | 2,596 |

The six kernel calls each cover 64 four-way groups, totaling 384. Every other
callback bucket and pass field stays identical. All five other fixed routes
retain the baseline histogram in both modes. At source level, removing a
four-row first-layer scratch store and reread accounts for
`384 * 4 * 32768 * 2 = 100663296` bytes (96 MiB) per target encode. This is
logical scratch materialization avoided, **not measured memory traffic,
bandwidth, elapsed time or a speedup**.

- All 24 Release/sanitizer × OFF/ON × six-cell workload records match the
  original. Twelve full Release parity files match standalone exact Leopard1
  `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`, totaling 122,028,032 bytes.
  Sanitizer traces are byte-identical to their Release counterparts.
- Each build passes 512 kernel cases against explicit scalar inverse layers
  followed by XOR: distances 1/4 and selected distance64 cases, fifteen lengths
  from zero to 65,598, exact allocation ends, offsets zero/17, all eight
  zero-skew masks, readonly inputs and repeat-accumulation cancellation.
  Two additional null-pointer empty operations pass.
- Each build passes fourteen directed public cases covering unaligned
  exact-end buffers, repeated immutable inputs, partial/omitted parity,
  two-byte tails, route/shape neighbors, native odd rejection, short scratch
  and overlap rejection. Outputs, omitted rows, inputs and canaries stay exact.
- Seventeen negative predicates, both modes, descriptor identity, reset and
  atomic sixteen-pass overflow rejection pass. Fourteen malformed/timing
  requests refuse before entering the workload.
- Ten pure structural/resource model tests and retained-only replay pass
  normally and with Python `-O`. Replay independently reconstructs the
  original traversal after undoing only the terminal substitution; it imports
  no experimental C++ code and executes no codec or collector.
- The 54 final native matrix/kernel/public scopes peak at 158,613,504 bytes
  under 256 MiB. Final archive/driver builds peak at 278,343,680 bytes and
  unit builds at 89,927,680 bytes under 512 MiB, one compiler at a time.
  All memory-event counters and swap are zero. No broad oracle was run.

## ISA findings and retained failures

The first candidate passed the functional matrix but failed the unchanged
Release ISA scanner: GCC packed local pointer arrays using `vpinsrq`, outside
the reviewed GFNI mnemonic allowlist. Named scalar pointers and explicit
load/accumulate calls remove that incidental packing. No ISA flag, allowlist
or field arithmetic was changed. Both full archives and all native checks were
rebuilt/rerun. The entire rejected first lane is retained in `rejected-v1/`.

The final full Release archive passes the existing scanner. Its all-nonzero
terminal hot loop uses VEX GFNI affine arithmetic, stores only the four
accumulators, and has no vector spills in the inspected loop body. Stack
metadata and compact-tail storage still exist. Assembly is not timing evidence.

An additional scan of the full sanitizer archive fails because the **unchanged**
`Leopard2BackendAVX512.cpp.o` contains ZMM stores. The original pinned sanitizer
archive fails identically; their offending member is byte-identical
(`07c42bb27aae1d67244c2841d2833a8a2217dbf25b0128536010dc97c758c57d`).
Observed `vmovdqu8` stores target ASan shadow offsets around `0x7fff8000`,
consistent with sanitizer instrumentation; this is an inference, not a completed
compiler diagnosis. The modified FF16/GFNI sanitizer members pass separately.
That narrow pass does not make the whole archive pass. The initial milestone
treated this as an unresolved gate; follow-up `leopard-79h.38.5.4.15` verified
that existing project policy deliberately separates sanitizer correctness
from unsanitized Release ISA auditing. See the correction report above.
No scanner or width policy was relaxed, and the raw failures remain retained.
The largest ancillary scan scope is 188,465,152 bytes with zero memory events
and swap, including retained nonzero ISA exits.

These are archive instruction scans only: compiler-metadata qualification and
scanner self-tests were not rerun. Existing v19 requirements are not satisfied
by these checks. The first provenance-script run also caught an overly broad
Git pathspec that included historical experiment additions; its source/log are
retained. Explicit top-level globs then verified the actual production boundary.

## Frozen identities and retention

| Artifact | SHA-256 |
| --- | --- |
| Added GFNI kernel source | `e53fffd31283cd96b73f832bb070a8bbb0fbbcbee41259db7bca1bd1753a0596` |
| Release archive | `bf189e85b309eda5dea55c8a47e5b78e38779822d3556f86f85c457e081be127` |
| ASan+UBSan archive | `ee18608af8c0d32d1ad20e35d4670bc03572e7cd7e206d77d3fb00a6feef342e` |
| Release check driver | `9219486c1ca4b0d781cd828b20020cc64a451bb4b588c231bddd1aa6881eb587` |
| Sanitizer check driver | `69c600f5352a0fc319777e5805760d6cfea306dcc26495b4dfd5794d42aca3a8` |
| Release directed test | `3e4722270607a3dedbfec4522fc0ecedddb47e98a1c40abf769081431ac85ccc` |
| Sanitizer directed test | `3dbdc1338571b4df5aaa124334392b35b4df42c6582309f215f685989b9fd91e` |

The read-only bundle `.research/leopard-79h/gfni-terminal.nEDGNa` is retained
locally and on ripper with both lanes, complete source/artifact provenance,
raw parity/traces, build/test/ISA failure logs, recipes and replay dependencies
authenticated by `SHA256SUMS`. The referenced baseline bundle's outer manifest
is `e83a217e680a23088e53208d9ba7747742e7918a1377e74efe9c35399b4ba342`.
Replay its frozen `replay/verify_gfni_terminal.py` using the experiment bundle
and baseline `preflight` paths in the established 256-MiB/no-swap scope.

## Next gate

No timing has been attempted. A separate driver without per-callback wrappers,
qualified capacity for its warmup/measurement loop, and a newly preregistered
frozen-binary OFF/ON screen are required. Keep first-stage forwarding OFF for
this isolated contrast. Retain the 5% target and 2% control gates, same-OFF
layout controls and zero-sibling requirement. Do not retry the exhausted `.10`
Leopard1 attempt or `.13` screen, multiply old ratios, or infer a production
speedup from the callback/traffic model. Any later combined experiment needs
its own evidence. Production integration additionally needs broader correctness,
unsanitized Release ISA/metadata, v19 and exact-Leopard1 qualification gates.

Review provenance is Codex self-review plus deterministic/adversarial/sanitizer
checks under the user's Claude opt-out, not independent-model `CONVERGED`.
No Claude/API, subagent, unrelated workload, kernel setting or affinity change
was used. The parent performance objective remains active.
