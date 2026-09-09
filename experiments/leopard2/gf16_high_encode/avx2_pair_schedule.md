# AVX2 inverse-pair scheduling prototype

Tracking: `leopard-79h.38.5.4.18.2`. This is an **untimed experimental
correctness milestone**, not a production optimization or measured speedup.
Production codec sources are unchanged.

The preceding [ISA comparison](avx2_isa_screen.md) found Leopard2 explicit
AVX2 ahead of AVX2-only Leopard1 by 16.6–22.8%, but behind original native
Leopard1 at K1000/R200/32KiB and K1000/R199/64KiB. Native remains the product
comparator. These results do not establish generic API overhead as the cause.

## Candidate and actual object code

The production inverse two-way GF16 butterfly reloads four broadcast table
vectors from its stack on each 64-byte iteration. The prototype preserves
the eight existing table vectors, AVX2 instruction ceiling, split transform,
tiling, source-copy policy, scalar tail, forward butterfly, and public API.
It changes only the normal AVX2 variant of `AVX2FF16Butterfly2Prepared<true>`.
Its callers include the existing split/range kernels; there is no new fusion.

Three compile-time variants were inspected with the exact production GCC
13.3 Release recipe, including its existing compiler-GC settings:

| Mode | Change | Inspected loop instructions | Stack references per iteration |
| --- | --- | ---: | ---: |
| 0 | Original implementation | 40 | 4 table loads |
| 1 | Early inverse-y stores, direct x accumulation, streamed nibble expressions | 40 | 4 table loads |
| 2 | Mode 1 plus compiler-only data/accumulator dependency boundaries | 36 | 0 |

All three loops retain eight `vpshufb` operations. Mode 1 is a codegen-negative
result for the spill-removal hypothesis, not a measured performance rejection.
Mode 2 keeps all table vectors in registers. Its empty GNU extended-assembly
statements emit no instructions and are not hardware memory fences. They
constrain the compiler's simultaneous live values and XOR reassociation.
They may nevertheless increase execution latency; instruction counts are
not dynamic work shares, cycles, or throughput estimates.

Early stores are legal because the documented `Butterfly2` contract requires
disjoint x/y shard buffers. Both x vectors and inverse-y values are loaded
before stores; the local vectors are then consumed by the product helper.
The candidate jointly changes store order, accumulation and scheduling, so
even a later gain cannot be attributed solely to table loads.

The overlay is [avx2_pair_schedule.patch](avx2_pair_schedule.patch) with
[avx2_pair_schedule.h](avx2_pair_schedule.h). It is applied only to a private
source copy. Mode 0 produced the **entire byte-identical production object**
`bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d`.
Mode 2's object is
`e0dc9278e24569c3f5112a7e6d8b811aa2a7043d5261bf05a8a4e2ec6b24f8a3`.
All three actual objects passed the no-EVEX/high-YMM/ZMM/ternary/GFNI scan.
No claim about other compilers or a portable MSVC implementation is made.

## Correctness evidence

Mode 2 replaced exactly one member of each pinned, both-field archive.
All other 23 members, including GFNI and AVX-512, were re-compared byte for
byte with their original archives. The sanitizer member was rebuilt with
the existing full ASan/UBSan instrumentation; leak detection stayed enabled.

Release and sanitizer profiles each passed:

- 66,147 pair cases: all 65,535 ordinary multiplier logs at an unaligned
  64-byte vector, plus 612 scalar-tail/vector/large-size boundary cases in
  both directions. Scalar table arithmetic is the independent kernel oracle.
  This does not exhaust every possible input value for every multiplier.
- 256 inverse split-range cases against scalar, including all eight
  zero-skew masks, distances 1/4, and both existing fusion-policy values.
- Eight public-API shapes, six full/prefix/sparse/no-output masks per shape,
  source and output guards, scratch bounds, unaligned buffers, even GF16 tails,
  odd GF16 rejection, and a GF8 case. The candidate AVX2 path executes the
  subset calls; unchanged GFNI supplies the full-output reference.
- Three full-parity/three-loss decode round trips: the actual AUTO-AVX2
  K1000/R199/32KiB neighbor, small GF16 with a vector tail, and odd-byte GF8.
  A separate four-thread test shares initialized tables while independently
  encoding/decoding both fields. This is not a race-detector qualification.
- Four public-encode parity files compared in full to the retained,
  independently linked original Leopard1 oracle: K1000/R200 at 32/64KiB,
  K1000/R199/64KiB and K4096/R512/4KiB. Across both profiles: **69,599,232
  bytes in eight full comparisons**. No new independent Leopard1 decode test
  is claimed.

There were **32 positive native records**, two deliberate benchmark-clock
guard aborts (exit 86 before a clock read), and eight malformed-CLI
rejections. The reused parity driver was linked with the clock-abort wrapper;
only its single-encode `--check` path generated parity evidence. No benchmark
samples were collected. A future timing front end must be separately
qualified, including its full warmup/sample sequence and same-binary switch.

All jobs ran locally under the canonical lock, one capped job at a time,
with swap disabled. Native children had a 30-second CPU limit.

| Job | Peak bytes | Limit |
| --- | ---: | ---: |
| Three serial codegen builds | 175,226,880 | 512 MiB |
| Release/sanitizer correctness builds | 263,520,256 | 512 MiB |
| Native checks | 180,895,744 | 256 MiB |
| Read-only final audit | 64,561,152 | 256 MiB |
| Optimized-Python audit | 66,904,064 | 256 MiB |

All six memory-event counters were zero. The native checks took 16.66 seconds
of wall time; that is a resource observation, not an encode timing result.

The first preparation rejected a source-identity mismatch before compilation:
the retained production bundle had the old two-line comment in `leopard2.cpp`.
The comment-only difference was inspected, that private input was replaced
with the exact committed source, and the complete source inventory was
checked against `3a2f064`. The original failed log is retained. No original
bundle or production file was modified.

[audit_avx2_pair_schedule.py](audit_avx2_pair_schedule.py) executes no codec.
It regenerates disassembly from each actual object, checks the full raw
instruction listing and AVX2 ceiling, identifies the inverse loop, replays
42 native records, re-compares all archive members and parity bytes, and
checks resource envelopes. It passed in normal and optimized Python.
Four focused parser tests passed in both modes. Review provenance is Codex
self-review plus deterministic checks; the user's Claude opt-out applies,
and no independent-model `CONVERGED` gate is claimed.

## Remaining gate

Keep the prototype isolated. Next implement and qualify a same-binary
original/mode-2 control with stable selection before codec execution, then
freeze the actual executables and **commit and push a fresh preregistration
before any timings**. Include native Leopard1, same-path controls, unchanged
neighbors and the actual AUTO-AVX2 neighbor. Preserve the 5% target gain,
2% controls/neighbors and zero-sibling requirements. A longer dependency chain
may cancel the saved loads. No retry of consumed experiments, pooling of old
ratios, change to explicit AVX2 requirements, remote workers or unrelated
process affinity changes is authorized by this result. Production integration
and broader performance qualification remain separate, open work.

Raw workspace: `/tmp/leopard-avx2-pair-schedule.o8wU0I`. The read-only retained
copy and outer manifest are recorded in Beads. This prototype does not close
the parent deficit issue or the overall Leopard2 performance goal.
