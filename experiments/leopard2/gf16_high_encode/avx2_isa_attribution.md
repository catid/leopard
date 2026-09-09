# Explicit AVX2 GF16 attribution: native Leopard1 is not AVX2-only

Bead: `leopard-79h.38.5.4.18` (open). Local host `work`, 2026-09-09.

## Finding and scope

The completed [post-Slipgate diagnostic](current_route_post_slipgate.md)
measured Leopard2 explicit AVX2 / original native Leopard1 throughput at
`0.9697681558960904` for K1000/R200/65536 bytes. That is about 3.12% higher
throughput for the native Leopard1 configuration. It is a valid product/build
comparison, **not an ISA-matched comparison**.

The retained native Leopard1 GF16 object actually contains EVEX instructions,
YMM registers above 15, and ternary logic. Its `IFFT_DIT4` AVX2-intrinsic branch
loads `%ymm30` at object offset `0x11f9` and uses `vpternlogq` at `0x139f`.
The original source's `CpuHasAVX2` branch is compiled with `-march=native`;
the name of that source branch does not impose an AVX2 compiler ceiling.
Leopard2's explicit AVX2 object is compiled with `-mavx2 -mno-avx512f` and
contains none of those instructions. Its object bytes are identical between
the measured `36dc0c8` codec and production integration `3a2f064`.

Static counts below cover the entire named function, including alternate
branches, **not executed instruction counts or time shares**:

| Leopard1 `IFFT_DIT4` build | EVEX instructions | Instructions using YMM16–31 | Ternary-logic instructions |
| --- | ---: | ---: | ---: |
| Original native | 390 | 140 | 88 |
| New AVX2-only attribution | 0 | 0 | 0 |

The native GF16 object also contains ZMM code elsewhere; `IFFT_DIT4` itself
has no ZMM instructions. Do not attribute those other instructions to the
measured hot loop or equate total object sizes/instruction counts with speed.
The restriction is part of Leopard2's existing portable ISA contract; it must
not be weakened or silently replaced with GFNI/AVX-512 for an explicit AVX2
request.

This evidence establishes an ISA/build-policy confound. It does **not** measure
its contribution to the observed 3.12% difference, prove API overhead, or
establish a new Leopard2 improvement.

## Attribution comparator and correctness

The existing, default-off `LEO_MAIN_PURE_AVX2=ON` comparison profile was used;
no Leopard1 or Leopard2 codec source was edited. All eight original source
and header files were checked against Git revision
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` and the retained native source copy.
The four-member, both-field archive uses GCC 13.3 and the existing effective
Release flags `-g -O0 -O3`, with
`-march=x86-64 -mtune=generic -mavx2 -mno-avx512f -Wall -Wextra -fopenmp`.
Compiler-only GC parameters remain `10/4096`; no field or optimization pass
was disabled. Native-to-pure comparison consequently changes compiler target
and tuning policy, not just one instruction family in isolation.

All four pure archive objects were disassembled and had zero EVEX,
YMM16–31, ZMM, ternary-logic, or GFNI occurrences. This is the reported
instruction-family scan, not a new general-purpose ISA certification tool.

The new `avx2_isa_check.cpp` has no clocks or measurement mode. Eleven native
executions checked source preservation and outer allocation guards:

- Pure AVX2: all six original diagnostic cells, including their three
  repeated K1000/R200/64KiB cases, plus GF8 K100/R20/64B and K12/R3/4096B.
- Original native: the target and the two additional GF8 cases.
- Nine full file comparisons: **74,134,784 bytes**, including all six pure
  GF16 outputs against their retained original native parity files. The
  retained file hashes were verified against the original pinned manifest.
- Fourteen malformed/timing CLI invocations were rejected. The disassembly
  parser's five tests passed under normal Python and `python -O`.

Archive/source/driver pins were checked before and after native executions.
These are focused Release parity checks, not full sanitizer validation of
the new pure-L1 build, decode coverage, or a production promotion gate.
Original-source GCC array-bounds warnings around the legacy skew-pointer
pattern remain visible in the build log; the exact baseline was not patched
or described as warning-free.

## Artifact identities

Raw local workspace: `/tmp/leopard-avx2-attribution.HbZGYl`.
The complete retained copy and outer manifest identity are recorded in Beads.

| Artifact | SHA-256 |
| --- | --- |
| Original native L1 archive | `3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1` |
| Original native L1 GF16 object | `42079dfa8d3d8ae10e9a545d40e54aa9ae44f18ea4dcb5318868e6108f2eecaf` |
| Measured and latest L2 AVX2 object | `bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d` |
| New pure-AVX2 L1 archive | `eaf35fac2953161df251a7fc837bf1c94409435ff127ce90811fad5b4f6524bc` |
| New pure-L1 check executable | `45311b6fb79bfd2d122c15d276c8b82164a10c840b59885c84a8aa1fa131452d` |
| New native-L1 check executable | `7b4621ac59f8d5d9473b65bbd333b4d2e75da32850b029a1081d3ad4a7e7068e` |

`native-audit-v2/report.json` includes per-function counts, raw examples,
exact compile recipes and archive/member hashes. `pure-checks-v2/report.json`
includes all four pure object inventories, eight source pins, native argv and
results, full parity comparisons and CLI rejections. Full disassemblies,
archives, parity files, build recipes and resource logs are retained.
An independent shell/awk replay rechecked the archive/member bytes, original
manifest identity, nine full parity comparisons, raw EVEX counts and fresh
four-object disassembly without importing the collector or executing codecs.
It passed at 7,286,784 bytes under 256MiB, with all resource counters zero.

All work was serialized through `/tmp/leopard-gf8-authoritative.lock`.
The pure library build peaked at 152,805,376 bytes under 512MiB, serial `-j1`;
driver builds at 132,116,480 bytes under 512MiB. Native checks peaked at
207,626,240 bytes under 256MiB. These completed scopes reported all six
`memory.events` counters zero and swap zero; native children had a 30-second
CPU limit. No unrelated process affinities or kernel settings changed.

Two harness-only failures are preserved: the initial static audit rejected
ambiguous production/test-hook compile recipes (fixed by selecting the exact
production target), and the initial pure check collector failed on a `./`
manifest-name prefix (normalized before lookup). Both stopped before any new
codec executions; neither was a timing attempt or a codec failure.

## Next decision

Keep native Leopard1 as the original performance target. A fresh, pushed
preregistration can compare native L1, pure-AVX2 L1, and unchanged explicit
L2 AVX2 with controls, zero-sibling and immutable-artifact gates. No timing
was performed here, and no old attempt may be reused.

If the ISA-matched comparison still shows a gap, use that evidence to choose
a genuinely new AVX2 kernel/scheduling candidate. Do not repeat rejected
copy removal, smaller tiling, cache blocking or broad four-way fusion, or
infer dispatch costs from the existing native-L1 ratio. Candidate promotion
still requires focused Release and ASan/UBSan/LSan correctness plus the
preregistered 5% target / 2% controls-and-neighbors gates and independent
exact-L1 validation. The broader goal and production-delivery blocker remain
open. Claude review is explicitly waived by the user; review here is Codex
self-review and deterministic checks, not independent-model `CONVERGED`.
