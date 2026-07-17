# Leopard2 selected-context R=1 XOR

## Result

Leopard2 now performs R=1 source pairing through the immutable backend table.
Encoding and one-loss recovery use one byte pass for each source pair on scalar,
SSSE3, AVX2, and the existing AArch64 execution path.  The change does not alter
the field, profile, coordinate map, parity bytes, scratch contract, or public
ABI.

The previous explicit-context path called the selected backend's ordinary XOR
twice for every source pair.  The process-default path had a special fused loop,
but that loop lived in the baseline translation unit.  Under the production ISA
isolation build it could therefore execute only baseline x86 instructions even
when an AVX2 context was selected.  A private `Ops::xor_memory_2to1` slot removes
both problems while preserving context tracing and runtime qualification.

## Contract and implementation

The internal primitive computes:

    destination[i] ^= source0[i] ^ source1[i]

The destination range must be disjoint from both source ranges.  Read-only
sources may alias exactly or partially because public shard inputs may alias.
A zero byte count accesses no memory.  Every implementation handles arbitrary
unaligned lengths without reading beyond the requested range:

- scalar uses unaligned-safe 64-bit `memcpy` chunks and byte tails;
- SSSE3 uses four 128-bit operations per main loop, then vector and byte tails;
- AVX2 uses four 256-bit operations per main loop, then vector and byte tails;
- AArch64 uses native 128-bit NEON loads, XORs, and stores before the portable
  tail.  The existing SSE2NEON build remains compile-tested.

Startup backend qualification independently checks zero and boundary sizes,
guards, exact source aliases, and output identity before publishing a table.
The dedicated test extends this through every byte count 0 through 257, selected
boundaries through 1,048,579 bytes, offsets, partially overlapping read-only
inputs, and concurrent immutable-table execution.

For legacy-high R=1 encoding, Leopard2 copies original zero to parity and then
pairs the remaining originals.  For direct one-loss recovery it copies parity
to the missing-original output, pairs all surviving originals, and handles one
unpaired survivor with the existing single-source XOR.  The latter pairing is
new; it avoids rereading and rewriting the restored shard once per survivor.

Correctness covers explicit GF8 and GF16 codecs.  GF16 cases include native
partial tiles, a padded-odd 33-byte payload in a 34-byte physical shard, and the
K=256 boundary whose 512-coordinate parent genuinely requires GF16.  Because
R=1 parity is bytewise XOR, field selection does not
change the arithmetic, but both codec identities are tested explicitly.

## Pinned crossover

The authoritative run compared baseline commit
`694881245bb07a58cb69b97fdffe5772d39f65e1` with runtime candidate commit
`800e0367baa06d8adb3304615b914a51698052aa`.  Both benchmarks were clean,
tests-disabled Release/Ninja/GCC 13.3 builds.  The executable hashes and other
recorded observations are frozen in
`experiments/leopard2/r1_xor/results/checkpoint.json`.

The SHA-256 values bind the exact executable and CMake-cache bytes used by the
run.  The source commit and dirty flag came from the source path recorded in
that cache at manifest-creation time.  Because the benchmark does not embed a
commit identity, those source fields are observational metadata, not a
cryptographic proof that a particular source tree produced the executable.  The
performance conclusion is therefore scoped to the exact executable hashes;
source reproducibility additionally depends on the documented clean-build
procedure.

The runner executed 51 cells: ten R=1 targets and seven neighbors on each of
scalar, SSSE3, and AVX2.  It spans K=3 through K=240, one byte through one MiB,
and batches 1, 8, and 64.  K=2 R=1 cells exercise the unaffected single-source
tail; R=2 cells exercise a neighboring non-fused transform path.  Every cell
used five baseline/candidate/candidate/baseline rounds, nine retained inner
samples, stable data seeds, and adaptive reuse targeting 256 MiB of offered
work.  Runs were pinned to logical CPU 14 while SMT sibling 30 stayed idle.

Improvement is `100 * (baseline_time / candidate_time - 1)`.  The reported 95%
interval is a deterministic percentile bootstrap over the five ABBA round
improvements.

| Backend | R=1 encode median range | Credible >=5% cells | R=1 decode median range | Credible >=5% cells | Worst neighbor median |
| --- | ---: | ---: | ---: | ---: | ---: |
| scalar | 6.54% to 31.45% | 10/10 | 4.70% to 31.37% | 8/10 | -0.61% |
| SSSE3 | 7.98% to 54.36% | 10/10 | 4.97% to 56.04% | 9/10 | -1.40% |
| AVX2 | 0.61% to 61.65% | 8/10 | 3.02% to 43.01% | 8/10 | -1.72% |

All affected backend/metric groups contain credible improvements above the 5%
promotion threshold.  No target or neighbor had an observed regression beyond
2%, and no credible regression beyond 2% was found.  The primitive and R=1
pairing are therefore enabled without a size threshold.  Tiny AVX2 encode cells
that do not clear 5% remain non-regressing, while the same change materially
improves their decode path; an extra dispatcher branch would add complexity
without protecting an observed losing region.

Raw artifacts remain ignored at
`.research/leopard2/r1-xor/pinned-800e036`: 51 job records, 1,020 benchmark
JSON files, and 2,040 stdout/stderr logs.  The compact checkpoint records hashes
for the manifest, analysis, runner, binaries, and the complete 3,113-file
artifact tree.

Exact timing command:

    tools/leopard2_r1_xor_crossover.py run \
      --baseline /home/catid/leopard-wt-r1-xor-baseline/build/r1-baseline/bench_leopard2 \
      --candidate /home/catid/leopard-wt-r1-xor/build/r1-benchmark/bench_leopard2 \
      --cpu 14 \
      --result-dir /home/catid/leopard-wt-r1-xor/.research/leopard2/r1-xor/pinned-800e036 \
      --backends scalar,ssse3,avx2 --grid compact \
      --abba-rounds 5 --iterations 9 --warmup 2 \
      --target-work-bytes 268435456 --maximum-reuse 100000 --timeout 180

Exact analysis command:

    tools/leopard2_r1_xor_crossover.py analyze \
      --result-dir /home/catid/leopard-wt-r1-xor/.research/leopard2/r1-xor/pinned-800e036

## Validation

The final implementation is covered by:

- Release CTest, including the backend startup KAT, explicit context tracing,
  public alias checks, fuzz smoke, compatibility, and portable ISA classifier;
- strict GCC 13.3 and Clang 18 builds with `-Wall -Wextra -Wpedantic -Werror`;
- Clang ASan+UBSan CTest (the assembly classifier is excluded because sanitizer
  instrumentation intentionally introduces instructions outside its baseline
  policy); focused Clang TSan context, R=1, and concurrency tests;
- scalar, SSSE3, and AVX2 production-table execution in one binary;
- AArch64/SSE2NEON cross-compilation at submodule commit
  `cad518a93b326f0f644b7972d488d04eaa2b0475`, with native NEON `eor v*.16b`
  instructions verified in the scalar-table object used by the existing ARM
  execution path.

The benchmark candidate predates test and documentation commits.  The only
later runtime-source edit narrows the inline NEON guard to AArch64, where ASIMD
is mandatory; this prevents optional-NEON ARMv7 from executing NEON through the
scalar backend before feature qualification.  That preprocessor-only ARM change
does not alter the benchmarked x86 executable path, so the checkpoint remains
exact evidence for scalar, SSSE3, and AVX2 on the recorded x86 machine.  The
AArch64 path is correctness- and cross-build-tested but has no performance claim
from this run.
