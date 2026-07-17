# Leopard2 context-routed butterfly tier

Status: production implementation landed. The correctness and 16-round
non-regression results retained below are historical evidence for the pre-R1
candidate `ae2b8566ad08cae02869c3fcf41de046b0b652ed`; they are not final-source
evidence for the later integrated tree. R1 changed the backend table, sources,
tests, and build closure, so an integrated-source backend matrix and pinned
campaign must be collected before this tier is evidence-qualified again. The
first final-source AVX2 v6 campaign was retained as failed evidence and caused
the GF16 size policy correction documented below; a fresh campaign on the
corrected source is still required.
Native NEON, AVX-512, GFNI, and wider fusion schedules remain separate work.

The current v6 collector requires two distinct final-source campaigns, one
with `--backend ssse3` and one with `--backend avx2`. Each campaign binds the
requested and resolved backend reported by every benchmark process to the
same passing forced-backend matrix variant. The v6 geometry has 26 cells,
1,664 process invocations, and 52 encode/decode gates per backend. In addition
to the historical cells, it covers high- and low-profile GF8 at 1,025 bytes
and GF16 at 130 bytes, immediately above the production fusion thresholds.

## Production architecture

Leopard2's immutable private backend table now owns every fixed-multiply, XOR,
two-way butterfly, and four-way butterfly used by a transform. The completed
four-way boundary contains:

- four independent XOR pairs;
- forward and inverse GF8 fused-two-layer butterflies; and
- forward and inverse GF16 fused-two-layer butterflies over the legacy ALTMAP
  representation.

The three four-way multipliers describe pairs `(0,1)`, `(2,3)`, and
`(0,2)/(1,3)`. The GF8 and GF16 modulus sentinels retain the XOR edge while
suppressing its multiplication. Scalar, SSSE3, and AVX2 backends implement the
same equations. Startup known-answer tests compare each implementation with an
independent reference across multiplier boundaries, vector tails, unaligned
buffers, and guard bytes before an immutable table may be published.

Production field translation units are compiled with their legacy whole-file
SSSE3/AVX2 branches disabled. Every grouped fallback therefore receives the
calling context's `Ops` table. Explicit scalar and SSSE3 contexts can coexist
with AVX2 and AUTO contexts in one process without falling through to the
process-default ISA. Tracing tests separately count forward and inverse
four-way dispatch for each high/low GF8 and GF16 encode/decode callsite. They
also exercise generic decode, shared immutable codecs/plans, mixed contexts,
and grouped-XOR entries; the concurrency cases include both profiles.

The selected schedules are deliberately size-specific:

- GF8 inverse butterflies remain fused for every shard size.
- GF8 forward butterflies are fused through 1,024 bytes and use the qualified
  split two-way schedule above that size.
- GF16 forward and inverse butterflies fuse only exact 64- and 128-byte
  transform working sets. This includes a public 66-byte compact tail after
  staging rounds it to 128 bytes. The field callsite classifies every radix-
  four invocation before backend dispatch; a public 130-byte shard stages to
  192 bytes and therefore uses the two-layer split schedule, as do all larger
  working sets. Test-only counters independently prove total, fused, and split
  inverse/forward selections for both profiles at 64, 66, 128, 130, and 1,026
  public bytes.
- Scalar kernels retain the straightforward reference-equivalent operation
  order.

Grouped XOR also fixed a latent count-tail bug: after processing complete
groups of four shards, both threaded and non-threaded loops now advance by the
number actually consumed before processing one to three residual shards. A
regression test covers counts 1 through 13 in both modes.

No public API, field representation, coordinate map, profile identifier, or
wire byte changed.

## Rejected variants

The production cutoffs came from whole-codec measurements, not microkernel
instruction counts.

- Keeping GF8 forward fusion for large shards lost to the split schedule, so
  fusion stops at 1,024 bytes.
- Pre-broadcasting all three GF16 tables and fusing through 1,024 bytes caused
  roughly a ten-percent low-GF16 encode regression from table/register
  pressure. The code was reverted.
- A paired split helper sharing the cross-pair table and a two-XOR-coalesced
  variant both remained marginally outside the `-2%` confidence floor. A
  sequential-nibble rewrite showed no useful gain. All were reverted.
- A field-side shortcut was rejected by the tracing test because it bypassed
  the selected context's four-way entry.

Only exact 64-/128-byte GF16 fused tiles and the measured GF8 cutover remain in
production.

## Failed final-source AVX2 campaign and policy correction

The first v6 final-source AVX2 campaign compared candidate
`ca52c4e97472dd3af5a7544fbedd795b4ea724b5` with matched split control
`f3e8dfbc6cb1943c8fd8728801e38cf0a1ab3816`. It completed all 1,664
isolated invocations but correctly published a failed manifest. Five gates in
four GF16 cells fell below the fixed `-2%` lower-confidence floor:

| Cell and metric | Point estimate | One-sided 95% lower bound |
| --- | ---: | ---: |
| high GF16, 130 B encode | -4.231% | -4.515% |
| high GF16, 130 B decode | -3.888% | -4.079% |
| high GF16, 1 KiB decode | -1.591% | -2.122% |
| low GF16, 130 B encode | -2.732% | -3.078% |
| low GF16, 16 KiB decode | -0.145% | -2.052% |

The 130-byte result exposed the implementation error most directly: staging
rounded it to a 192-byte transform buffer, but the field callsite sent every
size through the fused four-way backend. The production correction now sends
only 64 and 128 transform bytes through `ff16_*_butterfly4`; every larger size
executes the historical two `ff16_*_butterfly2` layers with identical modulus
XOR semantics. This is a measured size policy, not a field or wire change.

The ignored failed artifacts are under
`.research/leopard2/context-xor4/final-ca52c4e-v6-avx2-clean1/`:

- `abba_manifest.json` SHA-256:
  `b03c1553ba02dbd22baa943e912e52e48713de92ca0cff774427d4c64c55fea9`
- `abba_raw.json` SHA-256:
  `22265cf4f27c21401edd0e42bfaf4c9eef6d7f5b104c4631b286dc8724b62f1e`

These results are retained rather than resampled or excluded. They do not
qualify either the old candidate or the corrected source.

## Historical pre-R1 correctness, safety, and portability evidence

The retained pre-R1 validation on the 2026-07-17 x86 host used every one of the
30 logical CPUs granted to that worker (`0-14,16-30`) for parallel builds and
tests. These counts and hashes describe the named historical candidate, not the
current integrated source:

- Release: 44/44 CTests passed.
- GCC 13.3 strict warnings (`-Wall -Wextra -Wpedantic -Werror`): 44/44 passed.
- Clang 18 strict warnings: 44/44 passed.
- ASan plus UBSan: all 42 applicable runtime tests passed. The separate static
  ISA classifier rejects sanitizer-inserted `sahf` in a baseline member; the
  ordinary Release archive passes that classifier, so this is recorded as a
  sanitizer-code-generation limitation rather than a codec failure.
- TSan initially hit this host's `unexpected memory mapping` startup failure.
  Direct runs under `setarch x86_64 -R` passed the 16-thread backend-ops test
  (1,024 executions), mixed-context test, and four-profile concurrent encoder
  test (528 executions) without a TSan report.
- That AUTO/scalar/SSSE3/AVX2 matrix passed every deterministic comparison,
  backend-failure subprocess, portable-ISA gate, and the default build's
  CUDA-optional test. Its source fingerprint is
  `8b63b1e8ea56d8a924f6087bda46fad43619fe07d9ec899baaf0783ad7ebb3cb`;
  `matrix.json` SHA-256 is
  `8053489efe9bd04c847e669166ee5ab88a449fc89676616d63e786e037b13bbe`.
- AArch64/SSE2NEON cross-compilation completed after those field changes.
  Existing ignored-OpenMP and unused-x86-parameter warnings remain. This is
  compile-only preservation evidence, not native ARM correctness or
  performance evidence.

The retained v5 runner revision passed path-independent replay plus 52
adversarial mutations. These include missing/reordered ABBA slots, raw-output
edits, compiler and build-graph substitutions, source-closure changes,
topology and reservation changes, confidence-summary edits, and a high-
variance fixture. The current v6 self-test additionally runs successful AVX2
and SSSE3 mock campaigns and rejects independent changes to raw requested or
resolved backend identity and to manifest requested or resolved backend
identity.

## Historical pre-R1 non-regression evidence

The retained campaign compares baseline commit
`4cbe17c4374739ae087dfae9568949a17a15b2f2` with candidate commit
`ae2b8566ad08cae02869c3fcf41de046b0b652ed`. It contains 22 cells, 16 paired
A-B-B-A rounds per cell, seven measured samples and three warmups per
invocation, reuse eight, and one thread: 1,408 isolated process invocations and
44 encode/decode gates. Fresh baseline and candidate builds used all 30 allowed
CPUs. Measurements were pinned to CPU 14 while its SMT sibling, CPU 30, was
reserved and idle.

This tier qualifies a correctness-driven context-routing refactor. Every target
and neighbor therefore uses a one-sided 95% lower-confidence non-regression
floor of `-2%`; it does not require every cell to gain five percent. Statistics
are paired log ratios with 15 degrees of freedom and include conservative
within-invocation uncertainty derived from each median, MAD, minimum, and
maximum.

Representative point estimates and one-sided 95% lower bounds are percentages
relative to the pre-tier baseline:

| Cell | Encode point | Encode lower | Decode point | Decode lower |
| --- | ---: | ---: | ---: | ---: |
| high GF8, K=240 R=16, 64 KiB, L=4 | +2.422% | +1.248% | +4.840% | +3.754% |
| low GF8, K=32 R=224, 64 KiB, L=16 | +2.437% | +1.923% | +3.916% | +3.175% |
| high GF16, K=1000 R=200, 16 KiB, L=8 | +1.548% | +1.132% | +1.024% | +0.566% |
| low GF16, K=128 R=1024, 16 KiB, L=16 | +0.359% | -0.664% | +2.547% | +0.830% |
| balanced GF8, K=128 R=128, 64 KiB | +2.120% | +1.756% | +7.630% | +6.630% |
| balanced GF16, K=512 R=512, 16 KiB | +0.161% | -0.039% | +2.130% | +1.545% |
| low GF16 neighbor, K=128 R=1024, 1 KiB | -1.131% | -1.482% | +2.462% | +2.286% |

The low-GF16 1-KiB encode cell is the weakest lower bound and still clears the
fixed floor. All 44 gates passed.

The ignored machine-readable historical artifacts are under
`.research/leopard2/context-xor4/abba-ae2b856-v5-16round/`:

- `abba_manifest.json` SHA-256:
  `f25eeafc9857052e3c333b5b095fbd18eb08e80b86a2ac2d46fee45aca8051c8`
- `abba_raw.json` SHA-256:
  `0ebd33990052d717aced183a68224888a2e1a5c9315b3e39798adbf671d02d85`

Path-independent replay of the manifest, raw bundle, embedded matrix, source
closure, build graph, and all statistics passed. Caller template binaries are
not interchangeable with the exact fresh-build executables because build paths
affect their byte hashes; the collector verified the exact executables before
publishing the passed manifest.

## Retained failed evidence

Failures were not discarded or converted into exceptions:

- The eight-round v4 campaign failed exactly one of 44 gates: low-GF16 1-KiB
  encode had a `-1.619%` point estimate and `-2.114%` lower bound. Its manifest
  and raw SHA-256 values are
  `c8c6fa06898d1c09811c3154da80d6316009ef77c783d07da8aedbae4ae59796`
  and
  `13fbc5456c0bf622da8b4925403f3d9c1bfcf96c679bdaa739a28d199c411dea`.
- A predetermined 16-round confirmation at commit `2b7a4ae` cleared all 44
  performance gates, but publication failed because GCC reports a different
  first-line driver name when the same resolved binary is invoked through
  `cc/c++`. The failed manifest/raw hashes are
  `3057db060a236fe2be646dda65eb51e744f3a738d275d202dcd0c55a0b6eda63`
  and
  `ad3358f8c8c7ee040506df13d929ac0f74ce53f29125969bff5e6b88e9c6adcd`.
  Replay after normalizing only that argv0-dependent label proved it was the
  sole failure. The executable basename, binary digest, package/version text,
  and remaining version output are still checked.

That historical campaign was the one clean replacement after that provenance
fix; there was no repeated sampling after its result.

These results remain useful causal history for the four-way implementation but
must not be cited as qualification of a later source tree. Final-source closure
requires fresh source-bound matrix and ABBA artifacts after R1 integration.

## Reproduction

Use every allowed CPU for compilation and correctness work, capped at 128:

    JOBS="$(nproc)"
    if [ "$JOBS" -gt 128 ]; then JOBS=128; fi
    cmake -S . -B build/release \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON \
      -DLEO2_BUILD_BENCHMARKS=ON \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build build/release -j "$JOBS"
    ctest --test-dir build/release -j "$JOBS" --output-on-failure

Run the forced-backend differential matrix and runner mutation tests:

    python3 tools/leopard2_backend_matrix.py run \
      --source . \
      --build-root build/backend-matrix \
      --result-dir build/backend-matrix-results \
      --variants auto,scalar,ssse3,avx2 \
      --jobs "$JOBS" --variant-workers 4 --timeout 900 --no-resume
    python3 experiments/leopard2/backend_butterfly/run_abba.py self-test

The final-source control must be constructed from the final candidate tree by
reverting only the production four-way runtime selections. It retains the
final R1, locator, Windows, test, matrix, and collector changes, so the
experiment has a single production-runtime treatment. The GF8 control retains
the documented two-layer replacements. For GF16, leave the new >128-byte
split schedule byte-identical and change only `UseFusedButterfly4()` to select
no sizes in the control. Candidate and control therefore differ only for the
still-enabled 64-/128-byte GF16 fused regime; 130-byte and larger cells execute
the same split schedule in both. Configure clean tests-disabled Unix Makefiles
Release builds for both trees with `LEO2_BACKEND_VARIANT=auto`,
`LEO2_BUILD_BENCHMARKS=ON`, `LEO2_BUILD_TESTS=OFF`, and
`CMAKE_EXPORT_COMPILE_COMMANDS=ON`.

The full campaign requires those clean baseline and candidate worktrees, a
fresh final-source four-variant matrix, a canonical no-newline reservation
file, and a physical core whose allowed SMT sibling is idle. Run the collector
twice, using distinct output directories and otherwise identical arguments:

    python3 experiments/leopard2/backend_butterfly/run_abba.py run \
      --backend avx2 \
      --baseline CONTROL_BUILD/bench_leopard2 \
      --candidate CANDIDATE_BUILD/bench_leopard2 \
      --baseline-commit CONTROL_SHA --candidate-commit CANDIDATE_SHA \
      --baseline-source-root CONTROL_TREE \
      --candidate-source-root CANDIDATE_TREE \
      --baseline-compile-commands CONTROL_BUILD/compile_commands.json \
      --candidate-compile-commands CANDIDATE_BUILD/compile_commands.json \
      --baseline-cmake-cache CONTROL_BUILD/CMakeCache.txt \
      --candidate-cmake-cache CANDIDATE_BUILD/CMakeCache.txt \
      --baseline-library CONTROL_BUILD/liblibleopard.a \
      --candidate-library CANDIDATE_BUILD/liblibleopard.a \
      --matrix FINAL_MATRIX/matrix.json --output ABBA_AVX2 \
      --cpu CPU --reserved-sibling SIBLING \
      --reservation-file RESERVATION_JSON --build-jobs "$JOBS" --timeout 120

Repeat with `--backend ssse3 --output ABBA_SSSE3`. Do not run other memory-
intensive jobs during either pinned phase. Each manifest can then be replayed
with `run_abba.py verify`; v6 checks the embedded raw evidence and matrix and
will not accept the historical v5 schema.

Replay the historical evidence with the historical v5 collector from the
named candidate commit, without substituting a different build-path binary.
The current v6 collector intentionally rejects a v5 manifest. Successful replay
authenticates that historical campaign only:

    git worktree add --detach ../leopard-abba-v5 \
      ae2b8566ad08cae02869c3fcf41de046b0b652ed
    python3 ../leopard-abba-v5/experiments/leopard2/backend_butterfly/run_abba.py verify \
      --manifest .research/leopard2/context-xor4/abba-ae2b856-v5-16round/abba_manifest.json \
      --raw-bundle .research/leopard2/context-xor4/abba-ae2b856-v5-16round/abba_raw.json \
      --matrix .research/leopard2/context-xor4/backend-matrix-results-final/matrix.json
