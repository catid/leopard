# Leopard2 context-routed butterfly tier

Status: production implementation and the conservative GF16 size policy are
landed at `04a8ba35fc9b1eb068cd748358eb04b79bc5e63b`. Final-source Release,
GCC/Clang strict, ASan/UBSan, focused TSan, portable-ISA, and forced-backend
correctness gates pass, with the TSan test-link limitation documented below.
The retained timing evidence decisively supports every measured SSSE3/AVX2
GF16 treatment but does not pass the predeclared whole-campaign AVX2
statistical gate: two large-shard decode controls whose executed code is
identical crossed the
one-sided `-2%` confidence floor. That failure is retained and reported below;
it was not resampled or relabeled as a passing campaign. The production policy
remains because exact linked-code and data-layout comparison rules out a
size-policy change in those failing controls, while the actually changed
64/128-byte paths have large positive lower bounds.

The older correctness and 16-round non-regression results retained below are
historical evidence for the pre-R1 candidate
`ae2b8566ad08cae02869c3fcf41de046b0b652ed`; they are not substituted for the
final-source evidence.
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
- GF16 forward and inverse butterflies dispatch exact 64-byte transform
  working sets through the fused operation on every backend. AVX2 also fuses
  exact 128-byte working sets. SSSE3 keeps the two-layer split schedule at
  128 bytes because
  final-source whole-codec evidence found regressions for public 66-byte compact
  tails after staging rounded them to 128 bytes. A public 130-byte shard stages
  to 192 bytes and therefore uses the split schedule on every backend, as do all
  larger working sets. Test-only counters independently prove total, fused, and
  split inverse/forward selections for each traced x86 context and both profiles
  at 64, 66, 128, 130, and 1,026 public bytes.
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

Only exact 64-byte GF16 fused tiles on every backend, exact 128-byte tiles on
AVX2, and the measured GF8 cutover remain in production. The 64-byte route is
covered for correctness across the available x86 backends; the final-source
performance campaigns below cover explicit SSSE3 and AVX2 on this x86 host.

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
64 transform bytes on every backend and 128 bytes only on AVX2 through
`ff16_*_butterfly4`; every larger size executes the historical two
`ff16_*_butterfly2` layers with identical modulus XOR semantics. This is a
measured size/backend policy, not a field or wire change.

The ignored failed artifacts are under
`.research/leopard2/context-xor4/final-ca52c4e-v6-avx2-clean1/`:

- `abba_manifest.json` SHA-256:
  `b03c1553ba02dbd22baa943e912e52e48713de92ca0cff774427d4c64c55fea9`
- `abba_raw.json` SHA-256:
  `22265cf4f27c21401edd0e42bfaf4c9eef6d7f5b104c4631b286dc8724b62f1e`

These results are retained rather than resampled or excluded. They do not
qualify either the old candidate or the corrected source.

## Corrected final-source campaigns

The corrected size cutoff was measured from commit `5433f88` against a matched
control that disabled the remaining GF16 fusion and restored the historical
GF8 split choices. Both explicit-backend campaigns completed all 1,664 planned
ABBA invocations pinned to CPU 14 while the runner reserved SMT sibling 30.
Neither failed campaign was resampled.

The AVX2 candidate showed strong GF16 gains at the intended sizes: roughly
29 percent for high-rate 64-byte encode/decode, 20/5 percent for low-rate
64-byte encode/decode, and 10/8 percent for a high-rate public 66-byte tail.
Its manifest nevertheless failed the declared whole-campaign gate because two
large GF16 decode cells, whose candidate and control both selected the split
schedule, had near-zero point estimates but one-sided confidence bounds just
below the fixed -2 percent neighboring floor. The retained manifest SHA-256 is
`8124004310308ad9d232ff0c56f66a4b87366497c1a382975ac7a17d83c67532`;
the raw bundle SHA-256 is
`0b4c2a6fe17b1466a90a2812220b3720af45bcff37e44e152ddb5b5a252f9e23`.

The SSSE3 campaign exposed a real size-specific loss. Its exact 64-byte GF16
cells improved, but the staged 128-byte treatment regressed high-rate decode by
2.40 percent and low-rate encode by 5.92 percent for public 66-byte shards. The
retained manifest SHA-256 is
`86ac9698c9085b37db540796252dc90e909c1e068f0a24dac20d85bcf0353e29`;
the raw bundle SHA-256 is
`cc243b91e323405917cab0d969ac90111ce5850565915f26bd0171a6c2445b42`.
Production therefore keeps SSSE3 fusion at exactly 64 transform bytes and uses
the split schedule at 128 bytes. AVX2 retains both qualified sizes. The final
source-bound matrix and sensitivity evidence are recorded next.

## Final exact-source qualification and disposition

Commit `04a8ba35fc9b1eb068cd748358eb04b79bc5e63b` passed the final CPU
correctness and safety gates available on this host:

- Release: 50/50 CTests passed.
- GCC 13 strict warnings: 48/48 CTests passed.
- Clang 18 strict warnings: all 58 compile commands contained
  `-Wall -Wextra -Wpedantic -Werror`, the build emitted no warning or error,
  and 50/50 CTests passed.
- Clang 18 ASan plus UBSan: all 47 applicable CTests passed.
- Clang 18 TSan with OpenMP disabled passed seven focused executables without
  a report: backend operations, mixed contexts, encode concurrency, R=1, high
  and low shared-plan decode, and shared-context public batches. The plan-
  schedule test's intentional strong global `new`/`delete` overrides conflict
  with compiler-rt TSan at link time; that exact test passes in Release, while
  its concurrent shared-plan behavior is covered by the passing high/low and
  context TSan targets.
- The portable-ISA archive classifier passed. The AArch64/SSE2NEON compile-
  only gate also passed at exact gitlink `cad518a`, including native vector XOR
  instructions; this is not native-ARM performance evidence.
- The fresh AUTO/scalar/SSSE3/AVX2 matrix passed all four variants, 145/145
  semantic gates, and 156/156 normalized cross-backend comparisons. Its
  `matrix.json` SHA-256 is
  `21d5a32d492759fd68026c4ed0e971a2c9bf50d474d7ccd60a68a55adac4cc38`;
  the 58-file source fingerprint is
  `878afd91fb09d99eebfb6e7d5673ee218ba7479ab9d6622469a9ca01d48fad88`.

The exact-source policy-matched SSSE3 campaign passed all 26 cells and 1,664
ABBA entries. Exact 64-byte GF16 high-rate encode/decode improved by
12.708/10.699 percent (one-sided 95% lower bounds 12.332/10.390 percent), and
low-rate encode/decode improved by 4.052/1.866 percent (lower bounds
3.671/1.529 percent). Public 66-byte shards, which stage to 128 transform
bytes, remained on the split SSSE3 policy and cleared the neighbor floor. The
manifest/raw SHA-256 values are
`45623d45a183bb8e5bacbd2a9c3ddce318645f76111b4443de77395917383d23`
and
`465ad82e4affc346e49e85eeef263d78eb8fcbc050415304c8aeb340b3dee812`.

The corresponding exact-production AVX2 campaign retained strong intended
gains but failed two large LOW_V1 GF16 decode confidence bounds: the loss-16
cell had a -0.151 percent center and -2.226 percent lower bound; the one-loss
cell had a -0.666 percent center and -2.608 percent lower bound. Its retained
manifest/raw SHA-256 values are
`5c326327d335b464a071f9546620225e5de1c71438ada79aef81f1593e690219`
and
`76be32c6455b47d684352b8bc2d486c3c8b8a1e34ebe20580c802f4831a8a899`.
To separate a real treatment effect from layout or eligibility-check changes,
a final predeclared sensitivity pair was built from the production tree. The
false commit `4817f471424e3dfed177b98d50421f46446a09c3` and true commit
`ac87bc208625398a6de85c92c1ab06b2ca8d0c7a` have identical linked code,
symbols, relocations, allocated runtime section layout, and normalized build
recipes. Their non-loadable debug sections reflect different worktree-path
lengths and are not claimed identical. The binaries differ by one initialized
volatile gate byte at executable offset `0x4b030`; the true arm selects the
production schedule. The true arm's fresh four-backend matrix passed with zero
mismatch (`matrix.json` SHA-256
`c95d574c82c620652586289103652de04ce971212e670682ebafa6c802bac158`).

The AVX2 sensitivity campaign completed its one allowed set of 1,664 pinned
invocations. It did **not** pass the declared whole-campaign gate:

| Cell and metric | Point estimate | One-sided 95% lower bound |
| --- | ---: | ---: |
| low GF16, 16 KiB, 16 losses, decode | -0.3951% | -2.3600% |
| low GF16, 16 KiB, one loss, decode | -0.3537% | -2.4995% |

No other cell failed. Both failures execute the same large-shard split code in
the two layout-identical binaries. The paths actually changed by the gate were
strongly positive:

| GF16 treatment | Encode point/lower | Decode point/lower |
| --- | ---: | ---: |
| high, 64 B | +27.8547% / +27.0404% | +28.0752% / +27.4430% |
| high, public 66 B -> 128 B | +8.8682% / +8.3740% | +8.4342% / +8.1551% |
| low, 64 B | +19.6951% / +18.8706% | +4.8058% / +4.4769% |
| low, public 66 B -> 128 B | +5.0400% / +4.1726% | +1.4794% / +0.8448% |

The failed sensitivity manifest/raw SHA-256 values are
`ff443934cd38f3d24f21906503977c1da1255e5ce9bf038dd48e9d982c6071ab`
and
`09577453d6f36f8809763c23bb255b2f4fbdb23b045ef5fff78c40efb1df1830`;
the stable raw-evidence digest is
`540d2dc65826c7032c3a74d9ff1e1a6b97dd96463d546ca9b3667171ff4da2b8`.
An independent replay authenticated all provenance and all 1,664 raw records,
then reproduced the two policy failures. The current v6 `verify` command
rejects a deliberately failed status before full replay, and its CPU
reservation is advisory rather than a measurement of sibling activity. The
current collector additionally holds the same pair-wide filesystem and Linux
abstract-socket lease used by exact-main, low-copy, and Jerasure collectors, so
distinct reservation files cannot cause overlapping measurements on one
physical pair; those
tooling limitations are tracked separately and are not a reason to resample
this pair.

Final disposition: retain exact 64-byte GF16 fusion on all production backends,
retain exact 128-byte fusion only on AVX2, and keep every larger GF16 transform
on the split schedule. This is not reported as a passed AVX2 whole-campaign.
The decision combines exact-path identity for the two failing controls,
decisive positive lower bounds for every measured x86-SIMD changed path, the
separate exact-production matrix, and the passed SSSE3 campaign. Scalar
correctness is qualified, but no scalar performance promotion is claimed. No
wire or API identity depends on this kernel choice.

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

The ABBA collector bounds benchmark stdout, stderr, and elapsed time. Linux
collection uses temporary child-subreaper ownership, procfs ancestry snapshots,
pidfd signaling, and bounded adopted-child reaping, so a benchmark cannot escape
cleanup with `setsid()` and a double fork. The self-test exercises that exact
case and proves the daemon PID and delayed marker are gone. Authoritative
collection fails closed before spawning when these Linux facilities are not
available; structural replay remains portable. A separately implemented
emergency procfs/pidfd/prctl path runs after every post-spawn exception, so a
fault in normal attachment, signaling, or teardown cannot strand the process
tree or leave the runner in a changed child-subreaper state. Run
`experiments/leopard2/test_process_containment.py` to exercise those faults
against the butterfly, low-copy, and exact-main runners together.

The final-source control must be constructed from the final candidate tree by
reverting only the production four-way runtime selections. It retains the
final R1, locator, Windows, test, matrix, and collector changes, so the
experiment has a single production-runtime treatment. The GF8 control retains
the documented two-layer replacements. For GF16, leave the new >128-byte
split schedule byte-identical and change only `UseFusedButterfly4()` to select
no sizes in the control. Candidate and control therefore differ only for the
still-enabled backend-specific GF16 fused regime: 64/128 bytes for AVX2 and
64 bytes for SSSE3. The 130-byte and larger cells execute the same split
schedule in both. Configure clean tests-disabled Unix Makefiles Release builds
for both trees with `LEO2_BACKEND_VARIANT=auto`,
`LEO2_BUILD_BENCHMARKS=ON`, `LEO2_BUILD_TESTS=OFF`, and
`CMAKE_EXPORT_COMPILE_COMMANDS=ON`.

The full campaign requires those clean baseline and candidate worktrees, a
fresh final-source four-variant matrix, a canonical no-newline reservation
file, and a physical core whose allowed SMT sibling is idle. Run the collector
twice, using distinct output directories and otherwise identical arguments:

The current collector first holds a nonblocking exclusive flock on the owned
`/run/user/UID` directory inode and retains it until the reservation handle is
closed.  C7 and the exact-main collector use the same stable anchor.  Therefore
replacement of the reservation inode or its containing directory cannot let a
current peer acquire a disjoint file lock while a campaign is active.  The v6
evidence schema and historical portable replay remain unchanged because the
anchor is execution authority, not retained benchmark data.  This stable layer
conservatively serializes all current Leopard2 evidence campaigns for the UID,
including campaigns on disjoint pairs.

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
with `run_abba.py verify` when its campaign status is `passed`; v6 checks the
embedded raw evidence and matrix and will not accept the historical v5 schema.
The current command deliberately stops at `campaign status` for retained
failed manifests instead of replaying them to the expected statistical
failure. That limitation, and fail-closed measurement of reserved-sibling
activity, are tracked by the evidence-runner follow-up. Do not rerun the final
layout-matched AVX2 pair to work around this verifier behavior.

For AVX2, candidate and control differ at the 64- and 128-byte GF16 schedule
choices. For SSSE3, they differ only at 64 bytes; both select the split schedule
at 128 bytes. Scalar and native-ARM performance remain separate evidence gates.

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
