# 41 — Sparse-high production-reference attribution

## Result

The sign reversal between the sparse-high forced-path discovery and the
production-AUTO qualification is explained by a reference-path mismatch, not
by table preparation. The discovery benchmark's `force_transform` arm is a
test-only exact sparse-schedule transform. An ordinary production archive
cannot select that hook path. For the measured legacy-high GF8 cells, no
separately qualified production sparse schedule applies, so the archive builds
no exact schedule and executes the mature prefix transform. The production
qualification correctly compared AUTO direct against that ordinary prefix
path.

This distinction matters most at padded side two (`R=2`). The AVX2 T=2 prefix
encoder produces both transform rows together. Its K=2/K=3 generated terminal
uses one or two products shared across both outputs; larger K values begin with
a two-output butterfly and accumulate subsequent two-source blocks into both
outputs. The sparse-Q1 candidate instead uses the generic row-major direct
executor: for the one requested row it walks all K sources and performs one
copy/multiply followed by K-1 XOR or multiply-add operations. Skipping the
unrequested row therefore does not compensate for the mature T=2 transform's
shared arithmetic and memory traffic at the qualified production sizes.

The evidence does **not** authorize a narrower production selector after the
fact. Sparse-high AUTO remains compiled-default `OFF`. The next candidate is a
mechanism-derived boundary that excludes padded side two and retains the
predeclared R=4/8/16 region only for a fresh, independently versioned no-hook
production campaign.

## Control-flow proof

The source establishes four independent facts:

1. `TestForcesSparseEncodeSchedule` returns true for
   `LEO2_TEST_ENCODE_FORCE_TRANSFORM` only when `LEO2_ENABLE_TEST_HOOKS` is
   compiled. In an ordinary archive it always returns false.
2. At padded side two or greater, `EncodeLayout` reserves bounded
   exact-schedule storage for a hook build, and `PrepareSparseEncodePlans`
   compiles that schedule when the force-transform hook is active. If neither
   that hook nor the separately qualified sparse LOW-GF16 production case
   applies, `PrepareSparseEncodePlans` returns an empty view before schedule
   compilation.
3. `bench_leopard2_direct_encode` maps its transform mode to
   `LEO2_TEST_ENCODE_FORCE_TRANSFORM` and links the test-hooks archive. In
   contrast, `bench_leopard2_high_sparse_auto` rejects
   `LEO2_ENABLE_TEST_HOOKS` at compile time and requires an explicit
   `LEO2_HIGH_SPARSE_AUTO_LIBRARY_TEST_HOOKS=0` marker.
4. `ExecuteDirectEncodeRows` is the production sparse-Q1 candidate.
   `LeopardFF8.cpp::ReedSolomonEncode` treats the empty production plan as a
   dense schedule, admitting `ff8_high_encode_small`; a nonempty forced sparse
   plan bypasses that fast path and reaches `ExecuteSparseForwardPlan`. The
   AVX2 implementation behind `ff8_high_encode_small` routes side two through
   `AVX2FF8HighEncodeT2Direct<2>` or `AVX2FF8HighEncodeT2Direct<3>`, or through
   the two-output butterfly path. Both the AVX2 backend and `LeopardFF8.cpp`
   source blobs are byte-identical at the discovery source
   `c462f0ddef75979749c8cdd8e7199326ff1be892` and production source
   `b455108e61017c711e238ccf159d45da50e77ca2` (Git blobs
   `9b592eb30dba82f473d1c0a9b8d9e6b270dfa05f` and
   `9f95f143feaaa6d67f6c860e2fdac7a8e47619b6`). The
   `ExecuteDirectEncodeRows` body is also unchanged. The source therefore
   attributes the reversal to selection of a different reference mechanism,
   not to a changed direct executor or transform implementation; the campaigns
   did not perform a same-binary three-arm timing isolation.

The relevant implementation sites are `IsHighSparseDirectEncodeShape`,
`TestForcesSparseEncodeSchedule`, `TestBuildReservesSparseEncodeSchedule`,
`EncodeLayout`, `PrepareSparseEncodePlans`, `ExecuteTransformEncodePass`, and
`ExecuteDirectEncodeRows` in `leopard2.cpp`; `AVX2FF8HighEncodeT2Vector` and
the `side == 2` branch in `Leopard2BackendAVX2.cpp`; `ReedSolomonEncode` and
its `dense_schedule`/`ExecuteSparseForwardPlan` branches in `LeopardFF8.cpp`;
`TestMode` in `bench/leopard2/direct_encode_benchmark.cpp`; and the no-hook
compile-time guards in `bench/leopard2/high_sparse_auto_benchmark.cpp`.

## Quantitative attribution

The table below compares only the paired gain *within* each campaign. Absolute
microseconds are not compared or subtracted across campaigns: discovery ran on
`ripper.lan` and production qualification ran on `foureyes.lan`, using
different frozen binaries and CPU pairs even though both hosts report the same
Threadripper PRO 9985WX model.

| Shape/API pairing | Forced exact-schedule contrast | Production AUTO vs ordinary prefix | Production table point estimate (95% interval) |
| --- | ---: | ---: | ---: |
| `K=3,R=2`, 4096 B: diagnostic batch 1 / production one-shot | +55.81% (42.39% to 78.37%) | -27.06% (-31.88% to -21.90%) | +2.65% (-9.39% to +16.30%) |
| `K=16,R=2`, 4096 B: diagnostic batch 1 / production one-shot | +45.62% (40.10% to 52.37%) | -22.92% (-26.03% to -19.67%) | +1.32% (-5.12% to +8.20%) |
| `K=16,R=2`, 4096 B: diagnostic batch 1 / production batch 1 | +45.62% (40.10% to 52.37%) | -24.65% (-30.37% to -18.44%) | -0.70% (-1.43% to +0.03%) |
| batch 16 `K=16,R=2`, 4096 B | +40.74% (32.03% to 52.29%) | -10.94% (-12.38% to -9.48%) | -0.06% (-1.17% to +1.07%) |
| `K=16,R=4`, 4096 B: diagnostic batch 1 / production one-shot | +105.73% (104.43% to 107.26%) | +21.74% (14.73% to 29.19%) | -0.29% (-2.77% to +2.27%) |

All twelve production one-shot failures are R=2, including every 4096-byte
K in `{2,3,4,8,12,16}` and the six additional K=16/R=2 byte boundaries. All
six public batch/binding failures are also R=2. Conversely, all 56 old
candidate cells outside R=2 passed their original gate; the weakest R=4 route
lower bound was 14.73 percent. That separation supports the T=2 mechanism as a
candidate-design boundary, but those already-observed cells remain exploratory
for the redesign and cannot be reused as promotion evidence.

Moving generator-table setup cannot explain the loss: the production route
contrast holds tables enabled on both arms. The near-zero table contrasts in
the decisive K=16/R=2 public-API cells support that conclusion while the route
loss persists from batch 1 through batch 16. Binding creation already prepares
and validates reusable state before its timed execute path, yet the same R=2
loss remains.

## Evidence lineage and disposition

The forced-path checkpoint is
`experiments/leopard2/direct_encode/results/sparse_high_avx2_checkpoint_20260831.json`
(SHA-256
`47cb26cdbcd41fc29765c5bdb049c83a54f6e3307e838e74fb2926aac6ba616e`).
Its retained analysis SHA-256 is
`3f1c60ebc9f9bc780daccef8d22dcc558e7d150db749f22a26fff2cd0bb513fc`,
and its raw result root is
`ripper.lan:/home/catid/leopard-sparse-high-dddd577/results/leopard2/direct-encode-crossover/sparse-high-avx2-ripper-c462f0d-quiet-cpu14`.

The production checkpoint is
`experiments/leopard2/direct_encode/results/production_auto_avx2_checkpoint_20260831.json`
(SHA-256
`1d79118faccad907c80d578d07e42d382b4ae0514480aa1f96173099847fc666`).
Its retained analysis SHA-256 is
`a6da9452f530e3ac3af183a268822094c0715da84b06b19a01ad26f8eff6a853`,
and its raw result root is
`foureyes.lan:/tmp/leopard-production-auto-qualification-b455108-foureyes-v2`.

Report 39 remains a valid forced exact-schedule discovery result, and report 40
remains the authoritative negative result for the original 36-tuple production
policy. No measurements are pooled, relabeled, or deleted. The follow-up work
is tracked under Bead `leopard-79h.42.10`: first change the mechanism boundary,
then version and run a fresh no-hook campaign with R=2 structural negative
controls before making any promotion decision.
