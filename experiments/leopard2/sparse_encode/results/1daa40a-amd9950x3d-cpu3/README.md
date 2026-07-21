# LOW GF16 AVX2 sparse-output decision gate

This directory retains the accepted pinned kernel evidence for deciding whether
the sparse-output encoder needs an immutable encode-plan API to realize its
measured gain.  It does not promote a production dispatcher.  Production
`AUTO` remains unchanged until neighboring K/R cells and the public encoder
pass a separate end-to-end gate.

## Result

Ordinary call-local schedule construction is fast enough in the measured
LOW-GF16 AVX2 target.  Both sparse masks clear the 5% lower-confidence-bound
rule at 1 KiB and 64 KiB even though schedule setup is included in every call.
An immutable encode plan is therefore not a prerequisite for production
plumbing in this bounded region.  The 64-byte rows reject a blanket rule and
the dense-prefix controls reject applying the exact scheduler to dense output.

The primary values below are geometric call-local speedups derived from three
independent paired ABBA rounds.  Intervals are two-sided 95% Student-t
intervals with two degrees of freedom.  Prepared gain is a descriptive median
that excludes schedule setup.

| Shard bytes | Mask | Role | Call-local gain | 95% interval | Prepared gain |
|---:|---|---|---:|---:|---:|
| 64 | edge sparse | candidate | -29.160% | [-31.060%, -27.208%] | +6.161% |
| 64 | scattered sparse | candidate | -34.040% | [-34.621%, -33.453%] | -0.400% |
| 64 | dense prefix | control | -37.577% | [-38.585%, -36.552%] | -12.445% |
| 1,024 | edge sparse | candidate | +45.415% | [+44.977%, +45.855%] | +66.041% |
| 1,024 | scattered sparse | candidate | +33.806% | [+33.697%, +33.915%] | +52.798% |
| 1,024 | dense prefix | control | -14.364% | [-16.259%, -12.426%] | -8.293% |
| 65,536 | edge sparse | candidate | +61.198% | [+61.031%, +61.365%] | +61.762% |
| 65,536 | scattered sparse | candidate | +50.987% | [+50.867%, +51.107%] | +51.354% |
| 65,536 | dense prefix | control | -6.447% | [-6.484%, -6.410%] | -6.289% |

The aggregate preliminary gate is deliberately false: it sees the expected
64-byte and dense-control regressions.  This is not a contradiction.  A later
dispatcher must exclude those cells.  The retained summary also keeps
`production_promotion_sufficient=false` because this matrix does not cover
neighboring K/R values or ordinary public-API validation and staging costs.

## Identity and isolation

- Source: clean `1daa40a343bc14b648fba6b9c3225b62b358c56c`.
- Compiler: GCC 13.3.0, Release, C++11, test hooks disabled.
- Host: AMD Ryzen 9 9950X3D; allowed CPUs `0-14,16-30`.
- Benchmark CPU: 3; reserved SMT sibling: 19.
- Isolation: accepted; sibling non-idle delta exactly 0 jiffies over 1,124
  total jiffies.  The canonical abstract-socket plus `flock` pair lease was
  held for CPUs 3 and 19.
- Jobs: 9/9 passed; matrix `authoritative=true`.
- Iterations: 64; warmups: 4; independent ABBA rounds: 3; setup iterations:
  64; one worker; reuse models: 1, 8, and 64.
- Source fingerprint:
  `66e8d014162b6530178ba902eb14464dcfc34347e35ddcf5babd92d403d2e3e9`.
- Benchmark executable SHA-256:
  `369fe03c22501cedd6f8bb37c89eebcd94675ec9eafb0f34f9049c92866d09a6`.
- Link sidecar SHA-256:
  `eb3ba9ab6a4bf919ad0a08f036c2cb864fbebda95aa95ca55bc50897f7eead79`.
- Cell-manifest SHA-256:
  `50b6e2d196f59c4768777dc045f2aca65997ac9f90a6f51d8d3a6d7fb5cebc7a`.

Retained artifact hashes:

- `manifest.json`:
  `96f9d556329d6e61494339487aa838f5ed0fa5d44b45c0a332affe74126c6044`.
- `matrix.json`:
  `e3b94bab00875fc89c6ee4fa4106a41988379258dfdceeb61fee72f0393e3edf`.
- `summary.json`:
  `ee8c65484624173cbb34744b742312b36a57eb08f31e0b3028cfd7aa4bf78b3f`.

The manifest and matrix preserve the exact absolute paths observed during the
run.  Those paths are provenance, not portable build instructions.  The matrix
contains every raw benchmark result and ABBA observation; the separately
retained per-job stdout files would duplicate those embedded JSON objects.

## Reproduction

Start from the exact measured commit and configure after detaching so the
benchmark embeds the clean source identity:

    git switch --detach 1daa40a343bc14b648fba6b9c3225b62b358c56c
    cmake -S . -B build/sparse-pinned -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=OFF \
      -DLEO2_BUILD_BENCHMARKS=ON
    cmake --build build/sparse-pinned \
      --target bench_leopard2_sparse_encode -j 128
    python3 tools/leopard2_sparse_encode_benchmark_json_test.py \
      build/sparse-pinned/bench_leopard2_sparse_encode
    python3 tools/leopard2_sparse_encode_crossover.py self-test

Create a fresh isolation attestation naming an allowed physical CPU and its
complete SMT sibling set.  After every competing worker explicitly confirms
quiet, run:

    OMP_DYNAMIC=FALSE OMP_NUM_THREADS=1 \
    python3 tools/leopard2_sparse_encode_crossover.py pinned \
      --source . \
      --executable build/sparse-pinned/bench_leopard2_sparse_encode \
      --result-dir results/leopard2/sparse-encode-pinned \
      --cell-manifest \
        experiments/leopard2/sparse_encode/low_gf16_avx2_call_local_cells.json \
      --cpu CPU --workers 1 --no-resume \
      --isolation-attestation /absolute/path/isolation.json

Revalidate retained evidence while the exact source and executable still
exist:

    python3 tools/leopard2_sparse_encode_crossover.py analyze \
      --result-dir results/leopard2/sparse-encode-pinned
