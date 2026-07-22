# Expanded GF8 Leopard1/Leopard2 final map

`tools/leopard2_gf8_expanded_final_map.py` produces a deterministic, resumable
comparison of the exact Leopard `main` codec and the pure-AVX2 Leopard2 legacy
high profile.  It deliberately separates a broad parallel diagnostic screen
from publication-quality serial ABBA confirmation.

The baseline is fixed in the tool to commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` and executable SHA-256
`a43d7f43ff2e887ebcd47a1e94f806847a5d8b858a4e383e6c8d5e528a7dd910`.
The runner rejects another baseline.  It also rejects a dirty or mismatched
candidate source tree, a build directory configured from another tree, an
unexpected runtime backend/profile/field, digest disagreement, and any EVEX
instruction in either whole executable.  Both executables must contain AVX2
YMM instructions.  Before timing, the runner invokes the same candidate binary
once in schema-v5 source-attestation mode and requires the exact clean commit
and tree.  This benchmark-only instrumentation is inactive during the timed
schema-v3 decode-path calls unless `--attest-source` is explicitly requested.

## Prepare without benchmarking

Run the built-in tests and write the canonical matrix and host dry-run:

    python3 tools/leopard2_gf8_expanded_final_map.py self-test
    python3 tools/leopard2_gf8_expanded_final_map.py make-matrix \
      --output /tmp/leopard-gf8-expanded-map/matrix.json
    python3 tools/leopard2_gf8_expanded_final_map.py dry-run \
      --matrix /tmp/leopard-gf8-expanded-map/matrix.json \
      --output /tmp/leopard-gf8-expanded-map/dry-run.json

The matrix covers 64 B, 256 B, 1 KiB, 2 KiB, 4 KiB, 8 KiB, 16 KiB,
64 KiB, 1 MiB, and 16 MiB selected cells; T=1 through T=128; loss counts
0, 1, 2, 4, half, and maximum where valid; reuse counts 1, 8, and 64; and
batch counts 1 and 8.  It contains explicit R=32, K=33/34/35 at 16 KiB and
64 KiB to check that the promoted K=T+2 selector does not leak into the
rejected K=T+3 region.

## Stage 1: parallel diagnostic

Freeze and build the candidate first.  Substitute its full 40-character SHA,
source path, build path, and benchmark path below.  The runner uses every
allowed physical core pair by default and reserves each SMT sibling.  The
exclusive lock prevents every external build, test, diagnostic, or timing
phase using the canonical lock from overlapping the screen.  Only the runner's
own child jobs parallelize inside its lease.  Ratios from this stage remain
diagnostic only.

    python3 tools/leopard2_gf8_expanded_final_map.py diagnostic \
      --matrix /tmp/leopard-gf8-expanded-map/matrix.json \
      --baseline /tmp/leopard-gf8-final-map-9b1d439-20260722T030328Z/main-build/leopard_main_benchmark \
      --candidate /ABSOLUTE/CANDIDATE/BUILD/bench_leopard2 \
      --candidate-commit FULL_40_CHARACTER_CANDIDATE_SHA \
      --candidate-source /ABSOLUTE/FROZEN/CANDIDATE/SOURCE \
      --candidate-build /ABSOLUTE/CANDIDATE/BUILD \
      --run-dir /tmp/leopard-gf8-expanded-map/diagnostic \
      --memory-budget-gib 24 \
      --lock /tmp/leopard-gf8-authoritative.lock

Each completed cell is atomically stored under `cells/`; rerunning the exact
command resumes it.  Failed child stdout and stderr are retained under
`failures/`.  Changing the matrix, binary, source tree, CMake cache, CPU set,
or run options requires a new run directory.

## Stage 2: serial authoritative ABBA

The second stage selects every diagnostic near/slower cell (any primary ratio
at or below 1.10) plus all loss, reuse, batch, large-tiling, and
selector-isolation cells.  It holds the lock exclusively and runs three
baseline/candidate/candidate/baseline rounds on one physical core while
reserving its SMT sibling.  Explicit CPU numbers are recommended after
inspecting the dry-run CPU-pair list.

    python3 tools/leopard2_gf8_expanded_final_map.py abba \
      --matrix /tmp/leopard-gf8-expanded-map/matrix.json \
      --diagnostic-summary /tmp/leopard-gf8-expanded-map/diagnostic/summary.json \
      --baseline /tmp/leopard-gf8-final-map-9b1d439-20260722T030328Z/main-build/leopard_main_benchmark \
      --candidate /ABSOLUTE/CANDIDATE/BUILD/bench_leopard2 \
      --candidate-commit FULL_40_CHARACTER_CANDIDATE_SHA \
      --candidate-source /ABSOLUTE/FROZEN/CANDIDATE/SOURCE \
      --candidate-build /ABSOLUTE/CANDIDATE/BUILD \
      --run-dir /tmp/leopard-gf8-expanded-map/abba \
      --cpu ALLOWED_PHYSICAL_CPU --sibling ITS_RESERVED_SMT_SIBLING \
      --rounds 3 --near-ratio 1.10 \
      --lock /tmp/leopard-gf8-authoritative.lock

Finally, merge both summaries without changing their evidence status:

    python3 tools/leopard2_gf8_expanded_final_map.py merge \
      --diagnostic-summary /tmp/leopard-gf8-expanded-map/diagnostic/summary.json \
      --abba-summary /tmp/leopard-gf8-expanded-map/abba/summary.json \
      --output /tmp/leopard-gf8-expanded-map/evidence.json

The merged file labels the parallel screen non-authoritative and the serial
ABBA confidence intervals authoritative.  Ratios always use
`Leopard1 time / Leopard2 time`, so values above one favor Leopard2.
