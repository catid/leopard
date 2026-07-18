# GF16 AVX2 neighboring-regime evidence

This artifact closes the same-source neighboring-regime screen for production
candidate `7921998` against its direct parent `a725bf4`.  The only codec-source
difference is the GF16 AVX2 range-table hoist.  The integrated equivalent is
`201f6b9`; later commit `39cc427` adds tests but does not change the measured
kernel.

## Decision

The screen passed its predeclared non-regression rule.  None of the 13 encode
point estimates and none of the 13 decode point estimates exceeded a 2%
candidate slowdown.  The largest observed candidate slowdown was 0.60% for
encode and 0.065% for decode.  There was no statistically credible greater-than-
2% regression: no 95% interval had its upper bound below the 0.980392
control/candidate floor.

This is a neighboring-regime non-regression result, not a claim that every cell
is faster.  Nine metric intervals still cross the 2% floor because three
independent ABBA rounds provide limited precision in noisy cells.  Their point
estimates remain within the predeclared gate and their full intervals are
retained in `summary.json`; further rounds would be appropriate before making a
positive speedup claim in those cells.

All 156 processes used explicit GF16, explicit AVX2, one context thread, batch
one, fixed per-cell seeds, retained samples, and maximum valid original loss.
Every process passed the public round trip.  Original, transmitted-parity, and
recovered-original digests matched exactly between control and candidate in
each cell.  AUTO dispatch and Leopard-main are outside this same-source screen.

## Results

Speedup is control execution time divided by candidate execution time.  Values
above one favor the candidate.  Intervals are two-sided 95% Student-t intervals
over one paired-log contrast from each of three independent ABBA rounds.

| Cell | Encode speedup [95% CI] | Decode speedup [95% CI] |
|---|---:|---:|
| high 1000/200, 4 KiB | 1.0016 [0.9444, 1.0622] | 1.0053 [0.9827, 1.0283] |
| high 1000/200, 64 KiB | 1.0066 [0.9844, 1.0292] | 1.0031 [0.9988, 1.0074] |
| high 4096/512, 4 KiB | 1.0658 [1.0270, 1.1061] | 1.0073 [0.9707, 1.0451] |
| high 4096/512, 4098 B | 1.0543 [1.0440, 1.0647] | 0.9994 [0.9960, 1.0027] |
| low 512/1536, 4 KiB | 1.0374 [1.0246, 1.0504] | 1.0344 [1.0205, 1.0485] |
| low 1600/2496, 4 KiB | 1.0134 [0.9593, 1.0704] | 1.0122 [0.9660, 1.0605] |
| low 2032/2064, 4 KiB | 1.0035 [0.9782, 1.0294] | 1.0026 [0.9716, 1.0347] |
| high 1000/200, 64 B | 0.9940 [0.9831, 1.0051] | 1.0024 [0.9890, 1.0159] |
| high 1000/200, 66 B | 1.0076 [0.9917, 1.0237] | 1.0068 [1.0016, 1.0120] |
| high 1000/200, 128 B | 0.9943 [0.9704, 1.0188] | 1.0009 [0.9909, 1.0109] |
| high 1000/200, 130 B | 0.9971 [0.9890, 1.0053] | 0.9998 [0.9928, 1.0068] |
| low 512/1536, 66 B | 1.0196 [0.9867, 1.0536] | 1.0047 [0.9878, 1.0219] |
| low 512/1536, 130 B | 0.9950 [0.9670, 1.0238] | 1.0003 [0.9718, 1.0296] |

All decode cells resolved to the production tiled transform decoder.  The
non-fused cells therefore exercise the changed GF16 range kernel in both encode
and decode; an encode win was not allowed to hide a decode regression.  The
64-byte and 128-byte fused cells do not execute the changed branch and serve as
unaffected boundary controls.

The 64/66 and 128/130 cells are native-layout, complete-symbol GF16 physical
boundaries.  They are not odd-byte application-payload performance claims.
The versioned padded-odd layout maps payloads 65 and 129 to these physical
sizes, but its pack/scatter overhead was not measured here.  Commit `39cc427`
and `test_gf16_padded_odd` provide correctness-only coverage for that API.

## Isolation and provenance

Timing ran on CPU 12 of an AMD Ryzen 9 9950X3D while its topology-verified SMT
sibling CPU 28 was reserved and idle.  The runner held the canonical stable
anchor and pair-wide filesystem plus Linux abstract-socket leases.  A five-
second 0/0 non-idle quiet-pair presample passed before timing.  The retained
generic isolation record has `accepted=false` because that generic predicate
requires positive work on the measured CPU; the runner separately and
explicitly required zero non-idle jiffies on *both* CPUs for this quiet
presample.  Every complete timed ABBA round recorded positive work on CPU 12
and exactly zero non-idle scheduler jiffies on CPU 28; the minimum sibling
observation interval was 38 jiffies.

The runner itself was singleton-pinned to CPU 12, and every child was launched
with `taskset -c 12`.  The fixed environment forced one OpenMP thread.  Exact
source commits, clean tracked worktrees, executable/archive/AVX2-object hashes,
SHA-256 hashes of the link recipe and CMake cache, command, environment, raw
samples, resolved path, digests, and per-round `/proc/stat` snapshots are
retained.  The link-recipe and CMake-cache contents themselves are not retained.

The uncompressed raw JSON SHA-256 is
`a5c2ae2b377799558da48fe23c81affd3c8e58037aef173205f3d471b6ced4d4`.
It is stored deterministically as `raw.json.gz` (`gzip -n -9`) with SHA-256
`2fa051ff1530474f2b4f19ae76142cbceb2dd609fac03dbddad15c4e9cac08e5`.

## Reproduction and inspection

Build the exact clean pair in detached worktrees, then run only after reserving
one quiet physical core and its complete SMT sibling:

```sh
git worktree add --detach /tmp/leopard-gf16-control a725bf4
git worktree add --detach /tmp/leopard-gf16-candidate 7921998
cmake -S /tmp/leopard-gf16-control -B /tmp/leopard-gf16-control/build/neighbor-release -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON -DLEO2_ENABLE_CUDA=OFF
cmake --build /tmp/leopard-gf16-control/build/neighbor-release --target bench_leopard2 -j "$(nproc)"
cmake -S /tmp/leopard-gf16-candidate -B /tmp/leopard-gf16-candidate/build/neighbor-release -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=OFF -DLEO2_BUILD_BENCHMARKS=ON -DLEO2_ENABLE_CUDA=OFF
cmake --build /tmp/leopard-gf16-candidate/build/neighbor-release --target bench_leopard2 -j "$(nproc)"
taskset -c CPU python3 experiments/leopard2/gf16_high_encode/run_neighbor_abba.py --control-root /tmp/leopard-gf16-control --candidate-root /tmp/leopard-gf16-candidate --output /tmp/gf16-neighbor-evidence --cpu CPU --sibling SIBLING
```

The runner creates the predeclared `manifest.json` before measurement and
fails closed on source dirtiness, binary mutation, wrong profile/field/backend,
digest mismatch, incomplete ABBA order, invalid topology, or one non-idle
sibling jiffy.  The adjacent adversarial replay test validates all 156 retained
results and proves that mutations to requested profile, selected decode path,
or selected decode rule are rejected:

```sh
python3 -m unittest -v experiments/leopard2/gf16_high_encode/test_run_neighbor_abba.py
```

To inspect the retained raw data without modifying it:

```sh
gzip -cd experiments/leopard2/gf16_high_encode/results/7921998-neighbor-amd9950x3d/raw.json.gz | jq '.cells | length'
jq -r '.cells[] | [.cell.id, .encode.speedup_control_over_candidate, .decode.speedup_control_over_candidate] | @tsv' experiments/leopard2/gf16_high_encode/results/7921998-neighbor-amd9950x3d/summary.json
```

Exact Leopard-main comparison remains a separate gate and is not inferred from
this candidate-versus-parent result.
