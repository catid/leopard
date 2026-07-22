# General GF8/AVX2 one-loss direct-repair experiment

This experiment asks whether the existing exact direct-repair algebra should
cover two additional GF8/AVX2 regions: unequal legacy-high parents and LOW
profiles with `R > K`.  It is intentionally limited to one missing original.
The CMake option `LEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1` is OFF by
default, and its definition is source-scoped to `leopard2.cpp` (plus the
focused API test when tests are built).

For one loss, plan construction generates one exact systematic generator row,
selects one received parity equation, and stores at most `K` execution terms.
The immutable codec retains `K` barycentric weights.  Row generation costs
`O(D+K)` field operations for unequal HIGH and `O(P+K)` for LOW; unequal HIGH
codec construction also computes `K` denominators across `D` parent-message
coordinates once.  Execution initializes
the missing output from one selected parity term, then streams the `K-1`
surviving originals through fixed-coefficient AVX2 multiply-adds.  It uses no
full-shard scratch.  The candidate therefore has approximately `K` source
streams and one output stream, but its output is repeatedly read/modified in
cache; this is not a claim of `(K+1) * bytes` DRAM traffic.  Above 64 KiB the
experimental executor keeps one 64-KiB output tile resident while streaming
all source terms through it.  Focused correctness cases bracket that traversal
change at 65,536 and 65,537 bytes.

For these nontrivial shapes Leopard1 exposes `N` decode-work shard buffers.
At a 1 MiB shard that is 64 MiB for `(17,31)` and 256 MiB for every `N=256`
row in the table.  The direct candidate removes that shard-data workspace;
only bounded plan/call metadata remains.  The runner records the exact scratch
bytes from both same-source executables rather than inferring them from this
upper-level comparison.

The following operation counts are structural proxies, not interchangeable
cycle estimates.  “Specialized butterflies” is `N/2 * log2(T)` for HIGH and
`N/2 * log2(P)` for LOW.  “Full-parent butterflies” is one `N`-point radix-2
transform, `N/2 * log2(N)`; Leopard1 decoding performs more than one transform
and additional weighting/memory passes.  A direct coefficient application is
also not identical in cost to a butterfly.  The table is useful only for
prioritizing measurements.

| Profile | K | R | P/T | N | Direct source terms | Specialized butterfly proxy | One full-parent transform |
|---|---:|---:|---:|---:|---:|---:|---:|
| HIGH | 240 | 16 | T=16 | 256 | 240 | 512 | 1024 |
| HIGH | 224 | 32 | T=32 | 256 | 224 | 640 | 1024 |
| HIGH | 192 | 64 | T=64 | 256 | 192 | 768 | 1024 |
| HIGH | 200 | 30 | T=32 | 256 | 200 | 640 | 1024 |
| LOW | 17 | 31 | P=32 | 64 | 17 | 160 | 192 |
| LOW | 31 | 200 | P=32 | 256 | 31 | 640 | 1024 |
| LOW | 32 | 224 | P=32 | 256 | 32 | 640 | 1024 |
| LOW | 127 | 128 | P=128 | 256 | 127 | 896 | 1024 |

The direct path should win most clearly for redundancy-dominant LOW shapes and
large enough shards to amortize plan setup.  It can lose for tiny shards due
to generator-row/term setup, for large `K` when repeated output updates or
coefficient-table pressure dominate, and anywhere a mature transform has
better cache reuse than the simple operation proxy predicts.  `R == K` LOW
cells remain negative selector controls; `L=2` is probed explicitly and must
fall back to a transform.

## Prioritized timing matrix

The dedicated runner `run_general_l1_abba.py` uses these tiers:

1. `core`: the eight rows above at 64 B, 2 KiB, 4 KiB, 64 KiB, and 1 MiB,
each at plan reuse 1, 8, and 64 (120 cells).
2. `balanced-neighbor`: LOW `(17,17)`, `(31,31)`, `(32,32)`, `(127,127)`,
   and `(128,128)` over the full byte/reuse grid (75 selector-control cells).
3. `shape-neighbor`: valid GF8 `K`/`R` neighbors plus explicit legacy-high
   `(17,33)` and `(17,128)` at 4 KiB and 64 KiB with reuse 8 (44 cells).
4. `byte-neighbor`: each core shape at +/-1 around 64 B, 2 KiB, 4 KiB,
   64 KiB, and 1 MiB with reuse 8 (80 cells).

Before timing, eight 65-byte `L=2` probes must select a non-direct path.  Each
timed cell uses five serial ABBA rounds by default, pins one physical core,
requires an idle SMT sibling, and reports codec setup, plan setup, execution,
cold codec+plan+execution, and plan setup amortized at the declared reuse
count.  Codec setup is not hidden inside plan reuse.  A candidate clears the default gate
only when the lower 95% confidence bound is at least 1.05x.  A credible
neighbor regression below 0.98x is retained separately.

The runner configures same-source OFF and ON builds, permits exactly one
compile/object delta (`leopard2.cpp`), excludes AVX-512 at configuration and
disassembly, freezes and hashes both executables, and rechecks hashes on
resume.  An optional exact Leopard-main benchmark runs only after a valid
legacy-high cell first wins the same-source comparison.  LOW and non-winning
cells cannot invoke it.

Example commands after acquiring the canonical campaign lease:

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py self-test
    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py make-matrix \
        --output /tmp/leopard2-general-l1-matrix.json
    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-build \
        --run-dir /tmp/leopard2-general-l1-results

Pass `--main-benchmark`, the exact 40-character `--main-commit`, and its
pre-recorded 64-character `--main-sha256` only when a verified Leopard-main
adapter is available.  Build and result roots must be outside the source
checkout.
