# General GF8/AVX2 one-loss direct-repair experiment

This experiment asks whether the existing exact direct-repair algebra should
cover two additional GF8/AVX2 regions: unequal legacy-high parents and LOW
profiles with `R > K`.  It is intentionally limited to one missing original.
The CMake option `LEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1` is OFF by
default, and its definition is source-scoped to `leopard2.cpp` (plus the
focused API test when tests are built).

The nested option
`LEO2_EXPERIMENTAL_GF8_AVX2_DIRECT_L1_PAIR_FUSION` is also OFF by default and
requires the general experiment.  It qualifies an AVX2 operation that forms
one output from two independently scaled sources in one pass.  The direct
executor then initializes from its first source pair and accumulates later
pairs, halving hot-output loads/stores for the paired terms.  The original
runner defaults this nested option OFF so its results continue to isolate the
simple direct candidate.  `--pair-fusion` measures the complete fused candidate
against the transform control, while `--pair-attribution` compiles simple
direct repair into both binaries and permits only the qualification and AVX2
backend objects to differ.  The latter comparison determines whether pair
fusion needs a byte threshold instead of routing every shard size through the
paired traversal.

The pair experiment now uses a deterministic, exact-shape selector derived
from the immutable pair-only campaign.  Ordinary builds still compile with
the option OFF.  An enabled experimental build selects pair fusion only at or
above these byte counts; every other K/R shape retains the simple executor:

| Profile | K | R | Minimum shard bytes |
|---|---:|---:|---:|
| HIGH | 192 | 64 | 64 |
| HIGH | 200 | 30 | 64 |
| HIGH | 224 | 32 | 64 |
| HIGH | 240 | 16 | 65,536 |
| LOW | 17 | 31 | 1,048,576 |
| LOW | 31 | 200 | 64 |
| LOW | 32 | 224 | 2,048 |
| LOW | 127 | 128 | 65,536 |

The focused API suite executes each source shape at threshold minus one and
at threshold, and executes every valid immediate K/R neighbor at the source
threshold.  The benchmark JSON reports the actual selector decision, so a
campaign cannot attribute a no-op control cell to the paired kernel.  This
uses the explicit `--report-pair-fusion` schema v6 mode; the older strict
path-report schema v3 remains byte-for-byte structurally unchanged.

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
   64 KiB, and 1 MiB with reuse 8.  The eight 63-byte coordinates are owned by
   the tiny-byte tier below, leaving 72 cells in this tier.
5. `pair-selector-neighbor`: immediate K/R neighbors evaluated at the exact
   source-shape threshold when that coordinate is not already covered by the
   balanced or shape-neighbor tiers.
6. `tiny-byte`: all eight measured HIGH/LOW shapes at
   1,2,3,7,8,15,16,17,31,32,33,63 bytes and reuse 1 and 8 (192 cells).
   Run this tier in both the default transform-versus-simple mode
   and `--pair-fusion` mode.  Pair fusion itself must remain unselected below
   64 bytes; these cells gate the byte-independent general direct-plan choice
   and coordinate with the separate XMM-tail experiment.

The grouped-output XMM-tail work does not modify
`AVX2FF8LinearCombination2`, whose paired loop currently begins at 32 bytes
and uses scalar cleanup.  If forced-pair attribution confirms a 15/31/63-byte
cliff, a separate default-OFF follow-up will evaluate 16-byte and 8-byte VEX
tails for its generic and identity-specialized variants.  The conservative
selector remains unchanged until that complete-codec evidence exists.

Before timing, eight 65-byte `L=2` probes must select a non-direct path.  Each
timed cell uses five serial ABBA rounds by default, pins one physical core,
requires an idle SMT sibling, and reports codec setup, plan setup, execution,
cold codec+plan+execution, and plan setup amortized at the declared reuse
count.  Codec setup is not hidden inside plan reuse.  A candidate clears the default gate
only when the lower 95% confidence bound is at least 1.05x.  A credible
neighbor regression below 0.98x is retained separately.

By default the runner configures same-source OFF and ON builds and permits
exactly one compile/object delta (`leopard2.cpp`).  In `--pair-fusion` mode the
two pair-backend objects are the only additional deltas.  In
`--pair-attribution` mode `leopard2.cpp` must be byte-identical and only those
two backend objects may differ.  Every mode excludes AVX-512 at configuration
and disassembly, freezes and hashes both executables, and rechecks hashes on
resume.  An optional exact Leopard-main benchmark runs only after a valid
legacy-high cell first wins the same-source comparison.  LOW and non-winning
cells cannot invoke it, and pair-only attribution never invokes it.

Example commands after acquiring the canonical campaign lease:

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py self-test
    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py make-matrix \
        --output /tmp/leopard2-general-l1-matrix.json
    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-build \
        --run-dir /tmp/leopard2-general-l1-results

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --pair-attribution \
        --tiers core,byte-neighbor,shape-neighbor,balanced-neighbor,pair-selector-neighbor \
        --cell-regex='-u8$' \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-pair-build \
        --run-dir /tmp/leopard2-general-l1-pair-results

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --pair-fusion --tiers core \
        --cell-regex='(high-core-k(192-r64-b64|200-r30-b64|224-r32-b64|240-r16-b65536)-u8)$' \
        --main-benchmark /absolute/path/to/frozen-main-benchmark \
        --main-commit FULL_40_HEX_MAIN_COMMIT \
        --main-sha256 FULL_64_HEX_MAIN_EXECUTABLE_SHA256 \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-main-build \
        --run-dir /tmp/leopard2-general-l1-main-results

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --tiers tiny-byte \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-tiny-simple-build \
        --run-dir /tmp/leopard2-general-l1-tiny-simple-results

    python3 experiments/leopard2/direct_repair/run_general_l1_abba.py run \
        --pair-fusion --tiers tiny-byte \
        --source "$PWD" --commit "$(git rev-parse HEAD)" \
        --build-root /tmp/leopard2-general-l1-tiny-pair-build \
        --run-dir /tmp/leopard2-general-l1-tiny-pair-results

Pass `--main-benchmark`, the exact 40-character `--main-commit`, and its
pre-recorded 64-character `--main-sha256` only when a verified Leopard-main
adapter is available.  Build and result roots must be outside the source
checkout.
