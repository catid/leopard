# Algorithm 5 copy/no-copy evidence

This directory closes the attribution gap left by the Algorithm 5 evaluator
copy removal.  It does not change the installed library.  A private benchmark
links the non-installed `leopard_test_hooks` archive and selects either the
immutable-source evaluator or its retained whole-coefficient-block copy
fallback.  Both roles execute one binary compiled from one clean source
commit.  The ordinary `leopard` archive contains neither the selector symbols
nor the selector branch; `leopard2_test_hook_isolation` checks that boundary.

The benchmark emits schema `leopard2-benchmark-v4`.  In addition to ordinary
parity/recovery digests and selected-decoder-path metadata, it records output
block, mature-versus-pruned evaluator, out-of-place butterfly, and
compatibility-copy counters.  It fails
unless no-copy uses at least one out-of-place layer and zero fallbacks, or copy
mode uses exactly one fallback per output block and zero out-of-place layers.
The contract test covers materialized and tiled GF8/GF16, aligned shards,
GF8/GF16 compact tails, and balanced fully erased message blocks that must use
the nonpruned mature evaluator.  Adversarial mode, counter, path, and digest
mutations must be rejected.

## Build

Use a clean committed source tree.  The attribution build is diagnostic only:

    cmake -S . -B /tmp/leo2-high-copy-hooks \
      -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=ON \
      -DLEO2_BUILD_FUZZERS=OFF -DLEO2_ENABLE_CUDA=OFF \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build /tmp/leo2-high-copy-hooks -j "$(nproc)" \
      --target bench_leopard2_high_decode_copy_attribution

With compile-command export enabled, CTest also runs
`leopard2_high_decode_copy_build_identity`.  The same check is repeated before
every collected campaign: it resolves the actual link recipe, requires the
benchmark object and exactly one `leopard_test_hooks` archive, rejects the
production archive, and proves that both GF8 and GF16 selector symbols survive
in the private archive and executable.  The recursive
`leopard2_test_hook_isolation` test independently proves that neither selector
is present in the production archive.  Its benchmarks-enabled recursive mode
first builds the attribution executable and then proves that neither it nor the
private hook archive is installed.

Build an independent production candidate from the identical commit with test
hooks disabled, following `experiments/leopard2/main_compare/README.md`.  Build
the detached exact Leopard main adapter at commit
`6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198` as described there as well.

## Same-source ABBA

Create the canonical CPU reservation described by the exact-main runner, then
run the bounded 16-cell matrix.  It contains target, neighbor, nonpruned
full-block, and compact-tail cells for each field and workspace, and uses
ABBA/BAAB/ABBA order.  The collector
holds the coordinator reservation and the pair-wide lease, binds the clean
source/build/runtime closure, runs children with a strict environment, and
rejects any non-idle SMT-sibling jiffy.

    python3 experiments/leopard2/high_decode_copy/run_abba.py run \
      --source-root "$PWD" --source-commit "$(git rev-parse HEAD)" \
      --binary /tmp/leo2-high-copy-hooks/bench_leopard2_high_decode_copy_attribution \
      --hook-archive /tmp/leo2-high-copy-hooks/libleopard_test_hooks.a \
      --reservation-file build/leopard2-main-reservation.json \
      --output /tmp/leo2-high-copy-abba \
      --cpu 15 --reserved-sibling 31

The runner intentionally does not claim an exact-main result.  It is private
same-source attribution only.

## Independent exact-main campaigns

Run the production exact-main collector twice at the same Leopard2 source
commit, once per forced workspace.  Use only the six aligned cells for that
workspace from `matrix.json`.  For example, the materialized campaign uses:

    --candidate-mode materialized \
    --cell gf8-mat-target-one:240:16:65536:1:1701 \
    --cell gf8-mat-neighbor-full:224:32:4096:32:1703 \
    --cell gf8-mat-full-block:128:128:4096:128:1801 \
    --cell gf16-mat-target-one:1000:200:4096:1:1741 \
    --cell gf16-mat-neighbor-eight:4096:512:4096:8:1747 \
    --cell gf16-mat-full-block:256:256:4096:256:1811

The tiled campaign uses the corresponding `gf8-tiled-*` and `gf16-tiled-*`
cells and their seeds from `matrix.json`.  Supply all ordinary exact-main
runner arguments: detached baseline executable/archive/build/source,
production candidate executable/archive/build/source, reservation, CPU pair,
and a distinct output directory.

Finally bind the three independently verified campaigns:

    python3 experiments/leopard2/high_decode_copy/compose.py compose \
      --ab-manifest /tmp/leo2-high-copy-abba/manifest.json \
      --exact-main-materialized /tmp/leo2-main-materialized/manifest.json \
      --exact-main-tiled /tmp/leo2-main-tiled/manifest.json \
      --output /tmp/leo2-high-copy-composite.json

Composition reruns the current-input verifiers, requires both production
exact-main manifests to name the same source commit as the hook A/B campaign,
requires their complete K/R/shard-bytes/loss/seed cell metadata to match the
checked-in eligible subset, and compares parity, recovery, loss-set, mode,
backend, field, and path identities cell by cell.  The checked-in matrix
selects AVX2, so a production
exact-main campaign that resolves a different backend is valid on its own but
cannot be combined with this A/B campaign.  Non-64-byte tail cells are
explicitly classified as
`same-source-only-tail`; they can never acquire an inferred Leopard-main
result because the old API rejects those shard sizes.

All JSON digests in these runners are unkeyed integrity checks.  For
adversarial provenance, preserve the final bundles in an independently
authenticated store or transparency log.
