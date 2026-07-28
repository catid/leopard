# L3-aware GF16 tiling evidence runner

`run_abba.py` compares two caller-named Leopard2 lanes on either measured CPU
pair (CPU 4 with reserved SMT sibling 20, or CPU 14 with sibling 30). It uses
explicit GF16, AVX2, one thread, retained samples, and the fixed L3-boundary
matrix. Three independent ABBA rounds are mandatory.

Every lane names an exact clean Git worktree and its full commit. An existing
binary may be supplied:

    python3 experiments/leopard2/l3_tiling/run_abba.py \
      --lane fixed32 /tmp/leopard-fixed32 COMMIT build/release/bench_leopard2 \
      --lane candidate /tmp/leopard-candidate COMMIT build/release/bench_leopard2 \
      --comparison fixed32 candidate \
      --cpu 4 --sibling 20 \
      --output /tmp/leopard-l3-cpu4

Alternatively, replace either `--lane` with `--build-lane LABEL ROOT COMMIT`.
The runner configures a Release Ninja build with CUDA and tests disabled, builds
`bench_leopard2`, and then freezes the binary in the campaign output.

For a targeted contrast, add `--cells` followed by one or more IDs from the
canonical matrix, for example:

    --cells low-side-256-loss9 high-side-512-maxloss

The manifest retains the exact selected subset in canonical matrix order.

Run only when an authoritative campaign is intended. The runner acquires
`/tmp/leopard-gf8-authoritative.lock` across builds, checks, and timing. It also
holds the shared evidence-runner lease and a physical-pair lease, validates the
benchmark's compiled-in source identity, verifies the frozen executable before
and after every child, rejects sibling contamination, and never summarizes a
failed campaign.

The reported speedup is left-lane time divided by right-lane time, so a value
above one means the right lane is faster. Lane labels and comparison direction
are explicit to avoid assigning baseline/candidate semantics implicitly.
