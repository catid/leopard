# Balanced forced-decoder evidence protocol

The current protocol compares the generic active-parent decoder, materialized
Algorithm 5, and tiled Algorithm 5 in one production binary. It is a
measurement tool; it does not change production dispatch.

This forced-path comparison answers only which Leopard2 decoder kernel is
cheapest. Exact Leopard main remains the primary external baseline for release
claims and dispatcher promotion. A forced-path winner must therefore be joined
to current-source exact-main evidence; an internal win cannot conceal a
Leopard2 regression against main.

`run_balanced_abba.py` accepts an external, versioned JSON matrix. Every case
is balanced legacy-high GF8 full-original recovery and names two distinct
members of `generic`, `materialized`, and `tiled`. Those names map to exact
force tuples:

| Mode | Required benchmark flags | Required four-Boolean output tuple |
| --- | --- | --- |
| generic | `--force-generic` | generic true; specialized, tiled, and materialized false |
| materialized | `--force-specialized --force-materialized` | specialized and materialized true; generic and tiled false |
| tiled | `--force-specialized --force-tiled` | specialized and tiled true; generic and materialized false |

AUTO is not a matrix mode. In particular, a row with all force fields false
cannot be called specialized merely because `force_generic_decode` is false.
The analyzer rejects that mutation.

Each case runs three independent clustered rounds in ABBA, BAAB, ABBA order.
Every child retains all timing observations. Its raw record envelope binds:

- clean Git commit and tree;
- benchmark, decoder, and dispatch source bytes;
- benchmark and decoder object bytes, production archive, executable, CMake
  cache, and build graph;
- the external matrix, runner, common validation module, and CPU-lease source;
- exact argv and raw benchmark JSON;
- Python, kernel, CPU topology, runtime-library, affinity, and environment
  identity; and
- an exclusive CPU-pair lease plus `/proc/stat` evidence requiring at least one
  timed-CPU non-idle jiffy and exactly zero non-idle jiffies on the reserved SMT
  sibling for every invocation.

The runner refreshes `bench_leopard2` with one build worker before measuring.
That intentional serial phase is outside timing and prevents a stale object or
archive from being labeled as current source. The timed child is pinned to one
physical core with one OpenMP thread. Broad correctness and experiment work
should return to the full allowed CPU set outside this isolated phase.

An interrupted run may be resumed with `--resume`. Existing record envelopes
are authenticated before they are reused; partial or failed-isolation records
are never overwritten or accepted. The final analyzer takes no source-commit,
binary-digest, or mode-identity arguments. It obtains those only from signed
raw evidence, validates them against live files, recomputes every raw statistic
and derived rate, and reports one mean log contrast per independent round with
a Student-t 95 percent interval.

The committed smoke matrix exercises all three pairings but is not performance
evidence. A smoke collection, after building a Release benchmark, is:

    python3 experiments/leopard2/decoder_dispatch/run_balanced_abba.py \
        --source-root "$PWD" --source-commit "$(git rev-parse HEAD)" \
        --binary build/release/bench_leopard2 \
        --matrix experiments/leopard2/decoder_dispatch/balanced_forced_smoke_matrix.json \
        --output /tmp/leopard2-balanced-smoke --cpu PHYSICAL_CPU \
        --sibling ITS_IDLE_SMT_SIBLING

    python3 experiments/leopard2/decoder_dispatch/analyze_balanced.py \
        --manifest /tmp/leopard2-balanced-smoke/manifest.json \
        --output /tmp/leopard2-balanced-smoke-summary.json

The runner fails if the pair is not the exact two-thread sibling set reported
by sysfs, if either CPU is outside process affinity, if no housekeeping CPU is
available, if another timing/build process is active, or if the sibling does
work. These failures are expected safeguards, not reasons to weaken the
evidence contract.

The older `balanced_amd_9950x3d.json` result remains historical provenance for
the narrow existing dispatch rule. It used two reversed runs, separate backend
binaries, and inferred specialized from a false generic flag. It cannot be
reclassified as current three-mode evidence or used to widen the dispatch
table.
