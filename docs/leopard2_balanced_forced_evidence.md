# Balanced forced-decoder evidence protocol

The current protocol compares the generic active-parent decoder, materialized
Algorithm 5, and tiled Algorithm 5 in one production binary. It is a
measurement tool; it does not change production dispatch.

## Staged promotion matrix

`plan_balanced_promotion.py` makes the current-source promotion funnel
deterministic without treating a large blind sweep as authoritative evidence.
Its first stage contains only exact Leopard main versus current-source generic
decoder cells, sharded by transform side. A `(K, shard_bytes)` cell advances
only when the 95-percent lower confidence bounds for both first-use and
reuse-amortized decode speedup are at least 1.05. The selector authenticates the
schema-v4 exact-main bundles and derives the survivor set from their retained
analysis; an edited or merely re-signed plan, threshold, matrix, or survivor
file is rejected by canonical regeneration.

The initial boundary set is
`5,7,8,9,14,15,16,17,29,30,31,32,33,62,63,64,65,96,112,120,124,125,126,127,128`.
The interior `T=128` points prevent an apparent win at `K=128` from being
generalized across the whole transform bucket. The exact-main gate uses 256 B,
4 KiB, and 64 KiB aligned shards, for 75 cells and 900 timed child invocations
in five independently rerunnable transform-side shards. Exact-main runs three
repeated ABBA rounds. The same-binary forced runner instead uses ABBA, BAAB,
ABBA; the plan records these distinct order contracts explicitly.

If adjacent measured outcomes within one transform side and byte size disagree,
the boundary is resolved. If opposite outcomes still have unmeasured `K` values
between them, selection emits exact `runner_cell` refinements for every missing
integer and `advance` refuses to create later work. Once resolved, forced-mode
matrices contain only exact surviving `(K, shard_bytes)` cells; one survivor
does not activate its whole transform group.

Generate and validate the immutable-input plan outside the source tree:

    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        generate --output /tmp/leopard2-balanced-promotion-plan
    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        validate --input /tmp/leopard2-balanced-promotion-plan

Each generated gate JSON supplies the `runner_cell` strings for one invocation
of `main_compare/run_abba.py run` with `--candidate-mode generic`; all the other
source, binary, archive, build, CPU-pair, and reservation arguments remain the
runner's required authoritative inputs. After the five gate bundles finish,
derive the authenticated survivor set by repeating `--gate-manifest` once per
manifest:

    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        select --plan /tmp/leopard2-balanced-promotion-plan \
        --gate-manifest /results/t8/manifest.json \
        --gate-manifest /results/t16/manifest.json \
        --gate-manifest /results/t32/manifest.json \
        --gate-manifest /results/t64/manifest.json \
        --gate-manifest /results/t128/manifest.json \
        --output /tmp/leopard2-balanced-survivors.json

Run any `required_refinement_cells` through the same exact-main protocol and
repeat selection with those additional manifests. A zero refinement count can
advance to the conditional stage:

    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        advance --plan /tmp/leopard2-balanced-promotion-plan \
        --survivors /tmp/leopard2-balanced-survivors.json \
        --output /tmp/leopard2-balanced-confirmation-stage
    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        validate --plan /tmp/leopard2-balanced-promotion-plan \
        --input /tmp/leopard2-balanced-confirmation-stage

The conditional stage separates aligned current-source confirmation from true
tails. Generic, AUTO, materialized, and tiled exact-main confirmation uses
192 B, 4032 B, 4096 B, 64 KiB, and 1 MiB for every surviving `K`. The forced
same-binary tail matrices use only the non-aligned sizes 193 B, 4033 B, 4097 B,
and 1 MiB+1 on scalar, SSSE3, and AVX2. Exact Leopard main cannot execute those
tails, so they are correctness and internal-kernel evidence rather than an
external speedup claim.

AUTO rejection timing covers `L=1,4,floor(K/2),K-1`, `R=K-1`, and
`R=floor(K/2)` at 4 KiB and 64 KiB for every surviving `K`. Timing and workload
digests do not prove which AUTO rule ran. `path-attestation.json` is therefore
an exact-cell policy, not a command sketch or a transform-bucket inference. It
contains every evaluated full-loss gate cell, the loss/rate controls above,
the additional aligned confirmation sizes, and every true tail. Only an exact
`(K, shard_bytes)` gate cell whose two external confidence bounds passed may
select `generic/balanced_generic`. Rejected full-loss gates, all loss/rate
neighbors, the aligned sizes not present in the original gate, and every true
tail must select neither. A later promotion may change those expectations only
by generating a new canonically authenticated stage from evidence for those
exact cells.

The path result is collected rather than inferred from timing. The collector
refreshes `bench_leopard2`, requires a clean source at the survivor candidate
commit, records the benchmark source/object/archive/executable identity, runs
every canonical AUTO command, and validates exact K, R, seed, profile, field,
requested and resolved backend, byte count, loss coordinates, and selected
path/rule from schema-v3 benchmark JSON. The retained result binds the candidate
commit plus the plan, stage, survivor, and attestation digests. Verification
requires the same live source/build identity, exact raw file set, and a
deterministic reconstruction of every record; missing, duplicate, extra,
stale, edited, or merely re-signed outputs fail closed.

After the conditional timing work, collect and independently verify AUTO
selection from the exact candidate commit named by the stage. A dispatcher edit
changes that identity, so it requires regenerated exact-main evidence and a new
stage rather than reusing evidence for its predecessor:

    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        collect-attestation \
        --plan /tmp/leopard2-balanced-promotion-plan \
        --stage /tmp/leopard2-balanced-confirmation-stage \
        --source-root "$PWD" --binary build-release/bench_leopard2 \
        --output /tmp/leopard2-balanced-path-result
    python3 experiments/leopard2/decoder_dispatch/plan_balanced_promotion.py \
        verify-attestation \
        --plan /tmp/leopard2-balanced-promotion-plan \
        --stage /tmp/leopard2-balanced-confirmation-stage \
        --source-root "$PWD" --binary build-release/bench_leopard2 \
        --manifest /tmp/leopard2-balanced-path-result/manifest.json

An output directory placed inside the checkout must already be ignored. The
collector intentionally has no resume or online-tuning mode: path selection is
cheap and deterministic, so a partial directory is discarded and recollected
after fixing the cause. Timing still uses the isolated runners described below.

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
