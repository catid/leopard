# Leopard2 / Jerasure comparison protocol

> **Evidence status:** no correctness or timing artifact is committed.  A
> correctness campaign generated from the final integrated source belongs at
> the ignored `results/leopard2/jerasure/correctness.json` path and is not a
> performance result.  No Jerasure versus Leopard2 throughput claim is
> accepted until the separate pinned ABBA timing command completes and its
> artifact passes both portable and trusted-cache validation.

This experiment compares equivalent public erasure workloads, not wire
formats.  Jerasure and Leopard2 receive the same `K`, `R`, shard byte count,
deterministic source bytes, missing-original indices, batch size, reuse count,
and single-thread constraint.  Their fields/bases, coordinate sets, generator
matrices, and parity bytes can differ.  Only application-payload,
generated-output, and repaired-output rates are comparable.

## Entirely optional production dependency

The production CMake graph does not descend into
`experiments/leopard2/jerasure_compare`.  Its only normal-build integrations
are dependency-free Python policy/parser tests; they neither find, compile,
link, install, nor load either third-party library.  The full mutation campaign
is a manual evidence gate because it consumes the ignored correctness artifact.
The standalone project must be configured explicitly and builds private static
archives in the ignored research cache.  The full BSD-3-Clause notices are
retained beside the adapter.

The pinned sources are:

| Component | Source identity | License identity |
|---|---|---|
| Jerasure 2.0 | https://github.com/ceph/jerasure at `de1739cc8483696506829b52e7fda4f6bb195e6a`, tree `fb98f85c548038a5ff294141f89603dda70dd423` | BSD-3-Clause, `License.txt` SHA-256 `83b6b3ff237848fbccfa889bb52cfb13c331c1f83544b907617c7f8f31eb1769` |
| GF-Complete | https://github.com/ceph/gf-complete at `a6862d10c9db467148f20eef2c6445ac9afd94d8`, tree `5a13169b93b6e517184fbdf39033098b329d68a6` | BSD-3-Clause, `License.txt` SHA-256 `cb9790699b4a3d56a43bba1dd859f7f41361cd224e8745a24eef933ea134a280` |

Bootstrap rejects dirty or mismatched source trees.  It creates a detached
worktree at the exact clean Leopard commit and builds both benchmark
executables from that materialization.  Provenance records normalized compile
commands, actual CMake link commands, every object/archive input consumed by
those commands, the three private static archive hashes, compiler and build
tool identities, each executable hash, and the transitively closed ELF
interpreter/`DT_NEEDED` graph.  A trusted-cache replay reconstructs and exact
compares that closure.  Every replay also re-verifies the detached Leopard
HEAD, tree, clean status, tracked-tree listing, and detached-HEAD state.  The
source identities are included in the build identity whose object/archive and
executable closure is checked, so a coordinated source relabel cannot stand in
for the code that produced an executable.

Both standalone builds explicitly select the `Unix Makefiles` generator.  The
build identity validates exact Release configuration definitions and exact
adapter/Leopard2 targets, rejects extra recipe claims, and records nonempty
absolute paths plus hashes and reported versions for every build/link
inspection tool.  Each private static archive name is bound to one exact,
distinct normalized path and to a `static-archive` input actually consumed by
the adapter link.

## Codec and fairness semantics

The adapter uses Jerasure's documented
`reed_sol_vandermonde_coding_matrix` and `jerasure_matrix_encode` path.  The
field is GF(2^8) when `K+R <= 256` and GF(2^16) otherwise.  This comparison is
deliberately bounded to `K <= 4096`, `R <= 4096`, at most `2^24` coefficients
in any dense oracle/decode matrix, at most 8 GiB of application buffers, and a
paired Leopard2 parent no larger than 65,536 coordinates.  Both the executable
and matrix classifier reject cells outside that domain.  GF-Complete's
standard polynomials are
`0x11d` and `0x1100b`, respectively.

The independent oracle does not call Jerasure or GF-Complete arithmetic.  It
constructs the extended Vandermonde evaluations, independently inverts the
systematic prefix, applies Jerasure's documented generalized-RS row
normalizations, and checks all `K*R` coding coefficients.  It then recomputes
every requested parity symbol using scalar polynomial-basis multiplication.
The full correctness campaign checks every parity byte; timing children use a
deterministic boundary/random projection and are gated by the separate full
artifact.

Decode-plan setup calls `jerasure_make_decoding_matrix` once, with surviving
systematic rows followed by the lowest parity rows until exactly `K` rows are
selected.  Execution calls `jerasure_matrix_dotprod` only for requested
missing originals.  Parity rebuilding is excluded.  This separates matrix
construction/inversion from byte-heavy execution in the same way Leopard2
separates codec/plan setup from execution.  Offered received bytes and the
selected `K`-row byte count are recorded separately.

A no-loss decode constructs no decode matrix, selects and offers zero bytes,
and executes no repair dot product.  Its timing is latency/status evidence
only: input and output throughput fields are null, so the harness cannot claim
an offered-byte rate for a no-op.

Every application shard begins on a 64-byte boundary and shard length must be
a multiple of eight bytes, satisfying Jerasure's documented longword region
contract deterministically on the intended 64-bit comparison hosts.  Invalid
lengths fail closed; they are not silently padded.  Both codecs operate
directly on application buffers, and the adapter records zero staging or
format-conversion bytes.  The benchmark emits exact FNV-1a-64 digests for
original input, transmitted parity, and recovered originals.  Paired timing
children must match original and recovered digests and their independently
regenerated loss lists; parity digests are retained but intentionally need not
match across the two different codes.

Ambient loader, allocator, OpenMP, sanitizer, math-runtime, and
`GF_COMPLETE_*` dispatch controls are rejected before collection.  Children
receive a fixed one-thread environment.  Authoritative timing additionally
requires an explicitly reserved SMT pair plus at least one additional CPU in
the runner's launch affinity for housekeeping.  A pair-only launch affinity is
rejected before collection.  The runner then requires singleton CPU affinity
for every timed child and an explicitly reserved SMT sibling.  Both
logical CPUs must mutually report exactly that two-CPU pair through Linux
`thread_siblings_list`, with matching core and package identities.  The runner
holds both its cache-local coordinator locks and the cross-comparator,
pair-wide lock at
`/run/user/UID/leopard2-cpu-leases/leopard2-cpu-pair-UID-A-B.lock`.  The
per-user runtime root and child directory must be owned by that UID at mode
0700; the canonical lease file must be an owned, single-link mode-0600 regular
file.  Its device and inode are bound into evidence and the pathname, open
descriptor, permissions, ownership, link count, and contents are rechecked
before every accepted pair.  This closes the unlink-and-replace gap of a flat
predictable advisory lock.  ABBA provider order, four independent repetitions,
two warmups, and nine retained timing samples are still mandatory.  Setup,
execution, and setup amortized at the declared reuse count remain separate.

GF-Complete is compiled at portable compiler defaults.  The standalone build
does not apply whole-archive `-msse*`, `-mavx*`, `-mpclmul`, or native-target
flags merely because the compiler accepts them; doing so would not prove that
the runtime CPU can execute instructions inserted before library dispatch.
The normalized compile closure fails validation if these broad flags reappear.

Timing state is durable and exact-bound to the output path, correctness
artifact, source/build/executable identities, cells, CPU pair, environment,
and method.  Each child document is validated and atomically persisted.  A
repetition becomes resumable only after both providers, their matching loss
and workload identities, per-pair pre/post frequency snapshots, singleton
affinity evidence, exact SMT topology, lease inode, and CPU-activity audit are
atomically bound.  Immediately around every two-child repetition the runner
records the first eight non-double-counted fields for both logical CPUs from
`/proc/stat`.  It accepts the pair only when the benchmark CPU accumulated
non-idle work, the interval spans at least one sibling scheduler jiffy, and the
reserved sibling accumulated zero non-idle jiffies.  Every delta and acceptance
decision is recomputed during portable validation.  A failed isolation check
is retained under the durable state's `failures/` directory and aborts rather
than resampling.  A partial two-provider repetition is discarded and rerun; a
complete but tampered one fails closed.  Affinity is rechecked after every
timing child.  Final and correctness JSON files also use atomic replacement.

These controls are cooperative and unprivileged.  They cannot stop root-owned
tasks or kernel work from being scheduled on the pair, and `/proc/stat` has
scheduler-jiffy resolution.  Observed non-idle sibling activity therefore
fails the campaign instead of being described as OS-exclusive isolation.
The additional allowed CPU supplies a place outside the reserved core for
ambient housekeeping but does not move or constrain kernel work by itself.
Activity on the benchmark logical CPU cannot be separated from the intended
child's work by these counters, so an external host reservation and an idle
machine remain required.  Version-2 timing artifacts carry these semantics.
Retained version-1 checkpoints remain
portable-validation compatible, and version-1 durable manifests remain
structurally parseable for diagnostics, but version 1 does not acquire the new
isolation guarantee retroactively and cannot be resumed as a version-2
campaign.

## Correctness evidence lifecycle

Generated correctness evidence is never written into the source experiment
directory and is never committed.  It belongs under the ignored
`results/leopard2/jerasure/` tree; the correctness command defaults to
`results/leopard2/jerasure/correctness.json`.  Evidence is regenerated whenever
adapter executable semantics or the build/source-closure schema changes.  The
authoritative campaign contains 128 deterministic cases, uses ten parallel
subprocess workers, includes no-loss and maximum-original-loss cases, exercises
GF8 and GF16, and checks every parity byte against the independent scalar
oracle.  Its canonical, serialized, Leopard-source, executable, and provenance
hashes are reported with the final evidence; hashes from an earlier adapter
revision must never be reused.

## Reproduction

Normal-build optionality and synthetic policy checks require no external
source:

    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
      -DLEO2_BUILD_TESTS=ON -DLEO2_BUILD_BENCHMARKS=OFF \
      -DLEO2_ENABLE_CUDA=OFF
    cmake --build build/release -j 10
    ctest --test-dir build/release \
      -R 'leopard2_(external_comparison_self|jerasure_(comparison_self|default_optionality))_test' \
      --output-on-failure

These default tests are dependency-free: they do not read a generated
correctness artifact and do not discover, compile, link, or load Jerasure or
GF-Complete.  The full adversarial mutation campaign is intentionally not a
default CTest because it validates a generated, ignored 128-case correctness
artifact.  Run it manually after generating that artifact, using the command
below.

Bootstrap and replay use only the ignored cache for third-party source and
build products:

    python3 tools/leopard2_jerasure_compare.py bootstrap \
      --cache .research/leopard2 --jobs 10
    python3 tools/leopard2_jerasure_compare.py correctness \
      --cache .research/leopard2 --cases 128 --workers 10 \
      --output results/leopard2/jerasure/correctness.json
    python3 tools/leopard2_jerasure_compare.py validate-correctness \
      results/leopard2/jerasure/correctness.json
    python3 tools/leopard2_jerasure_compare.py validate-correctness \
      results/leopard2/jerasure/correctness.json \
      --cache .research/leopard2 --require-local-build-match
    python3 tools/leopard2_jerasure_compare.py mutation-test \
      results/leopard2/jerasure/correctness.json

Only after reserving one physical core and its SMT sibling may timing run.  The
runner must be launched with both logical CPUs and at least one additional
housekeeping CPU in its allowed affinity.  On Linux, verify that
`/run/user/$(id -u)` exists, is owned by the current UID, and has mode 0700; the
runner creates its protected child lease directory itself:

    stat -c '%u %a %n' /run/user/$(id -u)

    python3 tools/leopard2_jerasure_compare.py run \
      --cache .research/leopard2 \
      --cpu <isolated-allowed-cpu> \
      --reserved-idle-cpu <its-idle-smt-sibling> \
      --correctness-artifact results/leopard2/jerasure/correctness.json \
      --output results/leopard2/jerasure/checkpoint.json
    python3 tools/leopard2_jerasure_compare.py validate \
      results/leopard2/jerasure/checkpoint.json \
      --correctness-artifact results/leopard2/jerasure/correctness.json
    python3 tools/leopard2_jerasure_compare.py validate \
      results/leopard2/jerasure/checkpoint.json \
      --correctness-artifact results/leopard2/jerasure/correctness.json \
      --cache .research/leopard2 --require-local-build-match

The first checkpoint validation is portable and proves internal closure and
derivation, including mutual topology records, retained lease UID/path/inodes,
every scheduler-counter delta, and the zero-non-idle-sibling decision; it does
not require the collector's UID or lease inode to exist on the replay host.
The second is the required trusted-cache replay: it re-hashes the local pinned
sources, detached Leopard materialization, toolchain/build/link closure, and
executables.  Both must pass before any timing value is cited.

No performance number should be quoted from the correctness artifact or from
an unpinned direct invocation of the standalone executable.
