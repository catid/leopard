# Leopard2 all-K performance atlas harness

This standalone experiment generates the data and SVGs consumed by the
repository's performance README. It compares:

- Leopard2 from a clean source-attested build;
- the exact Leopard `main` baseline at
  `6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198`; and
- Wirehair's shipping C API at
  `067ca7cdb66aed424ec23f97557429bf791c6f0c`.

The runner deliberately does not add Wirehair to Leopard's normal build. The
adapter is compiled as a separate process against a pinned checkout, preserving
both projects' license boundaries and keeping Wirehair entirely optional.

## Matrix

`R=32`; `K` is every odd value through 223 plus every power of two and endpoint
224. Shard sizes are 64 B, 1 KiB, 4 KiB, and 1 MiB. The loss patterns are one
and two deterministic-random missing originals, 10% rounded up and clipped to
R, and maximum source loss `min(K,R)`. Equivalent patterns are timed once and
tagged with all applicable labels.

This dense boundary stops at K=224 so the high-profile construction remains a
GF8 parent of N=256. It is an intentional GF8 atlas, not a claim to benchmark
the much larger GF16 maximum with 1 MiB shards.

## Build the standalone Wirehair adapter

Clone or create a detached Wirehair checkout in the ignored research cache,
verify its commit, and build the adapter plus static library through the
fail-closed AVX2 compact-path experiment project:

    git clone https://github.com/catid/wirehair \
      .research/leopard2/wirehair-v1
    git -C .research/leopard2/wirehair-v1 checkout --detach \
      067ca7cdb66aed424ec23f97557429bf791c6f0c
    cmake -S experiments/leopard2/performance_atlas \
      -B .research/leopard2/performance-atlas-build/wirehair \
      -DCMAKE_BUILD_TYPE=Release \
      -DWIREHAIR_SOURCE_DIR="$PWD/.research/leopard2/wirehair-v1"
    cmake --build .research/leopard2/performance-atlas-build/wirehair -j 1

The pinned shipping artifact contains target-qualified AVX-512 helpers even
when its ordinary translation units use `-mno-avx512*`; this revision has no
supported option to compile those helpers out. They are reachable only through
Wirehair's opt-in thread-wide XOR facility. The adapter forces that facility
off and verifies it remains off across every measured phase, so it attests the
measured v1 compact path as AVX2 without falsely calling the whole artifact
AVX2-only.

The standalone build registers three bounded, serial CTest checks: the atlas
Python unit suite, the generator self-test, and a small Wirehair round trip
that also verifies the adapter and library compile commands. Run them after
building:

    ctest --test-dir \
      .research/leopard2/performance-atlas-build/wirehair \
      --output-on-failure -j 1

Build Leopard2's `bench_leopard2_allk` from a clean detached worktree with the
project's production pure-AVX2 recipe (AVX-512 feature probes forced false).
Build the exact main adapter using `../main_compare/CMakeLists.txt`, an exact
detached main source, and `-DLEO_MAIN_PURE_AVX2=ON`, as documented in
`../main_compare/README.md`. The runner requires Leopard2 AUTO to resolve to
AVX2 in every result.

## Validate and run

    python3 experiments/leopard2/performance_atlas/test_generate_atlas.py -v
    python3 experiments/leopard2/performance_atlas/generate_atlas.py self-test
    python3 experiments/leopard2/performance_atlas/test_wirehair_adapter.py \
      --binary .research/leopard2/performance-atlas-build/wirehair/wirehair_v1_benchmark \
      --compile-commands .research/leopard2/performance-atlas-build/wirehair/compile_commands.json \
      --wirehair-source .research/leopard2/wirehair-v1 \
      --expected-commit 067ca7cdb66aed424ec23f97557429bf791c6f0c

Use `generate_atlas.py all --help` for the complete run contract. Every binary
path requires a caller-supplied SHA-256. The runner copies those bytes into its
own read-only directory under the canonical benchmark lock, checks Leopard2's
embedded clean source identity, pins every child to one allowed CPU, caps its
address space, rotates codec order deterministically, writes each result
atomically, and revalidates all evidence during resume and summary generation.

For a quick early low-K Leopard2/Wirehair validation before the full campaign,
add `--max-cells 4`. It does not exercise exact Leopard main or every shard
size, and any partial run is smoke evidence only; it cannot generate the final
performance README. The standalone CTest smoke covers adapter/build mechanics,
while only the complete matrix is release evidence.

The SVG generator depends only on the Python standard library.
