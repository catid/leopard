# Leopard2 C1 parent-preserving dependency pruning

Status: scalar experiment complete; promote the flat schedule concept to a
bounded C++/SIMD prototype, but do not enable this Python implementation or its
thresholds in production dispatch.

The implementation is
`experiments/leopard2/non_power_of_two/c1/dependency_pruning.py`. It answers the
C1 question left open by the symbolic C0 model: can Leopard omit work caused by
shortened inputs and punctured outputs while computing the exact existing
dyadic parent transform? The answer is yes for the tested scalar model. This
experiment does not define an exact-length code and does not change a wire
profile.

## Preserved mathematics and wire identity

For a transform of length `N = 2^n` in `GF(2^m)`, `n <= m`, the evaluator keeps:

- Leopard's public Cantor-coordinate representation;
- the normalized LCH basis `Xbar_i`;
- the enclosing power-of-two additive subspace or aligned coset;
- the existing coordinate order and skew factors;
- every padded input as an explicit mathematical zero; and
- every unrequested output as present in the parent code, merely not computed.

The scalar skew generator reconstructs the field multipliers represented by
`FFTSkew` in `LeopardFF8.cpp`. Leopard stores their logarithms for table
multiplication; C1 stores the equivalent field elements so multipliers zero and
one are visible to schedule construction. The implementation then uses the
same forward butterfly

    x' = x + m*y
    y' = y + x'

and inverse butterfly

    y' = y + x
    x' = x + m*y'

in Leopard's existing stage order. A shift is an aligned additive-coset
coordinate, not a renumbered block. Child shifts are `shift` and
`shift + N/2`, where addition is XOR in the stored Cantor coordinates.

An independent oracle constructs every subspace polynomial directly,
normalizes every LCH basis polynomial, converts LCH coefficients to a monomial
polynomial, and evaluates it by Horner's rule. It shares neither the skew
generator nor butterfly graph with the candidate. This follows the active
subspace definitions in `docs/leopard2_math_and_sources.md` and the LCH basis in
R16/R17. Existing Leopard source is used only to freeze its schedule and wire
conventions.

## Exact pruning analysis

Plan setup performs two passes over the actual scalar butterfly DAG.

1. A forward structural pass marks whether each intermediate can be nonzero
   from the supplied active-input mask. Coefficient zero is honored exactly.
2. A backward pass starts at the requested-output mask and marks which incoming
   intermediates are required. It uses each butterfly's four exact matrix
   coefficients rather than assuming that both outputs depend on both inputs.

The intersection is a proof-carrying execution schedule: an operation is kept
only when a structurally live input contributes to a requested result. Identity
writes, all-zero outputs, coefficient-zero products, and unused single outputs
are removed during setup. The byte loop contains no mask branches. Plan setup
and byte execution remain separate in all forms.

The model supports arbitrary sparse input and output masks, although Leopard's
current encoding profiles primarily produce prefixes. A plan rejects input
data that violates its declared known-zero mask.

## Execution forms compared

| Form | Representation | Result |
| --- | --- | --- |
| Recursive | Pruned transform tree with a branch per visited node | Correct; small scalar branch overhead at tiny payloads. |
| Flat | Precompiled required-butterfly list | Correct; simplest C++ prototype candidate and usually the fastest or tied form. |
| Hybrid | Complete dyadic subtransforms plus scalar boundary operations | Correct; no clear Python advantage, but exposes regular regions to existing radix-4/SIMD kernels. |
| Generated | Experiment-only Python source with indices and constants embedded | Correct for selected common schedules; 6.9-8.2x faster than the generic scalar-symbol interpreter when called on prevalidated input, but source grows from 902 bytes at N=16 to 25,086 bytes at N=256. |

The generated callable validates the known-zero input mask by default, like the
other execution forms. Its isolated micro-row explicitly disables duplicate
validation after the benchmark has constructed valid input. No generated code
is used by the library or default build. A production study
should use templates, constexpr data, or offline generation with an explicit
code-size budget; it must not add runtime executable memory.

## Correctness evidence

The retained deterministic result is
`experiments/leopard2/non_power_of_two/c1/results/self_test.json`. Its SHA-256
is:

    6632eeef83cbb308ab5cc206b1b319aa334ed396ef6c125c1d929b687caa0aea

The 2026-07-16 run passed:

| Check | Count |
| --- | ---: |
| Direct-polynomial versus full LCH transform vectors | 134,404 |
| Full forward/inverse round trips | 134,404 |
| Every GF4 prefix geometry, forward and inverse | 816 plans |
| Every GF4 sparse output mask through N=8 for every input prefix | 2,400 plans |
| Every sparse input/output mask through GF4 N=4, plus all N=8 input masks against ten irregular output masks | 5,672 plans |
| GF8 dyadic boundary/prefix sweeps through N=256 | 856 plans |
| GF8 arbitrary sparse-input sweeps through N=256 | 144 plans |
| Recursive/flat/hybrid/generated requested-output comparisons | 29,687 executions and 187,666 symbols |
| Arbitrary-tail byte comparisons at 1, 7, 31, 65, 257, and 1,025 bytes | 253,686 bytes |
| Generated-kernel known-zero-mask rejection regression | 1 |

GF4 covers every aligned coset for each supported parent. GF8 covers the first
and last aligned cosets for sub-field-size parents and the full N=256 field.
For N=1 and N=2, the oracle sweep includes all field values; consequently the
GF8 N=2 portion alone covers all 65,536 two-coefficient messages. Inputs include
all-zero plans, empty output masks, prefixes on both sides of powers of two,
alternating sparse masks, edge-only masks, and deterministic irregular masks.

The deterministic self-test was run twice after the final execution changes
and reproduced the same JSON hash. Python development-mode checks and bytecode
compilation also passed. This isolated experiment allocates managed Python
objects, so C/C++ ASan and UBSan are not applicable; production translation
still requires the repository sanitizer matrix.

## Operation and memory results

Counts below are exact for the scalar schedules, per payload byte. Padded counts
already specialize Leopard's zero multiplier and identity destination rather
than charging a deliberately naive butterfly.

| N | Active inputs | Requested outputs | Kept / padded butterflies | Fixed multiplies pruned / padded | XORs pruned / padded | Loads+stores pruned / padded | Flat schedule |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 9 | 7 | 12 / 32 | 5 / 17 | 16 / 49 | 40 / 113 | 192 B |
| 64 | 33 | 17 | 63 / 192 | 32 / 129 | 80 / 321 | 206 / 705 | 1,008 B |
| 256 | 129 | 65 | 319 / 1,024 | 192 / 769 | 448 / 1,793 | 1,086 / 3,841 | 5,104 B |

These cases remove 62.5-68.8% of butterflies and 64.6-71.7% of modeled shard
loads and stores. The hybrid representation finds one complete regular subtree
in each case: 4, 32, and 192 butterflies respectively, leaving 8, 31, and 127
boundary butterflies. This supports calling existing fused kernels for complete
subtrees instead of interpreting every operation individually.

## Scalar timing evidence

The retained host-dependent result is
`experiments/leopard2/non_power_of_two/c1/results/benchmark_amd_9950x3d.json`.
Its SHA-256 is:

    d92695c50efe9a5a5f835066f0b815bfeea1f3c2a9ffdecb766a860a9b12ac2d

The run used Python 3.12.3 on an AMD Ryzen 9 9950X3D. CPU 0 was isolated with
`taskset`; the process affinity contained one CPU while the host exposed 32.
Other agents paused CPU- and memory-heavy work for this short phase. Every
form received two warmups, nine interleaved repetitions through 1 KiB and seven
above it. Execution order rotated each repetition. The maximum median absolute
deviation was 1.33% of the median; every other cell was below 0.7%. Setup is
excluded. Timed byte transforms are forward transforms; inverse transforms are
covered by the exhaustive correctness program.

| N | Active / requested | Shard bytes | Padded | Recursive | Flat | Hybrid | Best speedup |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 9 / 7 | 64 | 101.1 us | 36.9 us | 34.6 us | 36.1 us | 2.92x |
| 16 | 9 / 7 | 1,024 | 1.533 ms | 0.546 ms | 0.545 ms | 0.546 ms | 2.81x |
| 16 | 9 / 7 | 16,384 | 25.717 ms | 9.199 ms | 9.197 ms | 9.174 ms | 2.80x |
| 16 | 9 / 7 | 65,536 | 103.437 ms | 37.319 ms | 36.980 ms | 37.081 ms | 2.80x |
| 64 | 33 / 17 | 64 | 624.4 us | 177.2 us | 167.5 us | 181.1 us | 3.73x |
| 64 | 33 / 17 | 1,024 | 9.816 ms | 2.856 ms | 2.848 ms | 2.858 ms | 3.45x |
| 64 | 33 / 17 | 16,384 | 165.357 ms | 48.630 ms | 48.599 ms | 48.465 ms | 3.41x |
| 256 | 129 / 65 | 64 | 3.406 ms | 0.921 ms | 0.882 ms | 0.956 ms | 3.86x |
| 256 | 129 / 65 | 1,024 | 54.481 ms | 15.347 ms | 15.300 ms | 15.375 ms | 3.56x |

These are deliberately scalar Python timings. They confirm that the reduced
memory and arithmetic work survives schedule traversal overhead, but they do
not predict a C++ SIMD speedup, cache behavior, or production crossover. The
large gap between generated and generic single-symbol interpreters is mostly
Python object/property overhead; it justifies specialization research, not a
7-8x production claim.

## Disposition

Promote the following into a bounded C++ prototype behind tests and explicit
dispatch:

- exact forward-live/backward-needed plan construction;
- flat boundary operation lists;
- complete-subtransform descriptors that call existing fused kernels; and
- zero/one multiplier plus identity-write specialization.

Do not promote the Python executor, measured timing thresholds, or generated
source. Keep recursive execution as an oracle/debug form. The hybrid form is
algebraically useful but not independently faster in this scalar environment;
its production case depends on replacing complete groups with real radix-4 or
SIMD kernels. Generated common-pair kernels remain a separate code-size-limited
experiment.

This clears the C1 scalar correctness and 10% experiment thresholds in the
tested region, but production promotion remains blocked on:

- C++ GF8 and GF16 translation against existing scalar/SSSE3/AVX2 outputs;
- encode/decode integration using real profile masks and shifted blocks;
- arbitrary byte tails, aliasing, scratch, and immutable-plan concurrency tests;
- sanitizer, fuzz, and backend-determinism gates;
- end-to-end codec benchmarks with plan setup and reuse reported separately;
- code/table footprint and instruction-cache measurement; and
- a dispatcher study covering neighboring regions where little work prunes.

## Reproduction

Run the deterministic correctness program:

    PYTHONHASHSEED=0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c1/dependency_pruning.py \
        self-test --output /tmp/c1-self-test.json

Compare it with the retained artifact:

    sha256sum /tmp/c1-self-test.json

Run the scalar timing program on one otherwise idle allowed CPU:

    PYTHONHASHSEED=0 taskset -c 0 python3 -X dev \
        experiments/leopard2/non_power_of_two/c1/dependency_pruning.py \
        benchmark --output /tmp/c1-benchmark.json

CPU 0 is valid only on the recorded host. On another machine, select one CPU
from the process's actual affinity mask. Timing JSON is expected to vary by
host; mathematical counts and the self-test JSON are deterministic.

No default CMake target imports this experiment, and it requires only the
Python standard library.
