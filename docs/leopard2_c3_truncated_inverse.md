# Leopard2 C3: truncated inverse and basis conversion

## Result

C3 does not promote a production transform or wire profile.

The isolated scalar study proves two different inverse problems and keeps them
separate in every result row:

- The wire-compatible parent problem inverts a dyadic parent after fixing the
  unprovided evaluation suffix to zero.  A dense Lagrange-to-LCH column model
  matches the padded LCH inverse, but its deterministic operation model finds
  at most a 3.23% advantage.  Pinned Python timings expose a possible tiny-q
  direct region, so this candidate is **inconclusive for tiny q**, not promoted
  and not globally rejected.
- The exact-q problem interpolates the unique polynomial of degree less than q
  through exactly q prefix points.  Dense, Newton-to-LCH, and generic Lagrange
  constructions agree, and the dense path reaches a 52.94% best modeled
  advantage.  The sweep found 12,370 inputs for which it differs from the
  padded parent.  It is therefore **new-profile research**, not a
  wire-compatible optimization.
- The transparent Newton and generic Lagrange executors are slower than the
  dense executor in the pinned scalar cells.  They remain correctness
  references rather than production candidates.
- No fast Coxon recursive conversion or Tang-Han epsilon inverse was
  implemented.  Numerical agreement with the dense models is not used to
  claim either published algorithm.

The implementation is isolated in
`experiments/leopard2/non_power_of_two/c3/study.py`.  It is not referenced by
CMake or production headers and does not alter runtime dispatch.

## Mathematical boundary

Let N be the enclosing power of two, q be the active prefix length, B_N be the
N by N matrix that evaluates Leopard's normalized LCH basis on the shifted
parent coset, and y contain the q supplied values.

The parent-preserving inverse is

    c_parent = inverse(B_N) * [y, 0_(N-q)].

The dense parent model stores only the first q columns of inverse(B_N).  It is
therefore a direct Lagrange-to-LCH basis-conversion model and must match the
padded butterfly inverse exactly.  It does not change the code or parity.

The exact construction instead forms the q by q matrix from the first q LCH
basis polynomials and the q prefix points:

    c_exact = inverse(B_q) * y.

It represents a degree-less-than-q polynomial.  Unless q=N, it does not impose
zero values on the parent suffix and is not the same generator matrix.  C3
observed 12,370 nonzero witnesses where the exact polynomial evaluated
nonzero on the omitted suffix; the same 12,370 vectors distinguished the
parent and exact coefficient results.  Any use of this construction belongs
to the separately versioned C7/C8 exact-profile work and still requires its
own systematic/MDS/wire proof.

The Newton path computes divided differences at points shift xor i, accumulates
precomputed Newton-basis polynomials into monomials, and solves the triangular
monomial-to-LCH conversion.  Denominator inverses, Newton-basis polynomials,
and LCH leading-coefficient inverses are all plan setup, not execution.  The
generic path independently precomputes Lagrange cardinal polynomials and
performs the same explicit triangular basis conversion.  The dense exact path
inverts the evaluation matrix directly.  Direct monomial evaluation verifies
all three.

This scope follows the basis-conversion direction in Nicholas Coxon's
[Fast Transforms over Finite Fields of Characteristic Two](https://arxiv.org/pdf/1807.07785),
but the dense column model is not Coxon's fast recursive algorithm.  Tang and
Han's arbitrary-epsilon inverse in
[New Decoding of Reed-Solomon Codes Based on FFT and Modular Approach](https://arxiv.org/pdf/2207.11079)
was not transcribed in this bounded checkpoint.  The profile warning in
`docs/leopard2_math_and_sources.md` continues to apply.

## Correctness evidence

The deterministic analysis used 32 worker processes on the 32 CPUs allowed to
this checkout.  A one-worker rerun produced the identical result file.
The `verify` command independently recomputes that one-worker result and
compares the complete object, so a valid source hash cannot conceal edited
summary, model, or disposition evidence.

| Evidence | Count |
|---|---:|
| Correctness jobs | 76 |
| Input vectors | 14,327 |
| Field/result comparisons | 403,444 |
| Parent-versus-exact difference witnesses | 12,370 |
| Exact suffix-nonzero witnesses | 12,370 |
| Modeled byte/batch/reuse cells | 3,780 |

Coverage includes:

- GF(2^4) q=1 through 16, including every input vector for q at most 3;
- legacy-coordinate GF8 q=1,2,3,4,5,7,8,9,15,16,17,31;
- legacy-coordinate GF16 q=3,5,9,17;
- zero, middle, and final aligned shifted cosets where distinct;
- basis vectors, structured vectors, and stable SHA-256-derived random seeds;
- direct polynomial evaluation of both the padded parent and exact result;
- equality of padded butterflies and the independent dense parent columns;
- equality of dense, Newton, and generic-Lagrange exact coefficients;
- equality of parent and exact inverses when q is the full power of two.

GF4/GF8 multiplication uses independently generated coordinate-to-polynomial
tables.  GF16 uses the legacy Cantor basis and polynomial 0x1002d but performs
carryless multiplication and reduction directly rather than calling Leopard.

The committed deterministic artifact is intentionally retained in full.  Its
981,681 bytes contain every operation row and all 3,780 setup/reuse cells, so a
future dispatcher study can audit the result without rerunning this process.
It is evidence, not a generated production table.

## Setup and execution model

The deterministic score is an explicitly heuristic comparison:

    3 * fixed_multiplications + XORs + loads + 2 * stores
      + 12 * inversions.

It is evaluated separately for setup and per-symbol execution, then amortized
over GF8/GF16 symbol counts, shard bytes 64/1,024/65,536, batches 1/8/64, and
reuse 1/8/64/1,024.  Counts through parent 32 come from instrumenting the
plan-specific scalar schedule, including evaluation-matrix construction and
inversion.  Shared immutable field setup (normalized LCH bases and the skew
table) is excluded from every method.  Larger dense/Newton/Lagrange rows use
labeled algebraic upper models.  This score is useful for rejecting broad
claims, not for predicting a SIMD codec.

The best wire-compatible modeled cell was GF8 q=3, N=4, 65,536-byte shards,
batch 64, reuse 1,024, at 1.032258 times the padded score.  It misses the 10%
promotion threshold.  The best exact cell was GF8 q=5, N=8 at 1.529412 times
the padded score, but it changes the code definition.

Storage also differs: the parent dense executor retains N*q field constants,
the exact dense executor q*q, and the Newton plan retains q*(q-1)/2
denominator inverses plus q*(q+1)/2 Newton-basis coefficients and q LCH
leading inverses.  The padded baseline uses one shared field skew table.  Dense
setup uses Gaussian inversion and grows cubically.

## Pinned scalar timing

Timing was intentionally pinned to allowed CPU 16 with `OMP_NUM_THREADS=1`.
The host is an AMD Ryzen 9 9950X3D with 32 logical/16 physical cores and one
NUMA node.  Python 3.12.3 ran on Linux 6.8.0-134.  Each number is median
nanoseconds for one byte/symbol position across q input shards; it is not shard
throughput or a C++/SIMD performance claim.

| Field | q/N | Padded | Parent dense | Exact dense | Newton | Generic Lagrange | Parent setup |
|---|---:|---:|---:|---:|---:|---:|---:|
| GF8 | 3/4 | 2,144 | 1,159 | 984 | 2,914 | 2,439 | 19,140 |
| GF8 | 5/8 | 5,067 | 2,742 | 1,860 | 5,067 | 4,283 | 67,391 |
| GF8 | 9/16 | 11,277 | 7,630 | 4,484 | 11,725 | 10,258 | 357,994 |
| GF8 | 17/32 | 24,355 | 24,770 | 13,362 | 32,788 | 30,482 | 2,314,006 |
| GF16 | 3/4 | 4,230 | 3,346 | 5,237 | 15,246 | 17,486 | 156,902 |
| GF16 | 5/8 | 13,210 | 17,076 | 15,030 | 40,033 | 39,317 | 781,958 |
| GF16 | 9/16 | 39,954 | 89,112 | 55,202 | 135,183 | 135,541 | 5,721,014 |

The GF8 tiny-q parent result explains the cautious inconclusive disposition,
but production Leopard butterflies use mature vector kernels while this Python
dense loop uses table-backed scalar multiplication.  GF16 reverses the result
by q=5.  The evidence cannot justify a production dispatcher.

## Reproduction

From the repository root:

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3/study.py self-test

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev experiments/leopard2/non_power_of_two/c3/study.py \
      analyze --workers 32 \
      --output experiments/leopard2/non_power_of_two/c3/result.json

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3/study.py verify \
      experiments/leopard2/non_power_of_two/c3/result.json

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev experiments/leopard2/non_power_of_two/c3/study.py \
      analyze --workers 1 --output /tmp/c3-result-repeat.json

    cmp experiments/leopard2/non_power_of_two/c3/result.json \
      /tmp/c3-result-repeat.json

For timing, replace CPU 16 only if it is not in the allowed affinity set:

    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
      PYTHONMALLOC=debug taskset -c 16 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3/study.py benchmark \
      --field gf8 --counts 3 5 9 17 --samples 11 --iterations 64 \
      --output experiments/leopard2/non_power_of_two/c3/timing_gf8.json

    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
      PYTHONMALLOC=debug taskset -c 16 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3/study.py benchmark \
      --field gf16 --counts 3 5 9 --samples 7 --iterations 16 \
      --output experiments/leopard2/non_power_of_two/c3/timing_gf16.json

Exact SHA-256 values for this checkpoint are:

| File | SHA-256 |
|---|---|
| `study.py` | `c0030150633f133d0b04dda762096e3d6164173e99a60888b214476cb2d1b795` |
| `result.json` | `d2c602cd3ee8f704b21cf181d822be81e66356b09552f9d03fe2f72585a4b43d` |
| `timing_gf8.json` | `832fa4da05fec3f578b61b3423a9bbecdd6dca1f2d22f1273fe77635fb81b322` |
| `timing_gf16.json` | `6c9322ca1c987c352408ccd6eb7875a8b47587ac10fc829fe72586778cbff767` |

## Disposition and remaining boundary

- **Promoted:** nothing.
- **Wire-compatible dense parent columns:** inconclusive only for tiny q;
  route any production attempt through the existing direct-path and whole-codec
  SIMD gates.  The broad deterministic model misses the 10% threshold.
- **Exact dense/Newton/Lagrange:** mathematically correct but a new code
  definition.  Retain as C7/C8 input; do not decode or serialize it as a legacy
  profile.
- **Transparent Newton and generic executors:** killed as scalar production
  candidates in the measured cells, retained as independent oracles.
- **Fast Coxon/Tang-Han inverse:** not implemented and still inconclusive.

A production promotion would still require a derived fast conversion,
exhaustive exact-profile MDS proof, serialized identity, C++/SIMD implementation,
whole-codec setup/execution benchmarks, neighboring-regression checks, and a
deterministic dispatcher.  C3 deliberately stops before that integration.
