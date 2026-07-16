# Leopard2 C3b: executable fast inverse candidates

## Result

C3b implements the previously missing published algorithms, but promotes no
production code and no wire profile.

The default build is unchanged.  The implementation and its result files live
only under `experiments/leopard2/non_power_of_two/c3b/`; CMake does not compile,
install, or dispatch to them.

The principal result is negative:

- Coxon's Algorithm 3 Lagrange-to-LCH executor, Coxon's Algorithm 1
  Newton-to-LCH executor, and Tang--Han Appendix-A Algorithm 8 epsilon-point
  inverse all agree with independent exact interpolation, including arbitrary
  prefix lengths and non-aligned shifted cosets.
- None of those three scalar executors reaches the required 10% measured
  improvement over the padded LCH inverse.  Their best pinned gains are 0.828,
  0.812, and 0.862 respectively; a value below one is slower.
- A deterministic operation model predicts a 1.158 best gain for Coxon
  Algorithm 3, but the actual Python recursion overhead reverses that result.
  The measured result governs disposition.
- The small dense exact inverse reaches a 2.219 best pinned gain and exceeds
  10% in 143 of 288 amortization cells.  It remains input to C7/C8's direct
  small-K dispatcher study, not a production promotion: it changes the code
  definition, has quadratic execution storage/work, and has no serialized
  exact-profile identity or end-to-end MDS/codec result yet.
- Generic interpolation reaches only 1.066 and misses the 10% rule.

Thus the fast scalar candidates are rejected as executors for the current
production path and retained as executable algebra/oracles.  The dense exact
method remains a bounded new-profile candidate.  No result is described as
legacy compatible.

## Algorithms implemented

The experiment is a clean implementation from the papers rather than copied
research code.

### Coxon reduction tree

For a subspace basis `beta=(beta_0,...,beta_(n-1))`, the tree builder chooses
the largest power-of-two split `d<n` for which

    (beta_i / beta_0)^(2^d) = beta_i / beta_0,  0 <= i < d.

It constructs Coxon's child bases

    alpha = beta[0:d]
    delta_i = (beta_(d+i) / beta_0)^(2^d)
              + beta_(d+i) / beta_0.

This matters for the test-only GF4 polynomial basis, where assuming that every
child is a Cantor prefix would be wrong.  The legacy GF8/GF16 coordinate bases
admit the larger Cantor reductions.  The plan precomputes Coxon's `phi` values
for the requested shift and the `phi(sigma_(v,i))` update rows.

`coxon_l2x` is Algorithm 3, including all mixed `(c, ell, b)` states and its
strided row/column subproblems.  The top-level exact inverse calls it with
`c=ell=K` and `b=0`; only the first K evaluation values are inputs, and the
result has LCH coefficients K through N-1 fixed to zero.

`coxon_n2x` is Algorithm 1.  Input evaluations are first converted to
coefficients of Coxon's *normalized* Newton basis.  Ordinary divided
differences alone are not those coefficients: divided difference `d_i` is
multiplied by

    product_(j<i) ((shift + omega_i) + (shift + omega_j))

before Algorithm 1.  In characteristic two the shift cancels from this
normalizer.  Direct polynomial evaluation independently verifies the result.

Source: Nicholas Coxon, *Fast Transforms over Finite Fields of Characteristic
Two*, definitions and Lemma 2.1, Algorithms 1 and 3, Theorems 3.3 and 3.10:
https://arxiv.org/abs/1807.07785 and
https://arxiv.org/pdf/1807.07785 .  The paper was retrieved on 2026-07-16.

### Tang--Han epsilon inverse

`tang_han_inverse` implements Appendix-A Algorithm 8.  For `K<=N/2` it
recurses on the known lower evaluations, sets the upper coefficients to zero,
and evaluates the completed upper half.  For `K>N/2` it:

1. inverse-transforms the complete lower half to `w`;
2. evaluates `w` on the upper coset;
3. XORs that evaluation into the known upper prefix;
4. recursively solves the upper epsilon subproblem; and
5. recovers the lower coefficients with the normalized subspace-polynomial
   multiplier.

The implementation returns both the K-coefficient polynomial and Algorithm
8's completed N-evaluation vector.  The supplied prefix must remain unchanged,
and a direct full forward transform must reproduce the completion.

Source: Nianqi Tang and Yunghsiang S. Han, *New Decoding of Reed-Solomon Codes
Based on FFT and Modular Approach*, Appendix A, Lemmas 8--9 and Algorithm 8:
https://arxiv.org/abs/2207.11079 and
https://arxiv.org/pdf/2207.11079 .  The paper was retrieved on 2026-07-16.

## Exact-prefix versus parent-wire mathematics

Let `N=ceil_pow2(K)`, let `B_N(shift)` be the N-point matrix that evaluates
Leopard's normalized LCH basis at `shift + omega_i`, and let `y` contain the K
given values.

The existing parent-wire inverse is

    c_parent = inverse(B_N(shift)) * [y, 0_(N-K)].

It fixes missing *evaluations* to zero and generally returns N nonzero
coefficients.  C3b's Coxon, Tang--Han, dense, and generic methods instead solve

    B_K(shift) * c_exact = y,
    c_exact[K:N] = 0.

They interpolate the unique degree-less-than-K polynomial, and Algorithm 8
also computes its omitted evaluations.  For non-power-of-two K this is a new
exact profile.  At K=N, both definitions reduce to the same full LCH inverse;
all exact methods are required to match the parent byte-for-byte in that case.

This is why an exact algorithm can be faster and still be unusable under a
legacy wire identifier.  C3b records the method class in every deterministic
and timing cell.

## Correctness evidence

The deterministic checkpoint contains:

| Evidence | Count |
| --- | ---: |
| Correctness jobs | 165 |
| Input vectors | 24,605 |
| Algebra/result comparisons | 1,170,240 |
| Parent-versus-exact difference witnesses | 21,151 |
| Exact omitted-suffix nonzero witnesses | 21,151 |
| Setup/execution operation-model cells | 2,376 |

Coverage includes:

- GF(2^4) K=1 through 16, with all 16^K input vectors for K at most 3;
- legacy-coordinate GF8 K=1,2,3,4,5,7,8,9,15,16,17,31;
- legacy-coordinate GF16 K=1,3,5,9,17;
- shifts 0, 1, 3, a non-aligned half-field-plus-one shift, and the final field
  element;
- zero, all-one, every basis input, structured data, and stable SHA-256-seeded
  random vectors;
- equality of dense, generic, Coxon Algorithm 3, Coxon Algorithm 1, and
  Tang--Han Algorithm 8 exact coefficients;
- direct monomial evaluation at every supplied point;
- equality of Algorithm 8's completed vector with a direct full transform;
- direct evaluation of the parent inverse at all N points; and
- equality of exact and parent results for every full-parent case.

The field oracle maps explicit public coordinates to a polynomial basis and
performs carryless multiply/reduce.  GF8 uses polynomial `0x11d` and GF16 uses
`0x1002d` with the legacy Cantor-coordinate maps.  It does not call Leopard's
tables or transforms.  The general Coxon child-basis construction independently
exercises the non-Cantor GF4 basis.

## Setup and execution evidence

Plan-dependent tree, phi, divided-difference, dense-matrix, cardinal, and
epsilon metadata construction is measured as setup.  This is a Python scalar
experiment: execution allocates ordinary lists and tuples, and that interpreter
and allocation cost is included in the reported timing.  It is not an
allocation-free production-kernel claim.  Normalized subspace constants are
cached as the field-level analogue of Leopard's immutable skew table and
excluded equally from every plan.

`TangHanPlan.split_multipliers` retains the planned epsilon split constants for
validation, but this scalar executor does not consume that tuple: execution
looks up the equivalent value in the warmed field-level
`subspace_multiplier` cache.  The recorded Tang--Han plan setup is therefore
conservative and the execution includes the cache lookup; neither number is
evidence for a production split schedule.

The deterministic model uses

    score = 3*multiplications + XORs + loads + 2*stores + 12*inversions.

It applies per-symbol execution to shard bytes 64, 1,024, and 65,536; batches
1, 8, and 64; and reuse counts 1, 8, 64, and 1,024.  Its operation-count result
is useful for explaining algorithms, but not for promotion.

Pinned Python timing used 15 samples and 64 iterations per sample.  Entries
below are median nanoseconds for one symbol position spanning K shards.  MADs
are retained in `results/timing.json`; the largest execution MAD is 0.8383% of
its median.

| Field | K/N | Padded | Dense | Generic | Coxon L2X | Coxon N2X | Tang--Han |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GF8 | 3/4 | 2,094 | 1,266 | 2,633 | 3,243 | 2,863 | 2,702 |
| GF8 | 5/8 | 4,977 | 2,243 | 4,670 | 8,093 | 6,133 | 6,273 |
| GF8 | 9/16 | 11,137 | 5,464 | 11,962 | 21,474 | 15,100 | 14,173 |
| GF8 | 17/32 | 24,789 | 16,620 | 37,049 | 43,375 | 36,203 | 31,399 |
| GF16 | 3/4 | 4,029 | 5,446 | 16,396 | 5,071 | 8,055 | 4,674 |
| GF16 | 5/8 | 13,271 | 14,012 | 33,784 | 16,620 | 23,341 | 16,924 |
| GF16 | 9/16 | 35,426 | 48,983 | 127,897 | 45,557 | 85,311 | 49,056 |
| GF16 | 17/32 | 90,168 | 217,698 | 490,988 | 108,916 | 315,901 | 131,504 |

Plan-setup medians range from 11.4 microseconds for a GF8 K=3 generic plan to
8.31 milliseconds for a GF16 K=17 Coxon Newton plan.  The Lagrange-only Coxon
plan avoids the divided-difference constants (11.9 microseconds through 1.06
milliseconds).  Warm field-table Tang--Han plan creation is below 1.2
microseconds in these cells, but execution still loses.

An immediate second complete pinned run differed by at most 1.98% in any
execution median.  It is retained as an ignored `/tmp` repeat rather than
silently averaged into the committed observation.

`results/timing_summary.json` separately amortizes setup over the same byte,
batch, and reuse axes, producing 1,728 cells.  Those totals are extrapolations
from scalar-position timing, not claims of actual shard-memory throughput.

The host exposed 32 logical CPUs (16 physical cores, one NUMA node), not the
mission's requested 128.  Correctness analysis is capped at the assigned eight
workers.  Timing intentionally uses one process pinned to logical CPU 15 with
`OMP_NUM_THREADS=1`; using eight workers for a cache-sensitive latency result
would invalidate it.  The host is an AMD Ryzen 9 9950X3D and reported the
`powersave` governor.  These Python timings are research evidence, not a claim
about mature C++ SIMD kernels.

## Disposition under the 10% rule

| Candidate | Best measured gain over padded | Disposition |
| --- | ---: | --- |
| Coxon Algorithm 3 Lagrange-to-LCH | 0.828 | reject scalar executor; retain exact-profile algebra |
| Coxon Algorithm 1 Newton-to-LCH | 0.812 | reject scalar executor/oracle path |
| Tang--Han Algorithm 8 epsilon IFFT | 0.862 | reject general scalar executor; retain completion oracle |
| Generic interpolation | 1.066 | below threshold; retain oracle |
| Dense exact inverse | 2.219 | retain for bounded C7/C8 direct-path study only |

No wire-compatible candidate exists in this checkpoint: the exact methods
only match the parent code at K=N.  Even the dense candidate cannot be promoted
until an exact profile has a versioned coordinate map and serialized identity,
exhaustive MDS proof, whole-codec encode/decode integration, C++/SIMD evidence,
and a dispatcher that avoids its quadratic crossover.

## Reproduction

From the repository root:

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3b/fast_inverse.py self-test

    PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 PYTHONMALLOC=debug \
      python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3b/fast_inverse.py \
      analyze --workers 8 \
      --output experiments/leopard2/non_power_of_two/c3b/results/checkpoint.json

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3b/fast_inverse.py verify \
      experiments/leopard2/non_power_of_two/c3b/results/checkpoint.json

Use another allowed physical core if CPU 15 is unavailable:

    OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONHASHSEED=0 \
      PYTHONMALLOC=debug taskset -c 15 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3b/fast_inverse.py benchmark \
      --fields gf8 gf16 --counts 3 5 9 17 --shift 3 \
      --samples 15 --iterations 64 \
      --output experiments/leopard2/non_power_of_two/c3b/results/timing.json

    PYTHONDONTWRITEBYTECODE=1 python3 -B -X dev \
      experiments/leopard2/non_power_of_two/c3b/fast_inverse.py \
      summarize-timing \
      experiments/leopard2/non_power_of_two/c3b/results/timing.json \
      --output \
      experiments/leopard2/non_power_of_two/c3b/results/timing_summary.json

The committed JSON files retain every correctness digest, operation counter,
setup/execution median and MAD, and modeled or timed amortization cell.  Their
source hash makes stale evidence fail verification.
`summarize-timing` rejects a stale source hash and binds the raw timing file by
SHA-256, but there is no separate summary-verification subcommand; reproduce
the summarization command above and compare the committed hash when auditing
that derived file.

SHA-256 values for this checkpoint are:

| File | SHA-256 |
| --- | --- |
| `fast_inverse.py` | `bf63d8f228c7ca19b49f1971affc72ebfeebda870111111662c8b32ede0c9e84` |
| `results/checkpoint.json` | `5e71fa9ea2e72d9fec10f4e123a1acf85bd3de2c74738f523d59b80e1660e0bc` |
| `results/timing.json` | `018841e30afc5a646a1c8639ef56b39ad87112fa975957099140154d0d8d1ca3` |
| `results/timing_summary.json` | `9324311521e9c501b2caee96e1381cc37357d415c0fca93089a5f2e0970c191b` |
