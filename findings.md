# Findings: exact GFNI/LCH kernels, batching, setup, and application fusion

Evidence snapshot: 2026-08-19  
Writeup revision: 2026-08-26

## Executive finding

This episode began with a substantial implementation improvement inside the
known LCH additive-transform family. For the fixed 256-point transform over
`GF(2^8)`/`0x11b`, the forward kernel holds the complete state in four AVX-512
registers and fuses all eight stages; an independently searched masked-ternary
inverse later reduced the integrated inverse body to 860 bytes and 135 decoded
instructions. Their exact warm, direct-I/O table comparisons reached 34.7630x
forward and 46.2686x inverse on the recorded Threadripper PRO 9985WX.

The bounded follow-up did not stop at N=256. Frozen N=16/32/64/128 kernels,
shared-gamma batch-two and terminal batch-four kernels, and a separately
contracted fused degree-at-most-127 polynomial product all survived fresh
fixed-pair confirmations. The compact constructor was reduced to 84 products
for arbitrary bases, while the standard-basis time-memory frontier retained a
charged 2,048-byte gamma image. Under the final setup-plus-one-transform
contract, constructing the flat plan versus this compact plan measured
166.4703x forward and 170.4831x inverse; this large setup ratio is distinct from
the hot-transform ratios.

The investigation also preserved decisive negatives: setup-charged nonzero
offset transport lost, forward batch-four policy selection did not reproduce,
and bounded compiler, packed-upper, automorphism, tiny-transform, batch-eight,
tower-bitslice, synthesized-XOR, Frobenius, Cantor, and multiplicative-DFT
paths failed exact gates. The canonical ledger contains 41 branches, 171
hash-bound evidence files, and no pending indexed local experiment in that
snapshot.

None of these ratios may be stacked or generalized across contracts. They are
host-local repository implementation results with `claim_eligible=false`, not
proof of a new transform, portable speed, asymptotic improvement, mathematical
novelty, or a win over an independently qualified strongest LCH baseline.

### How to read this report

Every speed ratio is **comparator divided by candidate**. A value above one
therefore favors the candidate. Results are grouped by exact workload,
`freedom_mode`, setup policy, batch size, thread count, CPU, and comparator;
ratios from different rows are not multiplicative. A `confirmed-fixed-contract`
result means that a frozen selection survived a fresh paired campaign and all
preregistered correctness, order, layout, confidence-interval, and harness-floor
gates. A `diagnostic-selection-only` result chose what to confirm next but is
not itself the final performance evidence. A `closed-negative` result failed a
declared gate and is retained rather than rerun.

The machine-readable source of truth is
[`authoritative-index.json`](artifacts/lch-postfrontier/authoritative-index.json),
SHA-256
`c51b93ac0263637b8151815c4fb3299594232046510334c05b5a90cd2e32a50c`.
It indexes 41 material branches and 171 unique evidence files: 10 confirmed
fixed contracts, 3 correctness-only integrations, 3 diagnostic selections, 16
bounded negatives, 5 superseded-but-valid lineages, 3 zero-measurement failed
prefixes, and 1 external-authority boundary. No local branch is pending in that
snapshot.

### Confirmed result matrix

All intervals below are hierarchical 95% intervals over retained paired
measurements. `Setup excluded` rows are hot steady-state contracts; the compact
row includes plan construction and must not be compared numerically with them.

| Workload | Freedom/setup | Comparator/candidate | 95% interval | Evidence |
|---|---|---:|---:|---|
| N=256 forward direct I/O | `STRICT_ISOMORPHIC`, setup excluded | 34.763022x | [34.636236, 34.849303] | [forward bundle](artifacts/lch-direct-frontier/authenticated-forward-confirmation-v5-20x31-bootstrap50000/report.json) |
| N=256 inverse policy 40 | `STRICT_ISOMORPHIC`, setup excluded | 46.268605x | [46.137977, 46.376552] | [policy-40 report](artifacts/lch-postfrontier/authenticated-inverse-policy40-confirmation-v4-20x31-bootstrap50000/report.json) |
| N=16 forward / inverse | `STRICT_ISOMORPHIC`, setup excluded | 2.510703x / 2.093623x | [2.508568, 2.513653] / [2.024460, 2.097091] | [N=16/32 index](artifacts/lch-neighbor-frontier/small16-32-masked-fixed-pair-confirmation-20x31-cpu0-v1/index.json) |
| N=32 forward / inverse | `STRICT_ISOMORPHIC`, setup excluded | 16.233668x / 5.537158x | [16.126465, 16.303935] / [5.531838, 5.546492] | [N=16/32 index](artifacts/lch-neighbor-frontier/small16-32-masked-fixed-pair-confirmation-20x31-cpu0-v1/index.json) |
| N=64 forward / inverse | `STRICT_ISOMORPHIC`, setup excluded | 28.251709x / 18.235771x | [28.183304, 28.482910] / [18.211458, 18.243526] | [N=64/128 index](artifacts/lch-neighbor-frontier/neighbor64-128-id15-fixed-pair-confirmation-20x31-cpu0-v2/index.json) |
| N=128 forward / inverse | `STRICT_ISOMORPHIC`, setup excluded | 32.880183x / 23.257695x | [32.631987, 33.232124] / [23.241886, 23.275259] | [N=64/128 index](artifacts/lch-neighbor-frontier/neighbor64-128-id15-fixed-pair-confirmation-20x31-cpu0-v2/index.json) |
| Batch two forward / inverse | `STRICT_ISOMORPHIC`, setup excluded | 1.255705x / 1.185693x | [1.253849, 1.258319] / [1.183598, 1.188872] | [forward](artifacts/lch-postfrontier/batch2-forward-p0-confirmation-ready-v1/runs/one-shot-claimed/success/report.json), [inverse](artifacts/lch-postfrontier/batch2-inverse-p3-confirmation-ready-v1/runs/one-shot-claimed/success/report.json) |
| Batch four forward / inverse versus two batch-two calls | `STRICT_ISOMORPHIC`, setup excluded | 1.181308x / 1.354821x | [1.176661, 1.185979] / [1.340274, 1.360994] | [forward audit](artifacts/lch-postfrontier/batch4-forward-p0-vs-two-batch2-terminal-confirmation-postrun-audit-v2/index.json), [inverse audit](artifacts/lch-postfrontier/batch4-inverse-p3-vs-two-batch2-terminal-confirmation-postrun-audit-v2/index.json) |
| Compact gamma setup plus one transform, forward / inverse | `APPLICATION_CODESIGN`, setup included | 166.470313x / 170.483093x | [165.750000, 169.839228] / [169.998392, 170.677419] | [clean-room audit](artifacts/lch-postfrontier/compact-gamma-direct-setup-confirmation-v2-postrun-audit-v1/index.json) |
| Degree-at-most-127 polynomial product | `APPLICATION_CODESIGN`, setup excluded-fixed | 1.136470x | [1.135113, 1.137685] | [independent audit](artifacts/lch-postfrontier/degree127-polynomial-product-policy0-confirmation-postrun-audit-v1/index.json) |

Performance reproduction is deliberately **retained-only and no-fresh-timing**.
Fresh correctness builds are allowed, but the one-shot timing campaigns are
sealed. Section 11 shows how to rebuild correctness checks and recompute
hashes, schedules, checksums, statistics, and bootstrap intervals without
collecting new timings.

## 1. Initial exhaustive-census workload

The workload is one forward transform with the following fixed semantics:

| Property | Value |
|---|---|
| Field | `GF(2^8) = GF(2)[x]/(x^8+x^4+x^3+x+1)` (`0x11b`) |
| Input | 256 coefficients in the repository's LCH novel-polynomial basis |
| Output | Evaluations on the 256-point additive domain |
| Additive basis | Byte values `1, 2, 4, ..., 128` |
| Domain offset | `0` |
| Ordering | Additive-mask order: output index bits select basis vectors |
| Direction | Forward only |
| Threads | One |
| Performance mode | Warm-cache latency |
| Setup | Fixed precomputation excluded from timing but counted |
| Precomputation cap | 131,072 bytes |
| Workload hash | `4ef2ea225f5b8e5756c21b5d2b79de9d06029ec47be77cc5f88fb87b1b37799e` |

Field addition is bytewise XOR. The scalar field definition and the direct
novel-basis evaluation are implemented in
[`gf256.hpp`](cpp/include/addfft/gf256.hpp) and
[`lch.hpp`](cpp/include/addfft/lch.hpp). The external workload descriptor is
constructed in [`cli.py`](addfft_lab/cli.py), while the matched timing loop is
in [`benchmark.cpp`](cpp/src/benchmark.cpp).

The performance comparator is not the original recursive implementation. It is
an allocation-free flat plan using the same fused butterfly as the candidate,
the same input buffers, the same output-copy convention, and a full 256-by-256
multiplication table. This matters: the approximately `16x` ratio is not merely
the removal of recursive allocations or the algebraic butterfly repair. It is
the measured difference between the flat table path and the specialized GFNI
path under this contract.

## 2. The algebraic reduction

At one forward LCH butterfly, let `a` be the lower value, `b` the upper value,
and `gamma` the stage-and-block constant. The conventional outputs are

```text
L = a + gamma*b
R = a + (gamma+1)*b
```

Distributivity gives

```text
(gamma+1)*b = gamma*b + b,
```

and therefore

```text
P = gamma*b
L = a + P
R = L + b.
```

Because the field has characteristic two, both additions are XORs. The
butterfly consequently needs one field multiplication and two XORs rather than
two field multiplications and two XORs; when `gamma` is zero, the multiplication
can also be skipped. The portable implementation is the three-line inner loop
in [`flat_lch.hpp`](cpp/include/addfft/flat_lch.hpp).

This identity was checked by executable exhaustion over all
`256^3 = 16,777,216` triples `(gamma, a, b)`. Each multiplication backend was
also checked over all `256^2 = 65,536` operand pairs. It is an exact arithmetic
reduction for the known butterfly, not a claim of new mathematics.

## 3. Register-resident GFNI realization

The winning kernel maps one 256-byte state to four 64-byte ZMM registers. Its
implementation is in [`gfni_lch.cpp`](cpp/src/gfni_lch.cpp).

For levels 0 through 5, each register contains complete butterfly pairs. A
`VPERMB` permutation maps lane `i` to its partner `i xor (1 << level)`, masks
identify high lanes, and `VGF2P8MULB` performs bytewise multiplication in the
`0x11b` field. In simplified form, a stage is

```text
partner = permute(values, i xor (1 << level))
upper   = select_upper(values, partner)
product = gfni_multiply(upper, gamma)
low     = values xor product
high    = low xor partner
result  = blend(low, high)
```

Level 6 couples register pairs `(0,1)` and `(2,3)`. Level 7 couples `(0,2)` and
`(1,3)`. The implementation instantiates all 36 possible contiguous stage
ranges at compile time. A `forward_fixed_range<Highest,Lowest>` specialization
loads four registers once, inlines the requested descending stages, and stores
the four registers once.

The selected `7-to-0` body has the following retained disassembly properties:

- 1,245-byte hot body;
- no calls, stack frame, stack references, or ZMM spills;
- 24 `VPERMB` sites and 28 static GFNI-multiply sites;
- 25 GFNI multiplies on the measured offset-zero path because three high-stage
  zero-gamma cases skip their multiplication;
- four input vector loads and four output vector stores.

The runtime dispatcher requires GFNI plus AVX-512F, AVX-512BW, and AVX-512VBMI
and caches the feature test. A portable table fallback remains available, but
the measured winning path is ISA-specific. The retained disassembly is
[`confirmation-selected-range7-0-disassembly.txt`](artifacts/lch-schedule-search/confirmation-selected-range7-0-disassembly.txt),
SHA-256 `6b42ddd529677b993a2e7feb8a5f0f2d6bc33f39c854104cd8cbae0fcd203113`.

## 4. What was exhaustively searched

The search grammar deliberately froze four implementations for each of eight
stages:

| Digit | Backend |
|---:|---|
| 0 | Fused scalar multiplication |
| 1 | Fused full-table multiplication |
| 2 | Fused two-nibble-table multiplication |
| 3 | Fused GFNI/AVX-512 stage |

A policy is the resulting eight-digit base-4 number, so the census contains
exactly

```text
4^8 = 65,536 policies.
```

There are seven potential boundaries between adjacent stages. A fusion bit is
effective only when the stages on both sides use GFNI. If a policy has `a`
adjacent GFNI/GFNI boundaries, it has `2^a` distinct effective fusion variants.
Summing that count over every policy gives

```text
sum over policies of 2^a = 107,621 canonical variants.
```

That total has a short exact derivation. Let `N_n` and `G_n` be the weighted
counts of length-`n` backend strings ending in a non-GFNI stage and a GFNI
stage, where a string with `a` GFNI/GFNI adjacencies has weight `2^a`. For one
stage, `(N_1, G_1) = (3, 1)`. Appending one of the three non-GFNI backends or
one GFNI backend gives

```text
N_(n+1) = 3*(N_n + G_n)
G_(n+1) = N_n + 2*G_n.
```

The factor two in the second recurrence is the new effective fusion choice
when GFNI follows GFNI. Iterating to eight stages gives
`(N_8, G_8) = (75,036, 32,585)`, whose sum is `107,621`.

Each canonical variant represents the raw masks that differ only at ineffective
boundaries. For a policy with `a` effective boundaries, each canonical mask has
exactly `2^(7-a)` raw representatives because the other `7-a` bits are free
no-ops. The `2^a` canonical masks therefore represent
`2^a * 2^(7-a) = 128` raw masks for every policy, so the complete census
quotient-covers

```text
65,536 policies * 128 masks = 8,388,608 raw encodings.
```

The 8,388,608 raw encodings were not all timed separately; redundant encodings
were covered by an exact no-op equivalence. Every canonical variant was
screened. The best 256 were retimed, and the declared deterministic winner rule
ordered them by median latency, median absolute deviation, policy ID, and fusion
mask. The selected finalist was policy `65535`, mask `127`, with a search-retime
median of `27.3828125 ns`. Screening selected the implementation; it was not
used as the final performance claim.

The executable census is in
[`schedule_search.cpp`](cpp/src/schedule_search.cpp). Its authenticated runner
and validator are in
[`certify_lch_schedule_search.py`](tools/certify_lch_schedule_search.py).

## 5. Why the finite correctness argument covers every input

Correctness was established in layers before performance ranking:

1. **Field primitives.** Scalar, full-table, nibble-table, and GFNI
   multiplication agreed for all 65,536 input pairs.
2. **Butterfly identity.** The one-multiplication butterfly agreed with the
   conventional two-multiplication form for all 16,777,216 triples.
3. **Complete coordinate checks.** Every one of the 65,536 policies and every
   one of the 107,621 canonical variants was checked on all 256 coordinate unit
   vectors. This produced 16,777,216 policy-coordinate cases and 27,550,976
   variant-coordinate cases.
4. **Independent oracle witnesses.** The C++ executable emitted its complete
   256-column by 256-output transform matrix and all eight levels of gamma
   lanes. The certifier compared those bytes with the independent Python LCH
   oracle.
5. **Executable and source binding.** The certifier launched a receipt-bound
   search binary, checked its receipt before and after execution, recomputed
   timing summaries and finalist selection, and bound the output to frozen
   source hashes.

The coordinate test is complete, rather than a random sample, because every
admitted path is a fixed-gamma `GF(256)`-linear map. The schedules, lane
permutations, masks, and control paths do not depend on input data. Equality on
the 256 coordinate basis vectors therefore implies equality for every vector in
the 256-dimensional `GF(256)` input space, once multiplication equivalence has
been exhausted.

The two principal oracle witnesses are:

| Witness | Size | SHA-256 |
|---|---:|---|
| External transform matrix | 65,536 bytes | `5d916b1a117baaa62fa72047ca7dd101cee68aefa94c83a8d05be1a002843376` |
| Gamma lanes, levels 0–7 | 2,048 bytes | `74eaf720aad6545aa4f24345f41d27ab78e3f8121b1265995cb2d4194c8b3935` |

The final census contains these authenticated counts and witnesses in
[`authenticated-register-v3c-seed-20260817-manifest.json`](artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-manifest.json).
The raw record stream is retained both uncompressed and as
[`authenticated-register-v3c-seed-20260817-results.jsonl.gz`](artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-results.jsonl.gz).

## 6. Benchmark contract and results

Both measured paths use the flat fused plan, the same prepared input buffers,
256 warm-up transforms, and 2,048 transforms per timing sample. Within each
launch, baseline-first and candidate-first pairs are balanced and randomized.
The controller ran 20 independent process launches with 31 pairs each and used
a 50,000-resample hierarchical bootstrap over launches and pairs.

The recorded host was an AMD Ryzen Threadripper PRO 9985WX running Linux
`6.8.0-137-generic`, with GCC 13.3.0 and
`-O3 -DNDEBUG -std=c++20 -Wall -Wextra -Wconversion -Wshadow`. The CPU
fingerprint was
`7a08c75795edbce317aff77368fb93082a0d43acf351c12b6daa815f1b751d79`.

### Confirmation results

| Evidence | Original fresh confirmation | Snapshot-bound repeat |
|---|---:|---:|
| Launches × pairs | 20 × 31 | 20 × 31 |
| Retained/rejected pairs | 620 / 0 | 620 / 0 |
| Baseline median | 431.306 ns | 433.262 ns |
| Candidate median | 26.9865 ns | 27.1725 ns |
| Median paired ratio | 16.0801695× | 15.9864604× |
| Ratio MAD | 0.1603379 | 0.2517607 |
| Hierarchical-bootstrap 95% interval | [15.9819623, 16.1255484] | [15.8055742, 16.1991502] |
| Report SHA-256 | `1de944d56ac6ad024bc9e7f8106625c050652fa20cfe8d8c848711104f949ee5` | `2ae4ddd7c4ead90674add9f8449599a1e5e95c6c7db4f76b3de602accf488f73` |

The snapshot-bound repeat is the strongest long-term provenance chain. Its
build was configured directly from the retained source snapshot, its receipt
verifies against that path, and it produced the same benchmark binary SHA-256
as the original confirmation:
`b7efc937b7025ea6338e5bfdf6d8a0020eeb86d2bca04be842158102609562ac`.
Two of its 20 launch medians were about `12x`; all constituent pairs were kept
and no cause is assigned because hardware counters were unavailable.

The original report is
[`authenticated-confirmation-20x31-bootstrap50000.json`](artifacts/lch-schedule-search/authenticated-confirmation-20x31-bootstrap50000.json).
The snapshot-bound report and its explicit resource/command record are
[`authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000.json`](artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000.json)
and
[`authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000-context.json`](artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000-context.json).

### Resource accounting

| Resource | Baseline | Candidate |
|---|---:|---:|
| Declared precomputation | 75,785 bytes | 76,217 bytes |
| Contract cap | 131,072 bytes | 131,072 bytes |

The candidate's additional 432 bytes are the GFNI partner-index and lane-mask
tables. Specializing every possible contiguous range also creates a code-size
tradeoff:

- selected hot body: 1,245 bytes;
- all 36 specialized bodies: 20,303 bytes;
- kernel object executable text: 26,730 bytes;
- benchmark executable `.text`: 166,893 bytes;
- complete benchmark file: 214,744 bytes.

The selected body is small; the aggregate specialization footprint is the
tradeoff. Application-level instruction-cache stress, cycles, instructions,
branches, cache misses, stalls, and energy were not measured and are recorded
as unavailable rather than zero.

## 7. Other notable findings and repaired failures

### 7.1 A fair flat comparator was necessary

The original recursive forward transform allocated temporary vectors at every
internal recursion node and computed both field products in each butterfly.
Comparing a vectorized candidate only against that path would conflate
allocation removal, arithmetic repair, and SIMD. Building the allocation-free,
one-multiplication flat table path first made the final comparison materially
stronger.

### 7.2 Policy-first fusion pruning was not valid

An early trajectory screened all policies without fusion and explored fusion
only for selected policies. It chose policy `30207`, mask `47`, but fusion
changed the timing order enough that the combined search was not complete. That
trajectory was preserved and rejected; the repaired search enumerated every
canonical effective fusion subset for every policy.

### 7.3 Intrinsics and a cost model hid the real lowering

The first complete GFNI census selected mixed policy `62975`, mask `71`, and a
fresh confirmation showed only `1.1106608x`. Disassembly then revealed a scalar
64-iteration mask-building loop and a `0x140`-byte stack frame that materialized
all four logical ZMM values. Constant masks and compile-time range
specialization removed both defects. This is a useful methodological result:
the source-level intrinsic count was not a reliable description of the emitted
machine path.

### 7.4 The first cross-register shortcut was too narrow

An intermediate GFNI implementation hardcoded zero gamma at level 7 and in the
first level-6 block. That was correct for the measured offset-zero workload but
incorrect for the plan's public nonzero-offset constructor. The repair loads
the actual level-7 row and both level-6 gamma halves, skips only genuinely zero
blocks, and tests offsets `1`, `0x80`, and `0xff` against the recursive
reference. Those regressions are in
[`test_baselines.cpp`](cpp/tests/test_baselines.cpp).

### 7.5 Search certificates needed adversarial validation

Earlier certificate code checked enumeration shape but could accept fabricated
reported medians or an arbitrary finalist set. The repaired certifier rejects
non-finite timings, recomputes every median and MAD from full-precision samples,
reconstructs the exact top-256 selection and winner, validates footer totals,
compares matrix and gamma bytes, and authenticates the producing binary and
source. Tamper tests are retained in
[`test_schedule_search_certificate.py`](tests/test_schedule_search_certificate.py).

### 7.6 Build provenance is not baseline qualification

A valid local build receipt proves which source and binary were measured; it
does not prove that the comparator is an independently qualified Baseline 2.
The benchmark controller was hardened to reject missing, malformed, or
self-asserted qualification. Consequently both final reports correctly contain
`claim_passed: false`, regardless of the size of the measured ratio. The
contract tests are in
[`test_math_benchmark_contract.py`](tests/test_math_benchmark_contract.py).

### 7.7 Source-path binding mattered

The first confirmation's source-manifest bytes matched the retained snapshot,
but its CMake cache named the live repository path. Verifying that receipt
directly against the snapshot path therefore failed closed. Reconfiguring and
sealing from the snapshot produced byte-identical binary output and a receipt
that verifies directly against the retained source. The failed prefix and
repair are preserved in [`experiments.jsonl`](experiments.jsonl).

### 7.8 Artifact names can mislead

The generic `artifacts/lch-schedule-search/results.jsonl` and `manifest.json`
are old intermediate files and are not coherent with each other: the generic
results hash is
`b6b4036f7652df37c3d2f0d6c570689ed6236f577fbb17abc2f470650ff08ef1`,
while the generic manifest expects
`5782a403658b4c47937d5ecf1a0b115b56a9010488a94d6a4c44674ff05f4610`,
the stream retained as `generic-mask-results.jsonl`. The
`partial-policy-first-manifest.json` file likewise retains a stale generic
results path while its hash binds the descriptively named partial result. The
canonical evidence is the schema-v3c census, the original fresh confirmation,
and the snapshot-v4-final repeat listed in
[`artifacts/lch-schedule-search/README.md`](artifacts/lch-schedule-search/README.md).
The v4 repeat has the strongest path-self-contained provenance; it does not
invalidate or replace the original fresh measurement. Intermediate files remain
useful because they preserve false assumptions and repairs rather than
rewriting the episode as a clean success story.

## 8. Direct-I/O forward follow-up

The first census timed an API that copied a 256-byte input before invoking an
in-place kernel. Once the selected transform reached roughly 27 ns, that copy
and generic range dispatch became material. The follow-up therefore defined a
new, hashed host-local contract instead of silently changing the old workload:

| Property | Direct-forward value |
|---|---|
| Semantic parent hash | `4ef2ea225f5b8e5756c21b5d2b79de9d06029ec47be77cc5f88fb87b1b37799e` |
| Direct-I/O contract hash | `ccbc603128e5f5f5d7ce4029fe11300b86c174fab7cf6ff5ee7f370f2a7644f5` |
| Input/output | Immutable, distinct, nonoverlapping 256-byte buffers |
| Alignment | Both addresses congruent to 1 modulo 64 |
| Output | All 256 bytes overwritten by each arm |
| Dispatch | ISA/gamma validation and selection excluded; cached indirect call included |
| Allocation | Fixed stack output storage reserved by the caller before timing |
| Cache/threads/batch | Warm cache, one thread, one transform |
| Comparator | Copy input, then allocation-free full-table forward transform |
| Precomputation | 75,785 comparator bytes; 76,473 candidate bytes |
| Working memory | 2,880 bytes per arm |

The apparently natural direct table comparator was not selected merely because
its source avoided a copy. In the long diagnostic retime it measured 480.634 ns
versus 431.816 ns for copy-then-transform, so the faster copy path was frozen as
the comparator. This reversal is useful evidence that source-level operation
count and real code generation can disagree.

The candidate search replaced general `VPERMB` partner permutations with exact
width-specific instructions: `VPSHUFB` for levels 0/1, `VPSHUFD` for levels
2/3, and `VSHUFI64X2` for levels 4/5. A complete four-switch grammar covered all
`2^4 = 16` choices of those three stage-pair families and byte-versus-dword
upper/merge handling. Every member passed 256 direct-oracle coordinates, all
65,536 coordinate-byte singleton inputs, and 1,024 dense inputs. Lowering 15,
which enables all four specializations, won; a subsequent complete four-member
residual census found that a rotate or ternary-logic finish did not survive a
longer retime.

The first direct screen is deliberately retained as a failed prefix. Its C++
output object was value-initialized inside the timed loop, generating a hidden
320-byte `REP STOS`, and its checksum cancelled to zero. After both defects were
removed, disassembly showed no timed initialization, allocation, input copy, or
runtime arm branch. The fixed candidate still includes one cached indirect call
per transform.

The corrected self-contained 20-launch by 31-pair confirmation produced:

| Forward direct statistic | Value |
|---|---:|
| Valid/rejected pairs | 620 / 0 |
| Comparator median | 419.150146 ns |
| Candidate median | 12.059814 ns |
| Median paired ratio | 34.763022x |
| Ratio MAD | 0.397781 |
| Hierarchical-bootstrap 95% interval | [34.636236, 34.849303] |
| Baseline-first stratum | 309 pairs; 34.742288x median |
| Candidate-first stratum | 311 pairs; 34.787924x median |

The selected leaf is 1,053 bytes and 169 decoded instructions, including 24
`VGF2P8MULB`, eight `VPSHUFB`, eight `VPSHUFD`, and eight `VSHUFI64X2`
instructions. It contains no `VPERMB`, call, branch, stack reference, or spill.
Five 64-byte source operands cover the four input blocks because the second
block is reloaded once; four stores write the output. Operand counts are not a
claim about physical cache traffic.

The authoritative bundle is
[`authenticated-forward-confirmation-v5-20x31-bootstrap50000`](artifacts/lch-direct-frontier/authenticated-forward-confirmation-v5-20x31-bootstrap50000).
Its report SHA-256 is
`a2e6b001d8bffe729c7b72abe3710877ceb28e3c11873adadab509fc36b14e33`;
its context SHA-256 is
`f95ebaf9a5e224249500e5f175773b76b4d869fd4d0a41ea3ad0163f6461b63b`.
The v1 and v2 predecessors retain valid timing samples, but v1 has a misleading
workload-provenance label and no archived exact runner bytes, while v2 did not
bind the actual Bubblewrap path or source-bound test commands. They are
preserved along with v3 and v4, but superseded by v5: v3 bound the absolute
sandbox path and source-bound tests, v4 retained the sealed correctness builds,
and v5 archived the complete Python controller dependency closure.

## 9. Exact direct inverse

The forward butterfly

```text
L = a + gamma*u
R = L + u
```

has determinant one in characteristic two. Its inverse requires no field
division:

```text
u = L + R
a = L + gamma*u.
```

Inverse stages execute in ascending order. The implementation reuses the exact
stage-specific partner shuffles and the plan's gamma table. The known zero
blocks eliminate four of the nominal 28 vector field multiplies, leaving the
same 24 `VGF2P8MULB` instructions as the forward candidate, but the inverse
combine dataflow is shorter.

The inverse has its own semantic and performance identities:

| Property | Direct-inverse value |
|---|---|
| Semantic relation | Return the unique coefficients whose standard offset-zero LCH forward transform equals the input |
| Semantic hash | `b26fb61f370dccbead9b39cba980370766000cec75f01e47057ba7772c9236ee` |
| Direct-I/O contract hash | `7b9e0451adf4f161d940cbacd1b1eaa45a7aa40b4ca1d0a258763909ef1748d5` |
| I/O, alignment, setup, cache | Same direct-I/O discipline as Section 8 |
| Comparator | Direct-first-stage, allocation-free full-table inverse |
| Precomputation/working memory | 75,785 vs 76,473 bytes; 2,880 bytes per arm |

Correctness was gated by 256 independently recursive inverse coordinate cases,
all 65,536 coordinate-byte singleton cases, 1,024 dense forward/inverse round
trips, nonzero-offset compatibility/rejection checks, immutable unaligned input,
poisoned full output overwrite, canaries, and Release plus ASan/UBSan suites.
Because every fixed operation is `GF(256)`-linear, the complete coordinate
basis is also an all-input equivalence certificate; the larger singleton and
dense sets independently stress the implementation and API contract.

The fixed inverse leaf is 926 bytes and 148 decoded instructions. It contains
24 `VGF2P8MULB`, eight each of `VPSHUFB`, `VPSHUFD`, and `VSHUFI64X2`, four
source loads, four output stores, and no `VPERMB`, call, branch, stack reference,
or spill.

The self-contained inverse confirmation measured:

| Inverse direct statistic | Value |
|---|---:|
| Valid/rejected pairs | 620 / 0 |
| Comparator median | 408.360229 ns |
| Candidate median | 10.120239 ns |
| Median paired ratio | 40.435169x |
| Ratio MAD | 0.323327 |
| Hierarchical-bootstrap 95% interval | [40.310665, 40.574944] |
| Baseline-first stratum | 309 pairs; 40.395777x median |
| Candidate-first stratum | 311 pairs; 40.479694x median |

The authoritative bundle is
[`authenticated-inverse-confirmation-v3-20x31-bootstrap50000`](artifacts/lch-direct-frontier/authenticated-inverse-confirmation-v3-20x31-bootstrap50000).
Its report SHA-256 is
`4e965ee2752656487539a275b5de92c688da8e59b1eaeb1cc561ad4ea610f55c`;
its context SHA-256 is
`9ef201661ac9cd95eb227c798881f9be7ebffeeeece1677b5577601ea370f1f2`.

Both direct bundles are locally authority-verified and ran their measured
binary inside Bubblewrap with no home directory, source tree, controller key,
or network available. Each bundle archives the exact C++ source, executable,
receipt, complete trusted Python controller dependency closure, selection
inputs, source-bound correctness evidence and its receipt-bearing build
archive, and instruction-bound disassembly. The HMAC key is intentionally not
archived, so
the receipt is local provenance evidence rather than a portable signature or
independent performance authority.

## 10. Historical direct-I/O endpoint

At the end of the first direct-I/O phase, the supported statement was:

> On the recorded Threadripper PRO 9985WX/GCC 13.3 host, for the exact
> single-thread, warm-cache, 256-point, offset-zero workloads, a receipt-bound
> exhaustive census selected a fully fused all-GFNI LCH implementation. Its
> original forward contract observed about 16x versus the repository's flat
> table path. Under separately hashed immutable/distinct direct-I/O contracts,
> a complete 16-member forward lowering search confirmed 34.7630x, and an exact
> fixed inverse confirmed 40.4352x, each over 620 retained pairs with lower
> hierarchical 95% bounds of 34.6362x and 40.3107x respectively.

The following stronger interpretations are not supported:

- **Not “beats LCH.”** The comparator is an unqualified Baseline 2 candidate,
  not an independently recognized strongest implementation.
- **Not a new transform or asymptotic result.** The arithmetic and transform
  remain in the known LCH family.
- **Not a portability result.** The winning path requires GFNI and AVX-512F,
  BW, and VBMI. The source has a portable fallback, but it was not the measured
  winner.
- **Not a size-independent result.** This historical phase timed only `N=256`;
  later sections report separately contracted N=16/32/64/128 confirmations and
  must not be projected to unmeasured sizes.
- **Not an end-to-end setup or application result.** This historical contract
  is an isolated warm-cache transform microbenchmark. Later setup-inclusive and
  polynomial-product results use different workload hashes and freedom modes.
- **Not a complete search of all possible implementations.** Exhaustiveness is
  only over the frozen four-backend and adjacent-fusion grammar.
- **Not a public native-candidate handoff.** This laboratory episode has no root
  `task.json` or `candidate.json`; its `N=256`, multi-file, setup-excluded path
  is outside native-v1's `N<=128`, single-source, setup-included contract.
- **Not a novelty or human-review result.** No mathematical-novelty,
  strongest-known, or Claude review was available.

The 40.4352x inverse row is also superseded as the current inverse
implementation endpoint by the masked-ternary policy-40 leaf in Section 12.
It remains valid lineage for its frozen 926-byte implementation and contract.

## 11. Reproduction and retention

The reproduction boundary is intentionally stricter than ordinary benchmark
reruns. All successful one-shot timing directories are immutable. Reproduction
means checking their bytes and independently recomputing the recorded
arithmetic; it does not mean drawing a new sample. In particular, do not run
`run_snapshot_lch_confirmation.py`, `run_direct_lch_confirmation.py`, any
`runner/run_once.sh` authorization path, any `runner/measure_once.py`, or an
archived benchmark binary outside an explicit `--check`, `--verify`,
`--self-test`, or `--print-schedule` mode. A command containing
`--execute-one-shot` is always a fresh acquisition command and is forbidden for
reproduction.

### 11.1 Prerequisites

The full repository check expects Linux/x86-64, Bash, Python 3.12, `pip`, CMake,
Ninja, a C++20 compiler, `jq`, `sha256sum`, and the audited
`/usr/bin/bwrap`. Replaying the archived native correctness binaries also
requires GFNI, AVX-512F/BW, and AVX-512VBMI. The retained timing host used a
Threadripper PRO 9985WX, GCC 13.3.0, and Linux `6.8.0-137-generic`; a different
host can reproduce portable correctness and retained statistics but cannot turn
those historical numbers into a new matched performance result.

There is no root `task.json`. This is a multi-file laboratory episode rather
than a public native-v1 attempt, so no `candidate.json` was fabricated.

### 11.2 Fast schema and hash verification

This first command executes no candidate, compiler, or benchmark. It validates
the append-only record schemas, all referenced paths, and all 171 evidence
hashes:

```bash
python3.12 tools/validate_research_handoff.py --repo-root .
```

Expected summary:

```json
{"authoritative_index_entries":41,"evidence_hashes_checked":171,"experiments":85,"hypotheses":69,"passed":true,"schema_version":"addfft-research-handoff-validation-v1","timing_performed":false}
```

The underlying ledger can also be checked using only `jq` and `sha256sum`:

```bash
set -euo pipefail
jq -e . artifacts/lch-postfrontier/authoritative-index.json >/dev/null
while IFS=$'\t' read -r expected path; do
  test -f "$path"
  test "$(sha256sum "$path" | cut -d' ' -f1)" = "$expected"
done < <(
  jq -r '.entries[].evidence[] | [.sha256, .path] | @tsv' \
    artifacts/lch-postfrontier/authoritative-index.json
)
```

### 11.3 Full finite reproduction

From the repository root, run:

```bash
./reproduction.sh
```

The script creates or reuses `.venv`, installs the local package and declared
development dependencies, runs Ruff formatting and lint checks, mypy, all 397
Python tests, and the handoff validator. It then copies `cpp/` into a temporary
workspace, builds fresh Release and ASan/UBSan test binaries, and executes both
inside Bubblewrap with no home directory or network. Finally it replays six
independent postrun auditors from retained observations. The last verified run
created `/tmp/addfft-retained-reproduction.Pm76ps` and performed no timing.

This is finite but not a fully pinned binary environment. Dependency download
may require a configured package cache or network, and fresh compiler metadata
can differ across toolchains. The sealed timing bytes, their SHA-256 values,
and retained-audit outputs must match exactly.

### 11.4 Recompute the promoted statistics without timing

The following commands are narrower alternatives to the root script. They read
retained observations only.

For the N=256 direct forward and policy-40 inverse, each bundle's archived
controller implementation can recompute the complete 50,000-resample
hierarchical bootstrap directly from its 620 retained pairs. Changing into the
archived `controller_sources` directory prevents Python from importing the live
repository package by accident:

```bash
set -euo pipefail
for bundle in \
  artifacts/lch-direct-frontier/\
authenticated-forward-confirmation-v5-20x31-bootstrap50000 \
  artifacts/lch-postfrontier/\
authenticated-inverse-policy40-confirmation-v4-20x31-bootstrap50000; do
  (
    cd "$bundle/controller_sources"
    /usr/bin/python3 -B - ../report.json <<'PY'
import json, sys
from pathlib import Path
from addfft_lab.benchmarking import hierarchical_paired_statistics

report = json.loads(Path(sys.argv[1]).read_text())
baseline, candidate = [], []
for launch in range(report["launches"]):
    rows = [row for row in report["samples"] if row["launch"] == launch]
    baseline.append([row["baseline_ns"] for row in rows])
    candidate.append([row["candidate_ns"] for row in rows])
stats = hierarchical_paired_statistics(
    baseline, candidate, resamples=50_000, seed=0xADDFF7
)
expected = report["statistics"]
assert stats.sample_count == expected["sample_count"] == 620
assert stats.median_speedup == expected["median_speedup"]
assert stats.mad_speedup == expected["mad_speedup"]
assert stats.confidence_low == expected["confidence_low"]
assert stats.confidence_high == expected["confidence_high"]
assert list(stats.per_launch_median_speedup) == expected["per_launch_median_speedup"]
print(stats.median_speedup, stats.confidence_low, stats.confidence_high)
PY
  )
done
```

Expected output is:

```text
34.76302225004469 34.63623642502966 34.84930258922606
46.26860477370302 46.13797685890612 46.3765520020566
```

For N=16 and N=32, the independent `jq` auditor recomputes raw ratios,
order/floor balance, point estimates, and gates for all four workloads:

```bash
set -euo pipefail
ready=artifacts/lch-neighbor-frontier/\
small16-32-masked-fixed-pair-confirmation-ready-v1
result=artifacts/lch-neighbor-frontier/\
small16-32-masked-fixed-pair-confirmation-20x31-cpu0-v1
audit=artifacts/lch-neighbor-frontier/\
small16-32-masked-fixed-pair-confirmation-cpu0-audit-v1
tmpdir="$(mktemp -d --tmpdir addfft-small-audit.XXXXXX)"
for direction in forward inverse; do
  for n in 16 32; do
    contract="$ready/contracts/small_neighbor_confirmation_${direction}_${n}.json"
    contract_sha256="$(sha256sum "$contract" | cut -d' ' -f1)"
    jq -s \
      --slurpfile analysis "$result/analysis/${direction}-${n}.json" \
      --slurpfile contract "$contract" \
      --arg contract_sha256 "$contract_sha256" \
      --argjson expected_cpu 0 \
      -f "$audit/source/audit_small_neighbor_confirmation_result.jq" \
      "$result/results/${direction}-${n}"/launch-*.json \
      > "$tmpdir/${direction}-${n}.json"
    cmp "$tmpdir/${direction}-${n}.json" \
      "$audit/results/${direction}-${n}.json"
  done
done
```

For N=64 and N=128, the archived producer analyzer regenerates each retained
analysis, including its 50,000-resample bootstrap:

```bash
set -euo pipefail
ready=artifacts/lch-neighbor-frontier/\
neighbor64-128-id15-fixed-pair-confirmation-ready-v2
result=artifacts/lch-neighbor-frontier/\
neighbor64-128-id15-fixed-pair-confirmation-20x31-cpu0-v2
tmpdir="$(mktemp -d --tmpdir addfft-neighbor-audit.XXXXXX)"
for direction in forward inverse; do
  for n in 64 128; do
    python3 "$ready/source/analyze_neighbor_confirmation.py" \
      --contract "$ready/contracts/neighbor_confirmation_${direction}_${n}.json" \
      --input-dir "$result/results/${direction}-${n}" \
      --expected-cpu 0 > "$tmpdir/${direction}-${n}.json"
    cmp "$tmpdir/${direction}-${n}.json" \
      "$result/analysis/${direction}-${n}.json"
  done
done
```

There is a provenance limitation in the *separate* N=64/N=128 independent
audit bundle: its retained JSON names auditor source SHA-256 `bc5decda...`, but
those exact source bytes were not archived; the current repository tool is a
different hash. The audit output and index are still hash-bound evidence, but a
third party cannot byte-for-byte regenerate that independent audit from its
original source. The archived producer-analyzer route above is therefore the
reproducible statistics path and is not mislabeled as an independent audit.

The batch-two forward/inverse reports and the narrower batch-four inverse
policy-3-versus-policy-0 report can be regenerated from all retained launches
with their frozen validators:

```bash
set -euo pipefail
tmpdir="$(mktemp -d --tmpdir addfft-fixed-pair-audit.XXXXXX)"
for spec in \
  'forward 0 batch2-forward-p0-confirmation-ready-v1' \
  'inverse 3 batch2-inverse-p3-confirmation-ready-v1' \
  'inverse 3 batch4-inverse-p3-vs-p0-confirmation-ready-v1'; do
  read -r direction policy leaf <<< "$spec"
  bundle="artifacts/lch-postfrontier/$leaf"
  "$bundle/runner/validate_confirmation.py" \
    --direction "$direction" --policy "$policy" \
    --schedule "$bundle/evidence/timing-schedule.json" \
    --result-jsonl "$bundle/runs/one-shot-claimed/success/launches.jsonl" \
    --expected-cpu 0 > "$tmpdir/${leaf}-report.json"
  cmp "$tmpdir/${leaf}-report.json" \
    "$bundle/runs/one-shot-claimed/success/report.json"
done
```

The terminal batch-four, final compact-setup, affine-`ff`, gamma-representation,
and polynomial-product results have separate retained-only auditors:

```bash
set -euo pipefail
for script in \
  artifacts/lch-postfrontier/\
batch4-forward-p0-vs-two-batch2-terminal-confirmation-postrun-audit-v2/reproduction.sh \
  artifacts/lch-postfrontier/\
batch4-inverse-p3-vs-two-batch2-terminal-confirmation-postrun-audit-v2/reproduction.sh \
  artifacts/lch-postfrontier/\
compact-gamma-direct-setup-confirmation-v2-postrun-audit-v1/reproduction.sh \
  artifacts/lch-postfrontier/\
checked-affine-offset-inverse-ff-setup-confirmation-postrun-audit-v1/reproduction.sh; do
  env -i PATH=/usr/bin:/bin LANG=C "$PWD/$script"
done

(
  cd artifacts/lch-postfrontier/\
compact-gamma-standard-compression-screen-postrun-audit-v1
  env -i PATH=/usr/bin:/bin LANG=C ./reproduction.sh
)

product=artifacts/lch-postfrontier/\
degree127-polynomial-product-policy0-confirmation-ready-v1
audit=artifacts/lch-postfrontier/\
degree127-polynomial-product-policy0-confirmation-postrun-audit-v1
tmpdir="$(mktemp -d --tmpdir addfft-product-audit.XXXXXX)"
.venv/bin/python "$audit/source/audit_confirmation_result.py" \
  --bundle "$product" --output "$tmpdir/audit.json"
cmp "$tmpdir/audit.json" "$audit/audit.json"
```

### 11.5 Rebuild representative correctness and negative certificates

These bundle-local scripts compile or audit semantics and machine code, but do
not enter a timing mode:

```bash
set -euo pipefail
for bundle in \
  artifacts/lch-postfrontier/checked-affine-offset-plan-v1 \
  artifacts/lch-postfrontier/checked-affine-offset-plan-direct-gamma-compat-v1 \
  artifacts/lch-postfrontier/exact-alias-screen-ready-v1 \
  artifacts/lch-postfrontier/gf256-automorphism-edge-fusion-static-v1 \
  artifacts/lch-postfrontier/tower-bitslice-level01-static-falsifier-v1 \
  artifacts/lch-neighbor-frontier/n8-xor-circuit-static-falsifier-v1; do
  (cd "$bundle" && env -i PATH=/usr/bin:/bin LANG=C ./reproduction.sh)
done
```

The two algebraic obstruction certificates are deterministic and can be
regenerated byte-for-byte:

```bash
taskset -c 127 .venv/bin/python tools/certify_algorithmic_frontier.py |
  cmp - artifacts/lch-postfrontier/algorithmic-frontier-certificate-v1.json

PYTHONPATH=. taskset -c 127 .venv/bin/python \
  tools/certify_binary_frobenius_frontier.py |
  cmp - artifacts/lch-postfrontier/binary-novel-frobenius-frontier-v1.json
```

The selected source snapshots, raw observations, manifests, reports, and
auditors are still local repository artifacts. A clean or destructive
operation could remove uncommitted evidence. Local controller HMAC receipts
prove the recorded local provenance graph but are not portable public
signatures or external performance authority.

## 12. Postfrontier results

The direct-I/O endpoint triggered a bounded search over inverse lowering,
neighboring sizes, fixed-batch sharing, gamma representation, checked affine
offsets, exact aliasing, compiler/layout effects, restricted circuits, and one
application-level polynomial product. Each search declared its finite grammar
and stopping rule before its final measurement. The canonical status and hash
ledger is
[`authoritative-index.json`](artifacts/lch-postfrontier/authoritative-index.json);
superseded and failed prefixes remain in that ledger rather than disappearing.

### 12.1 Fixed single-transform kernels

The inverse butterfly was improved by a complete 81-policy merge grammar.
Policy 40 replaces 15 live merge operations with masked ternary logic. Its
integrated leaf is 860 bytes and 135 decoded instructions, with 24 GFNI
multiplies, 15 masked ternary operations, and no call, branch, or stack
reference. A fresh 620-pair table/candidate confirmation measured 46.268605x
with a hierarchical 95% interval `[46.137977, 46.376552]`. This supersedes the
926-byte/148-instruction inverse implementation as the current N=256 endpoint;
it does not invalidate the historical 40.4352x artifact for that older body.

Separately contracted neighboring-size confirmations all passed their strict
lower-bound, arm-order, and harness-floor gates:

| Workload | Frozen comparator | Ratio | Hierarchical 95% interval | Wins | Minimum candidate/floor |
|---|---|---:|---:|---:|---:|
| N=16 forward | direct-first-stage full table | 2.510703x | [2.508568, 2.513653] | 620/620 | 1.906385x |
| N=16 inverse | direct-first-stage full table | 2.093623x | [2.024460, 2.097091] | 617/620 | 2.182824x |
| N=32 forward | direct-first-stage full table | 16.233668x | [16.126465, 16.303935] | 620/620 | 1.999760x |
| N=32 inverse | direct-first-stage full table | 5.537158x | [5.531838, 5.546492] | 620/620 | 2.144503x |
| N=64 forward | direct-first-stage full table | 28.251709x | [28.183304, 28.482910] | 620/620 | 1.680509x |
| N=64 inverse | copy-full table | 18.235771x | [18.211458, 18.243526] | 620/620 | 1.537117x |
| N=128 forward | copy-full table | 32.880183x | [32.631987, 33.232124] | 620/620 | 3.041047x |
| N=128 inverse | direct-first-stage full table | 23.257695x | [23.241886, 23.275259] | 620/620 | 2.344007x |

These are not samples from one universal size law. N=16/32 use masked-ZMM
leaves, N=64/128 use frozen lowering 15, and each row has its own workload hash
and strongest admitted repository table formulation. The complete retained
results are indexed under
[`artifacts/lch-neighbor-frontier`](artifacts/lch-neighbor-frontier).
The N=16 inverse candidate lost three individual pairs while still clearing all
hierarchical, order-stratum, and floor gates; the other seven candidates won all
620 pairs.

### 12.2 Fixed-batch and application results

Sharing gamma loads across two transforms produced the following confirmations
over two sequential frozen batch-one calls:

| Direction | Frozen candidate | Comparator/candidate | Hierarchical 95% interval | Wins | Candidate/comparator ns per transform |
|---|---|---:|---:|---:|---:|
| Forward | policy 0 | 1.255705x | [1.253849, 1.258319] | 620/620 | 9.0732 / 11.3909 |
| Inverse | policy 3 | 1.185693x | [1.183598, 1.188872] | 616/620 | 7.4912 / 8.8948 |

The batch-four lowering census itself had mixed results: forward policy 2
failed fresh confirmation at 1.005000x, while inverse policy 3 confirmed a
narrow 1.014997x improvement over policy 0, interval
`[1.012014, 1.017417]`, with 571/620 wins. This internal lowering comparison is
distinct from the terminal batch-four throughput comparison below.

The terminal mechanism comparison used the stronger comparator of two already
confirmed batch-two calls, not merely four batch-one calls. The selected
batch-four forward policy 0 measured 1.181308x with interval
`[1.176661, 1.185979]` and 620/620 wins. Its first validator incorrectly
required equal final checksums even though alternating floor order deliberately
leaves either transformed output or a copied input in the destination. The raw
measurements were not repeated or changed: an independent scalar reconstruction
matched all 20 role-specific checksums and a repaired audit promoted the
preserved prefix. The inverse policy 3 terminal comparison measured 1.354821x
with interval `[1.340274, 1.360994]`, 620/620 wins, and a weakest
layout-by-order median of 1.351772x. Batch eight was not timed because its exact
monolithic bodies used all 32 ZMM registers and spilled heavily.

Under a separate `APPLICATION_CODESIGN` contract, a fused leaf computes the
ordinary product of two 128-coefficient normalized-LCH-novel polynomials into
256 coefficients. The pruned policy-0 versus modular-truncated confirmation
measured 1.136470x, interval `[1.135113, 1.137685]`, and 615/620 wins. A
standalone auditor regenerated every statistic and checksum from an independent
AES-field, monomial-product, and triangular novel-basis implementation. This is
an application implementation result, not evidence that the leaf is a native
LCH transform, globally strongest polynomial multiplier, or new algorithm.

### 12.3 Setup, memory, affine offsets, and aliasing

The arbitrary-basis compact gamma constructor was reduced to exactly 84 field
products (27 triangular, 18 batch-inversion, an 11-product inverse chain, and
28 normalizations) plus 247 subset XORs. For the actual standard basis, a
complete byte-aware screen retained the existing 2,048-byte gamma image: the
28-to-247-byte candidate failed one layout/order gate, while the larger staged
representations passed. The production object retains 2,048 mutable bytes and
charges 2,048 static gamma bytes plus 688 GFNI bytes, 4,784 total versus 75,785
for the flat table plan.

The final setup-inclusive contract allocates and constructs one plan and then
executes one identical N=256 transform in each fresh arm process. Across 620
pairs per direction, flat/compact was 166.470313x forward (95% interval
`[165.750000, 169.839228]`, medians 53,020/320 ns) and 170.483093x inverse
(`[169.998392, 170.677419]`, medians 53,045.5/311 ns), with 620/620 wins in
both directions. This large ratio is specifically the cost of constructing the
65,536-entry flat multiplication table versus copying the selected gamma image
and running one transform; it is not the hot-kernel ratio and excludes process
startup from the internal timer.

For a full independent basis, every affine offset is semantically an XOR-index
permutation. All 256 offsets pass matrix, checked-dispatch, overlap-rejection,
Release, and sanitizer gates. Four setup-excluded representatives initially
screened positively, but the weakest inverse coordinate `ff` failed the final
setup-plus-dispatch contract: general/candidate was 0.647796x, interval
`[0.646152, 0.650117]`, with only 1/620 candidate wins. Thus all-offset
correctness survives, while the nonzero-offset performance branch closes
negative.

Exact `src == dst` operation is also certified: each leaf reads all four source
blocks before the first complete destination store, and the checked API rejects
every partial overlap. The durable benefit is a 256-byte mutable payload instead
of the staged comparator's 512 bytes. Its claim-free timing had no advance gate
and is not promoted as a speed result.

### 12.4 Negative closures

The following are negative only for their declared bounded grammar and exact
contract. The last column states the smallest material change that could reopen
the question; another identical same-host rerun is not such a change.

| Branch | Observation and stopping reason | What could reopen it |
|---|---|---|
| Cached indirect versus direct entry | After repairing a devirtualized false prefix, forward ratios were 0.997217 and 1.034000 in opposite layouts; inverse ratios were 0.999281 and 0.825118. The sign changed with layout and direction, so cached indirect dispatch remained. | A materially different complete caller/binary placement. |
| Gamma alignment, broadcasts, affine substitution | Alignment reached only 1.007380x forward and 0.998993x inverse. Four legal high-stage broadcasts were roughly 0.982--1.003x, and no policy in the complete 32-policy GFNI/affine level grammar advanced. | A cold-gamma contract or different cache/ISA. |
| GF(256) automorphism edge fusion | All eight automorphisms were exact, but each nonidentity leaf was larger. The strongest reduced 24 multiplies to 23 while adding nine affine operations, so the static gate prohibited timing. | A smaller whole-boundary circuit or application-native nonidentity representation. |
| Packed-upper butterflies | Packing reduced 24 GFNI multiplies to 13 but added 12 two-source packs and 24 scatters. The best packed policy was 1.79% slower; the all-packed leaf was 1.474x slower. | A cheaper pack/scatter mechanism or different ISA. |
| Cantor and multiplicative routes | The public LCH graph has 769 useful scalar products. Canonical novel-to-monomial conversion already needs 4,480 nontrivial sites; every one of 128 enumerated Cantor chains needs an 11,120-entry nontrivial public-basis conversion. | A fused exact boundary converter materially smaller than these named graphs; these counts are not a universal lower bound. |
| Compiler switch | A separate-binary Clang build appeared 1.075975x faster. Linking exact GCC/Clang leaves behind one opaque caller produced only 0.99435--0.99655 across layouts and 0.99526--0.99605 across orders; every 1.02 gate failed. | Another compiler, target, CPU, or complete program. |
| N=64 lowering 11 | A 1.012621x single-layout selection became 1.004160 and 0.986946 in the opposite-layout fixed pair. Lowering 15 remained frozen. | A new lowering grammar or another CPU/compiler. |
| N=4/N=8 single transforms | The strongest candidate/table ratios were at most 0.6680 at N=4 and 0.9389 at N=8. N=2 has no nonzero-GFNI mechanism under this witness. | A different public-I/O contract or transform mechanism. |
| N=4/N=8 packed-many | Table-relative speed was large, but candidates were only 0.994--1.093x the exact I/O floor and failed the required 1.25 floor margin. | An application that owns an amortized or different public layout, with transport charged. |
| Synthesized N=8 XOR circuit | Exact mutually inverse circuits used 215/194 XOR gates, but forward compiled to 4,484 bytes, 997 instructions, and 354 stack/spill references; inverse to 4,333/954/342. Timing and conditional N=4 descent stopped. | A radically smaller synthesis or already-bitsliced application contract. |
| Batch-four forward policy 2 | The 1.02718x screen winner confirmed at only 1.005000x, interval `[1.003159, 1.007373]`, with 448/620 wins and both order medians below 1.01. | A new policy grammar, not a rerun of policy 2. |
| Monolithic batch eight | Exact leaves consumed all 32 ZMMs and incurred 111/116 stack references, including 47 spill loads per direction. The preregistered static gate forbade timing and larger all-live batches. | A zero-spill liveness schedule or more architectural vector registers. |
| Affine offset `ff`, setup included | Exact XOR-index transport became 1.5437x slower after allocation, coordinate solving, validation, resolution, and dispatch were charged: general/candidate 0.647796, interval `[0.646152, 0.650117]`, 1/620 wins. | A different setup/dispatch design, CPU, or amortization contract. |
| Standard-gamma compression | The 28-to-247-byte arm had a 1.052952x aggregate ratio but failed layout/order consistency. Subsequent staged gates selected 1,538, 1,792, and finally the existing 2,048-byte image. | A different host or retained-byte/setup-latency policy. |
| Level-0/1 tower bitslicing | All four leaves were exact and spill-free, but level 0 added 28 instructions, 186 bytes, and three memory operands; level 1 added 39--41 instructions, 262--274 bytes, and five memory operands over one memory-source GFNI stage. | A complete charged circuit statically smaller than GFNI, or another ISA. |
| Binary-novel Frobenius compression | The premise is false in the public novel basis: `X_2(2)^2 = 1`, but `X_2(4) = 6`. The violation map has rank 252 and only a four-dimensional compatible kernel. | A separately declared binary-monomial/orbit workload with all basis transport charged. |

Primary evidence for every row, including exact source and result SHA-256
values, is in the `closed-negative` entries of
[`authoritative-index.json`](artifacts/lch-postfrontier/authoritative-index.json).
Representative bundle-local regeneration commands are in Section 11.5.

### 12.5 Diagnostic, superseded, and failed trajectories

Three positive screens deliberately remain diagnostic rather than promoted:

- The first batch-four representative measured 1.679765x forward and 1.791920x
  inverse versus four batch-one calls. The terminal comparisons in Section
  12.2 use the stronger comparator of two confirmed batch-two calls.
- All 256 affine offsets are exact XOR-index permutations. Four setup-excluded
  representatives screened at 1.257--1.323x forward and 1.084--1.161x inverse,
  but the setup-included `ff` result controls the performance decision.
- The final constructor reached 84 products plus 247 XORs and selected the
  2,048-byte standard image. Only the later full setup-plus-one-transform
  campaign promotes the end-to-end effect.

Superseded evidence is not discarded. The separate-binary compiler selection,
short batch-two selection screens, old 2,057-byte compact plan, intermediate
constructor designs, and incomplete initial final-header build remain valid for
their exact old bytes and contracts. They do not control the final endpoint.

Three zero-measurement failures are also retained as failures, not performance
negatives:

- compiler-grid v1 stopped at an unsupported `-march=znver5` and then a
  runner/binary metadata mismatch;
- standard-gamma compression v1 failed working-directory reproduction, while
  v2 later stopped before candidate execution because `runs/` was mode `0555`;
- compact setup v1 failed cross-working-directory byte reproduction because
  GCC embedded the caller directory in `DW_AT_comp_dir`; v2 added deterministic
  path mapping before the sole timing campaign.

These prefixes matter because they identify the earliest false assumption and
the minimal repair. They must not be rewritten as if the successful bundle had
worked on its first attempt.

## 13. Maximum honest claim and remaining boundary

The episode establishes several host-local implementation advantages inside
their exact contracts: fixed N=16/32/64/128/256 LCH kernels, batch-two and
batch-four shared-gamma throughput, final compact-plan setup, and one fused
degree-127 polynomial-product workload. It does **not** combine those ratios,
promote one as a universal FFT speedup, or infer an asymptotic improvement.

Every performance artifact remains `claim_eligible=false`. The comparators are
strongest admitted repository implementations, not independently qualified
Baseline 2; the CPU is one Threadripper PRO 9985WX; GFNI/AVX-512 paths are not
portable evidence; and no human novelty or strongest-known review occurred.
The remaining high-value questions require new authority or a materially new
contract: an independently qualified comparator, a second microarchitecture,
privileged hardware-counter attribution, or human mathematical and novelty
review. Within the declared same-host grammars, no cheap untested mechanism
remains.

## References

1. [`cpp/include/addfft/gf256.hpp`](cpp/include/addfft/gf256.hpp) — finite-field
   definition and scalar multiplication.
2. [`cpp/include/addfft/lch.hpp`](cpp/include/addfft/lch.hpp) — recursive LCH,
   additive points, normalized subspace polynomials, and direct basis
   evaluation.
3. [`cpp/include/addfft/flat_lch.hpp`](cpp/include/addfft/flat_lch.hpp) — flat
   plan, fused butterfly, four backends, gamma construction, and resource
   accounting.
4. [`cpp/src/gfni_lch.cpp`](cpp/src/gfni_lch.cpp) — GFNI/AVX-512 stages,
   register-resident range specializations, and ISA dispatch.
5. [`cpp/src/schedule_search.cpp`](cpp/src/schedule_search.cpp) — exhaustive
   primitive, policy, variant, coordinate, and timing census.
6. [`tools/certify_lch_schedule_search.py`](tools/certify_lch_schedule_search.py)
   — authenticated execution, exact summary/selection validation, and Python
   oracle comparison.
7. [`cpp/src/benchmark.cpp`](cpp/src/benchmark.cpp) — matched warm-cache paired
   timing harness and embedded selected schedule.
8. [`cpp/tests/test_baselines.cpp`](cpp/tests/test_baselines.cpp) — scalar,
   table, nibble, GFNI, coordinate, inverse, and nonzero-offset tests.
9. [`tests/test_schedule_search_certificate.py`](tests/test_schedule_search_certificate.py)
   — certificate tamper and selection-integrity tests.
10. [`tests/test_math_benchmark_contract.py`](tests/test_math_benchmark_contract.py)
    — fairness, portability, receipt, and fail-closed qualification tests.
11. [`artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-manifest.json`](artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-manifest.json)
    — authoritative census manifest, SHA-256
    `0481457fb88344749ec33bad6025262fdc7c9618506f72b7f9a527bfefd6adb1`.
12. [`artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-results.jsonl.gz`](artifacts/lch-schedule-search/authenticated-register-v3c-seed-20260817-results.jsonl.gz)
    — deterministic gzip of the raw census, SHA-256
    `1e11e29c6d3d2eae41abc50a19564e9ad5a1881f20e7973a1cb1d6ae6235646e`;
    decompressed bytes hash to
    `7ea7815c7359f9b0f802cb062a1f00923e9324f56a110cdc7f982737a750726a`.
13. [`artifacts/lch-schedule-search/authenticated-confirmation-20x31-bootstrap50000.json`](artifacts/lch-schedule-search/authenticated-confirmation-20x31-bootstrap50000.json)
    — original fresh confirmation, SHA-256
    `1de944d56ac6ad024bc9e7f8106625c050652fa20cfe8d8c848711104f949ee5`.
14. [`artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000.json`](artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000.json)
    — snapshot-bound repeat, SHA-256
    `2ae4ddd7c4ead90674add9f8449599a1e5e95c6c7db4f76b3de602accf488f73`.
15. [`artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000-context.json`](artifacts/lch-schedule-search/authenticated-snapshot-confirmation-v4-final-20x31-bootstrap50000-context.json)
    — exact command, receipt, code size, seeds, unavailable counters, and source
    snapshot binding; SHA-256
    `9042e0c66b7ba9ba9d2b7f37f4a27c319bf24ef749ee6ff3e0c64a831944fac7`.
16. [`artifacts/lch-schedule-search/README.md`](artifacts/lch-schedule-search/README.md)
    — authoritative evidence index and stale/intermediate artifact warning.
17. [`experiments.jsonl`](experiments.jsonl) and
    [`hypotheses.jsonl`](hypotheses.jsonl) — append-oriented successful, failed,
    and repaired trajectories.
18. [`episode_summary.json`](episode_summary.json) — final episode contract,
    statistics, claim category, unresolved assumptions, and reproduction record.
19. [`reproduction.sh`](reproduction.sh) — finite repository-level reproduction
    entrypoint.
20. [`docs/native-candidate.md`](docs/native-candidate.md) — fail-closed public
    native-v1 size, source-layout, protocol, and setup contract that this
    laboratory experiment does not enter.
21. [`cpp/src/direct_kernel_search.cpp`](cpp/src/direct_kernel_search.cpp) —
    direct-I/O controls, complete 16-member lowering grammar, four residual
    descendants, deterministic correctness gates, and diagnostic timing.
22. [`cpp/src/inverse_kernel_search.cpp`](cpp/src/inverse_kernel_search.cpp) —
    recursive inverse oracle, complete singleton and dense gates, and both
    table-inverse comparators.
23. [`cpp/src/direct_benchmark.cpp`](cpp/src/direct_benchmark.cpp) — frozen
    forward/inverse fixed-pair benchmark, balanced order, I/O preflight, and
    machine-readable contract metadata.
24. [`tools/run_direct_lch_confirmation.py`](tools/run_direct_lch_confirmation.py)
    — receipt verification, absolute-path Bubblewrap isolation, paired
    confirmation, artifact archiving, and fail-closed metadata validation.
25. [`tools/capture_direct_lch_test_evidence.py`](tools/capture_direct_lch_test_evidence.py)
    — source-bound Release and ASan/UBSan correctness evidence.
26. [`artifacts/lch-direct-frontier/authenticated-forward-confirmation-v5-20x31-bootstrap50000`](artifacts/lch-direct-frontier/authenticated-forward-confirmation-v5-20x31-bootstrap50000)
    — authoritative self-contained forward direct confirmation bundle; report
    SHA-256
    `a2e6b001d8bffe729c7b72abe3710877ceb28e3c11873adadab509fc36b14e33`.
27. [`artifacts/lch-direct-frontier/authenticated-inverse-confirmation-v3-20x31-bootstrap50000`](artifacts/lch-direct-frontier/authenticated-inverse-confirmation-v3-20x31-bootstrap50000)
    — authoritative self-contained inverse direct confirmation bundle; report
    SHA-256
    `4e965ee2752656487539a275b5de92c688da8e59b1eaeb1cc561ad4ea610f55c`.
28. [`artifacts/lch-direct-frontier/final-bidirectional-correctness-evidence-v4.json`](artifacts/lch-direct-frontier/final-bidirectional-correctness-evidence-v4.json)
    — exact source manifest, CMake configurations, sandbox commands, test
    binary hashes, exit codes, and outputs; SHA-256
    `4b7383db396a6d5303d30b7b04294e47fd35eee4dcb8cb1873352fadf4c6747b`.
29. [`artifacts/lch-direct-frontier/final-correctness-builds-v4.tar.gz`](artifacts/lch-direct-frontier/final-correctness-builds-v4.tar.gz)
    — complete receipt-bearing Release and ASan/UBSan CMake build trees used by
    reference 28, without the local controller key; SHA-256
    `ab43f9185c5a97258ba17d1334bcfaaece081e402a9c0cd3dca0019cc492b1fb`.
30. [`artifacts/lch-postfrontier/README.md`](artifacts/lch-postfrontier/README.md)
    and
    [`authoritative-index.json`](artifacts/lch-postfrontier/authoritative-index.json)
    — canonical status vocabulary, claim boundaries, paths, and SHA-256 ledger
    for every postfrontier branch.
31. [`artifacts/lch-postfrontier/authenticated-inverse-policy40-confirmation-v4-20x31-bootstrap50000/report.json`](artifacts/lch-postfrontier/authenticated-inverse-policy40-confirmation-v4-20x31-bootstrap50000/report.json)
    — integrated masked-ternary inverse confirmation.
32. [`artifacts/lch-neighbor-frontier/small16-32-masked-fixed-pair-confirmation-20x31-cpu0-v1/index.json`](artifacts/lch-neighbor-frontier/small16-32-masked-fixed-pair-confirmation-20x31-cpu0-v1/index.json)
    and
    [`neighbor64-128 index`](artifacts/lch-neighbor-frontier/neighbor64-128-id15-fixed-pair-confirmation-20x31-cpu0-v2/index.json)
    — N=16/32/64/128 fixed-pair measurements and bootstrap intervals.
33. [`batch-two forward report`](artifacts/lch-postfrontier/batch2-forward-p0-confirmation-ready-v1/runs/one-shot-claimed/success/report.json)
    and
    [`batch-two inverse report`](artifacts/lch-postfrontier/batch2-inverse-p3-confirmation-ready-v1/runs/one-shot-claimed/success/report.json)
    — fixed shared-gamma batch-two confirmations.
34. [`batch-four forward terminal audit`](artifacts/lch-postfrontier/batch4-forward-p0-vs-two-batch2-terminal-confirmation-postrun-audit-v2/evidence/independent-postrun-audit.json)
    and
    [`batch-four inverse terminal report`](artifacts/lch-postfrontier/batch4-inverse-p3-vs-two-batch2-terminal-confirmation-ready-v2/runs/one-shot-claimed/success/report.json)
    — comparisons against two confirmed batch-two calls.
35. [`degree-127 polynomial-product independent audit`](artifacts/lch-postfrontier/degree127-polynomial-product-policy0-confirmation-postrun-audit-v1/audit.json)
    — fixed-pair statistics plus an independent field/product/basis oracle.
36. [`standard-gamma representation audit`](artifacts/lch-postfrontier/compact-gamma-standard-compression-screen-postrun-audit-v1/audit.json)
    and
    [`final constructor index`](artifacts/lch-postfrontier/compact-gamma-constructor-frontier-final-v1/index.json)
    — exact time-memory selection and general-constructor arithmetic census.
37. [`final setup-inclusive result`](artifacts/lch-postfrontier/compact-gamma-direct-setup-confirmation-ready-v2/runs/one-shot-claimed/success/result.json)
    and
    [`clean-room postrun audit`](artifacts/lch-postfrontier/compact-gamma-direct-setup-confirmation-v2-postrun-audit-v1/evidence/audit.json)
    — 1,240 retained flat/compact pairs, independent outputs, schedules, and
    100,000 total bootstrap resamples.
38. [`checked affine-offset audit`](artifacts/lch-postfrontier/checked-affine-offset-inverse-ff-setup-confirmation-postrun-audit-v1/audit.json)
    — retained all-offset lineage and the setup-charged coordinate-`ff`
    closed-negative result.
39. S.-J. Lin, T. Y. Al-Naffouri, and Y. S. Han, “FFT Algorithm for Binary
    Extension Finite Fields and Its Application to Reed-Solomon Codes,”
    [arXiv:1503.05761](https://arxiv.org/abs/1503.05761) — the principal
    additive-FFT/LCH construction used as background here.
40. S.-J. Lin, W.-H. Chung, and Y. S. Han, “Novel Polynomial Basis and Its
    Application to Reed-Solomon Erasure Codes,”
    [arXiv:1404.3458](https://arxiv.org/abs/1404.3458) — background for the
    normalized novel-basis transform and polynomial-product workload.

