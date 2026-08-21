# T=128 whole-transform GFNI experiment

## Decision

The transferable idea from the 256-point LCH findings works, but the fastest
Leopard realization does not keep the state in AES representation.  For the
exact legacy-high GF8 `K=65, R=65, B=64` packed encode terminal, a row-major
T=128 state and 512-bit direct-affine fixed multiplication are faster than both
the existing AVX2 transform and the same schedule kept in AES coordinates.

The direct-affine leaf is integrated as an optional, operation-specific AUTO
path.  It is not a new public backend and does not change Leopard's wire
representation.  Explicit backends and every forced diagnostic build retain
their exact implementation.

## What “AES representation” means

Leopard and AES use isomorphic copies of `GF(256)`:

- Leopard stores bytes in a Cantor basis for the field with polynomial
  `0x11d`.
- `VGF2P8MULB` interprets bytes in the AES polynomial basis, modulo `0x11b`.

An invertible GF(2)-linear map converts a Leopard byte to the corresponding
AES-field byte.  If the state and every multiplier are converted consistently,
addition remains XOR and multiplication produces the same abstract field
result.  The experiment uses the previously exhaustively checked maps:

    Leopard Cantor -> AES: 0xefd0aaca822e5cae
    AES -> Leopard Cantor: 0xffb08466dc727ea0

Consequently the transform mathematics works entirely in AES coordinates.
The question is only whether two boundary conversions and the chosen layout
cost less than Leopard's native representation.

## Candidates tested

`aes_gfni_t128.cpp` implements three exact T=128 inverse-plus-forward schedules
and compares all of them with Leopard's existing packed transform:

1. Coordinate-axis AES: transpose/gather four payload columns, convert once,
   keep each column in two ZMM registers through all stages, convert back, and
   scatter.
2. Row-major AES: copy the 65 input rows into one 8 KiB state, convert the state
   once, run row-major `VGF2P8MULB` butterflies, convert back, and copy 65 rows.
3. Row-major affine: keep Leopard's Cantor bytes and express each fixed
   multiplier as an 8-by-8 GF(2) matrix for `VGF2P8AFFINEQB`.

The coordinate-axis version most closely follows the supplied register-
resident kernel, but its gather/scatter transpose dominates this packed-shard
layout.  Row-major AES amortizes its conversion and is fast.  Native Cantor
affine avoids the conversions and is slightly faster again.

On the local AMD Threadripper 9980X, CPU 13, GCC 13.3, a 31-sample prototype
screen reported:

| Candidate | Median latency | Paired speedup vs current core |
|---|---:|---:|
| Current AVX2 core | 853.849 ns | 1.000x |
| Coordinate-axis AES/GFNI | 1211.077 ns | 0.717x |
| Row-major AES/GFNI | 314.688 ns | 2.719x |
| Row-major Cantor/GFNI affine | 307.360 ns | 2.782x |

All three candidates matched 65 systematic-basis cases and 129 deterministic
random cases.  These prototype timings select the implementation; the
production claim uses the ordinary benchmark's same-binary whole-call route
control, including public validation and copies.

## Larger-shard follow-up

`larger_shards_screen.cpp` extends the selected native-Cantor affine schedule
without changing the B64 body.  A first candidate transformed one 64-byte
payload vector at a time.  It remained faster through B2048 but fell to about
`0.88x` at B4096, so it was rejected.  The selected candidate transforms up to
four independent vectors together, reusing each skew and affine-matrix load
across a maximum 256-byte row microtile.

The post-hardening, selection-only screen on CPU 13 reported:

| Bytes | One-item batch | One-shot |
|---:|---:|---:|
| 128 | 2.556x | 2.544x |
| 256 | 2.231x | 2.233x |
| 512 | 2.185x | 2.150x |
| 1024 | 2.077x | 1.971x |
| 2048 | 1.944x | 2.004x |
| 4096 | 1.637x | 1.531x |

These numbers select the 256-byte microtile; they are not promotion evidence.
The retained v33 runner performs 25-round ABBA confirmation at every target,
16-times-work checks at B128 and B4096, and exact selector/call/tile/digest
checks at inactive byte and K/R neighbors.  Neighbor timing is descriptive and
cannot reject an otherwise exact no-op route.

### Frozen larger-shard confirmation

The authoritative v4 campaign at source commit
`1e5de49cde826f461ccefa035e5781dd6c0d8918` completed on the same
Threadripper 9980X class host, pinned to CPU 52 with SMT sibling 116 reserved.
It compared the enabled and disabled leaves in one frozen binary.  Each target
used 25 ABBA/BAAB rounds, 31 retained samples per launch, and both ordinary
one-item batch encode and one-shot encode as co-primary metrics.

| Bytes | One-item batch speedup (95% CI) | One-shot speedup (95% CI) |
|---:|---:|---:|
| 128 | 2.542x [2.523, 2.560] | 2.540x [2.525, 2.555] |
| 256 | 2.298x [2.270, 2.327] | 2.259x [2.230, 2.289] |
| 512 | 2.131x [2.118, 2.145] | 2.129x [2.110, 2.147] |
| 1024 | 1.958x [1.937, 1.979] | 1.957x [1.936, 1.979] |
| 2048 | 1.917x [1.893, 1.941] | 1.915x [1.889, 1.941] |
| 4096 | 1.587x [1.572, 1.601] | 1.569x [1.553, 1.586] |

The 16-times-work endpoints also passed: B128 measured 2.519x
[2.498, 2.540] and 2.517x [2.494, 2.542], while B4096 measured 1.580x
[1.551, 1.609] and 1.574x [1.544, 1.604] for batch and one-shot respectively.
Both all-or-none promotion regions therefore passed their predeclared 1.05
lower-confidence-bound gate.  All 16 inactive byte and K/R neighbors reported
selector false, zero candidate calls, zero vector tiles, and matching workload
digests.  The controller accepted 1,248 launches and retained but excluded 56
launches from 14 SMT-contaminated attempts.

The sealed local evidence is
`.research/leopard-8bd/1e5de49c-final-v4`.  Its principal hashes are:

- campaign report: `3f858cc1777100a48a073e2a916c4a621ce65274a6fee5489b04895b7845b11d`;
- event journal: `549c153e917597d1856bbce5e46bfaf2c07722ee689a499b95c28f91a4c1054d`;
- frozen benchmark: `83208fb63c0461faa385e36ccaa4a82ff2af69c467594bd8037ca06a434504ba`;
- canonical source archive: `d9c3506d1a2fa7eaab395b6ed36dd33cff1a8e3fb1ca601d325325ab060f6f00`;
- complete checksum index: `10603847fe8ec49cd3e56b6cb7f74e423a545cd2314638b4348da59239975cd2`.

Earlier v1 and v2 attempts failed before collecting timings.  V3 retained 20
launches, rejected every one for SMT-sibling activity, and accepted no round;
none contributes to this claim.  The confirmation is a host-local warm-cache
result against the current mature same-binary operation leaf, not an
exact-main, other-CPU, decode, multi-item, or neighboring-size claim.

## Production scope and gates

`Leopard2BackendAVX512GFNIT128.cpp` contains the selected row-major affine leaf
in its own translation unit.  It is compiled only with `-mavx2 -mavx512f
-mavx512bw -mavx512vl -mgfni`.  Baseline code calls its initializer only after
CPUID and XCR0 establish AVX2, AVX-512F/BW/VL, GFNI, and OS-managed ZMM state,
and only on the calibrated processor model.
Initialization regenerates all multiplier matrices from Leopard's scalar field
tables, exhaustively checks all 65,536 fixed-multiplier byte results, and runs
two complete scalar-versus-vector transform KATs.  Failure leaves the pointer
null and preserves the established path.

The original B64 leaf additionally requires:

- an AUTO context on the calibrated AMD family `0x1a`, model `0x08` host;
- legacy-high GF8, native layout, exactly `K=65, R=65, B=64`;
- the existing packed-terminal and public overlap/scratch validation; and
- a startup-KAT-qualified function pointer.

Unknown CPUs, neighboring sizes, explicit backends, forced builds, decode, and
all non-packed calls use their previous implementation.

The larger-shard callback is independently KAT-qualified at B128 and the
non-promoted B320 tail shape.  It uses four aligned 8 KiB states and is selected
only for ordinary one-shot or one-item batch calls with dense packed full
output at B128, B256, B512, B1024, B2048, or B4096.  Its diagnostic mode and
route counters are independent from the older B64 schema, so either campaign
cannot accidentally activate the other leaf.

## Reproduction

Build the standalone comparison against a Release library:

```sh
g++ -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Werror -I. \
  -mavx2 -mavx512f -mavx512bw -mavx512vbmi -mgfni \
  experiments/leopard2/aes_gfni_t128/aes_gfni_t128.cpp \
  <build>/libleopard.a -fopenmp -pthread -o <output>
taskset -c <isolated-cpu> <output>
```

The production benchmark route is schema v32 in `bench_leopard2`, selected by
`--k65r65-b64-avx512-gfni-mode 0|1`.  Its untimed route probe verifies that both
one-item batch and one-shot encode actually enter the candidate when expected;
the probe is normalized away before timing.

The larger route is schema v33 and uses
`--k65r65-t128-avx512-gfni-mode 0|1`.  The frozen confirmation controller is
`run_larger_shards_abba.py`; it refuses writable or hard-linked benchmark
artifacts, holds the canonical benchmark lock, verifies singleton CPU/SMT
isolation, exact-compares a canonical Git archive with the declared clean
commit and tree, executes the frozen JSON validator directly from its hashed
source bytes, journals every raw launch before applying gates, and writes an
exclusive atomic checkpoint report after every cell.
