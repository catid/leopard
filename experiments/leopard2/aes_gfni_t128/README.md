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

Selection additionally requires:

- an AUTO context on the calibrated AMD family `0x1a`, model `0x08` host;
- legacy-high GF8, native layout, exactly `K=65, R=65, B=64`;
- the existing packed-terminal and public overlap/scratch validation; and
- a startup-KAT-qualified function pointer.

Unknown CPUs, neighboring sizes, explicit backends, forced builds, decode, and
all non-packed calls use their previous implementation.

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
