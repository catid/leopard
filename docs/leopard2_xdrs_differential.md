# Leopard2 and XDRS differential study

## Result

XDRS commit `ae05a779e7f44be13c3d34e14d15b08b4bc02404` was built from
<https://github.com/fastecc/xdrs> on 2026-07-16.  The standalone comparison
explicitly calls `Algorithm2::ReedSolomondecodeL` for the low-rate cases and
`Algorithm3::ReedSolomondecodeH` for the high-rate cases.  This is important:
the checked-in XDRS benchmark selects Algorithm 1 by default.

Both XDRS specialized decoders and the corresponding Leopard2 specialized
decoders recovered four missing originals in every tested case.  XDRS parity
did not match Leopard2 parity after converting the stored symbols between
XDRS's polynomial representation and Leopard's Cantor representation:

| Algorithm | K | R | XDRS rows different from Leopard2 | Recovered originals |
|---|---:|---:|---:|---:|
| XDRS A2 / Leopard2 low | 8 | 248 | 248 | 4 |
| XDRS A2 / Leopard2 low | 16 | 240 | 240 | 4 |
| XDRS A2 / Leopard2 low | 32 | 224 | 224 | 4 |
| XDRS A2 / Leopard2 low | 64 | 192 | 192 | 4 |
| XDRS A2 / Leopard2 low | 128 | 128 | 128 | 4 |
| XDRS A3 / Leopard2 high | 192 | 64 | 64 | 4 |
| XDRS A3 / Leopard2 high | 224 | 32 | 32 | 4 |
| XDRS A3 / Leopard2 high | 240 | 16 | 16 | 4 |
| XDRS A3 / Leopard2 high | 248 | 8 | 8 | 4 |

This is an expected code-definition difference, not a decoder failure.  The
active XDRS initialization uses the binary polynomial-coordinate prefix as its
evaluation order.  Its alternative Cantor-basis initialization is commented
out.  Leopard uses the legacy Cantor-coordinate prefix.  A symbol-basis
conversion maps field values but does not map one systematic evaluation set to
the other, so systematic generator matrices and parity bytes need not agree.
Activating the commented XDRS initialization mechanically was also tested; it
did not pass XDRS's own specialized-decode oracle and is not used as evidence.

The useful differential result is therefore arithmetic and algorithmic:
independently implemented full-length low Algorithm 2 and high Algorithm 3
recover the same source values for their declared coordinate systems.  It is
not evidence that XDRS has the Leopard wire profile.  Legacy wire compatibility
continues to be gated by the old Leopard encoder and frozen golden vectors.

## Reproduction

XDRS remains in the ignored research cache; no Apache-2.0 implementation file
is copied into Leopard.  From the Leopard repository root:

    git clone https://github.com/fastecc/xdrs .research/leopard2/xdrs
    git -C .research/leopard2/xdrs checkout ae05a779e7f44be13c3d34e14d15b08b4bc02404
    cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DLEO2_BUILD_TESTS=ON
    cmake --build build/release -j "$(nproc)"
    c++ -std=c++11 -O2 -mavx2 -include cstring \
      -I. -I.research/leopard2/xdrs/src \
      .research/leopard2/xdrs/src/P_function.cpp \
      experiments/leopard2/xdrs_compare.cpp \
      build/release/libleopard.a -fopenmp -pthread \
      -o .research/leopard2/xdrs_compare
    taskset -c 0-15 .research/leopard2/xdrs_compare

`-include cstring` works around a missing direct include in the research
source.  The comparison allocates 256 skew entries because `init_dec()` touches
index 255, although the XDRS benchmark allocates only 255 entries.  It also
allocates the extra high-rate parity accumulator slots required by the XDRS
encoder.  XDRS uses process-global mutable parameters, so the harness runs
cases sequentially.

## Sources

- XDRS repository and Apache-2.0 license: <https://github.com/fastecc/xdrs>
- XDRS README paper link: <https://i4ai.org/hanyunghsiang/IT2026.pdf>
- Leopard legacy coordinate implementation: <https://github.com/catid/leopard>
