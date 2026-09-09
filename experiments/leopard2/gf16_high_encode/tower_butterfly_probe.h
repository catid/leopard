// Experiment-only, leopard-79h.18.20.3. Not a published backend interface.
#pragma once
#include "tower_avx2_probe.h"

struct TowerLowMapTables { uint8_t row[2][16]; };

// Even byte counts only: full 64-byte ALTMAP blocks plus a compact final
// q-low/q-high tile. Zero bytes touches no pointer. Conversion permits exact
// in-place operation or disjoint buffers, but no partial overlap.
// Actual GF16 high(u*b)=b makes this the same involution in both directions.
extern "C" void tower_convert_involution(const uint8_t* input, uint8_t* output,
                                        size_t bytes, const TowerLowMapTables* map);

// All values/accumulators are already in tower coordinates. The multiplier
// remains a canonical Leopard logarithm; log0 means one, log65535 means the
// butterfly zero-skew sentinel. For that sentinel tables may be null; for
// ordinary logs they must encode exactly the corresponding fixed multiplier.
// Pair buffers must be disjoint. Out/accum require all four buffers pairwise
// disjoint, matching the conservative existing accumulating contract.
extern "C" void tower_ifft_pair(uint8_t* x, uint8_t* y, size_t bytes,
                               uint16_t log, const TowerProductTables* tables);
extern "C" void tower_fft_pair(uint8_t* x, uint8_t* y, size_t bytes,
                              uint16_t log, const TowerProductTables* tables);
extern "C" void tower_fft_out(const uint8_t* x, const uint8_t* y,
                             uint8_t* u, uint8_t* v, size_t bytes,
                             uint16_t log, const TowerProductTables* tables);
extern "C" void tower_ifft_accumulate(const uint8_t* x, const uint8_t* y,
                                     uint8_t* u, uint8_t* v, size_t bytes,
                                     uint16_t log, const TowerProductTables* tables);
