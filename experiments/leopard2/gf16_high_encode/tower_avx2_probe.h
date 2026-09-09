// Experiment-only; leopard-79h.18.20.2. Not a codec/backend interface.
#pragma once
#include <cstddef>
#include <cstdint>

// Three GF256 multipliers, two 16-byte nibble tables each: 96 bytes.
struct TowerProductTables { uint8_t row[6][16]; };
// A byte-wide linear map into a pair of bytes: 64 bytes.
struct TowerConvertTables { uint8_t row[4][16]; };

// Complete 64-byte ALTMAP blocks only. Zero blocks is allowed; compact tails
// and partially overlapping buffers are outside this isolated probe contract.
extern "C" void tower_product_blocks(const uint8_t* source, uint8_t* destination,
                                    size_t blocks, const TowerProductTables* tables);
extern "C" void tower_convert_blocks(const uint8_t* source, uint8_t* destination,
                                    size_t blocks, const TowerConvertTables* tables);
