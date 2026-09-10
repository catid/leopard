// Link-specific tower diagnostics; shared Release driver for all L2 archives.
// leopard-79h.38.5.4.18.4.2. Count order documented below and in the verifier.
#pragma once
#include <cstdint>
struct LeoTowerCounts {
    // passes, source rows/bytes, output rows/bytes, inverse/forward/accumulating pairs
    uint64_t values[8];
    unsigned initializations; // process lifetime, NOT reset by Reset()
};
const char* LeoTowerPublicIdentity();
bool LeoTowerPublicScheduleValid(const char* schedule);
bool LeoTowerPublicSelect(char state);
unsigned LeoTowerPublicState();
bool LeoTowerPublicDefault();
bool LeoTowerPublicMeasurable();
bool LeoTowerPublicTraced();
void LeoTowerPublicReset();
LeoTowerCounts LeoTowerPublicCounts();
