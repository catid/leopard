#include "tower_public_link.h"
#include <cstring>
#ifndef LEO_TOWER_IDENTITY
#error Bind archive/profile identity in this shim, never the shared driver.
#endif
const char* LeoTowerPublicIdentity() { return LEO_TOWER_IDENTITY; }
#if defined(LEO_PAIRED_NATIVE) || defined(LEO_TOWER_ORIGINAL)
namespace {
#ifdef LEO_PAIRED_NATIVE
const char kState = 'N';
const char* const kSchedule = "NNNN";
#else
const char kState = 'P';
const char* const kSchedule = "PPPP";
#endif
}
bool LeoTowerPublicScheduleValid(const char* s) { return !std::strcmp(s, kSchedule); }
bool LeoTowerPublicSelect(char s) { return s == kState; }
unsigned LeoTowerPublicState() { return 2; }
bool LeoTowerPublicDefault() { return true; }
bool LeoTowerPublicTraced() { return false; }
void LeoTowerPublicReset() {}
LeoTowerCounts LeoTowerPublicCounts() { return {}; }
#else
#include "tower_encoder.h"
namespace { bool selected = false; }
bool LeoTowerPublicScheduleValid(const char* s)
{
    return !std::strcmp(s,"0110") || !std::strcmp(s,"1001") ||
        !std::strcmp(s,"0000") || !std::strcmp(s,"1111");
}
bool LeoTowerPublicSelect(char s)
{
    if (s != '0' && s != '1') return false;
    tower_encoder::SetEnabled(s == '1');
    selected = s == '1';
    return true;
}
unsigned LeoTowerPublicState() { return selected ? 1 : 0; }
bool LeoTowerPublicDefault()
{
    // Probe the ACTUAL default selector, not just this shim's local state.
    // At OFF this never increments selected-pass counts or initializes tables.
    leopard::backend::Ops ops = {};
    ops.kind = LEO2_BACKEND_AVX2;
    return !selected && !tower_encoder::Select(ops,32768,32768,256,nullptr);
}
bool LeoTowerPublicTraced() { return tower_encoder::TraceAvailable(); }
void LeoTowerPublicReset() { tower_encoder::ResetCounts(); }
LeoTowerCounts LeoTowerPublicCounts()
{
    const tower_encoder::Counts c = tower_encoder::GetCounts();
    return {{c.selected_passes,c.source_rows,c.source_bytes,c.output_rows,c.output_bytes,
        c.inverse_pairs,c.forward_pairs,c.accumulating_pairs},tower_encoder::InitializationCount()};
}
#endif
bool LeoTowerPublicMeasurable()
{
#ifdef __SANITIZE_ADDRESS__
    return false;
#else
    return !LeoTowerPublicTraced();
#endif
}
