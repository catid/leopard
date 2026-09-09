#include "avx2_adjacent_public_link.h"
#include <cstring>
#ifndef LEO_ADJACENT_IDENTITY
#error Bind archive and profile identity in this shim, never in the driver.
#endif
const char* LeoAdjacentPublicIdentity() { return LEO_ADJACENT_IDENTITY; }
#if defined(LEO_PAIRED_NATIVE) || defined(LEO_ADJACENT_ORIGINAL)
namespace {
#ifdef LEO_PAIRED_NATIVE
const char kState = 'N';
const char* const kSchedule = "NNNN";
#else
const char kState = 'P';
const char* const kSchedule = "PPPP";
#endif
}
bool LeoAdjacentPublicScheduleValid(const char* s) { return !std::strcmp(s, kSchedule); }
bool LeoAdjacentPublicSelect(char s) { return s == kState; }
unsigned LeoAdjacentPublicState() { return 2; }
bool LeoAdjacentPublicDefault() { return true; }
bool LeoAdjacentPublicTraced() { return false; }
void LeoAdjacentPublicReset() {}
LeoAdjacentCounts LeoAdjacentPublicCounts() { return {}; }
#else
bool LeoAdjacentPublicScheduleValid(const char* s)
{
    return !std::strcmp(s, "0110") || !std::strcmp(s, "1001") ||
        !std::strcmp(s, "0000") || !std::strcmp(s, "1111");
}
bool LeoAdjacentPublicSelect(char s)
{
    return (s == '0' || s == '1') && LeoAdjacentSetMode(static_cast<unsigned>(s - '0')) &&
        leo_adjacent_schedule_enabled == (s == '1');
}
unsigned LeoAdjacentPublicState() { return leo_adjacent_schedule_enabled ? 1 : 0; }
bool LeoAdjacentPublicDefault() { return !leo_adjacent_schedule_enabled; }
bool LeoAdjacentPublicTraced() { return LeoAdjacentTraceAvailable(); }
void LeoAdjacentPublicReset() { LeoAdjacentTraceReset(); }
LeoAdjacentCounts LeoAdjacentPublicCounts() { return LeoAdjacentTraceGet(); }
#endif
bool LeoAdjacentPublicMeasurable()
{
#ifdef __SANITIZE_ADDRESS__
    return false;
#else
    return !LeoAdjacentPublicTraced();
#endif
}
