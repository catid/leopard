#include "avx2_adjacent_control.h"
bool leo_adjacent_schedule_enabled = false;
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
static thread_local LeoAdjacentCounts leo_adjacent_counts = {};
LeoAdjacentCounts& LeoAdjacentTraceState() { return leo_adjacent_counts; }
#endif
bool LeoAdjacentSetMode(unsigned mode)
{
    if (mode > 1) return false;
    leo_adjacent_schedule_enabled = mode != 0;
    return true;
}
bool LeoAdjacentTraceAvailable()
{
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
    return true;
#else
    return false;
#endif
}
void LeoAdjacentTraceReset()
{
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
    leo_adjacent_counts = {};
#endif
}
LeoAdjacentCounts LeoAdjacentTraceGet()
{
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
    return leo_adjacent_counts;
#else
    return {};
#endif
}
