#include "avx2_pair_control.h"
bool leo_pair_schedule_enabled = false;
#if defined(LEO_PAIR_TRACE) && LEO_PAIR_TRACE
thread_local LeoPairCounts leo_pair_counts = {};
#endif
bool LeoPairSetMode(unsigned mode)
{
    if (mode > 1) return false;
    leo_pair_schedule_enabled = mode != 0;
    return true;
}
bool LeoPairTraceAvailable()
{
#if defined(LEO_PAIR_TRACE) && LEO_PAIR_TRACE
    return true;
#else
    return false;
#endif
}
void LeoPairTraceReset()
{
#if defined(LEO_PAIR_TRACE) && LEO_PAIR_TRACE
    leo_pair_counts = {};
#endif
}
LeoPairCounts LeoPairTraceGet()
{
#if defined(LEO_PAIR_TRACE) && LEO_PAIR_TRACE
    return leo_pair_counts;
#else
    return {};
#endif
}
