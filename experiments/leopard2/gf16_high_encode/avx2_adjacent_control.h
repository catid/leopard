// Experiment-only quiescent switch; leopard-79h.38.5.4.18.3.1.
#ifndef LEOPARD_EXPERIMENT_ADJACENT_CONTROL_H
#define LEOPARD_EXPERIMENT_ADJACENT_CONTROL_H
#include <cstdint>

// Set only while all codec operations are quiescent. Read-only during concurrent
// execution, as in the preceding inverse-only experiment. Not a production API.
extern bool leo_adjacent_schedule_enabled;
enum LeoAdjacentFamily { LeoAdjacentForward = 0, LeoAdjacentAccumulating = 1 };
struct LeoAdjacentCounts
{
    uint64_t calls[2][2];  // family, OFF/ON
    uint64_t blocks[2][2];
    bool overflow;
};
bool LeoAdjacentSetMode(unsigned mode);
bool LeoAdjacentTraceAvailable();
void LeoAdjacentTraceReset();
LeoAdjacentCounts LeoAdjacentTraceGet();
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
// Keep TLS ownership in its defining translation unit. Clients use a real
// accessor, including overflow tests; no extern-TLS initialization thunk.
LeoAdjacentCounts& LeoAdjacentTraceState();
inline void LeoAdjacentRecord(LeoAdjacentFamily family, bool enabled, uint64_t bytes)
{
    LeoAdjacentCounts& state = LeoAdjacentTraceState();
    uint64_t& calls = state.calls[family][enabled];
    uint64_t& blocks = state.blocks[family][enabled];
    if (calls == UINT64_MAX || blocks > UINT64_MAX - bytes / 64)
        state.overflow = true;
    else { ++calls; blocks += bytes / 64; }
}
#endif
#endif
