// Experiment-only same-binary control, leopard-79h.38.5.4.18.2.
#ifndef LEOPARD_EXPERIMENT_AVX2_PAIR_CONTROL_H
#define LEOPARD_EXPERIMENT_AVX2_PAIR_CONTROL_H
#include <cstdint>

// Set only with all codec operations quiescent, before the measured sequence.
// The ordinary bool is immutable during concurrent codec execution. No state
// is installed in the production library by this experimental interface.
extern bool leo_pair_schedule_enabled;
struct LeoPairCounts
{
    uint64_t off_calls, on_calls, off_blocks, on_blocks;
    bool overflow;
};
bool LeoPairSetMode(unsigned mode);
bool LeoPairTraceAvailable();
void LeoPairTraceReset();
LeoPairCounts LeoPairTraceGet();

#if defined(LEO_PAIR_TRACE) && LEO_PAIR_TRACE
extern thread_local LeoPairCounts leo_pair_counts;
inline void LeoPairRecord(bool enabled, uint64_t bytes)
{
    uint64_t& calls = enabled ? leo_pair_counts.on_calls : leo_pair_counts.off_calls;
    uint64_t& blocks = enabled ? leo_pair_counts.on_blocks : leo_pair_counts.off_blocks;
    if (calls == UINT64_MAX || blocks > UINT64_MAX - bytes / 64)
        leo_pair_counts.overflow = true;
    else { ++calls; blocks += bytes / 64; }
}
#endif
#endif
