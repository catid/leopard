#include "avx2_adjacent_control.h"
#include <cstdio>
#include <initializer_list>
#include <stdexcept>
#include <thread>

static void Check(bool value) { if (!value) throw std::runtime_error("control check"); }
static bool Empty(const LeoAdjacentCounts& counts)
{
    for (unsigned f=0; f<2; ++f) for (unsigned s=0; s<2; ++s)
        if (counts.calls[f][s] || counts.blocks[f][s]) return false;
    return !counts.overflow;
}
int main()
{
    try {
        Check(!leo_adjacent_schedule_enabled);
        for (unsigned mode : {0U,1U,0U,1U}) {
            Check(LeoAdjacentSetMode(mode) && leo_adjacent_schedule_enabled==(mode!=0));
            Check(!LeoAdjacentSetMode(2) && !LeoAdjacentSetMode(UINT32_MAX));
            Check(leo_adjacent_schedule_enabled==(mode!=0));
        }
        LeoAdjacentTraceReset(); Check(Empty(LeoAdjacentTraceGet()));
#if defined(LEO_ADJACENT_TRACE) && LEO_ADJACENT_TRACE
        Check(LeoAdjacentTraceAvailable());
        for (unsigned f=0; f<2; ++f) for (unsigned s=0; s<2; ++s) {
            LeoAdjacentTraceReset();
            for (uint64_t bytes : {0U,2U,62U,64U,66U,128U})
                LeoAdjacentRecord(static_cast<LeoAdjacentFamily>(f),s!=0,bytes);
            auto counts = LeoAdjacentTraceGet();
            Check(counts.calls[f][s]==6 && counts.blocks[f][s]==4 && !counts.overflow);
            counts.calls[f][s]=0; counts.blocks[f][s]=0; Check(Empty(counts));
            std::thread worker([] { Check(Empty(LeoAdjacentTraceGet()));
                LeoAdjacentRecord(LeoAdjacentForward,false,128);
                Check(LeoAdjacentTraceGet().blocks[0][0]==2); });
            worker.join(); Check(LeoAdjacentTraceGet().calls[f][s]==6);
            LeoAdjacentTraceState().calls[f][s]=UINT64_MAX;
            LeoAdjacentRecord(static_cast<LeoAdjacentFamily>(f),s!=0,64);
            Check(LeoAdjacentTraceGet().overflow && LeoAdjacentTraceGet().blocks[f][s]==4);
            LeoAdjacentTraceReset(); LeoAdjacentTraceState().blocks[f][s]=UINT64_MAX;
            LeoAdjacentRecord(static_cast<LeoAdjacentFamily>(f),s!=0,64);
            Check(LeoAdjacentTraceGet().overflow && LeoAdjacentTraceGet().calls[f][s]==0);
        }
#else
        Check(!LeoAdjacentTraceAvailable());
#endif
        LeoAdjacentTraceReset(); Check(Empty(LeoAdjacentTraceGet()));
        Check(LeoAdjacentSetMode(0) && !leo_adjacent_schedule_enabled);
        std::printf("{\"schema\":\"adjacent-control-unit/v1\",\"trace\":%s,\"timed\":false}\n",
                    LeoAdjacentTraceAvailable() ? "true" : "false");
        return 0;
    } catch (const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
