// Both states of the same experimental archive; no benchmark-clock interface.
// Generated from the unchanged qualified harness with only its final entry
// point renamed; nested legacy main-renaming macros remain untouched.
#include "adjacent_focused.inc"
#include "avx2_adjacent_control.h"

int main(int argc, char** argv)
{
    try {
        Require(argc==3 && (!std::strcmp(argv[1],"off") || !std::strcmp(argv[1],"on")),
                "usage: off|on focused_selector");
        Require(!leo_adjacent_schedule_enabled,"fresh process must default OFF");
        const bool enabled = !std::strcmp(argv[1],"on");
        Require(LeoAdjacentSetMode(enabled ? 1 : 0),"mode setter");
        Require(!LeoAdjacentSetMode(2) && leo_adjacent_schedule_enabled==enabled,"invalid mode");
        LeoAdjacentTraceReset();
        const int code = LeoAdjacentFocusedMain(argc-1,argv+1);
        Require(leo_adjacent_schedule_enabled==enabled && !LeoAdjacentTraceGet().overflow,"state/trace drift");
        const auto counts = LeoAdjacentTraceGet();
        std::printf("{\"schema\":\"adjacent-runtime-focused/v1\",\"mode\":%u,\"trace\":%s,"
                    "\"calls\":[[%llu,%llu],[%llu,%llu]],\"blocks\":[[%llu,%llu],[%llu,%llu]],\"timed\":false}\n",
                    unsigned(enabled),LeoAdjacentTraceAvailable() ? "true" : "false",
                    (unsigned long long)counts.calls[0][0],(unsigned long long)counts.calls[0][1],
                    (unsigned long long)counts.calls[1][0],(unsigned long long)counts.calls[1][1],
                    (unsigned long long)counts.blocks[0][0],(unsigned long long)counts.blocks[0][1],
                    (unsigned long long)counts.blocks[1][0],(unsigned long long)counts.blocks[1][1]);
        return code;
    } catch (const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
