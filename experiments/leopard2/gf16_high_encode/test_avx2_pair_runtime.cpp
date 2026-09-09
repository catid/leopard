// Reuse the deterministic kernel/API matrix with both same-binary states.
#define main LeoPairFocusedMain
#include "test_avx2_pair_schedule.cpp"
#undef main
#include "avx2_pair_control.h"

int main(int argc, char** argv)
{
    try
    {
        Require(argc == 3 && (!std::strcmp(argv[1], "off") || !std::strcmp(argv[1], "on")),
            "usage: off|on focused_selector");
        const bool enabled = !std::strcmp(argv[1], "on");
        Require(LeoPairSetMode(enabled ? 1 : 0), "mode setter");
        Require(!LeoPairSetMode(2) && leo_pair_schedule_enabled == enabled, "setter validation");
        LeoPairTraceReset();
        const int result = LeoPairFocusedMain(argc - 1, argv + 1);
        Require(leo_pair_schedule_enabled == enabled && !LeoPairTraceGet().overflow, "mode/trace drift");
        return result;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
