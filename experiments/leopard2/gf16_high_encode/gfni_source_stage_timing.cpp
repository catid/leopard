// Driver-only .38.5.4.11: one executable for both modes, no libc wrappers.
// Reuse the exact six public workloads and the existing 1+4+21 encode loop.
#define main SourceStageWorkloadMain
#include "current_route_screen.cpp"
#undef main
#include "gfni_source_stage_probe.h"
#include <omp.h>

static_assert(gfni_source_stage_probe::kCapacity == 64,
              "timing driver requires the checked 64-record trace");

int main(int argc, char** argv)
{
    try
    {
        using namespace gfni_source_stage_probe;
        Require((argc == 4 || argc == 5) &&
                (std::strcmp(argv[1], "--check") == 0 ||
                 std::strcmp(argv[1], "--measure") == 0 ||
                 std::strcmp(argv[1], "--exercise") == 0),
                "usage: gfni_source_stage_timing --check|--measure|--exercise cell --stage=0|--stage=1 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' &&
                argv[2][0] <= '5', "invalid cell");
        Require(std::strcmp(argv[3], "--stage=0") == 0 ||
                std::strcmp(argv[3], "--stage=1") == 0, "invalid stage mode");
        const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
        const bool enabled = argv[3][8] == '1';
        const bool measured = std::strcmp(argv[1], "--measure") == 0;
        const bool exercise = std::strcmp(argv[1], "--exercise") == 0;
        Require(argc == 4 || (!measured && !exercise), "parity dump is check-only");
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        Reset(enabled);
        char check_mode[] = "--check";
        char* workload_argv[] = {argv[0], exercise ? check_mode : argv[1],
            argv[2], argc == 5 ? argv[4] : NULL, NULL};
        // --exercise executes 26 independent check-mode workloads, never a
        // clock. It tests trace capacity/counts, not the measurement loop's
        // allocation or warmup behavior. --measure calls that loop unchanged.
        for (unsigned i = 0; i < (exercise ? 26U : 1U); ++i)
        {
            const int result = SourceStageWorkloadMain(argc == 5 ? 4 : 3, workload_argv);
            Require(result == 0, "public workload failed");
        }
        const unsigned encodes = measured || exercise ? 26 : 1;
        const unsigned passes[] = {2, 2, 1, 1, 2, 1};
        const unsigned kinds[] = {6, 3, 5, 3, 3, 3};
        const size_t tiles[] = {32768, 32768, 65536, 32768, 32768, 4096};
        const State& state = Get();
        Require(state.calls == encodes * passes[cell] &&
                state.matches == (cell == 0 ? encodes * 2 : 0) &&
                state.changed == (enabled ? state.matches : 0), "trace totals");
        for (unsigned i = 0; i < state.calls; ++i)
        {
            const Call& call = state.records[i];
            Require(call.kind == kinds[cell] && call.k == kCells[cell].k &&
                    call.r == kCells[cell].r && call.requested == kCells[cell].r &&
                    call.side == (cell == 5 ? 512U : 256U) && call.sparse_blocks == 0 &&
                    call.bytes == tiles[cell] && call.source_policy == kCells[cell].bytes &&
                    call.effective_policy == (cell == 0 && enabled ? 16384 : kCells[cell].bytes),
                    "changed internal pass identity");
        }
        std::fprintf(stderr,
            "{\"schema\":\"gfni-source-stage-timing/v1\",\"cell\":%u,\"enabled\":%s,"
            "\"encodes\":%u,\"calls\":%u,\"matches\":%u,\"changed\":%u,"
            "\"timed\":%s,\"exercise\":%s}\n", cell, enabled ? "true" : "false",
            encodes, state.calls, state.matches, state.changed,
            measured ? "true" : "false", exercise ? "true" : "false");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "source-stage timing driver: %s\n", error.what());
        return 1;
    }
}
