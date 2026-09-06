// Same executable for overlay-OFF and fused-ON; no Ops/libc/linker wrappers.
// The only hook is the once-per-pass hook in the frozen experimental FF16 object.
#define main Distance1WorkloadMain
#include "current_route_screen.cpp"
#undef main
#include "gfni_inverse_distance1.h"
#include <omp.h>

static_assert(gfni_inverse_distance1::kCapacity == 64, "timing requires 64-pass trace");

int main(int argc, char** argv)
{
    try
    {
        using namespace gfni_inverse_distance1;
        Require((argc == 4 || argc == 5) &&
                (std::strcmp(argv[1], "--check") == 0 ||
                 std::strcmp(argv[1], "--measure") == 0 ||
                 std::strcmp(argv[1], "--exercise") == 0),
                "usage: --check|--measure|--exercise cell[0..5] --fuse=0|--fuse=1 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '5', "invalid cell");
        Require(std::strcmp(argv[3], "--fuse=0") == 0 || std::strcmp(argv[3], "--fuse=1") == 0,
                "invalid fusion mode");
        const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
        const bool enabled = argv[3][7] == '1';
        const bool measured = std::strcmp(argv[1], "--measure") == 0;
        const bool exercise = std::strcmp(argv[1], "--exercise") == 0;
        Require(argc == 4 || (!measured && !exercise), "parity dump is check-only");
        omp_set_dynamic(0); omp_set_num_threads(1);
        Reset(enabled);
        char check_mode[] = "--check";
        char* workload[] = {argv[0],exercise ? check_mode : argv[1],argv[2],argc == 5 ? argv[4] : NULL,NULL};
        // Exercise checks 52-pass capacity without reading a clock. It does not
        // claim to reproduce the timing loop's allocation/warmup behavior.
        for (unsigned i=0; i<(exercise ? 26U : 1U); ++i)
            Require(Distance1WorkloadMain(argc == 5 ? 4 : 3, workload) == 0, "public workload failed");
        const unsigned encodes = measured || exercise ? 26 : 1;
        const unsigned passes[] = {2,2,1,1,2,1};
        const unsigned kinds[] = {6,3,5,3,3,3};
        const uint64_t tiles[] = {32768,32768,65536,32768,32768,4096};
        const State& state = Get();
        Require(state.calls == encodes*passes[cell] && state.matches == (cell == 0 ? encodes*2 : 0) &&
            state.changed == (enabled ? state.matches : 0), "pass totals");
        for (unsigned i=0; i<state.calls; ++i)
        {
            const Call& c = state.records[i];
            Require(c.kind == kinds[cell] && c.k == kCells[cell].k && c.r == kCells[cell].r &&
                c.requested == kCells[cell].r && c.side == (cell == 5 ? 512U : 256U) &&
                c.sparse_blocks == 0 && c.sparse_present && c.bytes == tiles[cell] &&
                c.source_policy == kCells[cell].bytes, "pass identity");
        }
        std::fprintf(stderr,"{\"schema\":\"gfni-inverse-distance1-timing/v1\",\"cell\":%u,\"enabled\":%s,"
            "\"encodes\":%u,\"calls\":%u,\"matches\":%u,\"changed\":%u,\"timed\":%s,\"exercise\":%s}\n",
            cell,enabled ? "true" : "false",encodes,state.calls,state.matches,state.changed,
            measured ? "true" : "false",exercise ? "true" : "false");
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,"distance-one timing driver: %s\n",error.what());
        return 1;
    }
}
