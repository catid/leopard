// Check-only public workload plus the already validated callback observer.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#include "gf16_callback_probe.cpp"
#include "gfni_inverse_distance1.h"

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 4 || argc == 5) && std::strcmp(argv[1], "--check") == 0,
            "usage: --check cell[0..5] --fuse=0|--fuse=1 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '5', "invalid cell");
        Require(std::strcmp(argv[3], "--fuse=0") == 0 || std::strcmp(argv[3], "--fuse=1") == 0,
            "invalid experiment mode");
        const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
        const bool enabled = argv[3][7] == '1';
        omp_set_dynamic(0); omp_set_num_threads(1);
        gfni_inverse_distance1::Reset(enabled);
        callback_probe::Reset();
        char* workload[] = {argv[0],argv[1],argv[2],argc == 5 ? argv[4] : NULL,NULL};
        Require(CallbackWorkloadMain(argc == 5 ? 4 : 3, workload) == 0, "public check failed");
        const auto& state = gfni_inverse_distance1::Get();
        gfni_inverse_distance1::Print();
        const unsigned passes[] = {2,2,1,1,2,1};
        Require(state.calls == passes[cell] && state.matches == (cell == 0 ? 2U : 0U) &&
            state.changed == (enabled ? state.matches : 0), "pass selection counts");
        Require(callback_probe::state.pass_count == state.calls, "observer pass count");
        callback_probe::Print();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "distance-one check: %s\n", error.what());
        return 1;
    }
}
