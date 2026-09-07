// Check-only public workload plus the already validated callback observer.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#include "gf16_callback_probe.cpp"
#include "gfni_terminal.h"

static unsigned kernel_calls=0, kernel_groups=0;
extern "C" void __real_LeoGFNIFinalAccumulate(const void* const*,void* const*,unsigned,uint16_t,uint16_t,uint16_t,uint64_t);
extern "C" void __wrap_LeoGFNIFinalAccumulate(const void* const* input,void* const* sums,unsigned distance,
    uint16_t a,uint16_t b,uint16_t c,uint64_t bytes)
{
    Require(callback_probe::active && callback_probe::active->kind == LEO2_BACKEND_GFNI &&
        kernel_calls<16 && distance==64 && bytes==32768 && a!=65535 && b!=65535 && c!=65535,
        "unexpected terminal kernel invocation");
    ++kernel_calls; kernel_groups+=distance;
    __real_LeoGFNIFinalAccumulate(input,sums,distance,a,b,c,bytes);
}

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 4 || argc == 5) && std::strcmp(argv[1], "--check") == 0,
            "usage: --check cell[0..5] --terminal=0|--terminal=1 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '5', "invalid cell");
        Require(std::strcmp(argv[3], "--terminal=0") == 0 || std::strcmp(argv[3], "--terminal=1") == 0,
            "invalid experiment mode");
        const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
        const bool enabled = std::strcmp(argv[3], "--terminal=1") == 0;
        omp_set_dynamic(0); omp_set_num_threads(1);
        gfni_terminal::Reset(enabled);
        callback_probe::Reset();
        char* workload[] = {argv[0],argv[1],argv[2],argc == 5 ? argv[4] : NULL,NULL};
        Require(CallbackWorkloadMain(argc == 5 ? 4 : 3, workload) == 0, "public check failed");
        const auto& state = gfni_terminal::Get();
        gfni_terminal::Print();
        const unsigned passes[] = {2,2,1,1,2,1};
        Require(state.calls == passes[cell] && state.matches == (cell == 0 ? 2U : 0U) &&
            state.changed == (enabled ? state.matches : 0), "pass selection counts");
        Require(callback_probe::state.pass_count == state.calls, "observer pass count");
        callback_probe::Print();
        Require(kernel_calls == (cell == 0 && enabled ? 6U : 0U), "terminal kernel call count");
        std::fprintf(stderr, "{\"schema\":\"gfni-terminal-kernel-counts/v1\",\"calls\":%u,\"lane_groups\":%u,\"timed\":false}\n",kernel_calls,kernel_groups);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "terminal check: %s\n", error.what());
        return 1;
    }
}
