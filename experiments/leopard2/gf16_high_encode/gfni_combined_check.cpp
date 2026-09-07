// Check-only public workload plus the already validated callback observer.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#include "gf16_callback_probe.cpp"
#include "gfni_combined.h"

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
            "usage: --check cell[0..5] --mode=0|1|2|3 [parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '5', "invalid cell");
        const unsigned mode = gfni_combined::ParseMode(argv[3]);
        const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
        omp_set_dynamic(0); omp_set_num_threads(1);
        gfni_combined::Reset(mode);
        callback_probe::Reset();
        char* workload[] = {argv[0],argv[1],argv[2],argc == 5 ? argv[4] : NULL,NULL};
        Require(CallbackWorkloadMain(argc == 5 ? 4 : 3, workload) == 0, "public check failed");
        const auto& state = gfni_combined::Get();
        gfni_combined::Print();
        const unsigned passes[] = {2,2,1,1,2,1};
        Require(state.calls == passes[cell] && state.matches == (cell == 0 ? 2U : 0U) &&
            state.first == ((mode & 1U) ? state.matches : 0U) &&
            state.terminal == ((mode & 2U) ? state.matches : 0U), "pass selection counts");
        Require(callback_probe::state.pass_count == state.calls, "observer pass count");
        callback_probe::Print();
        Require(kernel_calls == (cell == 0 && (mode & 2U) ? 6U : 0U), "terminal kernel call count");
        std::fprintf(stderr, "{\"schema\":\"gfni-terminal-kernel-counts/v1\",\"calls\":%u,\"lane_groups\":%u,\"timed\":false}\n",kernel_calls,kernel_groups);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "combined check: %s\n", error.what());
        return 1;
    }
}
