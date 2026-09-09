// Current-production path counts only; leopard-79h.38.5.4.18.3.
// The original archive is unchanged. Pair-control functions are driver-only
// compatibility stubs; no experimental kernel is linked by this observer.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#define LEO_GF16_CALLBACK_WORKLOAD "avx2_pair_screen.cpp"
#include "gf16_callback_probe.cpp"

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 4 || argc == 5) && !std::strcmp(argv[1], "--check") &&
                std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '7' &&
                !std::strcmp(argv[3], "off"),
                "check-only production observer: --check cell[0..7] off [parity_file]");
        omp_set_dynamic(0); omp_set_num_threads(1);
        callback_probe::Reset();
        const int result = CallbackWorkloadMain(argc, argv);
        if (result) return result;
        const bool gf8 = argv[2][0] == '7';
        Require(gf8 ? callback_probe::state.pass_count == 0 && callback_probe::state.calls == 0
                    : callback_probe::state.pass_count > 0 && callback_probe::state.calls > 0,
                "unexpected field observation");
        callback_probe::Print();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "adjacent observer: %s\n", error.what());
        return 1;
    }
}
