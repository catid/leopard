// Driver-only current-production callback observation; leopard-79h.38.5.4.19.
#define LEO_GF16_CALLBACK_NO_MAIN 1
#define LEO_GF16_CALLBACK_WORKLOAD "r199_boundary_screen.cpp"
#include "gf16_callback_probe.cpp"

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) && !std::strcmp(argv[1], "--check") &&
                (!std::strcmp(argv[2], "auto") || !std::strcmp(argv[2], "gfni")),
                "check-only observer: --check auto|gfni [parity_file]");
        callback_probe::Reset();
        const int result = CallbackWorkloadMain(argc, argv);
        if (result) return result;
        Require(callback_probe::state.pass_count == 1 && callback_probe::state.calls > 0,
                "expected exactly one GF16 encode pass");
        callback_probe::Print();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "R199 observer: %s\n", error.what());
        return 1;
    }
}
