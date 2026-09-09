// Focused exact AUTO addition; leopard-79h.38.5.4.19.1. No benchmark clocks.
#define LEO_AUTO_GFNI_R19932_TARGET 1
#define main OriginalBoundaryTestMain
#include "test_auto_gfni_boundary.cpp"
#undef main

namespace {
void NewRoutes()
{
    Require(diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1 &&
            diag::AutoGF16GFNIBoundariesEnabledForDiagnostics(), "old defaults changed");
    Spec target = Target(0), old200 = target;
    old200.r = 200;
    Config disabled(target), established(old200);
    Require(!disabled.Selected(32768) && disabled.Selected(65536) &&
            established.Selected(32768) && established.Selected(65536), "default routes changed");
    Require(diag::AutoGF16GFNIEncodeAvailableForDiagnostics(disabled.codec.get()),
            "R199 cached table for existing 64KiB route disappeared");
    Require(diag::SetAutoGF16GFNIR19932EnabledForDiagnostics(true), "enable new cell");
    Config enabled(target);
    Require(enabled.Selected(32768) && disabled.Selected(32768) &&
            enabled.Selected(65536) && established.Selected(32768) && established.Selected(65536),
            "new or established route mismatch");
    for (size_t bytes : {size_t(0), size_t(64), size_t(32766), size_t(32770),
                        size_t(65534), size_t(65538), size_t(131072)})
        Require(!enabled.Selected(bytes), "byte neighbor widened");
    for (unsigned variant = 0; variant < 12; ++variant)
    {
        Spec s = target;
        switch (variant)
        {
        case 0: s.k = 999; break;
        case 1: s.k = 1001; break;
        case 2: s.r = 198; break;
        case 3: s.r = 201; break;
        case 4: s.threads = 2; break;
        case 5: s.flags = LEO2_CODEC_FORCE_SPECIALIZED_DECODE; break;
        case 6: s.layout = LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1; break;
        case 7: s.profile = LEO2_PROFILE_LOW_V1; break;
        case 8: s.backend = LEO2_BACKEND_AVX2; break;
        case 9: s.backend = LEO2_BACKEND_GFNI; break;
        case 10: s.backend = LEO2_BACKEND_AVX512; break;
        case 11: s.k = 17; s.r = 7; s.field = LEO2_FIELD_GF8; break;
        }
        Config negative(s);
        Require(!negative.Selected(s.bytes) &&
                !diag::AutoGF16GFNIEncodeAvailableForDiagnostics(negative.codec.get()),
                "context or codec exclusion bypassed");
    }
    Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(false), "old control off");
    Config no_table(target);
    Require(!enabled.Selected(32768) && !enabled.Selected(65536) &&
            !established.Selected(32768) && established.Selected(65536) &&
            !diag::AutoGF16GFNIEncodeAvailableForDiagnostics(no_table.codec.get()), "old control bypassed");
    Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(true), "old control on");
    Require(!no_table.Selected(32768) && !no_table.Selected(65536), "late enable bypassed cache");
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(false), "global off");
    Config globally_disabled(target);
    Require(!enabled.Selected(32768) && !enabled.Selected(65536) &&
            !diag::AutoGF16GFNIEncodeAvailableForDiagnostics(globally_disabled.codec.get()),
            "global control bypassed");
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish global off");
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "global on");
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish global on");
    Require(!globally_disabled.Selected(32768), "global late enable bypassed cache");
    Require(diag::SetAutoGF16GFNIR19932EnabledForDiagnostics(false), "restore new control off");
    Require(!enabled.Selected(32768) && enabled.Selected(65536) &&
            established.Selected(32768) && established.Selected(65536), "new control changed old cases");
    Count(0);
}

void NewGuards(unsigned cell, bool enabled)
{
    Require(cell < 8, "guarded cell");
    const size_t bytes[] = {32768, 32768, 32770, 32766, 64, 66, 65, 66};
    const size_t offsets[] = {0, 1, 2, 1, 1, 2, 1, 1};
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "arm guards");
    CheckShape(cell, cell < 6 ? 1000 : 17, cell < 6 ? 199 : 7,
               bytes[cell], offsets[cell], cell == 6 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16,
               LEO2_BACKEND_AUTO, LEO2_BACKEND_AVX2);
    Count(enabled && cell < 2 ? 1 : 0);
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish guards");
}

void MixedConcurrency()
{
    Spec spec = Target(0), avx2 = spec, small = spec;
    avx2.backend = LEO2_BACKEND_AVX2;
    small.k = 17; small.r = 7; small.bytes = 65; small.field = LEO2_FIELD_GF8;
    Config candidate(spec), baseline(avx2), gf8(small);
    Work reference(baseline), small_reference(gf8);
    reference.Encode(); small_reference.Encode();
    std::atomic<unsigned> failures(0);
    const auto run = [&](const Config& codec, const Work& expected) {
        try { Work work(codec); for (unsigned i = 0; i < 4; ++i) { work.Encode(); work.Equal(expected); } }
        catch (...) { ++failures; }
    };
    std::thread a([&]() { run(candidate, reference); }), b([&]() { run(candidate, reference); });
    std::thread c([&]() { run(gf8, small_reference); }), d([&]() { run(gf8, small_reference); });
    a.join(); b.join(); c.join(); d.join();
    Require(failures.load() == 0 && candidate.Selected(32768), "mixed-field immutable concurrency");
}
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc >= 2 && !diag::AutoGF16GFNIR19932EnabledForDiagnostics(), "new default must be off");
        if (argc == 2 && !std::strcmp(argv[1], "--routes")) NewRoutes();
        else
        {
            bool enabled = true;
            if (!std::strcmp(argv[1], "--guards"))
            {
                Require(argc == 4 && std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '7' &&
                        (!std::strcmp(argv[3], "0") || !std::strcmp(argv[3], "1")), "guard CLI");
                enabled = argv[3][0] == '1';
            }
            Require(diag::SetAutoGF16GFNIR19932EnabledForDiagnostics(enabled), "set new control");
            if (argc == 2 && !std::strcmp(argv[1], "--api")) Api(0);
            else if (argc == 2 && !std::strcmp(argv[1], "--concurrent")) MixedConcurrency();
            else if (argc == 3 && !std::strcmp(argv[1], "--fault")) Fault(0, argv[2]);
            else if (!std::strcmp(argv[1], "--guards")) NewGuards(argv[2][0]-'0', enabled);
            else Require(false, "usage: --routes | --api | --concurrent | --fault host|unavailable|oom|kat | --guards cell 0|1");
        }
        std::printf("{\"schema\":\"leopard-auto-r19932-check/v1\",\"case\":\"%s\",\"passed\":true,\"timed\":false}\n", argv[1]);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
