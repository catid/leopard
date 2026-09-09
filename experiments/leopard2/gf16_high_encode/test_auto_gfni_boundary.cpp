// leopard-79h.38.5.4.17.1: focused tests of the default-off AUTO extension.
#define LEO_BOUNDARY_GUARD_NO_MAIN 1
#include "test_gfni_boundary.cpp"
#include "Leopard2Backend.h"
#include <atomic>
#include <thread>

namespace {
bool reject_host = false;
unsigned qualification_fault = 0, gfni_requests = 0;
}

extern "C" bool RealHost()
    asm("__real__ZN7leopard7backend34IsCalibratedAutoGF16GFNIEncodeHostEv");
extern "C" bool WrappedHost()
    asm("__wrap__ZN7leopard7backend34IsCalibratedAutoGF16GFNIEncodeHostEv");
extern "C" bool WrappedHost() { return !reject_host && RealHost(); }
extern "C" const leopard::backend::Ops* RealQualified(leo2_backend,
    leopard::backend::QualificationStatus*)
    asm("__real__ZN7leopard7backend15GetQualifiedOpsE12leo2_backendPNS0_19QualificationStatusE");
extern "C" const leopard::backend::Ops* WrappedQualified(leo2_backend,
    leopard::backend::QualificationStatus*)
    asm("__wrap__ZN7leopard7backend15GetQualifiedOpsE12leo2_backendPNS0_19QualificationStatusE");
extern "C" const leopard::backend::Ops* WrappedQualified(leo2_backend backend,
    leopard::backend::QualificationStatus* status)
{
    if (backend == LEO2_BACKEND_GFNI)
    {
        ++gfni_requests;
        if (qualification_fault)
        {
            if (status) *status = static_cast<leopard::backend::QualificationStatus>(
                qualification_fault);
            return NULL;
        }
    }
    return RealQualified(backend, status);
}

namespace {
namespace diag = leopard2_internal;

struct Spec
{
    unsigned k = 1000, r = 200, threads = 1, flags = 0;
    size_t bytes = 32768;
    leo2_backend backend = LEO2_BACKEND_AUTO;
    leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1;
    leo2_field field = LEO2_FIELD_GF16;
    leo2_shard_layout layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
};

Spec Target(unsigned cell)
{
    Require(cell < 2, "target cell");
    Spec result;
    if (cell) { result.r = 199; result.bytes = 65536; }
    return result;
}

struct Config
{
    Spec spec;
    Codec::Context context;
    Codec::Handle codec;
    explicit Config(Spec s) : spec(s), context(NULL, leo2_context_destroy),
        codec(NULL, leo2_codec_destroy)
    {
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = s.backend;
        options.thread_count = s.threads;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "test context");
        context.reset(raw_context);
        Require(leo2_context_backend(context.get()) ==
            (s.backend == LEO2_BACKEND_AUTO ? LEO2_BACKEND_AVX2 : s.backend), "exact backend");
        leo2_codec_options codec_options = {};
        codec_options.struct_size = sizeof(codec_options);
        codec_options.flags = s.flags;
        codec_options.shard_layout = s.layout;
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), s.k, s.r, s.profile, s.field,
            &codec_options, &raw_codec) == LEO2_SUCCESS, "test codec");
        codec.reset(raw_codec);
    }
    bool Selected(size_t bytes) const
    {
        return diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), bytes);
    }
};

struct Work
{
    const Config& config;
    Guard source, output, scratch;
    std::vector<const void*> inputs;
    std::vector<void*> outputs;
    uint64_t input_hash;
    static size_t Scratch(const Config& c)
    {
        size_t bytes = 0;
        Require(leo2_encode_scratch_size(c.codec.get(), c.spec.bytes, &bytes) == LEO2_SUCCESS,
            "test scratch query");
        return bytes;
    }
    explicit Work(const Config& c) : config(c), source(37 * c.spec.bytes),
        output(c.spec.r * c.spec.bytes), scratch(Scratch(c)),
        inputs(c.spec.k), outputs(c.spec.r), input_hash(0)
    {
        uint32_t random = 20260909;
        Fill(source.data, source.bytes, random);
        input_hash = Hash(source.data, source.bytes);
        for (unsigned i = 0; i < c.spec.k; ++i)
            inputs[i] = source.data + (i % 37) * c.spec.bytes;
        for (unsigned i = 0; i < c.spec.r; ++i)
            outputs[i] = output.data + i * c.spec.bytes;
    }
    void Encode()
    {
        Require(leo2_encode(config.codec.get(), config.spec.bytes, inputs.data(), outputs.data(),
            scratch.data, scratch.bytes) == LEO2_SUCCESS, "public encode");
    }
    leo2_encode_batch_item Item()
    {
        leo2_encode_batch_item item = {};
        item.shard_bytes = config.spec.bytes;
        item.original = inputs.data(); item.recovery = outputs.data();
        item.scratch = scratch.data; item.scratch_bytes = scratch.bytes;
        return item;
    }
    void Equal(const Work& reference)
    {
        Require(output.bytes == reference.output.bytes &&
            std::memcmp(output.data, reference.output.data, output.bytes) == 0, "full parity");
        Require(Hash(source.data, source.bytes) == input_hash, "input mutation");
        source.Check(); output.Check(); scratch.Check();
    }
};

void Count(unsigned count)
{
    Require(diag::AutoGF16GFNIEncodeCallCountForDiagnostics() == count, "AUTO route call count");
}

void Routes()
{
    Require(diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1, "old production route disabled");
    Config old(Target(0)), disabled(Target(1));
    Require(old.Selected(65536) && !old.Selected(32768), "default-off changed old R200 route");
    Require(!disabled.Selected(65536) &&
        !diag::AutoGF16GFNIEncodeAvailableForDiagnostics(disabled.codec.get()),
        "disabled R199 qualified optional table");
    Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(true), "enable boundaries");
    Require(!disabled.Selected(65536), "late enable bypassed cached qualification");
    Config added(Target(1));
    Require(old.Selected(32768) && old.Selected(65536) && added.Selected(65536) &&
        !added.Selected(32768), "exact boundary selector");
    for (size_t bytes : {size_t(32766), size_t(32770), size_t(65534), size_t(65538), size_t(131072)})
        Require(!old.Selected(bytes) && !added.Selected(bytes), "byte neighbor widened");
    for (unsigned cell = 0; cell < 2; ++cell)
    {
        for (unsigned variant = 0; variant < 12; ++variant)
        {
            Spec s = Target(cell);
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
                "context/codec gate escaped");
        }
    }
    // Negative processor identities use the real unchanged calibration predicate.
    for (leopard::backend::X86ProcessorIdentity id : {
            leopard::backend::X86ProcessorIdentity{false, 0x1a, 8},
            leopard::backend::X86ProcessorIdentity{true, 0x19, 8},
            leopard::backend::X86ProcessorIdentity{true, 0x1a, 7},
            leopard::backend::X86ProcessorIdentity{true, 0x1a, 9}})
        Require(!leopard::backend::IsCalibratedAutoGF16GFNIEncodeProcessor(id), "host gate escaped");
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(false), "global off");
    Require(!old.Selected(32768) && !old.Selected(65536) && !added.Selected(65536),
        "boundary control bypassed global off");
    Config globally_disabled(Target(1));
    Require(!diag::AutoGF16GFNIEncodeAvailableForDiagnostics(globally_disabled.codec.get()),
        "global-off optional qualification");
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish global off");
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "restore global on");
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish global on");
    Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(false), "restore boundaries off");
    Require(!old.Selected(32768) && old.Selected(65536) && !added.Selected(65536),
        "boundary-off changed old target");
    Count(0);
}

void Api(unsigned cell)
{
    Spec spec = Target(cell), avx2_spec = spec;
    avx2_spec.backend = LEO2_BACKEND_AVX2;
    Config candidate(spec), baseline(avx2_spec);
    Require(candidate.Selected(spec.bytes), "candidate unavailable");
    Work reference(baseline), work(candidate);
    reference.Encode();
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "arm route accounting");
    work.Encode(); work.Equal(reference); Count(1);
    work.outputs.back() = NULL;
    std::memset(work.output.data, 0x5a, work.output.bytes);
    work.Encode(); Count(1);
    Require(std::memcmp(work.output.data, reference.output.data,
        (spec.r - 1) * spec.bytes) == 0, "partial parity");
    for (size_t i = (spec.r - 1) * spec.bytes; i < work.output.bytes; ++i)
        Require(work.output.data[i] == 0x5a, "omitted output changed");
    work.outputs.back() = work.output.data + (spec.r - 1) * spec.bytes;
    leo2_encode_batch_item item = work.Item();
    Require(leo2_encode_batch(candidate.codec.get(), &item, 1) == LEO2_SUCCESS, "one-item batch");
    work.Equal(reference); Count(2);
    Require(leo2_encode_batch_with_preflight_scratch(candidate.codec.get(), &item, 1,
        NULL, 0) == LEO2_SUCCESS, "scalable one-item alias");
    work.Equal(reference); Count(2);
    leo2_encode_batch_binding* raw_binding = NULL;
    Require(leo2_encode_batch_binding_create(candidate.codec.get(), &item, 1,
        &raw_binding) == LEO2_SUCCESS, "binding create");
    std::unique_ptr<leo2_encode_batch_binding, decltype(&leo2_encode_batch_binding_destroy)>
        binding(raw_binding, leo2_encode_batch_binding_destroy);
    Require(leo2_encode_batch_binding_execute(binding.get()) == LEO2_SUCCESS, "binding execute");
    work.Equal(reference); Count(2);
    {
        Work second(candidate);
        leo2_encode_batch_item items[] = {work.Item(), second.Item()};
        Require(leo2_encode_batch(candidate.codec.get(), items, 2) == LEO2_SUCCESS, "two-item batch");
        work.Equal(reference); second.Equal(reference); Count(2);
        size_t preflight_bytes = 0;
        Require(leo2_encode_batch_preflight_scratch_size(candidate.codec.get(), 2,
            &preflight_bytes) == LEO2_SUCCESS && preflight_bytes > 0, "scalable query");
        Guard preflight(preflight_bytes);
        Require(leo2_encode_batch_with_preflight_scratch(candidate.codec.get(), items, 2,
            preflight.data, preflight.bytes) == LEO2_SUCCESS, "scalable two-item batch");
        work.Equal(reference); second.Equal(reference); preflight.Check(); Count(2);
        items[1].scratch_bytes = 0;
        std::memset(work.output.data, 0x5a, work.output.bytes);
        Require(leo2_encode_batch(candidate.codec.get(), items, 2) == LEO2_SCRATCH_TOO_SMALL,
            "invalid batch accepted");
        for (size_t i = 0; i < work.output.bytes; ++i)
            Require(work.output.data[i] == 0x5a, "invalid batch not atomic");
        Count(2);
    }
    work.Encode(); work.Equal(reference); Count(3);
    std::vector<uint8_t> originals(spec.k, 1), recovery(spec.r, 1);
    originals[0] = 0;
    leo2_decode_plan* raw_plan = NULL;
    Require(leo2_decode_plan_create(candidate.codec.get(), originals.data(), recovery.data(),
        &raw_plan) == LEO2_SUCCESS, "decode plan");
    std::unique_ptr<leo2_decode_plan, decltype(&leo2_decode_plan_destroy)>
        plan(raw_plan, leo2_decode_plan_destroy);
    size_t decode_bytes = 0;
    Require(leo2_decode_plan_scratch_size(plan.get(), spec.bytes, &decode_bytes) == LEO2_SUCCESS,
        "decode scratch");
    Guard decode_scratch(decode_bytes), restored(spec.bytes);
    std::vector<const void*> decode_inputs = work.inputs;
    decode_inputs[0] = NULL;
    std::vector<const void*> parity(work.outputs.begin(), work.outputs.end());
    std::vector<void*> restored_outputs(spec.k, NULL);
    restored_outputs[0] = restored.data;
    Require(leo2_decode_plan_execute(plan.get(), spec.bytes, decode_inputs.data(), parity.data(),
        restored_outputs.data(), decode_scratch.data, decode_scratch.bytes) == LEO2_SUCCESS,
        "baseline decode");
    Require(std::memcmp(restored.data, work.source.data, spec.bytes) == 0, "decode round trip");
    decode_scratch.Check(); restored.Check(); Count(3);
    Require(leo2_context_backend(candidate.context.get()) == LEO2_BACKEND_AVX2, "baseline changed");
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "normalize route probe");
}

void Concurrent(unsigned cell)
{
    Spec spec = Target(cell), avx2_spec = spec;
    avx2_spec.backend = LEO2_BACKEND_AVX2;
    Config candidate(spec), baseline(avx2_spec);
    Work reference(baseline);
    reference.Encode();
    std::atomic<unsigned> failures(0);
    std::thread first([&]() { try { Work work(candidate); work.Encode(); work.Equal(reference); }
                            catch (...) { ++failures; } });
    std::thread second([&]() { try { Work work(candidate); work.Encode(); work.Equal(reference); }
                             catch (...) { ++failures; } });
    first.join(); second.join();
    Require(failures.load() == 0 && candidate.Selected(spec.bytes), "concurrent immutable codec");
}

void Fault(unsigned cell, const char* kind)
{
    if (!std::strcmp(kind, "host")) reject_host = true;
    else if (!std::strcmp(kind, "unavailable")) qualification_fault = 1;
    else if (!std::strcmp(kind, "oom")) qualification_fault = 2;
    else if (!std::strcmp(kind, "kat")) qualification_fault = 3;
    else Require(false, "fault kind");
    Spec spec = Target(cell), avx2_spec = spec;
    avx2_spec.backend = LEO2_BACKEND_AVX2;
    Config candidate(spec), baseline(avx2_spec);
    Require(!candidate.Selected(spec.bytes) &&
        !diag::AutoGF16GFNIEncodeAvailableForDiagnostics(candidate.codec.get()), "failed qualification widened");
    Require(gfni_requests == (reject_host ? 0U : 1U), "optional qualification request count");
    Work reference(baseline), work(candidate);
    reference.Encode();
    Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "fallback probe");
    work.Encode(); work.Equal(reference); Count(0);
    Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish fallback probe");
}
}

int main(int argc, char** argv)
{
    try
    {
        Require(argc >= 2 && !diag::AutoGF16GFNIBoundariesEnabledForDiagnostics(),
            "candidate must default off");
        if (argc == 2 && !std::strcmp(argv[1], "--routes")) Routes();
        else
        {
            Require(argc >= 3 && std::strlen(argv[2]) == 1 && argv[2][0] >= '0' &&
                argv[2][0] <= '7', "cell");
            const unsigned cell = static_cast<unsigned>(argv[2][0] - '0');
            bool enabled = true;
            if (!std::strcmp(argv[1], "--guards"))
            {
                Require(argc == 4 && (!std::strcmp(argv[3], "0") || !std::strcmp(argv[3], "1")),
                    "guard mode");
                enabled = argv[3][0] == '1';
            }
            Require(diag::SetAutoGF16GFNIBoundariesEnabledForDiagnostics(enabled), "set boundaries");
            if (!std::strcmp(argv[1], "--guards"))
            {
                Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "guard probe");
                Check(cell, LEO2_BACKEND_AUTO);
                Count(enabled && cell < 4 ? 1 : 0);
                Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(), "finish guard probe");
            }
            else if (argc == 3 && !std::strcmp(argv[1], "--api")) Api(cell);
            else if (argc == 3 && !std::strcmp(argv[1], "--concurrent")) Concurrent(cell);
            else if (argc == 4 && !std::strcmp(argv[1], "--fault")) Fault(cell, argv[3]);
            else Require(false, "usage: --routes | --api CELL | --guards CELL 0|1 | --concurrent CELL | --fault CELL host|unavailable|oom|kat");
        }
        std::printf("{\"schema\":\"leopard-auto-gfni-boundary-check/v1\",\"case\":\"%s\","
            "\"passed\":true,\"timed\":false}\n", argv[1]);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
