// Link-only contrast; neither codec archive is edited or rebuilt.
#include "gfni_source_stage_probe.h"
#include <cstdio>
#include <stdexcept>

namespace gfni_source_stage_probe {
static State state = {};

bool Matches(const Call& call)
{
    return call.kind == LEO2_BACKEND_GFNI && call.k == 1000 && call.r == 200 &&
        call.requested == 200 && call.side == 256 && call.sparse_blocks == 0 &&
        call.bytes == 32768 && call.source_policy == 65536;
}

void Reset(bool enabled)
{
    state = State{};
    state.enabled = enabled;
}

const State& Get() { return state; }
}

// Exact GCC/Itanium name observed in the pinned dual-field archive. Build
// checks require this symbol; no substring or guessed-symbol interception.
#define LEO_STAGE_SYMBOL "_ZN7leopard4ff1633ReedSolomonEncodeWithSourcePolicyERKNS_7backend3OpsEmmjjjjPKPKvPPvPKN17leopard2_internal26SparseForwardPlanBatchViewE"
#define LEO_STAGE_ARGS const leopard::backend::Ops& ops, uint64_t bytes, \
    uint64_t policy, unsigned k, unsigned r, unsigned requested, unsigned side, \
    const void* const* data, void** work, \
    const leopard2_internal::SparseForwardPlanBatchView* sparse

extern "C" void RealSourceStage(LEO_STAGE_ARGS) asm("__real_" LEO_STAGE_SYMBOL);
extern "C" void WrappedSourceStage(LEO_STAGE_ARGS) asm("__wrap_" LEO_STAGE_SYMBOL);

extern "C" void WrappedSourceStage(LEO_STAGE_ARGS)
{
    using namespace gfni_source_stage_probe;
    if (state.calls >= kCapacity)
        throw std::runtime_error("source-stage probe record limit");
    Call call = {static_cast<unsigned>(ops.kind), k, r, requested, side,
        sparse ? sparse->block_count : 0U, bytes, policy, policy};
    if (Matches(call))
    {
        ++state.matches;
        if (state.enabled)
        {
            call.effective_policy = 16384;
            ++state.changed;
        }
    }
    state.records[state.calls++] = call;
    RealSourceStage(ops, bytes, call.effective_policy, k, r, requested, side,
        data, work, sparse);
}

#undef LEO_STAGE_ARGS
#undef LEO_STAGE_SYMBOL

#ifdef LEO_GFNI_SOURCE_STAGE_MAIN
#ifndef LEO_GFNI_SOURCE_STAGE_ENABLE
#define LEO_GFNI_SOURCE_STAGE_ENABLE 0
#endif
static_assert(LEO_GFNI_SOURCE_STAGE_ENABLE == 0 ||
              LEO_GFNI_SOURCE_STAGE_ENABLE == 1, "invalid experiment mode");
extern "C" int __real_main(int, char**);
extern "C" int __wrap_main(int argc, char** argv)
{
    using namespace gfni_source_stage_probe;
    Reset(LEO_GFNI_SOURCE_STAGE_ENABLE != 0);
    const int result = __real_main(argc, argv);
    if (result != 0) return result;
    const State& observed = Get();
    if (observed.calls == 0)
    {
        std::fprintf(stderr, "source-stage wrapper was not reached\n");
        return 1;
    }
    std::fprintf(stderr,
        "{\"schema\":\"gfni-source-stage/v1\",\"enabled\":%s,\"matches\":%u,"
        "\"changed\":%u,\"calls\":[", observed.enabled ? "true" : "false",
        observed.matches, observed.changed);
    for (unsigned i = 0; i < observed.calls; ++i)
    {
        const Call& call = observed.records[i];
        std::fprintf(stderr,
            "%s{\"kind\":%u,\"k\":%u,\"r\":%u,\"requested\":%u,\"side\":%u,"
            "\"sparse_blocks\":%u,\"bytes\":%llu,\"source_policy\":%llu,"
            "\"effective_policy\":%llu}", i ? "," : "", call.kind, call.k,
            call.r, call.requested, call.side, call.sparse_blocks,
            static_cast<unsigned long long>(call.bytes),
            static_cast<unsigned long long>(call.source_policy),
            static_cast<unsigned long long>(call.effective_policy));
    }
    std::fprintf(stderr, "],\"timed\":false}\n");
    return 0;
}
#endif
