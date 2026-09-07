#include "gfni_terminal.h"
#include <cstdio>
#include <stdexcept>

namespace gfni_terminal {
static State state = {};
bool Matches(const Call& c)
{
    return c.kind == LEO2_BACKEND_GFNI && c.k == 1000 && c.r == 200 &&
        c.requested == 200 && c.side == 256 && c.sparse_blocks == 0 &&
        c.sparse_present && c.bytes == 32768 && c.source_policy == 65536;
}
void Reset(bool enabled) { state = State{}; state.enabled = enabled; }
const State& Get() { return state; }
void Print()
{
    std::fprintf(stderr, "{\"schema\":\"gfni-terminal/v1\",\"timed\":false,"
        "\"enabled\":%s,\"calls\":%u,\"matches\":%u,\"changed\":%u,\"records\":[",
        state.enabled ? "true" : "false",state.calls,state.matches,state.changed);
    for (unsigned i=0; i<state.calls; ++i)
    {
        const Call& c = state.records[i];
        std::fprintf(stderr, "%s{\"kind\":%u,\"k\":%u,\"r\":%u,\"requested\":%u,"
            "\"side\":%u,\"sparse_blocks\":%u,\"bytes\":%llu,\"source_policy\":%llu,\"sparse_present\":%s}",
            i ? "," : "",c.kind,c.k,c.r,c.requested,c.side,c.sparse_blocks,
            static_cast<unsigned long long>(c.bytes),static_cast<unsigned long long>(c.source_policy),
            c.sparse_present ? "true" : "false");
    }
    std::fprintf(stderr, "]}\n");
}
}

extern "C" bool LeoGFNITerminalExperiment(
    const leopard::backend::Ops& ops, uint64_t bytes, uint64_t policy,
    unsigned k, unsigned r, unsigned requested, unsigned side,
    const leopard2_internal::SparseForwardPlanBatchView* sparse)
{
    using namespace gfni_terminal;
    if (state.calls >= kCapacity) throw std::runtime_error("terminal pass limit");
    const Call call = {static_cast<unsigned>(ops.kind),k,r,requested,side,
        sparse ? sparse->block_count : 0U,bytes,policy,sparse != NULL};
    const bool matched = Matches(call);
    const bool changed = state.enabled && matched;
    state.records[state.calls++] = call;
    state.matches += matched;
    state.changed += changed;
    return changed;
}
