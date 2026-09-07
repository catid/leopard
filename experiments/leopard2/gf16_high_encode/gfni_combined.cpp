#include "gfni_combined.h"
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace gfni_combined {
static State state = {};
bool Matches(const Call& c)
{
    return c.kind == LEO2_BACKEND_GFNI && c.k == 1000 && c.r == 200 &&
        c.requested == 200 && c.side == 256 && c.sparse_blocks == 0 &&
        c.sparse_present && c.bytes == 32768 && c.source_policy == 65536;
}
void Reset(unsigned mode)
{
    if (mode > 3) throw std::runtime_error("combined mode out of range");
    state = State{}; state.mode = mode;
}
unsigned ParseMode(const char* argument)
{
    if (!argument || std::strlen(argument) != 8 || std::strncmp(argument,"--mode=",7) != 0 ||
        argument[7] < '0' || argument[7] > '3')
        throw std::runtime_error("expected --mode=0|1|2|3");
    return static_cast<unsigned>(argument[7]-'0');
}
const State& Get() { return state; }
void Print()
{
    std::fprintf(stderr,"{\"schema\":\"gfni-combined/v1\",\"timed\":false,"
        "\"mode\":%u,\"calls\":%u,\"matches\":%u,\"first\":%u,\"terminal\":%u,\"records\":[",
        state.mode,state.calls,state.matches,state.first,state.terminal);
    for (unsigned i=0; i<state.calls; ++i)
    {
        const Call& c = state.records[i];
        std::fprintf(stderr,"%s{\"kind\":%u,\"k\":%u,\"r\":%u,\"requested\":%u,"
            "\"side\":%u,\"sparse_blocks\":%u,\"bytes\":%llu,\"source_policy\":%llu,\"sparse_present\":%s}",
            i ? "," : "",c.kind,c.k,c.r,c.requested,c.side,c.sparse_blocks,
            static_cast<unsigned long long>(c.bytes),static_cast<unsigned long long>(c.source_policy),
            c.sparse_present ? "true" : "false");
    }
    std::fprintf(stderr,"]}\n");
}
}

extern "C" unsigned LeoGFNICombinedExperiment(
    const leopard::backend::Ops& ops,uint64_t bytes,uint64_t policy,
    unsigned k,unsigned r,unsigned requested,unsigned side,
    const leopard2_internal::SparseForwardPlanBatchView* sparse)
{
    using namespace gfni_combined;
    if (state.calls >= kCapacity) throw std::runtime_error("combined pass limit");
    const Call call = {static_cast<unsigned>(ops.kind),k,r,requested,side,
        sparse ? sparse->block_count : 0U,bytes,policy,sparse != NULL};
    const bool matched = Matches(call);
    const unsigned mask = matched ? state.mode : 0U;
    state.records[state.calls++] = call;
    state.matches += matched;
    state.first += (mask & kFirst) != 0;
    state.terminal += (mask & kTerminal) != 0;
    return mask;
}
