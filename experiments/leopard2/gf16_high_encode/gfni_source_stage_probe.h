#pragma once
// Diagnostic-only policy contrast for leopard-79h.38.5.4.11.
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include <cstdint>

namespace gfni_source_stage_probe {
struct Call
{
    unsigned kind, k, r, requested, side, sparse_blocks;
    uint64_t bytes, source_policy, effective_policy;
};
struct State
{
    bool enabled;
    unsigned calls, matches, changed;
    Call records[16];
};
bool Matches(const Call& call);
void Reset(bool enabled);
const State& Get();
}
