#pragma once
// Default-off experiment for leopard-79h.38.5.4.13, not a production interface.
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include <cstdint>

namespace gfni_inverse_distance1 {
struct Call {
    unsigned kind, k, r, requested, side, sparse_blocks;
    uint64_t bytes, source_policy;
    bool sparse_present;
};
struct State {
    bool enabled;
    unsigned calls, matches, changed;
    Call records[16];
};
bool Matches(const Call& call);
void Reset(bool enabled);
const State& Get();
void Print();
}
extern "C" bool LeoGFNIInverseDistance1Experiment(
    const leopard::backend::Ops&, uint64_t, uint64_t,
    unsigned, unsigned, unsigned, unsigned,
    const leopard2_internal::SparseForwardPlanBatchView*);
