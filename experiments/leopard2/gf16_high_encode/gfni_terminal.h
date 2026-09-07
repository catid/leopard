#pragma once
// Default-off experiment for leopard-79h.38.5.4.14, not a production interface.
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include <cstdint>

namespace gfni_terminal {
#ifndef LEO_GFNI_TERMINAL_CAPACITY
#define LEO_GFNI_TERMINAL_CAPACITY 16
#endif
static_assert(LEO_GFNI_TERMINAL_CAPACITY == 16 || LEO_GFNI_TERMINAL_CAPACITY == 64,
    "unsupported experiment trace capacity");
constexpr unsigned kCapacity = LEO_GFNI_TERMINAL_CAPACITY;
struct Call {
    unsigned kind, k, r, requested, side, sparse_blocks;
    uint64_t bytes, source_policy;
    bool sparse_present;
};
struct State {
    bool enabled;
    unsigned calls, matches, changed;
    Call records[kCapacity];
};
bool Matches(const Call& call);
void Reset(bool enabled);
const State& Get();
void Print();
}
extern "C" bool LeoGFNITerminalExperiment(
    const leopard::backend::Ops&, uint64_t, uint64_t,
    unsigned, unsigned, unsigned, unsigned,
    const leopard2_internal::SparseForwardPlanBatchView*);

// Sources/accumulators and all four-distance coordinates must be disjoint.
// Complete even GF16 symbols; zero bytes touch no pointers. GFNI must be qualified.
extern "C" void LeoGFNIFinalAccumulate(const void* const*,void* const*,unsigned,
    uint16_t,uint16_t,uint16_t,uint64_t);
