#pragma once
// Four-state, default-off experiment for leopard-79h.38.5.4.16.
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include <cstdint>

namespace gfni_combined {
#ifndef LEO_GFNI_COMBINED_CAPACITY
#define LEO_GFNI_COMBINED_CAPACITY 16
#endif
static_assert(LEO_GFNI_COMBINED_CAPACITY == 16 || LEO_GFNI_COMBINED_CAPACITY == 64,
    "unsupported experiment trace capacity");
constexpr unsigned kCapacity = LEO_GFNI_COMBINED_CAPACITY;
constexpr unsigned kFirst = 1U, kTerminal = 2U;
struct Call {
    unsigned kind, k, r, requested, side, sparse_blocks;
    uint64_t bytes, source_policy;
    bool sparse_present;
};
struct State {
    unsigned mode, calls, matches, first, terminal;
    Call records[kCapacity];
};
bool Matches(const Call& call);
void Reset(unsigned mode);
unsigned ParseMode(const char* argument);
const State& Get();
void Print();
}
extern "C" unsigned LeoGFNICombinedExperiment(
    const leopard::backend::Ops&, uint64_t, uint64_t,
    unsigned, unsigned, unsigned, unsigned,
    const leopard2_internal::SparseForwardPlanBatchView*);

// Unchanged .14 kernel: disjoint rows/accumulators, complete even symbols,
// qualified GFNI. Empty byte/distance operations do not dereference pointers.
extern "C" void LeoGFNIFinalAccumulate(const void* const*,void* const*,unsigned,
    uint16_t,uint16_t,uint16_t,uint64_t);
