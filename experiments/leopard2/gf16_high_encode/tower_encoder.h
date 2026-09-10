// Experiment-only large high-rate GF16 integration; leopard-79h.38.5.4.18.4.
#pragma once
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"

namespace tower_encoder {
struct Counts {
    uint64_t selected_passes, source_rows, source_bytes, output_rows, output_bytes;
    uint64_t inverse_pairs, forward_pairs, accumulating_pairs;
};
// Set only while encode/decode calls are quiescent. Default is OFF.
void SetEnabled(bool enabled);
bool TraceAvailable();
void ResetCounts();
Counts GetCounts();
unsigned InitializationCount();

// Keeps the public codec/ISA contract intact. Selected passes must use the
// existing large copy-first and non-fused, dense-prefix transform schedule.
bool Select(const leopard::backend::Ops& ops, uint64_t bytes,
            uint64_t source_policy_bytes, unsigned side,
            const leopard2_internal::SparseForwardPlanBatchView* sparse);
const leopard::backend::Ops& GetOps(const leopard::backend::Ops& original);
void CopySource(const leopard::backend::Ops& ops, void* destination,
                const void* source, uint64_t bytes);
void Finish(void** work, unsigned recovery_prefix, uint64_t bytes);
} // namespace tower_encoder
