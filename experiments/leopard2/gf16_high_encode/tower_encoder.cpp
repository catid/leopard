// Experiment-only. This TU uses portable instructions; AVX2 stays in the
// already-qualified tower_butterfly_probe TU. No public payload allocation.
#include "tower_encoder.h"
#include "tower_butterfly_probe.h"
#include "LeopardFF16.h"
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <mutex>

#ifndef LEO_TOWER_TRACE
#define LEO_TOWER_TRACE 0
#endif

namespace tower_encoder {
namespace {
std::atomic<bool> enabled{false};
std::once_flag initialization;
std::atomic<unsigned> initialization_count{0};
leopard::backend::Ops tower_ops{};
TowerProductTables tables[65536];
TowerLowMapTables map{};
#if LEO_TOWER_TRACE
static thread_local Counts counts{};
#define TOWER_COUNT(member, value) (counts.member += (value))
#else
#define TOWER_COUNT(member, value) ((void)0)
#endif

const TowerProductTables* table(uint16_t log)
{
    return log == 65535 ? nullptr : &tables[log];
}
void inverse(void* x, void* y, uint16_t log, uint64_t bytes)
{
    TOWER_COUNT(inverse_pairs, 1);
    tower_ifft_pair(static_cast<uint8_t*>(x), static_cast<uint8_t*>(y), bytes, log, table(log));
}
void forward(void* x, void* y, uint16_t log, uint64_t bytes)
{
    TOWER_COUNT(forward_pairs, 1);
    tower_fft_pair(static_cast<uint8_t*>(x), static_cast<uint8_t*>(y), bytes, log, table(log));
}
void accumulate(const void* x, const void* y, void* u, void* v,
                uint16_t log, uint64_t bytes)
{
    TOWER_COUNT(accumulating_pairs, 1);
    tower_ifft_accumulate(static_cast<const uint8_t*>(x), static_cast<const uint8_t*>(y),
                         static_cast<uint8_t*>(u), static_cast<uint8_t*>(v), bytes, log, table(log));
}
template<bool Inverse>
void range(void* const* work, unsigned distance, uint16_t log01, uint16_t log23,
           uint16_t log02, uint64_t bytes, bool prefer_fused)
{
    // The selection proof excludes compact/fused passes. Fail loudly if a
    // future scheduler change violates it, instead of using canonical math.
    if (prefer_fused) std::abort();
    for (unsigned i = 0; i < distance; ++i) {
        void* x0 = work[i]; void* x1 = work[i+distance];
        void* x2 = work[i+2*distance]; void* x3 = work[i+3*distance];
        if (Inverse) {
            inverse(x0, x1, log01, bytes); inverse(x2, x3, log23, bytes);
            inverse(x0, x2, log02, bytes); inverse(x1, x3, log02, bytes);
        } else {
            forward(x0, x2, log02, bytes); forward(x1, x3, log02, bytes);
            forward(x0, x1, log01, bytes); forward(x2, x3, log23, bytes);
        }
    }
}
void unexpected_multiply(void*, const void*, uint16_t, uint64_t) { std::abort(); }
void unexpected_four(void*, void*, void*, void*, uint16_t, uint16_t, uint16_t, uint64_t)
{
    std::abort();
}
void unexpected_two_out(const void*, const void*, void*, void*, uint16_t, uint64_t)
{
    std::abort();
}
void unexpected_four_out(const void*, const void*, const void*, const void*,
                         void*, void*, void*, void*, uint16_t, uint16_t, uint16_t, uint64_t)
{
    std::abort();
}
void initialize(const leopard::backend::Ops& original)
{
    if (original.kind != LEO2_BACKEND_AVX2) std::abort();
    uint8_t subfield[65536];
    uint16_t u_times[256];
    for (unsigned a = 0; a < 256; ++a) {
        u_times[a] = leopard::ff16::MultiplyElements(256, uint16_t(a));
        if ((u_times[a] >> 8) != a) std::abort();
        for (unsigned b = 0; b < 256; ++b) {
            const uint16_t value = leopard::ff16::MultiplyElements(uint16_t(a), uint16_t(b));
            if (value >= 256) std::abort();
            subfield[a*256+b] = uint8_t(value);
        }
    }
    if ((leopard::ff16::MultiplyElements(256, 256) ^ 256) != 128) std::abort();
    for (unsigned i = 0; i < 16; ++i) {
        map.row[0][i] = uint8_t(u_times[i]);
        map.row[1][i] = uint8_t(u_times[i << 4]);
    }
    for (unsigned log = 0; log < 65536; ++log) {
        // Preserve raw log65535 == multiplication by one in the table;
        // butterfly callbacks independently interpret it as zero skew.
        const uint16_t canonical = leopard::ff16::MultiplyLogElement(1, uint16_t(log));
        const uint8_t d = uint8_t(canonical >> 8);
        const uint8_t c = uint8_t(canonical) ^ uint8_t(u_times[d]);
        const uint8_t coefficients[3] = {c, subfield[128*256+d], uint8_t(c ^ d)};
        for (unsigned k = 0; k < 3; ++k)
            for (unsigned i = 0; i < 16; ++i) {
                tables[log].row[2*k][i] = subfield[i*256+coefficients[k]];
                tables[log].row[2*k+1][i] = subfield[(i << 4)*256+coefficients[k]];
            }
    }
    tower_ops = original;
    tower_ops.name = "experiment-only GF16 tower / AVX2";
    tower_ops.ff16_ifft_butterfly2 = inverse;
    tower_ops.ff16_fft_butterfly2 = forward;
    tower_ops.ff16_ifft_butterfly2_xor = accumulate;
    tower_ops.ff16_ifft_butterfly4_range = range<true>;
    tower_ops.ff16_fft_butterfly4_range = range<false>;
    // These are unreachable under Select, not supported fallback operations.
    tower_ops.ff16_multiply = unexpected_multiply;
    tower_ops.ff16_multiply_add = unexpected_multiply;
    tower_ops.ff16_fft_butterfly2_out = unexpected_two_out;
    tower_ops.ff16_ifft_butterfly4 = unexpected_four;
    tower_ops.ff16_fft_butterfly4 = unexpected_four;
    tower_ops.ff16_ifft_butterfly4_out = unexpected_four_out;
    tower_ops.ff16_fft_butterfly4_out = unexpected_four_out;
    initialization_count.fetch_add(1, std::memory_order_relaxed);
}
} // namespace

void SetEnabled(bool value) { enabled.store(value, std::memory_order_relaxed); }
bool TraceAvailable() { return LEO_TOWER_TRACE != 0; }
void ResetCounts()
{
#if LEO_TOWER_TRACE
    counts = {};
#endif
}
Counts GetCounts()
{
#if LEO_TOWER_TRACE
    return counts;
#else
    return {};
#endif
}
unsigned InitializationCount() { return initialization_count.load(std::memory_order_relaxed); }
bool Select(const leopard::backend::Ops& ops, uint64_t bytes,
            uint64_t source_policy_bytes, unsigned side,
            const leopard2_internal::SparseForwardPlanBatchView* sparse)
{
    const bool selected = enabled.load(std::memory_order_relaxed) &&
        ops.kind == LEO2_BACKEND_AVX2 && side >= 256 &&
        source_policy_bytes > 16U*1024U && bytes >= 256 && bytes % 64 == 0 &&
        (!sparse || sparse->block_count == 0);
    if (selected) TOWER_COUNT(selected_passes, 1);
    return selected;
}
const leopard::backend::Ops& GetOps(const leopard::backend::Ops& original)
{
    std::call_once(initialization, [&original] { initialize(original); });
    return tower_ops;
}
void CopySource(const leopard::backend::Ops& ops, void* destination,
                const void* source, uint64_t bytes)
{
    if (&ops != &tower_ops) {
        std::memcpy(destination, source, bytes);
        return;
    }
    TOWER_COUNT(source_rows, 1); TOWER_COUNT(source_bytes, bytes);
    tower_convert_involution(static_cast<const uint8_t*>(source), static_cast<uint8_t*>(destination), bytes, &map);
}
void Finish(void** work, unsigned recovery_prefix, uint64_t bytes)
{
    for (unsigned i = 0; i < recovery_prefix; ++i) {
        TOWER_COUNT(output_rows, 1); TOWER_COUNT(output_bytes, bytes);
        auto* p = static_cast<uint8_t*>(work[i]);
        tower_convert_involution(p, p, bytes, &map);
    }
}
} // namespace tower_encoder
