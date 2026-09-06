// Clock-free, driver-only callback attribution for leopard-79h.38.5.4.12.
// Link to the pinned portable archive; do not use with legacy in-field SIMD.
#define main CallbackWorkloadMain
#include "current_route_screen.cpp"
#undef main
#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include <omp.h>

namespace callback_probe {
using namespace leopard::backend;
enum Operation { Multiply, MultiplyAdd, Xor, Xor2To1, Xor4, Copy, Ifft2, Fft2,
    Fft2Out, Ifft2Xor, Ifft4, Fft4, Ifft4Out, Fft4Out, Ifft4Range, Fft4Range, Count };
const char* const names[] = {"multiply", "multiply_add", "xor", "xor_2to1", "xor4",
    "copy", "ifft2", "fft2", "fft2_out", "ifft2_xor", "ifft4", "fft4",
    "ifft4_out", "fft4_out", "ifft4_range", "fft4_range"};
struct Bucket { unsigned operation, distance, zero_mask; bool hint; uint64_t bytes, calls; };
struct Pass { unsigned kind, k, r, requested, side, sparse_blocks; uint64_t bytes, policy; };
struct State { unsigned bucket_count, pass_count; uint64_t calls; Bucket buckets[256]; Pass passes[16]; };
static State state = {};
static const Ops* active = NULL;

void Reset()
{
    Require(active == NULL, "reset inside callback scope");
    state = State{};
}

struct Scope
{
    explicit Scope(const Ops& original)
    {
        Require(active == NULL && omp_get_max_threads() == 1, "reentry or threaded probe");
        active = &original;
    }
    ~Scope() { active = NULL; }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
};

void Record(unsigned operation, uint64_t bytes, unsigned distance = 1,
            unsigned zero_mask = 0, bool hint = false)
{
    Require(active && operation < Count && distance > 0 && distance <= 65536 &&
            bytes <= 1048576 && zero_mask < 8 && state.calls < 1000000,
            "callback record bounds");
    unsigned index = 0;
    for (; index < state.bucket_count; ++index)
    {
        const Bucket& b = state.buckets[index];
        if (b.operation == operation && b.distance == distance && b.zero_mask == zero_mask &&
            b.hint == hint && b.bytes == bytes) break;
    }
    Require(index < 256, "callback bucket limit");
    if (index == state.bucket_count)
    {
        state.buckets[index] = Bucket{operation, distance, zero_mask, hint, bytes, 0};
        ++state.bucket_count;
    }
    ++state.buckets[index].calls;
    ++state.calls;
}

unsigned Mask(uint16_t a, uint16_t b, uint16_t c)
{
    return (a == 65535 ? 1U : 0U) | (b == 65535 ? 2U : 0U) | (c == 65535 ? 4U : 0U);
}
template<unsigned Id, FixedMultiply Ops::* Field>
void Fixed(void* d, const void* s, uint16_t log, uint64_t bytes)
{ Record(Id, bytes, 1, log == 65535); (active->*Field)(d, s, log, bytes); }
template<unsigned Id, XorMemory Ops::* Field>
void Memory(void* d, const void* s, uint64_t bytes)
{ Record(Id, bytes); (active->*Field)(d, s, bytes); }
void Memory2(void* d, const void* s0, const void* s1, uint64_t bytes)
{ Record(Xor2To1, bytes); active->xor_memory_2to1(d, s0, s1, bytes); }
void Memory4(void* d0, const void* s0, void* d1, const void* s1,
             void* d2, const void* s2, void* d3, const void* s3, uint64_t bytes)
{ Record(Xor4, bytes); active->xor_memory4(d0, s0, d1, s1, d2, s2, d3, s3, bytes); }
template<unsigned Id, Butterfly2 Ops::* Field>
void Pair(void* x, void* y, uint16_t log, uint64_t bytes)
{ Record(Id, bytes, 1, log == 65535); (active->*Field)(x, y, log, bytes); }
template<unsigned Id, FFTButterfly2Out Ops::* Field>
void PairOut(const void* x, const void* y, void* u, void* v, uint16_t log, uint64_t bytes)
{ Record(Id, bytes, 1, log == 65535); (active->*Field)(x, y, u, v, log, bytes); }
template<unsigned Id, Butterfly4 Ops::* Field>
void Quad(void* a, void* b, void* c, void* d, uint16_t l0, uint16_t l1, uint16_t l2, uint64_t bytes)
{ Record(Id, bytes, 1, Mask(l0,l1,l2)); (active->*Field)(a,b,c,d,l0,l1,l2,bytes); }
template<unsigned Id, FFTButterfly4Out Ops::* Field>
void QuadOut(const void* a, const void* b, const void* c, const void* d,
             void* e, void* f, void* g, void* h, uint16_t l0, uint16_t l1, uint16_t l2, uint64_t bytes)
{ Record(Id, bytes, 1, Mask(l0,l1,l2)); (active->*Field)(a,b,c,d,e,f,g,h,l0,l1,l2,bytes); }
template<unsigned Id, Butterfly4Range Ops::* Field>
void Range(void* const* work, unsigned distance, uint16_t l0, uint16_t l1, uint16_t l2,
           uint64_t bytes, bool hint)
{ Record(Id, bytes, distance, Mask(l0,l1,l2), hint); (active->*Field)(work,distance,l0,l1,l2,bytes,hint); }

// Only these sixteen entries change in a private copy; null stays null.
// Restoring the entries also checks that no other Ops byte changed.
#define CALLBACK_FIELDS(X) \
    X(ff16_multiply, (Fixed<Multiply, &Ops::ff16_multiply>)) \
    X(ff16_multiply_add, (Fixed<MultiplyAdd, &Ops::ff16_multiply_add>)) \
    X(xor_memory, (Memory<Xor, &Ops::xor_memory>)) \
    X(xor_memory_2to1, Memory2) \
    X(xor_memory4, Memory4) \
    X(copy_memory, (Memory<Copy, &Ops::copy_memory>)) \
    X(ff16_ifft_butterfly2, (Pair<Ifft2, &Ops::ff16_ifft_butterfly2>)) \
    X(ff16_fft_butterfly2, (Pair<Fft2, &Ops::ff16_fft_butterfly2>)) \
    X(ff16_fft_butterfly2_out, (PairOut<Fft2Out, &Ops::ff16_fft_butterfly2_out>)) \
    X(ff16_ifft_butterfly2_xor, (PairOut<Ifft2Xor, &Ops::ff16_ifft_butterfly2_xor>)) \
    X(ff16_ifft_butterfly4, (Quad<Ifft4, &Ops::ff16_ifft_butterfly4>)) \
    X(ff16_fft_butterfly4, (Quad<Fft4, &Ops::ff16_fft_butterfly4>)) \
    X(ff16_ifft_butterfly4_out, (QuadOut<Ifft4Out, &Ops::ff16_ifft_butterfly4_out>)) \
    X(ff16_fft_butterfly4_out, (QuadOut<Fft4Out, &Ops::ff16_fft_butterfly4_out>)) \
    X(ff16_ifft_butterfly4_range, (Range<Ifft4Range, &Ops::ff16_ifft_butterfly4_range>)) \
    X(ff16_fft_butterfly4_range, (Range<Fft4Range, &Ops::ff16_fft_butterfly4_range>))

Ops View(const Ops& original)
{
    Ops view;
    std::memcpy(&view, &original, sizeof(view));
#define INSTALL(field, function) if (original.field) view.field = function;
    CALLBACK_FIELDS(INSTALL)
#undef INSTALL
    return view;
}

void CheckView(const Ops& view, const Ops& original)
{
    Ops restored;
    std::memcpy(&restored, &view, sizeof(restored));
#define RESTORE(field, function) restored.field = original.field;
    CALLBACK_FIELDS(RESTORE)
#undef RESTORE
    Require(std::memcmp(&restored, &original, sizeof(restored)) == 0, "unrelated Ops byte changed");
}
#undef CALLBACK_FIELDS

void Print()
{
    std::fprintf(stderr, "{\"schema\":\"gf16-callback-counts/v1\",\"timed\":false,\"calls\":%llu,\"passes\":[",
        static_cast<unsigned long long>(state.calls));
    for (unsigned i = 0; i < state.pass_count; ++i)
    {
        const Pass& p = state.passes[i];
        std::fprintf(stderr, "%s{\"kind\":%u,\"k\":%u,\"r\":%u,\"requested\":%u,\"side\":%u,"
            "\"sparse_blocks\":%u,\"bytes\":%llu,\"source_policy\":%llu}", i ? "," : "", p.kind,p.k,p.r,
            p.requested,p.side,p.sparse_blocks,static_cast<unsigned long long>(p.bytes),
            static_cast<unsigned long long>(p.policy));
    }
    std::fprintf(stderr, "],\"buckets\":[");
    for (unsigned i = 0; i < state.bucket_count; ++i)
    {
        const Bucket& b = state.buckets[i];
        std::fprintf(stderr, "%s{\"op\":\"%s\",\"distance\":%u,\"zero_mask\":%u,\"prefer_fused\":%s,"
            "\"bytes\":%llu,\"calls\":%llu}", i ? "," : "", names[b.operation], b.distance,b.zero_mask,
            b.hint ? "true" : "false", static_cast<unsigned long long>(b.bytes),
            static_cast<unsigned long long>(b.calls));
    }
    std::fprintf(stderr, "]}\n");
}
}

#define CALLBACK_SYMBOL "_ZN7leopard4ff1633ReedSolomonEncodeWithSourcePolicyERKNS_7backend3OpsEmmjjjjPKPKvPPvPKN17leopard2_internal26SparseForwardPlanBatchViewE"
#define CALLBACK_ARGS const leopard::backend::Ops& ops, uint64_t bytes, uint64_t policy, \
    unsigned k, unsigned r, unsigned requested, unsigned side, const void* const* data, \
    void** work, const leopard2_internal::SparseForwardPlanBatchView* sparse
extern "C" void RealCallbackEncode(CALLBACK_ARGS) asm("__real_" CALLBACK_SYMBOL);
extern "C" void WrappedCallbackEncode(CALLBACK_ARGS) asm("__wrap_" CALLBACK_SYMBOL);
extern "C" void WrappedCallbackEncode(CALLBACK_ARGS)
{
    using namespace callback_probe;
    Require(state.pass_count < 16, "source-policy pass limit");
    Scope scope(ops);
    Ops original;
    std::memcpy(&original, &ops, sizeof(original));
    Ops view = View(ops);
    CheckView(view, ops);
    state.passes[state.pass_count++] = Pass{static_cast<unsigned>(ops.kind),k,r,requested,side,
        sparse ? sparse->block_count : 0U,bytes,policy};
    RealCallbackEncode(view, bytes, policy, k,r,requested,side,data,work,sparse);
    CheckView(view, ops);
    Require(std::memcmp(&original, &ops, sizeof(ops)) == 0, "published Ops changed");
}
#undef CALLBACK_ARGS
#undef CALLBACK_SYMBOL

#ifndef LEO_GF16_CALLBACK_NO_MAIN
int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) && std::strcmp(argv[1], "--check") == 0,
                "callback probe is check-only: --check cell[0..5] [parity_file]");
        omp_set_dynamic(0); omp_set_num_threads(1);
        callback_probe::Reset();
        const int result = CallbackWorkloadMain(argc, argv);
        if (result != 0) return result;
        Require(callback_probe::state.pass_count && callback_probe::state.calls, "empty callback probe");
        callback_probe::Print();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "callback probe: %s\n", error.what());
        return 1;
    }
}
#endif
