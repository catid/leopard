// Test-only external public-call witness; not a performance executable.
#ifdef LEO_PAIRED_NATIVE
#include "leopard.h"
#else
#include "leopard2.h"
#include "Leopard2Direct.h"
#endif
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
struct Witness {
    unsigned calls = 0, states[3] = {}, apis[3] = {};
    uint64_t order = UINT64_C(14695981039346656037), bytes = 0;
    size_t scratch_bytes = 0;
    const void* codec = NULL;
    const void* scratch = NULL;
    const void* input_array = NULL;
    const void* output_array = NULL;
    std::vector<const void*> inputs, outputs;
    void Observe(unsigned api, unsigned state, const void* c, uint64_t b,
        unsigned k, unsigned r, const void* const* in, void* const* out,
        const void* work, size_t size)
    {
        if (calls == 0) {
            codec = c; bytes = b; scratch = work; scratch_bytes = size;
            input_array = in; output_array = out;
            inputs.assign(in, in + k); outputs.assign(out, out + r);
        }
        bool same = c == codec && b == bytes && work == scratch && size == scratch_bytes &&
            in == input_array && out == output_array && inputs.size() == k && outputs.size() == r;
        for (unsigned i = 0; same && i < k; ++i) same &= inputs[i] == in[i];
        for (unsigned i = 0; same && i < r; ++i) same &= outputs[i] == out[i];
        if (!same) { std::fputs("public call buffers changed\n", stderr); std::exit(87); }
        ++calls; ++states[state]; ++apis[api];
        order = (order ^ (state + 4 * api)) * UINT64_C(1099511628211);
    }
    ~Witness() {
        std::printf("{\"schema\":\"leopard-paired-witness/v1\",\"calls\":%u,"
            "\"states\":[%u,%u,%u],\"apis\":[%u,%u,%u],\"order_hash\":\"%016llx\"}\n",
            calls, states[0], states[1], states[2], apis[0], apis[1], apis[2],
            static_cast<unsigned long long>(order));
    }
} witness;
}

#ifdef LEO_PAIRED_NATIVE
extern "C" LeopardResult __real_leo_encode(uint64_t, unsigned, unsigned, unsigned, const void* const*, void**);
extern "C" LeopardResult __wrap_leo_encode(uint64_t b, unsigned k, unsigned r, unsigned n,
    const void* const* in, void** out)
{
    witness.Observe(2, 2, NULL, b, k, n, in, out, out[0], static_cast<size_t>(n) * b);
    return __real_leo_encode(b, k, r, n, in, out);
}
#else
extern "C" leo2_result __real_leo2_encode(const leo2_codec*, uint64_t, const void* const*, void* const*, void*, size_t);
extern "C" leo2_result __real_leo2_encode_batch(const leo2_codec*, const leo2_encode_batch_item*, size_t);
namespace {
void Observe(unsigned api, const leo2_codec* codec, uint64_t bytes,
    const void* const* in, void* const* out, void* scratch, size_t size)
{
    witness.Observe(api, leopard2_internal::AutoGF16GFNIR19932EnabledForDiagnostics() ? 1 : 0,
        codec, bytes, leo2_codec_original_count(codec), leo2_codec_recovery_count(codec),
        in, out, scratch, size);
}
}
extern "C" leo2_result __wrap_leo2_encode(const leo2_codec* codec, uint64_t bytes,
    const void* const* in, void* const* out, void* scratch, size_t size)
{
    Observe(0, codec, bytes, in, out, scratch, size);
    return __real_leo2_encode(codec, bytes, in, out, scratch, size);
}
extern "C" leo2_result __wrap_leo2_encode_batch(const leo2_codec* codec,
    const leo2_encode_batch_item* items, size_t count)
{
    if (count != 1) { std::fputs("not a one-item public batch\n", stderr); std::exit(87); }
    const leo2_encode_batch_item& item = items[0];
    Observe(1, codec, item.shard_bytes, item.original, item.recovery, item.scratch, item.scratch_bytes);
    return __real_leo2_encode_batch(codec, items, count);
}
#endif
