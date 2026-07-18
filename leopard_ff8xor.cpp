/* Experimental native plane-sliced FF8 XOR-circuit public API. */

#include "leopard_ff8xor.h"
#include "LeopardCommon.h"
#include "LeopardFF8Xor.h"

#include <string.h>

extern "C" int leo_internal_is_initialized();

namespace {

static bool EnsureExperimentalInitialized()
{
    // C++11 function-local static initialization provides a one-time,
    // thread-safe lazy boundary.  Keeping this reference out of leopard.cpp
    // avoids pulling generated circuit code into packed-only static clients.
    static const bool initialized = leopard::ff8xor::Initialize();
    return initialized;
}

static bool IsSupportedFF8Transform(
    unsigned original_count,
    unsigned recovery_count)
{
    if (original_count == 0 || original_count >= leopard::ff8xor::kOrder ||
        recovery_count == 0 || recovery_count > original_count)
        return false;

    const unsigned m = recovery_count == 1
        ? 1
        : leopard::NextPow2(recovery_count);
    return static_cast<uint64_t>(m) + original_count <=
        leopard::ff8xor::kOrder;
}

static void EncodeM1(
    uint64_t buffer_bytes,
    unsigned original_count,
    const void* const* original_data,
    void* recovery_data)
{
    memcpy(recovery_data, original_data[0], static_cast<size_t>(buffer_bytes));

    for (unsigned i = 1; i < original_count; ++i)
    {
        leopard::ff8xor::XorBuffer(
            buffer_bytes, recovery_data, original_data[i]);
    }
}

static void DecodeM1(
    uint64_t buffer_bytes,
    unsigned original_count,
    const void* const* original_data,
    const void* recovery_data,
    void* work_data)
{
    memcpy(work_data, recovery_data, static_cast<size_t>(buffer_bytes));

    for (unsigned i = 0; i < original_count; ++i)
    {
        if (original_data[i])
        {
            leopard::ff8xor::XorBuffer(
                buffer_bytes, work_data, original_data[i]);
        }
    }
}

} // namespace

extern "C" {

LEO_EXPORT unsigned leo_ff8xor_encode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    if (!IsSupportedFF8Transform(original_count, recovery_count))
        return 0;
    return leo_encode_work_count(original_count, recovery_count);
}

LEO_EXPORT LeopardResult leo_ff8xor_encode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned work_count,
    const void* const* original_data,
    void** work_data)
{
    // Fast paths below do not enter ReedSolomonEncode.  Reset at the public
    // call boundary so per-thread diagnostic statistics can never describe a
    // stale prior transform after k=1, r=1, or a rejected request.
    leopard::ff8xor::ResetMaterializationStatistics();

    if (buffer_bytes == 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;
    if (recovery_count == 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;
    if (!original_data || !work_data)
        return Leopard_InvalidInput;
    if (!leo_internal_is_initialized())
        return Leopard_CallInitialize;
    if (!EnsureExperimentalInitialized())
        return Leopard_Platform;
    if (!IsSupportedFF8Transform(original_count, recovery_count))
        return Leopard_TooMuchData;
    if (work_count != leo_ff8xor_encode_work_count(
            original_count, recovery_count))
        return Leopard_InvalidCounts;

    if (original_count == 1)
    {
        for (unsigned i = 0; i < recovery_count; ++i)
            memcpy(work_data[i], original_data[i], static_cast<size_t>(buffer_bytes));
        return Leopard_Success;
    }

    if (recovery_count == 1)
    {
        EncodeM1(buffer_bytes, original_count, original_data, work_data[0]);
        return Leopard_Success;
    }

    const unsigned m = leopard::NextPow2(recovery_count);
    leopard::ff8xor::ReedSolomonEncode(
        buffer_bytes,
        original_count,
        recovery_count,
        m,
        original_data,
        work_data);
    return Leopard_Success;
}

LEO_EXPORT unsigned leo_ff8xor_decode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    if (!IsSupportedFF8Transform(original_count, recovery_count))
        return 0;
    return leo_decode_work_count(original_count, recovery_count);
}

LEO_EXPORT LeopardResult leo_ff8xor_decode(
    uint64_t buffer_bytes,
    unsigned original_count,
    unsigned recovery_count,
    unsigned work_count,
    const void* const* original_data,
    const void* const* recovery_data,
    void** work_data)
{
    // This also covers no-loss, k=1, and r=1 exits that bypass the transform.
    leopard::ff8xor::ResetMaterializationStatistics();

    if (buffer_bytes == 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;
    if (recovery_count == 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;
    if (!original_data || !recovery_data || !work_data)
        return Leopard_InvalidInput;
    if (!leo_internal_is_initialized())
        return Leopard_CallInitialize;
    if (!EnsureExperimentalInitialized())
        return Leopard_Platform;
    // Reject unbounded caller counts before walking either pointer array.
    // Besides matching the FF8-only contract, this guarantees that malformed
    // out-of-range counts return TooMuchData instead of causing an array read
    // proportional to an unsupported count.
    if (!IsSupportedFF8Transform(original_count, recovery_count))
        return Leopard_TooMuchData;

    unsigned original_loss_count = 0;
    unsigned original_loss_i = 0;
    for (unsigned i = 0; i < original_count; ++i)
    {
        if (!original_data[i])
        {
            ++original_loss_count;
            original_loss_i = i;
        }
    }

    unsigned recovery_got_count = 0;
    unsigned recovery_got_i = 0;
    for (unsigned i = 0; i < recovery_count; ++i)
    {
        if (recovery_data[i])
        {
            ++recovery_got_count;
            recovery_got_i = i;
        }
    }
    if (recovery_got_count < original_loss_count)
        return Leopard_NeedMoreData;
    if (work_count != leo_ff8xor_decode_work_count(
            original_count, recovery_count))
        return Leopard_InvalidCounts;

    if (original_loss_count == 0)
    {
        for (unsigned i = 0; i < original_count; ++i)
            memcpy(work_data[i], original_data[i],
                static_cast<size_t>(buffer_bytes));
        return Leopard_Success;
    }

    if (original_count == 1)
    {
        memcpy(work_data[0], recovery_data[recovery_got_i],
            static_cast<size_t>(buffer_bytes));
        return Leopard_Success;
    }

    if (recovery_count == 1)
    {
        DecodeM1(
            buffer_bytes,
            original_count,
            original_data,
            recovery_data[0],
            work_data[original_loss_i]);
        return Leopard_Success;
    }

    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);
    leopard::ff8xor::ReedSolomonDecode(
        buffer_bytes,
        original_count,
        recovery_count,
        m,
        n,
        original_data,
        recovery_data,
        work_data);
    return Leopard_Success;
}

} // extern "C"
