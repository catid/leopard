/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "leopard.h"
#include "LeopardCommon.h"
#include "Leopard2Backend.h"

#ifdef LEO_HAS_FF8
    #include "LeopardFF8.h"
#endif // LEO_HAS_FF8
#ifdef LEO_HAS_FF16
    #include "LeopardFF16.h"
#endif // LEO_HAS_FF16

#include <atomic>
#include <mutex>
#include <string.h>

//------------------------------------------------------------------------------
// Initialization API

static std::atomic<bool> m_Initialized(false);
static const unsigned kLegacyFF8Order = 256;
static const unsigned kLegacyFF16Order = 65536;

enum LegacyGeometryResult
{
    LegacyGeometryValid,
    LegacyGeometryInvalidCounts,
    LegacyGeometryTooMuchData
};

struct LegacyGeometry
{
    unsigned m;
    unsigned n;
    unsigned encode_work_count;
    unsigned decode_work_count;
};

static unsigned LegacyNextPow2(unsigned value)
{
    // All callers bound value to [2, 65536].  Keep this legacy helper local so
    // its historical geometry remains explicit and independently bounded.
    unsigned result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static LegacyGeometryResult GetLegacyGeometry(
    unsigned original_count,
    unsigned recovery_count,
    LegacyGeometry& geometry)
{
    geometry = LegacyGeometry();
    if (original_count == 0 || recovery_count == 0 ||
        recovery_count > original_count)
        return LegacyGeometryInvalidCounts;

    const uint64_t transmitted_count =
        static_cast<uint64_t>(original_count) + recovery_count;
    if (transmitted_count > kLegacyFF16Order)
        return LegacyGeometryTooMuchData;

    if (original_count == 1)
    {
        geometry.m = 1;
        geometry.n = 1;
        geometry.encode_work_count = recovery_count;
        geometry.decode_work_count = original_count;
        return LegacyGeometryValid;
    }
    if (recovery_count == 1)
    {
        geometry.m = 1;
        geometry.n = original_count;
        geometry.encode_work_count = 1;
        geometry.decode_work_count = original_count;
        return LegacyGeometryValid;
    }

    geometry.m = LegacyNextPow2(recovery_count);
    const uint64_t parent_input =
        static_cast<uint64_t>(geometry.m) + original_count;
    if (parent_input > kLegacyFF16Order)
        return LegacyGeometryTooMuchData;
    geometry.n = LegacyNextPow2(static_cast<unsigned>(parent_input));
    geometry.encode_work_count = geometry.m * 2;
    geometry.decode_work_count = geometry.n;
    return LegacyGeometryValid;
}

static bool LegacyPointerArraySpanRepresentable(
    const void* pointer_array,
    unsigned pointer_count)
{
    if (!pointer_array || pointer_count == 0)
        return false;

    // Counts are bounded by GetLegacyGeometry(), so this multiplication is
    // representable even on a 32-bit target.  This check deliberately proves
    // only that indexing the top-level pointer array cannot wrap uintptr_t; it
    // cannot prove that the pointer names an accessible C object.
    const uintptr_t address = reinterpret_cast<uintptr_t>(pointer_array);
    const uintptr_t span_minus_one =
        static_cast<uintptr_t>(pointer_count) * sizeof(void*) - 1;
    return address <= UINTPTR_MAX - span_minus_one;
}

namespace leopard {

int InitializeLibrary(
    backend::QualificationStatus* qualification_failure_out)
{
    if (qualification_failure_out)
        *qualification_failure_out = backend::QualificationAvailable;

    // Initialization publishes process-global CPU features, field tables, and
    // backend tables.  Serialize the complete transaction so legacy callers
    // and Leopard2 context creation may safely contend in a fresh process.
    // Do not cache failures: a later call may succeed after transient resource
    // pressure, while a successfully published state remains idempotent.
    static std::mutex initialize_mutex;
    std::lock_guard<std::mutex> lock(initialize_mutex);
    if (m_Initialized.load(std::memory_order_acquire))
        return Leopard_Success;

    leopard::InitializeCPUArch();

#ifdef LEO_HAS_FF8
    if (!leopard::ff8::Initialize())
        return Leopard_Platform;
#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
    if (!leopard::ff16::Initialize())
        return Leopard_Platform;
#endif // LEO_HAS_FF16

    leopard::backend::InitializeArgs backend_args = {
#ifdef LEO_HAS_FF8
        leopard::ff8::MultiplyLogElement,
#else
        NULL,
#endif
#ifdef LEO_HAS_FF16
        leopard::ff16::MultiplyLogElement,
#else
        NULL,
#endif
#ifdef LEO_HAS_FF8
        leopard::ff8::SkewLogTable()
#else
        NULL
#endif
    };
    if (!leopard::backend::Initialize(backend_args))
    {
        if (qualification_failure_out)
            *qualification_failure_out =
                leopard::backend::StartupQualificationFailure();
        return Leopard_Platform;
    }


    m_Initialized.store(true, std::memory_order_release);
    return Leopard_Success;
}

} // namespace leopard

extern "C" {

LEO_EXPORT int leo_init_(int version)
{
    if (version != LEO_VERSION)
        return Leopard_InvalidInput;
    return leopard::InitializeLibrary(NULL);
}

//------------------------------------------------------------------------------
// Result

LEO_EXPORT const char* leo_result_string(LeopardResult result)
{
    switch (result)
    {
    case Leopard_Success: return "Operation succeeded";
    case Leopard_NeedMoreData: return "Not enough recovery data received";
    case Leopard_TooMuchData: return "Buffer counts are too high";
    case Leopard_InvalidSize: return "Buffer size must be a multiple of 64 bytes";
    case Leopard_InvalidCounts: return "Invalid counts provided";
    case Leopard_InvalidInput: return "A function parameter was invalid";
    case Leopard_Platform: return "Platform is unsupported";
    case Leopard_CallInitialize: return "Call leo_init() first";
    }
    return "Unknown";
}


//------------------------------------------------------------------------------
// Encoder API

LEO_EXPORT unsigned leo_encode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    LegacyGeometry geometry;
    return GetLegacyGeometry(original_count, recovery_count, geometry) ==
            LegacyGeometryValid
        ? geometry.encode_work_count
        : 0;
}

// recovery_data = parity of original_data (xor sum)
static void EncodeM1(
    uint64_t buffer_bytes,
    unsigned original_count,
    const void* const * const original_data,
    void* recovery_data)
{
    memcpy(recovery_data, original_data[0], buffer_bytes);

    leopard::XORSummer summer;
    summer.Initialize(recovery_data);

    for (unsigned i = 1; i < original_count; ++i)
        summer.Add(original_data[i], buffer_bytes);

    summer.Finalize(buffer_bytes);
}

LEO_EXPORT LeopardResult leo_encode(
    uint64_t buffer_bytes,                    // Number of bytes in each data buffer
    unsigned original_count,                  // Number of original_data[] buffer pointers
    unsigned recovery_count,                  // Number of recovery_data[] buffer pointers
    unsigned work_count,                      // Number of work_data[] buffer pointers, from leo_encode_work_count()
    const void* const * const original_data,  // Array of pointers to original data buffers
    void** work_data)                         // Array of work buffers
{
    if (!leopard::LegacyBufferBytesValidForAddressLimit(
            buffer_bytes, static_cast<uint64_t>(SIZE_MAX)))
        return Leopard_InvalidSize;

    LegacyGeometry geometry;
    const LegacyGeometryResult geometry_result = GetLegacyGeometry(
        original_count, recovery_count, geometry);
    if (geometry_result == LegacyGeometryInvalidCounts)
        return Leopard_InvalidCounts;
    if (geometry_result == LegacyGeometryTooMuchData)
        return Leopard_TooMuchData;

    if (!original_data || !work_data)
        return Leopard_InvalidInput;

    if (!m_Initialized.load(std::memory_order_acquire))
        return Leopard_CallInitialize;

    if (work_count != geometry.encode_work_count)
        return Leopard_InvalidCounts;

    if (!LegacyPointerArraySpanRepresentable(
            original_data, original_count) ||
        !LegacyPointerArraySpanRepresentable(work_data, work_count))
        return Leopard_InvalidInput;

    // Handle k = 1 case
    if (original_count == 1)
    {
        for (unsigned i = 0; i < recovery_count; ++i)
            memcpy(work_data[i], original_data[i], buffer_bytes);
        return Leopard_Success;
    }

    // Handle m = 1 case
    if (recovery_count == 1)
    {
        EncodeM1(
            buffer_bytes,
            original_count,
            original_data,
            work_data[0]);
        return Leopard_Success;
    }

    const unsigned m = geometry.m;
    const unsigned n = geometry.n;

    if (n <= kLegacyFF8Order)
    {
#ifdef LEO_HAS_FF8
        leopard::ff8::ReedSolomonEncode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
        return Leopard_Success;
#else
        return Leopard_Platform;
#endif // LEO_HAS_FF8
    }
    if (n <= kLegacyFF16Order)
    {
#ifdef LEO_HAS_FF16
        leopard::ff16::ReedSolomonEncode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
        return Leopard_Success;
#else
        return Leopard_Platform;
#endif // LEO_HAS_FF16
    }
    return Leopard_TooMuchData;
}


//------------------------------------------------------------------------------
// Decoder API

LEO_EXPORT unsigned leo_decode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    LegacyGeometry geometry;
    return GetLegacyGeometry(original_count, recovery_count, geometry) ==
            LegacyGeometryValid
        ? geometry.decode_work_count
        : 0;
}

static void DecodeM1(
    uint64_t buffer_bytes,
    unsigned original_count,
    const void* const * original_data,
    const void* recovery_data,
    void* work_data)
{
    memcpy(work_data, recovery_data, buffer_bytes);

    leopard::XORSummer summer;
    summer.Initialize(work_data);

    for (unsigned i = 0; i < original_count; ++i)
        if (original_data[i])
            summer.Add(original_data[i], buffer_bytes);

    summer.Finalize(buffer_bytes);
}

LEO_EXPORT LeopardResult leo_decode(
    uint64_t buffer_bytes,                    // Number of bytes in each data buffer
    unsigned original_count,                  // Number of original_data[] buffer pointers
    unsigned recovery_count,                  // Number of recovery_data[] buffer pointers
    unsigned work_count,                      // Number of buffer pointers in work_data[]
    const void* const * const original_data,  // Array of original data buffers
    const void* const * const recovery_data,  // Array of recovery data buffers
    void** work_data)                         // Array of work data buffers
{
    if (!leopard::LegacyBufferBytesValidForAddressLimit(
            buffer_bytes, static_cast<uint64_t>(SIZE_MAX)))
        return Leopard_InvalidSize;

    LegacyGeometry geometry;
    const LegacyGeometryResult geometry_result = GetLegacyGeometry(
        original_count, recovery_count, geometry);
    if (geometry_result == LegacyGeometryInvalidCounts)
        return Leopard_InvalidCounts;
    if (geometry_result == LegacyGeometryTooMuchData)
        return Leopard_TooMuchData;

    if (!original_data || !recovery_data || !work_data)
        return Leopard_InvalidInput;

    if (!m_Initialized.load(std::memory_order_acquire))
        return Leopard_CallInitialize;

    // Validate arrays that are scanned below before reading their elements.
    // The work array is checked after its count, preserving the legacy error
    // precedence for an incorrect work_count.
    if (!LegacyPointerArraySpanRepresentable(
            original_data, original_count) ||
        !LegacyPointerArraySpanRepresentable(
            recovery_data, recovery_count))
        return Leopard_InvalidInput;

    // Check if not enough recovery data arrived
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

    if (work_count != geometry.decode_work_count)
        return Leopard_InvalidCounts;

    if (!LegacyPointerArraySpanRepresentable(work_data, work_count))
        return Leopard_InvalidInput;

    // Handle case original_loss_count = 0
    if (original_loss_count == 0)
    {
        for(unsigned i = 0; i < original_count; i++)
            memcpy(work_data[i], original_data[i], buffer_bytes);
        return Leopard_Success;
    }

    // Handle k = 1 case.  This must follow the no-loss path: no recovery
    // shard is required when the sole original is already available.
    if (original_count == 1)
    {
        memcpy(work_data[0], recovery_data[recovery_got_i], buffer_bytes);
        return Leopard_Success;
    }

    // Handle m = 1 case
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

    const unsigned m = geometry.m;
    const unsigned n = geometry.n;

    if (n <= kLegacyFF8Order)
    {
#ifdef LEO_HAS_FF8
        leopard::ff8::ReedSolomonDecode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            n,
            original_data,
            recovery_data,
            work_data);
        return Leopard_Success;
#else
        return Leopard_Platform;
#endif // LEO_HAS_FF8
    }
    if (n <= kLegacyFF16Order)
    {
#ifdef LEO_HAS_FF16
        leopard::ff16::ReedSolomonDecode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            n,
            original_data,
            recovery_data,
            work_data);
        return Leopard_Success;
#else
        return Leopard_Platform;
#endif // LEO_HAS_FF16
    }
    return Leopard_TooMuchData;
}


} // extern "C"
