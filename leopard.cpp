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

#ifdef LEO_HAS_FF8
    #include "LeopardFF8.h"
#endif // LEO_HAS_FF8
#ifdef LEO_HAS_FF16
    #include "LeopardFF16.h"
#endif // LEO_HAS_FF16

#include <string.h>

extern "C" {


//------------------------------------------------------------------------------
// Initialization API

static bool m_Initialized = false;

LEO_EXPORT int leo_init_(int version)
{
    if (version != LEO_VERSION)
        return Leopard_InvalidInput;

    leopard::InitializeCPUArch();

#ifdef LEO_HAS_FF8
    if (!leopard::ff8::Initialize())
        return Leopard_Platform;
#endif // LEO_HAS_FF8

#ifdef LEO_HAS_FF16
    if (!leopard::ff16::Initialize())
        return Leopard_Platform;
#endif // LEO_HAS_FF16


    m_Initialized = true;
    return Leopard_Success;
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
    if (original_count == 1)
        return recovery_count;
    if (recovery_count == 1)
        return 1;
    return leopard::NextPow2(recovery_count) * 2;
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
    if (buffer_bytes <= 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;

    if (recovery_count <= 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;

    if (!original_data || !work_data)
        return Leopard_InvalidInput;

    if (!m_Initialized)
        return Leopard_CallInitialize;

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

    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);

    if (work_count != m * 2)
        return Leopard_InvalidCounts;

#ifdef LEO_HAS_FF8
    if (n <= leopard::ff8::kOrder)
    {
        leopard::ff8::ReedSolomonEncode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
    }
    else
#endif // LEO_HAS_FF8
#ifdef LEO_HAS_FF16
    if (n <= leopard::ff16::kOrder)
    {
        leopard::ff16::ReedSolomonEncode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
    }
    else
#endif // LEO_HAS_FF16
        return Leopard_TooMuchData;

    return Leopard_Success;
}


//------------------------------------------------------------------------------
// Decoder API

LEO_EXPORT unsigned leo_decode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    if (original_count == 1 || recovery_count == 1)
        return original_count;
    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);
    return n;
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
    if (buffer_bytes <= 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;

    if (recovery_count <= 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;

    if (!original_data || !recovery_data || !work_data)
        return Leopard_InvalidInput;

    if (!m_Initialized)
        return Leopard_CallInitialize;

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

    // Handle k = 1 case
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

    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);

    if (work_count != n)
        return Leopard_InvalidCounts;

#ifdef LEO_HAS_FF8
    if (n <= leopard::ff8::kOrder)
    {
        leopard::ff8::ReedSolomonDecode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            n,
            original_data,
            recovery_data,
            work_data);
    }
    else
#endif // LEO_HAS_FF8
#ifdef LEO_HAS_FF16
    if (n <= leopard::ff16::kOrder)
    {
        leopard::ff16::ReedSolomonDecode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            n,
            original_data,
            recovery_data,
            work_data);
    }
    else
#endif // LEO_HAS_FF16
        return Leopard_TooMuchData;

    return Leopard_Success;
}


} // extern "C"
