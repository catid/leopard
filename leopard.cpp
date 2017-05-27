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
#include "LeopardFF8.h"
#include "LeopardFF16.h"

extern "C" {


//------------------------------------------------------------------------------
// Initialization API

static bool m_Initialized = false;

LEO_EXPORT int leo_init_(int version)
{
    if (version != LEO_VERSION)
        return Leopard_InvalidInput;

    if (!leopard::ff8::Initialize())
        return Leopard_Platform;

    if (!leopard::ff16::Initialize())
        return Leopard_Platform;

    m_Initialized = true;
    return Leopard_Success;
}


//------------------------------------------------------------------------------
// Encoder API

LEO_EXPORT unsigned leo_encode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    return leopard::NextPow2(recovery_count) * 2;
}

LEO_EXPORT LeopardResult leo_encode(
    uint64_t buffer_bytes,              // Number of bytes in each data buffer
    unsigned original_count,            // Number of original_data[] buffer pointers
    unsigned recovery_count,            // Number of recovery_data[] buffer pointers
    unsigned work_count,                // Number of work_data[] buffer pointers, from leo_encode_work_count()
    void* const * const original_data,  // Array of pointers to original data buffers
    void** work_data,                   // Array of work buffers
    unsigned flags)                     // Operation flags
{
    if (buffer_bytes <= 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;

    if (recovery_count <= 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;

    if (!original_data || !work_data)
        return Leopard_InvalidInput;

    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);

    if (work_count != m * 2)
        return Leopard_InvalidCounts;

    const bool mt = (flags & LeopardFlags_Multithreaded) != 0;

    if (n <= leopard::ff8::kOrder)
    {
        leopard::ff8::Encode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
    }
    else if (n <= leopard::ff16::kOrder)
    {
        leopard::ff16::Encode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            original_data,
            work_data);
    }
    else
        return Leopard_TooMuchData;

    return Leopard_Success;
}


//------------------------------------------------------------------------------
// Decoder API

LEO_EXPORT unsigned leo_decode_work_count(
    unsigned original_count,
    unsigned recovery_count)
{
    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);
    return n;
}

LEO_EXPORT LeopardResult leo_decode(
    uint64_t buffer_bytes,              // Number of bytes in each data buffer
    unsigned original_count,            // Number of original_data[] buffer pointers
    unsigned recovery_count,            // Number of recovery_data[] buffer pointers
    unsigned work_count,                // Number of buffer pointers in work_data[]
    void* const * const original_data,  // Array of original data buffers
    void* const * const recovery_data,  // Array of recovery data buffers
    void** work_data,                   // Array of work data buffers
    unsigned flags)                     // Operation flags
{
    if (buffer_bytes <= 0 || buffer_bytes % 64 != 0)
        return Leopard_InvalidSize;

    if (recovery_count <= 0 || recovery_count > original_count)
        return Leopard_InvalidCounts;

    if (!original_data || !recovery_data || !work_data)
        return Leopard_InvalidInput;

    const unsigned m = leopard::NextPow2(recovery_count);
    const unsigned n = leopard::NextPow2(m + original_count);

    if (work_count != n)
        return Leopard_InvalidCounts;

    const bool mt = (flags & LeopardFlags_Multithreaded) != 0;

    if (n <= leopard::ff8::kOrder)
    {
        leopard::ff8::Decode(
            buffer_bytes,
            original_count,
            recovery_count,
            m,
            n,
            original_data,
            recovery_data,
            work_data);
    }
    else if (n <= leopard::ff16::kOrder)
    {
        leopard::ff16::Decode(
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
        return Leopard_TooMuchData;

    return Leopard_Success;
}


} // extern "C"
