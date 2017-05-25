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

#ifndef CAT_LEOPARD_RS_H
#define CAT_LEOPARD_RS_H

/*
    Leopard-RS: Reed-Solomon Error Correction Coding for Extremely Large Data

    S.-J. Lin, T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.
    http://ct.ee.ntust.edu.tw/it2016-2.pdf
*/

// Library version
#define LEO_VERSION 1

// Tweak if the functions are exported or statically linked
//#define LEO_DLL /* Defined when building/linking as DLL */
//#define LEO_BUILDING /* Defined by the library makefile */

#if defined(LEO_BUILDING)
# if defined(LEO_DLL)
    #define LEO_EXPORT __declspec(dllexport)
# else
    #define LEO_EXPORT
# endif
#else
# if defined(LEO_DLL)
    #define LEO_EXPORT __declspec(dllimport)
# else
    #define LEO_EXPORT extern
# endif
#endif


#ifdef __cplusplus
extern "C" {
#endif


//------------------------------------------------------------------------------
// Initialization API

/*
    leo_init()

    Perform static initialization for the library, verifying that the platform
    is supported.

    Returns 0 on success and other values on failure.
*/

LEO_EXPORT int leo_init_(int version);
#define leo_init() leo_init_(LEO_VERSION)


//------------------------------------------------------------------------------
// Shared Constants / Datatypes

// Results
typedef enum LeopardResultT
{
    Leopard_Success           =  0, // Operation succeeded

    Leopard_TooMuchData       = -1, // Buffer counts are too high
    Leopard_InvalidBlockSize  = -2, // Buffer size must be a multiple of 64 bytes
    Leopard_InvalidInput      = -3, // A function parameter was invalid
    Leopard_Platform          = -4, // Platform is unsupported
    Leopard_OutOfMemory       = -5, // Out of memory error occurred
    Leopard_Unexpected        = -6, // Unexpected error - Software bug?
} LeopardResult;

// Results
typedef enum LeopardFlagsT
{
    LeopardFlags_Defaults      = 0, // Default settings

    LeopardFlags_Multithreaded = 1, // Enable multiple threads
} LeopardFlags;


//------------------------------------------------------------------------------
// Encoder API

/*
    leo_encode_work_count()

	Calculate the number of work_data buffers to provide to leo_encode().

    The sum of original_count + recovery_count must not exceed 65536.

	Returns the work_count value to pass into leo_encode().
    Returns 0 on invalid input.
*/

LEO_EXPORT unsigned leo_encode_work_count(
    unsigned original_count,
    unsigned recovery_count);

/*
    leo_encode()

    Generate recovery data.

    original_count: Number of original_data[] buffers provided.
    recovery_count: Number of desired recovery data buffers.
    buffer_bytes:   Number of bytes in each data buffer.
    original_data:  Array of pointers to original data buffers.
    work_count:     Number of work_data[] buffers, from leo_encode_work_count().
    work_data:      Array of pointers to work data buffers.
    flags:          Flags for encoding e.g. LeopardFlag_Multithreaded

    The sum of original_count + recovery_count must not exceed 65536.
    The buffer_bytes must be a multiple of 64.
    Each buffer should have the same number of bytes.
    Even the last piece must be rounded up to the block size.

    Let buffer_bytes = The number of bytes in each buffer:

        original_count = static_cast<unsigned>(
            ((uint64_t)total_bytes + buffer_bytes - 1) / buffer_bytes);

    Or if the number of pieces is known:

        buffer_bytes = static_cast<unsigned>(
            ((uint64_t)total_bytes + original_count - 1) / original_count);

    Returns Leopard_Success on success.
    The first set of recovery_count buffers in work_data will be the result.

    Returns Leopard_TooMuchData if the data is too large.
    Returns Leopard_InvalidBlockSize if the data is the wrong size.
    Returns Leopard_InvalidInput on invalid input.
    Returns other values on errors.
*/
LEO_EXPORT LeopardResult leo_encode(
    unsigned buffer_bytes,              // Number of bytes in each data buffer
    unsigned original_count,            // Number of original_data[] buffer pointers
    unsigned recovery_count,            // Number of recovery_data[] buffer pointers
    unsigned work_count,                // Number of work_data[] buffer pointers, from leo_encode_work_count()
    void* const * const original_data,  // Array of pointers to original data buffers
    void** work_data,                   // Array of work buffers
    unsigned flags);                    // Operation flags


//------------------------------------------------------------------------------
// Decoder API

/*
    leo_decode_work_count()

	Calculate the number of work_data buffers to provide to leo_decode().

    The sum of original_count + recovery_count must not exceed 65536.

	Returns the work_count value to pass into leo_encode().
    Returns 0 on invalid input.
*/

LEO_EXPORT unsigned leo_decode_work_count(
    unsigned original_count,
    unsigned recovery_count);

/*
    leo_decode()

    Decode original data from recovery data.

    buffer_bytes:   Number of bytes in each data buffer.
    original_count: Number of original_data[] buffers provided.
    original_data:  Array of pointers to original data buffers.
    recovery_count: Number of recovery_data[] buffers provided.
    recovery_data:  Array of pointers to recovery data buffers.
    work_count:     Number of work_data[] buffers, from leo_decode_work_count().
    work_data:      Array of pointers to recovery data buffers.
    flags:          Flags for encoding e.g. LeopardFlag_Multithreaded

    Lost original/recovery data should be set to NULL.

    The sum of recovery_count + the number of non-NULL original data must be at
    least original_count in order to perform recovery.

    Returns Leopard_Success on success.
    Returns other values on errors.
*/
LEO_EXPORT LeopardResult leo_decode(
    unsigned buffer_bytes,              // Number of bytes in each data buffer
    unsigned original_count,            // Number of original_data[] buffer pointers
    unsigned recovery_count,            // Number of recovery_data[] buffer pointers
    unsigned work_count,                // Number of buffer pointers in work_data[]
    void* const * const original_data,  // Array of original data buffers
    void* const * const recovery_data,  // Array of recovery data buffers
    void** work_data,                   // Array of work data buffers
    unsigned flags);                    // Operation flags


#ifdef __cplusplus
}
#endif


#endif // CAT_LEOPARD_RS_H
