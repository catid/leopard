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
    Leopard-RS
    MDS Reed-Solomon Erasure Correction Codes for Large Data in C

    Algorithms are described in LeopardCommon.h


    Inspired by discussion with:

    Sian-Jhen Lin <sjhenglin@gmail.com> : Author of {1} {3}, basis for Leopard
    Bulat Ziganshin <bulat.ziganshin@gmail.com> : Author of FastECC
    Yutaka Sawada <tenfon@outlook.jp> : Author of MultiPar


    References:

    {1} S.-J. Lin, T. Y. Al-Naffouri, Y. S. Han, and W.-H. Chung,
    "Novel Polynomial Basis with Fast Fourier Transform
    and Its Application to Reed-Solomon Erasure Codes"
    IEEE Trans. on Information Theory, pp. 6284-6299, November, 2016.

    {2} D. G. Cantor, "On arithmetical algorithms over finite fields",
    Journal of Combinatorial Theory, Series A, vol. 50, no. 2, pp. 285-300, 1989.

    {3} Sian-Jheng Lin, Wei-Ho Chung, "An Efficient (n, k) Information
    Dispersal Algorithm for High Code Rate System over Fermat Fields,"
    IEEE Commun. Lett., vol.16, no.12, pp. 2036-2039, Dec. 2012.

    {4} Plank, J. S., Greenan, K. M., Miller, E. L., "Screaming fast Galois Field
    arithmetic using Intel SIMD instructions."  In: FAST-2013: 11th Usenix
    Conference on File and Storage Technologies, San Jose, 2013
*/

// Library version
#define LEO_VERSION 2

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

#include <stdint.h>

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

    Leopard_NeedMoreData      = -1, // Not enough recovery data received
    Leopard_TooMuchData       = -2, // Buffer counts are too high
    Leopard_InvalidSize       = -3, // Buffer size must be a multiple of 64 bytes
    Leopard_InvalidCounts     = -4, // Invalid counts provided
    Leopard_InvalidInput      = -5, // A function parameter was invalid
    Leopard_Platform          = -6, // Platform is unsupported
    Leopard_CallInitialize    = -7, // Call leo_init() first
} LeopardResult;

// Convert Leopard result to string
LEO_EXPORT const char* leo_result_string(LeopardResult result);


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

    The sum of original_count + recovery_count must not exceed 65536.
    The recovery_count <= original_count.

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
    * The first set of recovery_count buffers in work_data will be the result.
    Returns other values on errors.
*/
LEO_EXPORT LeopardResult leo_encode(
    uint64_t buffer_bytes,                    // Number of bytes in each data buffer
    unsigned original_count,                  // Number of original_data[] buffer pointers
    unsigned recovery_count,                  // Number of recovery_data[] buffer pointers
    unsigned work_count,                      // Number of work_data[] buffer pointers, from leo_encode_work_count()
    const void* const * const original_data,  // Array of pointers to original data buffers
    void** work_data);                        // Array of work buffers


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

    Lost original/recovery data should be set to NULL.

    The sum of recovery_count + the number of non-NULL original data must be at
    least original_count in order to perform recovery.

    Returns Leopard_Success on success.
    Returns other values on errors.
*/
LEO_EXPORT LeopardResult leo_decode(
    uint64_t buffer_bytes,                    // Number of bytes in each data buffer
    unsigned original_count,                  // Number of original_data[] buffer pointers
    unsigned recovery_count,                  // Number of recovery_data[] buffer pointers
    unsigned work_count,                      // Number of buffer pointers in work_data[]
    const void* const * const original_data,  // Array of original data buffers
    const void* const * const recovery_data,  // Array of recovery data buffers
    void** work_data);                        // Array of work data buffers


#ifdef __cplusplus
}
#endif


#endif // CAT_LEOPARD_RS_H
