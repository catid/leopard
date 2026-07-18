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

#ifndef CAT_LEOPARD2_DIRECT_TEST_H
#define CAT_LEOPARD2_DIRECT_TEST_H

#include "leopard2.h"

/*
    These hooks are intentionally absent from production builds and the public
    Leopard2 ABI.  Tests set the per-codec mode before concurrent execution.
*/
#ifdef LEO2_ENABLE_TEST_HOOKS

#ifdef __cplusplus
extern "C" {
#endif

typedef enum leo2_test_encode_mode {
    LEO2_TEST_ENCODE_AUTO = 0,
    LEO2_TEST_ENCODE_FORCE_DIRECT = 1,
    LEO2_TEST_ENCODE_FORCE_TRANSFORM = 2
} leo2_test_encode_mode;

LEO2_EXPORT leo2_result leo2_test_codec_set_encode_mode(
    leo2_codec* codec,
    leo2_test_encode_mode mode);

LEO2_EXPORT int leo2_test_codec_direct_encode_capable(
    const leo2_codec* codec);

LEO2_EXPORT leo2_result leo2_test_codec_encode_path(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    /* Number of non-null requested parity outputs, not a parity prefix. */
    uint32_t requested_recovery_count,
    int* direct_out);

LEO2_EXPORT void leo2_test_reset_generic_reveal_counts(void);

LEO2_EXPORT uint64_t leo2_test_generic_direct_reveal_shards(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LEO2_ENABLE_TEST_HOOKS */

#endif /* CAT_LEOPARD2_DIRECT_TEST_H */
