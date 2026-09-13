/* Minimal Leopard2 encode/decode example. */
#define _POSIX_C_SOURCE 200112L
#include "leopard2.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int ok(leo2_result result, const char *where) {
    if (result != LEO2_SUCCESS) {
        fprintf(stderr, "%s: %s\n", where, leo2_result_string(result));
        return 0;
    }
    return 1;
}

int main(void) {
    enum { K = 4, R = 2, BYTES = 64 };
    uint8_t original_storage[K][BYTES] = {{0}};
    uint8_t recovery_storage[R][BYTES] = {{0}};
    uint8_t restored_storage[BYTES] = {0};
    const void *original[K];
    void *recovery[R];
    const void *received_original[K];
    const void *received_recovery[R];
    void *restored[K] = {0};
    uint8_t original_present[K] = {1, 0, 1, 1};
    uint8_t recovery_present[R] = {1, 1};
    leo2_context *context = NULL;
    leo2_codec *codec = NULL;
    leo2_decode_plan *plan = NULL;
    leo2_context_options context_options = {
        sizeof(context_options), LEO2_BACKEND_AUTO, 1, 0
    };
    leo2_codec_options codec_options = {
        sizeof(codec_options), 0, 0, LEO2_SHARD_LAYOUT_NATIVE_V1
    };
    size_t scratch_bytes;
    void *scratch = NULL;

    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < BYTES; ++j)
            original_storage[i][j] = (uint8_t)(i * BYTES + j);
        original[i] = original_storage[i];
    }
    for (int i = 0; i < R; ++i) recovery[i] = recovery_storage[i];

    if (!ok(leo2_context_create(&context_options, &context), "context")) return 1;
    if (!ok(leo2_codec_create(context, K, R, LEO2_PROFILE_LOW_V1,
                              LEO2_FIELD_GF8, &codec_options, &codec), "codec")) return 1;
    scratch_bytes = 0;
    if (!ok(leo2_encode_scratch_size(codec, BYTES, &scratch_bytes), "encode scratch")) return 1;
    if (scratch_bytes && posix_memalign(&scratch, leo2_scratch_alignment(), scratch_bytes)) return 1;
    if (!ok(leo2_encode(codec, BYTES, original, recovery, scratch, scratch_bytes), "encode")) return 1;
    free(scratch);
    scratch = NULL;

    for (int i = 0; i < K; ++i) received_original[i] = original_present[i] ? original_storage[i] : NULL;
    for (int i = 0; i < R; ++i) received_recovery[i] = recovery_storage[i];
    restored[1] = restored_storage;
    if (!ok(leo2_decode_plan_create(codec, original_present, recovery_present, &plan), "decode plan")) return 1;
    scratch_bytes = 0;
    if (!ok(leo2_decode_plan_scratch_size(plan, BYTES, &scratch_bytes), "decode scratch")) return 1;
    if (scratch_bytes && posix_memalign(&scratch, leo2_scratch_alignment(), scratch_bytes)) return 1;
    if (!ok(leo2_decode_plan_execute(plan, BYTES, received_original, received_recovery,
                                     restored, scratch, scratch_bytes), "decode")) return 1;
    if (memcmp(restored_storage, original_storage[1], BYTES) != 0) return 1;

    free(scratch);
    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    puts("Leopard2 encode/decode succeeded");
    return 0;
}
