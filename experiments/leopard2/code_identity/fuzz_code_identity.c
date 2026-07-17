/*
 * Coverage-guided boundary for the experimental Leopard2 code identity.
 *
 * This target is built only when LEO2_BUILD_FUZZERS is enabled.  It does not
 * add the experimental identity implementation to libleopard or its public
 * ABI.
 */
#include "code_identity.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static void check(int condition)
{
    if (!condition) {
        __builtin_trap();
    }
}

static uint16_t fuzz_load16(const uint8_t *input)
{
    return (uint16_t)(((uint16_t)input[0] << 8) | (uint16_t)input[1]);
}

static uint32_t fuzz_load32(const uint8_t *input)
{
    return ((uint32_t)input[0] << 24) |
        ((uint32_t)input[1] << 16) |
        ((uint32_t)input[2] << 8) |
        (uint32_t)input[3];
}

static int identities_equal(
    const leo2_code_identity *left,
    const leo2_code_identity *right)
{
    size_t index;
    if (left->profile != right->profile ||
        left->profile_version != right->profile_version ||
        left->field != right->field ||
        left->field_version != right->field_version ||
        left->original_count != right->original_count ||
        left->recovery_count != right->recovery_count ||
        left->parent_count != right->parent_count ||
        left->padded_side != right->padded_side ||
        left->coordinate_map_version != right->coordinate_map_version ||
        left->metadata_count != right->metadata_count) {
        return 0;
    }
    for (index = 0; index < left->metadata_count; ++index) {
        const leo2_code_id_metadata *a = &left->metadata[index];
        const leo2_code_id_metadata *b = &right->metadata[index];
        if (a->type != b->type || a->value_bytes != b->value_bytes ||
            (a->value_bytes != 0u &&
             memcmp(a->value, b->value, a->value_bytes) != 0)) {
            return 0;
        }
    }
    return 1;
}

static void roundtrip_identity(
    const leo2_code_identity *identity,
    const uint8_t *canonical,
    size_t canonical_bytes)
{
    leo2_code_identity decoded;
    leo2_code_id_status status;
    uint8_t *output;
    size_t required = 0;

    check(leo2_code_identity_validate(identity) == LEO2_CODE_ID_OK);
    status = leo2_code_identity_serialize(identity, NULL, 0, &required);
    check(status == LEO2_CODE_ID_OK);
    check(required >= LEO2_CODE_ID_HEADER_BYTES);
    check(required <= LEO2_CODE_ID_MAX_BYTES);
    if (canonical != NULL) {
        check(required == canonical_bytes);
    }

    output = (uint8_t *)malloc(required);
    check(output != NULL);
    memset(output, 0xa5, required);
    if (required > 0u) {
        size_t undersized_required = 0;
        status = leo2_code_identity_serialize(identity, output, required - 1u,
            &undersized_required);
        check(status == LEO2_CODE_ID_BUFFER_TOO_SMALL);
        check(undersized_required == required);
        check(output[0] == UINT8_C(0xa5));
        check(output[required - 1u] == UINT8_C(0xa5));
    }
    status = leo2_code_identity_serialize(identity, output, required,
        &required);
    check(status == LEO2_CODE_ID_OK);
    if (canonical != NULL) {
        check(memcmp(output, canonical, canonical_bytes) == 0);
    }
    status = leo2_code_identity_deserialize(output, required, &decoded);
    check(status == LEO2_CODE_ID_OK);
    check(leo2_code_identity_validate(&decoded) == LEO2_CODE_ID_OK);
    check(identities_equal(identity, &decoded));
    free(output);
}

static void exercise_builder(const uint8_t *data, size_t size)
{
    leo2_code_identity identity;
    leo2_code_id_status status;
    size_t offset = 10u;
    unsigned item_count;
    unsigned item;

    if (size < offset) {
        return;
    }
    status = leo2_code_identity_make(&identity,
        (uint8_t)(data[0] % 5u),
        (uint8_t)(data[1] % 4u),
        fuzz_load32(data + 2),
        fuzz_load32(data + 6));
    if (status != LEO2_CODE_ID_OK) {
        check(leo2_code_id_status_string(status) != NULL);
        return;
    }

    item_count = size > offset ? (unsigned)(data[offset++] & 7u) : 0u;
    for (item = 0; item < item_count && size - offset >= 3u; ++item) {
        uint16_t type = fuzz_load16(data + offset);
        size_t value_bytes;
        offset += 2u;
        value_bytes = (size_t)data[offset++];
        if (value_bytes > size - offset) {
            value_bytes = size - offset;
        }
        status = leo2_code_identity_add_metadata(
            &identity, type, data + offset, value_bytes);
        check(leo2_code_id_status_string(status) != NULL);
        offset += value_bytes;
    }
    if (leo2_code_identity_validate(&identity) == LEO2_CODE_ID_OK) {
        roundtrip_identity(&identity, NULL, 0);
    }
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    leo2_code_identity identity;
    leo2_code_id_status status;

    if (size > LEO2_CODE_ID_MAX_BYTES) {
        return 0;
    }
    status = leo2_code_identity_deserialize(data, size, &identity);
    check(leo2_code_id_status_string(status) != NULL);
    if (status == LEO2_CODE_ID_OK) {
        roundtrip_identity(&identity, data, size);
    }
    exercise_builder(data, size);
    return 0;
}
