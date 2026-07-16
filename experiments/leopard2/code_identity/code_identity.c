/* Independent C99 oracle for the experimental Leopard2 code identity. */
#include "code_identity.h"

#include <limits.h>
#include <string.h>

static const uint8_t k_magic[4] = { 'L', '2', 'I', 'D' };

static uint16_t load_be16(const uint8_t *input)
{
    return (uint16_t)(((uint16_t)input[0] << 8) | (uint16_t)input[1]);
}

static uint32_t load_be32(const uint8_t *input)
{
    return ((uint32_t)input[0] << 24) |
        ((uint32_t)input[1] << 16) |
        ((uint32_t)input[2] << 8) |
        (uint32_t)input[3];
}

static void store_be16(uint8_t *output, uint16_t value)
{
    output[0] = (uint8_t)(value >> 8);
    output[1] = (uint8_t)value;
}

static void store_be32(uint8_t *output, uint32_t value)
{
    output[0] = (uint8_t)(value >> 24);
    output[1] = (uint8_t)(value >> 16);
    output[2] = (uint8_t)(value >> 8);
    output[3] = (uint8_t)value;
}

static int known_critical(uint16_t type)
{
    return type >= LEO2_CODE_ID_META_PROFILE_PARAMETERS_SHA256 &&
        type <= LEO2_CODE_ID_META_SHARD_LAYOUT;
}

static int known_digest(uint16_t type)
{
    return type >= LEO2_CODE_ID_META_PROFILE_PARAMETERS_SHA256 &&
        type <= LEO2_CODE_ID_META_PUNCTURING_SET_SHA256;
}

static int aliases_known_critical(uint16_t type)
{
    return (type & UINT16_C(0x8000)) == 0u &&
        known_critical(type | UINT16_C(0x8000));
}

static int ceil_pow2_bounded(uint64_t value, uint32_t *result)
{
    uint64_t power = 1;
    if (value == 0 || result == NULL) {
        return 0;
    }
    while (power < value) {
        if (power > UINT64_MAX / 2u) {
            return 0;
        }
        power *= 2u;
    }
    if (power > UINT32_MAX) {
        return 0;
    }
    *result = (uint32_t)power;
    return 1;
}

static leo2_code_id_status derive_counts(
    uint8_t profile,
    uint32_t original_count,
    uint32_t recovery_count,
    uint32_t *parent_count,
    uint32_t *padded_side)
{
    uint64_t parent_input;
    uint32_t side;

    if (original_count == 0 || recovery_count == 0 ||
        parent_count == NULL || padded_side == NULL) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (profile == LEO2_CODE_ID_PROFILE_LEGACY_HIGH) {
        if (!ceil_pow2_bounded(recovery_count, &side)) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        parent_input = (uint64_t)original_count + (uint64_t)side;
    } else if (profile == LEO2_CODE_ID_PROFILE_LOW) {
        if (!ceil_pow2_bounded(original_count, &side)) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        parent_input = (uint64_t)side + (uint64_t)recovery_count;
    } else if (profile == LEO2_CODE_ID_PROFILE_EXACT_LOW) {
        parent_input = (uint64_t)original_count + (uint64_t)recovery_count;
        if (parent_input > UINT32_MAX) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        *parent_count = (uint32_t)parent_input;
        *padded_side = original_count;
        return LEO2_CODE_ID_OK;
    } else {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (!ceil_pow2_bounded(parent_input, parent_count)) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    *padded_side = side;
    return LEO2_CODE_ID_OK;
}

leo2_code_id_status leo2_code_identity_validate(
    const leo2_code_identity *identity)
{
    uint32_t derived_parent;
    uint32_t derived_side;
    uint16_t previous = 0;
    size_t index;
    leo2_code_id_status status;

    if (identity == NULL) {
        return LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    if (identity->profile != LEO2_CODE_ID_PROFILE_LEGACY_HIGH &&
        identity->profile != LEO2_CODE_ID_PROFILE_LOW &&
        identity->profile != LEO2_CODE_ID_PROFILE_EXACT_LOW) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (identity->profile_version != 1u) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (identity->field != LEO2_CODE_ID_FIELD_GF8 &&
        identity->field != LEO2_CODE_ID_FIELD_GF16) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (identity->field_version != 1u ||
        identity->coordinate_map_version != 1u) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }

    status = derive_counts(identity->profile, identity->original_count,
        identity->recovery_count, &derived_parent, &derived_side);
    if (status != LEO2_CODE_ID_OK) {
        return status;
    }
    if (derived_parent != identity->parent_count ||
        derived_side != identity->padded_side) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (derived_parent > 65536u) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (identity->field == LEO2_CODE_ID_FIELD_GF8 && derived_parent > 256u) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (identity->metadata_count > LEO2_CODE_ID_MAX_METADATA) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }

    for (index = 0; index < identity->metadata_count; ++index) {
        const leo2_code_id_metadata *item = &identity->metadata[index];
        size_t prior;
        if (item->type == 0u) {
            return LEO2_CODE_ID_NONCANONICAL;
        }
        if (index != 0u && item->type <= previous) {
            return LEO2_CODE_ID_NONCANONICAL;
        }
        previous = item->type;
        for (prior = 0; prior < index; ++prior) {
            if ((identity->metadata[prior].type & UINT16_C(0x7fff)) ==
                (item->type & UINT16_C(0x7fff))) {
                return LEO2_CODE_ID_NONCANONICAL;
            }
        }
        if (aliases_known_critical(item->type)) {
            return LEO2_CODE_ID_NONCANONICAL;
        }
        if ((item->type & UINT16_C(0x8000)) != 0u &&
            !known_critical(item->type)) {
            return LEO2_CODE_ID_UNSUPPORTED;
        }
        if (item->value_bytes > LEO2_CODE_ID_MAX_METADATA_VALUE) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        if (item->value_bytes != 0u && item->value == NULL) {
            return LEO2_CODE_ID_INVALID_ARGUMENT;
        }
        if (known_digest(item->type) && item->value_bytes != 32u) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        if (identity->profile == LEO2_CODE_ID_PROFILE_EXACT_LOW &&
            (item->type == LEO2_CODE_ID_META_SHORTENING_SET_SHA256 ||
             item->type == LEO2_CODE_ID_META_PUNCTURING_SET_SHA256)) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        if (item->type == LEO2_CODE_ID_META_SHARD_LAYOUT) {
            uint8_t layout;
            if (item->value_bytes != 1u || item->value == NULL) {
                return LEO2_CODE_ID_INVALID_IDENTITY;
            }
            layout = item->value[0];
            if (layout == LEO2_CODE_ID_SHARD_LAYOUT_NATIVE_V1) {
                return LEO2_CODE_ID_NONCANONICAL;
            }
            if (layout != LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1) {
                return LEO2_CODE_ID_UNSUPPORTED;
            }
            if (identity->field != LEO2_CODE_ID_FIELD_GF16) {
                return LEO2_CODE_ID_INVALID_IDENTITY;
            }
        }
    }
    return LEO2_CODE_ID_OK;
}

leo2_code_id_status leo2_code_identity_make(
    leo2_code_identity *identity,
    uint8_t profile,
    uint8_t field,
    uint32_t original_count,
    uint32_t recovery_count)
{
    leo2_code_id_status status;
    if (identity == NULL) {
        return LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    memset(identity, 0, sizeof(*identity));
    identity->profile = profile;
    identity->profile_version = 1u;
    identity->field = field;
    identity->field_version = 1u;
    identity->original_count = original_count;
    identity->recovery_count = recovery_count;
    identity->coordinate_map_version = 1u;
    status = derive_counts(profile, original_count, recovery_count,
        &identity->parent_count, &identity->padded_side);
    if (status != LEO2_CODE_ID_OK) {
        return status;
    }
    return leo2_code_identity_validate(identity);
}

leo2_code_id_status leo2_code_identity_add_metadata(
    leo2_code_identity *identity,
    uint16_t type,
    const uint8_t *value,
    size_t value_bytes)
{
    size_t position;
    size_t index;
    leo2_code_id_metadata item;

    if (identity == NULL || (value_bytes != 0u && value == NULL)) {
        return LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    if (type == 0u || value_bytes > LEO2_CODE_ID_MAX_METADATA_VALUE ||
        value_bytes > UINT16_MAX) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (identity->metadata_count >= LEO2_CODE_ID_MAX_METADATA) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (aliases_known_critical(type)) {
        return LEO2_CODE_ID_NONCANONICAL;
    }
    if ((type & UINT16_C(0x8000)) != 0u && !known_critical(type)) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (known_digest(type) && value_bytes != 32u) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (identity->profile == LEO2_CODE_ID_PROFILE_EXACT_LOW &&
        (type == LEO2_CODE_ID_META_SHORTENING_SET_SHA256 ||
         type == LEO2_CODE_ID_META_PUNCTURING_SET_SHA256)) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    if (type == LEO2_CODE_ID_META_SHARD_LAYOUT) {
        uint8_t layout;
        if (value_bytes != 1u) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
        layout = value[0];
        if (layout == LEO2_CODE_ID_SHARD_LAYOUT_NATIVE_V1) {
            return LEO2_CODE_ID_NONCANONICAL;
        }
        if (layout != LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1) {
            return LEO2_CODE_ID_UNSUPPORTED;
        }
        if (identity->field != LEO2_CODE_ID_FIELD_GF16) {
            return LEO2_CODE_ID_INVALID_IDENTITY;
        }
    }
    for (index = 0; index < identity->metadata_count; ++index) {
        if ((identity->metadata[index].type & UINT16_C(0x7fff)) ==
            (type & UINT16_C(0x7fff))) {
            return LEO2_CODE_ID_NONCANONICAL;
        }
    }

    position = 0;
    while (position < identity->metadata_count &&
        identity->metadata[position].type < type) {
        ++position;
    }
    for (index = identity->metadata_count; index > position; --index) {
        identity->metadata[index] = identity->metadata[index - 1u];
    }
    item.type = type;
    item.value_bytes = (uint16_t)value_bytes;
    item.value = value;
    identity->metadata[position] = item;
    ++identity->metadata_count;
    return leo2_code_identity_validate(identity);
}

leo2_code_id_status leo2_code_identity_serialize(
    const leo2_code_identity *identity,
    uint8_t *output,
    size_t output_capacity,
    size_t *required_bytes)
{
    size_t total = LEO2_CODE_ID_HEADER_BYTES;
    size_t offset;
    size_t index;
    leo2_code_id_status status;

    if (required_bytes == NULL) {
        return LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    *required_bytes = 0;
    status = leo2_code_identity_validate(identity);
    if (status != LEO2_CODE_ID_OK) {
        return status;
    }
    for (index = 0; index < identity->metadata_count; ++index) {
        total += 4u + identity->metadata[index].value_bytes;
    }
    if (total > LEO2_CODE_ID_MAX_BYTES) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }
    *required_bytes = total;
    if (output == NULL) {
        return output_capacity == 0u ? LEO2_CODE_ID_OK :
            LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    if (output_capacity < total) {
        return LEO2_CODE_ID_BUFFER_TOO_SMALL;
    }

    memcpy(output, k_magic, sizeof(k_magic));
    output[4] = 1u;
    output[5] = 0u;
    store_be16(output + 6, LEO2_CODE_ID_HEADER_BYTES);
    store_be32(output + 8, (uint32_t)total);
    output[12] = identity->profile;
    output[13] = identity->profile_version;
    output[14] = identity->field;
    output[15] = identity->field_version;
    store_be32(output + 16, identity->original_count);
    store_be32(output + 20, identity->recovery_count);
    store_be32(output + 24, identity->parent_count);
    store_be32(output + 28, identity->padded_side);
    store_be16(output + 32, identity->coordinate_map_version);
    store_be16(output + 34, identity->metadata_count);

    offset = LEO2_CODE_ID_HEADER_BYTES;
    for (index = 0; index < identity->metadata_count; ++index) {
        const leo2_code_id_metadata *item = &identity->metadata[index];
        store_be16(output + offset, item->type);
        store_be16(output + offset + 2u, item->value_bytes);
        offset += 4u;
        if (item->value_bytes != 0u) {
            memcpy(output + offset, item->value, item->value_bytes);
            offset += item->value_bytes;
        }
    }
    return LEO2_CODE_ID_OK;
}

leo2_code_id_status leo2_code_identity_deserialize(
    const uint8_t *input,
    size_t input_bytes,
    leo2_code_identity *identity)
{
    uint32_t total_bytes;
    uint16_t metadata_count;
    size_t offset;
    size_t index;
    leo2_code_id_status status;

    if (input == NULL || identity == NULL) {
        return LEO2_CODE_ID_INVALID_ARGUMENT;
    }
    memset(identity, 0, sizeof(*identity));
    if (input_bytes < LEO2_CODE_ID_HEADER_BYTES) {
        return LEO2_CODE_ID_TRUNCATED;
    }
    if (memcmp(input, k_magic, sizeof(k_magic)) != 0) {
        return LEO2_CODE_ID_NONCANONICAL;
    }
    if (input[4] != 1u) {
        return LEO2_CODE_ID_UNSUPPORTED;
    }
    if (input[5] != 0u || load_be16(input + 6) != LEO2_CODE_ID_HEADER_BYTES) {
        return LEO2_CODE_ID_NONCANONICAL;
    }
    total_bytes = load_be32(input + 8);
    if ((uint64_t)total_bytes != (uint64_t)input_bytes ||
        total_bytes > LEO2_CODE_ID_MAX_BYTES) {
        return LEO2_CODE_ID_NONCANONICAL;
    }
    metadata_count = load_be16(input + 34);
    if (metadata_count > LEO2_CODE_ID_MAX_METADATA) {
        return LEO2_CODE_ID_INVALID_IDENTITY;
    }

    identity->profile = input[12];
    identity->profile_version = input[13];
    identity->field = input[14];
    identity->field_version = input[15];
    identity->original_count = load_be32(input + 16);
    identity->recovery_count = load_be32(input + 20);
    identity->parent_count = load_be32(input + 24);
    identity->padded_side = load_be32(input + 28);
    identity->coordinate_map_version = load_be16(input + 32);
    identity->metadata_count = metadata_count;

    offset = LEO2_CODE_ID_HEADER_BYTES;
    for (index = 0; index < metadata_count; ++index) {
        uint16_t value_bytes;
        if (input_bytes - offset < 4u) {
            return LEO2_CODE_ID_TRUNCATED;
        }
        identity->metadata[index].type = load_be16(input + offset);
        value_bytes = load_be16(input + offset + 2u);
        identity->metadata[index].value_bytes = value_bytes;
        offset += 4u;
        if (value_bytes > LEO2_CODE_ID_MAX_METADATA_VALUE ||
            (size_t)value_bytes > input_bytes - offset) {
            return LEO2_CODE_ID_TRUNCATED;
        }
        identity->metadata[index].value = input + offset;
        offset += value_bytes;
    }
    if (offset != input_bytes) {
        return LEO2_CODE_ID_NONCANONICAL;
    }
    status = leo2_code_identity_validate(identity);
    if (status != LEO2_CODE_ID_OK) {
        return status;
    }
    return LEO2_CODE_ID_OK;
}

const char *leo2_code_id_status_string(leo2_code_id_status status)
{
    switch (status) {
    case LEO2_CODE_ID_OK: return "ok";
    case LEO2_CODE_ID_INVALID_ARGUMENT: return "invalid argument";
    case LEO2_CODE_ID_INVALID_IDENTITY: return "invalid identity";
    case LEO2_CODE_ID_UNSUPPORTED: return "unsupported semantic extension";
    case LEO2_CODE_ID_NONCANONICAL: return "noncanonical input";
    case LEO2_CODE_ID_TRUNCATED: return "truncated input";
    case LEO2_CODE_ID_BUFFER_TOO_SMALL: return "buffer too small";
    default: return "unknown status";
    }
}
