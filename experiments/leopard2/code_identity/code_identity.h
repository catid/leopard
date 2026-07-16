/*
 * Experimental Leopard2 serialized code identity.
 *
 * This is an independent C99 implementation of the Experiment W proposal.  It
 * is deliberately isolated from the public Leopard API and is not a frozen ABI
 * or wire-format promise.
 */
#ifndef LEOPARD2_EXPERIMENT_CODE_IDENTITY_H
#define LEOPARD2_EXPERIMENT_CODE_IDENTITY_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LEO2_CODE_ID_HEADER_BYTES 36u
#define LEO2_CODE_ID_MAX_BYTES 65535u
#define LEO2_CODE_ID_MAX_METADATA 64u
#define LEO2_CODE_ID_MAX_METADATA_VALUE 4096u

#define LEO2_CODE_ID_PROFILE_LEGACY_HIGH 1u
#define LEO2_CODE_ID_PROFILE_LOW 2u
#define LEO2_CODE_ID_FIELD_GF8 1u
#define LEO2_CODE_ID_FIELD_GF16 2u

#define LEO2_CODE_ID_META_PROFILE_PARAMETERS_SHA256 UINT16_C(0x8001)
#define LEO2_CODE_ID_META_COORDINATE_SET_SHA256 UINT16_C(0x8002)
#define LEO2_CODE_ID_META_SHORTENING_SET_SHA256 UINT16_C(0x8003)
#define LEO2_CODE_ID_META_PUNCTURING_SET_SHA256 UINT16_C(0x8004)
#define LEO2_CODE_ID_META_SHARD_LAYOUT UINT16_C(0x8005)

/* Absence of META_SHARD_LAYOUT is the one canonical spelling of native V1. */
#define LEO2_CODE_ID_SHARD_LAYOUT_NATIVE_V1 UINT8_C(0)
#define LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1 UINT8_C(1)

typedef enum leo2_code_id_status {
    LEO2_CODE_ID_OK = 0,
    LEO2_CODE_ID_INVALID_ARGUMENT,
    LEO2_CODE_ID_INVALID_IDENTITY,
    LEO2_CODE_ID_UNSUPPORTED,
    LEO2_CODE_ID_NONCANONICAL,
    LEO2_CODE_ID_TRUNCATED,
    LEO2_CODE_ID_BUFFER_TOO_SMALL
} leo2_code_id_status;

typedef struct leo2_code_id_metadata {
    uint16_t type;
    uint16_t value_bytes;
    const uint8_t *value;
} leo2_code_id_metadata;

typedef struct leo2_code_identity {
    uint8_t profile;
    uint8_t profile_version;
    uint8_t field;
    uint8_t field_version;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t parent_count;
    uint32_t padded_side;
    uint16_t coordinate_map_version;
    uint16_t metadata_count;
    leo2_code_id_metadata metadata[LEO2_CODE_ID_MAX_METADATA];
} leo2_code_identity;

/*
 * Initialize the current version of a profile and derive parent_count and
 * padded_side. Metadata starts empty and may be added with
 * leo2_code_identity_add_metadata().
 */
leo2_code_id_status leo2_code_identity_make(
    leo2_code_identity *identity,
    uint8_t profile,
    uint8_t field,
    uint32_t original_count,
    uint32_t recovery_count);

/*
 * Add a TLV while maintaining canonical type order. The value is borrowed and
 * must remain valid through serialization. Duplicate base types are rejected.
 */
leo2_code_id_status leo2_code_identity_add_metadata(
    leo2_code_identity *identity,
    uint16_t type,
    const uint8_t *value,
    size_t value_bytes);

/* Validate all semantic and canonical invariants. */
leo2_code_id_status leo2_code_identity_validate(
    const leo2_code_identity *identity);

/*
 * Serialize into caller-owned storage. required_bytes is always populated for
 * a valid identity. Passing output == NULL and output_capacity == 0 is a size
 * query. No allocation is performed. Output must not overlap identity or any
 * borrowed metadata value.
 */
leo2_code_id_status leo2_code_identity_serialize(
    const leo2_code_identity *identity,
    uint8_t *output,
    size_t output_capacity,
    size_t *required_bytes);

/*
 * Parse a canonical identity without allocation. Metadata values are borrowed
 * views into input, which must remain alive while identity is used. Input must
 * not overlap identity.
 */
leo2_code_id_status leo2_code_identity_deserialize(
    const uint8_t *input,
    size_t input_bytes,
    leo2_code_identity *identity);

const char *leo2_code_id_status_string(leo2_code_id_status status);

#ifdef __cplusplus
}
#endif

#endif /* LEOPARD2_EXPERIMENT_CODE_IDENTITY_H */
