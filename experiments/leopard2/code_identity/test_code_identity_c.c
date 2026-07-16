/* C99 unit tests and line-protocol driver for the independent C oracle. */
#include "code_identity.h"

#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PROTOCOL_LINE_BYTES (LEO2_CODE_ID_MAX_BYTES * 2u + 4096u)

static unsigned long g_checks;

#define CHECK(expression) do {                                                \
    ++g_checks;                                                               \
    if (!(expression)) {                                                      \
        (void)fprintf(stderr, "check failed at %s:%d: %s\n",                 \
            __FILE__, __LINE__, #expression);                                 \
        return 0;                                                             \
    }                                                                         \
} while (0)

static int hex_nibble(char character, uint8_t *value)
{
    if (character >= '0' && character <= '9') {
        *value = (uint8_t)(character - '0');
        return 1;
    }
    if (character >= 'a' && character <= 'f') {
        *value = (uint8_t)(character - 'a' + 10);
        return 1;
    }
    if (character >= 'A' && character <= 'F') {
        *value = (uint8_t)(character - 'A' + 10);
        return 1;
    }
    return 0;
}

static int decode_hex(const char *hex, uint8_t *output, size_t capacity,
    size_t *output_bytes)
{
    size_t chars;
    size_t index;
    if (hex == NULL || output == NULL || output_bytes == NULL) {
        return 0;
    }
    if (strcmp(hex, "-") == 0) {
        *output_bytes = 0;
        return 1;
    }
    chars = strlen(hex);
    if ((chars & 1u) != 0u || chars / 2u > capacity) {
        return 0;
    }
    for (index = 0; index < chars / 2u; ++index) {
        uint8_t high;
        uint8_t low;
        if (!hex_nibble(hex[index * 2u], &high) ||
            !hex_nibble(hex[index * 2u + 1u], &low)) {
            return 0;
        }
        output[index] = (uint8_t)((uint8_t)(high << 4) | low);
    }
    *output_bytes = chars / 2u;
    return 1;
}

static void print_hex(const uint8_t *data, size_t bytes)
{
    static const char digits[] = "0123456789abcdef";
    size_t index;
    for (index = 0; index < bytes; ++index) {
        (void)putchar(digits[data[index] >> 4]);
        (void)putchar(digits[data[index] & 15u]);
    }
}

static int parse_u32(const char *text, uint32_t *value)
{
    char *end = NULL;
    uintmax_t parsed;
    if (text == NULL || value == NULL || text[0] == '\0' || text[0] == '-') {
        return 0;
    }
    errno = 0;
    parsed = strtoumax(text, &end, 0);
    if (errno != 0 || end == NULL || *end != '\0' || parsed > UINT32_MAX) {
        return 0;
    }
    *value = (uint32_t)parsed;
    return 1;
}

static int parse_u16(const char *text, uint16_t *value)
{
    uint32_t parsed;
    if (!parse_u32(text, &parsed) || parsed > UINT16_MAX) {
        return 0;
    }
    *value = (uint16_t)parsed;
    return 1;
}

static int parse_u8(const char *text, uint8_t *value)
{
    uint32_t parsed;
    if (!parse_u32(text, &parsed) || parsed > UINT8_MAX) {
        return 0;
    }
    *value = (uint8_t)parsed;
    return 1;
}

static int serialize_into(
    const leo2_code_identity *identity,
    uint8_t *output,
    size_t *output_bytes)
{
    size_t required = 0;
    leo2_code_id_status status = leo2_code_identity_serialize(
        identity, NULL, 0, &required);
    if (status != LEO2_CODE_ID_OK || required > LEO2_CODE_ID_MAX_BYTES) {
        return 0;
    }
    status = leo2_code_identity_serialize(identity, output,
        LEO2_CODE_ID_MAX_BYTES, &required);
    if (status != LEO2_CODE_ID_OK) {
        return 0;
    }
    *output_bytes = required;
    return 1;
}

static int deserialize_exact_prefix_rejected(
    const uint8_t *encoded,
    size_t encoded_bytes)
{
    uint8_t *exact = (uint8_t *)malloc(encoded_bytes == 0u ? 1u : encoded_bytes);
    leo2_code_identity decoded;
    leo2_code_id_status status;
    if (exact == NULL) {
        return 0;
    }
    if (encoded_bytes != 0u) {
        memcpy(exact, encoded, encoded_bytes);
    }
    status = leo2_code_identity_deserialize(exact, encoded_bytes, &decoded);
    free(exact);
    return status != LEO2_CODE_ID_OK;
}

static int test_golden(
    uint8_t profile,
    uint8_t field,
    uint32_t original_count,
    uint32_t recovery_count,
    uint16_t metadata_type,
    const char *metadata_hex,
    const char *expected_hex)
{
    leo2_code_identity identity;
    leo2_code_identity decoded;
    uint8_t metadata[LEO2_CODE_ID_MAX_METADATA_VALUE];
    uint8_t expected[LEO2_CODE_ID_MAX_BYTES];
    uint8_t encoded[LEO2_CODE_ID_MAX_BYTES];
    size_t metadata_bytes = 0;
    size_t expected_bytes = 0;
    size_t encoded_bytes = 0;
    size_t required = 0;

    CHECK(decode_hex(expected_hex, expected, sizeof(expected), &expected_bytes));
    CHECK(leo2_code_identity_make(&identity, profile, field,
        original_count, recovery_count) == LEO2_CODE_ID_OK);
    if (metadata_type != 0u) {
        CHECK(decode_hex(metadata_hex, metadata, sizeof(metadata),
            &metadata_bytes));
        CHECK(leo2_code_identity_add_metadata(&identity, metadata_type,
            metadata, metadata_bytes) == LEO2_CODE_ID_OK);
    }
    CHECK(leo2_code_identity_serialize(&identity, NULL, 0, &required) ==
        LEO2_CODE_ID_OK);
    CHECK(required == expected_bytes);
    CHECK(required != 0u);
    memset(encoded, 0xa5, sizeof(encoded));
    CHECK(leo2_code_identity_serialize(&identity, encoded, required - 1u,
        &encoded_bytes) == LEO2_CODE_ID_BUFFER_TOO_SMALL);
    CHECK(encoded_bytes == required);
    CHECK(encoded[0] == UINT8_C(0xa5));
    CHECK(encoded[required - 1u] == UINT8_C(0xa5));
    CHECK(serialize_into(&identity, encoded, &encoded_bytes));
    CHECK(encoded_bytes == expected_bytes);
    CHECK(memcmp(encoded, expected, expected_bytes) == 0);
    CHECK(leo2_code_identity_deserialize(encoded, encoded_bytes, &decoded) ==
        LEO2_CODE_ID_OK);
    CHECK(decoded.profile == identity.profile);
    CHECK(decoded.field == identity.field);
    CHECK(decoded.original_count == identity.original_count);
    CHECK(decoded.recovery_count == identity.recovery_count);
    CHECK(decoded.parent_count == identity.parent_count);
    CHECK(decoded.padded_side == identity.padded_side);
    CHECK(decoded.metadata_count == identity.metadata_count);
    return 1;
}

static int run_unit_tests(void)
{
    static const char golden_1[] =
        "4c324944010000240000002401010101000000f0000000100000010000000010"
        "00010000";
    static const char golden_2[] =
        "4c324944010000240000003b0201020100000081000000640000020000000100"
        "0001000112340013776972652d76656e646f722d6e65757472616c";
    static const char golden_3[] =
        "4c324944010000240000004801010201000003e8000000c80000080000000100"
        "000100018002002094b6af12d26e0acba766ecae613414f76b95c477d9a8d22"
        "a792ded9802a1a572";
    static const char golden_4[] =
        "4c32494401000024000000290201020100000002000001000000020000000002"
        "000100018005000101";
    static const char golden_5[] =
        "4c32494401000024000000240301010100000003000000fd0000010000000003"
        "00010000";
    leo2_code_identity identity;
    leo2_code_identity decoded;
    uint8_t encoded[LEO2_CODE_ID_MAX_BYTES];
    uint8_t digest[32];
    uint8_t one = 1u;
    size_t encoded_bytes = 0;
    size_t index;

    CHECK(leo2_code_identity_validate(NULL) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_make(NULL, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 1, 1) == LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_add_metadata(NULL, 1, &one, 1) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_deserialize(NULL, 0, &decoded) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_deserialize(digest, 0, NULL) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_serialize(NULL, NULL, 0, &encoded_bytes) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(encoded_bytes == 0u);

    CHECK(test_golden(LEO2_CODE_ID_PROFILE_LEGACY_HIGH,
        LEO2_CODE_ID_FIELD_GF8, 240, 16, 0, "-", golden_1));
    CHECK(test_golden(LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 129, 100, UINT16_C(0x1234),
        "776972652d76656e646f722d6e65757472616c", golden_2));
    CHECK(test_golden(LEO2_CODE_ID_PROFILE_LEGACY_HIGH,
        LEO2_CODE_ID_FIELD_GF16, 1000, 200,
        LEO2_CODE_ID_META_COORDINATE_SET_SHA256,
        "94b6af12d26e0acba766ecae613414f76b95c477d9a8d22a792ded9802a1a572",
        golden_3));
    CHECK(test_golden(LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 2, 256,
        LEO2_CODE_ID_META_SHARD_LAYOUT, "01", golden_4));
    CHECK(test_golden(LEO2_CODE_ID_PROFILE_EXACT_LOW,
        LEO2_CODE_ID_FIELD_GF8, 3, 253, 0, "-", golden_5));

    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 129, 100) == LEO2_CODE_ID_OK);
    memset(encoded, 0xa5, sizeof(encoded));
    CHECK(leo2_code_identity_serialize(&identity, encoded, sizeof(encoded),
        NULL) == LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(encoded[0] == UINT8_C(0xa5));
    CHECK(leo2_code_identity_serialize(&identity, NULL, 1, &encoded_bytes) ==
        LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(encoded_bytes == LEO2_CODE_ID_HEADER_BYTES);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x1234),
        NULL, 1) == LEO2_CODE_ID_INVALID_ARGUMENT);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x1234),
        &one, SIZE_MAX) == LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x1234),
        NULL, 0) == LEO2_CODE_ID_OK);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x0010),
        &one, 1) == LEO2_CODE_ID_OK);
    CHECK(identity.metadata_count == 2u);
    CHECK(identity.metadata[0].type == UINT16_C(0x0010));
    CHECK(identity.metadata[1].type == UINT16_C(0x1234));
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x9234),
        &one, 1) == LEO2_CODE_ID_UNSUPPORTED);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x9234),
        NULL, 0) == LEO2_CODE_ID_UNSUPPORTED);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x8001),
        &one, 1) == LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity, UINT16_C(0x8010),
        digest, sizeof(digest)) == LEO2_CODE_ID_UNSUPPORTED);
    for (index = 1; index <= 5; ++index) {
        CHECK(leo2_code_identity_add_metadata(&identity, (uint16_t)index,
            &one, 1) == LEO2_CODE_ID_NONCANONICAL);
    }

    CHECK(serialize_into(&identity, encoded, &encoded_bytes));
    for (index = 0; index < encoded_bytes; ++index) {
        CHECK(deserialize_exact_prefix_rejected(encoded, index));
    }
    CHECK(leo2_code_identity_deserialize(encoded, encoded_bytes, &decoded) ==
        LEO2_CODE_ID_OK);
    encoded[5] = 1u;
    CHECK(leo2_code_identity_deserialize(encoded, encoded_bytes, &decoded) ==
        LEO2_CODE_ID_NONCANONICAL);

    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 2, 256) == LEO2_CODE_ID_OK);
    one = LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, &one, 1) == LEO2_CODE_ID_OK);
    CHECK(leo2_code_identity_validate(&identity) == LEO2_CODE_ID_OK);
    CHECK(identity.metadata_count == 1u);
    CHECK(identity.metadata[0].value[0] ==
        LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1);

    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 2, 256) == LEO2_CODE_ID_OK);
    one = LEO2_CODE_ID_SHARD_LAYOUT_NATIVE_V1;
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, &one, 1) ==
        LEO2_CODE_ID_NONCANONICAL);
    one = 2u;
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, &one, 1) ==
        LEO2_CODE_ID_UNSUPPORTED);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, &one, 0) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, digest, 2) ==
        LEO2_CODE_ID_INVALID_IDENTITY);

    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF8, 2, 5) == LEO2_CODE_ID_OK);
    one = LEO2_CODE_ID_SHARD_LAYOUT_GF16_PADDED_ODD_V1;
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHARD_LAYOUT, &one, 1) ==
        LEO2_CODE_ID_INVALID_IDENTITY);

    CHECK(leo2_code_identity_make(&identity,
        LEO2_CODE_ID_PROFILE_LEGACY_HIGH, LEO2_CODE_ID_FIELD_GF8,
        129, 1) == LEO2_CODE_ID_OK);
    CHECK(identity.parent_count == 256u);
    CHECK(identity.padded_side == 1u);
    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF8, 129, 100) == LEO2_CODE_ID_UNSUPPORTED);
    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_EXACT_LOW,
        LEO2_CODE_ID_FIELD_GF8, 3, 253) == LEO2_CODE_ID_OK);
    CHECK(identity.parent_count == 256u);
    CHECK(identity.padded_side == 3u);
    memset(digest, 0, sizeof(digest));
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_COORDINATE_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHORTENING_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_PUNCTURING_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    memset(digest, 0xff, sizeof(digest));
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_COORDINATE_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_SHORTENING_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_add_metadata(&identity,
        LEO2_CODE_ID_META_PUNCTURING_SET_SHA256, digest, sizeof(digest)) ==
        LEO2_CODE_ID_INVALID_IDENTITY);
    CHECK(leo2_code_identity_make(&identity,
        LEO2_CODE_ID_PROFILE_EXACT_HIGH_RESERVED, LEO2_CODE_ID_FIELD_GF16,
        253, 3) == LEO2_CODE_ID_UNSUPPORTED);
    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 32768, 32768) == LEO2_CODE_ID_OK);
    CHECK(identity.parent_count == 65536u);
    CHECK(leo2_code_identity_make(&identity, LEO2_CODE_ID_PROFILE_LOW,
        LEO2_CODE_ID_FIELD_GF16, 65536, 1) == LEO2_CODE_ID_UNSUPPORTED);
    CHECK(strcmp(leo2_code_id_status_string(LEO2_CODE_ID_OK), "ok") == 0);

    (void)printf("C code-identity tests passed: %lu checks\n", g_checks);
    return 1;
}

static void protocol_error(void)
{
    (void)puts("ERR");
    (void)fflush(stdout);
}

static void protocol_ok(const uint8_t *encoded, size_t encoded_bytes)
{
    (void)fputs("OK ", stdout);
    print_hex(encoded, encoded_bytes);
    (void)putchar('\n');
    (void)fflush(stdout);
}

static void protocol_decode(char *cursor, uint8_t *input, uint8_t *output)
{
    char *hex = strtok(cursor, " \t\r\n");
    char *extra;
    size_t input_bytes = 0;
    size_t output_bytes = 0;
    leo2_code_identity identity;
    if (hex == NULL) {
        protocol_error();
        return;
    }
    extra = strtok(NULL, " \t\r\n");
    if (extra != NULL || !decode_hex(hex, input, LEO2_CODE_ID_MAX_BYTES,
            &input_bytes) ||
        leo2_code_identity_deserialize(input, input_bytes, &identity) !=
            LEO2_CODE_ID_OK ||
        !serialize_into(&identity, output, &output_bytes)) {
        protocol_error();
        return;
    }
    protocol_ok(output, output_bytes);
}

static void protocol_encode(char *cursor, uint8_t *values, uint8_t *output)
{
    char *token;
    uint8_t profile;
    uint8_t field;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t metadata_count;
    uint32_t index;
    size_t value_offset = 0;
    size_t output_bytes = 0;
    leo2_code_identity identity;

    token = strtok(cursor, " \t\r\n");
    if (!parse_u8(token, &profile)) { protocol_error(); return; }
    token = strtok(NULL, " \t\r\n");
    if (!parse_u8(token, &field)) { protocol_error(); return; }
    token = strtok(NULL, " \t\r\n");
    if (!parse_u32(token, &original_count)) { protocol_error(); return; }
    token = strtok(NULL, " \t\r\n");
    if (!parse_u32(token, &recovery_count)) { protocol_error(); return; }
    token = strtok(NULL, " \t\r\n");
    if (!parse_u32(token, &metadata_count) ||
        metadata_count > LEO2_CODE_ID_MAX_METADATA ||
        leo2_code_identity_make(&identity, profile, field, original_count,
            recovery_count) != LEO2_CODE_ID_OK) {
        protocol_error();
        return;
    }
    for (index = 0; index < metadata_count; ++index) {
        uint16_t type;
        size_t value_bytes;
        token = strtok(NULL, " \t\r\n");
        if (!parse_u16(token, &type)) { protocol_error(); return; }
        token = strtok(NULL, " \t\r\n");
        if (token == NULL || !decode_hex(token, values + value_offset,
                LEO2_CODE_ID_MAX_BYTES - value_offset, &value_bytes) ||
            value_bytes > LEO2_CODE_ID_MAX_METADATA_VALUE ||
            leo2_code_identity_add_metadata(&identity, type,
                values + value_offset, value_bytes) != LEO2_CODE_ID_OK) {
            protocol_error();
            return;
        }
        value_offset += value_bytes;
    }
    if (strtok(NULL, " \t\r\n") != NULL ||
        !serialize_into(&identity, output, &output_bytes)) {
        protocol_error();
        return;
    }
    protocol_ok(output, output_bytes);
}

static int run_protocol(void)
{
    char *line = (char *)malloc(PROTOCOL_LINE_BYTES);
    uint8_t *input = (uint8_t *)malloc(LEO2_CODE_ID_MAX_BYTES);
    uint8_t *output = (uint8_t *)malloc(LEO2_CODE_ID_MAX_BYTES);
    uint8_t *values = (uint8_t *)malloc(LEO2_CODE_ID_MAX_BYTES);
    if (line == NULL || input == NULL || output == NULL || values == NULL) {
        free(values);
        free(output);
        free(input);
        free(line);
        return 0;
    }
    while (fgets(line, (int)PROTOCOL_LINE_BYTES, stdin) != NULL) {
        if (line[0] == 'D' && (line[1] == ' ' || line[1] == '\t')) {
            protocol_decode(line + 2, input, output);
        } else if (line[0] == 'E' && (line[1] == ' ' || line[1] == '\t')) {
            protocol_encode(line + 2, values, output);
        } else if (strcmp(line, "Q\n") == 0 || strcmp(line, "Q\r\n") == 0) {
            break;
        } else {
            protocol_error();
        }
    }
    free(values);
    free(output);
    free(input);
    free(line);
    return ferror(stdin) == 0;
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        return run_unit_tests() ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (argc == 2 && strcmp(argv[1], "--protocol") == 0) {
        return run_protocol() ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    (void)fprintf(stderr, "usage: %s [--protocol]\n", argv[0]);
    return EXIT_FAILURE;
}
