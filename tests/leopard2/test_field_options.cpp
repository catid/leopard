#include "leopard.h"
#include "leopard2.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

uint32_t expected_mask()
{
    uint32_t mask = 0;
#if LEO2_EXPECT_GF8
    mask |= LEO2_FIELD_MASK_GF8;
#endif
#if LEO2_EXPECT_GF16
    mask |= LEO2_FIELD_MASK_GF16;
#endif
    return mask;
}

struct Scratch
{
    std::vector<uint8_t> storage;
    void* pointer;

    explicit Scratch(size_t bytes)
        : storage(bytes == 0 ? 0 : bytes + leo2_scratch_alignment() - 1)
        , pointer(NULL)
    {
        if (bytes != 0)
        {
            const uintptr_t raw = reinterpret_cast<uintptr_t>(storage.data());
            const uintptr_t alignment = leo2_scratch_alignment();
            pointer = reinterpret_cast<void*>(
                (raw + alignment - 1) & ~(alignment - 1));
        }
    }
};

std::vector<uint8_t> exercise_field(
    leo2_context* context, leo2_field field, size_t bytes,
    unsigned original_count = 3, unsigned recovery_count = 2)
{
    leo2_codec* codec = NULL;
    require(leo2_codec_create(context, original_count, recovery_count,
        LEO2_PROFILE_LEGACY_HIGH_V1, field, NULL, &codec) == LEO2_SUCCESS,
        "enabled field codec creation");
    require(codec != NULL && leo2_codec_field(codec) == field,
        "enabled field identity");

    // Contiguous shards also exercise the packed-tail encode layout.
    std::vector<uint8_t> originals(original_count * bytes);
    std::vector<uint8_t> recovery(recovery_count * bytes);
    std::vector<const void*> original_pointers(original_count);
    std::vector<void*> recovery_pointers(recovery_count);
    for (size_t i = 0; i < original_count; ++i)
    {
        for (size_t j = 0; j < bytes; ++j)
            originals[i * bytes + j] =
                static_cast<uint8_t>(i * 71U + j * 29U + 3U);
        original_pointers[i] = originals.data() + i * bytes;
    }
    for (size_t i = 0; i < recovery_count; ++i)
        recovery_pointers[i] = recovery.data() + i * bytes;

    size_t encode_scratch_bytes = 0;
    require(leo2_encode_scratch_size(codec, bytes, &encode_scratch_bytes) ==
        LEO2_SUCCESS, "encode scratch query");
    Scratch encode_scratch(encode_scratch_bytes);
    require(leo2_encode(codec, bytes, original_pointers.data(),
        recovery_pointers.data(), encode_scratch.pointer,
        encode_scratch_bytes) == LEO2_SUCCESS, "enabled field encode");

    std::vector<uint8_t> original_present(original_count, 1);
    std::vector<uint8_t> recovery_present(recovery_count, 1);
    original_present[1] = 0;
    leo2_decode_plan* plan = NULL;
    require(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), &plan) == LEO2_SUCCESS,
        "enabled field decode plan");
    size_t decode_scratch_bytes = 0;
    require(leo2_decode_plan_scratch_size(plan, bytes,
        &decode_scratch_bytes) == LEO2_SUCCESS, "decode scratch query");
    Scratch decode_scratch(decode_scratch_bytes);
    std::vector<uint8_t> restored(bytes, 0);
    original_pointers[1] = NULL;
    std::vector<const void*> recovery_inputs(recovery_count);
    for (size_t i = 0; i < recovery_count; ++i)
        recovery_inputs[i] = recovery.data() + i * bytes;
    std::vector<void*> restored_pointers(original_count, NULL);
    restored_pointers[1] = restored.data();
    require(leo2_decode_plan_execute(plan, bytes, original_pointers.data(),
        recovery_inputs.data(), restored_pointers.data(), decode_scratch.pointer,
        decode_scratch_bytes) == LEO2_SUCCESS, "enabled field decode");
    require(std::memcmp(restored.data(), originals.data() + bytes, bytes) == 0,
        "enabled field recovery bytes");

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    return recovery;
}

void check_packed_tail_backends(leo2_context* auto_context)
{
#if LEO2_EXPECT_GF8
    leo2_context_options options = {};
    options.struct_size = sizeof(options);
    options.thread_count = 1;
    options.backend = LEO2_BACKEND_SCALAR;
    leo2_context* scalar = NULL;
    require(leo2_context_create(&options, &scalar) == LEO2_SUCCESS,
        "scalar context creation");

    options.backend = LEO2_BACKEND_AVX2;
    leo2_context* avx2 = NULL;
    const leo2_result result = leo2_context_create(&options, &avx2);
    require((result == LEO2_SUCCESS && avx2 != NULL) ||
            (result == LEO2_UNSUPPORTED && avx2 == NULL),
        "AVX2 context creation or unavailable backend");
    for (unsigned original_count = 10; original_count <= 11; ++original_count)
    {
        for (unsigned recovery_count = 5; recovery_count <= 8; ++recovery_count)
        {
            const std::vector<uint8_t> reference = exercise_field(
                scalar, LEO2_FIELD_GF8, 256, original_count, recovery_count);
            require(exercise_field(auto_context, LEO2_FIELD_GF8, 256,
                original_count, recovery_count) == reference,
                "AUTO packed-tail parity matches scalar");
            if (avx2 != NULL)
                require(exercise_field(avx2, LEO2_FIELD_GF8, 256,
                    original_count, recovery_count) == reference,
                    "AVX2 packed-tail parity matches scalar");
        }
    }
    std::printf("Packed-tail parity/recovery passed: 8 shapes, AVX2=%s\n",
        avx2 != NULL ? "tested" : "unavailable");
    if (avx2 != NULL)
        leo2_context_destroy(avx2);
    leo2_context_destroy(scalar);
#else
    (void)auto_context;
#endif
}

void check_codec_selection(leo2_context* context)
{
    leo2_codec* codec = NULL;
    leo2_result result = leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec);
#if LEO2_EXPECT_GF8
    require(result == LEO2_SUCCESS, "explicit GF8 should be enabled");
    leo2_codec_destroy(codec);
    exercise_field(context, LEO2_FIELD_GF8, 17);
#else
    require(result == LEO2_UNSUPPORTED && codec == NULL,
        "explicit omitted GF8 must be unsupported");
#endif

    codec = NULL;
    result = leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL, &codec);
#if LEO2_EXPECT_GF16
    require(result == LEO2_SUCCESS, "explicit GF16 should be enabled");
    leo2_codec_destroy(codec);
    exercise_field(context, LEO2_FIELD_GF16, 66);
#else
    require(result == LEO2_UNSUPPORTED && codec == NULL,
        "explicit omitted GF16 must be unsupported");
#endif

    codec = NULL;
    result = leo2_codec_create(context, 3, 2,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_AUTO, NULL, &codec);
#if LEO2_EXPECT_GF8
    require(result == LEO2_SUCCESS && leo2_codec_field(codec) == LEO2_FIELD_GF8,
        "AUTO small code must retain canonical GF8");
    leo2_codec_destroy(codec);
#else
    require(result == LEO2_UNSUPPORTED && codec == NULL,
        "AUTO must not substitute GF16 for omitted canonical GF8");
#endif

    codec = NULL;
    result = leo2_codec_create(context, 193, 65,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_AUTO, NULL, &codec);
#if LEO2_EXPECT_GF16
    require(result == LEO2_SUCCESS && leo2_codec_field(codec) == LEO2_FIELD_GF16,
        "AUTO large code must retain canonical GF16");
    leo2_codec_destroy(codec);
#else
    require(result == LEO2_UNSUPPORTED && codec == NULL,
        "AUTO must not substitute GF8 for omitted canonical GF16");
#endif

    codec = NULL;
    require(leo2_codec_create(context, 193, 65,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, NULL, &codec) ==
        LEO2_INVALID_COUNTS && codec == NULL,
        "capacity errors precede compile support");
}

LeopardResult legacy_encode(unsigned originals, unsigned recoveries)
{
    const unsigned work_count = leo_encode_work_count(originals, recoveries);
    require(work_count != 0, "legacy work count");
    std::vector<std::vector<uint8_t> > input_storage(
        originals, std::vector<uint8_t>(64, 0x5a));
    std::vector<std::vector<uint8_t> > work_storage(
        work_count, std::vector<uint8_t>(64));
    std::vector<const void*> inputs(originals);
    std::vector<void*> work(work_count);
    for (unsigned i = 0; i < originals; ++i)
        inputs[i] = input_storage[i].data();
    for (unsigned i = 0; i < work_count; ++i)
        work[i] = work_storage[i].data();
    return leo_encode(64, originals, recoveries, work_count,
        inputs.data(), work.data());
}

void check_legacy_no_loss()
{
    const unsigned originals = 3;
    const unsigned recoveries = 2;
    const unsigned work_count = leo_decode_work_count(originals, recoveries);
    std::vector<std::vector<uint8_t> > original_storage(
        originals, std::vector<uint8_t>(64));
    std::vector<std::vector<uint8_t> > work_storage(
        work_count, std::vector<uint8_t>(64));
    std::vector<const void*> original_pointers(originals);
    std::vector<const void*> recovery_pointers(recoveries, NULL);
    std::vector<void*> work_pointers(work_count);
    for (unsigned i = 0; i < originals; ++i)
    {
        for (unsigned j = 0; j < 64; ++j)
            original_storage[i][j] = static_cast<uint8_t>(i * 31U + j);
        original_pointers[i] = original_storage[i].data();
    }
    for (unsigned i = 0; i < work_count; ++i)
        work_pointers[i] = work_storage[i].data();
    require(leo_decode(64, originals, recoveries, work_count,
        original_pointers.data(), recovery_pointers.data(),
        work_pointers.data()) == Leopard_Success,
        "field-independent legacy no-loss decode");
    for (unsigned i = 0; i < originals; ++i)
        require(work_storage[i] == original_storage[i],
            "legacy no-loss output bytes");

    // K=1 has its own recovery fast path.  With no loss it must still use the
    // available original rather than dereferencing an absent recovery shard.
    std::vector<uint8_t> single_original(64);
    std::vector<uint8_t> single_work(64, 0xa5);
    for (unsigned i = 0; i < 64; ++i)
        single_original[i] = static_cast<uint8_t>(i * 17U + 3U);
    const void* single_original_pointer[1] = { single_original.data() };
    const void* single_recovery_pointer[1] = { NULL };
    void* single_work_pointer[1] = { single_work.data() };
    require(leo_decode(64, 1, 1, 1, single_original_pointer,
        single_recovery_pointer, single_work_pointer) == Leopard_Success,
        "field-independent one-original no-loss decode");
    require(single_work == single_original,
        "legacy one-original no-loss output bytes");
}

void check_legacy_selection()
{
    require(leo_init() == Leopard_Success, "legacy initialization");
    require(leo_encode_work_count(1, 1) == 1 &&
            leo_decode_work_count(1, 1) == 1 &&
            leo_encode_work_count(32768, 32768) == 65536 &&
            leo_decode_work_count(32768, 32768) == 65536,
        "legacy valid count boundaries");
    require(leo_encode_work_count(0, 1) == 0 &&
            leo_decode_work_count(2, 3) == 0 &&
            leo_encode_work_count(65536, 1) == 0 &&
            leo_decode_work_count(32769, 32767) == 0 &&
            leo_decode_work_count(
                std::numeric_limits<unsigned>::max(), 2) == 0,
        "legacy invalid count boundaries");
#if LEO2_EXPECT_GF8
    require(legacy_encode(3, 2) == Leopard_Success,
        "legacy canonical GF8 enabled");
#else
    require(legacy_encode(3, 2) == Leopard_Platform,
        "legacy must not substitute GF16 for omitted GF8");
#endif
#if LEO2_EXPECT_GF16
    require(legacy_encode(193, 65) == Leopard_Success,
        "legacy canonical GF16 enabled");
#else
    require(legacy_encode(193, 65) == Leopard_Platform,
        "legacy omitted GF16 must fail deterministically");
#endif
    require(legacy_encode(3, 1) == Leopard_Success,
        "field-independent one-parity path remains available");
    require(legacy_encode(1, 1) == Leopard_Success,
        "field-independent one-original path remains available");
    check_legacy_no_loss();
}

} // namespace

int main()
{
    try
    {
        require(leo2_context_field_mask(NULL) == 0,
            "null context field mask");
        leo2_context* context = NULL;
        require(leo2_context_create(NULL, &context) == LEO2_SUCCESS,
            "context creation");
        require(leo2_context_field_mask(context) == expected_mask(),
            "context field mask");
        check_codec_selection(context);
        check_packed_tail_backends(context);
        check_legacy_selection();
        leo2_context_destroy(context);
        std::printf("Leopard field options passed: mask=%u\n", expected_mask());
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "Leopard field options failed: %s\n", error.what());
        return 1;
    }
}
