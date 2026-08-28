/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "Leopard2Direct.h"
#include "leopard2.h"

#ifdef LEO2_ENABLE_TEST_HOOKS
#error "production AUTO GF16/GFNI test must link the hookless archive"
#endif

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void* aligned_scratch(std::vector<uint8_t>& storage)
{
    const uintptr_t unaligned = reinterpret_cast<uintptr_t>(storage.data());
    const uintptr_t aligned =
        (unaligned + leo2_scratch_alignment() - 1U) &
        ~(static_cast<uintptr_t>(leo2_scratch_alignment()) - 1U);
    return reinterpret_cast<void*>(aligned);
}

void run_production_archive_route()
{
    if (!leopard::backend::IsCalibratedAutoGF16GFNIEncodeHost())
    {
        std::printf("Production AUTO GF16 GFNI route skipped: host is not "
            "calibrated\n");
        return;
    }

    leo2_context_options automatic_options;
    std::memset(&automatic_options, 0, sizeof(automatic_options));
    automatic_options.struct_size = sizeof(automatic_options);
    automatic_options.backend = LEO2_BACKEND_AUTO;
    automatic_options.thread_count = 1;
    leo2_context* automatic = NULL;
    require(leo2_context_create(&automatic_options, &automatic) ==
            LEO2_SUCCESS && automatic,
        "production AUTO context creation failed");
    if (leo2_context_backend(automatic) != LEO2_BACKEND_AVX2)
    {
        leo2_context_destroy(automatic);
        std::printf("Production AUTO GF16 GFNI route skipped: build is "
            "forced below AVX2\n");
        return;
    }

    leo2_context_options avx2_options = automatic_options;
    avx2_options.backend = LEO2_BACKEND_AVX2;
    leo2_context* avx2 = NULL;
    require(leo2_context_create(&avx2_options, &avx2) == LEO2_SUCCESS &&
            avx2 && leo2_context_backend(avx2) == LEO2_BACKEND_AVX2,
        "production explicit AVX2 context creation failed");

    leo2_codec* candidate = NULL;
    leo2_codec* comparator = NULL;
    require(leo2_codec_create(automatic, 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL,
            &candidate) == LEO2_SUCCESS && candidate,
        "production AUTO GF16 codec creation failed");
    require(leo2_codec_create(avx2, 1000, 200,
            LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16, NULL,
            &comparator) == LEO2_SUCCESS && comparator,
        "production explicit AVX2 codec creation failed");

    static const uint64_t kBytes = 64U * 1024U;
    require(leopard2_internal::AutoGF16GFNIEncodeModeForDiagnostics() == 1U,
        "production archive AUTO GF16 GFNI mode is not default-on");
    require(leopard2_internal::AutoGF16GFNIEncodeAvailableForDiagnostics(
                candidate) &&
            leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                candidate, kBytes),
        "production archive exact codec did not select qualified GFNI ops");
    require(!leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                candidate, kBytes - 2U) &&
            !leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                candidate, kBytes + 2U),
        "production archive selector escaped its exact byte cell");

    static const unsigned kDistinctOriginals = 17;
    std::vector<uint8_t> original(
        static_cast<size_t>(kDistinctOriginals) * kBytes);
    for (unsigned shard = 0; shard < kDistinctOriginals; ++shard)
        for (size_t i = 0; i < kBytes; ++i)
            original[static_cast<size_t>(shard) * kBytes + i] =
                static_cast<uint8_t>(
                    i * 41U + (i >> 7U) + shard * 29U);
    std::vector<const void*> original_ptrs(1000);
    for (size_t shard = 0; shard < original_ptrs.size(); ++shard)
        original_ptrs[shard] = original.data() +
            (shard % kDistinctOriginals) * kBytes;
    std::vector<uint8_t> automatic_recovery(200U * kBytes);
    std::vector<uint8_t> avx2_recovery(200U * kBytes);
    std::vector<uint8_t> batch_recovery(200U * kBytes);
    std::vector<void*> automatic_ptrs(200);
    std::vector<void*> avx2_ptrs(200);
    std::vector<void*> batch_ptrs(200);
    for (unsigned shard = 0; shard < 200; ++shard)
    {
        automatic_ptrs[shard] =
            automatic_recovery.data() + shard * kBytes;
        avx2_ptrs[shard] = avx2_recovery.data() + shard * kBytes;
        batch_ptrs[shard] = batch_recovery.data() + shard * kBytes;
    }

    size_t candidate_scratch_bytes = 0;
    size_t comparator_scratch_bytes = 0;
    require(leo2_encode_scratch_size(
            candidate, kBytes, &candidate_scratch_bytes) == LEO2_SUCCESS &&
            leo2_encode_scratch_size(
                comparator, kBytes, &comparator_scratch_bytes) ==
                LEO2_SUCCESS &&
            candidate_scratch_bytes == comparator_scratch_bytes &&
            candidate_scratch_bytes != 0,
        "production archive route changed scratch geometry");
    std::vector<uint8_t> scratch_storage(
        candidate_scratch_bytes + leo2_scratch_alignment());
    void* const scratch = aligned_scratch(scratch_storage);

    require(leo2_encode(candidate, kBytes, original_ptrs.data(),
            automatic_ptrs.data(), scratch, candidate_scratch_bytes) ==
            LEO2_SUCCESS,
        "production archive AUTO GF16 encode failed");
    require(leo2_encode(comparator, kBytes, original_ptrs.data(),
            avx2_ptrs.data(), scratch, comparator_scratch_bytes) ==
            LEO2_SUCCESS,
        "production archive explicit AVX2 encode failed");
    require(automatic_recovery == avx2_recovery,
        "production archive AUTO GFNI parity differs from explicit AVX2");

    leo2_encode_batch_item item;
    std::memset(&item, 0, sizeof(item));
    item.shard_bytes = kBytes;
    item.original = original_ptrs.data();
    item.recovery = batch_ptrs.data();
    item.scratch = scratch;
    item.scratch_bytes = candidate_scratch_bytes;
    require(leo2_encode_batch(candidate, &item, 1) == LEO2_SUCCESS &&
            batch_recovery == automatic_recovery,
        "production archive ordinary one-item batch changed parity bytes");
    require(leopard2_internal::AutoGF16GFNIEncodeModeForDiagnostics() == 1U &&
            leopard2_internal::AutoGF16GFNIEncodeSelectedForDiagnostics(
                candidate, kBytes) &&
            leo2_context_backend(automatic) == LEO2_BACKEND_AVX2,
        "production archive route changed its immutable public baseline");

    leo2_codec_destroy(comparator);
    leo2_codec_destroy(candidate);
    leo2_context_destroy(avx2);
    leo2_context_destroy(automatic);
    std::printf("Production AUTO GF16 GFNI route passed\n");
}

} // namespace

int main()
{
    try
    {
        run_production_archive_route();
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr,
            "Production AUTO GF16 GFNI route test failed: %s\n",
            error.what());
        return 1;
    }
}
