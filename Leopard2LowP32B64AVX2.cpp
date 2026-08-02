/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"

#include <immintrin.h>
#include <stdint.h>

namespace leopard { namespace backend {

#if LEO2_EXPERIMENT_LOW_P32_B64_TERMINAL && !defined(NO_LEO_HAS_FF8)

namespace {

static const unsigned kFF8TableBytes = 32;
static const unsigned kFF8HighTableOffset = 16;

static inline __m256i BroadcastTable(const unsigned char* table)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static inline __m256i ProductVector(
    __m256i data, __m256i low_table, __m256i high_table)
{
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, nibble_mask));
    const __m256i high = _mm256_shuffle_epi8(high_table,
        _mm256_and_si256(_mm256_srli_epi64(data, 4), nibble_mask));
    return _mm256_xor_si256(low, high);
}

static inline __m256i ApplyWeight(
    const unsigned char* tables, __m256i data, uint16_t weight_log)
{
    if (weight_log == 0 || weight_log == 255)
        return data;
    const unsigned char* table = tables + weight_log * kFF8TableBytes;
    return ProductVector(
        data, BroadcastTable(table),
        BroadcastTable(table + kFF8HighTableOffset));
}

static inline uint8_t AddLogs(uint8_t first, uint8_t second)
{
    const unsigned sum = static_cast<unsigned>(first) + second;
    return static_cast<uint8_t>(sum + (sum >> 8));
}

static inline void XorRow64(
    void* destination_pointer, const void* source_pointer)
{
    uint8_t* destination = static_cast<uint8_t*>(destination_pointer);
    const uint8_t* source = static_cast<const uint8_t*>(source_pointer);
    const __m256i first = _mm256_xor_si256(
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(destination)),
        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(source)));
    const __m256i second = _mm256_xor_si256(
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination + 32)),
        _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + 32)));
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(destination), first);
    _mm256_storeu_si256(
        reinterpret_cast<__m256i*>(destination + 32), second);
}

static void FinishInverse(
    const Ops& ops,
    void* const* work,
    const uint8_t* inverse_skew,
    void* const* xor_output)
{
    // The weighted input primitive already executed the distance-one and
    // distance-two inverse layers.  Complete distance four and then the
    // distance-sixteen terminal layer of the exact P=32 transform.
    ops.ff8_ifft_butterfly4_range(
        work, 4, inverse_skew[4], inverse_skew[12], inverse_skew[8],
        64, true);
    ops.ff8_ifft_butterfly4_range(
        work + 16, 4,
        inverse_skew[20], inverse_skew[28], inverse_skew[24], 64, true);

    const uint16_t final_log = inverse_skew[16];
    if (final_log != 255)
    {
        ops.ff8_ifft_butterfly2_range(
            work, xor_output, 16, final_log, 64);
        return;
    }

    // A zero skew suppresses the multiplication but retains y ^= x.  Preserve
    // the accumulating form used for block one without materializing it.
    for (unsigned lane = 0; lane < 16; ++lane)
    {
        if (xor_output)
        {
            XorRow64(xor_output[lane], work[lane]);
            XorRow64(xor_output[lane + 16], work[lane]);
            XorRow64(xor_output[lane + 16], work[lane + 16]);
        }
        else
            XorRow64(work[lane + 16], work[lane]);
    }
}

static void FinalFFTReveal(
    const Ops& ops,
    const unsigned char* tables,
    void* const* work,
    const uint8_t* forward_skew,
    uint32_t requested_mask,
    const uint8_t* locator_logs,
    void* const* restored)
{
    ops.ff8_fft_butterfly4_range(
        work, 8,
        forward_skew[8], forward_skew[24], forward_skew[16], 64, true);
    for (unsigned row = 0; row < 32; row += 8)
    {
        ops.ff8_fft_butterfly4_range(
            work + row, 2,
            forward_skew[row + 2], forward_skew[row + 6],
            forward_skew[row + 4], 64, true);
    }

    // Fuse the last distance-one forward layer with inverse-locator
    // multiplication and scatter.  Unrequested evaluator rows never make a
    // final round trip through the work array.
    for (unsigned row = 0; row < 32; row += 2)
    {
        const uint8_t* x = static_cast<const uint8_t*>(work[row]);
        const uint8_t* y = static_cast<const uint8_t*>(work[row + 1]);
        __m256i x0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x));
        __m256i x1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + 32));
        __m256i y0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y));
        __m256i y1 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + 32));

        const uint16_t transform_log = forward_skew[row + 1];
        if (transform_log != 255)
        {
            const unsigned char* table =
                tables + transform_log * kFF8TableBytes;
            const __m256i low = BroadcastTable(table);
            const __m256i high = BroadcastTable(
                table + kFF8HighTableOffset);
            x0 = _mm256_xor_si256(x0, ProductVector(y0, low, high));
            x1 = _mm256_xor_si256(x1, ProductVector(y1, low, high));
        }
        y0 = _mm256_xor_si256(y0, x0);
        y1 = _mm256_xor_si256(y1, x1);

        if ((requested_mask & (UINT32_C(1) << row)) != 0)
        {
            const uint16_t reveal_log =
                static_cast<uint8_t>(255U - locator_logs[row]);
            x0 = ApplyWeight(tables, x0, reveal_log);
            x1 = ApplyWeight(tables, x1, reveal_log);
            uint8_t* output = static_cast<uint8_t*>(restored[row]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), x0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output + 32), x1);
        }
        if ((requested_mask & (UINT32_C(1) << (row + 1))) != 0)
        {
            const uint16_t reveal_log =
                static_cast<uint8_t>(255U - locator_logs[row + 1]);
            y0 = ApplyWeight(tables, y0, reveal_log);
            y1 = ApplyWeight(tables, y1, reveal_log);
            uint8_t* output = static_cast<uint8_t*>(restored[row + 1]);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(output), y0);
            _mm256_storeu_si256(
                reinterpret_cast<__m256i*>(output + 32), y1);
        }
    }
}

} // namespace

#if defined(_MSC_VER)
#define LEO2_LOW_P32_B64_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && !defined(__APPLE__)
#define LEO2_LOW_P32_B64_NOINLINE \
    __attribute__((noinline, section(".text.leo2_low_p32_b64_terminal")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_LOW_P32_B64_NOINLINE __attribute__((noinline))
#else
#define LEO2_LOW_P32_B64_NOINLINE
#endif

bool LEO2_LOW_P32_B64_NOINLINE AVX2FF8LowP32B64Terminal(
    const Ops& ops,
    const void* const* coordinate_data,
    const uint32_t* requested_coordinates,
    uint32_t requested_count,
    const uint8_t* locator_logs,
    uint8_t block_factor,
    const uint8_t* inverse_skew0,
    const uint8_t* inverse_skew1,
    const uint8_t* forward_skew,
    void* const* restored,
    void* const* work)
{
    // uint8_t may inspect any object's representation.  GetAVX2FF8Tables()
    // deliberately erases the mature translation unit's private table type;
    // consume only its documented 32-byte low/high object representation here
    // rather than inventing a second, alias-incompatible struct definition.
    const unsigned char* tables = static_cast<const unsigned char*>(
        GetAVX2FF8Tables());
    if (ops.kind != LEO2_BACKEND_AVX2 || !ops.ff8_weighted_ifft_butterfly4 ||
        !ops.ff8_ifft_butterfly4_range || !ops.ff8_fft_butterfly4_range ||
        !ops.ff8_ifft_butterfly2_range || !coordinate_data ||
        !requested_coordinates || requested_count == 0 ||
        requested_count > 32 || !locator_logs || !inverse_skew0 ||
        !inverse_skew1 || !forward_skew || !restored || !work || !tables)
        return false;

    uint32_t requested_mask = 0;
    for (uint32_t i = 0; i < requested_count; ++i)
    {
        const uint32_t coordinate = requested_coordinates[i];
        if (coordinate >= 32 || !restored[coordinate] ||
            (requested_mask & (UINT32_C(1) << coordinate)) != 0)
            return false;
        requested_mask |= UINT32_C(1) << coordinate;
    }

    for (unsigned block = 0; block < 2; ++block)
    {
        const unsigned base = block * 32;
        const uint8_t* inverse_skew =
            block == 0 ? inverse_skew0 : inverse_skew1;
        for (unsigned row = 0; row < 32; row += 4)
        {
            const unsigned coordinate = base + row;
            uint8_t live_mask = 0;
            uint16_t weights[4];
            for (unsigned lane = 0; lane < 4; ++lane)
            {
                const unsigned index = coordinate + lane;
                if (coordinate_data[index])
                    live_mask = static_cast<uint8_t>(
                        live_mask | (1U << lane));
                weights[lane] = block == 0
                    ? locator_logs[index]
                    : AddLogs(locator_logs[index], block_factor);
            }
            ops.ff8_weighted_ifft_butterfly4(
                coordinate_data[coordinate],
                coordinate_data[coordinate + 1],
                coordinate_data[coordinate + 2],
                coordinate_data[coordinate + 3],
                work[coordinate], work[coordinate + 1],
                work[coordinate + 2], work[coordinate + 3],
                weights[0], weights[1], weights[2], weights[3], live_mask,
                inverse_skew[row + 1], inverse_skew[row + 3],
                inverse_skew[row + 2], 64);
        }
    }

    FinishInverse(ops, work, inverse_skew0, NULL);

    // In the normalized LCH basis, Algorithm 4's block-zero contribution is
    // A + A'.  The fixed triangular XOR circuit below is the exact
    // AddFormalDerivative(P=32) schedule used by the mature decoder.
    for (unsigned i = 1; i < 32; ++i)
    {
        const unsigned width = ((i ^ (i - 1)) + 1) >> 1;
        for (unsigned lane = 0; lane < width; ++lane)
            XorRow64(work[i - width + lane], work[i + lane]);
    }

    // The block factor was folded into block one's locator weights.  Complete
    // its inverse transform and accumulate its last layer directly into A.
    FinishInverse(ops, work + 32, inverse_skew1, work);
    FinalFFTReveal(
        ops, tables, work, forward_skew,
        requested_mask, locator_logs, restored);
    return true;
}

#undef LEO2_LOW_P32_B64_NOINLINE

#endif


}} // namespace leopard::backend
