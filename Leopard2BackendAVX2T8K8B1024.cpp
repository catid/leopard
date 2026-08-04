/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See Leopard2Backend.h for the full BSD license notice.
*/

#include "Leopard2Backend.h"
#include "LeopardCommon.h"

#include <immintrin.h>

namespace leopard { namespace backend {

#if defined(LEO2_HAVE_AVX2_BACKEND) && !defined(NO_LEO_HAS_FF8) && \
    !defined(LEO2_GFNI_VARIANT)

namespace {

static const uint8_t* CachedTables = NULL;

static LEO_FORCE_INLINE const uint8_t* Tables()
{
    return CachedTables;
}

static LEO_FORCE_INLINE __m256i Broadcast(const uint8_t* table)
{
    return _mm256_broadcastsi128_si256(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(table)));
}

static LEO_FORCE_INLINE __m256i Product(
    __m256i data,
    __m256i low_table,
    __m256i high_table)
{
    const __m256i mask = _mm256_set1_epi8(15);
    const __m256i low = _mm256_shuffle_epi8(
        low_table, _mm256_and_si256(data, mask));
    const __m256i high = _mm256_shuffle_epi8(
        high_table, _mm256_and_si256(_mm256_srli_epi64(data, 4), mask));
    return _mm256_xor_si256(low, high);
}

static LEO_FORCE_INLINE const uint8_t* Table(
    const uint8_t* tables,
    uint16_t log)
{
    const uint8_t* table = tables + static_cast<unsigned>(log) * 32U;
#if defined(__GNUC__) || defined(__clang__)
    /* Keep invariant table loads inside the byte loop.  Hoisting all seven
       legacy factors exhausts AVX2's sixteen vector registers and spills the
       live transform values.  The pointer remains a compile-time displacement
       from the immutable table base, while this empty dependency prevents
       loop-invariant vector-table retention. */
    __asm__ __volatile__("" : "+r"(table));
#endif
    return table;
}

static LEO_FORCE_INLINE void MulAddPrepared(
    __m256i& destination,
    __m256i source,
    __m256i low_table,
    __m256i high_table)
{
    destination = _mm256_xor_si256(
        destination, Product(source, low_table, high_table));
}

static LEO_FORCE_INLINE void MulAdd(
    const uint8_t* tables,
    __m256i& destination,
    __m256i source,
    uint16_t log)
{
    const uint8_t* const table = Table(tables, log);
    MulAddPrepared(
        destination, source, Broadcast(table), Broadcast(table + 16));
}

static LEO_FORCE_INLINE __m256i Multiply(
    const uint8_t* tables,
    __m256i source,
    uint16_t log)
{
    const uint8_t* const table = Table(tables, log);
    return Product(source, Broadcast(table), Broadcast(table + 16));
}

static LEO_FORCE_INLINE void IFFT2Nonzero(
    const uint8_t* tables,
    __m256i& x,
    __m256i& y,
    uint16_t log)
{
    y = _mm256_xor_si256(y, x);
    MulAdd(tables, x, y, log);
}

static LEO_FORCE_INLINE void FFT2Nonzero(
    const uint8_t* tables,
    __m256i& x,
    __m256i& y,
    uint16_t log)
{
    MulAdd(tables, x, y, log);
    y = _mm256_xor_si256(y, x);
}

static LEO_FORCE_INLINE void FFT2Prepared(
    __m256i& x,
    __m256i& y,
    __m256i low_table,
    __m256i high_table)
{
    MulAddPrepared(x, y, low_table, high_table);
    y = _mm256_xor_si256(y, x);
}

static LEO_FORCE_INLINE void IFFT4Nonzero(
    const uint8_t* tables,
    __m256i& x0,
    __m256i& x1,
    __m256i& x2,
    __m256i& x3,
    uint16_t log01,
    uint16_t log23,
    uint16_t log02)
{
    IFFT2Nonzero(tables, x0, x1, log01);
    IFFT2Nonzero(tables, x2, x3, log23);
    x2 = _mm256_xor_si256(x2, x0);
    x3 = _mm256_xor_si256(x3, x1);
    const uint8_t* const table = Table(tables, log02);
    const __m256i low_table = Broadcast(table);
    const __m256i high_table = Broadcast(table + 16);
    MulAddPrepared(x0, x2, low_table, high_table);
    MulAddPrepared(x1, x3, low_table, high_table);
}

static LEO_FORCE_INLINE void IFFTDistance4Nonzero(
    __m256i& x0,
    __m256i& x1,
    __m256i& x2,
    __m256i& x3,
    __m256i& x4,
    __m256i& x5,
    __m256i& x6,
    __m256i& x7,
    __m256i low_table,
    __m256i high_table)
{
    x4 = _mm256_xor_si256(x4, x0);
    x5 = _mm256_xor_si256(x5, x1);
    x6 = _mm256_xor_si256(x6, x2);
    x7 = _mm256_xor_si256(x7, x3);
    MulAddPrepared(x0, x4, low_table, high_table);
    MulAddPrepared(x1, x5, low_table, high_table);
    MulAddPrepared(x2, x6, low_table, high_table);
    MulAddPrepared(x3, x7, low_table, high_table);
}

static LEO_FORCE_INLINE void ForwardKnownZeros(
    const uint8_t* tables,
    __m256i& x0,
    __m256i& x1,
    __m256i& x2,
    __m256i& x3,
    __m256i& x4,
    __m256i& x5,
    __m256i& x6,
    __m256i& x7,
    __m256i low_85,
    __m256i high_85)
{
    /* skew[4], skew[2], and skew[1] are the zero-product sentinel 255
       for this versioned legacy parity coset. */
    x4 = _mm256_xor_si256(x4, x0);
    x5 = _mm256_xor_si256(x5, x1);
    x6 = _mm256_xor_si256(x6, x2);
    x7 = _mm256_xor_si256(x7, x3);
    x2 = _mm256_xor_si256(x2, x0);
    x3 = _mm256_xor_si256(x3, x1);

    MulAddPrepared(x4, x6, low_85, high_85);
    MulAddPrepared(x5, x7, low_85, high_85);
    x6 = _mm256_xor_si256(x6, x4);
    x7 = _mm256_xor_si256(x7, x5);

    x1 = _mm256_xor_si256(x1, x0);
    FFT2Prepared(x2, x3, low_85, high_85);
    FFT2Nonzero(tables, x4, x5, 17);
    FFT2Nonzero(tables, x6, x7, 34);
}

/*
    The public API describes eight independent shard objects.  The packed
    terminal proves their numeric addresses are adjacent, but ISO C++ pointer
    arithmetic may not step from one allocation into another.  This x86-only
    leaf therefore performs cross-shard displacement arithmetic on uintptr_t
    values and converts only the final, validated address to a vector pointer.
    GCC, Clang, and MSVC all use the required flat-address interpretation on
    the x86 targets for which this translation unit is built.
*/
static LEO_FORCE_INLINE __m256i LoadAddress(uintptr_t address)
{
    return _mm256_loadu_si256(
        reinterpret_cast<const __m256i*>(address));
}

static LEO_FORCE_INLINE __m256i LoadAddressFresh(uintptr_t address)
{
#if defined(__GNUC__) || defined(__clang__)
    /* The rank circuit deliberately trades repeated L1-hot input reads for
       fewer table products.  Do not let the compiler retain all seven inputs
       and spill the seven output accumulators. */
    __asm__ __volatile__("" : "+r"(address) : : "memory");
#endif
    return LoadAddress(address);
}

static LEO_FORCE_INLINE __m256i Rank12Product(
    __m256i data,
    __m256i low_table,
    __m256i high_table,
    __m256i mask)
{
#if defined(__GNUC__) || defined(__clang__)
    __m256i high;
    __asm__(
        "vpsrlq $4, %[data], %[high]\n\t"
        "vpand %[mask], %[data], %[data]\n\t"
        "vpand %[mask], %[high], %[high]\n\t"
        "vpshufb %[data], %[low], %[data]\n\t"
        "vpshufb %[high], %[high_table], %[high]\n\t"
        "vpxor %[high], %[data], %[data]"
        : [data] "+x"(data), [high] "=&x"(high)
        : [low] "x"(low_table), [high_table] "x"(high_table),
          [mask] "x"(mask));
    return data;
#else
    static_cast<void>(mask);
    return Product(data, low_table, high_table);
#endif
}

static LEO_FORCE_INLINE void StoreAddress(uintptr_t address, __m256i value)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i*>(address), value);
}

} // namespace

bool InitializeAVX2FF8HighEncodeK8R8T8B1024(const void* tables)
{
    if (!tables)
        return false;
    if (CachedTables)
        return CachedTables == tables;
    CachedTables = static_cast<const uint8_t*>(tables);
    return true;
}

#if defined(_MSC_VER)
# define LEO2_AVX2_K8_B1024_ENTRY __declspec(noinline)
#elif defined(__GNUC__) && !defined(__clang__) && defined(__ELF__)
# define LEO2_AVX2_K8_B1024_ENTRY \
    __attribute__((noinline, noipa, \
        section(".text.leo2_t8_k8_b1024"), aligned(64)))
#elif defined(__clang__) && defined(__ELF__)
# define LEO2_AVX2_K8_B1024_ENTRY \
    __attribute__((noinline, section(".text.leo2_t8_k8_b1024"), aligned(64)))
#elif defined(__GNUC__) || defined(__clang__)
# define LEO2_AVX2_K8_B1024_ENTRY __attribute__((noinline, aligned(64)))
#else
# define LEO2_AVX2_K8_B1024_ENTRY
#endif

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__) && \
    !defined(_WIN32)
# define LEO2_AVX2_RANK12_PRODUCT_ASM \
    "vpsrlq $4, %%ymm14, %%ymm15\n\t" \
    "vpand %%ymm13, %%ymm14, %%ymm14\n\t" \
    "vpand %%ymm13, %%ymm15, %%ymm15\n\t" \
    "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t" \
    "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t" \
    "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
#endif

LEO2_AVX2_K8_B1024_ENTRY void AVX2FF8HighEncodeK7R7T8B1024(
    const void* const* data,
    void* const* work)
{
    const uint8_t* const tables = Tables();
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;
    const uintptr_t input_base = reinterpret_cast<uintptr_t>(data[0]);
    const uintptr_t output_base = reinterpret_cast<uintptr_t>(work[0]);
    for (unsigned lane = 0; lane < 7; ++lane)
    {
        LEO_DEBUG_ASSERT(reinterpret_cast<uintptr_t>(data[lane]) ==
            input_base + lane * 1024U);
        LEO_DEBUG_ASSERT(reinterpret_cast<uintptr_t>(work[lane]) ==
            output_base + lane * 1024U);
    }

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__) && \
    !defined(_WIN32)
    uintptr_t asm_input = input_base;
    uintptr_t asm_output = output_base;
    __asm__ __volatile__(
        "movabsq $0x0f0f0f0f0f0f0f0f, %%r11\n\t"
        "vmovq %%r11, %%xmm13\n\t"
        "vpbroadcastq %%xmm13, %%ymm13\n\t"
        "movl $32, %%r11d\n\t"
        ".p2align 4\n"
        "1:\n\t"
        "vmovdqu 0(%[input]), %%ymm7\n\t"
        "vmovdqu 2048(%[input]), %%ymm8\n\t"
        "vmovdqu 4096(%[input]), %%ymm9\n\t"
        "vmovdqu 5120(%[input]), %%ymm10\n\t"

        "vmovdqa %%ymm9, %%ymm0\n\t"
        "vpxor 6144(%[input]), %%ymm0, %%ymm0\n\t"
        "vmovdqa %%ymm0, %%ymm2\n\t"
        "vpxor 3072(%[input]), %%ymm2, %%ymm2\n\t"
        "vpxor 1024(%[input]), %%ymm0, %%ymm0\n\t"
        "vmovdqa %%ymm7, %%ymm1\n\t"
        "vpxor %%ymm10, %%ymm1, %%ymm1\n\t"
        "vmovdqa %%ymm8, %%ymm3\n\t"
        "vpxor %%ymm10, %%ymm3, %%ymm3\n\t"
        "vmovdqa %%ymm7, %%ymm6\n\t"
        "vpxor %%ymm8, %%ymm6, %%ymm6\n\t"
        "vmovdqa %%ymm6, %%ymm4\n\t"
        "vpxor %%ymm10, %%ymm4, %%ymm4\n\t"
        "vmovdqu 1024(%[input]), %%ymm5\n\t"
        "vpxor 3072(%[input]), %%ymm5, %%ymm5\n\t"
        "vpxor %%ymm9, %%ymm5, %%ymm5\n\t"

        "vbroadcasti128 0xaa0(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0xab0(%[tables]), %%ymm12\n\t"

        "vmovdqa %%ymm6, %%ymm14\n\t"
        "vpxor 3072(%[input]), %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm3, %%ymm3\n\t"

        "vmovdqa %%ymm1, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"

        "vmovdqu 1024(%[input]), %%ymm14\n\t"
        "vpxor %%ymm8, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm5, %%ymm5\n\t"

        "vmovdqa %%ymm8, %%ymm14\n\t"
        "vpxor %%ymm10, %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm5, %%ymm5\n\t"
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"

        "vbroadcasti128 0x220(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0x230(%[tables]), %%ymm12\n\t"

        "vmovdqa %%ymm7, %%ymm14\n\t"
        "vpxor %%ymm8, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm5, %%ymm5\n\t"

        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"

        "vmovdqu 1024(%[input]), %%ymm14\n\t"
        "vpxor 3072(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"
        "vpxor %%ymm14, %%ymm5, %%ymm5\n\t"
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"

        "vmovdqa %%ymm7, %%ymm14\n\t"
        "vpxor %%ymm8, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm10, %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm14, %%ymm3, %%ymm3\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"

        "vbroadcasti128 0x1320(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0x1330(%[tables]), %%ymm12\n\t"

        "vmovdqa %%ymm7, %%ymm14\n\t"
        "vpxor %%ymm8, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"

        "vmovdqu 1024(%[input]), %%ymm14\n\t"
        "vpxor 3072(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm14, %%ymm3, %%ymm3\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"

        "vmovdqu 1024(%[input]), %%ymm14\n\t"
        "vpxor %%ymm8, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm10, %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"
        "vpxor %%ymm14, %%ymm5, %%ymm5\n\t"
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"

        "vmovdqa %%ymm9, %%ymm14\n\t"
        "vpxor %%ymm10, %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        "vpsrlq $4, %%ymm14, %%ymm15\n\t"
        "vpand %%ymm13, %%ymm14, %%ymm14\n\t"
        "vpand %%ymm13, %%ymm15, %%ymm15\n\t"
        "vpshufb %%ymm14, %%ymm11, %%ymm14\n\t"
        "vpshufb %%ymm15, %%ymm12, %%ymm15\n\t"
        "vpxor %%ymm15, %%ymm14, %%ymm14\n\t"
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"

        "vmovdqu %%ymm0, 0(%[output])\n\t"
        "vmovdqu %%ymm1, 1024(%[output])\n\t"
        "vmovdqu %%ymm2, 2048(%[output])\n\t"
        "vmovdqu %%ymm3, 3072(%[output])\n\t"
        "vmovdqu %%ymm4, 4096(%[output])\n\t"
        "vmovdqu %%ymm5, 5120(%[output])\n\t"
        "vmovdqu %%ymm6, 6144(%[output])\n\t"
        "addq $32, %[input]\n\t"
        "addq $32, %[output]\n\t"
        "decl %%r11d\n\t"
        "jnz 1b"
        : [input] "+r"(asm_input), [output] "+r"(asm_output)
        : [tables] "r"(tables)
        : "cc", "memory", "r11",
          "xmm0", "xmm1", "xmm2", "xmm3",
          "xmm4", "xmm5", "xmm6", "xmm7",
          "xmm8", "xmm9", "xmm10", "xmm11",
          "xmm12", "xmm13", "xmm14", "xmm15");
    return;
#endif

    uintptr_t input = input_base;
    uintptr_t output = output_base;
    const uintptr_t input_end = input_base + 1024U;
    for (; input != input_end; input += 32U, output += 32U)
    {
        const __m256i c0 = LoadAddress(input);
        const __m256i c2 = LoadAddress(input + 2048U);
        const __m256i c4 = LoadAddress(input + 4096U);
        const __m256i mask = _mm256_set1_epi8(15);

        /* Exact rank-12 representation of the shortened K=R=7 legacy
           generator matrix.  Four source rows remain resident; every other
           source occurrence is intentionally reloaded from the L1-hot slab
           so seven output accumulators and one prepared table pair fit in
           AVX2's sixteen architectural vector registers. */
        __m256i y0 = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U), c4);
        y0 = _mm256_xor_si256(
            y0, LoadAddressFresh(input + 6144U));
        const __m256i initial5 = LoadAddressFresh(input + 5120U);
        __m256i y1 = _mm256_xor_si256(c0, initial5);
        __m256i y2 = _mm256_xor_si256(
            LoadAddressFresh(input + 3072U), c4);
        y2 = _mm256_xor_si256(
            y2, LoadAddressFresh(input + 6144U));
        __m256i y3 = _mm256_xor_si256(c2, initial5);
        __m256i y4 = _mm256_xor_si256(c0, c2);
        y4 = _mm256_xor_si256(y4, initial5);
        __m256i y5 = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U),
            LoadAddressFresh(input + 3072U));
        y5 = _mm256_xor_si256(y5, c4);
        __m256i y6 = _mm256_xor_si256(c0, c2);

        const uint8_t* table = Table(tables, 85);
        __m256i low_table = Broadcast(table);
        __m256i high_table = Broadcast(table + 16);
        __m256i source = _mm256_xor_si256(c0, c2);
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 3072U));
        __m256i product = Rank12Product(
            source, low_table, high_table, mask);
        y0 = _mm256_xor_si256(y0, product);
        y2 = _mm256_xor_si256(y2, product);
        y3 = _mm256_xor_si256(y3, product);

        source = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U), c2);
        source = _mm256_xor_si256(source, c4);
        product = Rank12Product(source, low_table, high_table, mask);
        y0 = _mm256_xor_si256(y0, product);
        y1 = _mm256_xor_si256(y1, product);
        y2 = _mm256_xor_si256(y2, product);
        y5 = _mm256_xor_si256(y5, product);

        source = _mm256_xor_si256(c0, c4);
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 5120U));
        product = Rank12Product(source, low_table, high_table, mask);
        y1 = _mm256_xor_si256(y1, product);
        y2 = _mm256_xor_si256(y2, product);
        y4 = _mm256_xor_si256(y4, product);

        source = _mm256_xor_si256(
            c2, LoadAddressFresh(input + 5120U));
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 6144U));
        product = Rank12Product(source, low_table, high_table, mask);
        y2 = _mm256_xor_si256(y2, product);
        y5 = _mm256_xor_si256(y5, product);
        y6 = _mm256_xor_si256(y6, product);

        table = Table(tables, 17);
        low_table = Broadcast(table);
        high_table = Broadcast(table + 16);
        source = _mm256_xor_si256(c0, c2);
        product = Rank12Product(source, low_table, high_table, mask);
        y5 = _mm256_xor_si256(y5, product);

        product = Rank12Product(
            LoadAddressFresh(input + 5120U),
            low_table, high_table, mask);
        y0 = _mm256_xor_si256(y0, product);
        y2 = _mm256_xor_si256(y2, product);

        source = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U),
            LoadAddressFresh(input + 3072U));
        source = _mm256_xor_si256(source, c4);
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 6144U));
        product = Rank12Product(source, low_table, high_table, mask);
        y0 = _mm256_xor_si256(y0, product);
        y2 = _mm256_xor_si256(y2, product);
        y4 = _mm256_xor_si256(y4, product);
        y5 = _mm256_xor_si256(y5, product);
        y6 = _mm256_xor_si256(y6, product);

        source = _mm256_xor_si256(c0, c2);
        source = _mm256_xor_si256(source, c4);
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 5120U));
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 6144U));
        product = Rank12Product(source, low_table, high_table, mask);
        y1 = _mm256_xor_si256(y1, product);
        y3 = _mm256_xor_si256(y3, product);
        y4 = _mm256_xor_si256(y4, product);
        y6 = _mm256_xor_si256(y6, product);

        table = Table(tables, 153);
        low_table = Broadcast(table);
        high_table = Broadcast(table + 16);
        source = _mm256_xor_si256(c0, c2);
        source = _mm256_xor_si256(source, c4);
        product = Rank12Product(source, low_table, high_table, mask);
        y0 = _mm256_xor_si256(y0, product);
        y2 = _mm256_xor_si256(y2, product);
        y4 = _mm256_xor_si256(y4, product);

        source = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U),
            LoadAddressFresh(input + 3072U));
        source = _mm256_xor_si256(source, c4);
        product = Rank12Product(source, low_table, high_table, mask);
        y1 = _mm256_xor_si256(y1, product);
        y3 = _mm256_xor_si256(y3, product);
        y4 = _mm256_xor_si256(y4, product);

        source = _mm256_xor_si256(
            LoadAddressFresh(input + 1024U), c2);
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 5120U));
        product = Rank12Product(source, low_table, high_table, mask);
        y4 = _mm256_xor_si256(y4, product);
        y5 = _mm256_xor_si256(y5, product);
        y6 = _mm256_xor_si256(y6, product);

        source = _mm256_xor_si256(
            c4, LoadAddressFresh(input + 5120U));
        source = _mm256_xor_si256(
            source, LoadAddressFresh(input + 6144U));
        product = Rank12Product(source, low_table, high_table, mask);
        y1 = _mm256_xor_si256(y1, product);
        y2 = _mm256_xor_si256(y2, product);
        y4 = _mm256_xor_si256(y4, product);
        y6 = _mm256_xor_si256(y6, product);

        StoreAddress(output, y0);
        StoreAddress(output + 1024U, y1);
        StoreAddress(output + 2048U, y2);
        StoreAddress(output + 3072U, y3);
        StoreAddress(output + 4096U, y4);
        StoreAddress(output + 5120U, y5);
        StoreAddress(output + 6144U, y6);
    }
}

LEO2_AVX2_K8_B1024_ENTRY void AVX2FF8HighEncodeK8R8T8B1024(
    const void* const* data,
    void* const* work)
{
    const uint8_t* const tables = Tables();
    LEO_DEBUG_ASSERT(tables != NULL);
    if (!tables)
        return;
    const uintptr_t input_base = reinterpret_cast<uintptr_t>(data[0]);
    const uintptr_t output_base = reinterpret_cast<uintptr_t>(work[0]);
    for (unsigned lane = 0; lane < 8; ++lane)
    {
        LEO_DEBUG_ASSERT(reinterpret_cast<uintptr_t>(data[lane]) ==
            input_base + lane * 1024U);
        LEO_DEBUG_ASSERT(reinterpret_cast<uintptr_t>(work[lane]) ==
            output_base + lane * 1024U);
    }

#if (defined(__GNUC__) || defined(__clang__)) && defined(__x86_64__) && \
    !defined(_WIN32)
    uintptr_t asm_input = input_base;
    uintptr_t asm_output = output_base;
    /*
        Exact rank-12 circuit for the fixed legacy K=R=8 map.  Let M2, M4,
        and M8 denote multiplication by the raw-basis table constants whose
        logs are 85, 17, and 153 respectively.  The twelve product sources
        are:

          M2:  1357, 0457, 3467, 2567
          M4:  13,   46,   0257, 134567
          M8:  1347, 1257, 0357, 4567

        where, for example, 1357 means I1 ^ I3 ^ I5 ^ I7.  The final
        product fanout is evaluated as a shared XOR DAG rather than updating
        every destination independently.  Together with the input-side CSE,
        this lowers the loop from 100 to 79 XOR instructions while retaining
        exactly 24 table shuffles and no spills.  The independent mature-
        transform startup KAT checks the complete byte map before this entry
        is used.
    */
    __asm__ __volatile__(
        "movabsq $0x0f0f0f0f0f0f0f0f, %%r11\n\t"
        "vmovq %%r11, %%xmm13\n\t"
        "vpbroadcastq %%xmm13, %%ymm13\n\t"
        "movl $32, %%r11d\n\t"
        ".p2align 4\n"
        "1:\n\t"
        /* u = I5 ^ I7 is reused by all three multiplier groups. */
        "vmovdqu 5120(%[input]), %%ymm10\n\t"
        "vpxor 7168(%[input]), %%ymm10, %%ymm10\n\t"
        "vmovdqu 1024(%[input]), %%ymm8\n\t"
        "vpxor 3072(%[input]), %%ymm8, %%ymm8\n\t" /* a = I1 ^ I3. */
        "vmovdqu 4096(%[input]), %%ymm9\n\t"
        "vpxor 6144(%[input]), %%ymm9, %%ymm9\n\t" /* b = I4 ^ I6. */

        /* raw2 products p0..p3. */
        "vbroadcasti128 0xaa0(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0xab0(%[tables]), %%ymm12\n\t"
        "vmovdqa %%ymm8, %%ymm7\n\t"
        "vpxor %%ymm10, %%ymm7, %%ymm7\n\t" /* s0 = a ^ u. */
        "vmovdqa %%ymm7, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm0\n\t"
        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 0(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 4096(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm1\n\t"
        "vmovdqa %%ymm9, %%ymm14\n\t"
        "vpxor 3072(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 7168(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm2\n\t"
        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 2048(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm3\n\t"

        /* raw4 products p4..p7 and their shared correction nodes. */
        "vbroadcasti128 0x220(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0x230(%[tables]), %%ymm12\n\t"
        "vmovdqa %%ymm8, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm4\n\t"
        "vmovdqa %%ymm9, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm14, %%ymm5\n\t"
        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 0(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 2048(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vpxor %%ymm14, %%ymm2, %%ymm2\n\t" /* q12 = p2 ^ p6. */
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t" /* q17 = p4 ^ p6. */
        "vpxor %%ymm2, %%ymm5, %%ymm5\n\t"  /* q18 = p5 ^ q12. */
        "vmovdqa %%ymm7, %%ymm14\n\t"
        "vpxor %%ymm9, %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm0, %%ymm6\n\t"       /* y1 starts with p0. */
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t" /* q13 = p0 ^ p7. */
        "vmovdqa %%ymm14, %%ymm7\n\t"      /* y7 starts with p7. */

        /* raw8 products p8..p11.  Fold the shared correction DAG while each
           product is live so no product or output correction spills. */
        "vbroadcasti128 0x1320(%[tables]), %%ymm11\n\t"
        "vbroadcasti128 0x1330(%[tables]), %%ymm12\n\t"
        "vmovdqa %%ymm8, %%ymm14\n\t"
        "vpxor 4096(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 7168(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vmovdqa %%ymm1, %%ymm8\n\t"       /* y4 starts with p1. */
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t" /* q14 = p1 ^ p8. */
        "vmovdqa %%ymm14, %%ymm9\n\t"      /* y2 starts with p8. */
        "vpxor %%ymm1, %%ymm5, %%ymm5\n\t" /* q20 = q14 ^ q18. */

        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 1024(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 2048(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vpxor %%ymm14, %%ymm0, %%ymm0\n\t" /* q16 = q13 ^ p9. */
        "vpxor %%ymm4, %%ymm8, %%ymm8\n\t"  /* Preserve q17 in y4. */
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t" /* y6 = q17 ^ p9. */

        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 0(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 3072(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vpxor %%ymm14, %%ymm8, %%ymm8\n\t"
        "vpxor %%ymm14, %%ymm7, %%ymm7\n\t"
        "vpxor %%ymm0, %%ymm14, %%ymm14\n\t" /* q19 = p10 ^ q16. */
        "vpxor %%ymm14, %%ymm1, %%ymm1\n\t" /* y0 = q14 ^ q19. */
        "vpxor %%ymm14, %%ymm9, %%ymm9\n\t"
        "vpxor %%ymm2, %%ymm0, %%ymm0\n\t"  /* y5 = q16 ^ q12. */
        "vpxor %%ymm2, %%ymm7, %%ymm7\n\t"

        "vmovdqa %%ymm10, %%ymm14\n\t"
        "vpxor 4096(%[input]), %%ymm14, %%ymm14\n\t"
        "vpxor 6144(%[input]), %%ymm14, %%ymm14\n\t"
        LEO2_AVX2_RANK12_PRODUCT_ASM
        "vpxor %%ymm14, %%ymm8, %%ymm8\n\t"
        "vpxor %%ymm3, %%ymm14, %%ymm14\n\t" /* q15 = p3 ^ p11. */
        "vpxor %%ymm14, %%ymm6, %%ymm6\n\t"
        "vpxor %%ymm14, %%ymm9, %%ymm9\n\t"
        "vpxor %%ymm14, %%ymm4, %%ymm4\n\t"
        "vpxor %%ymm5, %%ymm3, %%ymm3\n\t"  /* y3 = p3 ^ q20. */
        "vpxor %%ymm5, %%ymm6, %%ymm6\n\t"

        /* Add the eight three-input base terms.  Four pairwise common
           subexpressions cost 19 XORs because u is still live. */
        "vmovdqu 4096(%[input]), %%ymm2\n\t"
        "vpxor 6144(%[input]), %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm2, %%ymm1, %%ymm1\n\t"
        "vpxor %%ymm2, %%ymm9, %%ymm9\n\t"
        "vpxor 1024(%[input]), %%ymm1, %%ymm1\n\t"
        "vpxor 3072(%[input]), %%ymm9, %%ymm9\n\t"
        "vpxor %%ymm10, %%ymm6, %%ymm6\n\t"
        "vpxor %%ymm10, %%ymm3, %%ymm3\n\t"
        "vpxor 0(%[input]), %%ymm6, %%ymm6\n\t"
        "vpxor 2048(%[input]), %%ymm3, %%ymm3\n\t"
        "vmovdqu 0(%[input]), %%ymm2\n\t"
        "vpxor 2048(%[input]), %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm2, %%ymm8, %%ymm8\n\t"
        "vpxor %%ymm2, %%ymm4, %%ymm4\n\t"
        "vpxor 5120(%[input]), %%ymm8, %%ymm8\n\t"
        "vpxor 7168(%[input]), %%ymm4, %%ymm4\n\t"
        "vmovdqu 1024(%[input]), %%ymm2\n\t"
        "vpxor 3072(%[input]), %%ymm2, %%ymm2\n\t"
        "vpxor %%ymm2, %%ymm0, %%ymm0\n\t"
        "vpxor %%ymm2, %%ymm7, %%ymm7\n\t"
        "vpxor 4096(%[input]), %%ymm0, %%ymm0\n\t"
        "vpxor 6144(%[input]), %%ymm7, %%ymm7\n\t"

        "vmovdqu %%ymm1, 0(%[output])\n\t"
        "vmovdqu %%ymm6, 1024(%[output])\n\t"
        "vmovdqu %%ymm9, 2048(%[output])\n\t"
        "vmovdqu %%ymm3, 3072(%[output])\n\t"
        "vmovdqu %%ymm8, 4096(%[output])\n\t"
        "vmovdqu %%ymm0, 5120(%[output])\n\t"
        "vmovdqu %%ymm4, 6144(%[output])\n\t"
        "vmovdqu %%ymm7, 7168(%[output])\n\t"
        "addq $32, %[input]\n\t"
        "addq $32, %[output]\n\t"
        "decl %%r11d\n\t"
        "jnz 1b"
        : [input] "+r"(asm_input), [output] "+r"(asm_output)
        : [tables] "r"(tables)
        : "cc", "memory", "r11",
          "xmm0", "xmm1", "xmm2", "xmm3",
          "xmm4", "xmm5", "xmm6", "xmm7",
          "xmm8", "xmm9", "xmm10", "xmm11",
          "xmm12", "xmm13", "xmm14", "xmm15");
    return;
#endif

    uintptr_t input = input_base;
    uintptr_t output = output_base;
    const uintptr_t input_end = input_base + 1024U;
    for (; input != input_end; input += 32U, output += 32U)
    {
        __m256i x0 = LoadAddress(input);
        __m256i x1 = LoadAddress(input + 1024U);
        __m256i x2 = LoadAddress(input + 2048U);
        __m256i x3 = LoadAddress(input + 3072U);
        __m256i x4 = LoadAddress(input + 4096U);
        __m256i x5 = LoadAddress(input + 5120U);
        __m256i x6 = LoadAddress(input + 6144U);
        __m256i x7 = LoadAddress(input + 7168U);

        /*
            Exact 15-product circuit for the fixed legacy K=R=8 map.  The
            ordinary butterfly composition needs 17 table products.  Two
            independently checked linear identities jointly reconstruct the
            products conventionally named p7/p9 and p14:

              p0^p1^p4^p5^p7^p9 = input4^input5^input6^input7
              p3^p4^p5^p11^p13^p14 = input0^input1^input6^input7

            q carries the first identity only through the inverse join; b
            carries the second through the forward distance-two layer.
        */
        x1 = _mm256_xor_si256(x1, x0);
        __m256i product = Multiply(tables, x1, 153);
        x0 = _mm256_xor_si256(x0, product);
        __m256i q = product;

        x3 = _mm256_xor_si256(x3, x2);
        product = Multiply(tables, x3, 102);
        x2 = _mm256_xor_si256(x2, product);
        q = _mm256_xor_si256(q, product);

        x2 = _mm256_xor_si256(x2, x0);
        MulAdd(tables, x0, x2, 17);
        x3 = _mm256_xor_si256(x3, x1);
        MulAdd(tables, x1, x3, 17);

        x5 = _mm256_xor_si256(x5, x4);
        product = Multiply(tables, x5, 51);
        x4 = _mm256_xor_si256(x4, product);
        q = _mm256_xor_si256(q, product);
        __m256i b = product;

        x7 = _mm256_xor_si256(x7, x6);
        product = Multiply(tables, x7, 187);
        x6 = _mm256_xor_si256(x6, product);
        q = _mm256_xor_si256(q, product);
        b = _mm256_xor_si256(b, product);
        b = _mm256_xor_si256(b, x1);
        b = _mm256_xor_si256(b, x7);

        x6 = _mm256_xor_si256(x6, x4);
        MulAdd(tables, x4, x6, 34);

        __m256i joint = _mm256_xor_si256(x1, x7);
        joint = _mm256_xor_si256(joint, q);
        x7 = _mm256_xor_si256(x7, x5);

        const uint8_t* const table_85 = Table(tables, 85);
        const __m256i low_85 = Broadcast(table_85);
        const __m256i high_85 = Broadcast(table_85 + 16);
        /* p9 = M2(joint) ^ joint and the next state is joint ^ p9,
           hence x5 is exactly M2(joint).  Writing the product directly into
           the now-dead old x5 shortens two live ranges in this register-tight
           AVX2 leaf. */
        x5 = Product(joint, low_85, high_85);
        x1 = _mm256_xor_si256(x1, joint);
        x1 = _mm256_xor_si256(x1, x5);

        x4 = _mm256_xor_si256(x4, x0);
        x6 = _mm256_xor_si256(x6, x2);
        x7 = _mm256_xor_si256(x7, x3);
        product = Product(x4, low_85, high_85);
        x0 = _mm256_xor_si256(x0, product);
        product = Product(x6, low_85, high_85);
        x2 = _mm256_xor_si256(x2, product);
        product = Product(x7, low_85, high_85);
        x3 = _mm256_xor_si256(x3, product);
        b = _mm256_xor_si256(b, product);

        /* Forward transform: the first layer and the low pair are zero-skew
           XORs.  Capture p13 while its product is live, then b is p14. */
        x4 = _mm256_xor_si256(x4, x0);
        x5 = _mm256_xor_si256(x5, x1);
        x6 = _mm256_xor_si256(x6, x2);
        x7 = _mm256_xor_si256(x7, x3);
        x2 = _mm256_xor_si256(x2, x0);
        x3 = _mm256_xor_si256(x3, x1);

        product = Product(x6, low_85, high_85);
        x4 = _mm256_xor_si256(x4, product);
        product = Product(x7, low_85, high_85);
        x5 = _mm256_xor_si256(x5, product);
        b = _mm256_xor_si256(b, product);
        x6 = _mm256_xor_si256(x6, x4);
        x7 = _mm256_xor_si256(x7, x5);

        x1 = _mm256_xor_si256(x1, x0);
        x2 = _mm256_xor_si256(x2, b);
        x3 = _mm256_xor_si256(x3, x2);
        FFT2Nonzero(tables, x4, x5, 17);
        FFT2Nonzero(tables, x6, x7, 34);

        StoreAddress(output, x0);
        StoreAddress(output + 1024U, x1);
        StoreAddress(output + 2048U, x2);
        StoreAddress(output + 3072U, x3);
        StoreAddress(output + 4096U, x4);
        StoreAddress(output + 5120U, x5);
        StoreAddress(output + 6144U, x6);
        StoreAddress(output + 7168U, x7);
    }
}

#undef LEO2_AVX2_RANK12_PRODUCT_ASM
#undef LEO2_AVX2_K8_B1024_ENTRY

#else

bool InitializeAVX2FF8HighEncodeK8R8T8B1024(const void*)
{
    return false;
}

void AVX2FF8HighEncodeK8R8T8B1024(
    const void* const*,
    void* const*)
{
}

void AVX2FF8HighEncodeK7R7T8B1024(
    const void* const*,
    void* const*)
{
}

#endif

}} // namespace leopard::backend
