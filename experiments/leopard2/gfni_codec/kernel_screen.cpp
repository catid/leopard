// Kernel-level screen: Leopard nibble-table GF kernels vs GFNI affine kernels.
//
// The field content is irrelevant to throughput, so this uses an arbitrary but
// self-consistent GF(2)-linear byte map: the nibble tables are derived from the
// same 8x8 bit matrix that feeds VGF2P8AFFINEQB.  Equality of the two kernels
// on all inputs also empirically pins the GFNI operand bit order on this host.
#include <immintrin.h>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <vector>
#include <algorithm>

namespace {

typedef std::chrono::steady_clock Clock;

uint64_t Rng(uint64_t& s)
{
    s ^= s << 13; s ^= s >> 7; s ^= s << 17;
    return s;
}

struct LinearMap
{
    uint8_t column[8];      // column[i] = A * (1 << i)
    uint8_t low[16];        // nibble table over bits 0..3
    uint8_t high[16];       // nibble table over bits 4..7
    uint64_t affine;        // VGF2P8AFFINEQB operand
};

LinearMap MakeMap(uint64_t& seed)
{
    LinearMap map;
    for (unsigned i = 0; i < 8; ++i)
        map.column[i] = static_cast<uint8_t>(Rng(seed));
    for (unsigned n = 0; n < 16; ++n)
    {
        uint8_t lo = 0, hi = 0;
        for (unsigned b = 0; b < 4; ++b)
            if (n & (1u << b))
            {
                lo ^= map.column[b];
                hi ^= map.column[b + 4];
            }
        map.low[n] = lo;
        map.high[n] = hi;
    }
    map.affine = 0;
    for (unsigned out = 0; out < 8; ++out)
        for (unsigned in = 0; in < 8; ++in)
            if (map.column[in] & (1u << out))
                map.affine |= 1ULL << (8 * (7 - out) + in);
    return map;
}

uint8_t Apply(const LinearMap& map, uint8_t value)
{
    return static_cast<uint8_t>(map.low[value & 15] ^ map.high[value >> 4]);
}

// ---------------------------------------------------------------- GF8 kernels

__attribute__((target("avx2")))
void Gf8ButterflyNibble256(uint8_t* x, uint8_t* y, const LinearMap& map,
                           size_t bytes)
{
    const __m256i low_table = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.low)));
    const __m256i high_table = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.high)));
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    for (size_t offset = 0; offset < bytes; offset += 32)
    {
        __m256i xv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset));
        __m256i yv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset));
        yv = _mm256_xor_si256(yv, xv);
        const __m256i lo = _mm256_shuffle_epi8(
            low_table, _mm256_and_si256(yv, nibble_mask));
        const __m256i hi = _mm256_shuffle_epi8(high_table,
            _mm256_and_si256(_mm256_srli_epi64(yv, 4), nibble_mask));
        xv = _mm256_xor_si256(xv, _mm256_xor_si256(lo, hi));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + offset), xv);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + offset), yv);
    }
}

__attribute__((target("avx512f,avx512bw,avx512vl")))
void Gf8ButterflyNibble512Regs(uint8_t* x, uint8_t* y, const LinearMap& map,
                               size_t bytes)
{
    // Same 256-bit data path, EVEX encoding and 32 architectural registers.
    // This is what the shipped AVX-512VL backend variant does today.
    const __m256i low_table = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.low)));
    const __m256i high_table = _mm256_broadcastsi128_si256(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.high)));
    const __m256i nibble_mask = _mm256_set1_epi8(15);
    for (size_t offset = 0; offset < bytes; offset += 32)
    {
        __m256i xv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset));
        __m256i yv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset));
        yv = _mm256_xor_si256(yv, xv);
        const __m256i lo = _mm256_shuffle_epi8(
            low_table, _mm256_and_si256(yv, nibble_mask));
        const __m256i hi = _mm256_shuffle_epi8(high_table,
            _mm256_and_si256(_mm256_srli_epi64(yv, 4), nibble_mask));
        xv = _mm256_xor_si256(xv, _mm256_xor_si256(lo, hi));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + offset), xv);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + offset), yv);
    }
}

__attribute__((target("avx512f,avx512bw,avx512vl")))
void Gf8ButterflyNibble512Wide(uint8_t* x, uint8_t* y, const LinearMap& map,
                               size_t bytes)
{
    const __m512i low_table = _mm512_broadcast_i32x4(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.low)));
    const __m512i high_table = _mm512_broadcast_i32x4(
        _mm_loadu_si128(reinterpret_cast<const __m128i*>(map.high)));
    const __m512i nibble_mask = _mm512_set1_epi8(15);
    for (size_t offset = 0; offset < bytes; offset += 64)
    {
        __m512i xv = _mm512_loadu_si512(x + offset);
        __m512i yv = _mm512_loadu_si512(y + offset);
        yv = _mm512_xor_si512(yv, xv);
        const __m512i lo = _mm512_shuffle_epi8(
            low_table, _mm512_and_si512(yv, nibble_mask));
        const __m512i hi = _mm512_shuffle_epi8(high_table,
            _mm512_and_si512(_mm512_srli_epi64(yv, 4), nibble_mask));
        xv = _mm512_xor_si512(xv, _mm512_xor_si512(lo, hi));
        _mm512_storeu_si512(x + offset, xv);
        _mm512_storeu_si512(y + offset, yv);
    }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni")))
void Gf8ButterflyGfni512(uint8_t* x, uint8_t* y, const LinearMap& map,
                         size_t bytes)
{
    const __m512i matrix = _mm512_set1_epi64(
        static_cast<long long>(map.affine));
    for (size_t offset = 0; offset < bytes; offset += 64)
    {
        __m512i xv = _mm512_loadu_si512(x + offset);
        __m512i yv = _mm512_loadu_si512(y + offset);
        yv = _mm512_xor_si512(yv, xv);
        xv = _mm512_xor_si512(xv,
            _mm512_gf2p8affine_epi64_epi8(yv, matrix, 0));
        _mm512_storeu_si512(x + offset, xv);
        _mm512_storeu_si512(y + offset, yv);
    }
}

__attribute__((target("avx2,gfni")))
void Gf8ButterflyGfni256(uint8_t* x, uint8_t* y, const LinearMap& map,
                         size_t bytes)
{
    const __m256i matrix = _mm256_set1_epi64x(
        static_cast<long long>(map.affine));
    for (size_t offset = 0; offset < bytes; offset += 32)
    {
        __m256i xv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(x + offset));
        __m256i yv = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(y + offset));
        yv = _mm256_xor_si256(yv, xv);
        xv = _mm256_xor_si256(xv,
            _mm256_gf2p8affine_epi64_epi8(yv, matrix, 0));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(x + offset), xv);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(y + offset), yv);
    }
}

// --------------------------------------------------------------- GF16 kernels
// Leopard's GF16 layout splits every 64-byte block into 32 low bytes followed
// by 32 high bytes of the 16-bit symbols.

struct Linear16
{
    uint16_t column[16];      // column[i] = M * (1 << i)
    uint8_t low[4][16];       // nibble tables producing the low product byte
    uint8_t high[4][16];      // nibble tables producing the high product byte
    uint64_t affine[4];       // A00, A01, A10, A11
};

Linear16 MakeMap16(uint64_t& seed)
{
    Linear16 map;
    for (unsigned i = 0; i < 16; ++i)
        map.column[i] = static_cast<uint16_t>(Rng(seed));
    for (unsigned part = 0; part < 4; ++part)
        for (unsigned n = 0; n < 16; ++n)
        {
            uint16_t value = 0;
            for (unsigned b = 0; b < 4; ++b)
                if (n & (1u << b))
                    value ^= map.column[part * 4 + b];
            map.low[part][n] = static_cast<uint8_t>(value);
            map.high[part][n] = static_cast<uint8_t>(value >> 8);
        }
    // A[out_half][in_half] with in_half 0 = low input byte (bits 0..7).
    for (unsigned block = 0; block < 4; ++block)
    {
        const unsigned out_half = block >> 1;
        const unsigned in_half = block & 1;
        uint64_t matrix = 0;
        for (unsigned out = 0; out < 8; ++out)
            for (unsigned in = 0; in < 8; ++in)
            {
                const uint16_t product = map.column[in_half * 8 + in];
                const uint8_t byte = static_cast<uint8_t>(
                    out_half ? (product >> 8) : product);
                if (byte & (1u << out))
                    matrix |= 1ULL << (8 * (7 - out) + in);
            }
        map.affine[block] = matrix;
    }
    return map;
}

uint16_t Apply16(const Linear16& map, uint16_t value)
{
    uint16_t result = 0;
    for (unsigned part = 0; part < 4; ++part)
    {
        const unsigned n = (value >> (part * 4)) & 15;
        result ^= static_cast<uint16_t>(
            map.low[part][n] | (map.high[part][n] << 8));
    }
    return result;
}

__attribute__((target("avx2")))
void Gf16MulAddNibble256(uint8_t* out, const uint8_t* in, const Linear16& map,
                         size_t bytes)
{
    __m256i low_tables[4], high_tables[4];
    for (unsigned i = 0; i < 4; ++i)
    {
        low_tables[i] = _mm256_broadcastsi128_si256(_mm_loadu_si128(
            reinterpret_cast<const __m128i*>(map.low[i])));
        high_tables[i] = _mm256_broadcastsi128_si256(_mm_loadu_si128(
            reinterpret_cast<const __m128i*>(map.high[i])));
    }
    const __m256i mask = _mm256_set1_epi8(15);
    for (size_t offset = 0; offset < bytes; offset += 64)
    {
        const __m256i low_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(in + offset));
        const __m256i high_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(in + offset + 32));
        const __m256i nibbles[4] = {
            _mm256_and_si256(low_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(low_data, 4), mask),
            _mm256_and_si256(high_data, mask),
            _mm256_and_si256(_mm256_srli_epi64(high_data, 4), mask)
        };
        __m256i product_low = _mm256_setzero_si256();
        __m256i product_high = _mm256_setzero_si256();
        for (unsigned i = 0; i < 4; ++i)
        {
            product_low = _mm256_xor_si256(product_low,
                _mm256_shuffle_epi8(low_tables[i], nibbles[i]));
            product_high = _mm256_xor_si256(product_high,
                _mm256_shuffle_epi8(high_tables[i], nibbles[i]));
        }
        product_low = _mm256_xor_si256(product_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(out + offset)));
        product_high = _mm256_xor_si256(product_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(out + offset + 32)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(out + offset), product_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(out + offset + 32), product_high);
    }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni")))
void Gf16MulAddGfni512(uint8_t* out, const uint8_t* in, const Linear16& map,
                       size_t bytes)
{
    // Each 64-byte Leopard block holds 32 low bytes then 32 high bytes, so a
    // 512-bit register spans exactly one block's low half and high half of two
    // consecutive blocks.  Process two blocks at a time with 256-bit halves
    // gathered into one 512-bit operand pair.
    const __m256i a00 = _mm256_set1_epi64x(
        static_cast<long long>(map.affine[0]));
    const __m256i a01 = _mm256_set1_epi64x(
        static_cast<long long>(map.affine[1]));
    const __m256i a10 = _mm256_set1_epi64x(
        static_cast<long long>(map.affine[2]));
    const __m256i a11 = _mm256_set1_epi64x(
        static_cast<long long>(map.affine[3]));
    for (size_t offset = 0; offset < bytes; offset += 64)
    {
        const __m256i low_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(in + offset));
        const __m256i high_data = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(in + offset + 32));
        __m256i product_low = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(low_data, a00, 0),
            _mm256_gf2p8affine_epi64_epi8(high_data, a01, 0));
        __m256i product_high = _mm256_xor_si256(
            _mm256_gf2p8affine_epi64_epi8(low_data, a10, 0),
            _mm256_gf2p8affine_epi64_epi8(high_data, a11, 0));
        product_low = _mm256_xor_si256(product_low, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(out + offset)));
        product_high = _mm256_xor_si256(product_high, _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(out + offset + 32)));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(out + offset), product_low);
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(out + offset + 32), product_high);
    }
}

__attribute__((target("avx512f,avx512bw,avx512vl,gfni")))
void Gf16MulAddGfni512Wide(uint8_t* out, const uint8_t* in, const Linear16& map,
                           size_t bytes)
{
    // Two 64-byte blocks at a time: gather both low halves into one ZMM and
    // both high halves into another using 64-bit lane permutes.
    const __m512i a00 = _mm512_set1_epi64(
        static_cast<long long>(map.affine[0]));
    const __m512i a01 = _mm512_set1_epi64(
        static_cast<long long>(map.affine[1]));
    const __m512i a10 = _mm512_set1_epi64(
        static_cast<long long>(map.affine[2]));
    const __m512i a11 = _mm512_set1_epi64(
        static_cast<long long>(map.affine[3]));
    const __m512i gather_low = _mm512_setr_epi64(0, 1, 2, 3, 8, 9, 10, 11);
    const __m512i gather_high = _mm512_setr_epi64(4, 5, 6, 7, 12, 13, 14, 15);
    const __m512i scatter_a = _mm512_setr_epi64(0, 1, 2, 3, 8, 9, 10, 11);
    const __m512i scatter_b = _mm512_setr_epi64(4, 5, 6, 7, 12, 13, 14, 15);
    for (size_t offset = 0; offset + 128 <= bytes; offset += 128)
    {
        const __m512i block0 = _mm512_loadu_si512(in + offset);
        const __m512i block1 = _mm512_loadu_si512(in + offset + 64);
        const __m512i low_data =
            _mm512_permutex2var_epi64(block0, gather_low, block1);
        const __m512i high_data =
            _mm512_permutex2var_epi64(block0, gather_high, block1);
        __m512i product_low = _mm512_xor_si512(
            _mm512_gf2p8affine_epi64_epi8(low_data, a00, 0),
            _mm512_gf2p8affine_epi64_epi8(high_data, a01, 0));
        __m512i product_high = _mm512_xor_si512(
            _mm512_gf2p8affine_epi64_epi8(low_data, a10, 0),
            _mm512_gf2p8affine_epi64_epi8(high_data, a11, 0));
        __m512i out0 =
            _mm512_permutex2var_epi64(product_low, scatter_a, product_high);
        __m512i out1 =
            _mm512_permutex2var_epi64(product_low, scatter_b, product_high);
        out0 = _mm512_xor_si512(out0, _mm512_loadu_si512(out + offset));
        out1 = _mm512_xor_si512(out1, _mm512_loadu_si512(out + offset + 64));
        _mm512_storeu_si512(out + offset, out0);
        _mm512_storeu_si512(out + offset + 64, out1);
    }
}

// -------------------------------------------------------------------- harness

struct Buffers
{
    std::vector<uint8_t*> shard;
    size_t bytes;
    ~Buffers() { for (size_t i = 0; i < shard.size(); ++i) free(shard[i]); }
};

void Allocate(Buffers& buffers, size_t count, size_t bytes, uint64_t& seed)
{
    buffers.bytes = bytes;
    buffers.shard.resize(count);
    for (size_t i = 0; i < count; ++i)
    {
        void* pointer = NULL;
        if (posix_memalign(&pointer, 64, bytes) != 0) abort();
        uint8_t* data = static_cast<uint8_t*>(pointer);
        for (size_t j = 0; j < bytes; ++j)
            data[j] = static_cast<uint8_t>(Rng(seed));
        buffers.shard[i] = data;
    }
}

template<typename Fn>
double TimeSweep(Fn kernel, Buffers& buffers, unsigned layers, unsigned reps)
{
    const size_t count = buffers.shard.size();
    double best = 1e30;
    for (unsigned rep = 0; rep < reps; ++rep)
    {
        const Clock::time_point start = Clock::now();
        for (unsigned layer = 0; layer < layers; ++layer)
        {
            const size_t span = size_t(1) << layer;
            for (size_t base = 0; base + 2 * span <= count; base += 2 * span)
                for (size_t i = 0; i < span; ++i)
                    kernel(buffers.shard[base + i],
                           buffers.shard[base + span + i], buffers.bytes);
        }
        const double us = std::chrono::duration<double, std::micro>(
            Clock::now() - start).count();
        best = std::min(best, us);
    }
    return best;
}

bool VerifyGf8(const LinearMap& map)
{
    uint8_t x[64], y[64], xr[64], yr[64];
    for (unsigned i = 0; i < 64; ++i) { x[i] = uint8_t(i * 7 + 1); y[i] = uint8_t(i * 13 + 5); }
    memcpy(xr, x, 64); memcpy(yr, y, 64);
    for (unsigned i = 0; i < 64; ++i)
    {
        yr[i] ^= xr[i];
        xr[i] ^= Apply(map, yr[i]);
    }
    uint8_t a[64], b[64];
    memcpy(a, x, 64); memcpy(b, y, 64);
    Gf8ButterflyGfni512(a, b, map, 64);
    if (memcmp(a, xr, 64) != 0 || memcmp(b, yr, 64) != 0) return false;
    memcpy(a, x, 64); memcpy(b, y, 64);
    Gf8ButterflyNibble256(a, b, map, 64);
    if (memcmp(a, xr, 64) != 0 || memcmp(b, yr, 64) != 0) return false;
    memcpy(a, x, 64); memcpy(b, y, 64);
    Gf8ButterflyNibble512Wide(a, b, map, 64);
    if (memcmp(a, xr, 64) != 0 || memcmp(b, yr, 64) != 0) return false;
    memcpy(a, x, 64); memcpy(b, y, 64);
    Gf8ButterflyGfni256(a, b, map, 64);
    return memcmp(a, xr, 64) == 0 && memcmp(b, yr, 64) == 0;
}

bool VerifyGf16(const Linear16& map)
{
    uint8_t in[128], out[128], expect[128];
    for (unsigned i = 0; i < 128; ++i) { in[i] = uint8_t(i * 11 + 3); out[i] = uint8_t(i * 5 + 9); }
    memcpy(expect, out, 128);
    for (unsigned block = 0; block < 2; ++block)
        for (unsigned i = 0; i < 32; ++i)
        {
            const unsigned lo = block * 64 + i;
            const unsigned hi = block * 64 + 32 + i;
            const uint16_t value = uint16_t(in[lo] | (in[hi] << 8));
            const uint16_t product = Apply16(map, value);
            expect[lo] ^= uint8_t(product);
            expect[hi] ^= uint8_t(product >> 8);
        }
    uint8_t work[128];
    memcpy(work, out, 128);
    Gf16MulAddNibble256(work, in, map, 128);
    if (memcmp(work, expect, 128) != 0) { fprintf(stderr, "gf16 nibble mismatch\n"); return false; }
    memcpy(work, out, 128);
    Gf16MulAddGfni512(work, in, map, 128);
    if (memcmp(work, expect, 128) != 0) { fprintf(stderr, "gf16 gfni256 mismatch\n"); return false; }
    memcpy(work, out, 128);
    Gf16MulAddGfni512Wide(work, in, map, 128);
    if (memcmp(work, expect, 128) != 0) { fprintf(stderr, "gf16 gfni512 mismatch\n"); return false; }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    const unsigned shards = argc > 1 ? unsigned(atoi(argv[1])) : 32;
    const unsigned reps = argc > 2 ? unsigned(atoi(argv[2])) : 40;
    uint64_t seed = 0x243F6A8885A308D3ULL;
    const LinearMap map = MakeMap(seed);
    const Linear16 map16 = MakeMap16(seed);
    if (!VerifyGf8(map)) { fprintf(stderr, "GF8 kernel mismatch\n"); return 1; }
    if (!VerifyGf16(map16)) { fprintf(stderr, "GF16 kernel mismatch\n"); return 1; }
    printf("kernel equivalence: OK (GFNI operand order confirmed on host)\n");

    unsigned layers = 0;
    while ((1u << (layers + 1)) <= shards) ++layers;

    const size_t sizes[] = {256, 1024, 4096, 16384, 65536, 262144};
    printf("\nGF8 IFFT-style butterfly sweep, %u shards x %u layers\n",
           shards, layers);
    printf("%10s %12s %12s %12s %12s %12s\n", "bytes", "avx2-256",
           "avx512-256", "avx512-512", "gfni-256", "gfni-512");
    for (size_t index = 0; index < sizeof(sizes) / sizeof(sizes[0]); ++index)
    {
        const size_t bytes = sizes[index];
        Buffers buffers;
        Allocate(buffers, shards, bytes, seed);
        const double a = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf8ButterflyNibble256(x, y, map, n); }, buffers, layers, reps);
        const double b = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf8ButterflyNibble512Regs(x, y, map, n); }, buffers, layers, reps);
        const double c = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf8ButterflyNibble512Wide(x, y, map, n); }, buffers, layers, reps);
        const double d = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf8ButterflyGfni256(x, y, map, n); }, buffers, layers, reps);
        const double e = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf8ButterflyGfni512(x, y, map, n); }, buffers, layers, reps);
        printf("%10zu %12.2f %12.2f %12.2f %12.2f %12.2f   speedup gfni512/avx2=%.2fx\n",
               bytes, a, b, c, d, e, a / e);
    }

    printf("\nGF16 multiply-add sweep (%u shard pairs)\n", shards);
    printf("%10s %12s %12s %12s\n", "bytes", "avx2-nibble", "gfni-256",
           "gfni-512");
    for (size_t index = 0; index < sizeof(sizes) / sizeof(sizes[0]); ++index)
    {
        const size_t bytes = sizes[index];
        Buffers buffers;
        Allocate(buffers, shards, bytes, seed);
        const double a = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf16MulAddNibble256(x, y, map16, n); }, buffers, layers, reps);
        const double b = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf16MulAddGfni512(x, y, map16, n); }, buffers, layers, reps);
        const double c = TimeSweep([&](uint8_t* x, uint8_t* y, size_t n) {
            Gf16MulAddGfni512Wide(x, y, map16, n); }, buffers, layers, reps);
        printf("%10zu %12.2f %12.2f %12.2f   speedup gfni512/avx2=%.2fx\n",
               bytes, a, b, c, a / c);
    }
    return 0;
}
