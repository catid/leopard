// Untimed independent polynomial/original-field oracle, leopard-79h.18.20.2.
#include "tower_avx2_probe.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace {
const uint16_t basis[16] = {
    0x0001, 0xACCA, 0x3C0E, 0x163E, 0xC582, 0xED2E, 0x914C, 0x4012,
    0x6C98, 0x10D8, 0x6A72, 0xB900, 0xFDB8, 0xFB34, 0xFF38, 0x991E
};
std::array<uint16_t, 65536> polynomial, canonical;
std::array<uint8_t, 65536> subfield;
std::array<uint16_t, 256> u_times;
std::array<uint8_t, 256> high_inverse;
uint8_t delta;

void require(bool condition, const char* message)
{
    if (!condition) { std::fprintf(stderr, "FAIL: %s\n", message); std::exit(1); }
}

// Direct carryless polynomial convolution and reduction, not logarithm tables.
uint16_t polynomial_product(uint16_t a, uint16_t b)
{
    uint32_t product = 0;
    for (unsigned i = 0; i < 16; ++i)
        if ((b >> i) & 1) product ^= uint32_t(a) << i;
    for (int i = 30; i >= 16; --i)
        if ((product >> i) & 1) product ^= uint32_t(0x1002D) << (i - 16);
    return static_cast<uint16_t>(product);
}

uint16_t oracle(uint16_t a, uint16_t b)
{
    return canonical[polynomial_product(polynomial[a], polynomial[b])];
}

uint8_t small(uint8_t a, uint8_t b) { return subfield[unsigned(a) * 256 + b]; }

uint16_t to_tower(uint16_t value)
{
    const uint8_t b = high_inverse[value >> 8];
    return uint16_t(uint8_t(value) ^ uint8_t(u_times[b])) | (uint16_t(b) << 8);
}

uint16_t from_tower(uint16_t value) { return uint8_t(value) ^ u_times[value >> 8]; }

uint16_t product(uint16_t x, uint16_t y)
{
    const uint8_t a = uint8_t(x), b = uint8_t(x >> 8);
    const uint8_t c = uint8_t(y), d = uint8_t(y >> 8);
    const uint8_t ac = small(a, c);
    const uint8_t bd = small(b, small(delta, d));
    const uint8_t cross = small(a ^ b, c ^ d);
    return uint16_t(ac ^ bd) | (uint16_t(cross ^ ac) << 8);
}

uint8_t nibble(const uint8_t lo[16], const uint8_t hi[16], uint8_t value)
{
    return lo[value & 15] ^ hi[value >> 4];
}

void conversion_tables(TowerConvertTables& forward, TowerConvertTables& inverse)
{
    for (unsigned i = 0; i < 16; ++i) {
        for (unsigned half = 0; half < 2; ++half) {
            const unsigned value = i << (4 * half);
            const uint16_t f = to_tower(uint16_t(value << 8));
            const uint16_t r = u_times[value];
            forward.row[half][i] = uint8_t(f);
            forward.row[half + 2][i] = uint8_t(f >> 8);
            inverse.row[half][i] = uint8_t(r);
            inverse.row[half + 2][i] = uint8_t(r >> 8);
        }
    }
    for (unsigned i = 0; i < 256; ++i) {
        for (const TowerConvertTables* table : {&forward, &inverse}) {
            const uint16_t mapped = nibble(table->row[0], table->row[1], uint8_t(i)) |
                (uint16_t(nibble(table->row[2], table->row[3], uint8_t(i))) << 8);
            const uint16_t wanted = table == &forward ? to_tower(uint16_t(i << 8)) : u_times[i];
            require(mapped == wanted, "conversion nibble linearity");
        }
    }
}

void product_tables(uint16_t coefficient, TowerProductTables& tables)
{
    const uint16_t tower = to_tower(coefficient);
    const uint8_t c = uint8_t(tower), d = uint8_t(tower >> 8);
    const uint8_t multipliers[3] = {c, small(delta, d), uint8_t(c ^ d)};
    for (unsigned m = 0; m < 3; ++m)
        for (unsigned n = 0; n < 16; ++n) {
            tables.row[2*m][n] = small(uint8_t(n), multipliers[m]);
            tables.row[2*m+1][n] = small(uint8_t(n << 4), multipliers[m]);
        }
}

void store(uint8_t* p, unsigned i, uint16_t value)
{
    p[(i / 32) * 64 + (i % 32)] = uint8_t(value);
    p[(i / 32) * 64 + (i % 32) + 32] = uint8_t(value >> 8);
}
uint16_t load(const uint8_t* p, unsigned i)
{
    return p[(i / 32) * 64 + (i % 32)] |
        (uint16_t(p[(i / 32) * 64 + (i % 32) + 32]) << 8);
}
uint32_t random_state = 20260909;
uint16_t next_random()
{
    random_state ^= random_state << 13;
    random_state ^= random_state >> 17;
    random_state ^= random_state << 5;
    return uint16_t(random_state);
}
} // namespace

int main(int argc, char** argv)
{
    if (argc != 1) { std::fprintf(stderr, "No timing/CLI options supported: %s\n", argv[1]); return 2; }
    require(leo_init() == 0, "native library initialization");
    std::array<bool, 65536> seen{};
    for (unsigned x = 0; x < 65536; ++x) {
        uint16_t p = 0;
        for (unsigned bit = 0; bit < 16; ++bit) if ((x >> bit) & 1) p ^= basis[bit];
        require(!seen[p], "Cantor basis invertible");
        seen[p] = true; polynomial[x] = p; canonical[p] = uint16_t(x);
    }
    for (unsigned a = 0; a < 256; ++a) {
        for (unsigned b = 0; b < 256; ++b) {
            const uint16_t actual = oracle(uint16_t(a), uint16_t(b));
            require(actual < 256, "first eight coordinates form closed subfield");
            require(actual == leopard::ff16::MultiplyElements(uint16_t(a), uint16_t(b)), "subfield/native oracle");
            subfield[a*256+b] = uint8_t(actual);
        }
        uint16_t frobenius = uint16_t(a);
        for (unsigned i = 0; i < 8; ++i) frobenius = oracle(frobenius, frobenius);
        require(frobenius == a, "GF256 Frobenius");
    }
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b) {
            uint8_t linear = 0;
            for (unsigned bit = 0; bit < 8; ++bit) if ((a >> bit) & 1) linear ^= small(uint8_t(1U << bit), uint8_t(b));
            require(linear == small(uint8_t(a), uint8_t(b)), "GF256 multiplier linearity");
        }
    const uint16_t relation = oracle(256, 256) ^ 256;
    require(relation < 256 && relation != 0, "quadratic relation in subfield");
    delta = uint8_t(relation);
    for (unsigned a = 0; a < 256; ++a) require((small(uint8_t(a), uint8_t(a)) ^ a) != delta, "quadratic irreducible");
    std::array<bool, 256> high_seen{};
    for (unsigned b = 0; b < 256; ++b) {
        u_times[b] = oracle(256, uint16_t(b));
        const unsigned hi = u_times[b] >> 8;
        require(!high_seen[hi], "tower high map invertible");
        high_seen[hi] = true; high_inverse[hi] = uint8_t(b);
    }
    for (unsigned x = 0; x < 65536; ++x) {
        require(from_tower(to_tower(uint16_t(x))) == x, "canonical roundtrip");
        require(to_tower(from_tower(uint16_t(x))) == x, "tower roundtrip");
        uint16_t linear = 0;
        for (unsigned bit = 0; bit < 16; ++bit) if ((x >> bit) & 1) linear ^= to_tower(uint16_t(1U << bit));
        require(linear == to_tower(uint16_t(x)), "tower conversion GF2 linearity");
    }
    TowerConvertTables forward{}, inverse{};
    conversion_tables(forward, inverse);
    // Exhaustive vector conversion on all symbols, including in-place operation.
    std::array<uint8_t, 128> input{}, mapped{}, restored{};
    for (unsigned base = 0; base < 65536; base += 64) {
        for (unsigned i = 0; i < 64; ++i) store(input.data(), i, uint16_t(base + i));
        tower_convert_blocks(input.data(), mapped.data(), 2, &forward);
        for (unsigned i = 0; i < 64; ++i) require(load(mapped.data(), i) == to_tower(uint16_t(base+i)), "vector conversion");
        tower_convert_blocks(mapped.data(), restored.data(), 2, &inverse);
        require(input == restored, "vector conversion roundtrip");
        tower_convert_blocks(input.data(), input.data(), 2, &forward);
        require(input == mapped, "in-place forward conversion");
        tower_convert_blocks(input.data(), input.data(), 2, &inverse);
        require(input == restored, "in-place inverse conversion");
    }
    unsigned mutation_rejections = 0;
    TowerProductTables tables{};
    for (unsigned coefficient = 0; coefficient < 65536; ++coefficient) {
        product_tables(uint16_t(coefficient), tables);
        // Offset by one with prefix/suffix canaries: one exact unaligned block.
        std::array<uint8_t, 66> canonical_input{}, tower_input{}, tower_output{}, canonical_output{};
        canonical_input.fill(0xA7); tower_input.fill(0xA7);
        tower_output.fill(0xA7); canonical_output.fill(0xA7);
        for (unsigned i = 0; i < 32; ++i) {
            const uint16_t x = i < 16 ? uint16_t(1U << i) : (i == 16 ? 0 : (i == 17 ? 65535 : next_random()));
            store(canonical_input.data()+1, i, x);
        }
        tower_convert_blocks(canonical_input.data()+1, tower_input.data()+1, 1, &forward);
        tower_product_blocks(tower_input.data()+1, tower_output.data()+1, 1, &tables);
        tower_convert_blocks(tower_output.data()+1, canonical_output.data()+1, 1, &inverse);
        for (unsigned i = 0; i < 32; ++i) {
            const uint16_t x = load(canonical_input.data()+1, i);
            const uint16_t expected = oracle(x, uint16_t(coefficient));
            require(expected == leopard::ff16::MultiplyElements(x, uint16_t(coefficient)), "polynomial/native full product");
            require(from_tower(product(to_tower(x), to_tower(uint16_t(coefficient)))) == expected, "scalar tower product");
            require(load(canonical_output.data()+1, i) == expected, "vector tower product");
        }
        tower_product_blocks(tower_input.data()+1, tower_input.data()+1, 1, &tables);
        require(tower_input == tower_output, "in-place product");
        require(canonical_output.front() == 0xA7 && canonical_output.back() == 0xA7 &&
                tower_input.front() == 0xA7 && tower_input.back() == 0xA7 &&
                tower_output.front() == 0xA7 && tower_output.back() == 0xA7, "vector write canaries");
        if (coefficient == 257) {
            // Table/representation corruption must actually be detected.
            for (unsigned row = 0; row < 6; ++row) {
                TowerProductTables bad = tables; bad.row[row][1] ^= 1;
                std::array<uint8_t, 64> all_nibbles{}, good{}, corrupt{};
                for (unsigned i = 0; i < 32; ++i) {
                    const uint8_t value = uint8_t((i & 15) * 17);
                    all_nibbles[i] = i < 16 ? value : 0;
                    all_nibbles[i+32] = i < 16 ? 0 : value;
                }
                tower_product_blocks(all_nibbles.data(), good.data(), 1, &tables);
                tower_product_blocks(all_nibbles.data(), corrupt.data(), 1, &bad);
                require(good != corrupt, "corrupt product table rejected");
                ++mutation_rejections;
            }
        }
    }
    tower_product_blocks(nullptr, nullptr, 0, nullptr);
    tower_convert_blocks(nullptr, nullptr, 0, nullptr);
    std::printf("{\"tracker\":\"leopard-79h.18.20.2\",\"polynomial\":65581,\"u\":256,\"delta\":%u,"
                "\"subfield_products\":65536,\"conversion_symbols\":65536,\"constant_basis_products\":1048576,"
                "\"additional_products\":1048576,\"vector_product_checks\":2097152,\"mutated_table_rejections\":%u,"
                "\"table_bytes_per_multiplier\":96,\"conversion_table_bytes_each\":64,\"timed\":false,"
                "\"production_integrated\":false,\"u_times_basis\":[", unsigned(delta), mutation_rejections);
    for (unsigned i = 0; i < 8; ++i) std::printf("%s%u", i ? "," : "", unsigned(u_times[1U << i]));
    std::printf("]}\n");
    return 0;
}
