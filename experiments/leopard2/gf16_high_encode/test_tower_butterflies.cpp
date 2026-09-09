// Untimed full-pair qualification, leopard-79h.18.20.3.
// Reuse the already-qualified polynomial/Cantor oracle and table constructor,
// but do not execute the prior standalone qualification again.
#define main PreviousTowerAlgebraMain
#include "test_tower_algebra.cpp"
#undef main
#include "tower_butterfly_probe.h"
#include "Leopard2Backend.h"
#include <algorithm>
#include <vector>
#if defined(__SANITIZE_ADDRESS__)
#include <sanitizer/asan_interface.h>
#endif

namespace {
struct Guard {
    std::vector<uint8_t> storage;
    const size_t bytes, offset;
    uint8_t* data;
    Guard(size_t n, size_t misalignment) : storage(n+192, 0xD3), bytes(n), offset(64+misalignment), data(storage.data()+offset) { poison(true); }
    ~Guard() { poison(false); }
    void poison(bool enabled) {
#if defined(__SANITIZE_ADDRESS__)
        if (enabled) {
            __asan_poison_memory_region(storage.data(), offset);
            __asan_poison_memory_region(data+bytes, storage.size()-offset-bytes);
        } else __asan_unpoison_memory_region(storage.data(), storage.size());
#else
        (void)enabled;
#endif
    }
    void check() {
        poison(false);
        for (size_t i = 0; i < offset; ++i) require(storage[i] == 0xD3, "prefix guard");
        for (size_t i = offset+bytes; i < storage.size(); ++i) require(storage[i] == 0xD3, "suffix guard");
        poison(true);
    }
};

uint16_t symbol(const uint8_t* p, size_t bytes, size_t i)
{
    const size_t base = (i/32)*64, lane = i%32;
    const size_t q = std::min(size_t(32), (bytes-base)/2);
    return p[base+lane] | (uint16_t(p[base+q+lane]) << 8);
}
void put(uint8_t* p, size_t bytes, size_t i, uint16_t value)
{
    const size_t base = (i/32)*64, lane = i%32;
    const size_t q = std::min(size_t(32), (bytes-base)/2);
    p[base+lane] = uint8_t(value); p[base+q+lane] = uint8_t(value >> 8);
}

TowerLowMapTables low_map{};
void setup()
{
    require(leo_init() == 0, "original field initialized");
    for (unsigned x = 0; x < 65536; ++x) {
        uint16_t value = 0;
        for (unsigned bit = 0; bit < 16; ++bit) if ((x >> bit) & 1) value ^= basis[bit];
        polynomial[x] = value; canonical[value] = uint16_t(x);
    }
    for (unsigned a = 0; a < 256; ++a)
        for (unsigned b = 0; b < 256; ++b) subfield[a*256+b] = uint8_t(oracle(uint16_t(a), uint16_t(b)));
    delta = 128;
    for (unsigned b = 0; b < 256; ++b) {
        u_times[b] = oracle(256, uint16_t(b));
        require(u_times[b] >> 8 == b, "actual full high-byte map is identity");
        high_inverse[b] = uint8_t(b);
    }
    for (unsigned i = 0; i < 16; ++i) {
        low_map.row[0][i] = uint8_t(u_times[i]);
        low_map.row[1][i] = uint8_t(u_times[i << 4]);
    }
    for (unsigned i = 0; i < 256; ++i)
        require(nibble(low_map.row[0], low_map.row[1], uint8_t(i)) == uint8_t(u_times[i]), "low map nibble decomposition");
    require(leopard::ff16::MultiplyLogElement(1, 0) == 1 &&
            leopard::ff16::MultiplyLogElement(1, 65535) == 1, "raw log endpoint semantics");
}

struct Fixture {
    size_t bytes;
    Guard x, y, u, v, restored_x, restored_y;
    std::vector<uint8_t> cx, cy, cu, cv, expected_x, expected_y, bx, by, bu, bv;
    Fixture(size_t n, size_t offset) : bytes(n), x(n,offset), y(n,offset), u(n,offset), v(n,offset),
        restored_x(n,offset), restored_y(n,offset), cx(n), cy(n), cu(n), cv(n),
        expected_x(n), expected_y(n), bx(n), by(n), bu(n), bv(n) {}
    void guards() { for (Guard* p : {&x,&y,&u,&v,&restored_x,&restored_y}) p->check(); }
};

void baseline(const leopard::backend::Ops& ops, unsigned mode, Fixture& f, uint16_t log)
{
    const size_t bytes = f.bytes;
    f.bx = f.cx; f.by = f.cy; f.bu = f.cu; f.bv = f.cv;
    if (log == 65535 && mode != 2) {
        // Match the documented transform-caller XOR specialization, not
        // the raw pair table at log65535 (which would multiply by one).
        for (size_t i = 0; i < bytes; ++i) {
            if (mode == 3) { f.bu[i] ^= f.bx[i]; f.bv[i] ^= f.bx[i] ^ f.by[i]; }
            else f.by[i] ^= f.bx[i];
        }
    } else if (mode == 0) ops.ff16_ifft_butterfly2(f.bx.data(), f.by.data(), log, bytes);
    else if (mode == 1) ops.ff16_fft_butterfly2(f.bx.data(), f.by.data(), log, bytes);
    else if (mode == 2) ops.ff16_fft_butterfly2_out(f.bx.data(), f.by.data(), f.bu.data(), f.bv.data(), log, bytes);
    else ops.ff16_ifft_butterfly2_xor(f.bx.data(), f.by.data(), f.bu.data(), f.bv.data(), log, bytes);
    require((mode < 2 ? f.bx : f.bu) == f.expected_x &&
            (mode < 2 ? f.by : f.bv) == f.expected_y, "original backend/polynomial butterfly");
    if (mode >= 2) require(f.bx == f.cx && f.by == f.cy, "original out inputs preserved");
}

uint64_t check_case(Fixture& f, uint16_t log, bool basis_inputs,
                    const leopard::backend::Ops& scalar, const leopard::backend::Ops& avx2)
{
    const size_t bytes = f.bytes;
    for (size_t i = 0; i < bytes/2; ++i) {
        const unsigned lane = unsigned(i%32);
        put(f.cx.data(), bytes, i, basis_inputs ? (lane < 16 ? uint16_t(1U << lane) : 0) : next_random());
        put(f.cy.data(), bytes, i, basis_inputs ? (lane >= 16 ? uint16_t(1U << (lane-16)) : 0) : next_random());
        put(f.cu.data(), bytes, i, next_random()); put(f.cv.data(), bytes, i, next_random());
    }
    const uint16_t coefficient = log == 65535 ? 0 : leopard::ff16::MultiplyLogElement(1, log);
    TowerProductTables tables{};
    product_tables(coefficient, tables);
    const TowerProductTables* selected = log == 65535 ? nullptr : &tables;
    for (unsigned mode = 0; mode < 4; ++mode) {
        tower_convert_involution(f.cx.data(), f.x.data, bytes, &low_map);
        tower_convert_involution(f.cy.data(), f.y.data, bytes, &low_map);
        tower_convert_involution(f.cu.data(), f.u.data, bytes, &low_map);
        tower_convert_involution(f.cv.data(), f.v.data, bytes, &low_map);
        for (size_t i = 0; i < bytes/2; ++i) {
            uint16_t a = symbol(f.cx.data(), bytes, i), b = symbol(f.cy.data(), bytes, i);
            if (mode == 0 || mode == 3) b ^= a;
            a ^= oracle(b, coefficient);
            if (mode == 1 || mode == 2) b ^= a;
            if (mode == 3) { a ^= symbol(f.cu.data(), bytes, i); b ^= symbol(f.cv.data(), bytes, i); }
            put(f.expected_x.data(), bytes, i, a); put(f.expected_y.data(), bytes, i, b);
        }
        baseline(scalar, mode, f, log); baseline(avx2, mode, f, log);
        if (mode == 0) tower_ifft_pair(f.x.data, f.y.data, bytes, log, selected);
        else if (mode == 1) tower_fft_pair(f.x.data, f.y.data, bytes, log, selected);
        else if (mode == 2) tower_fft_out(f.x.data, f.y.data, f.u.data, f.v.data, bytes, log, selected);
        else tower_ifft_accumulate(f.x.data, f.y.data, f.u.data, f.v.data, bytes, log, selected);
        tower_convert_involution(mode < 2 ? f.x.data : f.u.data, f.restored_x.data, bytes, &low_map);
        tower_convert_involution(mode < 2 ? f.y.data : f.v.data, f.restored_y.data, bytes, &low_map);
        require(!std::memcmp(f.restored_x.data, f.expected_x.data(), bytes) &&
                !std::memcmp(f.restored_y.data, f.expected_y.data(), bytes), "tower full butterfly result");
        if (mode >= 2) {
            for (size_t i = 0; i < bytes/2; ++i) {
                require(symbol(f.x.data, bytes, i) == to_tower(symbol(f.cx.data(), bytes, i)), "out x input unchanged");
                require(symbol(f.y.data, bytes, i) == to_tower(symbol(f.cy.data(), bytes, i)), "out y input unchanged");
            }
        }
        if (mode == 3) {
            tower_ifft_accumulate(f.x.data, f.y.data, f.u.data, f.v.data, bytes, log, selected);
            tower_convert_involution(f.u.data, f.u.data, bytes, &low_map);
            tower_convert_involution(f.v.data, f.v.data, bytes, &low_map);
            require(!std::memcmp(f.u.data, f.cu.data(), bytes) && !std::memcmp(f.v.data, f.cv.data(), bytes), "accumulating twice restores outputs");
        }
        f.guards();
    }
    return 4;
}
} // namespace

int main(int argc, char**)
{
    if (argc != 1) { std::fprintf(stderr, "No options or timing mode supported\n"); return 2; }
    setup();
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    const auto* avx2 = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    require(scalar && avx2, "original scalar/AVX2 backends required");
    uint64_t exhaustive = 0, boundary = 0;
    Fixture one(64, 1);
    for (unsigned log = 0; log < 65536; ++log) exhaustive += check_case(one, uint16_t(log), true, *scalar, *avx2);
    for (size_t bytes : {size_t(2), size_t(30), size_t(32), size_t(62), size_t(66), size_t(126), size_t(128),
                         size_t(130), size_t(8190), size_t(8192), size_t(32768), size_t(32770), size_t(65536), size_t(65598)})
        for (size_t offset : {size_t(0), size_t(1), size_t(17)}) {
            Fixture f(bytes, offset);
            for (unsigned log : {0U, 1U, 255U, 256U, 32768U, 65534U, 65535U}) boundary += check_case(f, uint16_t(log), false, *scalar, *avx2);
        }
    // Every compact even tail, including zero, after two full blocks.
    for (size_t tail = 0; tail < 64; tail += 2) {
        Fixture f(128+tail, 17);
        for (unsigned log : {0U, 12345U, 65535U}) boundary += check_case(f, uint16_t(log), false, *scalar, *avx2);
    }
    // All symbols through the new two-shuffle conversion, disjoint and in-place.
    Fixture conversion(128, 1);
    for (unsigned base = 0; base < 65536; base += 64) {
        for (unsigned i = 0; i < 64; ++i) put(conversion.cx.data(), 128, i, uint16_t(base+i));
        tower_convert_involution(conversion.cx.data(), conversion.x.data, 128, &low_map);
        for (unsigned i = 0; i < 64; ++i) require(symbol(conversion.x.data, 128, i) == to_tower(uint16_t(base+i)), "new conversion all symbols");
        tower_convert_involution(conversion.x.data, conversion.x.data, 128, &low_map);
        require(!std::memcmp(conversion.x.data, conversion.cx.data(), 128), "new involution all symbols");
        conversion.guards();
    }
    for (uint16_t log : {uint16_t(0), uint16_t(65535)}) {
        tower_ifft_pair(nullptr, nullptr, 0, log, nullptr);
        tower_fft_pair(nullptr, nullptr, 0, log, nullptr);
        tower_fft_out(nullptr, nullptr, nullptr, nullptr, 0, log, nullptr);
        tower_ifft_accumulate(nullptr, nullptr, nullptr, nullptr, 0, log, nullptr);
    }
    tower_convert_involution(nullptr, nullptr, 0, nullptr);
    std::printf("{\"tracker\":\"leopard-79h.18.20.3\",\"exhaustive_butterfly_cases\":%llu,"
                "\"boundary_butterfly_cases\":%llu,\"input_pair_basis_symbols_per_log\":32,"
                "\"logs_including_zero_skew\":65536,\"conversion_symbols\":65536,"
                "\"high_identity_values\":256,\"zero_byte_null_checks\":9,\"timed\":false,"
                "\"public_codec_qualified\":false}\n", (unsigned long long)exhaustive, (unsigned long long)boundary);
}
