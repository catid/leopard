// Untimed integration-cache/dispatch oracle, leopard-79h.38.5.4.18.4.
#include "tower_encoder.h"
#include "LeopardFF16.h"
#include "leopard.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>

namespace {
void require(bool ok, const char* message) {
    if (!ok) { std::fprintf(stderr, "%s\n", message); std::exit(1); }
}
// Independent polynomial basis from the previously qualified algebra oracle.
const uint16_t basis[] = {0x0001,0xACCA,0x3C0E,0x163E,0xC582,0xED2E,0x914C,0x4012,
                         0x6C98,0x10D8,0x6A72,0xB900,0xFDB8,0xFB34,0xFF38,0x991E};
uint16_t polynomial[65536], canonical[65536];
uint16_t multiply(uint16_t a, uint16_t b) {
    uint32_t product = 0;
    for (unsigned bit = 0; bit < 16; ++bit)
        if ((polynomial[b] >> bit) & 1) product ^= uint32_t(polynomial[a]) << bit;
    for (int bit = 30; bit >= 16; --bit)
        if ((product >> bit) & 1) product ^= uint32_t(0x1002D) << (bit-16);
    return canonical[product];
}
void put(uint8_t* p, unsigned lane, uint16_t x) { p[lane] = uint8_t(x); p[32+lane] = uint8_t(x >> 8); }
uint16_t get(const uint8_t* p, unsigned lane) { return p[lane] | (uint16_t(p[32+lane]) << 8); }

unsigned selectors() {
    unsigned cases = 0;
    leopard::backend::Ops ops{};
    leopard2_internal::SparseForwardPlanBatchView empty{nullptr,0,0}, sparse{nullptr,0,1};
    for (bool enabled : {false,true}) {
        tower_encoder::SetEnabled(enabled);
        for (leo2_backend backend : {LEO2_BACKEND_AUTO,LEO2_BACKEND_SCALAR,LEO2_BACKEND_SSSE3,
                                     LEO2_BACKEND_AVX2,LEO2_BACKEND_AVX512,LEO2_BACKEND_GFNI}) {
            ops.kind = backend;
            for (uint64_t bytes : {0U,64U,128U,192U,255U,256U,257U,320U,16384U,16448U})
                for (uint64_t policy : {0U,64U,16384U,16385U,32768U})
                    for (unsigned side : {1U,128U,255U,256U,512U})
                        for (const auto* view : {static_cast<const leopard2_internal::SparseForwardPlanBatchView*>(nullptr),
                                                static_cast<const leopard2_internal::SparseForwardPlanBatchView*>(&empty),
                                                static_cast<const leopard2_internal::SparseForwardPlanBatchView*>(&sparse)}) {
                            const bool expected = enabled && backend == LEO2_BACKEND_AVX2 &&
                                (bytes == 256 || bytes == 320 || bytes == 16384 || bytes == 16448) &&
                                (policy == 16385 || policy == 32768) && (side == 256 || side == 512) && view != &sparse;
                            require(tower_encoder::Select(ops,bytes,policy,side,view) == expected, "selector boundary");
                            ++cases;
                        }
        }
    }
    require(tower_encoder::InitializationCount() == 0, "predicate initialized cache");
    return cases;
}
}

int main(int argc, char** argv)
{
    require(argc == 1 || (argc == 3 && !std::strcmp(argv[1], "--abort") &&
            std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '8'), "no timing mode");
    const unsigned selector_cases = selectors();
    require(leo_init() == 0, "field initialization");
    const auto* original = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX2);
    require(original != nullptr, "qualified AVX2 required");
    const auto& ops = tower_encoder::GetOps(*original);
    require(tower_encoder::InitializationCount() == 1 && &tower_encoder::GetOps(*original) == &ops,
            "immutable ops publication");
    if (argc == 3) {
        switch (argv[2][0]-'0') {
        case 0: ops.ff16_multiply(nullptr,nullptr,0,0); break;
        case 1: ops.ff16_multiply_add(nullptr,nullptr,0,0); break;
        case 2: ops.ff16_fft_butterfly2_out(nullptr,nullptr,nullptr,nullptr,0,0); break;
        case 3: ops.ff16_ifft_butterfly4(nullptr,nullptr,nullptr,nullptr,0,0,0,0); break;
        case 4: ops.ff16_fft_butterfly4(nullptr,nullptr,nullptr,nullptr,0,0,0,0); break;
        case 5: ops.ff16_ifft_butterfly4_out(nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,0,0,0,0); break;
        case 6: ops.ff16_fft_butterfly4_out(nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,0,0,0,0); break;
        case 7: ops.ff16_ifft_butterfly4_range(nullptr,0,0,0,0,0,true); break;
        case 8: ops.ff16_fft_butterfly4_range(nullptr,0,0,0,0,0,true); break;
        }
        require(false, "unsupported callback returned");
    }
    for (unsigned value = 0; value < 65536; ++value) {
        uint16_t p = 0;
        for (unsigned bit = 0; bit < 16; ++bit) if ((value >> bit) & 1) p ^= basis[bit];
        polynomial[value] = p; canonical[p] = uint16_t(value);
    }
    alignas(64) uint8_t x[64]{}, y[64]{}, tx[64]{}, ty[64]{};
    for (unsigned lane = 0; lane < 32; ++lane) {
        put(x,lane,lane < 16 ? uint16_t(1U << lane) : 0);
        put(y,lane,lane >= 16 ? uint16_t(1U << (lane-16)) : 0);
    }
    for (unsigned log = 0; log < 65536; ++log) {
        const uint16_t coefficient = log == 65535 ? 0 : leopard::ff16::MultiplyLogElement(1,uint16_t(log));
        tower_encoder::CopySource(ops,tx,x,64); tower_encoder::CopySource(ops,ty,y,64);
        ops.ff16_ifft_butterfly2(tx,ty,uint16_t(log),64);
        void* work[] = {tx,ty}; tower_encoder::Finish(work,2,64);
        for (unsigned lane = 0; lane < 32; ++lane) {
            const uint16_t a = get(x,lane), b = get(y,lane) ^ a;
            const uint16_t product = multiply(b,coefficient);
            require(product == leopard::ff16::MultiplyElements(b,coefficient), "polynomial/scalar oracle");
            require(get(tx,lane) == (a ^ product) && get(ty,lane) == b, "cached log table/basis");
        }
    }
    tower_encoder::CopySource(*original,tx,x,64);
    require(!std::memcmp(tx,x,64), "canonical copy remains canonical");
    std::printf("{\"tracker\":\"leopard-79h.38.5.4.18.4\",\"selector_cases\":%u,"
                "\"log_cases\":65536,\"basis_pairs\":2097152,\"initializations\":%u,\"timed\":false}\n",
                selector_cases,tower_encoder::InitializationCount());
}
