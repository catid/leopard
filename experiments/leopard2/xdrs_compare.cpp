/*
    Differential harness for the public XDRS API.  It intentionally contains
    no XDRS implementation code; provide P_function.cpp and P_function.h from
    the separately cloned Apache-2.0 research repository when compiling it.
*/

#include "leopard2.h"
#include "P_function.h"

#include <stdint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

typedef std::vector<std::vector<uint8_t> > Shards;

uint8_t cantor_to_polynomial(uint8_t value)
{
    static const uint8_t basis[8] = { 1, 214, 152, 146, 86, 200, 88, 230 };
    uint8_t result = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (value & (1u << bit))
            result ^= basis[bit];
    return result;
}

uint8_t polynomial_to_cantor(uint8_t value)
{
    for (unsigned candidate = 0; candidate < 256; ++candidate)
        if (cantor_to_polynomial(static_cast<uint8_t>(candidate)) == value)
            return static_cast<uint8_t>(candidate);
    throw std::runtime_error("Cantor map is not invertible");
}

class Aligned
{
public:
    explicit Aligned(size_t bytes) : p(NULL), n(bytes)
    {
        if (n && posix_memalign(&p, leo2_scratch_alignment(), n))
            p = NULL;
    }
    ~Aligned() { free(p); }
    void* p;
    size_t n;
};

void check(bool ok, const char* what)
{
    if (!ok)
        throw std::runtime_error(what);
}

uint64_t next(uint64_t* state)
{
    uint64_t z = (*state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

struct XdrsState
{
    std::vector<GFSymbol> index;
    std::vector<GFSymbol> alpha;
    std::vector<GFSymbol> skew;
    std::vector<GFSymbol> b;
    std::vector<GFSymbol> log_walsh;
    std::vector<GFSymbol> base;
    std::vector<GFSymbol> coef_l;
    std::vector<GFSymbol> coef_h;
    std::vector<int> s;

    explicit XdrsState(unsigned k)
        : index(Size), alpha(Size), skew(Size), b(Size), log_walsh(Size)
        , base(len - 1), coef_l(Size / k > 1 ? Size / k - 1 : 1)
        , coef_h(Size), s(len)
    {
        function::main_params(64, 1, k);
        function::array_params(&index[0], &alpha[0], &skew[0], &b[0],
            &log_walsh[0], &base[0], &coef_l[0], &coef_h[0], &s[0]);
        function::init();
        function::init_dec();
    }
};

void run_case(leo2_context* context, unsigned k)
{
    const bool low = k <= 128;
    const unsigned r = Size - k;
    XdrsState state(k);
    (void)state;

    Shards original(k, std::vector<uint8_t>(64));
    uint64_t random = UINT64_C(0x5844525300000000) ^ k;
    for (unsigned i = 0; i < k; ++i)
        for (unsigned j = 0; j < 64; ++j)
            original[i][j] = static_cast<uint8_t>(next(&random) >> 56);

    Shards xdrs_original_storage(k, std::vector<uint8_t>(64));
    for (unsigned i = 0; i < k; ++i)
        for (unsigned j = 0; j < 64; ++j)
            xdrs_original_storage[i][j] = cantor_to_polynomial(original[i][j]);

    const unsigned xdrs_parity_slots = low ? r : r * 2;
    Shards xdrs_parity(xdrs_parity_slots, std::vector<uint8_t>(64));
    std::vector<const GFSymbol*> xdrs_original(k);
    std::vector<GFSymbol*> xdrs_output(xdrs_parity_slots);
    std::vector<const GFSymbol*> xdrs_parity_input(xdrs_parity_slots);
    for (unsigned i = 0; i < k; ++i)
        xdrs_original[i] = &xdrs_original_storage[i][0];
    for (unsigned i = 0; i < xdrs_parity_slots; ++i)
    {
        xdrs_output[i] = &xdrs_parity[i][0];
        xdrs_parity_input[i] = &xdrs_parity[i][0];
    }
    if (low)
        function::ReedSolomonEncodeL(&xdrs_original[0], &xdrs_output[0]);
    else
        function::ReedSolomonEncodeH(&xdrs_original[0], &xdrs_output[0]);

    leo2_codec* codec = NULL;
    const leo2_profile profile = low
        ? LEO2_PROFILE_LOW_V1 : LEO2_PROFILE_LEGACY_HIGH_V1;
    check(leo2_codec_create(context, k, r, profile, LEO2_FIELD_GF8,
        NULL, &codec) == LEO2_SUCCESS, "Leopard2 codec create");
    Shards parity(r, std::vector<uint8_t>(64));
    std::vector<const void*> input(k);
    std::vector<void*> output(r);
    for (unsigned i = 0; i < k; ++i)
        input[i] = &original[i][0];
    for (unsigned i = 0; i < r; ++i)
        output[i] = &parity[i][0];
    size_t scratch_bytes = 0;
    check(leo2_encode_scratch_size(codec, 64, &scratch_bytes) == LEO2_SUCCESS,
        "Leopard2 encode scratch");
    Aligned encode_scratch(scratch_bytes);
    check(encode_scratch.n == 0 || encode_scratch.p,
        "Leopard2 encode scratch allocation");
    check(leo2_encode(codec, 64, &input[0], &output[0],
        encode_scratch.p, encode_scratch.n) == LEO2_SUCCESS,
        "Leopard2 encode");
    Shards xdrs_parity_cantor(r, std::vector<uint8_t>(64));
    for (unsigned i = 0; i < r; ++i)
        for (unsigned j = 0; j < 64; ++j)
            xdrs_parity_cantor[i][j] = polynomial_to_cantor(xdrs_parity[i][j]);
    unsigned parity_mismatches = 0;
    for (unsigned i = 0; i < r; ++i)
        if (parity[i] != xdrs_parity_cantor[i])
        {
            unsigned byte = 0;
            while (byte < 64 && parity[i][byte] == xdrs_parity_cantor[i][byte])
                ++byte;
            if (parity_mismatches == 0)
                std::cerr << "first parity mismatch K=" << k << " row=" << i
                          << " byte=" << byte << " leopard="
                          << static_cast<unsigned>(parity[i][byte]) << " xdrs="
                          << static_cast<unsigned>(xdrs_parity_cantor[i][byte]) << std::endl;
            ++parity_mismatches;
        }

    const unsigned losses = k < 4 ? k : 4;
    std::vector<bool> erased_vector(Size, false);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> parity_present(r, 0);
    std::vector<unsigned> missing(losses);
    for (unsigned i = 0; i < losses; ++i)
    {
        missing[i] = losses == 1 ? k / 2 : i * (k - 1) / (losses - 1);
        original_present[missing[i]] = 0;
        const unsigned coordinate = low ? missing[i] : r + missing[i];
        erased_vector[coordinate] = true;
        parity_present[i] = 1;
    }
    for (unsigned i = losses; i < r; ++i)
    {
        const unsigned coordinate = low ? k + i : i;
        erased_vector[coordinate] = true;
    }
    bool* erased = new bool[Size];
    for (unsigned i = 0; i < Size; ++i)
        erased[i] = erased_vector[i];
    Shards xdrs_codeword(Size, std::vector<uint8_t>(64));
    std::vector<GFSymbol*> xdrs_codeword_ptr(Size);
    for (unsigned i = 0; i < Size; ++i)
        xdrs_codeword_ptr[i] = &xdrs_codeword[i][0];
    std::vector<GFSymbol> locator(Size);
    if (low)
        function::Algorithm2::ReedSolomondecodeL(&xdrs_original[0],
            &xdrs_parity_input[0], &xdrs_codeword_ptr[0], erased, &locator[0]);
    else
        function::Algorithm3::ReedSolomondecodeH(&xdrs_original[0],
            &xdrs_parity_input[0], &xdrs_codeword_ptr[0], erased, &locator[0]);
    for (unsigned i = 0; i < losses; ++i)
    {
        const unsigned coordinate = low ? missing[i] : r + missing[i];
        check(xdrs_codeword[coordinate] == xdrs_original_storage[missing[i]],
            "XDRS specialized decode mismatch");
    }
    delete[] erased;

    leo2_decode_plan* plan = NULL;
    check(leo2_decode_plan_create(codec, &original_present[0],
        &parity_present[0], &plan) == LEO2_SUCCESS,
        "Leopard2 plan create");
    std::vector<const void*> received_original(k, NULL);
    std::vector<const void*> received_parity(r, NULL);
    Shards restored(k, std::vector<uint8_t>(64));
    std::vector<void*> restored_output(k, NULL);
    for (unsigned i = 0; i < k; ++i)
        if (original_present[i])
            received_original[i] = &original[i][0];
        else
            restored_output[i] = &restored[i][0];
    for (unsigned i = 0; i < r; ++i)
        if (parity_present[i])
            received_parity[i] = &parity[i][0];
    check(leo2_decode_plan_scratch_size(plan, 64, &scratch_bytes) ==
        LEO2_SUCCESS, "Leopard2 decode scratch");
    Aligned decode_scratch(scratch_bytes);
    check(decode_scratch.n == 0 || decode_scratch.p,
        "Leopard2 decode scratch allocation");
    check(leo2_decode_plan_execute(plan, 64, &received_original[0],
        &received_parity[0], &restored_output[0], decode_scratch.p,
        decode_scratch.n) == LEO2_SUCCESS, "Leopard2 decode");
    for (unsigned i = 0; i < losses; ++i)
        check(restored[missing[i]] == original[missing[i]],
            "Leopard2 specialized decode mismatch");

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    std::cout << "xdrs differential profile=" << (low ? "low-A2" : "high-A3")
              << " K=" << k << " R=" << r << " parity_bytes=" << r * 64
              << " parity_mismatches=" << parity_mismatches
              << " recovered=" << losses << std::endl;
}

} // namespace

int main()
{
    try
    {
        leo2_context* context = NULL;
        check(leo2_context_create(NULL, &context) == LEO2_SUCCESS,
            "Leopard2 context create");
        const unsigned cases[] = { 8, 16, 32, 64, 128, 192, 224, 240, 248 };
        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i)
            run_case(context, cases[i]);
        leo2_context_destroy(context);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "xdrs differential failed: " << error.what() << std::endl;
        return 1;
    }
}
