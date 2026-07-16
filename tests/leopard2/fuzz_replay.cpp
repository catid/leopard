/* Deterministic smoke driver for the libFuzzer entry point. */

#include <stddef.h>
#include <stdint.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

extern "C" int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size);

static uint64_t parse(const char* text)
{
    if (!text || !*text)
        return 0;
    errno = 0;
    char* end = NULL;
    const unsigned long long value = strtoull(text, &end, 0);
    return errno == 0 && end && *end == '\0' ? value : 0;
}

int main(int argc, char** argv)
{
    if (argc > 3)
        return 2;
    uint64_t state = argc > 1 ? parse(argv[1]) : UINT64_C(0x4c656f7061726432);
    const uint64_t iterations = argc > 2 ? parse(argv[2]) : 512;
    if (state == 0 || iterations == 0 || iterations > UINT64_C(10000000))
        return 2;
    uint8_t input[96];
    for (uint64_t iteration = 0; iteration < iterations; ++iteration)
    {
        for (size_t i = 0; i < sizeof(input); ++i)
        {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            input[i] = static_cast<uint8_t>(
                state * UINT64_C(2685821657736338717));
        }
        if (LLVMFuzzerTestOneInput(input, sizeof(input)) != 0)
            return 1;
    }
    printf("leopard2_fuzz_replay seed=%llu iterations=%llu passed\n",
        static_cast<unsigned long long>(argc > 1 ? parse(argv[1]) :
            UINT64_C(0x4c656f7061726432)),
        static_cast<unsigned long long>(iterations));
    return 0;
}
