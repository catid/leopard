/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.
    See leopard.h for the full BSD license notice.
*/

#include "LeopardCommon.h"
#include "leopard.h"

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace {

static int Fail(const char* message, int detail = 0)
{
    std::fprintf(stderr, "%s: %d\n", message, detail);
    return 1;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
        return Fail("usage: test_legacy_simd_init_failure ff8|ff16");

    unsigned failure_index = 0;
    unsigned expected_retry_allocations = 0;
    if (std::strcmp(argv[1], "ff8") == 0)
    {
        failure_index = 1;
        expected_retry_allocations = 2;
    }
    else if (std::strcmp(argv[1], "ff16") == 0)
    {
        failure_index = 2;
        expected_retry_allocations = 1;
    }
    else
        return Fail("unknown field argument");

    leopard::TestSetSIMDSafeAllocationFailure(failure_index);
    const int failed_result = leo_init();
    if (failed_result != Leopard_Platform)
        return Fail("injected allocation returned wrong result", failed_result);
    if (leopard::TestSIMDSafeAllocationAttempts() != failure_index)
        return Fail("allocation failure occurred at wrong attempt",
            static_cast<int>(leopard::TestSIMDSafeAllocationAttempts()));

    uint8_t original_storage[64];
    uint8_t recovery_storage[64];
    for (unsigned i = 0; i < sizeof(original_storage); ++i)
        original_storage[i] = static_cast<uint8_t>(i * 17u + 3u);
    const void* original[1] = { original_storage };
    void* work[1] = { recovery_storage };
    if (leo_encode(sizeof(original_storage), 1, 1, 1, original, work) !=
        Leopard_CallInitialize)
        return Fail("failed setup partially published public readiness");

    leopard::TestSetSIMDSafeAllocationFailure(0);
    const int retry_result = leo_init();
    if (retry_result != Leopard_Success)
        return Fail("initialization retry failed", retry_result);
    if (leopard::TestSIMDSafeAllocationAttempts() !=
        expected_retry_allocations)
        return Fail("retry initialized the wrong field-table set",
            static_cast<int>(leopard::TestSIMDSafeAllocationAttempts()));
    if (leo_init() != Leopard_Success)
        return Fail("successful retry was not idempotent");
    if (leopard::TestSIMDSafeAllocationAttempts() !=
        expected_retry_allocations)
        return Fail("idempotent initialization allocated again");

    if (leo_encode(sizeof(original_storage), 1, 1, 1, original, work) !=
            Leopard_Success ||
        std::memcmp(original_storage, recovery_storage,
            sizeof(original_storage)) != 0)
        return Fail("retried initialization did not publish usable state");

    return 0;
}
