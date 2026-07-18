#include "../leopard.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <vector>

namespace {

static bool ExpectUnsigned(
    const char* label,
    unsigned actual,
    unsigned expected)
{
    if (actual == expected)
        return true;

    fprintf(stderr, "%s: got %u, expected %u\n", label, actual, expected);
    return false;
}

static bool ExpectResult(
    const char* label,
    LeopardResult actual,
    LeopardResult expected)
{
    if (actual == expected)
        return true;

    fprintf(stderr,
        "%s: got %d (%s), expected %d (%s)\n",
        label,
        static_cast<int>(actual),
        leo_result_string(actual),
        static_cast<int>(expected),
        leo_result_string(expected));
    return false;
}

struct WorkCountCase
{
    unsigned OriginalCount;
    unsigned RecoveryCount;
    unsigned EncodeWorkCount;
    unsigned DecodeWorkCount;
    const char* Label;
};

static bool TestWorkCounts()
{
    static const WorkCountCase kValidCases[] = {
        { 1, 1, 1, 1, "k=1 r=1" },
        { 8, 1, 1, 8, "m=1" },
        { 2, 2, 4, 4, "small power of two" },
        { 32767, 32767, 65536, 65536, "below count boundary" },
        { 32768, 32768, 65536, 65536, "exact count boundary" },
        // The public count constraint is satisfied here, although padding the
        // transform requires 131072 decoder work buffers and the operation
        // subsequently reports Leopard_TooMuchData.
        { 32769, 32767, 65536, 131072, "padded transform boundary" },
        { 65535, 1, 1, 65535, "largest m=1 count boundary" },
    };

    bool ok = true;
    for (unsigned i = 0;
         i < sizeof(kValidCases) / sizeof(kValidCases[0]);
         ++i)
    {
        const WorkCountCase& test = kValidCases[i];
        ok = ExpectUnsigned(
                 test.Label,
                 leo_encode_work_count(
                     test.OriginalCount, test.RecoveryCount),
                 test.EncodeWorkCount) && ok;
        ok = ExpectUnsigned(
                 test.Label,
                 leo_decode_work_count(
                     test.OriginalCount, test.RecoveryCount),
                 test.DecodeWorkCount) && ok;
    }

    static const unsigned kInvalidCases[][2] = {
        { 0, 0 },
        { 0, 1 },
        { 1, 0 },
        { 1, 2 },
        { 65535, 2 },
        { UINT_MAX, 1 },
        { UINT_MAX, UINT_MAX },
        // Both of these sums wrap in unsigned arithmetic.
        { UINT_MAX - 1U, 2 },
        { UINT_MAX / 2U + 1U, UINT_MAX / 2U + 1U },
    };

    for (unsigned i = 0;
         i < sizeof(kInvalidCases) / sizeof(kInvalidCases[0]);
         ++i)
    {
        const unsigned original_count = kInvalidCases[i][0];
        const unsigned recovery_count = kInvalidCases[i][1];
        char label[96];
        snprintf(label, sizeof(label),
            "invalid encode work count k=%u r=%u",
            original_count, recovery_count);
        ok = ExpectUnsigned(label,
                 leo_encode_work_count(original_count, recovery_count),
                 0) && ok;
        snprintf(label, sizeof(label),
            "invalid decode work count k=%u r=%u",
            original_count, recovery_count);
        ok = ExpectUnsigned(label,
                 leo_decode_work_count(original_count, recovery_count),
                 0) && ok;
    }

    return ok;
}

struct ApiCountCase
{
    unsigned OriginalCount;
    unsigned RecoveryCount;
    LeopardResult Expected;
    const char* Label;
};

static bool TestApiCountValidation()
{
    uint8_t data[64] = {};
    uint8_t output[64] = {};
    const void* const bounded_original[1] = { data };
    const void* const bounded_recovery[1] = { data };
    void* bounded_work[1] = { output };

    static const ApiCountCase kCases[] = {
        { 0, 0, Leopard_InvalidCounts, "zero counts" },
        { 4, 0, Leopard_InvalidCounts, "zero recovery count" },
        { 4, 5, Leopard_InvalidCounts, "recovery exceeds original" },
        { 1, UINT_MAX, Leopard_InvalidCounts,
          "malformed relationship wins over capacity" },
        { 65536, 1, Leopard_TooMuchData, "one above count boundary" },
        { UINT_MAX, 1, Leopard_TooMuchData, "UINT_MAX m=1 count" },
        { UINT_MAX, UINT_MAX, Leopard_TooMuchData,
          "UINT_MAX equal counts" },
        { UINT_MAX - 1U, 2, Leopard_TooMuchData,
          "wrapped total" },
        { UINT_MAX / 2U + 1U, UINT_MAX / 2U + 1U,
          Leopard_TooMuchData, "wrapped equal-count total" },
    };

    bool ok = true;
    for (unsigned i = 0; i < sizeof(kCases) / sizeof(kCases[0]); ++i)
    {
        const ApiCountCase& test = kCases[i];
        char label[128];
        snprintf(label, sizeof(label), "encode %s", test.Label);
        ok = ExpectResult(label,
                 leo_encode(
                     64,
                     test.OriginalCount,
                     test.RecoveryCount,
                     0,
                     bounded_original,
                     bounded_work),
                 test.Expected) && ok;

        snprintf(label, sizeof(label), "decode %s", test.Label);
        ok = ExpectResult(label,
                 leo_decode(
                     64,
                     test.OriginalCount,
                     test.RecoveryCount,
                     0,
                     bounded_original,
                     bounded_recovery,
                     bounded_work),
                 test.Expected) && ok;
    }

    return ok;
}

static bool TestSpecialCases()
{
    uint8_t originals[4][64];
    uint8_t parity[64];
    uint8_t decoded[4][64];
    const void* original_ptrs[4];
    void* encode_work[1] = { parity };
    void* decode_work[4] = {
        decoded[0], decoded[1], decoded[2], decoded[3]
    };

    for (unsigned i = 0; i < 4; ++i)
    {
        original_ptrs[i] = originals[i];
        for (unsigned j = 0; j < 64; ++j)
            originals[i][j] = static_cast<uint8_t>(i * 61U + j * 17U + 3U);
    }

    bool ok = ExpectResult(
        "m=1 encode",
        leo_encode(64, 4, 1, 1, original_ptrs, encode_work),
        Leopard_Success);
    ok = ExpectResult(
             "general encode rejects zero work count",
             leo_encode(64, 4, 2, 0, original_ptrs, encode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "m=1 encode rejects zero work count",
             leo_encode(64, 4, 1, 0, original_ptrs, encode_work),
             Leopard_InvalidCounts) && ok;
    for (unsigned j = 0; j < 64; ++j)
    {
        const uint8_t expected = static_cast<uint8_t>(
            originals[0][j] ^ originals[1][j] ^
            originals[2][j] ^ originals[3][j]);
        if (parity[j] != expected)
        {
            fprintf(stderr,
                "m=1 encode byte %u: got %u, expected %u\n",
                j, parity[j], expected);
            ok = false;
            break;
        }
    }

    const void* received_originals[4] = {
        originals[0], originals[1], NULL, originals[3]
    };
    const void* recovery_ptrs[1] = { parity };
    const void* missing_recovery_ptrs[1] = { NULL };
    const void* missing_general_recovery[2] = { NULL, NULL };
    memset(decoded, 0, sizeof(decoded));
    ok = ExpectResult(
             "general no-loss decode rejects zero work count",
             leo_decode(
                 64, 4, 2, 0,
                 original_ptrs, missing_general_recovery, decode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "need-more-data precedes work-count validation",
             leo_decode(
                 64, 4, 1, 0,
                 received_originals, missing_recovery_ptrs, decode_work),
             Leopard_NeedMoreData) && ok;
    ok = ExpectResult(
             "m=1 decode rejects zero work count",
             leo_decode(
                 64, 4, 1, 0,
                 received_originals, recovery_ptrs, decode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "m=1 decode",
             leo_decode(
                 64, 4, 1, 4,
                 received_originals, recovery_ptrs, decode_work),
             Leopard_Success) && ok;
    if (memcmp(decoded[2], originals[2], 64) != 0)
    {
        fprintf(stderr, "m=1 decode recovered the wrong bytes\n");
        ok = false;
    }

    uint8_t singleton_recovery[64] = {};
    uint8_t singleton_decoded[64] = {};
    const void* singleton_original[1] = { originals[0] };
    const void* missing_singleton[1] = { NULL };
    const void* singleton_recovery_ptr[1] = { singleton_recovery };
    void* singleton_encode_work[1] = { singleton_recovery };
    void* singleton_decode_work[1] = { singleton_decoded };
    ok = ExpectResult(
             "k=1 encode rejects zero work count",
             leo_encode(
                 64, 1, 1, 0,
                 singleton_original, singleton_encode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "k=1 encode",
             leo_encode(
                 64, 1, 1, 1,
                 singleton_original, singleton_encode_work),
             Leopard_Success) && ok;
    ok = ExpectResult(
             "k=1 decode rejects zero work count",
             leo_decode(
                 64, 1, 1, 0,
                 missing_singleton,
                 singleton_recovery_ptr,
                 singleton_decode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "k=1 decode",
             leo_decode(
                 64, 1, 1, 1,
                 missing_singleton,
                 singleton_recovery_ptr,
                 singleton_decode_work),
             Leopard_Success) && ok;
    if (memcmp(singleton_decoded, originals[0], 64) != 0)
    {
        fprintf(stderr, "k=1 round trip recovered the wrong bytes\n");
        ok = false;
    }

    // A no-loss decode must use the original and must not require or prefer a
    // recovery pointer.  Historically the packed k=1 shortcut ran first and
    // dereferenced recovery_data[0] even when it was legitimately NULL.
    const void* missing_singleton_recovery[1] = { NULL };
    memset(singleton_decoded, 0, sizeof(singleton_decoded));
    ok = ExpectResult(
             "k=1 no-loss decode rejects zero work count",
             leo_decode(
                 64, 1, 1, 0,
                 singleton_original,
                 missing_singleton_recovery,
                 singleton_decode_work),
             Leopard_InvalidCounts) && ok;
    ok = ExpectResult(
             "k=1 no-loss decode without recovery",
             leo_decode(
                 64, 1, 1, 1,
                 singleton_original,
                 missing_singleton_recovery,
                 singleton_decode_work),
             Leopard_Success) && ok;
    if (memcmp(singleton_decoded, originals[0], 64) != 0)
    {
        fprintf(stderr, "k=1 no-loss decode did not copy the original\n");
        ok = false;
    }

    singleton_recovery[0] ^= UINT8_C(0xff);
    memset(singleton_decoded, 0, sizeof(singleton_decoded));
    ok = ExpectResult(
             "k=1 no-loss decode ignores stale recovery",
             leo_decode(
                 64, 1, 1, 1,
                 singleton_original,
                 singleton_recovery_ptr,
                 singleton_decode_work),
             Leopard_Success) && ok;
    if (memcmp(singleton_decoded, originals[0], 64) != 0)
    {
        fprintf(stderr, "k=1 no-loss decode preferred stale recovery\n");
        ok = false;
    }

    return ok;
}

static bool TestFF16PackedGolden()
{
    static const unsigned kOriginalCount = 300;
    static const unsigned kRecoveryCount = 100;
    static const unsigned kBufferBytes = 64;
    const unsigned work_count = leo_encode_work_count(
        kOriginalCount, kRecoveryCount);
    if (work_count != 256)
    {
        fprintf(stderr, "FF16 golden work count changed: %u\n", work_count);
        return false;
    }

    std::vector<std::vector<uint8_t> > originals(
        kOriginalCount, std::vector<uint8_t>(kBufferBytes));
    std::vector<std::vector<uint8_t> > work(
        work_count, std::vector<uint8_t>(kBufferBytes));
    std::vector<const void*> original_ptrs(kOriginalCount);
    std::vector<void*> work_ptrs(work_count);
    for (unsigned i = 0; i < kOriginalCount; ++i)
    {
        original_ptrs[i] = originals[i].data();
        for (unsigned j = 0; j < kBufferBytes; ++j)
        {
            originals[i][j] = static_cast<uint8_t>(
                i * 37U + j * 19U + (i >> 3) * 11U + 5U);
        }
    }
    for (unsigned i = 0; i < work_count; ++i)
        work_ptrs[i] = work[i].data();

    if (!ExpectResult(
            "FF16 packed golden encode",
            leo_encode(
                kBufferBytes,
                kOriginalCount,
                kRecoveryCount,
                work_count,
                original_ptrs.data(),
                work_ptrs.data()),
            Leopard_Success))
        return false;

    // This historical regression seed is deliberately retained verbatim; it
    // is FNV-1a-like but is not the canonical 64-bit FNV offset basis.
    uint64_t hash = UINT64_C(1469598103934665603);
    for (unsigned i = 0; i < kRecoveryCount; ++i)
    {
        for (unsigned j = 0; j < kBufferBytes; ++j)
        {
            hash ^= work[i][j];
            hash *= UINT64_C(1099511628211);
        }
    }
    if (hash != UINT64_C(0xe74c4664e37dac03))
    {
        fprintf(stderr,
            "FF16 packed golden hash changed: 0x%016llx\n",
            static_cast<unsigned long long>(hash));
        return false;
    }

    // Exercise the FF16 decoder's locator, IFFT, derivative, truncated FFT,
    // and final scaling paths with both original and recovery erasures.
    const unsigned decode_work_count = leo_decode_work_count(
        kOriginalCount, kRecoveryCount);
    if (decode_work_count != 512)
    {
        fprintf(stderr,
            "FF16 decode work count changed: %u\n", decode_work_count);
        return false;
    }
    std::vector<const void*> received_originals(original_ptrs);
    static const unsigned kLostOriginals[] = { 0, 17, 299 };
    for (unsigned i = 0;
         i < sizeof(kLostOriginals) / sizeof(kLostOriginals[0]); ++i)
    {
        received_originals[kLostOriginals[i]] = NULL;
    }
    std::vector<const void*> recovery_ptrs(kRecoveryCount);
    for (unsigned i = 0; i < kRecoveryCount; ++i)
        recovery_ptrs[i] = work[i].data();
    recovery_ptrs[3] = NULL;
    recovery_ptrs[57] = NULL;

    std::vector<std::vector<uint8_t> > decode_work(
        decode_work_count, std::vector<uint8_t>(kBufferBytes));
    std::vector<void*> decode_work_ptrs(decode_work_count);
    for (unsigned i = 0; i < decode_work_count; ++i)
        decode_work_ptrs[i] = decode_work[i].data();

    if (!ExpectResult(
            "FF16 packed decode",
            leo_decode(
                kBufferBytes,
                kOriginalCount,
                kRecoveryCount,
                decode_work_count,
                received_originals.data(),
                recovery_ptrs.data(),
                decode_work_ptrs.data()),
            Leopard_Success))
        return false;
    for (unsigned i = 0;
         i < sizeof(kLostOriginals) / sizeof(kLostOriginals[0]); ++i)
    {
        const unsigned index = kLostOriginals[i];
        if (memcmp(
                decode_work[index].data(),
                originals[index].data(),
                kBufferBytes) != 0)
        {
            fprintf(stderr,
                "FF16 packed decode recovered original %u incorrectly\n",
                index);
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    bool ok = TestWorkCounts();

    const LeopardResult init_result =
        static_cast<LeopardResult>(leo_init());
    if (!ExpectResult("leo_init", init_result, Leopard_Success))
        return 1;

    ok = TestApiCountValidation() && ok;
    ok = TestSpecialCases() && ok;
    ok = TestFF16PackedGolden() && ok;

    if (!ok)
        return 1;

    printf("Packed API count validation tests passed\n");
    return 0;
}
