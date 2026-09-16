#include "EncodeTiming.h"
#include <climits>
#include <cstdio>
#include <stdexcept>

static void Require(bool condition)
{
    if (!condition) throw std::runtime_error("group timing unit mismatch");
}

int main()
{
    unsigned tests = 0;
    for (unsigned group : {1U, 2U, 16U, 128U, 256U, 1048576U, 4194304U}) {
        for (bool sample : {false, true}) {
            uint64_t encodes = 0;
            unsigned clocks = 0;
            auto values = native_release::Run(group, sample, [&]() { ++encodes; }, [&]() {
                const unsigned endpoint = clocks++;
                Require(encodes == 4 + static_cast<uint64_t>((endpoint + 1) / 2) * group);
                return INT64_C(1000) + (endpoint / 2) * INT64_C(1000000) +
                    (endpoint % 2 ? 100 + endpoint / 2 : 0);
            });
            Require(values.size() == 9 && clocks == (sample ? 18U : 0U) &&
                encodes == 4 + UINT64_C(9) * group);
            for (unsigned i = 0; i < 9; ++i)
                Require(values[i] == (sample ? 100 + i : 0));
            ++tests;
        }
    }
    for (unsigned group : {0U, 3U, 4194305U, UINT_MAX}) {
        unsigned work = 0;
        bool rejected = false;
        try {
            native_release::Run(group, true, [&]() { ++work; }, [&]() {
                ++work; return INT64_C(0);
            });
        } catch (const std::runtime_error&) { rejected = true; }
        Require(rejected && work == 0);
        ++tests;
    }
    const int64_t bad[][2] = {
        {-1, 1}, {INT64_MIN, INT64_MAX}, {0, 0}, {2, 1},
        {0, INT64_C(9007199254740992)}, {0, INT64_MAX}
    };
    for (const auto& pair : bad) {
        bool rejected = false;
        try { native_release::Duration(pair[0], pair[1]); }
        catch (const std::runtime_error&) { rejected = true; }
        Require(rejected);
        ++tests;
    }
    Require(native_release::Duration(0, INT64_C(9007199254740991)) == UINT64_C(9007199254740991));
    Require(native_release::Duration(INT64_MAX - 1, INT64_MAX) == 1);
    tests += 2;
    std::printf("{\"schema\":\"native-encode-timing-unit/v1\",\"cases\":%u,\"real_clocks\":0}\n", tests);
}
