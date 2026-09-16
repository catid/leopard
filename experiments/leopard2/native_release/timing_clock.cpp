// Identical driver object links to steady, synthetic-witness, or abort clocks.
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

namespace {
unsigned calls = 0;
uint64_t witnesses[18] = {};
}
#ifdef LEO_NATIVE_CLOCK_SYNTHETIC
extern "C" uint64_t NativeWitnessCalls();
#endif

extern "C" const char* NativeClockKind()
{
#if defined(LEO_NATIVE_CLOCK_SYNTHETIC)
    return "synthetic";
#elif defined(LEO_NATIVE_CLOCK_ABORT)
    return "abort";
#else
    return "steady";
#endif
}
extern "C" unsigned NativeClockCalls() { return calls; }
extern "C" uint64_t NativeClockWitnessAt(unsigned index)
{
    if (index >= 18) throw std::runtime_error("clock trace overflow");
    return witnesses[index];
}
extern "C" int64_t NativeNow()
{
#if defined(LEO_NATIVE_CLOCK_ABORT)
    std::fputs("unexpected timing clock in clock-free probe\n", stderr);
    std::_Exit(86);
#else
    if (calls >= 18) throw std::runtime_error("too many clock calls");
    const unsigned index = calls++;
#if defined(LEO_NATIVE_CLOCK_SYNTHETIC)
    witnesses[index] = NativeWitnessCalls();
    const int64_t begin = INT64_C(1000000) + (index / 2) * INT64_C(1000000000);
    const char* fault = std::getenv("LEO_NATIVE_TEST_CLOCK_FAULT");
    if (fault) {
        if (std::strcmp(fault, "equal") == 0) return begin;
        if (std::strcmp(fault, "reverse") == 0) return begin - (index % 2);
        if (std::strcmp(fault, "negative") == 0) return -1;
        if (std::strcmp(fault, "huge") == 0) return index % 2 ? INT64_MAX : 0;
        throw std::runtime_error("unknown clock fault");
    }
    return begin + (index % 2 ? 31000000 + (index / 2) * 100 : 0);
#else
    (void)index;
    return std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
#endif
#endif
}
