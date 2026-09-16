// leopard-79h.57.12.2: grouped ordinary public encode, no buffer resets in-loop.
#ifndef LEOPARD_NATIVE_RELEASE_ENCODE_TIMING_H
#define LEOPARD_NATIVE_RELEASE_ENCODE_TIMING_H
#include <cstdint>
#include <stdexcept>
#include <vector>

extern "C" const char* NativeClockKind();
extern "C" int64_t NativeNow();
extern "C" unsigned NativeClockCalls();
extern "C" uint64_t NativeClockWitnessAt(unsigned index);

namespace native_release {
static const unsigned kGroups[] = {4194304, 1048576, 256, 256, 256, 16, 16, 128};
static const unsigned kSamples = 9;
static const unsigned kWarmup = 4;

inline uint64_t Duration(int64_t begin, int64_t end)
{
    if (begin < 0 || end <= begin)
        throw std::runtime_error("nonpositive or reversed clock interval");
    const uint64_t elapsed = static_cast<uint64_t>(end - begin);
    if (elapsed > UINT64_C(9007199254740991))
        throw std::runtime_error("clock interval exceeds exact binary64 range");
    return elapsed;
}

inline void CheckCount(unsigned count)
{
    if (!count || count > 4194304 || (count & (count - 1)))
        throw std::runtime_error("invalid power-of-two group size");
}

template<class Encode, class Clock>
uint64_t Group(unsigned count, bool sample, Encode&& encode, Clock&& now)
{
    CheckCount(count);
    const int64_t begin = sample ? now() : 0;
    for (unsigned i = 0; i < count; ++i) encode();
    const int64_t end = sample ? now() : 0;
    return sample ? Duration(begin, end) : 0;
}

template<class Encode, class Clock>
std::vector<uint64_t> Run(unsigned count, bool sample, Encode&& encode, Clock&& now)
{
    CheckCount(count);
    // Reserve and warm up before clocks; reporting, verification and
    // normalization remain outside groups. API result checks, loop and timer
    // bookkeeping overhead remain in the measured interval, without subtraction.
    std::vector<uint64_t> elapsed;
    elapsed.reserve(kSamples);
    for (unsigned i = 0; i < kWarmup; ++i) encode();
    for (unsigned i = 0; i < kSamples; ++i)
        elapsed.push_back(Group(count, sample, encode, now));
    return elapsed;
}
} // namespace native_release
#endif
