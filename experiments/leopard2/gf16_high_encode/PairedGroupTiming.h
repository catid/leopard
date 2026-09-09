// leopard-79h.38.5.4.19.1.1: exact grouped public-call boundary.
#ifndef LEOPARD_PAIRED_GROUP_TIMING_H
#define LEOPARD_PAIRED_GROUP_TIMING_H
#include <cstdint>
#include <stdexcept>

namespace leopard_paired {
struct Sample { uint64_t elapsed_ns; double ns_per_call; };

inline void CheckGroup(unsigned operations)
{
    if (operations != 1 && operations != 256)
        throw std::runtime_error("invalid grouped operation count");
}

inline Sample Normalize(int64_t begin, int64_t end, unsigned operations)
{
    CheckGroup(operations);
    // Validate before subtracting: malformed clocks must not cause signed UB.
    if (begin < 0 || end <= begin)
        throw std::runtime_error("nonpositive or reversed grouped clock interval");
    const int64_t elapsed = end - begin;
    if (elapsed > INT64_C(9007199254740991))
        throw std::runtime_error("group duration exceeds exact binary64 integer range");
    // 1 and 256 are powers of two. With the bound above both conversions are
    // exact; this is a group-average cost, never a single-call latency sample.
    return Sample{static_cast<uint64_t>(elapsed), static_cast<double>(elapsed) / operations};
}

template<class Encode, class Clock>
Sample RunGroup(unsigned operations, bool sample, Encode&& encode, Clock&& now)
{
    CheckGroup(operations);
    const int64_t begin = sample ? now() : 0;
    for (unsigned n = 0; n < operations; ++n) encode();
    const int64_t end = sample ? now() : 0;
    return sample ? Normalize(begin, end, operations) : Sample{0, 0.};
}
}
#endif
