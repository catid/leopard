// The driver object is identical with real, aborting and synthetic clocks.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#ifndef LEO_PAIRED_SYNTHETIC
#ifdef LEO_PAIRED_ABORT
extern "C" const char* LeoPairedClockKind() { return "abort"; }
#else
extern "C" const char* LeoPairedClockKind() { return "steady"; }
#endif
#else
extern "C" const char* LeoPairedClockKind() { return "synthetic"; }
extern "C" unsigned LeoPairedWitnessCalls();
namespace {
struct Trace {
    unsigned count = 0;
    unsigned calls[168] = {};
    ~Trace() {
        std::printf("{\"schema\":\"paired-synthetic-clock/v1\",\"clock_calls\":%u,\"public_calls_at_clock\":[",count);
        for (unsigned i=0; i<count; ++i) std::printf("%s%u",i ? "," : "",calls[i]);
        std::puts("],\"timed\":false}");
    }
} trace;
}
extern "C" std::chrono::steady_clock::time_point
__wrap__ZNSt6chrono3_V212steady_clock3nowEv()
{
    if (trace.count >= 168) { std::fputs("too many synthetic clocks\n",stderr); std::exit(89); }
    const unsigned index = trace.count++;
    trace.calls[index] = LeoPairedWitnessCalls();
    const unsigned sample = index/2;
    const int64_t start = INT64_C(1000000) + sample*INT64_C(100000);
    int64_t tick = start + (index%2 ? 257 + sample*17 : 0);
    const char* fault = std::getenv("LEO_PAIRED_TEST_CLOCK_FAULT");
    if (fault && *fault) {
        if (!std::strcmp(fault,"equal")) tick = start;
        else if (!std::strcmp(fault,"reverse")) tick = start - (index%2 ? 1 : 0);
        else if (!std::strcmp(fault,"negative")) tick = index%2 ? INT64_MAX : INT64_MIN;
        else if (!std::strcmp(fault,"huge")) tick = index%2 ? INT64_MAX : 0;
        else { std::fputs("unknown synthetic clock fault\n",stderr); std::exit(89); }
    }
    return std::chrono::steady_clock::time_point(std::chrono::nanoseconds(tick));
}
#endif
