// Link with --wrap for direct clock references from the static probe/codec.
// This does not interpose private calls inside dynamically linked libraries.
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>

static void ClockCalled()
{
    std::fputs("unexpected timing clock in clock-free probe\n", stderr);
    std::_Exit(86);
}

extern "C" int __wrap_clock_gettime(clockid_t, struct timespec*) { ClockCalled(); return -1; }
extern "C" int __wrap_gettimeofday(struct timeval*, void*) { ClockCalled(); return -1; }
extern "C" clock_t __wrap_clock() { ClockCalled(); return static_cast<clock_t>(-1); }
