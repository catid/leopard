// Exercise the actual generated synthetic wrapper; no codec or real clock.
#define LEO_PAIRED_SYNTHETIC 1
#include "paired_timer_clock.cpp"
extern "C" unsigned LeoPairedWitnessCalls() { return 0; }

int main(int argc, char** argv)
{
    if (argc != 2) return 98;
    if (!std::strcmp(argv[1],"unknown-fault")) setenv("LEO_PAIRED_TEST_CLOCK_FAULT","unknown",1);
    else if (!std::strcmp(argv[1],"epoch2-no-fault")) setenv("LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH","2",1);
    else if (std::strcmp(argv[1],"capacity") && std::strcmp(argv[1],"overflow"))
        setenv("LEO_PAIRED_TEST_CLOCK_FAULT_EPOCH",argv[1],1);
    for (unsigned i=0;i<504;++i) {
        const int64_t actual = __wrap__ZNSt6chrono3_V212steady_clock3nowEv().time_since_epoch().count();
        const int64_t sample = i/2;
        const int64_t expected = INT64_C(1000000)+sample*INT64_C(100000)+(i%2 ? 257+sample*17 : 0);
        if (actual != expected || trace.calls[i] != 0) return 97;
    }
    if (!std::strcmp(argv[1],"overflow")) { __wrap__ZNSt6chrono3_V212steady_clock3nowEv(); return 99; }
    std::puts("{\"schema\":\"paired-epoch-clock-unit/v1\",\"endpoints\":504,\"timed\":false}");
    return 0;
}
