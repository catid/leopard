// Pure clock-free tests for the actual timing boundary used by the frontend.
#include "PairedGroupTiming.h"
#include <cstdio>
#include <limits>
#include <string>
#include <utility>

namespace {
void Check(bool ok) { if (!ok) throw std::runtime_error("group timing test failed"); }
template<class Function> void Reject(Function&& function)
{
    bool failed = false;
    try { function(); } catch (const std::runtime_error&) { failed = true; }
    Check(failed);
}
}

int main()
{
    try {
        using namespace leopard_paired;
        unsigned cases = 0;
        for (unsigned group : {1U,256U}) {
            std::string trace;
            unsigned clocks = 0;
            const auto clock = [&]() -> int64_t {
                trace += clocks++ ? 'Z' : 'A';
                return clocks == 1 ? 1000 : 1257;
            };
            const auto encode = [&]() { trace += 'E'; };
            const Sample s = RunGroup(group,true,encode,clock);
            Check(trace == "A" + std::string(group,'E') + "Z" && clocks == 2);
            Check(s.elapsed_ns == 257 && s.ns_per_call == 257. / group); ++cases;
            trace.clear(); clocks = 0;
            const Sample untimed = RunGroup(group,false,encode,clock);
            Check(trace == std::string(group,'E') && clocks == 0);
            Check(untimed.elapsed_ns == 0 && untimed.ns_per_call == 0.); ++cases;
            for (int64_t elapsed : {INT64_C(1),INT64_C(255),INT64_C(256),INT64_C(257),
                                    INT64_C(9007199254740991)}) {
                const Sample value = Normalize(0,elapsed,group);
                Check(value.elapsed_ns == static_cast<uint64_t>(elapsed));
                Check(value.ns_per_call * group == static_cast<double>(elapsed)); ++cases;
            }
            for (auto endpoints : {std::make_pair(INT64_C(0),INT64_C(0)),
                                   std::make_pair(INT64_C(2),INT64_C(1)),
                                   std::make_pair(INT64_MIN,INT64_MAX),
                                   std::make_pair(INT64_C(0),INT64_MAX),
                                   std::make_pair(INT64_C(0),INT64_C(9007199254740992))}) {
                Reject([&]() { Normalize(endpoints.first,endpoints.second,group); }); ++cases;
            }
            trace.clear(); clocks = 0;
            Reject([&]() { RunGroup(group,true,[&]() { trace+='E'; throw std::runtime_error("encode"); },clock); });
            Check(trace == "AE" && clocks == 1); ++cases;
            trace.clear();
            Reject([&]() { RunGroup(group,true,encode,[&]() -> int64_t {
                trace+='A'; throw std::runtime_error("start clock"); }); });
            Check(trace == "A"); ++cases;
            trace.clear(); clocks = 0;
            Reject([&]() { RunGroup(group,true,encode,[&]() -> int64_t {
                if (++clocks == 2) { trace+='Z'; throw std::runtime_error("end clock"); }
                trace+='A'; return 1;
            }); });
            Check(trace == "A"+std::string(group,'E')+"Z" && clocks == 2); ++cases;
            Check(Normalize(INT64_MAX-1,INT64_MAX,group).elapsed_ns == 1); ++cases;
        }
        for (unsigned group : {0U,2U,255U,257U,std::numeric_limits<unsigned>::max()}) {
            unsigned calls = 0;
            Reject([&]() { RunGroup(group,true,[&]() { ++calls; },[&]() -> int64_t { ++calls; return 0; }); });
            Check(calls == 0); ++cases;
        }
        std::printf("{\"schema\":\"paired-group-unit/v1\",\"cases\":%u,\"timed\":false}\n",cases);
        return 0;
    } catch (const std::exception& e) { std::fprintf(stderr,"%s\n",e.what()); return 1; }
}
