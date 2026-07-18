/*
    Finite comparative benchmark for packed Leopard FF8 and the experimental
    native plane-sliced FF8 XOR-circuit backend.
*/

#include "../LeopardCommon.h"
#include "../LeopardFF8.h"
#include "../LeopardFF8Xor.h"
#include "../LeopardFF8XorTranspose.h"
#include "../leopard.h"
#include "../leopard_ff8xor.h"
#include "FF8XorCacheColoredBuffers.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <sstream>
#include <string>
#include <vector>

#if defined(__linux__)
    #include <linux/perf_event.h>
    #include <sched.h>
    #include <sys/ioctl.h>
    #include <sys/syscall.h>
    #include <sys/utsname.h>
    #include <unistd.h>
#endif

namespace {

static volatile uint64_t BenchmarkSink = 0;
static leopard::ff8xor::transpose::Mode BoundaryTransposeMode =
    leopard::ff8xor::transpose::Mode::Auto;

struct Options
{
    bool quick;
    bool csv;
    bool json;
    bool include_transpose;
    bool portable_transpose;
    bool counters;
    bool abba;
    bool pin;
    bool cache_color;
    leopard::ff8xor::KernelMode ff8xor_mode;
    std::string ff8xor_mode_name;
    leopard::ff8xor::FourBufferMode four_buffer_mode;
    std::string four_buffer_mode_name;
    int pin_cpu;
    unsigned warmups;
    unsigned iterations;
    double minimum_sample_usec;

    Options()
        : quick(false)
        , csv(false)
        , json(false)
        , include_transpose(false)
        , portable_transpose(false)
        , counters(false)
        , abba(false)
        , pin(true)
        , cache_color(false)
        , ff8xor_mode(leopard::ff8xor::KernelMode::Auto)
        , ff8xor_mode_name("auto")
        , four_buffer_mode(leopard::ff8xor::FourBufferMode::Disabled)
        , four_buffer_mode_name("disabled")
        , pin_cpu(-1)
        , warmups(2)
        , iterations(7)
        , minimum_sample_usec(1000.)
    {
    }
};

struct Environment
{
    std::string compiler;
    std::string build_type;
    std::string build_flags;
    std::string cpu;
    std::string simd;
    std::string ff8xor_mode_requested;
    std::string four_buffer_mode_requested;
    std::string operating_system;
    std::string affinity;
    std::string counter_backend;
};

struct CounterMetric
{
    std::string name;
    bool available;
    double median_per_call;
    std::string status;

    CounterMetric()
        : available(false)
        , median_per_call(0)
    {
    }
};

struct Timing
{
    double median_usec;
    double best_usec;
    unsigned calls_per_sample;
    unsigned sample_count;
    std::vector<CounterMetric> counters;
    double median_ipc;
    double median_effective_ghz;

    Timing()
        : median_usec(0)
        , best_usec(0)
        , calls_per_sample(0)
        , sample_count(0)
        , median_ipc(0)
        , median_effective_ghz(0)
    {
    }
};

struct Result
{
    std::string record;
    std::string backend;
    std::string operation;
    bool transpose_included;
    bool cache_coloring_applied;
    unsigned original_count;
    unsigned recovery_count;
    uint64_t buffer_bytes;
    unsigned loss_count;
    unsigned warmups;
    unsigned iterations;
    Timing timing;
    uint64_t input_bytes;
    uint64_t output_bytes;
    // The scheduled model is the pre-materialization transform traffic.  The
    // signed elision estimate comes from the most recent native codec call;
    // adjusted is scheduled - elided.  Keeping all three avoids presenting
    // the old static model as work that was actually performed.
    uint64_t modeled_payload_bytes;
    int64_t modeled_payload_bytes_elided;
    int64_t modeled_payload_bytes_adjusted;
    leopard::ff8xor::MaterializationStatistics materialization_statistics;
    leopard::ff8xor::FourBufferStatistics four_buffer_statistics;
    double ratio_vs_packed;
    unsigned locator_shift;
    std::string schedule_id;
    std::string measurement_order;
    std::string note;

    Result()
        : transpose_included(false)
        , cache_coloring_applied(false)
        , original_count(0)
        , recovery_count(0)
        , buffer_bytes(0)
        , loss_count(0)
        , warmups(0)
        , iterations(0)
        , input_bytes(0)
        , output_bytes(0)
        , modeled_payload_bytes(0)
        , modeled_payload_bytes_elided(0)
        , modeled_payload_bytes_adjusted(0)
        , materialization_statistics()
        , four_buffer_statistics()
        , ratio_vs_packed(0)
        , locator_shift(std::numeric_limits<unsigned>::max())
    {
    }
};

typedef leopard_ff8xor_test::BufferSet BufferSet;

class Random
{
public:
    explicit Random(uint64_t seed) : State(seed ? seed : UINT64_C(1)) {}

    uint64_t Next()
    {
        uint64_t value = State;
        value ^= value >> 12;
        value ^= value << 25;
        value ^= value >> 27;
        State = value;
        return value * UINT64_C(2685821657736338717);
    }

private:
    uint64_t State;
};

// packed[8*g + lane] bit plane_index ==
// plane[plane_index * (bytes / 8) + g] bit lane.
static void PackedToPlane(
    const uint8_t* packed,
    uint8_t* plane,
    uint64_t bytes)
{
    leopard::ff8xor::transpose::PackedToPlane(
        packed, plane, bytes, BoundaryTransposeMode);
}

static void PlaneToPacked(
    const uint8_t* plane,
    uint8_t* packed,
    uint64_t bytes)
{
    leopard::ff8xor::transpose::PlaneToPacked(
        plane, packed, bytes, BoundaryTransposeMode);
}

static const char* AutoBoundaryTransposeNote()
{
    const bool bitalg = leopard::ff8xor::transpose::IsModeAvailable(
        leopard::ff8xor::transpose::Mode::Avx512Bitalg);
    const bool vbmi = leopard::ff8xor::transpose::IsModeAvailable(
        leopard::ff8xor::transpose::Mode::Avx512Vbmi);
    if (bitalg && vbmi)
        return "auto: BITALG ZMM packed-to-plane; VBMI ZMM plane-to-packed; AVX2/portable tails";
    if (bitalg)
        return "auto: BITALG ZMM packed-to-plane; AVX2/portable inverse and tails";
    if (vbmi)
        return "auto: AVX2/portable forward; VBMI ZMM plane-to-packed; AVX2/portable tails";
    if (leopard::ff8xor::transpose::IsModeAvailable(
            leopard::ff8xor::transpose::Mode::Avx2))
        return "auto-dispatched blocked AVX2 boundary transpose with portable tails";
    return "portable boundary transpose fallback";
}

static void FillRandom(BufferSet& buffers, uint64_t bytes, uint64_t seed)
{
    Random random(seed);
    for (unsigned shard = 0; shard < buffers.Count(); ++shard)
    {
        for (uint64_t offset = 0; offset < bytes; offset += 8)
        {
            const uint64_t value = random.Next();
            memcpy(buffers[shard] + offset, &value, sizeof(value));
        }
    }
}

static bool Equal(const void* a, const void* b, uint64_t bytes)
{
    return memcmp(a, b, static_cast<size_t>(bytes)) == 0;
}

static double Median(std::vector<double> values)
{
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    if ((values.size() & 1) != 0)
        return values[middle];
    return (values[middle - 1] + values[middle]) * 0.5;
}

struct PerfDescriptor
{
    const char* name;
    uint32_t type;
    uint64_t config;
};

static std::vector<PerfDescriptor> PerfDescriptors()
{
    std::vector<PerfDescriptor> result;
#if defined(__linux__)
    const uint64_t read_miss =
        static_cast<uint64_t>(PERF_COUNT_HW_CACHE_OP_READ) << 8 |
        static_cast<uint64_t>(PERF_COUNT_HW_CACHE_RESULT_MISS) << 16;
    result.push_back({ "cycles", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_CPU_CYCLES });
    result.push_back({ "instructions", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_INSTRUCTIONS });
    result.push_back({ "reference_cycles", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_REF_CPU_CYCLES });
    result.push_back({ "cache_references", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_CACHE_REFERENCES });
    result.push_back({ "cache_misses", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_CACHE_MISSES });
    result.push_back({ "frontend_stalled_cycles", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_STALLED_CYCLES_FRONTEND });
    result.push_back({ "backend_stalled_cycles", PERF_TYPE_HARDWARE,
        PERF_COUNT_HW_STALLED_CYCLES_BACKEND });
    result.push_back({ "l1d_load_misses", PERF_TYPE_HW_CACHE,
        PERF_COUNT_HW_CACHE_L1D | read_miss });
    result.push_back({ "dtlb_load_misses", PERF_TYPE_HW_CACHE,
        PERF_COUNT_HW_CACHE_DTLB | read_miss });
    result.push_back({ "itlb_load_misses", PERF_TYPE_HW_CACHE,
        PERF_COUNT_HW_CACHE_ITLB | read_miss });
#else
    // Keep the schema identical on non-Linux hosts.
    static const char* names[] = {
        "cycles", "instructions", "reference_cycles", "cache_references",
        "cache_misses", "frontend_stalled_cycles", "backend_stalled_cycles",
        "l1d_load_misses", "dtlb_load_misses", "itlb_load_misses"
    };
    for (size_t i = 0; i < sizeof(names) / sizeof(names[0]); ++i)
        result.push_back({ names[i], 0, 0 });
#endif
    return result;
}

class PerfCounterSet
{
    struct Event
    {
        PerfDescriptor descriptor;
        int fd;
        bool available;
        std::string status;
        std::vector<double> samples;

        Event()
            : fd(-1)
            , available(false)
        {
        }
    };

    struct CounterGroup
    {
        int leader_fd;
        std::vector<size_t> event_indices;

        CounterGroup()
            : leader_fd(-1)
        {
        }
    };

public:
    explicit PerfCounterSet(bool requested)
    {
        const std::vector<PerfDescriptor> descriptors = PerfDescriptors();
        Events.resize(descriptors.size());
        for (size_t i = 0; i < descriptors.size(); ++i)
        {
            Events[i].descriptor = descriptors[i];
            if (!requested)
                Events[i].status = "disabled (pass --counters)";
        }
        if (!requested)
            return;
#if defined(__linux__)
        // A perf group is scheduled only when all of its members fit on the
        // PMU at once.  Typical x86 CPUs expose fewer than ten programmable
        // counters, so one group containing this complete descriptor set may
        // never run.  Use small groups that can multiplex as units.  The first
        // group deliberately keeps cycles, instructions, and reference cycles
        // together so IPC and frequency ratios share an exact payload window.
        static const size_t kMaximumEventsPerGroup = 3;
        for (size_t group_begin = 0;
             group_begin < descriptors.size();
             group_begin += kMaximumEventsPerGroup)
        {
            CounterGroup group;
            const size_t group_end = std::min(
                descriptors.size(), group_begin + kMaximumEventsPerGroup);
            for (size_t i = group_begin; i < group_end; ++i)
            {
                const bool leader = group.leader_fd < 0;
                struct perf_event_attr attributes;
                memset(&attributes, 0, sizeof(attributes));
                attributes.type = descriptors[i].type;
                attributes.size = sizeof(attributes);
                attributes.config = descriptors[i].config;
                attributes.disabled = leader ? 1 : 0;
                attributes.exclude_kernel = 1;
                attributes.exclude_hv = 1;
                if (leader)
                {
                    attributes.read_format = PERF_FORMAT_GROUP |
                        PERF_FORMAT_TOTAL_TIME_ENABLED |
                        PERF_FORMAT_TOTAL_TIME_RUNNING;
                }
                const int fd = static_cast<int>(syscall(
                    __NR_perf_event_open, &attributes, 0, -1,
                    leader ? -1 : group.leader_fd, 0));
                if (fd < 0)
                {
                    Events[i].status = std::string("unavailable: ") +
                        strerror(errno);
                    continue;
                }

                Events[i].fd = fd;
                Events[i].available = true;
                Events[i].status = "available";
                group.event_indices.push_back(i);
                if (leader)
                    group.leader_fd = fd;
            }
            if (group.leader_fd >= 0)
                Groups.push_back(group);
        }
#else
        for (size_t i = 0; i < descriptors.size(); ++i)
            Events[i].status = "unavailable: perf_event_open is Linux-only";
#endif
    }

    ~PerfCounterSet()
    {
#if defined(__linux__)
        for (size_t i = 0; i < Events.size(); ++i)
        {
            if (Events[i].fd >= 0)
                close(Events[i].fd);
        }
#endif
    }

    void Begin()
    {
#if defined(__linux__)
        // Reset every group before enabling any of them, then enable all
        // groups in a tight pass immediately before the payload region.
        for (size_t group_index = 0; group_index < Groups.size(); ++group_index)
        {
            CounterGroup& group = Groups[group_index];
            if (group.event_indices.empty() ||
                !Events[group.event_indices[0]].available)
                continue;
            if (ioctl(group.leader_fd, PERF_EVENT_IOC_RESET,
                    PERF_IOC_FLAG_GROUP) != 0)
            {
                MarkUnavailable(group,
                    std::string("unavailable while resetting PMU group: ") +
                    strerror(errno));
            }
        }
        // Enable the primary cycles/instructions/reference group last.  End()
        // disables that group first, so its counter window brackets the same
        // payload interval used for wall-clock frequency and IPC calculations
        // instead of including ioctls for the auxiliary groups.
        for (size_t group_index = Groups.size(); group_index-- > 0;)
        {
            CounterGroup& group = Groups[group_index];
            if (group.event_indices.empty() ||
                !Events[group.event_indices[0]].available)
                continue;
            if (ioctl(group.leader_fd, PERF_EVENT_IOC_ENABLE,
                    PERF_IOC_FLAG_GROUP) != 0)
            {
                MarkUnavailable(group,
                    std::string("unavailable while enabling PMU group: ") +
                    strerror(errno));
            }
        }
#endif
    }

    void End(unsigned calls_per_sample, double usec_per_call)
    {
        std::vector<double> values;
#if defined(__linux__)
        // Stop every group before allocating or reading result buffers.  This
        // keeps counter post-processing out of all payload measurement windows.
        for (size_t group_index = 0; group_index < Groups.size(); ++group_index)
        {
            CounterGroup& counter_group = Groups[group_index];
            if (counter_group.event_indices.empty() ||
                !Events[counter_group.event_indices[0]].available)
                continue;

            if (ioctl(counter_group.leader_fd, PERF_EVENT_IOC_DISABLE,
                    PERF_IOC_FLAG_GROUP) != 0)
            {
                MarkUnavailable(counter_group,
                    std::string("unavailable while disabling PMU group: ") +
                    strerror(errno));
            }
        }

        values.assign(Events.size(), 0);
        for (size_t group_index = 0; group_index < Groups.size(); ++group_index)
        {
            CounterGroup& counter_group = Groups[group_index];
            if (counter_group.event_indices.empty() ||
                !Events[counter_group.event_indices[0]].available)
                continue;
            // PERF_FORMAT_GROUP returns nr, enabled, running, then values in
            // group-open order.  Each group is small enough to schedule as a
            // unit; enabled/running scaling accounts for multiplexing between
            // the independent groups.
            std::vector<uint64_t> group(
                3 + counter_group.event_indices.size(), 0);
            const size_t expected_bytes = group.size() * sizeof(uint64_t);
            const ssize_t bytes = read(
                counter_group.leader_fd, &group[0], expected_bytes);
            if (bytes != static_cast<ssize_t>(expected_bytes) ||
                group[0] != counter_group.event_indices.size() ||
                group[2] == 0)
            {
                MarkUnavailable(counter_group,
                    "unavailable while reading PMU group");
                continue;
            }

            const double scale = group[1] == group[2]
                ? 1. : static_cast<double>(group[1]) /
                    static_cast<double>(group[2]);
            for (size_t i = 0;
                 i < counter_group.event_indices.size(); ++i)
            {
                const size_t event_index =
                    counter_group.event_indices[i];
                values[event_index] =
                    static_cast<double>(group[3 + i]) * scale /
                    calls_per_sample;
                Events[event_index].samples.push_back(values[event_index]);
            }
        }
#else
        (void)calls_per_sample;
        values.assign(Events.size(), 0);
#endif
        if (Events.size() >= 2 && Events[0].available && Events[1].available &&
            values[0] > 0)
            IPCSamples.push_back(values[1] / values[0]);
        if (!Events.empty() && Events[0].available && usec_per_call > 0)
            FrequencySamples.push_back(values[0] / usec_per_call / 1000.);
    }

    void Finish(Timing& timing) const
    {
        timing.counters.clear();
        timing.counters.reserve(Events.size());
        for (size_t i = 0; i < Events.size(); ++i)
        {
            CounterMetric metric;
            metric.name = Events[i].descriptor.name;
            metric.available = Events[i].available && !Events[i].samples.empty();
            metric.status = Events[i].status;
            if (metric.available)
                metric.median_per_call = Median(Events[i].samples);
            timing.counters.push_back(metric);
        }
        timing.median_ipc = Median(IPCSamples);
        timing.median_effective_ghz = Median(FrequencySamples);
    }

    bool AnyAvailable() const
    {
        for (size_t i = 0; i < Events.size(); ++i)
        {
            if (Events[i].available)
                return true;
        }
        return false;
    }

private:
    void MarkUnavailable(CounterGroup& group, const std::string& status)
    {
        for (size_t i = 0; i < group.event_indices.size(); ++i)
        {
            Event& event = Events[group.event_indices[i]];
            event.available = false;
            event.status = status;
        }
    }

    std::vector<Event> Events;
    std::vector<CounterGroup> Groups;
    std::vector<double> IPCSamples;
    std::vector<double> FrequencySamples;
};

template <typename Function>
static bool Measure(
    const Options& options,
    Function function,
    Timing& timing,
    LeopardResult& result)
{
    typedef std::chrono::steady_clock Clock;
    for (unsigned i = 0; i < options.warmups; ++i)
    {
        result = function();
        if (result != Leopard_Success)
            return false;
    }

    // Measure one untimed calibration call, then batch fast operations so a
    // measured sample is long enough that clock overhead does not dominate.
    // Reported times below are divided back down to one whole-buffer call.
    const Clock::time_point calibration_begin = Clock::now();
    result = function();
    const Clock::time_point calibration_end = Clock::now();
    if (result != Leopard_Success)
        return false;
    const double calibration_usec =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
            calibration_end - calibration_begin).count();

    static const unsigned kMaximumCallsPerSample = 1U << 20;
    unsigned calls_per_sample = 1;
    if (calibration_usec > 0. &&
        calibration_usec < options.minimum_sample_usec)
    {
        const double wanted =
            options.minimum_sample_usec / calibration_usec;
        if (wanted >= static_cast<double>(kMaximumCallsPerSample))
            calls_per_sample = kMaximumCallsPerSample;
        else
            calls_per_sample = std::max(1U,
                static_cast<unsigned>(wanted + 0.999999));
    }
    timing.calls_per_sample = calls_per_sample;
    timing.sample_count = options.iterations;
    PerfCounterSet counters(options.counters);

    std::vector<double> samples;
    samples.reserve(options.iterations);
    for (unsigned i = 0; i < options.iterations; ++i)
    {
        const Clock::time_point begin = Clock::now();
        for (unsigned call = 0; call < calls_per_sample; ++call)
        {
            result = function();
            if (result != Leopard_Success)
                return false;
        }
        const Clock::time_point end = Clock::now();
        const double usec =
            std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
            end - begin).count() / calls_per_sample;
        samples.push_back(usec);
    }

    std::sort(samples.begin(), samples.end());
    timing.best_usec = samples.front();
    const size_t middle = samples.size() / 2;
    if ((samples.size() & 1) != 0)
        timing.median_usec = samples[middle];
    else
        timing.median_usec = (samples[middle - 1] + samples[middle]) * 0.5;
    if (counters.AnyAvailable())
    {
        for (unsigned i = 0; i < options.iterations; ++i)
        {
            counters.Begin();
            const Clock::time_point begin = Clock::now();
            for (unsigned call = 0; call < calls_per_sample; ++call)
            {
                result = function();
                if (result != Leopard_Success)
                    return false;
            }
            const Clock::time_point end = Clock::now();
            const double usec =
                std::chrono::duration_cast<
                    std::chrono::duration<double, std::micro> >(
                        end - begin).count() / calls_per_sample;
            counters.End(calls_per_sample, usec);
        }
    }
    counters.Finish(timing);
    return true;
}

template <typename Function>
static bool CalibrateCalls(
    const Options& options,
    Function function,
    unsigned& calls_per_sample,
    LeopardResult& result)
{
    typedef std::chrono::steady_clock Clock;
    const Clock::time_point begin = Clock::now();
    result = function();
    const Clock::time_point end = Clock::now();
    if (result != Leopard_Success)
        return false;
    const double usec =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
            end - begin).count();
    static const unsigned kMaximumCallsPerSample = 1U << 20;
    calls_per_sample = 1;
    if (usec > 0. && usec < options.minimum_sample_usec)
    {
        const double wanted = options.minimum_sample_usec / usec;
        calls_per_sample = wanted >= kMaximumCallsPerSample
            ? kMaximumCallsPerSample
            : std::max(1U, static_cast<unsigned>(wanted + 0.999999));
    }
    return true;
}

template <typename Function>
static bool MeasureBlock(
    Function function,
    unsigned calls_per_sample,
    std::vector<double>& samples,
    LeopardResult& result)
{
    typedef std::chrono::steady_clock Clock;
    const Clock::time_point begin = Clock::now();
    for (unsigned call = 0; call < calls_per_sample; ++call)
    {
        result = function();
        if (result != Leopard_Success)
            return false;
    }
    const Clock::time_point end = Clock::now();
    const double usec =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
            end - begin).count() / calls_per_sample;
    samples.push_back(usec);
    return true;
}

template <typename Function>
static bool CountBlock(
    Function function,
    unsigned calls_per_sample,
    PerfCounterSet& counters,
    LeopardResult& result)
{
    typedef std::chrono::steady_clock Clock;
    counters.Begin();
    const Clock::time_point begin = Clock::now();
    for (unsigned call = 0; call < calls_per_sample; ++call)
    {
        result = function();
        if (result != Leopard_Success)
            return false;
    }
    const Clock::time_point end = Clock::now();
    const double usec =
        std::chrono::duration_cast<std::chrono::duration<double, std::micro> >(
            end - begin).count() / calls_per_sample;
    counters.End(calls_per_sample, usec);
    return true;
}

static void FinishTiming(
    std::vector<double> samples,
    unsigned calls_per_sample,
    const PerfCounterSet& counters,
    Timing& timing)
{
    std::sort(samples.begin(), samples.end());
    timing.best_usec = samples.front();
    timing.median_usec = Median(samples);
    timing.calls_per_sample = calls_per_sample;
    timing.sample_count = static_cast<unsigned>(samples.size());
    counters.Finish(timing);
}

template <typename FunctionA, typename FunctionB>
static bool MeasurePairABBA(
    const Options& options,
    FunctionA function_a,
    FunctionB function_b,
    Timing& timing_a,
    Timing& timing_b,
    LeopardResult& result)
{
    for (unsigned warmup = 0; warmup < options.warmups; ++warmup)
    {
        result = function_a();
        if (result != Leopard_Success)
            return false;
        result = function_b();
        if (result != Leopard_Success)
            return false;
        result = function_b();
        if (result != Leopard_Success)
            return false;
        result = function_a();
        if (result != Leopard_Success)
            return false;
    }

    unsigned calls_a = 1;
    unsigned calls_b = 1;
    if (!CalibrateCalls(options, function_a, calls_a, result) ||
        !CalibrateCalls(options, function_b, calls_b, result))
        return false;

    PerfCounterSet counters_a(options.counters);
    PerfCounterSet counters_b(options.counters);
    std::vector<double> samples_a;
    std::vector<double> samples_b;
    samples_a.reserve(options.iterations * 2);
    samples_b.reserve(options.iterations * 2);
    for (unsigned round = 0; round < options.iterations; ++round)
    {
        if (!MeasureBlock(function_a, calls_a, samples_a, result) ||
            !MeasureBlock(function_b, calls_b, samples_b, result) ||
            !MeasureBlock(function_b, calls_b, samples_b, result) ||
            !MeasureBlock(function_a, calls_a, samples_a, result))
            return false;
    }
    if (counters_a.AnyAvailable() || counters_b.AnyAvailable())
    {
        for (unsigned round = 0; round < options.iterations; ++round)
        {
            if (!CountBlock(function_a, calls_a, counters_a, result) ||
                !CountBlock(function_b, calls_b, counters_b, result) ||
                !CountBlock(function_b, calls_b, counters_b, result) ||
                !CountBlock(function_a, calls_a, counters_a, result))
                return false;
        }
    }
    FinishTiming(samples_a, calls_a, counters_a, timing_a);
    FinishTiming(samples_b, calls_b, counters_b, timing_b);
    return true;
}

static double Throughput(uint64_t bytes, double usec)
{
    if (usec <= 0)
        return 0;
    // Decimal MB/s: bytes/usec is numerically equal to MB/s.
    return static_cast<double>(bytes) / usec;
}

static std::string Csv(const std::string& text)
{
    std::string escaped;
    escaped.reserve(text.size() + 2);
    escaped.push_back('"');
    for (size_t i = 0; i < text.size(); ++i)
    {
        if (text[i] == '"')
            escaped.push_back('"');
        escaped.push_back(text[i]);
    }
    escaped.push_back('"');
    return escaped;
}

static std::string Json(const std::string& text)
{
    std::ostringstream output;
    output << '"';
    for (size_t i = 0; i < text.size(); ++i)
    {
        const unsigned char value = static_cast<unsigned char>(text[i]);
        switch (value)
        {
        case '"': output << "\\\""; break;
        case '\\': output << "\\\\"; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (value < 0x20)
            {
                output << "\\u" << std::hex << std::setw(4)
                    << std::setfill('0') << static_cast<unsigned>(value)
                    << std::dec << std::setfill(' ');
            }
            else
                output << static_cast<char>(value);
        }
    }
    output << '"';
    return output.str();
}

static const CounterMetric* FindCounter(
    const Timing& timing,
    const char* name)
{
    for (size_t i = 0; i < timing.counters.size(); ++i)
    {
        if (timing.counters[i].name == name)
            return &timing.counters[i];
    }
    return NULL;
}

static std::string CounterStatus(const Timing& timing)
{
    std::ostringstream status;
    for (size_t i = 0; i < timing.counters.size(); ++i)
    {
        if (i != 0)
            status << "; ";
        status << timing.counters[i].name << '=' << timing.counters[i].status;
    }
    return status.str();
}

class Reporter
{
public:
    Reporter(bool csv, bool json, const Environment& environment)
        : CSV(csv), Env(environment)
        , JSON(json)
    {
    }

    bool MachineReadable() const { return CSV || JSON; }

    void Begin(const Options& options)
    {
        const leopard::ff8xor::CircuitStatistics multiply =
            leopard::ff8xor::GetMultiplyCircuitStatistics();
        const leopard::ff8xor::CircuitStatistics fft =
            leopard::ff8xor::GetFFTCircuitStatistics();
        const leopard::ff8xor::CircuitStatistics ifft =
            leopard::ff8xor::GetIFFTCircuitStatistics();
        MeasurementOrder = options.abba
            ? "end-to-end pairs ABBA; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential"
            : "packed-then-ff8xor end-to-end; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential";
        if (CSV)
        {
            std::cout
                << "record,backend,operation,transpose_included,cache_coloring_applied,"
                << "k,r,buffer_bytes,loss_count,"
                << "warmups,iterations,calls_per_sample,median_us,best_us,median_input_MBps,"
                << "best_input_MBps,median_output_MBps,best_output_MBps,"
                << "speed_ratio_vs_packed,compiler,build_type,cpu,simd,"
                << "ff8xor_mode_requested,four_buffer_mode_requested,checksum,"
                << "cost_profile_id,cost_profile_checksum,gate_min,"
                << "gate_max,gate_average,note,schedule_id,"
                << "modeled_payload_bytes_scheduled,modeled_payload_bytes_elided,"
                << "modeled_payload_bytes_adjusted,"
                << "materialization_deferred_zero_fills,"
                << "materialization_added_zero_fills,"
                << "materialization_butterflies_skipped,"
                << "materialization_butterflies_reduced,"
                << "materialization_xors_skipped,"
                << "materialization_xors_replaced_by_copies,"
                << "materialization_identity_operations_elided,"
                << "four_buffer_fused_units,"
                << "four_buffer_payload_bytes_elided,"
                << "locator_shift,measurement_order,operating_system,affinity,build_flags,"
                << "counter_backend,cycles,instructions,reference_cycles,cache_references,"
                << "cache_misses,frontend_stalled_cycles,backend_stalled_cycles,"
                << "l1d_load_misses,dtlb_load_misses,itlb_load_misses,ipc,effective_ghz,"
                << "counter_status\n";
            GateCSV("multiply",
                multiply.MinimumGateCount,
                multiply.MaximumGateCount,
                multiply.AverageGateCount);
            GateCSV("fft",
                fft.MinimumGateCount,
                fft.MaximumGateCount,
                fft.AverageGateCount);
            GateCSV("ifft",
                ifft.MinimumGateCount,
                ifft.MaximumGateCount,
                ifft.AverageGateCount);
            return;
        }

        if (JSON)
        {
            std::cout << "{" << Json("record") << ':' << Json("metadata")
                << ',' << Json("schema") << ':'
                << Json("leopard.ff8xor.benchmark.jsonl.v2")
                << ',' << Json("compiler") << ':' << Json(Env.compiler)
                << ',' << Json("build_type") << ':' << Json(Env.build_type)
                << ',' << Json("build_flags") << ':' << Json(Env.build_flags)
                << ',' << Json("cpu") << ':' << Json(Env.cpu)
                << ',' << Json("simd") << ':' << Json(Env.simd)
                << ',' << Json("ff8xor_mode_requested") << ':'
                << Json(Env.ff8xor_mode_requested)
                << ',' << Json("four_buffer_mode_requested") << ':'
                << Json(Env.four_buffer_mode_requested)
                << ',' << Json("circuit_checksum") << ':'
                << Json(leopard::ff8xor::GetCircuitChecksum())
                << ',' << Json("circuit_cost_profile_id") << ':'
                << Json(leopard::ff8xor::GetCircuitCostProfileId())
                << ',' << Json("circuit_cost_profile_checksum") << ':'
                << Json(leopard::ff8xor::GetCircuitCostProfileChecksum())
                << ',' << Json("operating_system") << ':'
                << Json(Env.operating_system)
                << ',' << Json("affinity") << ':' << Json(Env.affinity)
                << ',' << Json("counter_backend") << ':'
                << Json(Env.counter_backend)
                << ',' << Json("quick") << ':'
                << (options.quick ? "true" : "false")
                << ',' << Json("warmups") << ':' << options.warmups
                << ',' << Json("iterations") << ':' << options.iterations
                << ',' << Json("minimum_sample_usec") << ':'
                << options.minimum_sample_usec
                << ',' << Json("transpose_included") << ':'
                << (options.include_transpose ? "true" : "false")
                << ',' << Json("pmu_requested") << ':'
                << (options.counters ? "true" : "false")
                << ',' << Json("cache_coloring_requested") << ':'
                << (options.cache_color ? "true" : "false")
                << ',' << Json("measurement_order") << ':'
                << Json(options.abba
                    ? "end-to-end pairs ABBA; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential"
                    : "packed-then-ff8xor end-to-end; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential")
                << "}\n";
            GateJSON("multiply", multiply);
            GateJSON("fft", fft);
            GateJSON("ifft", ifft);
            return;
        }

        std::cout << "FF8 XOR-circuit comparative benchmark\n"
            << "compiler: " << Env.compiler << '\n'
            << "build: " << Env.build_type << '\n'
            << "cpu: " << Env.cpu << '\n'
            << "simd: " << Env.simd << '\n'
            << "ff8xor mode requested: " << Env.ff8xor_mode_requested << '\n'
            << "four-buffer mode requested: "
            << Env.four_buffer_mode_requested << '\n'
            << "os: " << Env.operating_system << '\n'
            << "affinity: " << Env.affinity << '\n'
            << "build flags: " << Env.build_flags << '\n'
            << "counter backend: " << Env.counter_backend << '\n'
            << "mode: " << (options.quick ? "quick" : "full")
            << ", warmups=" << options.warmups
            << ", iterations=" << options.iterations
            << ", minimum sample=" << options.minimum_sample_usec << " us"
            << ", packed-boundary transpose measurements="
            << (options.include_transpose ? "included" : "disabled") << '\n'
            << "measurement order: "
            << (options.abba
                ? "end-to-end pairs ABBA; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential"
                : "packed then ff8xor end-to-end; XOR-batch pairs ABBA; locator-selector pairs ABBA; transpose pairs ABBA; other micro/boundary sequential") << '\n'
            << "PMU counters: "
            << (options.counters ? "requested" : "disabled") << '\n'
            << "ff8xor allocations: "
            << (options.cache_color
                ? "4-KiB-relative colors for fully aliasing plane strides"
                : "allocator default (uncolored)") << '\n'
            << "timing: calls/sample is auto-calibrated and each result is per call\n"
            << "circuit checksum: "
            << leopard::ff8xor::GetCircuitChecksum() << '\n'
            << "circuit cost profile: "
            << leopard::ff8xor::GetCircuitCostProfileId() << " ("
            << leopard::ff8xor::GetCircuitCostProfileChecksum() << ")\n"
            << std::fixed << std::setprecision(3)
            << "gates multiply="
            << multiply.MinimumGateCount << ".."
            << multiply.MaximumGateCount << " avg="
            << multiply.AverageGateCount
            << ", fft=" << fft.MinimumGateCount << ".."
            << fft.MaximumGateCount << " avg="
            << fft.AverageGateCount
            << ", ifft=" << ifft.MinimumGateCount << ".."
            << ifft.MaximumGateCount << " avg="
            << ifft.AverageGateCount << "\n\n";
    }

    void Print(const Result& result)
    {
        const double median_input = Throughput(
            result.input_bytes, result.timing.median_usec);
        const double best_input = Throughput(
            result.input_bytes, result.timing.best_usec);
        const double median_output = Throughput(
            result.output_bytes, result.timing.median_usec);
        const double best_output = Throughput(
            result.output_bytes, result.timing.best_usec);

        if (CSV)
        {
            std::cout << Csv(result.record) << ','
                << Csv(result.backend) << ',' << Csv(result.operation) << ','
                << (result.transpose_included ? 1 : 0) << ','
                << (result.cache_coloring_applied ? 1 : 0) << ','
                << result.original_count << ',' << result.recovery_count << ','
                << result.buffer_bytes << ',' << result.loss_count << ','
                << result.warmups << ',' << result.iterations << ','
                << result.timing.calls_per_sample << ','
                << std::fixed << std::setprecision(3)
                << result.timing.median_usec << ',' << result.timing.best_usec << ','
                << median_input << ',' << best_input << ','
                << median_output << ',' << best_output << ','
                << result.ratio_vs_packed << ','
                << Csv(Env.compiler) << ',' << Csv(Env.build_type) << ','
                << Csv(Env.cpu) << ',' << Csv(Env.simd) << ','
                << Csv(Env.ff8xor_mode_requested) << ','
                << Csv(Env.four_buffer_mode_requested) << ','
                << Csv("") << ",,,,,," << Csv(result.note) << ','
                << Csv(result.schedule_id) << ','
                << result.modeled_payload_bytes << ','
                << result.modeled_payload_bytes_elided << ','
                << result.modeled_payload_bytes_adjusted << ','
                << result.materialization_statistics.DeferredZeroFills << ','
                << result.materialization_statistics.AddedZeroFills << ','
                << result.materialization_statistics.ButterfliesSkipped << ','
                << result.materialization_statistics.ButterfliesReduced << ','
                << result.materialization_statistics.XorsSkipped << ','
                << result.materialization_statistics.XorsReplacedByCopies << ','
                << result.materialization_statistics.IdentityOperationsElided << ','
                << result.four_buffer_statistics.FusedUnits << ','
                << result.four_buffer_statistics.EstimatedPayloadBytesElided << ',';
            if (result.locator_shift != std::numeric_limits<unsigned>::max())
                std::cout << result.locator_shift;
            std::cout << ',' << Csv(result.measurement_order) << ','
                << Csv(Env.operating_system) << ',' << Csv(Env.affinity) << ','
                << Csv(Env.build_flags) << ',' << Csv(Env.counter_backend);
            static const char* counter_names[] = {
                "cycles", "instructions", "reference_cycles", "cache_references",
                "cache_misses", "frontend_stalled_cycles", "backend_stalled_cycles",
                "l1d_load_misses", "dtlb_load_misses", "itlb_load_misses"
            };
            for (size_t i = 0;
                 i < sizeof(counter_names) / sizeof(counter_names[0]); ++i)
            {
                std::cout << ',';
                const CounterMetric* metric =
                    FindCounter(result.timing, counter_names[i]);
                if (metric && metric->available)
                    std::cout << metric->median_per_call;
            }
            std::cout << ',';
            if (result.timing.median_ipc > 0)
                std::cout << result.timing.median_ipc;
            std::cout << ',';
            if (result.timing.median_effective_ghz > 0)
                std::cout << result.timing.median_effective_ghz;
            std::cout << ',' << Csv(CounterStatus(result.timing)) << '\n';
            return;
        }

        if (JSON)
        {
            std::cout << "{" << Json("record") << ':' << Json(result.record)
                << ',' << Json("backend") << ':' << Json(result.backend)
                << ',' << Json("operation") << ':' << Json(result.operation)
                << ',' << Json("transpose_included") << ':'
                << (result.transpose_included ? "true" : "false")
                << ',' << Json("cache_coloring_applied") << ':'
                << (result.cache_coloring_applied ? "true" : "false")
                << ',' << Json("k") << ':' << result.original_count
                << ',' << Json("r") << ':' << result.recovery_count
                << ',' << Json("buffer_bytes") << ':' << result.buffer_bytes
                << ',' << Json("loss_count") << ':' << result.loss_count
                << ',' << Json("schedule_id") << ':' << Json(result.schedule_id)
                << ',' << Json("measurement_order") << ':'
                << Json(result.measurement_order)
                << ',' << Json("warmups") << ':' << result.warmups
                << ',' << Json("iterations") << ':' << result.iterations
                << ',' << Json("calls_per_sample") << ':'
                << result.timing.calls_per_sample
                << ',' << Json("median_us") << ':'
                << result.timing.median_usec
                << ',' << Json("best_us") << ':' << result.timing.best_usec
                << ',' << Json("median_input_MBps") << ':' << median_input
                << ',' << Json("best_input_MBps") << ':' << best_input
                << ',' << Json("median_output_MBps") << ':' << median_output
                << ',' << Json("best_output_MBps") << ':' << best_output
                << ',' << Json("speed_ratio_vs_packed") << ':'
                << result.ratio_vs_packed
                << ',' << Json("modeled_payload_bytes_scheduled") << ':'
                << result.modeled_payload_bytes
                << ',' << Json("modeled_payload_bytes_elided") << ':'
                << result.modeled_payload_bytes_elided
                << ',' << Json("modeled_payload_bytes_adjusted") << ':'
                << result.modeled_payload_bytes_adjusted
                << ',' << Json("materialization_deferred_zero_fills") << ':'
                << result.materialization_statistics.DeferredZeroFills
                << ',' << Json("materialization_added_zero_fills") << ':'
                << result.materialization_statistics.AddedZeroFills
                << ',' << Json("materialization_butterflies_skipped") << ':'
                << result.materialization_statistics.ButterfliesSkipped
                << ',' << Json("materialization_butterflies_reduced") << ':'
                << result.materialization_statistics.ButterfliesReduced
                << ',' << Json("materialization_xors_skipped") << ':'
                << result.materialization_statistics.XorsSkipped
                << ',' << Json("materialization_xors_replaced_by_copies") << ':'
                << result.materialization_statistics.XorsReplacedByCopies
                << ',' << Json("materialization_identity_operations_elided") << ':'
                << result.materialization_statistics.IdentityOperationsElided
                << ',' << Json("four_buffer_fused_units") << ':'
                << result.four_buffer_statistics.FusedUnits
                << ',' << Json("four_buffer_payload_bytes_elided") << ':'
                << result.four_buffer_statistics.EstimatedPayloadBytesElided
                << ',' << Json("locator_shift") << ':';
            if (result.locator_shift == std::numeric_limits<unsigned>::max())
                std::cout << "null";
            else
                std::cout << result.locator_shift;
            std::cout << ',' << Json("counters") << ":{";
            for (size_t i = 0; i < result.timing.counters.size(); ++i)
            {
                const CounterMetric& metric = result.timing.counters[i];
                if (i != 0)
                    std::cout << ',';
                std::cout << Json(metric.name) << ":{" << Json("value_per_call")
                    << ':';
                if (metric.available)
                    std::cout << metric.median_per_call;
                else
                    std::cout << "null";
                std::cout << ',' << Json("status") << ':' << Json(metric.status)
                    << '}';
            }
            std::cout << "}," << Json("ipc") << ':';
            if (result.timing.median_ipc > 0)
                std::cout << result.timing.median_ipc;
            else
                std::cout << "null";
            std::cout << ',' << Json("effective_ghz") << ':';
            if (result.timing.median_effective_ghz > 0)
                std::cout << result.timing.median_effective_ghz;
            else
                std::cout << "null";
            std::cout << ',' << Json("note") << ':' << Json(result.note)
                << "}\n";
            return;
        }

        std::cout << std::left << std::setw(25) << result.backend
            << std::setw(10) << result.operation << std::right
            << " k=" << std::setw(3) << result.original_count
            << " r=" << std::setw(3) << result.recovery_count
            << " B=" << std::setw(8) << result.buffer_bytes
            << " loss=" << std::setw(3) << result.loss_count
            << " warmups=" << std::setw(2) << result.warmups
            << " samples=" << std::setw(2) << result.iterations
            << " calls/sample=" << std::setw(6) << result.timing.calls_per_sample
            << " transpose=" << (result.transpose_included ? "yes" : "no ")
            << " cache-color="
            << (result.cache_coloring_applied ? "yes" : "no ")
            << std::fixed << std::setprecision(3)
            << " median=" << std::setw(10) << result.timing.median_usec << " us"
            << " best=" << std::setw(10) << result.timing.best_usec << " us"
            << " input=" << std::setw(10) << median_input << " MB/s"
            << " output=" << std::setw(10) << median_output << " MB/s";
        if (result.modeled_payload_bytes != 0)
        {
            std::cout << " modeled_scheduled=" << result.modeled_payload_bytes
                << " B elided=" << result.modeled_payload_bytes_elided
                << " B adjusted=" << result.modeled_payload_bytes_adjusted
                << " B materialization_ops={deferred_zero:"
                << result.materialization_statistics.DeferredZeroFills
                << ",added_zero:"
                << result.materialization_statistics.AddedZeroFills
                << ",butterfly_skip:"
                << result.materialization_statistics.ButterfliesSkipped
                << ",butterfly_reduce:"
                << result.materialization_statistics.ButterfliesReduced
                << ",xor_skip:"
                << result.materialization_statistics.XorsSkipped
                << ",xor_copy:"
                << result.materialization_statistics.XorsReplacedByCopies
                << ",identity:"
                << result.materialization_statistics.IdentityOperationsElided
                << "} four_buffer={units:"
                << result.four_buffer_statistics.FusedUnits
                << ",bytes_elided:"
                << result.four_buffer_statistics.EstimatedPayloadBytesElided
                << '}';
        }
        if (result.timing.median_ipc > 0)
            std::cout << " IPC=" << result.timing.median_ipc;
        if (result.timing.median_effective_ghz > 0)
            std::cout << " GHz=" << result.timing.median_effective_ghz;
        if (result.ratio_vs_packed > 0)
            std::cout << " speed=" << result.ratio_vs_packed << "x packed";
        if (!result.note.empty())
            std::cout << " (" << result.note << ')';
        std::cout << '\n';
    }

    void Skip(
        unsigned original_count,
        unsigned recovery_count,
        uint64_t buffer_bytes,
        unsigned loss_count,
        const std::string& note)
    {
        Result result;
        result.record = "skip";
        result.backend = "all";
        result.operation = "skip";
        result.original_count = original_count;
        result.recovery_count = recovery_count;
        result.buffer_bytes = buffer_bytes;
        result.loss_count = loss_count;
        result.note = note;
        Print(result);
    }

private:
    void GateJSON(
        const char* operation,
        const leopard::ff8xor::CircuitStatistics& statistics)
    {
        std::cout << "{" << Json("record") << ':' << Json("circuit_metadata")
            << ',' << Json("family") << ':' << Json(operation)
            << ',' << Json("checksum") << ':'
            << Json(leopard::ff8xor::GetCircuitChecksum())
            << ',' << Json("cost_profile_id") << ':'
            << Json(leopard::ff8xor::GetCircuitCostProfileId())
            << ',' << Json("cost_profile_checksum") << ':'
            << Json(leopard::ff8xor::GetCircuitCostProfileChecksum())
            << ',' << Json("gate_min") << ':' << statistics.MinimumGateCount
            << ',' << Json("gate_max") << ':' << statistics.MaximumGateCount
            << ',' << Json("gate_average") << ':'
            << statistics.AverageGateCount
            << ',' << Json("depth_min") << ':' << statistics.MinimumDepth
            << ',' << Json("depth_max") << ':' << statistics.MaximumDepth
            << ',' << Json("depth_average") << ':' << statistics.AverageDepth
            << "}\n";
    }

    void GateCSV(const char* operation, unsigned minimum, unsigned maximum, double average)
    {
        std::cout << Csv("metadata") << ',' << Csv("generated") << ','
            << Csv(operation) << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
            << Csv(Env.compiler) << ',' << Csv(Env.build_type) << ','
            << Csv(Env.cpu) << ',' << Csv(Env.simd) << ','
            << Csv(Env.ff8xor_mode_requested) << ','
            << Csv(Env.four_buffer_mode_requested) << ','
            << Csv(leopard::ff8xor::GetCircuitChecksum()) << ','
            << Csv(leopard::ff8xor::GetCircuitCostProfileId()) << ','
            << Csv(leopard::ff8xor::GetCircuitCostProfileChecksum()) << ','
            << minimum << ',' << maximum << ',' << std::fixed
            << std::setprecision(6) << average << ',' << Csv("") << ','
            << Csv("") << ",0,0,0,0,0,0,0,0,0,0,0,0,,"
            << Csv(MeasurementOrder) << ','
            << Csv(Env.operating_system) << ',' << Csv(Env.affinity) << ','
            << Csv(Env.build_flags) << ',' << Csv(Env.counter_backend)
            ;
        // Ten PMU values plus IPC and effective GHz are unavailable for
        // circuit-metadata rows.  Keep the CSV schema rectangular.
        for (unsigned i = 0; i < 12; ++i)
            std::cout << ',';
        std::cout << ',' << Csv("") << '\n';
    }

    bool CSV;
    Environment Env;
    bool JSON;
    std::string MeasurementOrder;
};

static std::string CompilerName()
{
#if defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#elif defined(_MSC_VER)
    std::ostringstream text;
    text << "MSVC " << _MSC_VER;
    return text.str();
#else
    return "unknown";
#endif
}

static std::string BuildType()
{
#if defined(LEO_BENCH_CMAKE_BUILD_TYPE)
    std::string result = LEO_BENCH_CMAKE_BUILD_TYPE;
    result += " (";
    result +=
#if defined(NDEBUG)
        "NDEBUG";
#else
        "NDEBUG unset";
#endif
    result += ')';
    return result;
#elif defined(NDEBUG)
    return "optimized (NDEBUG)";
#elif defined(__OPTIMIZE__) || (defined(_MSC_VER) && !defined(_DEBUG))
    return "optimized (NDEBUG unset)";
#else
    return "Debug/unoptimized";
#endif
}

static std::string BuildFlags()
{
#if defined(LEO_BENCH_BUILD_FLAGS)
    return LEO_BENCH_BUILD_FLAGS;
#else
    return "unavailable (build system did not provide flags)";
#endif
}

static std::string Trim(const std::string& input)
{
    const size_t begin = input.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
        return "";
    const size_t end = input.find_last_not_of(" \t\r\n");
    return input.substr(begin, end - begin + 1);
}

static std::string CPUName()
{
#if defined(__linux__)
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line))
    {
        if (line.compare(0, 10, "model name") == 0)
        {
            const size_t colon = line.find(':');
            if (colon != std::string::npos)
                return Trim(line.substr(colon + 1));
        }
    }
#elif defined(_WIN32)
    const char* identifier = std::getenv("PROCESSOR_IDENTIFIER");
    if (identifier)
        return identifier;
#endif
    return "unavailable";
}

static std::string SIMDName()
{
    std::string packed = "portable";
#if defined(LEO_TRY_AVX2)
    if (leopard::CpuHasAVX2)
        packed = "AVX2";
    else
#endif
    if (leopard::CpuHasSSSE3)
        packed = "128-bit/SSSE3";
#if defined(LEO_TRY_NEON)
    else if (leopard::CpuHasNeon)
        packed = "NEON";
#endif
    return std::string("packed=") + packed + "; ff8xor=" +
        leopard::ff8xor::GetKernelBackendName();
}

static std::string OperatingSystemName()
{
#if defined(__linux__)
    struct utsname name;
    if (uname(&name) == 0)
    {
        std::ostringstream result;
        result << name.sysname << ' ' << name.release << ' ' << name.machine;
        return result.str();
    }
#elif defined(_WIN32)
    return "Windows";
#endif
    return "unavailable";
}

static std::string CounterBackendName()
{
#if defined(__linux__)
    std::ifstream input("/proc/sys/kernel/perf_event_paranoid");
    int paranoid = 0;
    if (input >> paranoid)
    {
        std::ostringstream result;
        result << "Linux perf_event_open; perf_event_paranoid=" << paranoid;
        return result.str();
    }
    return "Linux perf_event_open; perf_event_paranoid unavailable";
#else
    return "unavailable: perf_event_open is Linux-only";
#endif
}

static std::string ConfigureAffinity(const Options& options)
{
    if (!options.pin)
        return "not pinned (--no-pin)";
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0)
        return std::string("pinning unavailable: sched_getaffinity: ") +
            strerror(errno);

    int cpu = options.pin_cpu;
    if (cpu < 0)
    {
        for (int candidate = 0; candidate < CPU_SETSIZE; ++candidate)
        {
            if (CPU_ISSET(candidate, &allowed))
            {
                cpu = candidate;
                break;
            }
        }
    }
    if (cpu < 0 || cpu >= CPU_SETSIZE || !CPU_ISSET(cpu, &allowed))
    {
        std::ostringstream result;
        result << "pinning unavailable: requested CPU " << cpu
            << " is outside the allowed affinity set";
        return result.str();
    }

    cpu_set_t selected;
    CPU_ZERO(&selected);
    CPU_SET(cpu, &selected);
    if (sched_setaffinity(0, sizeof(selected), &selected) != 0)
        return std::string("pinning unavailable: sched_setaffinity: ") +
            strerror(errno);
    std::ostringstream result;
    result << "pinned to logical CPU " << cpu;
    return result.str();
#else
    return "pinning unavailable on this platform";
#endif
}

static Environment GetEnvironment(const std::string& affinity)
{
    Environment environment;
    environment.compiler = CompilerName();
    environment.build_type = BuildType();
    environment.build_flags = BuildFlags();
    environment.cpu = CPUName();
    environment.simd = SIMDName();
    environment.operating_system = OperatingSystemName();
    environment.affinity = affinity;
    environment.counter_backend = CounterBackendName();
    return environment;
}

static std::vector<unsigned> ShuffledIndices(unsigned count, uint64_t seed);

static unsigned NextPowerOfTwo(unsigned value)
{
    unsigned result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static bool IsPowerOfTwo(unsigned value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

static uint64_t TransformMemoryUnits(
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base)
{
    uint64_t units = 0;
    for (unsigned distance = 1; distance < count; distance <<= 1)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0; range < count_truncated; range += span)
        {
            // FFTSkewPadded is the sentinel at padded indices 1,2,4,...128.
            // The sentinel specialization reads x/y and writes y (3B), while
            // a general butterfly reads and writes both buffers (4B).
            const bool sentinel =
                IsPowerOfTwo(skew_base + range + distance);
            units += static_cast<uint64_t>(distance) * (sentinel ? 3 : 4);
        }
    }
    return units;
}

static uint64_t ErrorFFTMemoryUnits(
    unsigned count_truncated,
    unsigned count,
    unsigned skew_base,
    const std::vector<unsigned>& missing_locations,
    unsigned initial_distance)
{
    (void)count;
    uint64_t units = 0;
    for (unsigned distance = initial_distance;
        distance != 0;
        distance >>= 1)
    {
        const unsigned span = distance << 1;
        for (unsigned range = 0; range < count_truncated; range += span)
        {
            bool needed = false;
            for (size_t i = 0; i < missing_locations.size(); ++i)
            {
                if (missing_locations[i] >= range &&
                    missing_locations[i] < range + span)
                {
                    needed = true;
                    break;
                }
            }
            if (needed)
            {
                const bool sentinel =
                    IsPowerOfTwo(skew_base + range + distance);
                units += static_cast<uint64_t>(distance) *
                    (sentinel ? 3 : 4);
            }
        }
    }
    return units;
}

static uint64_t FusedDerivativeBoundaryMemoryUnits(unsigned count)
{
    const unsigned half_count = count >> 1;
    uint64_t units = 0;
    for (unsigned q = 0; q < half_count; ++q)
    {
        // The direct row reads A[q] and R[q], writes both outputs, and reads
        // one higher right-half source for every zero bit of q.  This models
        // memory operands folded into VPTERNLOG as reads too.
        units += 4;
        for (unsigned bit = 1; bit < half_count; bit <<= 1)
        {
            if ((q & bit) == 0)
                ++units;
        }
    }
    return units;
}

static uint64_t ModeledEncodePayloadBytes(
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    bool transpose)
{
    const unsigned m = NextPowerOfTwo(recovery_count);
    unsigned chunks = 0;
    unsigned zero_buffers = 0;
    uint64_t transform_units = 0;
    for (unsigned chunk_start = 0; chunk_start < original_count;
         chunk_start += m)
    {
        const unsigned chunk_count = std::min(m, original_count - chunk_start);
        ++chunks;
        zero_buffers += m - chunk_count;
        transform_units += TransformMemoryUnits(
            chunk_count, m, m + chunk_start);
    }
    transform_units += TransformMemoryUnits(recovery_count, m, 0);

    uint64_t units = 0;
    units += static_cast<uint64_t>(original_count) * 2; // copy
    units += zero_buffers; // zero store
    units += static_cast<uint64_t>(m) * (chunks - 1) * 3; // XOR load/load/store
    units += transform_units;
    if (transpose)
        units += static_cast<uint64_t>(original_count + recovery_count) * 2;
    return units * buffer_bytes;
}

static uint64_t ModeledDecodePayloadBytes(
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    unsigned loss_count,
    bool transpose)
{
    const unsigned m = NextPowerOfTwo(recovery_count);
    const unsigned n = NextPowerOfTwo(m + original_count);
    const uint64_t seed = UINT64_C(0xff8c000000000000) ^
        (static_cast<uint64_t>(original_count) << 40) ^
        (static_cast<uint64_t>(recovery_count) << 24) ^ buffer_bytes;
    const uint64_t loss_seed = seed ^
        (static_cast<uint64_t>(loss_count) << 8) ^ UINT64_C(0xd3c0de);
    const std::vector<unsigned> original_order =
        ShuffledIndices(original_count, loss_seed);
    std::vector<unsigned> missing_locations;
    for (unsigned i = 0; i < loss_count; ++i)
        missing_locations.push_back(m + original_order[i]);

    // The deterministic benchmark supplies one available recovery for every
    // missing original, so (k-loss) originals + loss recoveries == k inputs.
    const unsigned present_originals = original_count - loss_count;
    const unsigned available_recoveries = loss_count;
    const unsigned present_inputs = present_originals + available_recoveries;
    const unsigned zero_buffers = n - present_inputs;
    const unsigned half_count = n >> 1;
    uint64_t left_derivative_buffers = 0;
    for (unsigned index = 1; index < half_count; ++index)
        left_derivative_buffers += ((index ^ (index - 1)) + 1) >> 1;
    const uint64_t derivative_boundary_units =
        FusedDerivativeBoundaryMemoryUnits(n);
    const uint64_t transform_units =
        TransformMemoryUnits(m + original_count, n, 0) +
        ErrorFFTMemoryUnits(
            m + original_count,
            n,
            0,
            missing_locations,
            n >> 2);

    uint64_t units = 0;
    units += static_cast<uint64_t>(present_inputs + loss_count) * 2;
    units += zero_buffers;
    units += left_derivative_buffers * 3;
    units += derivative_boundary_units;
    units += transform_units;
    if (transpose)
    {
        // Present originals and recoveries are transposed into the native
        // layout; recovered originals are transposed out.
        units += static_cast<uint64_t>(original_count + loss_count) * 2;
    }
    return units * buffer_bytes;
}

static std::string ScheduleID(
    const char* operation,
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    unsigned loss_count)
{
    if (original_count == 0 || recovery_count == 0)
        return "";
    std::ostringstream id;
    id << operation << "-k" << original_count << "-r" << recovery_count
        << "-b" << buffer_bytes;
    if (strcmp(operation, "decode") == 0 ||
        strcmp(operation, "locator_shift_select") == 0)
        id << "-loss" << loss_count;
    return id.str();
}

static Result MakeResult(
    const Options& options,
    const char* backend,
    const char* operation,
    bool transpose,
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    unsigned loss_count,
    const Timing& timing,
    uint64_t input_bytes,
    uint64_t output_bytes,
    double ratio,
    const char* note = "",
    bool cache_coloring_applied = false)
{
    Result result;
    result.record = "benchmark";
    result.backend = backend;
    result.operation = operation;
    result.transpose_included = transpose;
    result.cache_coloring_applied = cache_coloring_applied;
    result.original_count = original_count;
    result.recovery_count = recovery_count;
    result.buffer_bytes = buffer_bytes;
    result.loss_count = loss_count;
    result.warmups = options.warmups;
    result.iterations = timing.sample_count != 0
        ? timing.sample_count : options.iterations;
    result.timing = timing;
    result.input_bytes = input_bytes;
    result.output_bytes = output_bytes;
    result.ratio_vs_packed = ratio;
    result.schedule_id = ScheduleID(operation, original_count, recovery_count,
        buffer_bytes, loss_count);
    result.measurement_order = options.abba && !transpose &&
        original_count != 0 && recovery_count != 0 &&
        (strcmp(backend, "packed_ff8") == 0 ||
         strcmp(backend, "ff8xor_native") == 0)
            ? "ABBA" : "sequential";
    if (strncmp(backend, "ff8xor", 6) == 0)
    {
        if (strcmp(operation, "encode") == 0)
        {
            result.modeled_payload_bytes = ModeledEncodePayloadBytes(
                original_count, recovery_count, buffer_bytes, transpose);
        }
        else if (strcmp(operation, "decode") == 0)
        {
            result.modeled_payload_bytes = ModeledDecodePayloadBytes(
                original_count, recovery_count, buffer_bytes, loss_count,
                transpose);
        }

        if (strcmp(operation, "encode") == 0 ||
            strcmp(operation, "decode") == 0)
        {
            result.materialization_statistics =
                leopard::ff8xor::GetLastMaterializationStatistics();
            result.four_buffer_statistics =
                leopard::ff8xor::GetLastFourBufferStatistics();
            result.modeled_payload_bytes_elided =
                result.materialization_statistics.EstimatedPayloadBytesElided +
                static_cast<int64_t>(result.four_buffer_statistics.
                    EstimatedPayloadBytesElided);
        }
        result.modeled_payload_bytes_adjusted =
            static_cast<int64_t>(result.modeled_payload_bytes) -
            result.modeled_payload_bytes_elided;
    }
    result.note = note;
    if (result.cache_coloring_applied)
    {
        if (!result.note.empty())
            result.note += "; ";
        result.note += "cache-colored transform buffers";
    }
    return result;
}

static void AddUnique(std::vector<unsigned>& values, unsigned value)
{
    if (std::find(values.begin(), values.end(), value) == values.end())
        values.push_back(value);
}

static std::vector<unsigned> LossCounts(unsigned recovery_count)
{
    std::vector<unsigned> losses;
    AddUnique(losses, 1);
    AddUnique(losses, std::min(4U, recovery_count));
    AddUnique(losses, std::max(1U, recovery_count / 2));
    AddUnique(losses, recovery_count);
    return losses;
}

static std::vector<unsigned> ShuffledIndices(unsigned count, uint64_t seed)
{
    std::vector<unsigned> indices(count);
    for (unsigned i = 0; i < count; ++i)
        indices[i] = i;
    Random random(seed);
    for (unsigned i = count; i > 1; --i)
    {
        const unsigned other = static_cast<unsigned>(random.Next() % i);
        std::swap(indices[i - 1], indices[other]);
    }
    return indices;
}

static bool CheckTransposeHelper()
{
    BufferSet packed;
    BufferSet plane;
    BufferSet roundtrip;
    if (!packed.Allocate(1, 64) || !plane.Allocate(1, 64) ||
        !roundtrip.Allocate(1, 64))
        return false;
    FillRandom(packed, 64, UINT64_C(0x746573745f747270));
    PackedToPlane(packed[0], plane[0], 64);
    PlaneToPacked(plane[0], roundtrip[0], 64);
    return Equal(packed[0], roundtrip[0], 64);
}

static bool CheckCacheColorHelper()
{
    using leopard_ff8xor_test::BufferSet;
    using leopard_ff8xor_test::SaltedTransformCacheColor;
    using leopard_ff8xor_test::TransformCacheColor;

    // Decode copies original i to work slot m + i.  Cover every valid FF8
    // padded recovery size and original count, including additions that carry
    // across transform-index bits.  A salt of 63, for example, collides at
    // m=64 and i=64; the selected salt must survive the complete domain.
    for (unsigned index = 0; index < 256; ++index)
    {
        if (TransformCacheColor(index) ==
            SaltedTransformCacheColor(
                index, leopard_ff8xor_test::kDecodeWorkCacheColorSalt))
        {
            std::cerr << "decode recovery/work cache-color collision: index="
                << index << '\n';
            return false;
        }
    }

    for (unsigned m = 2; m <= 128; m <<= 1)
    {
        for (unsigned original_count = 1;
             original_count + m <= 256;
             ++original_count)
        {
            for (unsigned index = 0; index < original_count; ++index)
            {
                if (TransformCacheColor(index) ==
                    SaltedTransformCacheColor(
                        m + index,
                        leopard_ff8xor_test::kDecodeWorkCacheColorSalt))
                {
                    std::cerr << "decode source/work cache-color collision: m="
                        << m << " k=" << original_count
                        << " index=" << index << '\n';
                    return false;
                }
                if (SaltedTransformCacheColor(
                        index,
                        leopard_ff8xor_test::kDecodeWorkCacheColorSalt) ==
                    SaltedTransformCacheColor(
                        m + index,
                        leopard_ff8xor_test::kDecodeWorkCacheColorSalt))
                {
                    std::cerr << "decode output/input work cache-color collision: m="
                        << m << " k=" << original_count
                        << " index=" << index << '\n';
                    return false;
                }
            }
        }
    }

    // Radix-2 partners differ in one transform-index bit.  Check every such
    // pair for every bounded prefix rather than only the power-of-two sizes
    // used by the codec, so a future truncated schedule cannot invalidate the
    // allocation invariant silently.
    for (unsigned count = 1; count <= 256; ++count)
    {
        for (unsigned distance = 1; distance < 256; distance <<= 1)
        {
            for (unsigned index = 0; index < count; ++index)
            {
                const unsigned partner = index ^ distance;
                if (partner < count &&
                    TransformCacheColor(index) == TransformCacheColor(partner))
                {
                    std::cerr << "cache-color collision: n=" << count
                        << " distance=" << distance
                        << " index=" << index
                        << " partner=" << partner << '\n';
                    return false;
                }
            }
        }
    }

    bool seen[leopard_ff8xor_test::kCacheColorCount] = {};
    for (unsigned index = 0; index < 64; ++index)
        seen[TransformCacheColor(index)] = true;
    for (unsigned color = 0;
         color < leopard_ff8xor_test::kCacheColorCount;
         ++color)
    {
        if (!seen[color])
        {
            std::cerr << "cache-color map does not cover color " << color << '\n';
            return false;
        }
    }

    static const size_t kSizes[] = { 64, 1024, 64 * 1024 };
    static const unsigned kSalts[] = {
        0, leopard_ff8xor_test::kDecodeWorkCacheColorSalt
    };
    for (size_t size_index = 0;
         size_index < sizeof(kSizes) / sizeof(kSizes[0]);
         ++size_index)
    {
        for (size_t salt_index = 0;
             salt_index < sizeof(kSalts) / sizeof(kSalts[0]);
             ++salt_index)
        {
            BufferSet buffers;
            if (!buffers.Allocate(256, kSizes[size_index], true, 0,
                    kSalts[salt_index]))
            {
                std::cerr << "cache-color allocation self-test failed\n";
                return false;
            }
            for (unsigned index = 0; index < 256; ++index)
            {
                if (!buffers.ValidateColoredAllocation(
                        index, kSizes[size_index], index, kSalts[salt_index]))
                {
                    std::cerr << "invalid cache-colored allocation: index="
                        << index << " bytes=" << kSizes[size_index]
                        << " salt=" << kSalts[salt_index] << '\n';
                    return false;
                }
                const unsigned expected_color = SaltedTransformCacheColor(
                    index, kSalts[salt_index]);
                const uintptr_t address =
                    reinterpret_cast<uintptr_t>(buffers[index]);
                if ((address / leopard_ff8xor_test::kCacheLineBytes) %
                        leopard_ff8xor_test::kCacheColorCount != expected_color)
                    return false;
                if (leopard_ff8xor_test::PlaneStartsFullyCacheAlias(
                        kSizes[size_index]))
                {
                    const size_t plane_bytes = kSizes[size_index] / 8;
                    for (unsigned plane = 0; plane < 8; ++plane)
                    {
                        if ((address + plane * plane_bytes) %
                                leopard_ff8xor_test::kCacheColorPeriod !=
                            expected_color *
                                leopard_ff8xor_test::kCacheLineBytes)
                        {
                            std::cerr << "plane-start color mismatch: index="
                                << index << " plane=" << plane << '\n';
                            return false;
                        }
                    }
                }
                buffers[index][0] = static_cast<uint8_t>(index);
                buffers[index][kSizes[size_index] - 1] =
                    static_cast<uint8_t>(index ^ 0xff);
            }
        }
    }
    return true;
}

static bool VerifyPackedRecovery(
    const BufferSet& packed_recovery,
    const BufferSet& plane_work,
    unsigned recovery_count,
    uint64_t buffer_bytes,
    uint8_t* scratch)
{
    for (unsigned i = 0; i < recovery_count; ++i)
    {
        PlaneToPacked(plane_work[i], scratch, buffer_bytes);
        if (!Equal(packed_recovery[i], scratch, buffer_bytes))
        {
            std::cerr << "encode equivalence mismatch in recovery shard "
                << i << '\n';
            return false;
        }
    }
    return true;
}

static bool RunParameter(
    const Options& options,
    Reporter& reporter,
    unsigned original_count,
    unsigned recovery_count,
    uint64_t buffer_bytes)
{
    const bool use_cache_color = options.cache_color &&
        leopard_ff8xor_test::PlaneStartsFullyCacheAlias(
            static_cast<size_t>(buffer_bytes));
    const unsigned packed_encode_count =
        leo_encode_work_count(original_count, recovery_count);
    const unsigned xor_encode_count =
        leo_ff8xor_encode_work_count(original_count, recovery_count);
    if (packed_encode_count == 0 || xor_encode_count != packed_encode_count)
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "invalid or inconsistent encode work count");
        return true;
    }

    BufferSet packed_original;
    BufferSet plane_original;
    BufferSet scratch;
    if (!packed_original.Allocate(original_count, static_cast<size_t>(buffer_bytes)) ||
        !plane_original.Allocate(original_count, static_cast<size_t>(buffer_bytes),
            use_cache_color, 0, 0) ||
        !scratch.Allocate(1, static_cast<size_t>(buffer_bytes)))
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "allocation failed for originals");
        return true;
    }

    const uint64_t seed = UINT64_C(0xff8c000000000000) ^
        (static_cast<uint64_t>(original_count) << 40) ^
        (static_cast<uint64_t>(recovery_count) << 24) ^ buffer_bytes;
    FillRandom(packed_original, buffer_bytes, seed);
    for (unsigned i = 0; i < original_count; ++i)
        PackedToPlane(packed_original[i], plane_original[i], buffer_bytes);

    std::vector<const void*> packed_original_ptrs(original_count);
    std::vector<const void*> plane_original_ptrs(original_count);
    for (unsigned i = 0; i < original_count; ++i)
    {
        packed_original_ptrs[i] = packed_original[i];
        plane_original_ptrs[i] = plane_original[i];
    }

    LeopardResult api_result = Leopard_Success;
    Timing packed_encode_timing;
    Timing xor_encode_timing;
    BufferSet packed_encode_work;
    BufferSet xor_encode_work;
    if (!packed_encode_work.Allocate(
            packed_encode_count, static_cast<size_t>(buffer_bytes)) ||
        !xor_encode_work.Allocate(
            xor_encode_count, static_cast<size_t>(buffer_bytes),
            use_cache_color, 0, 0))
    {
        reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
            "allocation failed for encode work");
        return true;
    }
    const auto packed_encode = [&]() -> LeopardResult {
            return leo_encode(buffer_bytes, original_count, recovery_count,
                packed_encode_count, &packed_original_ptrs[0],
                packed_encode_work.Data());
        };
    const auto xor_encode = [&]() -> LeopardResult {
            return leo_ff8xor_encode(buffer_bytes, original_count, recovery_count,
                xor_encode_count, &plane_original_ptrs[0],
                xor_encode_work.Data());
        };
    const bool encode_measured = options.abba
        ? MeasurePairABBA(options, packed_encode, xor_encode,
            packed_encode_timing, xor_encode_timing, api_result)
        : (Measure(options, packed_encode, packed_encode_timing, api_result) &&
           Measure(options, xor_encode, xor_encode_timing, api_result));
    if (!encode_measured)
    {
        if (api_result == Leopard_TooMuchData)
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
                "backend reports unsupported parameters");
            return true;
        }
        std::cerr << "encode failed: " << leo_result_string(api_result) << '\n';
        return false;
    }

    const uint64_t encode_input_bytes = buffer_bytes * original_count;
    const uint64_t encode_output_bytes = buffer_bytes * recovery_count;
    reporter.Print(MakeResult(options, "packed_ff8", "encode", false,
        original_count, recovery_count, buffer_bytes, 0,
        packed_encode_timing, encode_input_bytes, encode_output_bytes, 1.0));

    BufferSet packed_recovery;
    if (!packed_encode_work.TransferFirst(recovery_count, packed_recovery))
    {
        std::cerr << "unable to retain packed recovery buffers\n";
        return false;
    }

    if (!VerifyPackedRecovery(packed_recovery, xor_encode_work,
            recovery_count, buffer_bytes, scratch[0]))
        return false;

    reporter.Print(MakeResult(options, "ff8xor_native", "encode", false,
        original_count, recovery_count, buffer_bytes, 0,
        xor_encode_timing, encode_input_bytes, encode_output_bytes,
        packed_encode_timing.median_usec / xor_encode_timing.median_usec,
        "", use_cache_color));

    if (options.include_transpose)
    {
        BufferSet packed_output;
        if (!packed_output.Allocate(
                recovery_count, static_cast<size_t>(buffer_bytes)))
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, 0,
                "allocation failed for transpose-inclusive encode output");
        }
        else
        {
            Timing included_timing;
            if (!Measure(options,
                [&]() -> LeopardResult {
                    for (unsigned i = 0; i < original_count; ++i)
                        PackedToPlane(packed_original[i], plane_original[i], buffer_bytes);
                    LeopardResult result = leo_ff8xor_encode(
                        buffer_bytes, original_count, recovery_count,
                        xor_encode_count, &plane_original_ptrs[0],
                        xor_encode_work.Data());
                    if (result != Leopard_Success)
                        return result;
                    for (unsigned i = 0; i < recovery_count; ++i)
                        PlaneToPacked(xor_encode_work[i], packed_output[i], buffer_bytes);
                    return Leopard_Success;
                }, included_timing, api_result))
            {
                std::cerr << "transpose-inclusive ff8xor encode failed: "
                    << leo_result_string(api_result) << '\n';
                return false;
            }
            for (unsigned i = 0; i < recovery_count; ++i)
            {
                if (!Equal(packed_recovery[i], packed_output[i], buffer_bytes))
                {
                    std::cerr << "transpose-inclusive encode mismatch at shard "
                        << i << '\n';
                    return false;
                }
            }
            reporter.Print(MakeResult(options, "ff8xor_packed_boundary", "encode", true,
                original_count, recovery_count, buffer_bytes, 0,
                included_timing, encode_input_bytes, encode_output_bytes,
                packed_encode_timing.median_usec / included_timing.median_usec,
                options.portable_transpose ?
                    "portable boundary transpose control" :
                    AutoBoundaryTransposeNote(),
                use_cache_color));
        }
    }

    BufferSet plane_recovery;
    if (!xor_encode_work.TransferFirst(recovery_count, plane_recovery))
    {
        std::cerr << "unable to retain ff8xor recovery buffers\n";
        return false;
    }

    const std::vector<unsigned> loss_counts = LossCounts(recovery_count);
    for (size_t loss_case = 0; loss_case < loss_counts.size(); ++loss_case)
    {
        const unsigned loss_count = loss_counts[loss_case];
        const uint64_t loss_seed = seed ^
            (static_cast<uint64_t>(loss_count) << 8) ^ UINT64_C(0xd3c0de);
        const std::vector<unsigned> original_order =
            ShuffledIndices(original_count, loss_seed);
        const std::vector<unsigned> recovery_order =
            ShuffledIndices(recovery_count, loss_seed ^ UINT64_C(0x7265636f76657279));

        std::vector<bool> original_missing(original_count, false);
        std::vector<bool> recovery_available(recovery_count, false);
        for (unsigned i = 0; i < loss_count; ++i)
        {
            original_missing[original_order[i]] = true;
            recovery_available[recovery_order[i]] = true;
        }

        std::vector<const void*> packed_decode_original(original_count);
        std::vector<const void*> plane_decode_original(original_count);
        std::vector<const void*> packed_decode_recovery(recovery_count);
        std::vector<const void*> plane_decode_recovery(recovery_count);
        for (unsigned i = 0; i < original_count; ++i)
        {
            packed_decode_original[i] = original_missing[i] ? NULL : packed_original[i];
            plane_decode_original[i] = original_missing[i] ? NULL : plane_original[i];
        }
        for (unsigned i = 0; i < recovery_count; ++i)
        {
            packed_decode_recovery[i] = recovery_available[i]
                ? packed_recovery[i] : NULL;
            plane_decode_recovery[i] = recovery_available[i]
                ? plane_recovery[i] : NULL;
        }

        const unsigned packed_decode_count =
            leo_decode_work_count(original_count, recovery_count);
        const unsigned xor_decode_count =
            leo_ff8xor_decode_work_count(original_count, recovery_count);
        if (packed_decode_count == 0 || xor_decode_count != packed_decode_count)
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                "invalid or inconsistent decode work count");
            continue;
        }

        BufferSet packed_decode_work;
        BufferSet xor_decode_work;
        if (!packed_decode_work.Allocate(
                packed_decode_count, static_cast<size_t>(buffer_bytes)) ||
            !xor_decode_work.Allocate(
                xor_decode_count, static_cast<size_t>(buffer_bytes),
                use_cache_color, 0,
                leopard_ff8xor_test::kDecodeWorkCacheColorSalt))
        {
            reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                "allocation failed for decode work");
            continue;
        }
        Timing packed_decode_timing;
        Timing xor_decode_timing;
        const auto packed_decode = [&]() -> LeopardResult {
                return leo_decode(buffer_bytes, original_count, recovery_count,
                    packed_decode_count, &packed_decode_original[0],
                    &packed_decode_recovery[0], packed_decode_work.Data());
            };
        const auto xor_decode = [&]() -> LeopardResult {
                return leo_ff8xor_decode(buffer_bytes, original_count, recovery_count,
                    xor_decode_count, &plane_decode_original[0],
                    &plane_decode_recovery[0], xor_decode_work.Data());
            };
        const bool decode_measured = options.abba
            ? MeasurePairABBA(options, packed_decode, xor_decode,
                packed_decode_timing, xor_decode_timing, api_result)
            : (Measure(options, packed_decode, packed_decode_timing, api_result) &&
               Measure(options, xor_decode, xor_decode_timing, api_result));
        if (!decode_measured)
        {
            std::cerr << "decode failed: " << leo_result_string(api_result) << '\n';
            return false;
        }
        for (unsigned i = 0; i < loss_count; ++i)
        {
            const unsigned index = original_order[i];
            if (!Equal(packed_original[index], packed_decode_work[index], buffer_bytes))
            {
                std::cerr << "packed decode mismatch at original shard " << index << '\n';
                return false;
            }
        }

        const uint64_t decode_input_bytes = buffer_bytes * original_count;
        const uint64_t decode_output_bytes = buffer_bytes * loss_count;
        reporter.Print(MakeResult(options, "packed_ff8", "decode", false,
            original_count, recovery_count, buffer_bytes, loss_count,
            packed_decode_timing, decode_input_bytes, decode_output_bytes, 1.0));
        for (unsigned i = 0; i < loss_count; ++i)
        {
            const unsigned index = original_order[i];
            if (!Equal(plane_original[index], xor_decode_work[index], buffer_bytes))
            {
                std::cerr << "native ff8xor decode mismatch at original shard "
                    << index << '\n';
                return false;
            }
            PlaneToPacked(xor_decode_work[index], scratch[0], buffer_bytes);
            if (!Equal(packed_original[index], scratch[0], buffer_bytes))
            {
                std::cerr << "packed ff8xor decode equivalence mismatch at shard "
                    << index << '\n';
                return false;
            }
        }
        Result xor_decode_result = MakeResult(
            options, "ff8xor_native", "decode", false,
            original_count, recovery_count, buffer_bytes, loss_count,
            xor_decode_timing, decode_input_bytes, decode_output_bytes,
            packed_decode_timing.median_usec / xor_decode_timing.median_usec,
            "", use_cache_color);
        xor_decode_result.locator_shift =
            leopard::ff8xor::GetLastLocatorShiftForTesting();
        reporter.Print(xor_decode_result);

        if (options.include_transpose)
        {
            BufferSet packed_output;
            if (!packed_output.Allocate(
                    loss_count, static_cast<size_t>(buffer_bytes)))
            {
                reporter.Skip(original_count, recovery_count, buffer_bytes, loss_count,
                    "allocation failed for transpose-inclusive decode output");
            }
            else
            {
                Timing included_timing;
                if (!Measure(options,
                    [&]() -> LeopardResult {
                        for (unsigned i = 0; i < original_count; ++i)
                            if (!original_missing[i])
                                PackedToPlane(packed_original[i], plane_original[i],
                                    buffer_bytes);
                        for (unsigned i = 0; i < recovery_count; ++i)
                            if (recovery_available[i])
                                PackedToPlane(packed_recovery[i], plane_recovery[i],
                                    buffer_bytes);
                        LeopardResult result = leo_ff8xor_decode(
                            buffer_bytes, original_count, recovery_count,
                            xor_decode_count, &plane_decode_original[0],
                            &plane_decode_recovery[0], xor_decode_work.Data());
                        if (result != Leopard_Success)
                            return result;
                        for (unsigned i = 0; i < loss_count; ++i)
                        {
                            const unsigned index = original_order[i];
                            PlaneToPacked(xor_decode_work[index], packed_output[i],
                                buffer_bytes);
                        }
                        return Leopard_Success;
                    }, included_timing, api_result))
                {
                    std::cerr << "transpose-inclusive ff8xor decode failed: "
                        << leo_result_string(api_result) << '\n';
                    return false;
                }
                for (unsigned i = 0; i < loss_count; ++i)
                {
                    const unsigned index = original_order[i];
                    if (!Equal(packed_original[index], packed_output[i], buffer_bytes))
                    {
                        std::cerr << "transpose-inclusive decode mismatch at shard "
                            << index << '\n';
                        return false;
                    }
                }
                Result included_result = MakeResult(
                    options, "ff8xor_packed_boundary", "decode", true,
                    original_count, recovery_count, buffer_bytes, loss_count,
                    included_timing, decode_input_bytes, decode_output_bytes,
                    packed_decode_timing.median_usec / included_timing.median_usec,
                    options.portable_transpose ?
                        "portable boundary transpose control" :
                        AutoBoundaryTransposeNote(),
                    use_cache_color);
                included_result.locator_shift =
                    leopard::ff8xor::GetLastLocatorShiftForTesting();
                reporter.Print(included_result);
            }
        }
    }

    if (!options.csv && !options.json)
        std::cout << '\n';
    return true;
}

static void PrintMicro(
    const Options& options,
    Reporter& reporter,
    const char* backend,
    const char* operation,
    const std::string& note,
    const Timing& timing,
    uint64_t buffer_bytes,
    uint64_t input_bytes,
    uint64_t output_bytes,
    uint64_t modeled_payload_bytes,
    const char* measurement_order = "sequential",
    unsigned original_count = 0,
    unsigned recovery_count = 0,
    unsigned loss_count = 0,
    unsigned locator_shift = std::numeric_limits<unsigned>::max())
{
    Result result = MakeResult(options, backend, operation, false,
        original_count, recovery_count, buffer_bytes, loss_count,
        timing, input_bytes, output_bytes, 0, note.c_str());
    result.record = "microbenchmark";
    result.modeled_payload_bytes = modeled_payload_bytes;
    result.modeled_payload_bytes_adjusted =
        static_cast<int64_t>(modeled_payload_bytes);
    result.measurement_order = measurement_order;
    result.locator_shift = locator_shift;
    reporter.Print(result);
}

static bool RunLocatorSelectorMicrobenchmarks(
    const Options& options,
    Reporter& reporter)
{
    struct SelectorCase
    {
        unsigned OriginalCount;
        unsigned RecoveryCount;
        unsigned LossCount;
    };
    static const SelectorCase kCases[] = {
        { 8, 2, 1 },
        { 16, 4, 4 },
        { 32, 8, 4 },
        { 64, 16, 8 },
        { 128, 32, 16 },
        { 128, 128, 128 }
    };

    for (size_t case_index = 0;
         case_index < sizeof(kCases) / sizeof(kCases[0]);
         ++case_index)
    {
        const SelectorCase parameters = kCases[case_index];
        const unsigned term_count =
            parameters.OriginalCount + parameters.LossCount;
        leopard::ff8xor::ffe_t logarithms[leopard::ff8xor::kOrder];
        bool inverse[leopard::ff8xor::kOrder];
        uint32_t random = UINT32_C(0x51ec7001) ^
            (parameters.OriginalCount << 16) ^
            (parameters.RecoveryCount << 8) ^ parameters.LossCount;
        for (unsigned index = 0; index < term_count; ++index)
        {
            random ^= random << 13;
            random ^= random >> 17;
            random ^= random << 5;
            logarithms[index] = static_cast<leopard::ff8xor::ffe_t>(random);
            inverse[index] = index >= parameters.OriginalCount;
        }
        // Exercise the redundant multiplier-one logarithm in every timed case.
        logarithms[0] = leopard::ff8xor::kModulus;

        const unsigned reference_shift =
            leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
                logarithms, inverse, term_count);
        const unsigned rotated_shift =
            leopard::ff8xor::SelectLocatorShiftForTesting(
                logarithms, inverse, term_count);
        if (reference_shift != rotated_shift ||
            reference_shift >= leopard::ff8xor::kModulus)
        {
            std::cerr << "locator selector microbenchmark mismatch: k="
                << parameters.OriginalCount
                << " r=" << parameters.RecoveryCount
                << " loss=" << parameters.LossCount
                << " reference=" << reference_shift
                << " rotated=" << rotated_shift << '\n';
            return false;
        }

        typedef unsigned (*SelectorFunction)(
            const leopard::ff8xor::ffe_t*, const bool*, unsigned);
        SelectorFunction volatile reference_function =
            &leopard::ff8xor::SelectLocatorShiftReferenceForTesting;
        SelectorFunction volatile rotated_function =
            &leopard::ff8xor::SelectLocatorShiftForTesting;
        const auto reference = [&]() -> LeopardResult {
            BenchmarkSink ^=
                reference_function(logarithms, inverse, term_count);
            return Leopard_Success;
        };
        const auto rotated = [&]() -> LeopardResult {
            BenchmarkSink ^=
                rotated_function(logarithms, inverse, term_count);
            return Leopard_Success;
        };
        Timing reference_timing;
        Timing rotated_timing;
        LeopardResult result = Leopard_Success;
        if (!MeasurePairABBA(
                options,
                reference,
                rotated,
                reference_timing,
                rotated_timing,
                result))
        {
            return false;
        }

        std::ostringstream reference_note;
        reference_note << "old shift-major exact gate/depth scoring; terms="
            << term_count << "; selected_shift=" << reference_shift
            << "; paired ABBA control";
        std::ostringstream rotated_note;
        rotated_note << "generated rotated gate rows; depth only on gate ties; "
            << "terms=" << term_count
            << "; selected_shift=" << rotated_shift
            << "; paired ABBA";
        PrintMicro(
            options,
            reporter,
            "ff8xor_selector_reference",
            "locator_shift_select",
            reference_note.str(),
            reference_timing,
            0,
            0,
            0,
            0,
            "ABBA",
            parameters.OriginalCount,
            parameters.RecoveryCount,
            parameters.LossCount,
            reference_shift);
        PrintMicro(
            options,
            reporter,
            "ff8xor_selector_rotated",
            "locator_shift_select",
            rotated_note.str(),
            rotated_timing,
            0,
            0,
            0,
            0,
            "ABBA",
            parameters.OriginalCount,
            parameters.RecoveryCount,
            parameters.LossCount,
            rotated_shift);
    }
    return true;
}

struct CircuitCase
{
    unsigned coefficient;
    unsigned gates;
    unsigned depth;
    const char* description;
};

typedef leopard::ff8xor::CircuitCost (*CircuitCostGetter)(
    leopard::ff8xor::ffe_t coefficient);

static CircuitCase SelectCircuitCase(
    CircuitCostGetter get_cost,
    unsigned begin,
    unsigned end,
    double target_gate_count,
    bool select_minimum,
    bool select_maximum,
    const char* description)
{
    CircuitCase selected = { begin, 0, 0, description };
    double selected_distance = std::numeric_limits<double>::max();
    bool have_selected = false;

    for (unsigned coefficient = begin; coefficient < end; ++coefficient)
    {
        const leopard::ff8xor::CircuitCost cost = get_cost(
            static_cast<leopard::ff8xor::ffe_t>(coefficient));
        const double distance = cost.GateCount > target_gate_count ?
            cost.GateCount - target_gate_count :
            target_gate_count - cost.GateCount;
        const bool better = !have_selected ||
            (select_minimum &&
                (cost.GateCount < selected.gates ||
                 (cost.GateCount == selected.gates &&
                  (cost.Depth < selected.depth ||
                   (cost.Depth == selected.depth &&
                    coefficient < selected.coefficient))))) ||
            (select_maximum &&
                (cost.GateCount > selected.gates ||
                 (cost.GateCount == selected.gates &&
                  (cost.Depth > selected.depth ||
                   (cost.Depth == selected.depth &&
                    coefficient < selected.coefficient))))) ||
            (!select_minimum && !select_maximum &&
                (distance < selected_distance ||
                 (distance == selected_distance &&
                  (cost.Depth < selected.depth ||
                   (cost.Depth == selected.depth &&
                    coefficient < selected.coefficient)))));
        if (better)
        {
            selected.coefficient = coefficient;
            selected.gates = cost.GateCount;
            selected.depth = cost.Depth;
            selected_distance = distance;
            have_selected = true;
        }
    }
    return selected;
}

static double AverageCircuitGateCount(
    CircuitCostGetter get_cost,
    unsigned begin,
    unsigned end)
{
    uint64_t sum = 0;
    for (unsigned coefficient = begin; coefficient < end; ++coefficient)
    {
        sum += get_cost(static_cast<leopard::ff8xor::ffe_t>(coefficient))
            .GateCount;
    }
    return static_cast<double>(sum) / (end - begin);
}

static std::string CircuitCaseNote(
    const char* coefficient_name,
    const CircuitCase& circuit,
    const char* suffix)
{
    std::ostringstream note;
    note << coefficient_name << '=' << circuit.coefficient
        << "; gates=" << circuit.gates
        << "; depth=" << circuit.depth
        << "; " << circuit.description
        << "; " << suffix;
    return note.str();
}

static bool RunMicrobenchmarks(const Options& base_options, Reporter& reporter)
{
    Options options = base_options;
    options.warmups = base_options.quick ? 2 : 3;
    options.iterations = base_options.quick ? 15 : 31;
    const uint64_t bytes = base_options.quick ? 64 * 1024 : 1024 * 1024;

    if (!RunLocatorSelectorMicrobenchmarks(options, reporter))
        return false;

    BufferSet a;
    BufferSet b;
    BufferSet c;
    BufferSet d;
    BufferSet e;
    BufferSet xor_destination;
    BufferSet xor_source;
    BufferSet xor_initial;
    if (!a.Allocate(1, static_cast<size_t>(bytes)) ||
        !b.Allocate(1, static_cast<size_t>(bytes)) ||
        !c.Allocate(1, static_cast<size_t>(bytes)) ||
        !d.Allocate(1, static_cast<size_t>(bytes)) ||
        !e.Allocate(1, static_cast<size_t>(bytes)) ||
        !xor_destination.Allocate(4, static_cast<size_t>(bytes)) ||
        !xor_source.Allocate(4, static_cast<size_t>(bytes)) ||
        !xor_initial.Allocate(4, static_cast<size_t>(bytes)))
    {
        reporter.Skip(0, 0, bytes, 0, "allocation failed for microbenchmarks");
        return true;
    }

    FillRandom(
        xor_destination, bytes, UINT64_C(0x786f725f64657374));
    FillRandom(xor_source, bytes, UINT64_C(0x786f725f736f7572));
    for (unsigned stream = 0; stream < 4; ++stream)
    {
        memcpy(
            xor_initial[stream], xor_destination[stream],
            static_cast<size_t>(bytes));
    }

    LeopardResult result = Leopard_Success;
    const uint64_t xor_sizes[] = {
        64, 1024, 4096, 64 * 1024, 1024 * 1024
    };
    const size_t xor_size_count = base_options.quick ? 4 : 5;
    const unsigned xor_counts[] = { 2, 4 };
    for (size_t count_index = 0;
         count_index < sizeof(xor_counts) / sizeof(xor_counts[0]);
         ++count_index)
    {
        const unsigned xor_count = xor_counts[count_index];
        for (size_t size_index = 0; size_index < xor_size_count; ++size_index)
        {
            const uint64_t xor_bytes = xor_sizes[size_index];
            for (unsigned stream = 0; stream < xor_count; ++stream)
            {
                memcpy(
                    xor_destination[stream], xor_initial[stream],
                    static_cast<size_t>(xor_bytes));
            }

            // Validate every measured size/count outside timing.  The paired
            // A-B-B-A loop then uses the same buffers: each XOR toggles them,
            // and each complete round returns to its starting contents.
            leopard::ff8xor::XorBuffers(
                xor_bytes, xor_count,
                xor_destination.Data(), xor_source.Data());
            for (unsigned stream = 0; stream < xor_count; ++stream)
            {
                for (uint64_t offset = 0; offset < xor_bytes; ++offset)
                {
                    if (xor_destination[stream][offset] !=
                        static_cast<uint8_t>(
                            xor_initial[stream][offset] ^
                            xor_source[stream][offset]))
                    {
                        std::cerr
                            << "microbenchmark batched XOR mismatch: count="
                            << xor_count << " bytes=" << xor_bytes
                            << " stream=" << stream
                            << " offset=" << offset << '\n';
                        return false;
                    }
                }
                memcpy(
                    xor_destination[stream], xor_initial[stream],
                    static_cast<size_t>(xor_bytes));
            }

            const auto sequential_xor = [&]() -> LeopardResult {
                    for (unsigned stream = 0; stream < xor_count; ++stream)
                    {
                        leopard::ff8xor::XorBuffer(
                            xor_bytes,
                            xor_destination[stream],
                            xor_source[stream]);
                    }
                    BenchmarkSink ^= xor_destination[0][0];
                    return Leopard_Success;
                };
            const auto batched_xor = [&]() -> LeopardResult {
                    leopard::ff8xor::XorBuffers(
                        xor_bytes, xor_count,
                        xor_destination.Data(), xor_source.Data());
                    BenchmarkSink ^= xor_destination[0][0];
                    return Leopard_Success;
                };
            Timing sequential_timing;
            Timing batched_timing;
            if (!MeasurePairABBA(
                    options,
                    sequential_xor,
                    batched_xor,
                    sequential_timing,
                    batched_timing,
                    result))
                return false;

            const char* sequential_operation = xor_count == 2 ?
                "xor2_sequential" : "xor4_sequential";
            const char* batched_operation = xor_count == 2 ?
                "xor2_batched" : "xor4_batched";
            std::ostringstream sequential_note;
            sequential_note << xor_count
                << " independent whole-buffer XOR calls; paired ABBA control";
            std::ostringstream batched_note;
            batched_note << "one SIMD-mode resolution; ";
            if (xor_bytes > 1024)
                batched_note << "tuned sequential vector loops";
            else
                batched_note << xor_count << "-stream vector loop";
            batched_note << "; paired ABBA";
            PrintMicro(
                options, reporter, "ff8xor_native", sequential_operation,
                sequential_note.str(), sequential_timing, xor_bytes,
                xor_bytes * xor_count * 2, xor_bytes * xor_count,
                xor_bytes * xor_count * 3, "ABBA");
            PrintMicro(
                options, reporter, "ff8xor_native", batched_operation,
                batched_note.str(), batched_timing, xor_bytes,
                xor_bytes * xor_count * 2, xor_bytes * xor_count,
                xor_bytes * xor_count * 3, "ABBA");
        }
    }

    // Cover both logarithmic spellings of identity, the least expensive
    // non-identity circuit, a near-average circuit, and a maximum circuit.
    // This guards both logarithmic identity spellings and their out-of-place
    // copy fast path, while exposing
    // coefficient sensitivity rather than a favorable hand-picked case.
    const leopard::ff8xor::CircuitCost multiply_identity =
        leopard::ff8xor::GetMultiplyCircuitCost(0);
    const leopard::ff8xor::CircuitCost multiply_identity_255 =
        leopard::ff8xor::GetMultiplyCircuitCost(255);
    const double multiply_nonidentity_average = AverageCircuitGateCount(
        &leopard::ff8xor::GetMultiplyCircuitCost, 1, 255);
    const CircuitCase multiply_cases[] = {
        { 0, multiply_identity.GateCount, multiply_identity.Depth,
            "identity multiplier (canonical logarithm)" },
        SelectCircuitCase(&leopard::ff8xor::GetMultiplyCircuitCost,
            1, 255, 0., true, false,
            "minimum-gate non-identity multiplier"),
        SelectCircuitCase(&leopard::ff8xor::GetMultiplyCircuitCost,
            1, 255, multiply_nonidentity_average, false, false,
            "near-average-gate multiplier"),
        SelectCircuitCase(&leopard::ff8xor::GetMultiplyCircuitCost,
            1, 255, 0., false, true,
            "maximum-gate multiplier"),
        { 255, multiply_identity_255.GateCount,
            multiply_identity_255.Depth,
            "identity multiplier (redundant logarithm 255)" }
    };

    FillRandom(b, bytes, UINT64_C(0x6d6963726f5f6232));
    PackedToPlane(b[0], d[0], bytes);

    for (size_t case_index = 0;
         case_index < sizeof(multiply_cases) / sizeof(multiply_cases[0]);
         ++case_index)
    {
        const CircuitCase circuit = multiply_cases[case_index];

        // Verify packed and plane multiplication equivalence once, outside the
        // timed regions.  A zero packed destination turns multiply-add into a
        // multiplication-only reference result.
        memset(a[0], 0, static_cast<size_t>(bytes));
        leopard::ff8::ExperimentalPackedMulAdd(
            bytes, a[0], b[0], static_cast<uint8_t>(circuit.coefficient));
        leopard::ff8xor::MultiplyBuffer(
            bytes, e[0], d[0], static_cast<uint8_t>(circuit.coefficient));
        PlaneToPacked(e[0], c[0], bytes);
        if (!Equal(a[0], c[0], bytes))
        {
            std::cerr << "microbenchmark multiply equivalence failed for log="
                << circuit.coefficient << '\n';
            return false;
        }

        FillRandom(a, bytes, UINT64_C(0x6d6963726f5f7800) ^
            circuit.coefficient);
        Timing packed_timing;
        if (!Measure(options,
            [&]() -> LeopardResult {
                leopard::ff8::ExperimentalPackedMulAdd(
                    bytes, a[0], b[0],
                    static_cast<uint8_t>(circuit.coefficient));
                BenchmarkSink ^= a[0][0];
                return Leopard_Success;
            }, packed_timing, result))
            return false;
        PrintMicro(options, reporter, "packed_ff8", "multiply_add",
            CircuitCaseNote("log", circuit,
                "existing packed lookup multiply-add; gate metadata is for ff8xor"),
            packed_timing, bytes, bytes * 2, bytes, bytes * 3);

        Timing xor_timing;
        if (!Measure(options,
            [&]() -> LeopardResult {
                leopard::ff8xor::MultiplyBuffer(
                    bytes, e[0], d[0],
                    static_cast<uint8_t>(circuit.coefficient));
                BenchmarkSink ^= e[0][0];
                return Leopard_Success;
            }, xor_timing, result))
            return false;
        PrintMicro(options, reporter, "ff8xor_native", "multiply",
            CircuitCaseNote("log", circuit,
                circuit.coefficient == 0 || circuit.coefficient == 255 ?
                    "out-of-place identity memmove fast path" :
                    "generated multiplication-only whole-buffer kernel"),
            xor_timing, bytes, bytes, bytes, bytes * 2);
    }

    const leopard::ff8xor::CircuitCost sentinel_cost =
        leopard::ff8xor::GetFFTCircuitCost(255);
    const double fft_nonsentinel_average = AverageCircuitGateCount(
        &leopard::ff8xor::GetFFTCircuitCost, 0, 255);
    const CircuitCase butterfly_cases[] = {
        { 255, sentinel_cost.GateCount, sentinel_cost.Depth,
            "sentinel y^=x fast path" },
        SelectCircuitCase(&leopard::ff8xor::GetFFTCircuitCost,
            0, 255, 0., true, false,
            "minimum-gate non-sentinel butterfly"),
        SelectCircuitCase(&leopard::ff8xor::GetFFTCircuitCost,
            0, 255, fft_nonsentinel_average, false, false,
            "near-average-gate butterfly"),
        SelectCircuitCase(&leopard::ff8xor::GetFFTCircuitCost,
            0, 255, 0., false, true,
            "maximum-gate butterfly")
    };
    const leopard::ff8xor::CircuitCost inverse_sentinel_cost =
        leopard::ff8xor::GetIFFTCircuitCost(255);
    const double ifft_nonsentinel_average = AverageCircuitGateCount(
        &leopard::ff8xor::GetIFFTCircuitCost, 0, 255);
    const CircuitCase inverse_butterfly_cases[] = {
        { 255, inverse_sentinel_cost.GateCount,
            inverse_sentinel_cost.Depth, "sentinel y^=x fast path" },
        SelectCircuitCase(&leopard::ff8xor::GetIFFTCircuitCost,
            0, 255, 0., true, false,
            "minimum-gate non-sentinel inverse butterfly"),
        SelectCircuitCase(&leopard::ff8xor::GetIFFTCircuitCost,
            0, 255, ifft_nonsentinel_average, false, false,
            "near-average-gate inverse butterfly"),
        SelectCircuitCase(&leopard::ff8xor::GetIFFTCircuitCost,
            0, 255, 0., false, true,
            "maximum-gate inverse butterfly")
    };
    for (size_t case_index = 0;
         case_index < sizeof(butterfly_cases) / sizeof(butterfly_cases[0]);
         ++case_index)
    {
        const CircuitCase circuit = butterfly_cases[case_index];
        const CircuitCase inverse_circuit =
            inverse_butterfly_cases[case_index];
        FillRandom(a, bytes, UINT64_C(0x6d6963726f5f6100) ^
            circuit.coefficient);
        FillRandom(b, bytes, UINT64_C(0x6d6963726f5f6200) ^
            circuit.coefficient);
        memcpy(d[0], a[0], static_cast<size_t>(bytes));
        memcpy(e[0], b[0], static_cast<size_t>(bytes));

        // Verify both circuit orders are mutual inverses outside timing.
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, a[0], b[0], static_cast<uint8_t>(circuit.coefficient));
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, a[0], b[0], static_cast<uint8_t>(circuit.coefficient));
        if (!Equal(a[0], d[0], bytes) || !Equal(b[0], e[0], bytes))
        {
            std::cerr << "microbenchmark FFT/IFFT inverse check failed for skew="
                << circuit.coefficient << '\n';
            return false;
        }

        if (inverse_circuit.coefficient != circuit.coefficient)
        {
            leopard::ff8xor::IFFTButterflyBuffer(
                bytes, a[0], b[0],
                static_cast<uint8_t>(inverse_circuit.coefficient));
            leopard::ff8xor::FFTButterflyBuffer(
                bytes, a[0], b[0],
                static_cast<uint8_t>(inverse_circuit.coefficient));
            if (!Equal(a[0], d[0], bytes) || !Equal(b[0], e[0], bytes))
            {
                std::cerr << "microbenchmark selected IFFT/FFT inverse check "
                    "failed for skew=" << inverse_circuit.coefficient << '\n';
                return false;
            }
        }
        leopard::ff8xor::IFFTButterflyBuffer(
            bytes, a[0], b[0], static_cast<uint8_t>(circuit.coefficient));
        leopard::ff8xor::FFTButterflyBuffer(
            bytes, a[0], b[0], static_cast<uint8_t>(circuit.coefficient));
        if (!Equal(a[0], d[0], bytes) || !Equal(b[0], e[0], bytes))
        {
            std::cerr << "microbenchmark IFFT/FFT inverse check failed for skew="
                << circuit.coefficient << '\n';
            return false;
        }

        Timing fft_timing;
        if (!Measure(options,
            [&]() -> LeopardResult {
                leopard::ff8xor::FFTButterflyBuffer(
                    bytes, a[0], b[0],
                    static_cast<uint8_t>(circuit.coefficient));
                BenchmarkSink ^= a[0][0];
                return Leopard_Success;
            }, fft_timing, result))
            return false;
        PrintMicro(options, reporter, "ff8xor_native", "fft_butterfly",
            CircuitCaseNote("skew", circuit, "generated two-buffer kernel"),
            fft_timing, bytes, bytes * 2, bytes * 2,
            bytes * (circuit.coefficient == 255 ? 3 : 4));

        memcpy(a[0], d[0], static_cast<size_t>(bytes));
        memcpy(b[0], e[0], static_cast<size_t>(bytes));
        Timing ifft_timing;
        if (!Measure(options,
            [&]() -> LeopardResult {
                leopard::ff8xor::IFFTButterflyBuffer(
                    bytes, a[0], b[0],
                    static_cast<uint8_t>(inverse_circuit.coefficient));
                BenchmarkSink ^= b[0][0];
                return Leopard_Success;
            }, ifft_timing, result))
            return false;
        PrintMicro(options, reporter, "ff8xor_native", "ifft_butterfly",
            CircuitCaseNote("skew", inverse_circuit,
                "generated two-buffer kernel"),
            ifft_timing, bytes, bytes * 2, bytes * 2,
            bytes * (inverse_circuit.coefficient == 255 ? 3 : 4));
    }

    const uint64_t transpose_sizes[] = {
        512, 1024, 4096, 64 * 1024, 1024 * 1024
    };
    const size_t transpose_size_count = base_options.quick ? 4 : 5;
    const bool avx2_transpose =
        leopard::ff8xor::transpose::IsModeAvailable(
            leopard::ff8xor::transpose::Mode::Avx2);
    const bool bitalg_transpose =
        leopard::ff8xor::transpose::IsModeAvailable(
            leopard::ff8xor::transpose::Mode::Avx512Bitalg);
    const bool vbmi_transpose =
        leopard::ff8xor::transpose::IsModeAvailable(
            leopard::ff8xor::transpose::Mode::Avx512Vbmi);
    for (size_t size_index = 0;
         size_index < transpose_size_count; ++size_index)
    {
        const uint64_t transpose_bytes = transpose_sizes[size_index];
        FillRandom(a, transpose_bytes,
            UINT64_C(0x7472616e73706f73) ^ transpose_bytes);
        leopard::ff8xor::transpose::PackedToPlane(
            a[0], d[0], transpose_bytes,
            leopard::ff8xor::transpose::Mode::Portable);
        leopard::ff8xor::transpose::PackedToPlane(
            a[0], c[0], transpose_bytes,
            leopard::ff8xor::transpose::Mode::Auto);
        if (!Equal(c[0], d[0], transpose_bytes))
        {
            std::cerr << "microbenchmark transpose backend mismatch at bytes="
                << transpose_bytes << '\n';
            return false;
        }
        leopard::ff8xor::transpose::PlaneToPacked(
            c[0], b[0], transpose_bytes,
            leopard::ff8xor::transpose::Mode::Auto);
        if (!Equal(a[0], b[0], transpose_bytes))
        {
            std::cerr << "microbenchmark transpose round trip failed at bytes="
                << transpose_bytes << '\n';
            return false;
        }

        Timing portable_timing;
        Timing auto_timing;
        if (!MeasurePairABBA(options,
            [&]() -> LeopardResult {
                leopard::ff8xor::transpose::PackedToPlane(
                    a[0], d[0], transpose_bytes,
                    leopard::ff8xor::transpose::Mode::Portable);
                BenchmarkSink ^= d[0][0];
                return Leopard_Success;
            },
            [&]() -> LeopardResult {
                leopard::ff8xor::transpose::PackedToPlane(
                    a[0], c[0], transpose_bytes,
                    leopard::ff8xor::transpose::Mode::Auto);
                BenchmarkSink ^= c[0][0];
                return Leopard_Success;
            }, portable_timing, auto_timing, result))
            return false;
        PrintMicro(options, reporter, "transpose",
            "packed_to_plane_portable",
            "portable 8x8 word transpose; no allocation; paired ABBA control",
            portable_timing, transpose_bytes, transpose_bytes,
            transpose_bytes, transpose_bytes * 2, "ABBA");
        PrintMicro(options, reporter, "transpose", "packed_to_plane",
            bitalg_transpose ?
                "auto-dispatched AVX-512BITALG ZMM 64-byte blocks; no allocation; paired ABBA" :
            avx2_transpose ?
                "auto-dispatched blocked AVX2 8x8 transpose; portable tail; no allocation; paired ABBA" :
                "portable fallback (AVX2 transpose unavailable); no allocation; paired ABBA",
            auto_timing, transpose_bytes, transpose_bytes,
            transpose_bytes, transpose_bytes * 2, "ABBA");

        if (avx2_transpose && bitalg_transpose)
        {
            Timing avx2_timing;
            Timing bitalg_timing;
            if (!MeasurePairABBA(options,
                [&]() -> LeopardResult {
                    leopard::ff8xor::transpose::PackedToPlane(
                        a[0], d[0], transpose_bytes,
                        leopard::ff8xor::transpose::Mode::Avx2);
                    BenchmarkSink ^= d[0][0];
                    return Leopard_Success;
                },
                [&]() -> LeopardResult {
                    leopard::ff8xor::transpose::PackedToPlane(
                        a[0], c[0], transpose_bytes,
                        leopard::ff8xor::transpose::Mode::Avx512Bitalg);
                    BenchmarkSink ^= c[0][0];
                    return Leopard_Success;
                }, avx2_timing, bitalg_timing, result))
                return false;
            PrintMicro(options, reporter, "transpose",
                "packed_to_plane_avx2_control",
                "forced AVX2 control for retained-width comparison; paired ABBA",
                avx2_timing, transpose_bytes, transpose_bytes,
                transpose_bytes, transpose_bytes * 2, "ABBA");
            PrintMicro(options, reporter, "transpose",
                "packed_to_plane_bitalg_zmm",
                "forced AVX-512BITALG ZMM; paired against AVX2; no allocation",
                bitalg_timing, transpose_bytes, transpose_bytes,
                transpose_bytes, transpose_bytes * 2, "ABBA");
        }

        Timing inverse_portable_timing;
        Timing inverse_timing;
        if (!MeasurePairABBA(options,
            [&]() -> LeopardResult {
                leopard::ff8xor::transpose::PlaneToPacked(
                    c[0], e[0], transpose_bytes,
                    leopard::ff8xor::transpose::Mode::Portable);
                BenchmarkSink ^= e[0][0];
                return Leopard_Success;
            },
            [&]() -> LeopardResult {
                leopard::ff8xor::transpose::PlaneToPacked(
                    c[0], b[0], transpose_bytes,
                    leopard::ff8xor::transpose::Mode::Auto);
                BenchmarkSink ^= b[0][0];
                return Leopard_Success;
            }, inverse_portable_timing, inverse_timing, result))
            return false;
        PrintMicro(options, reporter, "transpose",
            "plane_to_packed_portable",
            "portable inverse 8x8 word transpose; no allocation; paired ABBA control",
            inverse_portable_timing, transpose_bytes, transpose_bytes,
            transpose_bytes, transpose_bytes * 2, "ABBA");
        PrintMicro(options, reporter, "transpose", "plane_to_packed",
            vbmi_transpose ?
                "auto-dispatched AVX-512VBMI ZMM 512-byte hierarchy; AVX2/portable tail; no allocation; paired ABBA" :
            avx2_transpose ?
                "auto-dispatched blocked AVX2 inverse 8x8 transpose; portable tail; no allocation; paired ABBA" :
                "portable inverse fallback (AVX2 transpose unavailable); no allocation; paired ABBA",
            inverse_timing, transpose_bytes, transpose_bytes,
            transpose_bytes, transpose_bytes * 2, "ABBA");

        if (avx2_transpose && vbmi_transpose)
        {
            Timing inverse_avx2_timing;
            Timing vbmi_timing;
            if (!MeasurePairABBA(options,
                [&]() -> LeopardResult {
                    leopard::ff8xor::transpose::PlaneToPacked(
                        c[0], e[0], transpose_bytes,
                        leopard::ff8xor::transpose::Mode::Avx2);
                    BenchmarkSink ^= e[0][0];
                    return Leopard_Success;
                },
                [&]() -> LeopardResult {
                    leopard::ff8xor::transpose::PlaneToPacked(
                        c[0], b[0], transpose_bytes,
                        leopard::ff8xor::transpose::Mode::Avx512Vbmi);
                    BenchmarkSink ^= b[0][0];
                    return Leopard_Success;
                }, inverse_avx2_timing, vbmi_timing, result))
                return false;
            PrintMicro(options, reporter, "transpose",
                "plane_to_packed_avx2_control",
                "forced AVX2 control for retained-width comparison; paired ABBA",
                inverse_avx2_timing, transpose_bytes, transpose_bytes,
                transpose_bytes, transpose_bytes * 2, "ABBA");
            PrintMicro(options, reporter, "transpose",
                "plane_to_packed_vbmi_zmm",
                "forced AVX-512VBMI ZMM; paired against AVX2; portable tail",
                vbmi_timing, transpose_bytes, transpose_bytes,
                transpose_bytes, transpose_bytes * 2, "ABBA");
        }
    }

    if (!base_options.csv && !base_options.json)
        std::cout << '\n';
    return true;
}

static bool ParseFF8XorMode(
    const std::string& text,
    leopard::ff8xor::KernelMode& mode)
{
    if (text == "auto")
        mode = leopard::ff8xor::KernelMode::Auto;
    else if (text == "portable")
        mode = leopard::ff8xor::KernelMode::Portable;
    else if (text == "simd128")
        mode = leopard::ff8xor::KernelMode::Simd128;
    else if (text == "avx2")
        mode = leopard::ff8xor::KernelMode::Avx2;
    else if (text == "avx512vl")
        mode = leopard::ff8xor::KernelMode::Avx512VL;
    else if (text == "avx512zmm")
        mode = leopard::ff8xor::KernelMode::Avx512Zmm;
    else
        return false;
    return true;
}

static bool ParseFourBufferMode(
    const std::string& text,
    leopard::ff8xor::FourBufferMode& mode)
{
    if (text == "disabled")
        mode = leopard::ff8xor::FourBufferMode::Disabled;
    else if (text == "xor2")
        mode = leopard::ff8xor::FourBufferMode::Xor2;
    else if (text == "xor3")
        mode = leopard::ff8xor::FourBufferMode::Xor3;
    else
        return false;
    return true;
}

static void Usage(const char* program)
{
    std::cout << "Usage: " << program
        << " [--quick] [--csv|--json] [--include-transpose] [--counters]"
        << " [--portable-transpose] [--abba] [--cache-color]"
        << " [--ff8xor-mode MODE]"
        << " [--four-buffer-mode MODE]"
        << " [--cpu N|--no-pin]\n"
        << "  --quick              Run a bounded development/CI subset.\n"
        << "  --csv                Emit machine-readable CSV.\n"
        << "  --json               Emit newline-delimited machine-readable JSON.\n"
        << "  --include-transpose  Also time packed boundary transposes around ff8xor.\n"
        << "  --portable-transpose Use the portable transpose control for packed-boundary rows.\n"
        << "  --counters           Request Linux perf_event PMU counters; unavailable"
        << " events remain null.\n"
        << "  --abba               Measure paired end-to-end rows in A-B-B-A order.\n"
        << "  --cache-color        Color native transform buffers by page-relative"
        << " L1D set when all eight plane starts alias.\n"
        << "  --ff8xor-mode MODE   Force auto, portable, simd128, avx2,"
        << " avx512vl, or avx512zmm; unavailable modes skip.\n"
        << "  --four-buffer-mode MODE  Select disabled, xor2, or xor3"
        << " AVX-512 radix-4 circuits (default: disabled).\n"
        << "  --cpu N              Pin the benchmark to allowed logical CPU N.\n"
        << "  --no-pin             Disable default pinning to the first allowed CPU.\n";
}

static bool ParseOptions(int argc, char** argv, Options& options)
{
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--quick")
            options.quick = true;
        else if (argument == "--csv")
            options.csv = true;
        else if (argument == "--json")
            options.json = true;
        else if (argument == "--include-transpose")
            options.include_transpose = true;
        else if (argument == "--portable-transpose")
            options.portable_transpose = true;
        else if (argument == "--counters")
            options.counters = true;
        else if (argument == "--abba")
            options.abba = true;
        else if (argument == "--cache-color")
            options.cache_color = true;
        else if (argument == "--ff8xor-mode")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "--ff8xor-mode requires a mode name\n";
                return false;
            }
            options.ff8xor_mode_name = argv[++i];
            if (!ParseFF8XorMode(
                    options.ff8xor_mode_name, options.ff8xor_mode))
            {
                std::cerr << "Invalid --ff8xor-mode value: "
                    << options.ff8xor_mode_name << '\n';
                return false;
            }
        }
        else if (argument == "--four-buffer-mode")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "--four-buffer-mode requires a mode name\n";
                return false;
            }
            options.four_buffer_mode_name = argv[++i];
            if (!ParseFourBufferMode(
                    options.four_buffer_mode_name,
                    options.four_buffer_mode))
            {
                std::cerr << "Invalid --four-buffer-mode value: "
                    << options.four_buffer_mode_name << '\n';
                return false;
            }
        }
        else if (argument == "--no-pin")
            options.pin = false;
        else if (argument == "--cpu")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "--cpu requires a logical CPU number\n";
                return false;
            }
            char* end = NULL;
            errno = 0;
            const long value = strtol(argv[++i], &end, 10);
            if (errno != 0 || end == argv[i] || *end != '\0' ||
                value < 0 || value > std::numeric_limits<int>::max())
            {
                std::cerr << "Invalid --cpu value: " << argv[i] << '\n';
                return false;
            }
            options.pin = true;
            options.pin_cpu = static_cast<int>(value);
        }
        else if (argument == "--help" || argument == "-h")
        {
            Usage(argv[0]);
            std::exit(0);
        }
        else
        {
            std::cerr << "Unknown option: " << argument << '\n';
            Usage(argv[0]);
            return false;
        }
    }
    if (options.csv && options.json)
    {
        std::cerr << "--csv and --json are mutually exclusive\n";
        return false;
    }
    if (options.quick)
    {
        options.warmups = 1;
        options.iterations = 3;
        options.minimum_sample_usec = 250.;
    }
    return true;
}

struct Counts
{
    unsigned original;
    unsigned recovery;
};

} // namespace

int main(int argc, char** argv)
{
    Options options;
    if (!ParseOptions(argc, argv, options))
        return 2;

    if (options.four_buffer_mode !=
            leopard::ff8xor::FourBufferMode::Disabled &&
        options.ff8xor_mode != leopard::ff8xor::KernelMode::Avx512Zmm)
    {
        std::cerr << "--four-buffer-mode requires --ff8xor-mode avx512zmm\n";
        return 2;
    }

    BoundaryTransposeMode = options.portable_transpose ?
        leopard::ff8xor::transpose::Mode::Portable :
        leopard::ff8xor::transpose::Mode::Auto;

    const std::string affinity = ConfigureAffinity(options);

    const LeopardResult init_result = static_cast<LeopardResult>(leo_init());
    if (init_result != Leopard_Success)
    {
        std::cerr << "leo_init failed: " << leo_result_string(init_result) << '\n';
        return 1;
    }

    if (options.ff8xor_mode != leopard::ff8xor::KernelMode::Auto &&
        !leopard::ff8xor::IsKernelModeAvailable(options.ff8xor_mode))
    {
        std::cerr << "Requested ff8xor mode is unavailable on this build/host: "
            << options.ff8xor_mode_name << '\n';
        return 77;
    }
    leopard::ff8xor::SetKernelMode(options.ff8xor_mode);
    leopard::ff8xor::SetFourBufferMode(options.four_buffer_mode);

    Environment environment = GetEnvironment(affinity);
    environment.ff8xor_mode_requested = options.ff8xor_mode_name;
    environment.four_buffer_mode_requested = options.four_buffer_mode_name;
    Reporter reporter(options.csv, options.json, environment);
    reporter.Begin(options);

    if (!CheckTransposeHelper())
    {
        std::cerr << "8x8 transpose self-test failed\n";
        return 1;
    }
    if (!CheckCacheColorHelper())
    {
        std::cerr << "cache-color allocator self-test failed\n";
        return 1;
    }

    if (!RunMicrobenchmarks(options, reporter))
        return 1;

    static const Counts full_counts[] = {
        { 8, 2 }, { 16, 4 }, { 32, 8 },
        { 64, 16 }, { 128, 32 }, { 128, 128 }
    };
    static const Counts quick_counts[] = {
        { 8, 2 }, { 32, 8 }
    };
    static const uint64_t full_bytes[] = {
        1024, 4096, 64 * 1024, 1024 * 1024
    };
    static const uint64_t quick_bytes[] = {
        1024, 64 * 1024
    };

    const Counts* counts = options.quick ? quick_counts : full_counts;
    const size_t counts_size = options.quick
        ? sizeof(quick_counts) / sizeof(quick_counts[0])
        : sizeof(full_counts) / sizeof(full_counts[0]);
    const uint64_t* sizes = options.quick ? quick_bytes : full_bytes;
    const size_t sizes_count = options.quick
        ? sizeof(quick_bytes) / sizeof(quick_bytes[0])
        : sizeof(full_bytes) / sizeof(full_bytes[0]);

    try
    {
        for (size_t pair = 0; pair < counts_size; ++pair)
        {
            for (size_t size = 0; size < sizes_count; ++size)
            {
                if (!RunParameter(options, reporter,
                        counts[pair].original, counts[pair].recovery,
                        sizes[size]))
                    return 1;
            }
        }
    }
    catch (const std::bad_alloc&)
    {
        std::cerr << "allocation failed while constructing benchmark metadata\n";
        return 1;
    }

    if (!options.csv && !options.json)
        std::cout << "Benchmark complete. sink=" << BenchmarkSink << '\n';
    return 0;
}
