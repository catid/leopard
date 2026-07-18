/*
    Independent metadata tests for the generated FF8 XOR circuits.

    The generated circuits are executed with symbolic wire values.  Every
    literal XorValue(destination, source) records one CNOT and advances the
    conservative dependency layer of both wires.  This derives counts and
    depths from the emitted circuit bodies rather than trusting their emitted
    metadata arrays.
*/

#include "../LeopardFF8Xor.h"

#include <algorithm>
#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

namespace leopard { namespace ff8xor {

struct StatisticsTraceValue
{
    unsigned Wire;
};

struct StatisticsTraceContext
{
    explicit StatisticsTraceContext(unsigned width)
        : Width(width)
        , GateCount(0)
        , MaximumDepth(0)
        , Valid(width <= 16)
    {
        memset(LastLayer, 0, sizeof(LastLayer));
    }

    unsigned Width;
    unsigned GateCount;
    unsigned MaximumDepth;
    unsigned LastLayer[16];
    bool Valid;
};

static StatisticsTraceContext* ActiveStatisticsTrace = NULL;

static LEO_FORCE_INLINE StatisticsTraceValue XorValue(
    StatisticsTraceValue destination,
    StatisticsTraceValue source)
{
    StatisticsTraceContext* const trace = ActiveStatisticsTrace;
    if (trace == NULL ||
        destination.Wire >= trace->Width ||
        source.Wire >= trace->Width ||
        destination.Wire == source.Wire)
    {
        if (trace != NULL)
            trace->Valid = false;
        return destination;
    }

    const unsigned layer = 1 + std::max(
        trace->LastLayer[destination.Wire],
        trace->LastLayer[source.Wire]);
    trace->LastLayer[destination.Wire] = layer;
    trace->LastLayer[source.Wire] = layer;
    trace->MaximumDepth = std::max(trace->MaximumDepth, layer);
    ++trace->GateCount;
    return destination;
}

}} // namespace leopard::ff8xor

#include "../generated/LeopardFF8XorCircuits.inl"

namespace {

struct CircuitObservation
{
    unsigned GateCount;
    unsigned Depth;
    bool Valid;
};

template <unsigned Coefficient>
static CircuitObservation ObserveMultiplyCircuit()
{
    leopard::ff8xor::StatisticsTraceContext trace(8);
    leopard::ff8xor::StatisticsTraceValue wires[8];
    for (unsigned wire = 0; wire < 8; ++wire)
        wires[wire].Wire = wire;

    leopard::ff8xor::ActiveStatisticsTrace = &trace;
    leopard::ff8xor::generated::MultiplyCircuit<Coefficient>::Apply(
        wires[0], wires[1], wires[2], wires[3],
        wires[4], wires[5], wires[6], wires[7]);
    leopard::ff8xor::ActiveStatisticsTrace = NULL;

    CircuitObservation result = {
        trace.GateCount, trace.MaximumDepth, trace.Valid
    };
    return result;
}

template <unsigned Skew, bool Inverse>
static CircuitObservation ObserveButterflyCircuit()
{
    leopard::ff8xor::StatisticsTraceContext trace(16);
    leopard::ff8xor::StatisticsTraceValue wires[16];
    for (unsigned wire = 0; wire < 16; ++wire)
        wires[wire].Wire = wire;

    leopard::ff8xor::ActiveStatisticsTrace = &trace;
    if (Inverse)
    {
        leopard::ff8xor::generated::IFFTCircuit<Skew>::Apply(
            wires[0], wires[1], wires[2], wires[3],
            wires[4], wires[5], wires[6], wires[7],
            wires[8], wires[9], wires[10], wires[11],
            wires[12], wires[13], wires[14], wires[15]);
    }
    else
    {
        leopard::ff8xor::generated::FFTCircuit<Skew>::Apply(
            wires[0], wires[1], wires[2], wires[3],
            wires[4], wires[5], wires[6], wires[7],
            wires[8], wires[9], wires[10], wires[11],
            wires[12], wires[13], wires[14], wires[15]);
    }
    leopard::ff8xor::ActiveStatisticsTrace = NULL;

    CircuitObservation result = {
        trace.GateCount, trace.MaximumDepth, trace.Valid
    };
    return result;
}

typedef CircuitObservation (*CircuitObserver)();

#define LEO_FF8XOR_MULTIPLY_OBSERVER(Coefficient) \
    &ObserveMultiplyCircuit<Coefficient>,
static const CircuitObserver MultiplyObservers[256] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_MULTIPLY_OBSERVER)
};
#undef LEO_FF8XOR_MULTIPLY_OBSERVER

#define LEO_FF8XOR_FFT_OBSERVER(Skew) \
    &ObserveButterflyCircuit<Skew, false>,
static const CircuitObserver FFTObservers[256] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_FFT_OBSERVER)
};
#undef LEO_FF8XOR_FFT_OBSERVER

#define LEO_FF8XOR_IFFT_OBSERVER(Skew) \
    &ObserveButterflyCircuit<Skew, true>,
static const CircuitObserver IFFTObservers[256] = {
    LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_FF8XOR_IFFT_OBSERVER)
};
#undef LEO_FF8XOR_IFFT_OBSERVER

struct FamilyFixture
{
    unsigned MinimumGateCount;
    unsigned MaximumGateCount;
    uint64_t GateCountSum;
    unsigned MinimumDepth;
    unsigned MaximumDepth;
    uint64_t DepthSum;
};

struct GeneratedStatistics
{
    unsigned MinimumGateCount;
    unsigned MaximumGateCount;
    double AverageGateCount;
    unsigned MinimumDepth;
    unsigned MaximumDepth;
    double AverageDepth;
};

typedef leopard::ff8xor::CircuitCost (*CircuitCostGetter)(
    leopard::ff8xor::ffe_t coefficient);

static bool CheckUnsigned(
    const char* family,
    const char* field,
    unsigned actual,
    unsigned expected)
{
    if (actual == expected)
        return true;
    fprintf(stderr, "%s %s mismatch: expected=%u actual=%u\n",
        family, field, expected, actual);
    return false;
}

static bool CheckSum(
    const char* family,
    const char* field,
    uint64_t actual,
    uint64_t expected)
{
    if (actual == expected)
        return true;
    fprintf(stderr, "%s %s mismatch: expected=%llu actual=%llu\n",
        family,
        field,
        static_cast<unsigned long long>(expected),
        static_cast<unsigned long long>(actual));
    return false;
}

static bool CheckAverage(
    const char* family,
    const char* field,
    double actual,
    double expected)
{
    if (actual == expected)
        return true;
    fprintf(stderr, "%s %s mismatch: expected=%.9f actual=%.9f\n",
        family, field, expected, actual);
    return false;
}

static bool CheckFamily(
    const char* family,
    const CircuitObserver* observers,
    const uint8_t* generated_gate_counts,
    const uint8_t* generated_depths,
    CircuitCostGetter published_cost,
    const GeneratedStatistics& generated,
    const leopard::ff8xor::CircuitStatistics& published,
    const FamilyFixture& fixture)
{
    unsigned minimum_gate_count = std::numeric_limits<unsigned>::max();
    unsigned maximum_gate_count = 0;
    uint64_t gate_count_sum = 0;
    unsigned minimum_depth = std::numeric_limits<unsigned>::max();
    unsigned maximum_depth = 0;
    uint64_t depth_sum = 0;
    bool ok = true;

    for (unsigned index = 0; index < 256; ++index)
    {
        const CircuitObservation observed = observers[index]();
        if (!observed.Valid)
        {
            fprintf(stderr, "%s circuit %u produced an invalid trace\n",
                family, index);
            ok = false;
            continue;
        }
        if (observed.GateCount != generated_gate_counts[index])
        {
            fprintf(stderr,
                "%s circuit %u gate count mismatch: generated=%u observed=%u\n",
                family,
                index,
                static_cast<unsigned>(generated_gate_counts[index]),
                observed.GateCount);
            ok = false;
        }
        if (observed.Depth != generated_depths[index])
        {
            fprintf(stderr,
                "%s circuit %u depth mismatch: generated=%u observed=%u\n",
                family,
                index,
                static_cast<unsigned>(generated_depths[index]),
                observed.Depth);
            ok = false;
        }

        const leopard::ff8xor::CircuitCost cost = published_cost(
            static_cast<leopard::ff8xor::ffe_t>(index));
        if (cost.GateCount != observed.GateCount ||
            cost.Depth != observed.Depth)
        {
            fprintf(stderr,
                "%s circuit %u published cost mismatch: "
                "published=%u/%u observed=%u/%u\n",
                family,
                index,
                cost.GateCount,
                cost.Depth,
                observed.GateCount,
                observed.Depth);
            ok = false;
        }

        minimum_gate_count = std::min(minimum_gate_count, observed.GateCount);
        maximum_gate_count = std::max(maximum_gate_count, observed.GateCount);
        gate_count_sum += observed.GateCount;
        minimum_depth = std::min(minimum_depth, observed.Depth);
        maximum_depth = std::max(maximum_depth, observed.Depth);
        depth_sum += observed.Depth;
    }

    const double average_gate_count = gate_count_sum / 256.0;
    const double average_depth = depth_sum / 256.0;

    ok = CheckUnsigned(family, "fixture minimum gate count",
        minimum_gate_count, fixture.MinimumGateCount) && ok;
    ok = CheckUnsigned(family, "fixture maximum gate count",
        maximum_gate_count, fixture.MaximumGateCount) && ok;
    ok = CheckSum(family, "fixture gate count sum",
        gate_count_sum, fixture.GateCountSum) && ok;
    ok = CheckUnsigned(family, "fixture minimum depth",
        minimum_depth, fixture.MinimumDepth) && ok;
    ok = CheckUnsigned(family, "fixture maximum depth",
        maximum_depth, fixture.MaximumDepth) && ok;
    ok = CheckSum(family, "fixture depth sum",
        depth_sum, fixture.DepthSum) && ok;

    ok = CheckUnsigned(family, "generated minimum gate count",
        generated.MinimumGateCount, minimum_gate_count) && ok;
    ok = CheckUnsigned(family, "generated maximum gate count",
        generated.MaximumGateCount, maximum_gate_count) && ok;
    ok = CheckAverage(family, "generated average gate count",
        generated.AverageGateCount, average_gate_count) && ok;
    ok = CheckUnsigned(family, "generated minimum depth",
        generated.MinimumDepth, minimum_depth) && ok;
    ok = CheckUnsigned(family, "generated maximum depth",
        generated.MaximumDepth, maximum_depth) && ok;
    ok = CheckAverage(family, "generated average depth",
        generated.AverageDepth, average_depth) && ok;

    ok = CheckUnsigned(family, "published minimum gate count",
        published.MinimumGateCount, minimum_gate_count) && ok;
    ok = CheckUnsigned(family, "published maximum gate count",
        published.MaximumGateCount, maximum_gate_count) && ok;
    ok = CheckAverage(family, "published average gate count",
        published.AverageGateCount, average_gate_count) && ok;
    ok = CheckUnsigned(family, "published minimum depth",
        published.MinimumDepth, minimum_depth) && ok;
    ok = CheckUnsigned(family, "published maximum depth",
        published.MaximumDepth, maximum_depth) && ok;
    ok = CheckAverage(family, "published average depth",
        published.AverageDepth, average_depth) && ok;

    if (ok)
    {
        printf("%s: gates=%u..%u avg=%.9f depth=%u..%u avg=%.9f\n",
            family,
            minimum_gate_count,
            maximum_gate_count,
            average_gate_count,
            minimum_depth,
            maximum_depth,
            average_depth);
    }
    return ok;
}

} // namespace

int main()
{
    using namespace leopard::ff8xor::generated;

    const GeneratedStatistics multiply_generated = {
        kMultiplyMinGateCount,
        kMultiplyMaxGateCount,
        kMultiplyAverageGateCount,
        kMultiplyMinDepth,
        kMultiplyMaxDepth,
        kMultiplyAverageDepth
    };
    const GeneratedStatistics fft_generated = {
        kFFTMinGateCount,
        kFFTMaxGateCount,
        kFFTAverageGateCount,
        kFFTMinDepth,
        kFFTMaxDepth,
        kFFTAverageDepth
    };
    const GeneratedStatistics ifft_generated = {
        kIFFTMinGateCount,
        kIFFTMaxGateCount,
        kIFFTAverageGateCount,
        kIFFTMinDepth,
        kIFFTMaxDepth,
        kIFFTAverageDepth
    };

    // These fixtures were computed independently from the generated summaries.
    // Integer sums make the averages exact binary fractions when divided by 256.
    const FamilyFixture multiply_fixture = { 0, 23, 4903, 0, 17, 2912 };
    const FamilyFixture fft_fixture = { 8, 51, 10240, 1, 14, 2780 };
    const FamilyFixture ifft_fixture = { 8, 51, 10240, 1, 14, 2780 };

    bool ok = true;
    const leopard::ff8xor::CircuitCost identity_zero =
        leopard::ff8xor::GetMultiplyCircuitCost(0);
    const leopard::ff8xor::CircuitCost identity_redundant =
        leopard::ff8xor::GetMultiplyCircuitCost(255);
    if (identity_zero.GateCount != 0 || identity_zero.Depth != 0 ||
        identity_redundant.GateCount != 0 ||
        identity_redundant.Depth != 0)
    {
        fprintf(stderr,
            "multiply identity costs must both be zero: log0=%u/%u "
            "log255=%u/%u\n",
            identity_zero.GateCount,
            identity_zero.Depth,
            identity_redundant.GateCount,
            identity_redundant.Depth);
        ok = false;
    }
    ok = CheckFamily(
        "multiply",
        MultiplyObservers,
        kMultiplyGateCounts,
        kMultiplyDepths,
        &leopard::ff8xor::GetMultiplyCircuitCost,
        multiply_generated,
        leopard::ff8xor::GetMultiplyCircuitStatistics(),
        multiply_fixture) && ok;
    ok = CheckFamily(
        "FFT",
        FFTObservers,
        kFFTGateCounts,
        kFFTDepths,
        &leopard::ff8xor::GetFFTCircuitCost,
        fft_generated,
        leopard::ff8xor::GetFFTCircuitStatistics(),
        fft_fixture) && ok;
    ok = CheckFamily(
        "IFFT",
        IFFTObservers,
        kIFFTGateCounts,
        kIFFTDepths,
        &leopard::ff8xor::GetIFFTCircuitCost,
        ifft_generated,
        leopard::ff8xor::GetIFFTCircuitStatistics(),
        ifft_fixture) && ok;

    if (!ok)
        return 1;
    puts("FF8 XOR circuit statistics passed");
    return 0;
}
