/*
    Exact differential tests for the experimental FF8 XOR error-locator
    shift selector.  The oracle below deliberately uses the original
    shift-major objective and the public multiplier-cost API, not the new
    generated rotation rows.
*/

#include "../LeopardFF8Xor.h"
#include "../generated/LeopardFF8XorLocatorRotations.inl"

#include <stdint.h>
#include <stdio.h>

namespace {

static const unsigned kOrder = 256;
static const unsigned kModulus = 255;

struct CostTable
{
    unsigned Gates[kOrder];
    unsigned Depths[kOrder];
};

struct ReferenceResult
{
    unsigned Shift;
    unsigned Gates;
    unsigned Depth;
    unsigned GateTieCount;
    unsigned ExactTieCount;
};

static CostTable Costs;
static uint64_t CheckedCases = 0;

static unsigned CanonicalLog(unsigned logarithm)
{
    return logarithm == kModulus ? 0 : logarithm;
}

static unsigned ShiftedLog(unsigned logarithm, unsigned shift)
{
    unsigned result = CanonicalLog(logarithm) + shift;
    if (result >= kModulus)
        result -= kModulus;
    return result;
}

static unsigned ScoredLog(
    unsigned logarithm,
    bool inverse,
    unsigned shift)
{
    const unsigned shifted = ShiftedLog(logarithm, shift);
    return inverse && shifted != 0 ? kModulus - shifted : shifted;
}

static void ScoreShift(
    const leopard::ff8xor::ffe_t* logarithms,
    const bool* inverse,
    unsigned count,
    unsigned shift,
    unsigned& gates,
    unsigned& depth)
{
    gates = 0;
    depth = 0;
    for (unsigned index = 0; index < count; ++index)
    {
        const unsigned log = ScoredLog(
            logarithms[index], inverse[index], shift);
        gates += Costs.Gates[log];
        depth += Costs.Depths[log];
    }
}

static ReferenceResult SelectReference(
    const leopard::ff8xor::ffe_t* logarithms,
    const bool* inverse,
    unsigned count)
{
    ReferenceResult result = { 0, UINT32_MAX, UINT32_MAX, 0, 0 };
    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        unsigned gates;
        unsigned depth;
        ScoreShift(logarithms, inverse, count, shift, gates, depth);
        if (gates < result.Gates ||
            (gates == result.Gates && depth < result.Depth))
        {
            result.Shift = shift;
            result.Gates = gates;
            result.Depth = depth;
        }
    }

    for (unsigned shift = 0; shift < kModulus; ++shift)
    {
        unsigned gates;
        unsigned depth;
        ScoreShift(logarithms, inverse, count, shift, gates, depth);
        if (gates == result.Gates)
            ++result.GateTieCount;
        if (gates == result.Gates && depth == result.Depth)
            ++result.ExactTieCount;
    }
    return result;
}

static bool CheckCase(
    const char* label,
    const leopard::ff8xor::ffe_t* logarithms,
    const bool* inverse,
    unsigned count,
    ReferenceResult* observation = NULL)
{
    const ReferenceResult expected = SelectReference(
        logarithms, inverse, count);
    const unsigned selected = leopard::ff8xor::SelectLocatorShiftForTesting(
        logarithms, inverse, count);
    const unsigned retained =
        leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
            logarithms, inverse, count);
    ++CheckedCases;

    if (selected != expected.Shift || retained != expected.Shift)
    {
        fprintf(stderr,
            "%s selector mismatch: count=%u expected=%u fast=%u old=%u "
            "gates=%u depth=%u\n",
            label,
            count,
            expected.Shift,
            selected,
            retained,
            expected.Gates,
            expected.Depth);
        return false;
    }
    if (observation != NULL)
        *observation = expected;
    return true;
}

static bool CheckGeneratedRotations()
{
    using namespace leopard::ff8xor::generated;
    for (unsigned raw_log = 0; raw_log < kOrder; ++raw_log)
    {
        const unsigned base = CanonicalLog(raw_log);
        for (unsigned shift = 0; shift < kOrder; ++shift)
        {
            const unsigned positive_log = ScoredLog(raw_log, false, shift);
            const unsigned inverse_log = ScoredLog(raw_log, true, shift);
            const unsigned offset = base + shift;
            if (kLocatorPositiveGateCosts[offset] !=
                    Costs.Gates[positive_log] ||
                kLocatorPositiveDepthCosts[offset] !=
                    Costs.Depths[positive_log] ||
                kLocatorInverseGateCosts[offset] !=
                    Costs.Gates[inverse_log] ||
                kLocatorInverseDepthCosts[offset] !=
                    Costs.Depths[inverse_log])
            {
                fprintf(stderr,
                    "locator rotation mismatch: raw_log=%u base=%u "
                    "shift=%u offset=%u\n",
                    raw_log, base, shift, offset);
                return false;
            }
        }
    }
    return true;
}

static bool CheckInvalidInputs()
{
    leopard::ff8xor::ffe_t logarithms[kOrder] = {};
    bool inverse[kOrder] = {};
    if (leopard::ff8xor::SelectLocatorShiftForTesting(NULL, NULL, 0) != 0 ||
        leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
            NULL, NULL, 0) != 0 ||
        leopard::ff8xor::SelectLocatorShiftForTesting(
            NULL, inverse, 1) != kModulus ||
        leopard::ff8xor::SelectLocatorShiftForTesting(
            logarithms, NULL, 1) != kModulus ||
        leopard::ff8xor::SelectLocatorShiftForTesting(
            logarithms, inverse, kOrder + 1) != kModulus ||
        leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
            NULL, inverse, 1) != kModulus ||
        leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
            logarithms, NULL, 1) != kModulus ||
        leopard::ff8xor::SelectLocatorShiftReferenceForTesting(
            logarithms, inverse, kOrder + 1) != kModulus)
    {
        fprintf(stderr, "locator selector invalid-input contract failed\n");
        return false;
    }
    return true;
}

static bool CheckSingletonsAndPairs()
{
    leopard::ff8xor::ffe_t logarithms[2];
    bool inverse[2];

    // Raw log 255 is included explicitly in both signs; it must be identical
    // to raw log zero for multiplication-only scoring.
    for (unsigned raw_log = 0; raw_log < kOrder; ++raw_log)
    {
        logarithms[0] = static_cast<leopard::ff8xor::ffe_t>(raw_log);
        for (unsigned sign = 0; sign < 2; ++sign)
        {
            inverse[0] = sign != 0;
            if (!CheckCase("signed singleton", logarithms, inverse, 1))
                return false;
        }
    }

    // Exhaust every ordered pair of canonical base logs and signs.
    for (unsigned left = 0; left < kModulus * 2; ++left)
    {
        logarithms[0] = static_cast<leopard::ff8xor::ffe_t>(left / 2);
        inverse[0] = (left & 1) != 0;
        for (unsigned right = 0; right < kModulus * 2; ++right)
        {
            logarithms[1] = static_cast<leopard::ff8xor::ffe_t>(right / 2);
            inverse[1] = (right & 1) != 0;
            if (!CheckCase("signed pair", logarithms, inverse, 2))
                return false;
        }
    }
    return true;
}

static bool CheckTieBreaks()
{
    leopard::ff8xor::ffe_t depth_logs[2] = { 0, 7 };
    bool depth_inverse[2] = { false, false };
    ReferenceResult depth_result = {};
    if (!CheckCase(
            "depth tie-break", depth_logs, depth_inverse, 2,
            &depth_result) ||
        depth_result.Shift != 248 ||
        depth_result.GateTieCount != 2 ||
        depth_result.ExactTieCount != 1)
    {
        fprintf(stderr,
            "expected {+0,+7} to select shift 248 via depth: "
            "shift=%u gate_ties=%u exact_ties=%u\n",
            depth_result.Shift,
            depth_result.GateTieCount,
            depth_result.ExactTieCount);
        return false;
    }

    leopard::ff8xor::ffe_t cycle_logs[kModulus];
    bool cycle_inverse[kModulus] = {};
    for (unsigned index = 0; index < kModulus; ++index)
        cycle_logs[index] = static_cast<leopard::ff8xor::ffe_t>(index);
    ReferenceResult cycle_result = {};
    if (!CheckCase(
            "numeric tie-break", cycle_logs, cycle_inverse, kModulus,
            &cycle_result) ||
        cycle_result.Shift != 0 ||
        cycle_result.ExactTieCount != kModulus)
    {
        fprintf(stderr,
            "full-cycle lowest-shift tie failed: shift=%u exact_ties=%u\n",
            cycle_result.Shift,
            cycle_result.ExactTieCount);
        return false;
    }
    return true;
}

static uint32_t NextRandom(uint32_t& state)
{
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

static bool CheckLongLists()
{
    leopard::ff8xor::ffe_t logarithms[kOrder];
    bool inverse[kOrder];

    // Shape the maximum term count exactly like k=r=128 with all recovery
    // shards available and a mixture of present/missing originals.
    for (unsigned index = 0; index < kOrder; ++index)
    {
        logarithms[index] = static_cast<leopard::ff8xor::ffe_t>(
            index % 31 == 0 ? kModulus : (index * 73 + 19) % kModulus);
        inverse[index] = index >= 128 && (index % 17 == 0);
    }
    if (!CheckCase("k128 r128 max terms", logarithms, inverse, kOrder))
        return false;

    // Also model the asymmetric k=254,r=2 edge: two available recoveries,
    // then 254 original decisions, while retaining the full 256-term bound.
    for (unsigned index = 0; index < kOrder; ++index)
    {
        logarithms[index] = static_cast<leopard::ff8xor::ffe_t>(
            index % 29 == 0 ? kModulus : (index * 41 + 3) % kModulus);
        inverse[index] = index >= 2 && (index == 2 || index == 173);
    }
    if (!CheckCase("k254 r2 max terms", logarithms, inverse, kOrder))
        return false;

    uint32_t state = UINT32_C(0x51ec70a5);
    for (unsigned trial = 0; trial < 512; ++trial)
    {
        const unsigned count = NextRandom(state) % (kOrder + 1);
        for (unsigned index = 0; index < count; ++index)
        {
            logarithms[index] = static_cast<leopard::ff8xor::ffe_t>(
                NextRandom(state) & 255U);
            inverse[index] = (NextRandom(state) & 1U) != 0;
        }
        if (!CheckCase("deterministic random list",
                logarithms, inverse, count))
            return false;
    }
    return true;
}

} // namespace

int main()
{
    for (unsigned log = 0; log < kOrder; ++log)
    {
        const leopard::ff8xor::CircuitCost cost =
            leopard::ff8xor::GetMultiplyCircuitCost(
                static_cast<leopard::ff8xor::ffe_t>(log));
        Costs.Gates[log] = cost.GateCount;
        Costs.Depths[log] = cost.Depth;
    }

    if (!CheckGeneratedRotations() ||
        !CheckInvalidInputs() ||
        !CheckCase("empty", NULL, NULL, 0) ||
        !CheckSingletonsAndPairs() ||
        !CheckTieBreaks() ||
        !CheckLongLists())
    {
        return 1;
    }

    printf("FF8 XOR locator selector passed: %llu differential cases\n",
        static_cast<unsigned long long>(CheckedCases));
    return 0;
}
