// Compiler-only schedule and allocation regression for leopard-79h.56.
#include "Leopard2Plan.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <vector>

namespace allocation_audit {

// This executable is deliberately single-threaded. Generations keep storage
// owned by an earlier plan outside the current compiler's allocation account.
struct alignas(std::max_align_t) Header { size_t bytes, generation; };
bool active = false;
size_t generation = 0, calls = 0, live = 0, peak = 0, fail_at = 0;
bool failed = false;

struct Measurement
{
    size_t calls, live, peak;
    bool failed;
};

void begin(size_t failure = 0)
{
    ++generation;
    calls = live = peak = 0;
    failed = false;
    fail_at = failure;
    active = true;
}

Measurement end()
{
    active = false;
    const Measurement result = { calls, live, peak, failed };
    return result;
}

} // namespace allocation_audit

#if defined(_MSC_VER)
#define COMPILER_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define COMPILER_TEST_NOINLINE __attribute__((noinline))
#else
#define COMPILER_TEST_NOINLINE
#endif

COMPILER_TEST_NOINLINE void* operator new(size_t bytes)
{
    using namespace allocation_audit;
    if (active && ++calls == fail_at)
    {
        failed = true;
        throw std::bad_alloc();
    }
    if (bytes == 0) bytes = 1;
    if (bytes > SIZE_MAX - sizeof(Header)) throw std::bad_alloc();
    Header* header = static_cast<Header*>(std::malloc(sizeof(Header) + bytes));
    if (!header) throw std::bad_alloc();
    header->bytes = bytes;
    header->generation = active ? generation : 0;
    if (active)
    {
        live += bytes;
        if (live > peak) peak = live;
    }
    return header + 1;
}

COMPILER_TEST_NOINLINE void operator delete(void* pointer) noexcept
{
    using namespace allocation_audit;
    if (!pointer) return;
    Header* header = static_cast<Header*>(pointer) - 1;
    if (active && header->generation == generation) live -= header->bytes;
    std::free(header);
}

COMPILER_TEST_NOINLINE void* operator new[](size_t bytes)
{
    return ::operator new(bytes);
}
COMPILER_TEST_NOINLINE void operator delete[](void* pointer) noexcept
{
    ::operator delete(pointer);
}
COMPILER_TEST_NOINLINE void* operator new(size_t bytes, const std::nothrow_t&) noexcept
{
    try { return ::operator new(bytes); } catch (...) { return NULL; }
}
COMPILER_TEST_NOINLINE void* operator new[](size_t bytes, const std::nothrow_t&) noexcept
{
    return ::operator new(bytes, std::nothrow);
}
COMPILER_TEST_NOINLINE void operator delete(void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}
COMPILER_TEST_NOINLINE void operator delete[](void* pointer, const std::nothrow_t&) noexcept
{
    ::operator delete(pointer);
}
#if defined(__cpp_sized_deallocation)
COMPILER_TEST_NOINLINE void operator delete(void* pointer, size_t) noexcept
{
    ::operator delete(pointer);
}
COMPILER_TEST_NOINLINE void operator delete[](void* pointer, size_t) noexcept
{
    ::operator delete(pointer);
}
#endif
#undef COMPILER_TEST_NOINLINE

namespace {

using leopard2_internal::PrunedTransformPlan;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

// Hash semantic values, never structure padding, host endianness or addresses.
struct Digest
{
    uint64_t value;
    Digest() : value(UINT64_C(14695981039346656037)) {}
    void add(uint64_t input)
    {
        for (unsigned byte = 0; byte < 8; ++byte)
        {
            value ^= (input >> (byte * 8)) & 255U;
            value *= UINT64_C(1099511628211);
        }
    }
    template<class T> void add_vector(const std::vector<T>& input)
    {
        add(input.size());
        for (size_t i = 0; i < input.size(); ++i) add(input[i]);
    }
    void add_plan(const PrunedTransformPlan& plan)
    {
        add(plan.size); add(plan.shift); add(plan.zero_multiplier_log);
        add(plan.first_layer_multiplier_log); add(plan.inverse);
        add_vector(plan.input_mask); add_vector(plan.output_mask);
        add(plan.operations.size());
        for (size_t i = 0; i < plan.operations.size(); ++i)
        {
            const leopard2_internal::PrunedTransformOperation& op = plan.operations[i];
            add(op.x); add(op.y); add(op.multiplier_log); add(op.flags);
        }
        add_vector(plan.inverse_accumulation_flags);
        add_vector(plan.fused_four_starts); add_vector(plan.zero_outputs);
        add(plan.full_butterfly_count); add(plan.one_output_butterflies);
        add(plan.input_zero_specializations);
        add(plan.zero_multiplier_butterflies); add(plan.one_multiplier_butterflies);
    }
};

struct Provider
{
    uint32_t order, mode;
    size_t fail_at;
    bool invalid;
    mutable size_t calls;
    mutable Digest digest;
    Provider(uint32_t field_order, uint32_t kind)
        : order(field_order), mode(kind), fail_at(0), invalid(false), calls(0) {}
};

uint16_t multiplier(const void* opaque, uint32_t index)
{
    const Provider& provider = *static_cast<const Provider*>(opaque);
    provider.digest.add(index);
    if (++provider.calls == provider.fail_at)
    {
        if (provider.invalid) return 256; // Only used with GF8.
        throw std::runtime_error("injected multiplier failure");
    }
    const uint32_t mixed = index * UINT32_C(2654435761) ^ (index >> 3);
    if (provider.mode == 0 || (provider.mode == 3 && (mixed & 7) == 0))
        return static_cast<uint16_t>(provider.order - 1);
    if (provider.mode == 1 || (provider.mode == 3 && (mixed & 7) == 1))
        return 0;
    return static_cast<uint16_t>(1 + mixed % (provider.order - 2));
}

bool compile(uint32_t size, uint32_t shift, bool inverse,
    const std::vector<uint8_t>& input, const std::vector<uint8_t>& output,
    const Provider& provider, PrunedTransformPlan& plan)
{
    return leopard2_internal::CompilePrunedTransformPlan(
        provider.order, static_cast<uint16_t>(provider.order - 1),
        size, shift, inverse, input.data(), output.data(), multiplier, &provider, plan);
}

std::vector<uint8_t> mask(uint32_t size, unsigned pattern)
{
    std::vector<uint8_t> result(size);
    for (uint32_t i = 0; i < size; ++i)
    {
        switch (pattern)
        {
        case 0: result[i] = 0; break;
        case 1: result[i] = 1; break;
        case 2: result[i] = i <= size / 3; break;
        case 3: result[i] = i >= size / 2; break;
        case 4: result[i] = ((i * 17U + (i >> 2)) % 5U) == 0; break;
        default: result[i] = i == size / 2; break;
        }
    }
    return result;
}

void record_case(Digest& digest, uint32_t order, uint32_t size,
    uint32_t shift, bool inverse, unsigned mode,
    const std::vector<uint8_t>& input, const std::vector<uint8_t>& output)
{
    Provider provider(order, mode);
    PrunedTransformPlan plan;
    require(compile(size, shift, inverse, input, output, provider, plan),
        "schedule baseline compilation failed");
    digest.add_plan(plan);
    digest.add(provider.calls); digest.add(provider.digest.value);
}

uint64_t schedule_grid()
{
    Digest digest;
    const unsigned pairs[][2] = {
        {0,0}, {0,1}, {1,0}, {1,1}, {2,1}, {1,2},
        {2,3}, {3,2}, {4,5}, {5,4}, {3,4}, {4,3}
    };
    const uint32_t orders[] = {256, 65536};
    for (unsigned field = 0; field < 2; ++field)
    {
        const uint32_t order = orders[field];
        for (uint32_t size = 2; size <= order; size *= 2)
        for (unsigned coset = 0; coset < (size == order ? 1U : 2U); ++coset)
        for (unsigned inverse = 0; inverse < 2; ++inverse)
        for (unsigned p = 0; p < sizeof(pairs) / sizeof(pairs[0]); ++p)
        {
            const std::vector<uint8_t> input = mask(size, pairs[p][0]);
            const std::vector<uint8_t> output = mask(size, pairs[p][1]);
            record_case(digest, order, size, coset ? order - size : 0,
                inverse != 0, 3, input, output);
            if (size == order || size <= 4)
                for (unsigned mode = 0; mode < 3; ++mode)
                    record_case(digest, order, size, coset ? order - size : 0,
                        inverse != 0, mode, input, output);
        }
    }
    return digest.value;
}

uint64_t exhaustive_small()
{
    Digest digest;
    const uint32_t orders[] = {256, 65536};
    for (unsigned field = 0; field < 2; ++field)
    for (uint32_t size = 2; size <= 8; size *= 2)
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    for (unsigned inputs = 0; inputs < (1U << size); ++inputs)
    for (unsigned outputs = 0; outputs < (1U << size); ++outputs)
    {
        std::vector<uint8_t> input(size), output(size);
        for (uint32_t i = 0; i < size; ++i)
        {
            input[i] = (inputs >> i) & 1U;
            output[i] = (outputs >> i) & 1U;
        }
        record_case(digest, orders[field], size, orders[field] - size,
            inverse != 0, 3, input, output);
    }
    return digest.value;
}

uint64_t plan_digest(const PrunedTransformPlan& plan)
{
    Digest digest;
    digest.add_plan(plan);
    return digest.value;
}

void allocation_failures()
{
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    for (unsigned pattern = 0; pattern < 6; ++pattern)
    {
        const std::vector<uint8_t> input = mask(64, pattern);
        const std::vector<uint8_t> output = mask(64, 1);
        Provider provider(256, 3);
        PrunedTransformPlan reference;
        allocation_audit::begin();
        const bool success = compile(64, 64, inverse != 0, input, output, provider, reference);
        const allocation_audit::Measurement measured = allocation_audit::end();
        require(success && measured.calls > 0 && measured.live > 0,
            "compiler allocation audit unavailable");
        const uint64_t expected = plan_digest(reference);
        for (size_t failure = 1; failure <= measured.calls; ++failure)
        {
            PrunedTransformPlan plan(reference);
            Provider failing(256, 3);
            allocation_audit::begin(failure);
            const bool result = compile(64, 64, inverse != 0, input, output, failing, plan);
            const allocation_audit::Measurement attempt = allocation_audit::end();
            require(attempt.failed && !result && attempt.live == 0 &&
                plan_digest(plan) == expected, "allocation failure was not atomic");
        }
        const size_t fail_calls[] = {1, provider.calls / 2, provider.calls};
        for (unsigned invalid = 0; invalid < 2; ++invalid)
        for (unsigned i = 0; i < 3; ++i)
        {
            Provider failing(256, 3);
            failing.fail_at = fail_calls[i];
            failing.invalid = invalid != 0;
            require(!compile(64, 64, inverse != 0, input, output, failing, reference) &&
                plan_digest(reference) == expected, "provider failure was not atomic");
        }
    }
}

void maximum_allocation()
{
    const std::vector<uint8_t> input = mask(65536, 1);
    const std::vector<uint8_t> output = mask(65536, 1);
    for (unsigned inverse = 0; inverse < 2; ++inverse)
    {
        Provider provider(65536, 3);
        PrunedTransformPlan plan;
        allocation_audit::begin();
        const bool result = compile(65536, 0, inverse != 0, input, output, provider, plan);
        const allocation_audit::Measurement measured = allocation_audit::end();
        require(result && measured.calls > 0 && measured.peak >= measured.live,
            "maximum compiler allocation audit failed");
        // One full mutable operation array plus bounded liveness and vector
        // growth overhead. The old raw+snapshot+planned representation exceeds
        // this budget by almost 7 MiB at the maximum GF16 side.
        const size_t transient_budget =
            (sizeof(leopard2_internal::PrunedTransformOperation) + 2U) *
            plan.full_butterfly_count;
        require(measured.peak - measured.live <= transient_budget,
            "compiler retained duplicate full transient representations");
        std::printf("N65536 inverse=%u calls=%zu peak=%zu retained=%zu transient=%zu\n",
            inverse, measured.calls, measured.peak, measured.live,
            measured.peak - measured.live);
    }
}

} // namespace

int main()
{
    try
    {
        const uint64_t grid = schedule_grid();
        const uint64_t exhaustive = exhaustive_small();
        std::printf("Schedule grid=%016llx exhaustive=%016llx\n",
            static_cast<unsigned long long>(grid),
            static_cast<unsigned long long>(exhaustive));
        allocation_failures();
        maximum_allocation();
        // Captured from the unchanged compiler at 29dd2e8, not regenerated
        // from an optimized compiler when its schedule changes.
        require(grid == UINT64_C(0x0d9fc94c6e9165d2) &&
            exhaustive == UINT64_C(0x1d25a40d577f1985),
            "compiler schedule/provider order differs from frozen baseline");
        return 0;
    }
    catch (const std::exception& error)
    {
        allocation_audit::active = false;
        std::fprintf(stderr, "Pruned compiler regression: %s\n", error.what());
        return 1;
    }
}
