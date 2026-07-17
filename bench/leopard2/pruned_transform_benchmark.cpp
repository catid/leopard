/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "pruned_transform_benchmark requires LEO2_ENABLE_TEST_HOOKS"
#endif

#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef LEO2_PRUNED_SOURCE_GIT_SHA
#define LEO2_PRUNED_SOURCE_GIT_SHA "unknown"
#endif

#ifndef LEO2_PRUNED_SOURCE_DIRTY
#define LEO2_PRUNED_SOURCE_DIRTY 1
#endif

namespace {

typedef std::chrono::steady_clock Clock;

enum MaskShape
{
    MaskPrefix,
    MaskHoley
};

struct Options
{
    unsigned field_bits;
    unsigned size;
    unsigned shift;
    bool inverse;
    uint64_t shard_bytes;
    unsigned input_active;
    unsigned output_requested;
    MaskShape input_shape;
    MaskShape output_shape;
    unsigned full_input_prefix;
    unsigned full_output_prefix;
    leo2_backend backend;
    unsigned iterations;
    unsigned samples;
    unsigned warmups;
    unsigned setup_iterations;
    uint64_t memory_mib;
    uint64_t seed;

    Options()
        : field_bits(8)
        , size(64)
        , shift(0)
        , inverse(false)
        , shard_bytes(1024)
        , input_active(33)
        , output_requested(17)
        , input_shape(MaskPrefix)
        , output_shape(MaskPrefix)
        , full_input_prefix(0)
        , full_output_prefix(0)
        , backend(LEO2_BACKEND_AUTO)
        , iterations(64)
        , samples(7)
        , warmups(2)
        , setup_iterations(16)
        , memory_mib(256)
        , seed(UINT64_C(0x433142454e43484d))
    {}
};

void fail(const std::string& message)
{
    throw std::runtime_error(message);
}

uint64_t parse_unsigned(const char* text, const char* name)
{
    if (!text || text[0] == '\0' || text[0] == '-')
        fail(std::string("invalid ") + name);
    errno = 0;
    char* end = NULL;
    const unsigned long long value = std::strtoull(text, &end, 0);
    if (errno == ERANGE || !end || *end != '\0')
        fail(std::string("invalid ") + name);
    return static_cast<uint64_t>(value);
}

unsigned parse_unsigned32(const char* text, const char* name)
{
    const uint64_t value = parse_unsigned(text, name);
    if (value > std::numeric_limits<unsigned>::max())
        fail(std::string(name) + " exceeds the 32-bit range");
    return static_cast<unsigned>(value);
}

const char* require_value(int& index, int argc, char** argv)
{
    if (++index >= argc)
        fail(std::string("missing value for ") + argv[index - 1]);
    return argv[index];
}

leo2_backend parse_backend(const std::string& value)
{
    if (value == "auto")
        return LEO2_BACKEND_AUTO;
    if (value == "scalar")
        return LEO2_BACKEND_SCALAR;
    if (value == "ssse3")
        return LEO2_BACKEND_SSSE3;
    if (value == "avx2")
        return LEO2_BACKEND_AVX2;
    fail("invalid backend");
    return LEO2_BACKEND_AUTO;
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_AUTO: return "auto";
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    default: return "unknown";
    }
}

MaskShape parse_mask_shape(const std::string& value)
{
    if (value == "prefix")
        return MaskPrefix;
    if (value == "holey")
        return MaskHoley;
    fail("invalid mask shape");
    return MaskPrefix;
}

const char* mask_shape_name(MaskShape shape)
{
    return shape == MaskPrefix ? "prefix" : "holey";
}

void usage(const char* executable)
{
    std::cerr
        << "Usage: " << executable << " [options]\n"
        << "  --field gf8|gf16\n"
        << "  --size N                  dyadic transform size\n"
        << "  --shift N                 aligned additive-coset shift\n"
        << "  --direction forward|inverse\n"
        << "  --bytes N                 physical bytes per shard\n"
        << "  --input-active N          active input prefix\n"
        << "  --output-requested N      requested output prefix\n"
        << "  --input-shape prefix|holey\n"
        << "  --output-shape prefix|holey\n"
        << "  --backend auto|scalar|ssse3|avx2\n"
        << "  --iterations N            executions per timing sample\n"
        << "  --samples N               retained timing samples\n"
        << "  --warmups N\n"
        << "  --setup-iterations N\n"
        << "  --memory-mib N            workspace cap\n"
        << "  --seed N\n";
}

Options parse_options(int argc, char** argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if (argument == "--help" || argument == "-h")
        {
            usage(argv[0]);
            std::exit(0);
        }
        const char* value = require_value(i, argc, argv);
        if (argument == "--field")
        {
            const std::string field(value);
            if (field == "gf8")
                options.field_bits = 8;
            else if (field == "gf16")
                options.field_bits = 16;
            else
                fail("invalid field");
        }
        else if (argument == "--size")
            options.size = parse_unsigned32(value, "size");
        else if (argument == "--shift")
            options.shift = parse_unsigned32(value, "shift");
        else if (argument == "--direction")
        {
            const std::string direction(value);
            if (direction == "forward")
                options.inverse = false;
            else if (direction == "inverse")
                options.inverse = true;
            else
                fail("invalid direction");
        }
        else if (argument == "--bytes")
            options.shard_bytes = parse_unsigned(value, "bytes");
        else if (argument == "--input-active")
            options.input_active = parse_unsigned32(value, "input-active");
        else if (argument == "--output-requested")
            options.output_requested = parse_unsigned32(value, "output-requested");
        else if (argument == "--input-shape")
            options.input_shape = parse_mask_shape(value);
        else if (argument == "--output-shape")
            options.output_shape = parse_mask_shape(value);
        else if (argument == "--backend")
            options.backend = parse_backend(value);
        else if (argument == "--iterations")
            options.iterations = parse_unsigned32(value, "iterations");
        else if (argument == "--samples")
            options.samples = parse_unsigned32(value, "samples");
        else if (argument == "--warmups")
            options.warmups = parse_unsigned32(value, "warmups");
        else if (argument == "--setup-iterations")
            options.setup_iterations = parse_unsigned32(value, "setup-iterations");
        else if (argument == "--memory-mib")
            options.memory_mib = parse_unsigned(value, "memory-mib");
        else if (argument == "--seed")
            options.seed = parse_unsigned(value, "seed");
        else
            fail(std::string("unknown option: ") + argument);
    }
    return options;
}

bool is_power_of_two(unsigned value)
{
    return value >= 2 && (value & (value - 1)) == 0;
}

void validate_options(const Options& options)
{
    const unsigned order = options.field_bits == 8 ? 256u : 65536u;
    if (!is_power_of_two(options.size) || options.size > order)
        fail("size must be a supported power of two >= 2");
    if ((options.shift & (options.size - 1)) != 0 ||
        options.shift > order - options.size)
        fail("shift must name an aligned in-field coset");
    if (options.shard_bytes == 0 ||
        options.shard_bytes > std::numeric_limits<size_t>::max())
        fail("bytes is out of range");
    if (options.field_bits == 16 && (options.shard_bytes & 1) != 0)
        fail("GF16 requires an even physical byte count");
    if (options.input_active == 0 || options.input_active > options.size)
        fail("input-active must be in 1..size");
    if (options.output_requested == 0 ||
        options.output_requested > options.size)
        fail("output-requested must be in 1..size");
    if (options.iterations == 0 || options.iterations > 1000000 ||
        options.samples < 3 || (options.samples & 1u) == 0 ||
        options.samples > 101 || options.warmups > 1000 ||
        options.setup_iterations == 0 ||
        options.setup_iterations > 1000000 || options.memory_mib == 0)
        fail("timing counts or memory cap are out of range");
}

uint64_t mix64(uint64_t value)
{
    value ^= value >> 30;
    value *= UINT64_C(0xbf58476d1ce4e5b9);
    value ^= value >> 27;
    value *= UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

uint64_t checked_workspace_bytes(const Options& options)
{
    if (options.shard_bytes >
        std::numeric_limits<uint64_t>::max() / options.size)
        fail("workspace size overflow");
    const uint64_t result = options.shard_bytes * options.size;
    if (result > std::numeric_limits<size_t>::max())
        fail("workspace exceeds addressable storage");
    return result;
}

std::vector<uint8_t> make_baseline(
    const Options& options,
    const std::vector<uint8_t>& input_mask)
{
    const uint64_t total = checked_workspace_bytes(options);
    std::vector<uint8_t> result(static_cast<size_t>(total), 0);
    for (unsigned shard = 0; shard < options.size; ++shard)
    {
        if (!input_mask[shard])
            continue;
        for (uint64_t offset = 0; offset < options.shard_bytes; ++offset)
        {
            result[static_cast<size_t>(shard * options.shard_bytes + offset)] =
                static_cast<uint8_t>(mix64(
                    options.seed ^ (static_cast<uint64_t>(shard) << 40) ^
                    offset));
        }
    }
    return result;
}

std::vector<uint8_t> make_mask(
    unsigned size,
    unsigned count,
    MaskShape shape)
{
    std::vector<uint8_t> mask(size, 0);
    if (shape == MaskPrefix)
    {
        for (unsigned i = 0; i < count; ++i)
            mask[i] = 1;
    }
    else
    {
        // Every supported size is dyadic and five is odd, so this is a
        // permutation rather than a sampling loop with duplicates.
        for (unsigned i = 0; i < count; ++i)
            mask[(i * 5u + 3u) & (size - 1u)] = 1;
    }
    return mask;
}

unsigned prefix_extent(const std::vector<uint8_t>& mask)
{
    unsigned extent = static_cast<unsigned>(mask.size());
    while (extent != 0 && !mask[extent - 1])
        --extent;
    return extent;
}

class Workspace
{
public:
    Workspace(unsigned size, uint64_t bytes)
        : bytes_(static_cast<size_t>(bytes))
        , storage_(static_cast<size_t>(size) * bytes_, 0)
        , pointers_(size, NULL)
    {
        for (unsigned i = 0; i < size; ++i)
            pointers_[i] = storage_.data() + static_cast<size_t>(i) * bytes_;
    }

    void Reset(const std::vector<uint8_t>& baseline)
    {
        std::memcpy(storage_.data(), baseline.data(), baseline.size());
    }

    void** Pointers() { return pointers_.data(); }
    const std::vector<uint8_t>& Storage() const { return storage_; }

private:
    size_t bytes_;
    std::vector<uint8_t> storage_;
    std::vector<void*> pointers_;
};

uint64_t plan_storage_bytes(const leopard2_internal::PrunedTransformPlan& plan)
{
    return sizeof(plan) +
        plan.input_mask.capacity() * sizeof(plan.input_mask[0]) +
        plan.output_mask.capacity() * sizeof(plan.output_mask[0]) +
        plan.operations.capacity() * sizeof(plan.operations[0]) +
        plan.fused_four_starts.capacity() *
            sizeof(plan.fused_four_starts[0]) +
        plan.zero_outputs.capacity() * sizeof(plan.zero_outputs[0]);
}

template<class Field>
leopard2_internal::PrunedTransformPlan compile_plan(
    const Options& options,
    const std::vector<uint8_t>& input_mask,
    const std::vector<uint8_t>& output_mask)
{
    leopard2_internal::PrunedTransformPlan plan;
    if (!Field::Prepare(
            options.size, options.shift, options.inverse,
            input_mask.data(), output_mask.data(), plan))
        fail("pruned plan construction failed");
    return plan;
}

struct GF8
{
    static bool Prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }

    static void Full(
        const leopard::backend::Ops& ops,
        const Options& options,
        void** work)
    {
        if (options.inverse)
            leopard::ff8::TestOnlyLchInverse(
                ops, options.shard_bytes, options.size, options.shift,
                options.full_input_prefix, work);
        else
            leopard::ff8::TestOnlyLchForward(
                ops, options.shard_bytes, options.size, options.shift,
                options.full_output_prefix, work);
    }
};

struct GF16
{
    static bool Prepare(
        unsigned size, unsigned shift, bool inverse,
        const uint8_t* input, const uint8_t* output,
        leopard2_internal::PrunedTransformPlan& plan)
    {
        return leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input, output, plan);
    }

    static void Full(
        const leopard::backend::Ops& ops,
        const Options& options,
        void** work)
    {
        if (options.inverse)
            leopard::ff16::TestOnlyLchInverse(
                ops, options.shard_bytes, options.size, options.shift,
                options.full_input_prefix, work);
        else
            leopard::ff16::TestOnlyLchForward(
                ops, options.shard_bytes, options.size, options.shift,
                options.full_output_prefix, work);
    }
};

enum Form
{
    FormFull = 0,
    FormFlat = 1,
    FormFused = 2,
    FormCount = 3
};

template<class Field>
void execute_form(
    Form form,
    const leopard::backend::Ops& ops,
    const Options& options,
    const leopard2_internal::PrunedTransformPlan& flat_plan,
    const leopard2_internal::PrunedTransformPlan& fused_plan,
    void** work)
{
    if (form == FormFull)
    {
        Field::Full(ops, options, work);
        return;
    }
    const leopard2_internal::PrunedTransformPlan& plan =
        form == FormFlat ? flat_plan : fused_plan;
    if (!leopard2_internal::ExecutePrunedTransformPlan(
            ops, options.shard_bytes, plan, work))
        fail("pruned transform execution failed");
}

uint64_t fnv1a64(uint64_t hash, const void* data, size_t bytes)
{
    const uint8_t* input = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        hash ^= input[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

template<class Field>
uint64_t correctness_gate(
    const leopard::backend::Ops& ops,
    const Options& options,
    const std::vector<uint8_t>& baseline,
    const std::vector<uint8_t>& output_mask,
    const leopard2_internal::PrunedTransformPlan& flat_plan,
    const leopard2_internal::PrunedTransformPlan& fused_plan)
{
    Workspace full(options.size, options.shard_bytes);
    Workspace flat(options.size, options.shard_bytes);
    Workspace fused(options.size, options.shard_bytes);
    full.Reset(baseline);
    flat.Reset(baseline);
    fused.Reset(baseline);
    execute_form<Field>(FormFull, ops, options, flat_plan, fused_plan,
        full.Pointers());
    execute_form<Field>(FormFlat, ops, options, flat_plan, fused_plan,
        flat.Pointers());
    execute_form<Field>(FormFused, ops, options, flat_plan, fused_plan,
        fused.Pointers());

    uint64_t digest = UINT64_C(1469598103934665603);
    const size_t bytes = static_cast<size_t>(options.shard_bytes);
    for (unsigned i = 0; i < options.size; ++i)
    {
        if (!output_mask[i])
            continue;
        const size_t offset = static_cast<size_t>(i) * bytes;
        if (std::memcmp(
                full.Storage().data() + offset,
                flat.Storage().data() + offset, bytes) != 0 ||
            std::memcmp(
                full.Storage().data() + offset,
                fused.Storage().data() + offset, bytes) != 0)
            fail("full, flat, and fused requested outputs differ");
        digest = fnv1a64(digest, full.Storage().data() + offset, bytes);
    }
    return digest;
}

template<class Field>
double time_form(
    Form form,
    const leopard::backend::Ops& ops,
    const Options& options,
    const std::vector<uint8_t>& baseline,
    const leopard2_internal::PrunedTransformPlan& flat_plan,
    const leopard2_internal::PrunedTransformPlan& fused_plan,
    std::vector<Workspace>& workspaces)
{
    unsigned remaining = options.iterations;
    uint64_t elapsed_ns = 0;
    while (remaining != 0)
    {
        const unsigned count = std::min<unsigned>(
            remaining, static_cast<unsigned>(workspaces.size()));
        for (unsigned i = 0; i < count; ++i)
            workspaces[i].Reset(baseline);
        const Clock::time_point begin = Clock::now();
        for (unsigned i = 0; i < count; ++i)
            execute_form<Field>(
                form, ops, options, flat_plan, fused_plan,
                workspaces[i].Pointers());
        const Clock::time_point end = Clock::now();
        elapsed_ns += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - begin).count());
        remaining -= count;
    }
    return static_cast<double>(elapsed_ns) / options.iterations;
}

double median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() & 1 ? values[middle]
        : (values[middle - 1] + values[middle]) * 0.5;
}

void write_array(const std::vector<double>& values)
{
    std::cout << '[';
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0)
            std::cout << ',';
        std::cout << std::fixed << std::setprecision(3) << values[i];
    }
    std::cout << ']';
}

template<class Field>
void run(
    Options options,
    const leopard::backend::Ops& ops)
{
    const std::vector<uint8_t> input_mask = make_mask(
        options.size, options.input_active, options.input_shape);
    const std::vector<uint8_t> output_mask = make_mask(
        options.size, options.output_requested, options.output_shape);
    options.full_input_prefix = prefix_extent(input_mask);
    options.full_output_prefix = prefix_extent(output_mask);

    leopard2_internal::PrunedTransformPlan fused_plan =
        compile_plan<Field>(options, input_mask, output_mask);
    leopard2_internal::PrunedTransformPlan flat_plan(fused_plan);
    std::vector<uint32_t>().swap(flat_plan.fused_four_starts);

    const std::vector<uint8_t> baseline =
        make_baseline(options, input_mask);
    const uint64_t digest = correctness_gate<Field>(
        ops, options, baseline, output_mask, flat_plan, fused_plan);

    const uint64_t bytes_per_workspace = checked_workspace_bytes(options);
    const uint64_t memory_bytes = options.memory_mib >
        std::numeric_limits<uint64_t>::max() / (1024u * 1024u)
        ? std::numeric_limits<uint64_t>::max()
        : options.memory_mib * 1024u * 1024u;
    uint64_t instance_count = memory_bytes / bytes_per_workspace;
    if (instance_count == 0)
        instance_count = 1;
    if (instance_count > options.iterations)
        instance_count = options.iterations;
    std::vector<Workspace> workspaces;
    workspaces.reserve(static_cast<size_t>(instance_count));
    for (uint64_t i = 0; i < instance_count; ++i)
        workspaces.push_back(Workspace(options.size, options.shard_bytes));

    for (unsigned warmup = 0; warmup < options.warmups; ++warmup)
    {
        Workspace& workspace = workspaces[warmup % workspaces.size()];
        for (unsigned form = 0; form < FormCount; ++form)
        {
            workspace.Reset(baseline);
            execute_form<Field>(
                static_cast<Form>(form), ops, options,
                flat_plan, fused_plan, workspace.Pointers());
        }
    }

    std::vector<double> timing[FormCount];
    for (unsigned sample = 0; sample < options.samples; ++sample)
    {
        for (unsigned offset = 0; offset < FormCount; ++offset)
        {
            const Form form = static_cast<Form>((sample + offset) % FormCount);
            timing[form].push_back(time_form<Field>(
                form, ops, options, baseline, flat_plan, fused_plan,
                workspaces));
        }
    }

    std::vector<double> setup_samples;
    for (unsigned sample = 0; sample < options.samples; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        uint64_t setup_digest = 0;
        for (unsigned i = 0; i < options.setup_iterations; ++i)
        {
            const leopard2_internal::PrunedTransformPlan plan =
                compile_plan<Field>(options, input_mask, output_mask);
            setup_digest ^= plan.operations.size() +
                (static_cast<uint64_t>(plan.fused_four_starts.size()) << 32);
        }
        const Clock::time_point end = Clock::now();
        if (setup_digest == UINT64_MAX)
            fail("unreachable setup digest");
        const double elapsed = static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                end - begin).count());
        setup_samples.push_back(elapsed / options.setup_iterations);
    }

    const double full_median = median(timing[FormFull]);
    const double flat_median = median(timing[FormFlat]);
    const double fused_median = median(timing[FormFused]);
    if (full_median <= 0 || flat_median <= 0 || fused_median <= 0)
        fail("clock resolution was insufficient for this timing cell");
    const size_t effective_steps = fused_plan.operations.size() -
        fused_plan.fused_four_starts.size() * 3u;

    std::ostringstream digest_text;
    digest_text << "0x" << std::hex << std::setw(16) << std::setfill('0')
                << digest;

    std::cout << "{\n"
        << "  \"schema\": \"leopard2-pruned-transform-benchmark-v1\",\n"
        << "  \"authoritative\": false,\n"
        << "  \"authority_note\": \"requires an external isolated runner and host-state attestation\",\n"
        << "  \"build\": {\"source_git_sha\": \""
        << LEO2_PRUNED_SOURCE_GIT_SHA << "\", \"source_dirty\": "
        << LEO2_PRUNED_SOURCE_DIRTY << "},\n"
        << "  \"parameters\": {\"field\": \"gf" << options.field_bits
        << "\", \"size\": " << options.size
        << ", \"shift\": " << options.shift
        << ", \"direction\": \"" << (options.inverse ? "inverse" : "forward")
        << "\", \"shard_bytes\": " << options.shard_bytes
        << ", \"input_active\": " << options.input_active
        << ", \"output_requested\": " << options.output_requested
        << ", \"input_shape\": \"" << mask_shape_name(options.input_shape)
        << "\", \"output_shape\": \"" << mask_shape_name(options.output_shape)
        << "\", \"full_input_prefix\": " << options.full_input_prefix
        << ", \"full_output_prefix\": " << options.full_output_prefix
        << ", \"requested_backend\": \"" << backend_name(options.backend)
        << "\", \"iterations\": " << options.iterations
        << ", \"samples\": " << options.samples
        << ", \"warmups\": " << options.warmups
        << ", \"setup_iterations\": " << options.setup_iterations
        << ", \"memory_mib\": " << options.memory_mib
        << ", \"seed\": " << options.seed << "},\n"
        << "  \"resolved\": {\"backend\": \"" << ops.name
        << "\", \"workspace_instances\": " << instance_count << "},\n"
        << "  \"plan\": {\"full_butterflies\": "
        << fused_plan.full_butterfly_count
        << ", \"operations\": " << fused_plan.operations.size()
        << ", \"fused_four_descriptors\": "
        << fused_plan.fused_four_starts.size()
        << ", \"effective_execution_steps\": " << effective_steps
        << ", \"flat_storage_bytes\": " << plan_storage_bytes(flat_plan)
        << ", \"fused_storage_bytes\": " << plan_storage_bytes(fused_plan)
        << "},\n"
        << "  \"correctness\": {\"matches_full\": true, \"digest_fnv1a64\": \""
        << digest_text.str() << "\"},\n"
        << "  \"metrics\": {\n"
        << "    \"setup_ns\": {\"median\": " << std::fixed
        << std::setprecision(3) << median(setup_samples) << ", \"samples\": ";
    write_array(setup_samples);
    std::cout << "},\n"
        << "    \"full_ns\": {\"median\": " << full_median
        << ", \"samples\": ";
    write_array(timing[FormFull]);
    std::cout << "},\n"
        << "    \"flat_ns\": {\"median\": " << flat_median
        << ", \"samples\": ";
    write_array(timing[FormFlat]);
    std::cout << "},\n"
        << "    \"fused_ns\": {\"median\": " << fused_median
        << ", \"samples\": ";
    write_array(timing[FormFused]);
    std::cout << "},\n"
        << "    \"full_over_flat\": " << (full_median / flat_median)
        << ", \"full_over_fused\": " << (full_median / fused_median)
        << "\n  }\n}\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        validate_options(options);
        if (leo_init() != Leopard_Success)
            fail("Leopard initialization failed");

        leopard::backend::QualificationStatus status =
            leopard::backend::QualificationAvailable;
        const leopard::backend::Ops* ops =
            leopard::backend::GetQualifiedOps(options.backend, &status);
        if (!ops)
        {
            std::ostringstream stream;
            stream << "requested backend is unavailable (status="
                   << static_cast<int>(status) << ')';
            fail(stream.str());
        }

        if (options.field_bits == 8)
            run<GF8>(options, *ops);
        else
            run<GF16>(options, *ops);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 pruned transform benchmark: "
                  << error.what() << std::endl;
        return 1;
    }
}
