/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "sparse_encode_benchmark requires LEO2_ENABLE_TEST_HOOKS"
#endif

#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

#ifndef LEO2_SPARSE_ENCODE_SOURCE_GIT_SHA
#define LEO2_SPARSE_ENCODE_SOURCE_GIT_SHA "unknown"
#endif

#ifndef LEO2_SPARSE_ENCODE_SOURCE_DIRTY
#define LEO2_SPARSE_ENCODE_SOURCE_DIRTY 1
#endif

#define LEO2_STRINGIZE_DETAIL(value) #value
#define LEO2_STRINGIZE(value) LEO2_STRINGIZE_DETAIL(value)

namespace {

typedef std::chrono::steady_clock Clock;

enum Profile
{
    ProfileHigh,
    ProfileLow
};

enum Form
{
    FormPrefix,
    FormExactPrepared,
    FormExactCallLocal,
    FormCount
};

struct Options
{
    Profile profile;
    unsigned field_bits;
    unsigned k;
    unsigned r;
    uint64_t shard_bytes;
    std::string requested_text;
    leo2_backend backend;
    unsigned iterations;
    unsigned samples;
    unsigned warmups;
    unsigned setup_iterations;
    uint64_t memory_mib;
    uint64_t seed;
    std::vector<unsigned> reuse;

    Options()
        : profile(ProfileHigh)
        , field_bits(8)
        , k(48)
        , r(16)
        , shard_bytes(1024)
        , requested_text("0,7,15")
        , backend(LEO2_BACKEND_AUTO)
        , iterations(32)
        , samples(7)
        , warmups(2)
        , setup_iterations(32)
        , memory_mib(512)
        , seed(UINT64_C(0x535041525345454e))
    {
        reuse.push_back(1);
        reuse.push_back(8);
        reuse.push_back(64);
    }
};

struct Summary
{
    double median;
    double mad;
    double minimum;
    double maximum;
    std::vector<double> samples;
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
    if (errno == ERANGE || !end || end == text || *end != '\0')
        fail(std::string("invalid ") + name);
    return static_cast<uint64_t>(value);
}

unsigned parse_unsigned32(const char* text, const char* name)
{
    const uint64_t value = parse_unsigned(text, name);
    if (value > std::numeric_limits<unsigned>::max())
        fail(std::string(name) + " exceeds uint32");
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
    if (value == "auto") return LEO2_BACKEND_AUTO;
    if (value == "scalar") return LEO2_BACKEND_SCALAR;
    if (value == "ssse3") return LEO2_BACKEND_SSSE3;
    if (value == "avx2") return LEO2_BACKEND_AVX2;
    fail("backend must be auto, scalar, ssse3, or avx2");
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
    }
    return "unknown";
}

const char* compiler_name()
{
#if defined(__clang__)
    return "clang";
#elif defined(__GNUC__)
    return "gcc";
#elif defined(_MSC_VER)
    return "msvc";
#else
    return "unknown";
#endif
}

const char* compiler_version()
{
#if defined(__clang__)
    return __clang_version__;
#elif defined(__GNUC__)
    return __VERSION__;
#elif defined(_MSC_FULL_VER)
    return LEO2_STRINGIZE(_MSC_FULL_VER);
#else
    return "unknown";
#endif
}

std::string json_escape(const char* text)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    for (const unsigned char* cursor =
             reinterpret_cast<const unsigned char*>(text);
         *cursor; ++cursor)
    {
        switch (*cursor)
        {
        case '\\': output << "\\\\"; break;
        case '"': output << "\\\""; break;
        case '\b': output << "\\b"; break;
        case '\f': output << "\\f"; break;
        case '\n': output << "\\n"; break;
        case '\r': output << "\\r"; break;
        case '\t': output << "\\t"; break;
        default:
            if (*cursor < 0x20)
                output << "\\u" << std::hex << std::setw(4)
                       << std::setfill('0') << static_cast<unsigned>(*cursor)
                       << std::dec << std::setfill(' ');
            else
                output << static_cast<char>(*cursor);
        }
    }
    return output.str();
}

void usage(const char* executable)
{
    std::cerr
        << "Usage: " << executable << " [options]\n"
        << "  --profile high|low\n"
        << "  --field gf8|gf16\n"
        << "  --k N --r N\n"
        << "  --bytes N\n"
        << "  --requested-parity MASK   comma-separated indices/ranges\n"
        << "  --backend auto|scalar|ssse3|avx2\n"
        << "  --iterations N            calls per timing sample\n"
        << "  --samples N               odd retained sample count\n"
        << "  --warmups N\n"
        << "  --setup-iterations N\n"
        << "  --reuse N[,N...]          setup amortization points\n"
        << "  --memory-mib N\n"
        << "  --seed N\n";
}

std::vector<unsigned> parse_reuse(const std::string& text)
{
    std::vector<unsigned> result;
    size_t begin = 0;
    while (begin < text.size())
    {
        const size_t comma = text.find(',', begin);
        const size_t end = comma == std::string::npos ? text.size() : comma;
        if (end == begin)
            fail("empty reuse item");
        const unsigned value = parse_unsigned32(
            text.substr(begin, end - begin).c_str(), "reuse");
        if (value == 0 || value > 1000000)
            fail("reuse values must be in 1..1000000");
        if (std::find(result.begin(), result.end(), value) != result.end())
            fail("duplicate reuse value");
        result.push_back(value);
        if (comma == std::string::npos)
            break;
        begin = comma + 1;
    }
    if (result.empty())
        fail("reuse list is empty");
    std::sort(result.begin(), result.end());
    return result;
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
        if (argument == "--profile")
        {
            const std::string profile(value);
            if (profile == "high") options.profile = ProfileHigh;
            else if (profile == "low") options.profile = ProfileLow;
            else fail("profile must be high or low");
        }
        else if (argument == "--field")
        {
            const std::string field(value);
            if (field == "gf8") options.field_bits = 8;
            else if (field == "gf16") options.field_bits = 16;
            else fail("field must be gf8 or gf16");
        }
        else if (argument == "--k")
            options.k = parse_unsigned32(value, "k");
        else if (argument == "--r")
            options.r = parse_unsigned32(value, "r");
        else if (argument == "--bytes")
            options.shard_bytes = parse_unsigned(value, "bytes");
        else if (argument == "--requested-parity" ||
                 argument == "--requested-mask")
            options.requested_text = value;
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
        else if (argument == "--reuse")
            options.reuse = parse_reuse(value);
        else if (argument == "--memory-mib")
            options.memory_mib = parse_unsigned(value, "memory-mib");
        else if (argument == "--seed")
            options.seed = parse_unsigned(value, "seed");
        else
            fail("unknown option: " + argument);
    }
    return options;
}

unsigned next_power_of_two(unsigned value)
{
    if (value == 0)
        return 0;
    --value;
    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    return value + 1;
}

size_t checked_size(uint64_t count, uint64_t bytes, const char* name)
{
    if (count && bytes > std::numeric_limits<size_t>::max() / count)
        fail(std::string(name) + " overflows size_t");
    return static_cast<size_t>(count * bytes);
}

uint64_t checked_add(uint64_t left, uint64_t right, const char* name)
{
    if (right > std::numeric_limits<uint64_t>::max() - left)
        fail(std::string(name) + " overflows uint64");
    return left + right;
}

struct Geometry
{
    unsigned order;
    unsigned side;
    unsigned block_count;
};

Geometry validate_options(const Options& options)
{
    if (options.k == 0 || options.r == 0)
        fail("K and R must be positive");
    if (options.shard_bytes == 0 ||
        options.shard_bytes > std::numeric_limits<size_t>::max())
        fail("bytes must be positive and fit size_t");
    if (options.field_bits == 16 && (options.shard_bytes & 1) != 0)
        fail("raw GF16 transform measurements require an even byte count");
    if (options.iterations == 0 || options.iterations > 1000000 ||
        options.samples < 3 || options.samples > 101 ||
        (options.samples & 1) == 0 || options.warmups > 1000 ||
        options.setup_iterations == 0 || options.setup_iterations > 1000000 ||
        options.memory_mib == 0 || options.memory_mib > 1048576)
        fail("invalid timing or memory bounds");

    Geometry geometry;
    geometry.order = options.field_bits == 8 ? 256u : 65536u;
    geometry.side = next_power_of_two(
        options.profile == ProfileHigh ? options.r : options.k);
    if (geometry.side < 2 || geometry.side > geometry.order)
        fail("padded transform side is outside the field");
    const uint64_t parent_span = options.profile == ProfileHigh
        ? static_cast<uint64_t>(options.k) + geometry.side
        : static_cast<uint64_t>(geometry.side) + options.r;
    if (parent_span > geometry.order)
        fail("profile parent does not fit the selected field");
    if (options.profile == ProfileLow && geometry.side > geometry.order / 2)
        fail("low profile requires room for at least one parity coordinate");
    geometry.block_count = options.profile == ProfileHigh
        ? 1u : (options.r + geometry.side - 1u) / geometry.side;
    return geometry;
}

std::vector<uint8_t> parse_requested(const Options& options)
{
    std::vector<uint8_t> requested(options.r, 0);
    size_t begin = 0;
    while (begin < options.requested_text.size())
    {
        const size_t comma = options.requested_text.find(',', begin);
        const size_t end = comma == std::string::npos
            ? options.requested_text.size() : comma;
        if (end == begin)
            fail("empty requested-parity item");
        const std::string item = options.requested_text.substr(begin, end - begin);
        const size_t dash = item.find('-');
        unsigned first = 0;
        unsigned last = 0;
        if (dash == std::string::npos)
            first = last = parse_unsigned32(item.c_str(), "requested parity");
        else
        {
            if (dash == 0 || dash + 1 == item.size() ||
                item.find('-', dash + 1) != std::string::npos)
                fail("invalid requested-parity range");
            first = parse_unsigned32(item.substr(0, dash).c_str(), "requested parity");
            last = parse_unsigned32(item.substr(dash + 1).c_str(), "requested parity");
            if (last < first)
                fail("descending requested-parity range");
        }
        if (last >= options.r)
            fail("requested parity is outside [0,R)");
        for (unsigned index = first; index <= last; ++index)
        {
            if (requested[index])
                fail("duplicate requested parity");
            requested[index] = 1;
        }
        if (comma == std::string::npos)
            break;
        begin = comma + 1;
    }
    if (std::count(requested.begin(), requested.end(), static_cast<uint8_t>(1)) == 0)
        fail("requested parity mask is empty");
    return requested;
}

void aligned_free(void* pointer)
{
#if defined(_MSC_VER)
    _aligned_free(pointer);
#else
    std::free(pointer);
#endif
}

void* aligned_allocate(size_t bytes)
{
#if defined(_MSC_VER)
    return _aligned_malloc(bytes, 64);
#else
    void* pointer = NULL;
    return posix_memalign(&pointer, 64, bytes) == 0 ? pointer : NULL;
#endif
}

class AlignedBuffer
{
public:
    AlignedBuffer() : data_(NULL), size_(0) {}
    ~AlignedBuffer() { aligned_free(data_); }

    void reset(size_t bytes)
    {
        aligned_free(data_);
        data_ = NULL;
        size_ = 0;
        if (bytes == 0)
            return;
        data_ = aligned_allocate(bytes);
        if (!data_)
            throw std::bad_alloc();
        size_ = bytes;
        std::memset(data_, 0, bytes);
    }

    uint8_t* bytes() { return static_cast<uint8_t*>(data_); }
    const uint8_t* bytes() const { return static_cast<const uint8_t*>(data_); }
    size_t size() const { return size_; }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* data_;
    size_t size_;
};

class Random
{
public:
    explicit Random(uint64_t seed)
        : state_(seed ? seed : UINT64_C(0x9e3779b97f4a7c15)) {}
    uint64_t next()
    {
        uint64_t value = state_;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        state_ = value;
        return value;
    }
private:
    uint64_t state_;
};

struct PlanStorage
{
    std::vector<uint8_t> masks;
    std::vector<uint8_t> workspace;
    leopard2_internal::SparseForwardPlanBatchView view;
    size_t full_butterflies;
    size_t prefix_butterflies;
    size_t retained_butterflies;
    size_t one_output_butterflies;
    size_t fused_four_groups;

    PlanStorage(const Geometry& geometry)
        : masks(geometry.block_count *
            leopard2_internal::SparseForwardRetainedBytes(geometry.side), 0)
        , workspace(
            leopard2_internal::SparseForwardDependencyBytes(geometry.side), 0)
        , full_butterflies(0)
        , prefix_butterflies(0)
        , retained_butterflies(0)
        , one_output_butterflies(0)
        , fused_four_groups(0)
    {
        view.operation_masks = masks.empty() ? NULL : &masks[0];
        view.operation_stride =
            leopard2_internal::SparseForwardRetainedBytes(geometry.side);
        view.block_count = geometry.block_count;
        if (masks.size() > 65536u || workspace.size() > 65536u - masks.size())
            fail("cell exceeds the production call-local schedule budget");
    }
};

struct GF8
{
    static bool prepare(unsigned size, unsigned shift, uint8_t* workspace,
        size_t workspace_bytes, uint8_t* masks, size_t mask_bytes,
        leopard2_internal::SparseForwardPlanStats& stats)
    {
        return leopard::ff8::PrepareSparseForwardPlan(size, shift, workspace,
            workspace_bytes, masks, mask_bytes, stats);
    }
    static void high(const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned k, unsigned prefix, unsigned requested, unsigned side,
        const void* const* data, void** work,
        const leopard2_internal::SparseForwardPlanBatchView* plan)
    {
        leopard::ff8::ReedSolomonEncode(
            ops, bytes, k, prefix, requested, side, data, work, plan);
    }
    static void low(const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned k, unsigned r, unsigned side, const void* const* data,
        void* const* recovery, void** work,
        const leopard2_internal::SparseForwardPlanBatchView* plan)
    {
        leopard::ff8::ReedSolomonEncodeLow(
            ops, bytes, k, r, side, data, recovery, work, plan);
    }
};

struct GF16
{
    static bool prepare(unsigned size, unsigned shift, uint8_t* workspace,
        size_t workspace_bytes, uint8_t* masks, size_t mask_bytes,
        leopard2_internal::SparseForwardPlanStats& stats)
    {
        return leopard::ff16::PrepareSparseForwardPlan(size, shift, workspace,
            workspace_bytes, masks, mask_bytes, stats);
    }
    static void high(const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned k, unsigned prefix, unsigned requested, unsigned side,
        const void* const* data, void** work,
        const leopard2_internal::SparseForwardPlanBatchView* plan)
    {
        leopard::ff16::ReedSolomonEncode(
            ops, bytes, k, prefix, requested, side, data, work, plan);
    }
    static void low(const leopard::backend::Ops& ops, uint64_t bytes,
        unsigned k, unsigned r, unsigned side, const void* const* data,
        void* const* recovery, void** work,
        const leopard2_internal::SparseForwardPlanBatchView* plan)
    {
        leopard::ff16::ReedSolomonEncodeLow(
            ops, bytes, k, r, side, data, recovery, work, plan);
    }
};

template<class Field>
void compile_plan(const Options& options, const Geometry& geometry,
    const std::vector<uint8_t>& requested, PlanStorage& storage,
    bool capture_structure)
{
    if (capture_structure)
    {
        storage.full_butterflies = 0;
        storage.prefix_butterflies = 0;
        storage.retained_butterflies = 0;
        storage.one_output_butterflies = 0;
        storage.fused_four_groups = 0;
    }
    for (unsigned block = 0; block < geometry.block_count; ++block)
    {
        const unsigned offset = options.profile == ProfileHigh
            ? 0u : block * geometry.side;
        const unsigned count = std::min<unsigned>(geometry.side, options.r - offset);
        unsigned prefix = 0;
        std::fill(storage.workspace.begin(), storage.workspace.end(), 0);
        for (unsigned i = 0; i < count; ++i)
        {
            if (!requested[offset + i])
                continue;
            storage.workspace[i >> 3] |= static_cast<uint8_t>(1u << (i & 7u));
            prefix = i + 1;
        }
        if (prefix == 0)
            continue;
        leopard2_internal::SparseForwardPlanStats stats;
        uint8_t* masks = &storage.masks[0] +
            static_cast<size_t>(block) * storage.view.operation_stride;
        const unsigned shift = options.profile == ProfileHigh
            ? 0u : geometry.side + offset;
        if (!Field::prepare(geometry.side, shift, &storage.workspace[0],
                storage.workspace.size(), masks,
                storage.view.operation_stride, stats))
            fail("exact sparse schedule compilation failed");
        if (capture_structure)
        {
            storage.full_butterflies += stats.full_butterfly_count;
            storage.prefix_butterflies +=
                leopard2_internal::PrefixForwardButterflyCount(
                    geometry.side, prefix);
            storage.retained_butterflies += stats.retained_butterfly_count;
            storage.one_output_butterflies += stats.one_output_butterflies;
            storage.fused_four_groups += stats.fused_four_groups;
        }
    }
    if (storage.retained_butterflies == 0)
        fail("compiled schedule retained no work");
}

struct Buffers
{
    AlignedBuffer original_storage;
    AlignedBuffer prefix_work_storage;
    AlignedBuffer exact_work_storage;
    AlignedBuffer prefix_output_storage;
    AlignedBuffer exact_output_storage;
    std::vector<const void*> original;
    std::vector<void*> prefix_work;
    std::vector<void*> exact_work;
    std::vector<void*> prefix_recovery;
    std::vector<void*> exact_recovery;

    Buffers(const Options& options, const Geometry& geometry,
        const std::vector<uint8_t>& requested)
    {
        const size_t shard_bytes = static_cast<size_t>(options.shard_bytes);
        const size_t original_bytes = checked_size(options.k, shard_bytes, "original");
        const size_t work_bytes = checked_size(
            static_cast<uint64_t>(geometry.side) * 2u, shard_bytes, "work");
        const size_t requested_count = static_cast<size_t>(std::count(
            requested.begin(), requested.end(), static_cast<uint8_t>(1)));
        const size_t output_bytes = options.profile == ProfileLow
            ? checked_size(requested_count, shard_bytes, "requested output") : 0;
        uint64_t total = static_cast<uint64_t>(original_bytes);
        total = checked_add(total, static_cast<uint64_t>(work_bytes), "memory total");
        total = checked_add(total, static_cast<uint64_t>(work_bytes), "memory total");
        total = checked_add(total, static_cast<uint64_t>(output_bytes), "memory total");
        total = checked_add(total, static_cast<uint64_t>(output_bytes), "memory total");
        const uint64_t limit = options.memory_mib * UINT64_C(1024) * UINT64_C(1024);
        if (total > limit)
            fail("cell exceeds --memory-mib");

        original_storage.reset(original_bytes);
        prefix_work_storage.reset(work_bytes);
        exact_work_storage.reset(work_bytes);
        prefix_output_storage.reset(output_bytes);
        exact_output_storage.reset(output_bytes);
        original.resize(options.k);
        prefix_work.resize(static_cast<size_t>(geometry.side) * 2u);
        exact_work.resize(static_cast<size_t>(geometry.side) * 2u);
        prefix_recovery.assign(options.r, static_cast<void*>(NULL));
        exact_recovery.assign(options.r, static_cast<void*>(NULL));

        Random random(options.seed);
        for (size_t i = 0; i < original_bytes; ++i)
            original_storage.bytes()[i] = static_cast<uint8_t>(random.next() >> 56);
        for (unsigned i = 0; i < options.k; ++i)
            original[i] = original_storage.bytes() + static_cast<size_t>(i) * shard_bytes;
        for (unsigned i = 0; i < geometry.side * 2u; ++i)
        {
            prefix_work[i] = prefix_work_storage.bytes() +
                static_cast<size_t>(i) * shard_bytes;
            exact_work[i] = exact_work_storage.bytes() +
                static_cast<size_t>(i) * shard_bytes;
        }
        size_t output_index = 0;
        for (unsigned i = 0; i < options.r; ++i)
        {
            if (!requested[i])
                continue;
            if (options.profile == ProfileLow)
            {
                prefix_recovery[i] = prefix_output_storage.bytes() +
                    output_index * shard_bytes;
                exact_recovery[i] = exact_output_storage.bytes() +
                    output_index * shard_bytes;
            }
            ++output_index;
        }
    }
};

unsigned requested_count(const std::vector<uint8_t>& requested)
{
    return static_cast<unsigned>(std::count(
        requested.begin(), requested.end(), static_cast<uint8_t>(1)));
}

unsigned requested_prefix(const std::vector<uint8_t>& requested)
{
    for (size_t i = requested.size(); i > 0; --i)
        if (requested[i - 1])
            return static_cast<unsigned>(i);
    return 0;
}

template<class Field>
void execute(const leopard::backend::Ops& ops, const Options& options,
    const Geometry& geometry, const std::vector<uint8_t>& requested,
    Buffers& buffers, const PlanStorage& plan, bool exact)
{
    if (options.profile == ProfileHigh)
    {
        Field::high(ops, options.shard_bytes, options.k,
            requested_prefix(requested), requested_count(requested), geometry.side,
            &buffers.original[0],
            exact ? &buffers.exact_work[0] : &buffers.prefix_work[0],
            exact ? &plan.view : NULL);
    }
    else
    {
        Field::low(ops, options.shard_bytes, options.k, options.r, geometry.side,
            &buffers.original[0],
            exact ? &buffers.exact_recovery[0] : &buffers.prefix_recovery[0],
            exact ? &buffers.exact_work[0] : &buffers.prefix_work[0],
            exact ? &plan.view : NULL);
    }
}

const uint8_t* output_pointer(const Options& options, const Geometry& geometry,
    const std::vector<uint8_t>& requested, const Buffers& buffers,
    unsigned recovery, bool exact)
{
    if (options.profile == ProfileHigh)
    {
        const std::vector<void*>& work = exact
            ? buffers.exact_work : buffers.prefix_work;
        return static_cast<const uint8_t*>(work[recovery]);
    }
    size_t output_index = 0;
    for (unsigned i = 0; i < recovery; ++i)
        output_index += requested[i] != 0;
    const AlignedBuffer& storage = exact
        ? buffers.exact_output_storage : buffers.prefix_output_storage;
    (void)geometry;
    return storage.bytes() + output_index * static_cast<size_t>(options.shard_bytes);
}

uint64_t verify_and_digest(const Options& options, const Geometry& geometry,
    const std::vector<uint8_t>& requested, const Buffers& buffers)
{
    uint64_t digest = UINT64_C(14695981039346656037);
    const size_t bytes = static_cast<size_t>(options.shard_bytes);
    for (unsigned recovery = 0; recovery < options.r; ++recovery)
    {
        if (!requested[recovery])
            continue;
        const uint8_t* prefix = output_pointer(
            options, geometry, requested, buffers, recovery, false);
        const uint8_t* exact = output_pointer(
            options, geometry, requested, buffers, recovery, true);
        if (std::memcmp(prefix, exact, bytes) != 0)
            fail("exact and mature-prefix parity differ");
        for (size_t i = 0; i < bytes; ++i)
        {
            digest ^= exact[i];
            digest *= UINT64_C(1099511628211);
        }
    }
    return digest;
}

double median(std::vector<double> values)
{
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    return values[middle];
}

Summary summarize(const std::vector<double>& samples)
{
    if (samples.empty())
        fail("cannot summarize empty samples");
    Summary result;
    result.samples = samples;
    result.median = median(samples);
    std::vector<double> deviations(samples.size());
    for (size_t i = 0; i < samples.size(); ++i)
        deviations[i] = std::fabs(samples[i] - result.median);
    result.mad = median(deviations);
    result.minimum = *std::min_element(samples.begin(), samples.end());
    result.maximum = *std::max_element(samples.begin(), samples.end());
    return result;
}

template<class Field>
double time_form(Form form, const leopard::backend::Ops& ops,
    const Options& options, const Geometry& geometry,
    const std::vector<uint8_t>& requested, Buffers& buffers, PlanStorage& plan)
{
    const Clock::time_point begin = Clock::now();
    for (unsigned i = 0; i < options.iterations; ++i)
    {
        if (form == FormExactCallLocal)
            compile_plan<Field>(options, geometry, requested, plan, false);
        execute<Field>(ops, options, geometry, requested, buffers, plan,
            form != FormPrefix);
    }
    const Clock::time_point end = Clock::now();
    const double elapsed = static_cast<double>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
    return elapsed / options.iterations;
}

template<class Field>
void run(const leopard::backend::Ops& ops, const Options& options,
    const Geometry& geometry, const std::vector<uint8_t>& requested)
{
    PlanStorage plan(geometry);
    compile_plan<Field>(options, geometry, requested, plan, true);
    Buffers buffers(options, geometry, requested);
    execute<Field>(ops, options, geometry, requested, buffers, plan, false);
    execute<Field>(ops, options, geometry, requested, buffers, plan, true);
    const uint64_t digest = verify_and_digest(
        options, geometry, requested, buffers);

    for (unsigned warmup = 0; warmup < options.warmups; ++warmup)
    {
        for (unsigned form = 0; form < FormCount; ++form)
        {
            if (form == FormExactCallLocal)
                compile_plan<Field>(options, geometry, requested, plan, false);
            execute<Field>(ops, options, geometry, requested, buffers, plan,
                form != FormPrefix);
        }
    }

    std::vector<double> form_samples[FormCount];
    for (unsigned sample = 0; sample < options.samples; ++sample)
    {
        for (unsigned offset = 0; offset < FormCount; ++offset)
        {
            const Form form = static_cast<Form>((sample + offset) % FormCount);
            form_samples[form].push_back(time_form<Field>(
                form, ops, options, geometry, requested, buffers, plan));
        }
    }

    std::vector<double> setup_samples;
    uint64_t setup_digest = 0;
    for (unsigned sample = 0; sample < options.samples; ++sample)
    {
        const Clock::time_point begin = Clock::now();
        for (unsigned i = 0; i < options.setup_iterations; ++i)
        {
            compile_plan<Field>(options, geometry, requested, plan, false);
            setup_digest ^= plan.retained_butterflies +
                (static_cast<uint64_t>(plan.one_output_butterflies) << 32);
        }
        const Clock::time_point end = Clock::now();
        const double elapsed = static_cast<double>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count());
        setup_samples.push_back(elapsed / options.setup_iterations);
    }
    if (setup_digest == UINT64_MAX)
        fail("unreachable setup digest");

    execute<Field>(ops, options, geometry, requested, buffers, plan, false);
    execute<Field>(ops, options, geometry, requested, buffers, plan, true);
    if (verify_and_digest(options, geometry, requested, buffers) != digest)
        fail("parity changed during timing");

    const Summary prefix = summarize(form_samples[FormPrefix]);
    const Summary exact = summarize(form_samples[FormExactPrepared]);
    const Summary call_local = summarize(form_samples[FormExactCallLocal]);
    const Summary setup = summarize(setup_samples);
    if (prefix.median <= 0 || exact.median <= 0 ||
        call_local.median <= 0 || setup.median <= 0)
        fail("clock resolution is insufficient");

    std::ostringstream digest_text;
    digest_text << "0x" << std::hex << std::setw(16) << std::setfill('0')
                << digest;
    std::cout.imbue(std::locale::classic());
    std::cout << std::fixed << std::setprecision(3)
        << "{\n"
        << "  \"schema\": \"leopard2-sparse-encode-benchmark-v1\",\n"
        << "  \"authoritative\": false,\n"
        << "  \"authority_note\": \"requires the fail-closed pinned runner "
        << "and host isolation attestation\",\n"
        << "  \"build\": {\"source_git_sha\": \""
        << LEO2_SPARSE_ENCODE_SOURCE_GIT_SHA << "\", \"source_dirty\": "
        << LEO2_SPARSE_ENCODE_SOURCE_DIRTY << ", \"compiler\": \""
        << compiler_name() << "\", \"compiler_version\": \""
        << json_escape(compiler_version()) << "\", \"cplusplus\": " << __cplusplus
        << "},\n"
        << "  \"parameters\": {\"profile\": \""
        << (options.profile == ProfileHigh ? "high" : "low")
        << "\", \"field\": \"gf" << options.field_bits
        << "\", \"K\": " << options.k << ", \"R\": " << options.r
        << ", \"shard_bytes\": " << options.shard_bytes
        << ", \"requested_parity\": [";
    bool first = true;
    for (unsigned i = 0; i < options.r; ++i)
    {
        if (!requested[i]) continue;
        if (!first) std::cout << ',';
        std::cout << i;
        first = false;
    }
    std::cout << "], \"requested_backend\": \"" << backend_name(options.backend)
        << "\", \"iterations\": " << options.iterations
        << ", \"samples\": " << options.samples
        << ", \"warmups\": " << options.warmups
        << ", \"setup_iterations\": " << options.setup_iterations
        << ", \"reuse\": [";
    for (size_t i = 0; i < options.reuse.size(); ++i)
    {
        if (i) std::cout << ',';
        std::cout << options.reuse[i];
    }
    std::cout << "], \"memory_mib\": " << options.memory_mib
        << ", \"seed\": " << options.seed << "},\n"
        << "  \"resolved\": {\"backend\": \"" << ops.name
        << "\", \"padded_side\": " << geometry.side
        << ", \"block_count\": " << geometry.block_count << "},\n"
        << "  \"plan\": {\"schedule_bytes\": " << plan.masks.size()
        << ", \"dependency_workspace_bytes\": " << plan.workspace.size()
        << ", \"full_butterflies\": " << plan.full_butterflies
        << ", \"prefix_butterflies\": " << plan.prefix_butterflies
        << ", \"retained_butterflies\": " << plan.retained_butterflies
        << ", \"one_output_butterflies\": " << plan.one_output_butterflies
        << ", \"fused_four_groups\": " << plan.fused_four_groups << "},\n"
        << "  \"correctness\": {\"exact_prefix_parity_match\": true, "
        << "\"digest_fnv1a64\": \"" << digest_text.str() << "\"},\n"
        << "  \"metrics\": {\n";

    const Summary summaries[] = { setup, prefix, exact, call_local };
    const char* names[] = {
        "schedule_setup_ns", "prefix_execution_ns",
        "exact_prepared_execution_ns", "exact_call_local_total_ns"
    };
    for (unsigned metric = 0; metric < 4; ++metric)
    {
        const Summary& value = summaries[metric];
        std::cout << "    \"" << names[metric] << "\": {\"median\": "
            << value.median << ", \"mad\": " << value.mad
            << ", \"minimum\": " << value.minimum
            << ", \"maximum\": " << value.maximum << ", \"samples\": [";
        for (size_t i = 0; i < value.samples.size(); ++i)
        {
            if (i) std::cout << ',';
            std::cout << value.samples[i];
        }
        std::cout << "]},\n";
    }
    std::cout << "    \"prefix_over_exact_prepared\": "
        << prefix.median / exact.median
        << ", \"prefix_over_exact_call_local\": "
        << prefix.median / call_local.median << ",\n"
        << "    \"amortized_exact\": [\n";
    for (size_t i = 0; i < options.reuse.size(); ++i)
    {
        const double amortized = exact.median +
            setup.median / static_cast<double>(options.reuse[i]);
        std::cout << "      {\"reuse\": " << options.reuse[i]
            << ", \"modeled_ns\": " << amortized
            << ", \"prefix_over_modeled_exact\": "
            << prefix.median / amortized << '}';
        std::cout << (i + 1 == options.reuse.size() ? "\n" : ",\n");
    }
    std::cout << "    ]\n  }\n}\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Options options = parse_options(argc, argv);
        const Geometry geometry = validate_options(options);
        const std::vector<uint8_t> requested = parse_requested(options);
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
            run<GF8>(*ops, options, geometry, requested);
        else
            run<GF16>(*ops, options, geometry, requested);
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 sparse encode benchmark: "
                  << error.what() << std::endl;
        return 1;
    }
}
