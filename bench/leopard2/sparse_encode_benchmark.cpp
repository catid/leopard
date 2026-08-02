/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

#ifndef LEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS
// Fail closed for manually compiled binaries.  The production CMake target
// sets this to zero only while linking the uninstrumented leopard archive.
#define LEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS 1
#endif

#include "Leopard2Backend.h"
#include "Leopard2Plan.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "direct_oracle.h"
#include "leopard.h"
#include "leopard2.h"

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

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#endif

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
    unsigned rounds;
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
        , rounds(3)
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

struct AbbaRound
{
    Form order[4];
    double observations[4];
    double prefix_median;
    double call_local_median;
    double log_contrast;
    double candidate_gain_percent;
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
    if (value == "avx512" || value == "avx512vl")
        return LEO2_BACKEND_AVX512;
    if (value == "gfni" || value == "avx2-gfni")
        return LEO2_BACKEND_GFNI;
    fail("backend must be auto, scalar, ssse3, avx2, avx512, or gfni");
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
    case LEO2_BACKEND_AVX512: return "avx512";
    case LEO2_BACKEND_GFNI: return "avx2-gfni";
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

std::vector<unsigned> runtime_affinity()
{
    std::vector<unsigned> result;
#if defined(__linux__)
    cpu_set_t set;
    CPU_ZERO(&set);
    if (sched_getaffinity(0, sizeof(set), &set) != 0)
        fail("cannot read runtime CPU affinity");
    for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
        if (CPU_ISSET(cpu, &set))
            result.push_back(cpu);
#endif
    return result;
}

std::string runtime_executable_path()
{
#if defined(__linux__)
    std::vector<char> buffer(4096, 0);
    const ssize_t length = readlink("/proc/self/exe", &buffer[0], buffer.size() - 1);
    if (length <= 0 || static_cast<size_t>(length) >= buffer.size())
        fail("cannot resolve runtime executable");
    buffer[static_cast<size_t>(length)] = '\0';
    return std::string(&buffer[0]);
#else
    return std::string();
#endif
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
        << "  --backend auto|scalar|ssse3|avx2|avx512|gfni\n"
        << "  --iterations N            calls per timing sample\n"
        << "  --rounds 3                independent ABBA/BAAB/ABBA rounds\n"
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
        else if (argument == "--rounds" || argument == "--samples")
            options.rounds = parse_unsigned32(value, "rounds");
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
        options.rounds != 3 ||
        options.warmups > 1000 ||
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
    AlignedBuffer() : allocation_(NULL), data_(NULL), size_(0) {}
    ~AlignedBuffer() { aligned_free(allocation_); }

    void reset(size_t bytes)
    {
        aligned_free(allocation_);
        allocation_ = NULL;
        data_ = NULL;
        size_ = 0;
        if (bytes == 0)
            return;
        if (bytes > std::numeric_limits<size_t>::max() - 128u)
            fail("guarded allocation overflows size_t");
        allocation_ = aligned_allocate(bytes + 128u);
        if (!allocation_)
            throw std::bad_alloc();
        data_ = static_cast<uint8_t*>(allocation_) + 64u;
        size_ = bytes;
        std::memset(allocation_, 0xd3, 64u);
        std::memset(data_, 0, bytes);
        std::memset(static_cast<uint8_t*>(data_) + bytes, 0xd3, 64u);
    }

    uint8_t* bytes() { return static_cast<uint8_t*>(data_); }
    const uint8_t* bytes() const { return static_cast<const uint8_t*>(data_); }
    size_t size() const { return size_; }

    void verify_guards(const char* name) const
    {
        if (!allocation_)
            return;
        const uint8_t* prefix = static_cast<const uint8_t*>(allocation_);
        const uint8_t* suffix = static_cast<const uint8_t*>(data_) + size_;
        for (size_t i = 0; i < 64u; ++i)
            if (prefix[i] != 0xd3 || suffix[i] != 0xd3)
                fail(std::string(name) + " allocation guard changed");
    }

private:
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);
    void* allocation_;
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

class Sha256
{
public:
    Sha256() : total_bytes_(0), buffered_(0)
    {
        state_[0] = UINT32_C(0x6a09e667);
        state_[1] = UINT32_C(0xbb67ae85);
        state_[2] = UINT32_C(0x3c6ef372);
        state_[3] = UINT32_C(0xa54ff53a);
        state_[4] = UINT32_C(0x510e527f);
        state_[5] = UINT32_C(0x9b05688c);
        state_[6] = UINT32_C(0x1f83d9ab);
        state_[7] = UINT32_C(0x5be0cd19);
    }

    void update(const void* input, size_t bytes)
    {
        const uint8_t* data = static_cast<const uint8_t*>(input);
        total_bytes_ += bytes;
        while (bytes != 0)
        {
            const size_t count = std::min<size_t>(bytes, 64u - buffered_);
            std::memcpy(block_ + buffered_, data, count);
            buffered_ += count;
            data += count;
            bytes -= count;
            if (buffered_ == 64u)
            {
                transform(block_);
                buffered_ = 0;
            }
        }
    }

    std::string finish()
    {
        const uint64_t bit_count = total_bytes_ * UINT64_C(8);
        block_[buffered_++] = 0x80;
        if (buffered_ > 56u)
        {
            std::memset(block_ + buffered_, 0, 64u - buffered_);
            transform(block_);
            buffered_ = 0;
        }
        std::memset(block_ + buffered_, 0, 56u - buffered_);
        for (unsigned i = 0; i < 8; ++i)
            block_[56u + i] = static_cast<uint8_t>(bit_count >> (56u - 8u * i));
        transform(block_);
        std::ostringstream text;
        text << std::hex << std::setfill('0');
        for (unsigned i = 0; i < 8; ++i)
            text << std::setw(8) << state_[i];
        return text.str();
    }

private:
    static uint32_t rotate_right(uint32_t value, unsigned count)
    {
        return (value >> count) | (value << (32u - count));
    }

    void transform(const uint8_t* input)
    {
        static const uint32_t constants[64] = {
            0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,
            0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
            0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
            0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
            0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,
            0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
            0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,
            0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
            0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
            0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
            0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,
            0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
            0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,
            0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
            0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
            0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
        };
        uint32_t words[64];
        for (unsigned i = 0; i < 16; ++i)
            words[i] = (static_cast<uint32_t>(input[i * 4]) << 24) |
                (static_cast<uint32_t>(input[i * 4 + 1]) << 16) |
                (static_cast<uint32_t>(input[i * 4 + 2]) << 8) |
                static_cast<uint32_t>(input[i * 4 + 3]);
        for (unsigned i = 16; i < 64; ++i)
        {
            const uint32_t s0 = rotate_right(words[i - 15], 7) ^
                rotate_right(words[i - 15], 18) ^ (words[i - 15] >> 3);
            const uint32_t s1 = rotate_right(words[i - 2], 17) ^
                rotate_right(words[i - 2], 19) ^ (words[i - 2] >> 10);
            words[i] = words[i - 16] + s0 + words[i - 7] + s1;
        }
        uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
        uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
        for (unsigned i = 0; i < 64; ++i)
        {
            const uint32_t upper = rotate_right(e, 6) ^ rotate_right(e, 11) ^
                rotate_right(e, 25);
            const uint32_t choose = (e & f) ^ ((~e) & g);
            const uint32_t first = h + upper + choose + constants[i] + words[i];
            const uint32_t lower = rotate_right(a, 2) ^ rotate_right(a, 13) ^
                rotate_right(a, 22);
            const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t second = lower + majority;
            h = g; g = f; f = e; e = d + first;
            d = c; c = b; b = a; a = first + second;
        }
        state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
        state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
    }

    uint32_t state_[8];
    uint64_t total_bytes_;
    uint8_t block_[64];
    size_t buffered_;
};

std::string sha256(const void* data, size_t bytes)
{
    Sha256 hash;
    hash.update(data, bytes);
    return hash.finish();
}

void verify_sha256_implementation()
{
    static const char empty[] = "";
    static const char abc[] = "abc";
    if (sha256(empty, 0) !=
            "e3b0c44298fc1c149afbf4c8996fb924"
            "27ae41e4649b934ca495991b7852b855" ||
        sha256(abc, 3) !=
            "ba7816bf8f01cfea414140de5dae2223"
            "b00361a396177a9cb410ff61f20015ad")
        fail("internal SHA-256 known-answer test failed");
}

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
            ops, bytes, k, prefix, requested, side, data, work, plan,
            bytes == 64U);
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
    std::vector<uint8_t> original_snapshot;
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
        const size_t output_bytes = options.profile == ProfileLow
            ? checked_size(options.r, shard_bytes, "parity output") : 0;
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
        if (output_bytes != 0)
        {
            std::memset(prefix_output_storage.bytes(), 0xa5, output_bytes);
            std::memset(exact_output_storage.bytes(), 0xa5, output_bytes);
        }
        original.resize(options.k);
        prefix_work.resize(static_cast<size_t>(geometry.side) * 2u);
        exact_work.resize(static_cast<size_t>(geometry.side) * 2u);
        prefix_recovery.assign(options.r, static_cast<void*>(NULL));
        exact_recovery.assign(options.r, static_cast<void*>(NULL));

        Random random(options.seed);
        for (size_t i = 0; i < original_bytes; ++i)
            original_storage.bytes()[i] = static_cast<uint8_t>(random.next() >> 56);
        original_snapshot.assign(
            original_storage.bytes(), original_storage.bytes() + original_bytes);
        for (unsigned i = 0; i < options.k; ++i)
            original[i] = original_storage.bytes() + static_cast<size_t>(i) * shard_bytes;
        for (unsigned i = 0; i < geometry.side * 2u; ++i)
        {
            prefix_work[i] = prefix_work_storage.bytes() +
                static_cast<size_t>(i) * shard_bytes;
            exact_work[i] = exact_work_storage.bytes() +
                static_cast<size_t>(i) * shard_bytes;
        }
        for (unsigned i = 0; i < options.r; ++i)
        {
            if (!requested[i])
                continue;
            if (options.profile == ProfileLow)
            {
                prefix_recovery[i] = prefix_output_storage.bytes() +
                    static_cast<size_t>(i) * shard_bytes;
                exact_recovery[i] = exact_output_storage.bytes() +
                    static_cast<size_t>(i) * shard_bytes;
            }
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
    const AlignedBuffer& storage = exact
        ? buffers.exact_output_storage : buffers.prefix_output_storage;
    (void)geometry;
    (void)requested;
    return storage.bytes() + static_cast<size_t>(recovery) *
        static_cast<size_t>(options.shard_bytes);
}

void verify_buffer_integrity(const Options& options,
    const std::vector<uint8_t>& requested, const Buffers& buffers)
{
    if (buffers.original_snapshot.size() != buffers.original_storage.size() ||
        std::memcmp(&buffers.original_snapshot[0], buffers.original_storage.bytes(),
            buffers.original_snapshot.size()) != 0)
        fail("encoder modified immutable source input");
    buffers.original_storage.verify_guards("original");
    buffers.prefix_work_storage.verify_guards("prefix work");
    buffers.exact_work_storage.verify_guards("exact work");
    buffers.prefix_output_storage.verify_guards("prefix output");
    buffers.exact_output_storage.verify_guards("exact output");
    if (options.profile != ProfileLow)
        return;
    const size_t bytes = static_cast<size_t>(options.shard_bytes);
    for (unsigned recovery = 0; recovery < options.r; ++recovery)
    {
        if (requested[recovery])
            continue;
        const size_t offset = static_cast<size_t>(recovery) * bytes;
        for (size_t i = 0; i < bytes; ++i)
            if (buffers.prefix_output_storage.bytes()[offset + i] != 0xa5 ||
                buffers.exact_output_storage.bytes()[offset + i] != 0xa5)
                fail("encoder modified an unrequested low-profile parity canary");
    }
}

typedef leopard2_test::Element OracleElement;
typedef std::vector<std::vector<OracleElement> > OracleRows;

OracleElement read_symbol(
    const uint8_t* shard, unsigned field_bits, size_t bytes, size_t symbol)
{
    if (field_bits == 8)
        return shard[symbol];
    const size_t total_symbols = bytes / 2u;
    const size_t tile_symbol = symbol & 31u;
    const size_t tile_first = symbol - tile_symbol;
    const size_t tile_symbols = std::min<size_t>(32u, total_symbols - tile_first);
    const size_t tile_byte = tile_first * 2u;
    return static_cast<OracleElement>(shard[tile_byte + tile_symbol] |
        (static_cast<unsigned>(shard[tile_byte + tile_symbols + tile_symbol]) << 8));
}

OracleRows independent_generator_rows(const Options& options,
    const std::vector<uint8_t>& requested,
    const leopard2_test::BinaryField& field)
{
    const leopard2_test::ProfileLayout layout =
        leopard2_test::make_profile_layout(
            options.profile == ProfileHigh
                ? leopard2_test::kLegacyHigh : leopard2_test::kLow,
            options.k, options.r);
    if (layout.parent_size > field.order())
        fail("oracle profile exceeds field order");
    std::vector<OracleElement> weights(options.k, 0);
    for (unsigned original = 0; original < options.k; ++original)
    {
        const OracleElement x = static_cast<OracleElement>(
            layout.systematic_coordinates[original]);
        OracleElement denominator = 1;
        for (unsigned other = 0; other < layout.parent_dimension; ++other)
        {
            const OracleElement y = static_cast<OracleElement>(
                layout.systematic_coordinates[other]);
            if (x != y)
                denominator = field.multiply(denominator, field.add(x, y));
        }
        weights[original] = field.inverse(denominator);
    }
    OracleRows rows;
    for (unsigned parity = 0; parity < options.r; ++parity)
    {
        if (!requested[parity])
            continue;
        const OracleElement point = static_cast<OracleElement>(
            layout.parity_coordinates[parity]);
        OracleElement vanishing = 1;
        for (unsigned systematic = 0;
             systematic < layout.parent_dimension; ++systematic)
        {
            vanishing = field.multiply(vanishing, field.add(point,
                static_cast<OracleElement>(
                    layout.systematic_coordinates[systematic])));
        }
        rows.push_back(std::vector<OracleElement>(options.k, 0));
        for (unsigned original = 0; original < options.k; ++original)
        {
            const OracleElement difference = field.add(point,
                static_cast<OracleElement>(
                    layout.systematic_coordinates[original]));
            rows.back()[original] = field.multiply(
                field.multiply(vanishing, field.inverse(difference)),
                weights[original]);
        }
    }
    return rows;
}

struct CorrectnessResult
{
    unsigned direct_generator_symbols;
    std::string direct_generator_sha256;
    std::string input_sha256;
    std::string parity_sha256;
    std::string recovery_source_sha256;
    std::string recovered_original_sha256;
};

void require_leo2(leo2_result result, const char* operation)
{
    if (result == LEO2_SUCCESS)
        return;
    std::ostringstream stream;
    stream << operation << ": " << leo2_result_string(result);
    fail(stream.str());
}

CorrectnessResult verify_independent_generator_and_recovery(
    const Options& options, const Geometry& geometry,
    const std::vector<uint8_t>& requested, Buffers& buffers)
{
    const leopard2_test::BinaryField field = options.field_bits == 8
        ? leopard2_test::make_legacy_gf8()
        : leopard2_test::make_legacy_gf16();
    const OracleRows rows = independent_generator_rows(options, requested, field);
    const size_t bytes = static_cast<size_t>(options.shard_bytes);
    const size_t symbols = options.field_bits == 8 ? bytes : bytes / 2u;
    const size_t sample_count = std::min<size_t>(symbols, 8u);
    Sha256 direct_digest;
    unsigned row_index = 0;
    unsigned checked = 0;
    for (unsigned parity = 0; parity < options.r; ++parity)
    {
        if (!requested[parity])
            continue;
        const uint8_t* actual = output_pointer(
            options, geometry, requested, buffers, parity, true);
        for (size_t sample = 0; sample < sample_count; ++sample)
        {
            const size_t symbol = sample_count == 1 ? 0 :
                sample * (symbols - 1u) / (sample_count - 1u);
            OracleElement expected = 0;
            for (unsigned original = 0; original < options.k; ++original)
            {
                expected ^= field.multiply(rows[row_index][original],
                    read_symbol(static_cast<const uint8_t*>(buffers.original[original]),
                        options.field_bits, bytes, symbol));
            }
            if (read_symbol(actual, options.field_bits, bytes, symbol) != expected)
                fail("exact parity differs from independent direct generator");
            const uint8_t encoded[2] = {
                static_cast<uint8_t>(expected),
                static_cast<uint8_t>(expected >> 8)
            };
            direct_digest.update(encoded, options.field_bits == 8 ? 1u : 2u);
            ++checked;
        }
        ++row_index;
    }

    Sha256 input_digest;
    input_digest.update(buffers.original_storage.bytes(),
        buffers.original_storage.size());
    Sha256 parity_digest;
    for (unsigned parity = 0; parity < options.r; ++parity)
        if (requested[parity])
            parity_digest.update(output_pointer(
                options, geometry, requested, buffers, parity, true), bytes);

    unsigned selected_parity = 0;
    while (selected_parity < options.r && !requested[selected_parity])
        ++selected_parity;
    if (selected_parity == options.r)
        fail("no requested parity is available for recovery");
    const unsigned missing = options.k / 2u;
    leo2_context_options context_options = {};
    context_options.struct_size = sizeof(context_options);
    context_options.backend = options.backend;
    context_options.thread_count = 1;
    leo2_context* context = NULL;
    require_leo2(leo2_context_create(&context_options, &context), "context create");
    leo2_codec_options codec_options = {};
    codec_options.struct_size = sizeof(codec_options);
    codec_options.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    leo2_codec* codec = NULL;
    require_leo2(leo2_codec_create(context, options.k, options.r,
        options.profile == ProfileHigh
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1,
        options.field_bits == 8 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16,
        &codec_options, &codec), "codec create");
    std::vector<uint8_t> original_present(options.k, 1);
    std::vector<uint8_t> recovery_present(options.r, 0);
    original_present[missing] = 0;
    recovery_present[selected_parity] = 1;
    leo2_decode_plan* decode_plan = NULL;
    require_leo2(leo2_decode_plan_create(codec, &original_present[0],
        &recovery_present[0], &decode_plan), "decode plan create");
    size_t scratch_bytes = 0;
    require_leo2(leo2_decode_plan_scratch_size(
        decode_plan, options.shard_bytes, &scratch_bytes), "decode scratch query");
    AlignedBuffer scratch;
    AlignedBuffer restored;
    scratch.reset(scratch_bytes);
    restored.reset(bytes);
    std::memset(restored.bytes(), 0x6c, bytes);
    std::vector<const void*> decode_original = buffers.original;
    std::vector<const void*> decode_recovery(options.r, static_cast<const void*>(NULL));
    std::vector<void*> restored_original(options.k, static_cast<void*>(NULL));
    decode_original[missing] = NULL;
    decode_recovery[selected_parity] = output_pointer(
        options, geometry, requested, buffers, selected_parity, true);
    restored_original[missing] = restored.bytes();
    require_leo2(leo2_decode_plan_execute(decode_plan, options.shard_bytes,
        &decode_original[0], &decode_recovery[0], &restored_original[0],
        scratch.bytes(), scratch.size()), "decode execute");
    if (std::memcmp(restored.bytes(), buffers.original[missing], bytes) != 0)
        fail("candidate parity did not recover the missing original");
    scratch.verify_guards("decode scratch");
    restored.verify_guards("restored original");
    const std::string restored_digest = sha256(restored.bytes(), bytes);
    const std::string source_digest = sha256(buffers.original[missing], bytes);
    if (restored_digest != source_digest)
        fail("recovered-original SHA-256 differs from source");
    leo2_decode_plan_destroy(decode_plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    verify_buffer_integrity(options, requested, buffers);

    CorrectnessResult result;
    result.direct_generator_symbols = checked;
    result.direct_generator_sha256 = direct_digest.finish();
    result.input_sha256 = input_digest.finish();
    result.parity_sha256 = parity_digest.finish();
    result.recovery_source_sha256 = source_digest;
    result.recovered_original_sha256 = restored_digest;
    return result;
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
    const double upper = values[middle];
    if ((values.size() & 1u) != 0)
        return upper;
    std::nth_element(
        values.begin(), values.begin() + middle - 1, values.begin() + middle);
    return (values[middle - 1] + upper) * 0.5;
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
    verify_buffer_integrity(options, requested, buffers);
    const CorrectnessResult correctness =
        verify_independent_generator_and_recovery(
            options, geometry, requested, buffers);

    for (unsigned warmup = 0; warmup < options.warmups; ++warmup)
    {
        time_form<Field>(FormPrefix, ops, options, geometry,
            requested, buffers, plan);
        time_form<Field>(FormExactCallLocal, ops, options, geometry,
            requested, buffers, plan);
        time_form<Field>(FormExactPrepared, ops, options, geometry,
            requested, buffers, plan);
    }

    // Promotion inference uses only the setup-inclusive pair.  Each independent
    // round is drift-balanced, alternating A/B/B/A with B/A/A/B; prepared
    // execution is measured later in a separate diagnostic block and cannot
    // influence the primary comparison.
    std::vector<AbbaRound> rounds;
    std::vector<double> prefix_round_medians;
    std::vector<double> call_local_round_medians;
    std::vector<double> round_gains;
    rounds.reserve(options.rounds);
    for (unsigned round_index = 0; round_index < options.rounds; ++round_index)
    {
        AbbaRound round;
        if ((round_index & 1u) == 0)
        {
            round.order[0] = FormPrefix;
            round.order[1] = FormExactCallLocal;
            round.order[2] = FormExactCallLocal;
            round.order[3] = FormPrefix;
        }
        else
        {
            round.order[0] = FormExactCallLocal;
            round.order[1] = FormPrefix;
            round.order[2] = FormPrefix;
            round.order[3] = FormExactCallLocal;
        }
        std::vector<double> prefix_values;
        std::vector<double> candidate_values;
        for (unsigned position = 0; position < 4; ++position)
        {
            round.observations[position] = time_form<Field>(
                round.order[position], ops, options, geometry,
                requested, buffers, plan);
            if (round.order[position] == FormPrefix)
                prefix_values.push_back(round.observations[position]);
            else
                candidate_values.push_back(round.observations[position]);
        }
        round.prefix_median = median(prefix_values);
        round.call_local_median = median(candidate_values);
        round.log_contrast = 0.5 * (
            std::log(prefix_values[0]) + std::log(prefix_values[1]) -
            std::log(candidate_values[0]) - std::log(candidate_values[1]));
        round.candidate_gain_percent =
            (std::exp(round.log_contrast) - 1.0) * 100.0;
        rounds.push_back(round);
        prefix_round_medians.push_back(round.prefix_median);
        call_local_round_medians.push_back(round.call_local_median);
        round_gains.push_back(round.candidate_gain_percent);
    }

    std::vector<double> prepared_samples;
    for (unsigned round_index = 0; round_index < options.rounds; ++round_index)
        prepared_samples.push_back(time_form<Field>(FormExactPrepared,
            ops, options, geometry, requested, buffers, plan));

    std::vector<double> setup_samples;
    uint64_t setup_digest = 0;
    for (unsigned round_index = 0; round_index < options.rounds; ++round_index)
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
    verify_buffer_integrity(options, requested, buffers);
    const CorrectnessResult after_correctness =
        verify_independent_generator_and_recovery(
            options, geometry, requested, buffers);
    if (after_correctness.direct_generator_symbols !=
            correctness.direct_generator_symbols ||
        after_correctness.direct_generator_sha256 !=
            correctness.direct_generator_sha256 ||
        after_correctness.input_sha256 != correctness.input_sha256 ||
        after_correctness.parity_sha256 != correctness.parity_sha256 ||
        after_correctness.recovery_source_sha256 !=
            correctness.recovery_source_sha256 ||
        after_correctness.recovered_original_sha256 !=
            correctness.recovered_original_sha256)
        fail("correctness digests changed during timing");

    const Summary prefix = summarize(prefix_round_medians);
    const Summary exact = summarize(prepared_samples);
    const Summary call_local = summarize(call_local_round_medians);
    const Summary gains = summarize(round_gains);
    const Summary setup = summarize(setup_samples);
    if (prefix.median <= 0 || exact.median <= 0 ||
        call_local.median <= 0 || setup.median <= 0)
        fail("clock resolution is insufficient");

    std::ostringstream digest_text;
    digest_text << "0x" << std::hex << std::setw(16) << std::setfill('0')
                << digest;
    const std::vector<unsigned> affinity = runtime_affinity();
    const std::string executable_path = runtime_executable_path();
    std::cout.imbue(std::locale::classic());
    std::cout << std::fixed << std::setprecision(6)
        << "{\n"
        << "  \"schema\": \"leopard2-sparse-encode-benchmark-v4\",\n"
        << "  \"authoritative\": false,\n"
        << "  \"authority_note\": \"requires the fail-closed pinned runner "
        << "and host isolation attestation\",\n"
        << "  \"build\": {\"source_git_sha\": \""
        << LEO2_SPARSE_ENCODE_SOURCE_GIT_SHA << "\", \"source_dirty\": "
        << LEO2_SPARSE_ENCODE_SOURCE_DIRTY
        << ", \"library_test_hooks\": "
        << (LEO2_SPARSE_ENCODE_LIBRARY_TEST_HOOKS ? "true" : "false")
        << ", \"compiler\": \""
        << compiler_name() << "\", \"compiler_version\": \""
        << json_escape(compiler_version()) << "\", \"cplusplus\": " << __cplusplus
        << "},\n"
        << "  \"runtime\": {\"linux_procfs_affinity_attested\": "
#if defined(__linux__)
        << "true"
#else
        << "false"
#endif
        << ", \"executable_path\": \"" << json_escape(executable_path.c_str())
        << "\", \"allowed_cpus\": [";
    for (size_t cpu_index = 0; cpu_index < affinity.size(); ++cpu_index)
    {
        if (cpu_index) std::cout << ',';
        std::cout << affinity[cpu_index];
    }
    std::cout << "]},\n"
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
        << ", \"rounds\": " << options.rounds
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
        << "\"direct_generator_parity_match\": true, "
        << "\"direct_generator_symbols_checked\": "
        << correctness.direct_generator_symbols
        << ", \"direct_generator_sample_sha256\": \""
        << correctness.direct_generator_sha256
        << "\", \"encode_decode_recovery_match\": true, "
        << "\"input_immutable\": true, \"allocation_guards_match\": true, "
        << "\"unrequested_output_canary\": {\"applicable\": "
        << (options.profile == ProfileLow ? "true" : "false")
        << ", \"match\": true}, \"post_timing_recheck_match\": true, "
        << "\"input_sha256\": \"" << correctness.input_sha256
        << "\", \"parity_sha256\": \"" << correctness.parity_sha256
        << "\", \"recovery_source_sha256\": \""
        << correctness.recovery_source_sha256
        << "\", \"recovered_original_sha256\": \""
        << correctness.recovered_original_sha256
        << "\", \"digest_fnv1a64\": \"" << digest_text.str() << "\"},\n"
        << "  \"primary_abba\": {\n"
        << "    \"design\": \"pairwise_prefix_call_local_ABBA\",\n"
        << "    \"order_policy\": \"alternating_ABBA_BAAB\",\n"
        << "    \"rounds\": [\n";
    for (size_t round_index = 0; round_index < rounds.size(); ++round_index)
    {
        const AbbaRound& round = rounds[round_index];
        std::cout << "      {\"round\": " << round_index
            << ", \"order\": [";
        for (unsigned position = 0; position < 4; ++position)
        {
            if (position) std::cout << ',';
            std::cout << (round.order[position] == FormPrefix
                ? "\"prefix\"" : "\"call_local\"");
        }
        std::cout << "], \"observations_ns\": ["
            << round.observations[0] << ',' << round.observations[1] << ','
            << round.observations[2] << ',' << round.observations[3]
            << "], \"prefix_median_ns\": " << round.prefix_median
            << ", \"call_local_median_ns\": " << round.call_local_median
            << ", \"log_contrast\": " << round.log_contrast
            << ", \"candidate_gain_percent\": "
            << round.candidate_gain_percent << '}';
        std::cout << (round_index + 1 == rounds.size() ? "\n" : ",\n");
    }
    std::cout << "    ],\n"
        << "    \"round_gain_percent\": {\"median\": " << gains.median
        << ", \"mad\": " << gains.mad
        << ", \"minimum\": " << gains.minimum
        << ", \"maximum\": " << gains.maximum << ", \"samples\": [";
    for (size_t i = 0; i < gains.samples.size(); ++i)
    {
        if (i) std::cout << ',';
        std::cout << gains.samples[i];
    }
    std::cout << "]}\n  },\n"
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
        verify_sha256_implementation();
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
