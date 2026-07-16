/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

// C5 research checkpoint: canonical binary dyadic-block execution.
//
// This standalone, test-hook-only program is deliberately absent from CMake
// and from the production dispatcher.  It evaluates a nonzero prefix of the
// existing padded parent's normalized LCH coefficients by splitting the prefix
// into aligned dyadic blocks.  Complete blocks use Leopard's production LCH
// kernels; block factors are derived independently in the legacy field.  Every
// result is compared with the complete padded production transform, while
// sampled symbols are also compared with direct normalized-LCH evaluation.

#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "leopard2.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <sys/utsname.h>
#include <unistd.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#if !defined(LEO2_C5_SOURCE_SHA256)
#define LEO2_C5_SOURCE_SHA256 "unbound"
#endif
#if !defined(LEO2_C5_CORE_GIT_SHA)
#define LEO2_C5_CORE_GIT_SHA "unbound"
#endif
#if !defined(LEO2_C5_LIBRARY_SHA256)
#define LEO2_C5_LIBRARY_SHA256 "unbound"
#endif
#if !defined(LEO2_C5_SANITIZER_MODE)
#define LEO2_C5_SANITIZER_MODE "unbound"
#endif

namespace {

typedef uint16_t Element;
typedef std::chrono::steady_clock Clock;

static const size_t kKernelQuantum = 64;
static const size_t kGuardBytes = 64;
static const uint8_t kGuardValue = 0xa5;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

size_t round_up_64(size_t value)
{
    require(value != 0, "shard byte count must be positive");
    require(value <= std::numeric_limits<size_t>::max() - 63,
        "shard byte count overflow");
    return (value + 63) & ~static_cast<size_t>(63);
}

unsigned ceil_log2(unsigned value)
{
    unsigned result = 0;
    unsigned power = 1;
    while (power < value)
    {
        require(power <= std::numeric_limits<unsigned>::max() / 2,
            "ceil_log2 overflow");
        power <<= 1;
        ++result;
    }
    return result;
}

unsigned ceil_pow2(unsigned value)
{
    return 1u << ceil_log2(value);
}

unsigned floor_pow2(unsigned value)
{
    require(value != 0, "floor_pow2 requires a positive value");
    unsigned result = 1;
    while (result <= value / 2)
        result <<= 1;
    return result;
}

unsigned integer_log2(unsigned value)
{
    require(value != 0 && (value & (value - 1)) == 0,
        "integer_log2 requires a power of two");
    unsigned result = 0;
    while ((1u << result) != value)
        ++result;
    return result;
}

class BinaryField
{
public:
    BinaryField(unsigned bits, uint32_t polynomial, const Element* basis)
        : bits_(bits)
        , order_(1u << bits)
        , polynomial_(polynomial)
        , coordinate_to_polynomial_(order_, 0)
        , polynomial_to_coordinate_(order_, UINT32_MAX)
        , subspace_normalizer_(bits, 0)
    {
        require(bits == 8 || bits == 16, "C5 supports legacy GF8/GF16 only");
        require((polynomial & (1u << bits)) != 0,
            "field polynomial has the wrong degree");
        for (unsigned coordinate = 0; coordinate < order_; ++coordinate)
        {
            Element polynomial_value = 0;
            for (unsigned bit = 0; bit < bits_; ++bit)
                if ((coordinate & (1u << bit)) != 0)
                    polynomial_value ^= basis[bit];
            require(polynomial_value < order_, "basis value outside field");
            require(polynomial_to_coordinate_[polynomial_value] == UINT32_MAX,
                "coordinate basis is not independent");
            coordinate_to_polynomial_[coordinate] = polynomial_value;
            polynomial_to_coordinate_[polynomial_value] = coordinate;
        }

        // c_j = s_j(omega_(2^j)).  The recurrence follows from the
        // linearized subspace-polynomial identity
        // s_(j+1)(x) = s_j(x) * (s_j(x) + c_j).
        for (unsigned bit = 0; bit < bits_; ++bit)
        {
            Element value = static_cast<Element>(1u << bit);
            for (unsigned level = 0; level < bit; ++level)
                value = multiply(value,
                    static_cast<Element>(value ^ subspace_normalizer_[level]));
            require(value != 0, "subspace normalizer is zero");
            subspace_normalizer_[bit] = value;
        }
    }

    unsigned bits() const { return bits_; }
    unsigned order() const { return order_; }

    Element multiply(Element left, Element right) const
    {
        uint32_t a = coordinate_to_polynomial_[left];
        uint32_t b = coordinate_to_polynomial_[right];
        uint32_t product = 0;
        while (b != 0)
        {
            if ((b & 1u) != 0)
                product ^= a;
            b >>= 1;
            a <<= 1;
        }
        for (int bit = static_cast<int>(2 * bits_ - 2);
             bit >= static_cast<int>(bits_); --bit)
        {
            if ((product & (1u << bit)) != 0)
                product ^= polynomial_ << (bit - bits_);
        }
        require(product < order_ &&
                polynomial_to_coordinate_[product] != UINT32_MAX,
            "field reduction escaped the coordinate map");
        return static_cast<Element>(polynomial_to_coordinate_[product]);
    }

    Element power(Element value, uint32_t exponent) const
    {
        Element result = 1;
        while (exponent != 0)
        {
            if ((exponent & 1u) != 0)
                result = multiply(result, value);
            exponent >>= 1;
            if (exponent != 0)
                value = multiply(value, value);
        }
        return result;
    }

    Element inverse(Element value) const
    {
        require(value != 0, "zero has no inverse");
        return power(value, order_ - 2);
    }

    Element subspace_at(unsigned bit, Element point) const
    {
        require(bit < bits_, "subspace bit outside field");
        Element value = point;
        for (unsigned level = 0; level < bit; ++level)
            value = multiply(value,
                static_cast<Element>(value ^ subspace_normalizer_[level]));
        return value;
    }

    Element lch_normalizer(unsigned index) const
    {
        require(index < order_, "LCH index outside field");
        Element result = 1;
        for (unsigned bit = 0; bit < bits_; ++bit)
            if ((index & (1u << bit)) != 0)
                result = multiply(result, subspace_normalizer_[bit]);
        return result;
    }

    Element lch_basis(unsigned index, Element point) const
    {
        require(index < order_, "LCH index outside field");
        Element numerator = 1;
        for (unsigned bit = 0; bit < bits_; ++bit)
            if ((index & (1u << bit)) != 0)
                numerator = multiply(numerator, subspace_at(bit, point));
        if (numerator == 0)
            return 0;
        return multiply(numerator, inverse(lch_normalizer(index)));
    }

    Element load_symbol(const uint8_t* shard, size_t chunk, unsigned lane) const
    {
        const uint8_t* base = shard + chunk * kKernelQuantum;
        if (bits_ == 8)
        {
            require(lane < 64, "GF8 lane outside chunk");
            return base[lane];
        }
        require(lane < 32, "GF16 lane outside chunk");
        return static_cast<Element>(base[lane] |
            (static_cast<unsigned>(base[32 + lane]) << 8));
    }

    void store_symbol(
        uint8_t* shard, size_t chunk, unsigned lane, Element value) const
    {
        uint8_t* base = shard + chunk * kKernelQuantum;
        if (bits_ == 8)
        {
            require(lane < 64, "GF8 lane outside chunk");
            base[lane] = static_cast<uint8_t>(value);
            return;
        }
        require(lane < 32, "GF16 lane outside chunk");
        base[lane] = static_cast<uint8_t>(value);
        base[32 + lane] = static_cast<uint8_t>(value >> 8);
    }

private:
    unsigned bits_;
    unsigned order_;
    uint32_t polynomial_;
    std::vector<Element> coordinate_to_polynomial_;
    std::vector<uint32_t> polynomial_to_coordinate_;
    std::vector<Element> subspace_normalizer_;
};

BinaryField make_gf8()
{
    static const Element basis[] = {
        1, 214, 152, 146, 86, 200, 88, 230
    };
    return BinaryField(8, 0x11d, basis);
}

BinaryField make_gf16()
{
    static const Element basis[] = {
        0x0001, 0xacca, 0x3c0e, 0x163e,
        0xc582, 0xed2e, 0x914c, 0x4012,
        0x6c98, 0x10d8, 0x6a72, 0xb900,
        0xfdb8, 0xfb34, 0xff38, 0x991e
    };
    return BinaryField(16, 0x1002d, basis);
}

Element production_normalizer(const BinaryField& field, unsigned index)
{
    return field.bits() == 8
        ? static_cast<Element>(leopard::ff8::TestOnlyLchNormalizer(index))
        : static_cast<Element>(leopard::ff16::TestOnlyLchNormalizer(index));
}

Element production_subspace_at(
    const BinaryField& field, unsigned size, unsigned shift)
{
    return field.bits() == 8
        ? static_cast<Element>(leopard::ff8::TestOnlySubspaceAt(size, shift))
        : static_cast<Element>(leopard::ff16::TestOnlySubspaceAt(size, shift));
}

void production_forward(
    const BinaryField& field,
    size_t processed_bytes,
    unsigned size,
    unsigned shift,
    std::vector<void*>* pointers)
{
    require(pointers->size() >= size, "pointer workspace is too small");
    if (field.bits() == 8)
        leopard::ff8::TestOnlyLchForward(
            processed_bytes, size, shift, size, &(*pointers)[0]);
    else
        leopard::ff16::TestOnlyLchForward(
            processed_bytes, size, shift, size, &(*pointers)[0]);
}

class ShardSet
{
public:
    ShardSet() : count_(0), valid_bytes_(0), processed_bytes_(0), pitch_(0) {}

    ShardSet(unsigned count, size_t valid_bytes)
        : count_(count)
        , valid_bytes_(valid_bytes)
        , processed_bytes_(round_up_64(valid_bytes))
        , pitch_(processed_bytes_ + kGuardBytes)
        , memory_(static_cast<size_t>(count) * pitch_, kGuardValue)
    {
        require(count != 0, "ShardSet requires at least one shard");
        for (unsigned i = 0; i < count_; ++i)
            std::memset(ptr(i), 0, processed_bytes_);
    }

    unsigned count() const { return count_; }
    size_t valid_bytes() const { return valid_bytes_; }
    size_t processed_bytes() const { return processed_bytes_; }
    size_t resident_bytes() const { return memory_.size(); }

    uint8_t* ptr(unsigned index)
    {
        require(index < count_, "ShardSet index out of range");
        return &memory_[static_cast<size_t>(index) * pitch_];
    }

    const uint8_t* ptr(unsigned index) const
    {
        require(index < count_, "ShardSet index out of range");
        return &memory_[static_cast<size_t>(index) * pitch_];
    }

    bool guards_intact() const
    {
        for (unsigned shard = 0; shard < count_; ++shard)
        {
            const uint8_t* guard = ptr(shard) + processed_bytes_;
            for (size_t i = 0; i < kGuardBytes; ++i)
                if (guard[i] != kGuardValue)
                    return false;
        }
        return true;
    }

private:
    unsigned count_;
    size_t valid_bytes_;
    size_t processed_bytes_;
    size_t pitch_;
    std::vector<uint8_t> memory_;
};

uint64_t splitmix64(uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

void fill_inputs(ShardSet* shards, uint64_t seed)
{
    for (unsigned shard = 0; shard < shards->count(); ++shard)
    {
        uint8_t* data = shards->ptr(shard);
        for (size_t byte = 0; byte < shards->valid_bytes(); ++byte)
        {
            const uint64_t value = seed ^
                (static_cast<uint64_t>(shard + 1) * UINT64_C(0xd6e8feb86659fd93)) ^
                (static_cast<uint64_t>(byte + 1) * UINT64_C(0xa0761d6478bd642f));
            data[byte] = static_cast<uint8_t>(splitmix64(value));
        }
        // The old kernels operate in 64-byte quanta.  Zero padding defines
        // the exact tail experiment without reading outside allocated memory.
        if (shards->processed_bytes() > shards->valid_bytes())
            std::memset(data + shards->valid_bytes(), 0,
                shards->processed_bytes() - shards->valid_bytes());
    }
}

struct Block
{
    unsigned offset;
    unsigned size;
};

std::vector<Block> binary_prefix_blocks(unsigned q)
{
    require(q != 0, "binary prefix must be nonempty");
    std::vector<Block> result;
    unsigned offset = 0;
    for (unsigned size = floor_pow2(q); size != 0; size >>= 1)
    {
        if ((q & size) == 0)
            continue;
        require((offset & (size - 1)) == 0,
            "canonical binary block is not aligned");
        Block block = { offset, size };
        result.push_back(block);
        offset += size;
    }
    require(offset == q, "binary blocks do not cover prefix");
    return result;
}

struct Job
{
    unsigned offset;
    unsigned size;
    unsigned shift;
    Element factor;
    bool base;
};

struct Traffic
{
    uint64_t input_read_bytes;
    uint64_t output_read_bytes;
    uint64_t output_write_bytes;
    uint64_t scratch_read_bytes;
    uint64_t scratch_write_bytes;
    uint64_t butterfly_equivalents;
    uint64_t factor_symbols;
    uint64_t accumulation_symbols;

    Traffic()
        : input_read_bytes(0)
        , output_read_bytes(0)
        , output_write_bytes(0)
        , scratch_read_bytes(0)
        , scratch_write_bytes(0)
        , butterfly_equivalents(0)
        , factor_symbols(0)
        , accumulation_symbols(0)
    {}

    uint64_t total_bytes() const
    {
        return input_read_bytes + output_read_bytes + output_write_bytes +
            scratch_read_bytes + scratch_write_bytes;
    }
};

class Plan
{
public:
    Plan(const BinaryField& field, unsigned q, size_t processed_bytes)
        : field_(&field)
        , q_(q)
        , parent_(ceil_pow2(q))
        , processed_bytes_(processed_bytes)
        , largest_block_(floor_pow2(q))
        , largest_tail_block_(0)
        , fused_power_plus_one_(q == floor_pow2(q) + 1)
        , factor_checks_(0)
    {
        require(q <= field.order(), "prefix exceeds field");
        require(processed_bytes != 0 && processed_bytes % 64 == 0,
            "processed bytes must be a nonzero multiple of 64");
        const std::vector<Block> blocks = binary_prefix_blocks(q);
        blocks_ = blocks;
        for (size_t block_i = 0; block_i < blocks.size(); ++block_i)
        {
            const Block& block = blocks[block_i];
            if (block_i != 0)
                largest_tail_block_ = std::max(largest_tail_block_, block.size);
            require(production_normalizer(field, block.offset) ==
                    field.lch_normalizer(block.offset),
                "independent LCH normalizer differs from production");
            ++factor_checks_;
            for (unsigned shift = 0; shift < parent_; shift += block.size)
            {
                const Element factor = field.lch_basis(
                    block.offset, static_cast<Element>(shift));
                for (unsigned bit = 0; bit < field.bits(); ++bit)
                {
                    if ((block.offset & (1u << bit)) == 0)
                        continue;
                    const unsigned size = 1u << bit;
                    const Element expected = field.subspace_at(
                        bit, static_cast<Element>(shift));
                    const Element actual = production_subspace_at(
                        field, size, shift);
                    require(actual == expected,
                        "independent subspace factor differs from production");
                    ++factor_checks_;
                }
                if (factor == 0)
                    continue;
                const Job job = {
                    block.offset, block.size, shift, factor, block_i == 0
                };
                jobs_.push_back(job);
            }
        }
        require(!jobs_.empty(), "plan has no active jobs");
        for (unsigned shift = 0; shift < parent_; shift += largest_block_)
        {
            bool found = false;
            for (size_t i = 0; i < jobs_.size(); ++i)
                found = found || (jobs_[i].base && jobs_[i].shift == shift);
            require(found, "base block does not initialize every output");
        }
        candidate_traffic_ = model_candidate();
        padded_traffic_ = model_padded();
    }

    const BinaryField& field() const { return *field_; }
    unsigned q() const { return q_; }
    unsigned parent() const { return parent_; }
    size_t processed_bytes() const { return processed_bytes_; }
    unsigned largest_block() const { return largest_block_; }
    unsigned largest_tail_block() const { return largest_tail_block_; }
    bool fused_power_plus_one() const { return fused_power_plus_one_; }
    const std::vector<Block>& blocks() const { return blocks_; }
    const std::vector<Job>& jobs() const { return jobs_; }
    uint64_t factor_checks() const { return factor_checks_; }
    const Traffic& candidate_traffic() const { return candidate_traffic_; }
    const Traffic& padded_traffic() const { return padded_traffic_; }

    size_t resident_plan_bytes() const
    {
        return sizeof(*this) + blocks_.capacity() * sizeof(Block) +
            jobs_.capacity() * sizeof(Job);
    }

private:
    Traffic model_candidate() const
    {
        Traffic result;
        const uint64_t symbols_per_shard = processed_bytes_ /
            (field_->bits() == 8 ? 1 : 2);
        for (size_t i = 0; i < jobs_.size(); ++i)
        {
            const Job& job = jobs_[i];
            if (fused_power_plus_one_ && !job.base)
                continue;
            const uint64_t shard_bytes =
                static_cast<uint64_t>(job.size) * processed_bytes_;
            result.input_read_bytes += shard_bytes;
            result.butterfly_equivalents += static_cast<uint64_t>(job.size / 2) *
                integer_log2(job.size);
            const uint64_t butterfly_bytes =
                static_cast<uint64_t>(4) * processed_bytes_ *
                (static_cast<uint64_t>(job.size / 2) * integer_log2(job.size));
            if (job.base)
            {
                result.output_write_bytes += shard_bytes;
                result.output_read_bytes += butterfly_bytes / 2;
                result.output_write_bytes += butterfly_bytes / 2;
                if (fused_power_plus_one_ &&
                    job.shift == largest_block_)
                {
                    // The one-element tail is X_largest.  It is zero on the
                    // lower coset and one on the upper coset, so inject it
                    // into coefficient zero before the shifted transform.
                    result.input_read_bytes += processed_bytes_;
                    result.output_read_bytes += processed_bytes_;
                    result.output_write_bytes += processed_bytes_;
                    result.accumulation_symbols += symbols_per_shard;
                }
            }
            else
            {
                if (job.size > 1)
                {
                    result.scratch_write_bytes += shard_bytes;
                    result.scratch_read_bytes += butterfly_bytes / 2;
                    result.scratch_write_bytes += butterfly_bytes / 2;
                    result.scratch_read_bytes += shard_bytes;
                }
                result.output_read_bytes += shard_bytes;
                result.output_write_bytes += shard_bytes;
                result.accumulation_symbols +=
                    static_cast<uint64_t>(job.size) * symbols_per_shard;
                if (job.factor != 1)
                    result.factor_symbols +=
                        static_cast<uint64_t>(job.size) * symbols_per_shard;
            }
        }
        return result;
    }

    Traffic model_padded() const
    {
        Traffic result;
        result.input_read_bytes = static_cast<uint64_t>(q_) * processed_bytes_;
        result.output_write_bytes = static_cast<uint64_t>(parent_) * processed_bytes_;
        result.butterfly_equivalents = static_cast<uint64_t>(parent_ / 2) *
            integer_log2(parent_);
        const uint64_t butterfly_bytes = static_cast<uint64_t>(4) *
            processed_bytes_ * result.butterfly_equivalents;
        result.output_read_bytes += butterfly_bytes / 2;
        result.output_write_bytes += butterfly_bytes / 2;
        return result;
    }

    const BinaryField* field_;
    unsigned q_;
    unsigned parent_;
    size_t processed_bytes_;
    unsigned largest_block_;
    unsigned largest_tail_block_;
    bool fused_power_plus_one_;
    uint64_t factor_checks_;
    std::vector<Block> blocks_;
    std::vector<Job> jobs_;
    Traffic candidate_traffic_;
    Traffic padded_traffic_;
};

void xor_bytes(uint8_t* destination, const uint8_t* source, size_t bytes)
{
    for (size_t i = 0; i < bytes; ++i)
        destination[i] ^= source[i];
}

void multiply_xor(
    const BinaryField& field,
    Element factor,
    uint8_t* destination,
    const uint8_t* source,
    size_t bytes)
{
    require(bytes % kKernelQuantum == 0,
        "multiply_xor requires complete kernel chunks");
    if (factor == 0)
        return;
    if (factor == 1)
    {
        xor_bytes(destination, source, bytes);
        return;
    }
    const size_t chunks = bytes / kKernelQuantum;
    const unsigned lanes = field.bits() == 8 ? 64 : 32;
    for (size_t chunk = 0; chunk < chunks; ++chunk)
    {
        for (unsigned lane = 0; lane < lanes; ++lane)
        {
            const Element product = field.multiply(
                factor, field.load_symbol(source, chunk, lane));
            const Element value = static_cast<Element>(
                field.load_symbol(destination, chunk, lane) ^ product);
            field.store_symbol(destination, chunk, lane, value);
        }
    }
}

class CandidateWorkspace
{
public:
    explicit CandidateWorkspace(const Plan& plan)
        : pointers_(plan.largest_block(), NULL)
    {
        if (!plan.fused_power_plus_one() && plan.largest_tail_block() > 1)
            scratch_.reset(new ShardSet(
                plan.largest_tail_block(), plan.processed_bytes()));
    }

    size_t resident_bytes() const
    {
        return (scratch_.get() ? scratch_->resident_bytes() : 0) +
            pointers_.capacity() * sizeof(void*);
    }

    void execute(const Plan& plan, const ShardSet& input, ShardSet* output)
    {
        require(input.count() == plan.q(), "candidate input count mismatch");
        require(output->count() == plan.parent(),
            "candidate output count mismatch");
        require(input.processed_bytes() == plan.processed_bytes() &&
                output->processed_bytes() == plan.processed_bytes(),
            "candidate byte count mismatch");
        for (size_t job_i = 0; job_i < plan.jobs().size(); ++job_i)
        {
            const Job& job = plan.jobs()[job_i];
            if (plan.fused_power_plus_one() && !job.base)
                continue;
            if (job.base)
            {
                require(job.factor == 1, "base job factor is not one");
                for (unsigned i = 0; i < job.size; ++i)
                {
                    std::memcpy(output->ptr(job.shift + i), input.ptr(i),
                        plan.processed_bytes());
                    pointers_[i] = output->ptr(job.shift + i);
                }
                if (plan.fused_power_plus_one() &&
                    job.shift == plan.largest_block())
                {
                    xor_bytes(output->ptr(job.shift),
                        input.ptr(plan.largest_block()),
                        plan.processed_bytes());
                }
                if (job.size > 1)
                    production_forward(plan.field(), plan.processed_bytes(),
                        job.size, job.shift, &pointers_);
                continue;
            }

            if (job.size == 1)
            {
                multiply_xor(plan.field(), job.factor,
                    output->ptr(job.shift), input.ptr(job.offset),
                    plan.processed_bytes());
                continue;
            }
            require(scratch_.get() && job.size <= scratch_->count(),
                "candidate scratch is too small");
            for (unsigned i = 0; i < job.size; ++i)
            {
                std::memcpy(scratch_->ptr(i), input.ptr(job.offset + i),
                    plan.processed_bytes());
                pointers_[i] = scratch_->ptr(i);
            }
            if (job.size > 1)
                production_forward(plan.field(), plan.processed_bytes(),
                    job.size, job.shift, &pointers_);
            for (unsigned i = 0; i < job.size; ++i)
                multiply_xor(plan.field(), job.factor,
                    output->ptr(job.shift + i), scratch_->ptr(i),
                    plan.processed_bytes());
        }
        require(!scratch_.get() || scratch_->guards_intact(),
            "candidate scratch guard changed");
        require(input.guards_intact(), "candidate input guard changed");
        require(output->guards_intact(), "candidate output guard changed");
    }

private:
    std::unique_ptr<ShardSet> scratch_;
    std::vector<void*> pointers_;
};

class PaddedWorkspace
{
public:
    explicit PaddedWorkspace(const Plan& plan)
        : pointers_(plan.parent(), NULL)
    {}

    size_t resident_bytes() const
    {
        return pointers_.capacity() * sizeof(void*);
    }

    void execute(const Plan& plan, const ShardSet& input, ShardSet* output)
    {
        require(input.count() == plan.q(), "padded input count mismatch");
        require(output->count() == plan.parent(), "padded output count mismatch");
        for (unsigned i = 0; i < plan.q(); ++i)
            std::memcpy(output->ptr(i), input.ptr(i), plan.processed_bytes());
        for (unsigned i = plan.q(); i < plan.parent(); ++i)
            std::memset(output->ptr(i), 0, plan.processed_bytes());
        for (unsigned i = 0; i < plan.parent(); ++i)
            pointers_[i] = output->ptr(i);
        if (plan.parent() > 1)
            production_forward(plan.field(), plan.processed_bytes(),
                plan.parent(), 0, &pointers_);
        require(input.guards_intact(), "padded input guard changed");
        require(output->guards_intact(), "padded output guard changed");
    }

private:
    std::vector<void*> pointers_;
};

bool equal_shards(const ShardSet& left, const ShardSet& right)
{
    require(left.count() == right.count() &&
            left.processed_bytes() == right.processed_bytes(),
        "cannot compare unlike shard sets");
    for (unsigned i = 0; i < left.count(); ++i)
        if (std::memcmp(left.ptr(i), right.ptr(i),
                left.processed_bytes()) != 0)
            return false;
    return true;
}

uint64_t fnv_update(uint64_t digest, const void* data, size_t bytes)
{
    const uint8_t* source = static_cast<const uint8_t*>(data);
    for (size_t i = 0; i < bytes; ++i)
    {
        digest ^= source[i];
        digest *= UINT64_C(1099511628211);
    }
    return digest;
}

uint64_t digest_shards(uint64_t digest, const ShardSet& shards)
{
    for (unsigned i = 0; i < shards.count(); ++i)
        digest = fnv_update(digest, shards.ptr(i), shards.valid_bytes());
    return digest;
}

uint64_t verify_direct_samples(
    const Plan& plan,
    const ShardSet& input,
    const ShardSet& output,
    uint64_t* digest)
{
    std::vector<unsigned> points;
    points.push_back(0);
    if (plan.parent() > 1)
        points.push_back(1);
    if (plan.parent() > 2)
        points.push_back(plan.parent() / 2);
    if (plan.parent() > 3)
        points.push_back(plan.parent() - 1);
    const unsigned lanes = plan.field().bits() == 8 ? 64 : 32;
    const unsigned sampled_lanes = std::min(4u, lanes);
    uint64_t checks = 0;
    for (size_t point_i = 0; point_i < points.size(); ++point_i)
    {
        const unsigned point = points[point_i];
        for (unsigned lane = 0; lane < sampled_lanes; ++lane)
        {
            Element expected = 0;
            for (unsigned coefficient = 0; coefficient < plan.q(); ++coefficient)
            {
                const Element value = plan.field().load_symbol(
                    input.ptr(coefficient), 0, lane);
                const Element basis = plan.field().lch_basis(
                    coefficient, static_cast<Element>(point));
                if (value != 0 && basis != 0)
                    expected ^= plan.field().multiply(value, basis);
            }
            const Element actual = plan.field().load_symbol(
                output.ptr(point), 0, lane);
            require(actual == expected,
                "direct normalized-LCH sample differs from production output");
            *digest = fnv_update(*digest, &actual, sizeof(actual));
            ++checks;
        }
    }
    return checks;
}

struct CorrectnessCase
{
    unsigned field_bits;
    unsigned q;
    size_t valid_bytes;
};

struct CorrectnessResult
{
    unsigned field_bits;
    unsigned q;
    unsigned parent;
    size_t valid_bytes;
    size_t processed_bytes;
    unsigned blocks;
    unsigned jobs;
    uint64_t factor_checks;
    uint64_t compared_bytes;
    uint64_t direct_symbol_checks;
    size_t candidate_scratch_bytes;
    size_t padded_scratch_bytes;
    uint64_t candidate_traffic_bytes;
    uint64_t padded_traffic_bytes;
};

CorrectnessResult verify_case(
    const BinaryField& field,
    const CorrectnessCase& test_case,
    bool run_direct,
    uint64_t* digest)
{
    const size_t processed_bytes = round_up_64(test_case.valid_bytes);
    const Plan plan(field, test_case.q, processed_bytes);
    ShardSet input(plan.q(), test_case.valid_bytes);
    ShardSet candidate(plan.parent(), test_case.valid_bytes);
    ShardSet padded(plan.parent(), test_case.valid_bytes);
    fill_inputs(&input,
        UINT64_C(0xc5d1ad1c5eed0000) ^
        (static_cast<uint64_t>(test_case.field_bits) << 48) ^
        (static_cast<uint64_t>(test_case.q) << 16) ^ test_case.valid_bytes);
    CandidateWorkspace candidate_workspace(plan);
    PaddedWorkspace padded_workspace(plan);
    candidate_workspace.execute(plan, input, &candidate);
    padded_workspace.execute(plan, input, &padded);
    require(equal_shards(candidate, padded),
        "dyadic candidate differs from padded parent");
    *digest = digest_shards(*digest, candidate);
    const uint64_t direct_checks = run_direct
        ? verify_direct_samples(plan, input, padded, digest)
        : 0;

    CorrectnessResult result;
    result.field_bits = field.bits();
    result.q = plan.q();
    result.parent = plan.parent();
    result.valid_bytes = test_case.valid_bytes;
    result.processed_bytes = processed_bytes;
    result.blocks = static_cast<unsigned>(plan.blocks().size());
    result.jobs = static_cast<unsigned>(plan.jobs().size());
    result.factor_checks = plan.factor_checks();
    result.compared_bytes = static_cast<uint64_t>(plan.parent()) * processed_bytes;
    result.direct_symbol_checks = direct_checks;
    result.candidate_scratch_bytes = candidate_workspace.resident_bytes();
    result.padded_scratch_bytes = padded_workspace.resident_bytes();
    result.candidate_traffic_bytes = plan.candidate_traffic().total_bytes();
    result.padded_traffic_bytes = plan.padded_traffic().total_bytes();
    return result;
}

double median(std::vector<double> values)
{
    require(!values.empty(), "median requires samples");
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2;
    return values.size() % 2 != 0
        ? values[middle]
        : (values[middle - 1] + values[middle]) / 2.;
}

double mad(const std::vector<double>& values, double center)
{
    std::vector<double> deviations;
    deviations.reserve(values.size());
    for (size_t i = 0; i < values.size(); ++i)
        deviations.push_back(std::fabs(values[i] - center));
    return median(deviations);
}

double elapsed_us(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration_cast<
        std::chrono::duration<double, std::micro> >(end - begin).count();
}

struct BenchmarkSpec
{
    unsigned field_bits;
    unsigned q;
    size_t valid_bytes;
    unsigned batch;
    unsigned reuse;
};

struct Stripe
{
    ShardSet input;
    ShardSet output;
    CandidateWorkspace candidate;
    PaddedWorkspace padded;

    Stripe(const Plan& plan, size_t valid_bytes, uint64_t seed)
        : input(plan.q(), valid_bytes)
        , output(plan.parent(), valid_bytes)
        , candidate(plan)
        , padded(plan)
    {
        fill_inputs(&input, seed);
    }
};

struct BenchmarkResult
{
    unsigned field_bits;
    unsigned q;
    unsigned parent;
    size_t valid_bytes;
    size_t processed_bytes;
    unsigned batch;
    unsigned reuse;
    unsigned samples;
    double setup_median_us;
    double setup_mad_us;
    double candidate_median_us;
    double candidate_mad_us;
    double padded_median_us;
    double padded_mad_us;
    double padded_over_candidate;
    double credible_gain_percent;
    uint64_t candidate_scratch_bytes_per_stripe;
    uint64_t padded_scratch_bytes_per_stripe;
    uint64_t resident_batch_bytes;
    uint64_t candidate_traffic_bytes_per_execution;
    uint64_t padded_traffic_bytes_per_execution;
};

BenchmarkResult benchmark_case(
    const BinaryField& field,
    const BenchmarkSpec& spec,
    unsigned scale,
    volatile uint64_t* sink)
{
    const size_t processed_bytes = round_up_64(spec.valid_bytes);
    std::vector<double> setup_samples;
    const unsigned setup_repetitions = 31 * std::max(1u, scale);
    for (unsigned i = 0; i < setup_repetitions; ++i)
    {
        const Clock::time_point begin = Clock::now();
        const Plan temporary(field, spec.q, processed_bytes);
        const Clock::time_point end = Clock::now();
        setup_samples.push_back(elapsed_us(begin, end));
        *sink ^= temporary.jobs().size();
    }
    const Plan plan(field, spec.q, processed_bytes);

    std::vector<std::unique_ptr<Stripe> > stripes;
    stripes.reserve(spec.batch);
    uint64_t resident_bytes = 0;
    for (unsigned i = 0; i < spec.batch; ++i)
    {
        stripes.push_back(std::unique_ptr<Stripe>(new Stripe(
            plan, spec.valid_bytes,
            UINT64_C(0xc5ba7c4000000000) ^
                (static_cast<uint64_t>(spec.q) << 32) ^
                (static_cast<uint64_t>(spec.valid_bytes) << 8) ^ i)));
        resident_bytes += stripes.back()->input.resident_bytes();
        resident_bytes += stripes.back()->output.resident_bytes();
        resident_bytes += stripes.back()->candidate.resident_bytes();
        resident_bytes += stripes.back()->padded.resident_bytes();
    }

    for (unsigned warmup = 0; warmup < 3; ++warmup)
    {
        for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
            stripes[stripe]->candidate.execute(
                plan, stripes[stripe]->input, &stripes[stripe]->output);
        for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
            stripes[stripe]->padded.execute(
                plan, stripes[stripe]->input, &stripes[stripe]->output);
    }

    unsigned sample_count = spec.valid_bytes >= 65536 ? 7 : 11;
    sample_count *= std::max(1u, scale);
    if ((sample_count & 1u) == 0)
        ++sample_count;
    std::vector<double> candidate_samples;
    std::vector<double> padded_samples;
    const double executions = static_cast<double>(spec.batch) * spec.reuse;
    for (unsigned sample = 0; sample < sample_count; ++sample)
    {
        const bool candidate_first = (sample & 1u) == 0;
        for (unsigned order = 0; order < 2; ++order)
        {
            const bool run_candidate = (order == 0) == candidate_first;
            const Clock::time_point begin = Clock::now();
            for (unsigned reuse = 0; reuse < spec.reuse; ++reuse)
                for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
                {
                    Stripe& item = *stripes[stripe];
                    if (run_candidate)
                        item.candidate.execute(plan, item.input, &item.output);
                    else
                        item.padded.execute(plan, item.input, &item.output);
                }
            const Clock::time_point end = Clock::now();
            const double per_execution = elapsed_us(begin, end) / executions;
            if (run_candidate)
                candidate_samples.push_back(per_execution);
            else
                padded_samples.push_back(per_execution);
            Stripe& selected = *stripes[sample % stripes.size()];
            *sink ^= selected.output.ptr(sample % plan.parent())[
                sample % selected.output.valid_bytes()];
        }
    }

    // Untimed final identity gate for every distinct batch member.
    for (unsigned stripe = 0; stripe < spec.batch; ++stripe)
    {
        Stripe& item = *stripes[stripe];
        ShardSet expected(plan.parent(), spec.valid_bytes);
        item.candidate.execute(plan, item.input, &item.output);
        item.padded.execute(plan, item.input, &expected);
        require(equal_shards(item.output, expected),
            "benchmark candidate differs from padded parent");
    }

    BenchmarkResult result;
    result.field_bits = field.bits();
    result.q = spec.q;
    result.parent = plan.parent();
    result.valid_bytes = spec.valid_bytes;
    result.processed_bytes = processed_bytes;
    result.batch = spec.batch;
    result.reuse = spec.reuse;
    result.samples = sample_count;
    result.setup_median_us = median(setup_samples);
    result.setup_mad_us = mad(setup_samples, result.setup_median_us);
    result.candidate_median_us = median(candidate_samples);
    result.candidate_mad_us = mad(candidate_samples, result.candidate_median_us);
    result.padded_median_us = median(padded_samples);
    result.padded_mad_us = mad(padded_samples, result.padded_median_us);
    result.padded_over_candidate =
        result.padded_median_us / result.candidate_median_us;
    result.credible_gain_percent = 100. * (
        (result.padded_median_us - result.padded_mad_us) /
        (result.candidate_median_us + result.candidate_mad_us) - 1.);
    result.candidate_scratch_bytes_per_stripe =
        stripes[0]->candidate.resident_bytes();
    result.padded_scratch_bytes_per_stripe =
        stripes[0]->padded.resident_bytes();
    result.resident_batch_bytes = resident_bytes;
    result.candidate_traffic_bytes_per_execution =
        plan.candidate_traffic().total_bytes();
    result.padded_traffic_bytes_per_execution =
        plan.padded_traffic().total_bytes();
    return result;
}

const char* backend_name(leo2_backend backend)
{
    switch (backend)
    {
    case LEO2_BACKEND_SCALAR: return "scalar";
    case LEO2_BACKEND_SSSE3: return "ssse3";
    case LEO2_BACKEND_AVX2: return "avx2";
    case LEO2_BACKEND_NEON: return "neon";
    default: return "auto";
    }
}

std::string json_escape(const std::string& value)
{
    std::ostringstream output;
    for (size_t i = 0; i < value.size(); ++i)
    {
        const unsigned char ch = static_cast<unsigned char>(value[i]);
        if (ch == '\\' || ch == '"')
            output << '\\' << ch;
        else if (ch >= 0x20)
            output << ch;
    }
    return output.str();
}

struct RuntimeEnvironment
{
    std::vector<unsigned> affinity_cpus;
    std::string omp_num_threads;
    unsigned openmp_max_threads;
    std::string hostname;
    std::string machine;
    std::string cpu_model;
};

std::string trim(const std::string& value)
{
    const size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    const size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

RuntimeEnvironment read_runtime_environment()
{
    RuntimeEnvironment result;
    const char* omp_threads = std::getenv("OMP_NUM_THREADS");
    result.omp_num_threads = omp_threads ? omp_threads : "unset";
#if defined(_OPENMP)
    result.openmp_max_threads = static_cast<unsigned>(omp_get_max_threads());
#else
    result.openmp_max_threads = 0;
#endif
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) == 0)
        for (unsigned cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &allowed))
                result.affinity_cpus.push_back(cpu);
    std::array<char, 256> hostname = {{ 0 }};
    result.hostname = gethostname(hostname.data(), hostname.size()) == 0
        ? hostname.data() : "unavailable";
    struct utsname name;
    result.machine = uname(&name) == 0
        ? std::string(name.sysname) + " " + name.release + " " + name.machine
        : "unavailable";
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line))
    {
        const size_t colon = line.find(':');
        if (colon != std::string::npos &&
            trim(line.substr(0, colon)) == "model name")
        {
            result.cpu_model = trim(line.substr(colon + 1));
            break;
        }
    }
#else
    result.hostname = result.machine = "unavailable";
#endif
    if (result.cpu_model.empty())
        result.cpu_model = "unavailable";
    return result;
}

struct Arguments
{
    std::string output_path;
    std::string backend_label;
    std::string mode;
    unsigned benchmark_scale;

    Arguments()
        : backend_label("unspecified"), mode("all"), benchmark_scale(1)
    {}
};

Arguments parse_arguments(int argc, char** argv)
{
    Arguments result;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument(argv[i]);
        if ((argument == "--output" || argument == "--backend-label" ||
             argument == "--mode" || argument == "--benchmark-scale") &&
            i + 1 >= argc)
            throw std::invalid_argument(argument + " requires a value");
        if (argument == "--output")
            result.output_path = argv[++i];
        else if (argument == "--backend-label")
            result.backend_label = argv[++i];
        else if (argument == "--mode")
            result.mode = argv[++i];
        else if (argument == "--benchmark-scale")
            result.benchmark_scale = static_cast<unsigned>(
                std::strtoul(argv[++i], NULL, 10));
        else
            throw std::invalid_argument("unknown argument: " + argument);
    }
    require(result.mode == "all" || result.mode == "correctness" ||
            result.mode == "backend",
        "mode must be all, correctness, or backend");
    require(result.benchmark_scale >= 1, "benchmark scale must be positive");
    return result;
}

std::vector<CorrectnessCase> make_correctness_cases()
{
    static const unsigned gf8_q[] = {
        1, 3, 5, 7, 9, 15, 17, 31, 33, 63, 65, 127, 129, 191, 255
    };
    static const unsigned gf16_q[] = {
        257, 259, 511, 513, 1000, 1025, 4097, 8191
    };
    static const size_t tails[] = { 1, 17, 63, 64, 65, 129 };
    std::vector<CorrectnessCase> result;
    for (size_t q_i = 0; q_i < sizeof(gf8_q) / sizeof(gf8_q[0]); ++q_i)
        for (size_t b_i = 0; b_i < sizeof(tails) / sizeof(tails[0]); ++b_i)
        {
            CorrectnessCase item = { 8, gf8_q[q_i], tails[b_i] };
            result.push_back(item);
        }
    for (size_t q_i = 0; q_i < sizeof(gf16_q) / sizeof(gf16_q[0]); ++q_i)
        for (size_t b_i = 0; b_i < sizeof(tails) / sizeof(tails[0]); ++b_i)
        {
            CorrectnessCase item = { 16, gf16_q[q_i], tails[b_i] };
            result.push_back(item);
        }
    static const size_t boundary_bytes[] = { 4095, 4096, 4097 };
    for (size_t i = 0;
         i < sizeof(boundary_bytes) / sizeof(boundary_bytes[0]); ++i)
    {
        CorrectnessCase gf8 = { 8, 129, boundary_bytes[i] };
        CorrectnessCase gf16 = { 16, 257, boundary_bytes[i] };
        result.push_back(gf8);
        result.push_back(gf16);
    }
    return result;
}

std::vector<BenchmarkSpec> make_benchmark_specs(bool backend_only)
{
    std::vector<BenchmarkSpec> result;
    const unsigned fields[] = { 8, 8, 8, 16, 16 };
    const unsigned qs[] = { 33, 65, 129, 257, 513 };
    const unsigned geometry_count =
        static_cast<unsigned>(sizeof(qs) / sizeof(qs[0]));
    if (backend_only)
    {
        for (unsigned i = 0; i < geometry_count; ++i)
        {
            BenchmarkSpec item = { fields[i], qs[i], 1024, 1, 8 };
            result.push_back(item);
        }
        return result;
    }
    for (unsigned geometry = 0; geometry < geometry_count; ++geometry)
    {
        const size_t small_bytes[] = { 64, 1024 };
        for (size_t b = 0; b < 2; ++b)
            for (unsigned batch_i = 0; batch_i < 2; ++batch_i)
                for (unsigned reuse_i = 0; reuse_i < 2; ++reuse_i)
                {
                    BenchmarkSpec item = {
                        fields[geometry], qs[geometry], small_bytes[b],
                        batch_i == 0 ? 1u : 8u,
                        reuse_i == 0 ? 1u : 8u
                    };
                    result.push_back(item);
                }
        for (unsigned reuse_i = 0; reuse_i < 2; ++reuse_i)
        {
            BenchmarkSpec item = {
                fields[geometry], qs[geometry], 65536, 1,
                reuse_i == 0 ? 1u : 8u
            };
            result.push_back(item);
        }
    }
    const unsigned large_q[] = { 33, 65, 129, 257 };
    for (unsigned i = 0; i < 4; ++i)
    {
        BenchmarkSpec item = {
            i < 3 ? 8u : 16u, large_q[i], 1024 * 1024, 1, 1
        };
        result.push_back(item);
    }
    return result;
}

void write_results(
    std::ostream& output,
    const Arguments& arguments,
    leo2_backend runtime_backend,
    const RuntimeEnvironment& environment,
    uint64_t digest,
    uint64_t direct_checks,
    uint64_t factor_checks,
    const std::vector<CorrectnessResult>& correctness,
    const std::vector<BenchmarkResult>& benchmarks)
{
    output << std::fixed << std::setprecision(6);
    output << "{\n"
           << "  \"schema_version\": \"leopard2-c5-cpp-v1\",\n"
           << "  \"status\": \"pass\",\n"
           << "  \"wire_identity\": \"existing padded dyadic parent\",\n"
           << "  \"exact_profile_implemented\": false,\n"
           << "  \"default_build_integration\": false,\n"
           << "  \"mode\": \"" << json_escape(arguments.mode) << "\",\n"
           << "  \"requested_backend\": \""
           << json_escape(arguments.backend_label) << "\",\n"
           << "  \"runtime_backend\": \"" << backend_name(runtime_backend)
           << "\",\n"
           << "  \"kernel_quantum_bytes\": " << kKernelQuantum << ",\n"
           << "  \"tail_policy\": \"zero-pad to kernel quantum; compare valid and padded bytes\",\n"
           << "  \"correctness_digest_fnv1a64\": \"0x"
           << std::hex << std::setw(16) << std::setfill('0') << digest
           << std::dec << std::setfill(' ') << "\",\n"
           << "  \"direct_symbol_checks\": " << direct_checks << ",\n"
           << "  \"factor_checks\": " << factor_checks << ",\n"
           << "  \"runtime_environment\": {\n"
           << "    \"process_affinity_cpus\": [";
    for (size_t i = 0; i < environment.affinity_cpus.size(); ++i)
    {
        if (i != 0)
            output << ", ";
        output << environment.affinity_cpus[i];
    }
    output << "],\n"
           << "    \"omp_num_threads_env\": \""
           << json_escape(environment.omp_num_threads) << "\",\n"
           << "    \"openmp_max_threads\": "
           << environment.openmp_max_threads << ",\n"
           << "    \"hostname\": \"" << json_escape(environment.hostname)
           << "\",\n"
           << "    \"machine\": \"" << json_escape(environment.machine)
           << "\",\n"
           << "    \"cpu_model\": \"" << json_escape(environment.cpu_model)
           << "\"\n"
           << "  },\n"
           << "  \"build_binding\": {\n"
           << "    \"source_sha256\": \"" << LEO2_C5_SOURCE_SHA256 << "\",\n"
           << "    \"core_git_sha\": \"" << LEO2_C5_CORE_GIT_SHA << "\",\n"
           << "    \"linked_library_sha256\": \""
           << LEO2_C5_LIBRARY_SHA256 << "\",\n"
           << "    \"sanitizer_mode\": \"" << LEO2_C5_SANITIZER_MODE << "\",\n"
           << "    \"compiler\": \"" << json_escape(__VERSION__) << "\"\n"
           << "  },\n"
           << "  \"correctness_cases\": [\n";
    for (size_t i = 0; i < correctness.size(); ++i)
    {
        const CorrectnessResult& item = correctness[i];
        output << "    {\"field_bits\": " << item.field_bits
               << ", \"q\": " << item.q
               << ", \"parent\": " << item.parent
               << ", \"valid_bytes\": " << item.valid_bytes
               << ", \"processed_bytes\": " << item.processed_bytes
               << ", \"blocks\": " << item.blocks
               << ", \"jobs\": " << item.jobs
               << ", \"factor_checks\": " << item.factor_checks
               << ", \"compared_bytes\": " << item.compared_bytes
               << ", \"direct_symbol_checks\": " << item.direct_symbol_checks
               << ", \"candidate_scratch_bytes\": "
               << item.candidate_scratch_bytes
               << ", \"padded_scratch_bytes\": "
               << item.padded_scratch_bytes
               << ", \"candidate_traffic_bytes\": "
               << item.candidate_traffic_bytes
               << ", \"padded_traffic_bytes\": "
               << item.padded_traffic_bytes << "}"
               << (i + 1 == correctness.size() ? "\n" : ",\n");
    }
    output << "  ],\n  \"benchmarks\": [\n";
    for (size_t i = 0; i < benchmarks.size(); ++i)
    {
        const BenchmarkResult& item = benchmarks[i];
        output << "    {\"field_bits\": " << item.field_bits
               << ", \"q\": " << item.q
               << ", \"parent\": " << item.parent
               << ", \"valid_bytes\": " << item.valid_bytes
               << ", \"processed_bytes\": " << item.processed_bytes
               << ", \"batch\": " << item.batch
               << ", \"reuse\": " << item.reuse
               << ", \"samples\": " << item.samples
               << ", \"setup_median_us\": " << item.setup_median_us
               << ", \"setup_mad_us\": " << item.setup_mad_us
               << ", \"candidate_median_us\": " << item.candidate_median_us
               << ", \"candidate_mad_us\": " << item.candidate_mad_us
               << ", \"padded_median_us\": " << item.padded_median_us
               << ", \"padded_mad_us\": " << item.padded_mad_us
               << ", \"padded_over_candidate\": "
               << item.padded_over_candidate
               << ", \"credible_gain_percent\": "
               << item.credible_gain_percent
               << ", \"candidate_scratch_bytes_per_stripe\": "
               << item.candidate_scratch_bytes_per_stripe
               << ", \"padded_scratch_bytes_per_stripe\": "
               << item.padded_scratch_bytes_per_stripe
               << ", \"resident_batch_bytes\": " << item.resident_batch_bytes
               << ", \"candidate_traffic_bytes_per_execution\": "
               << item.candidate_traffic_bytes_per_execution
               << ", \"padded_traffic_bytes_per_execution\": "
               << item.padded_traffic_bytes_per_execution << "}"
               << (i + 1 == benchmarks.size() ? "\n" : ",\n");
    }
    output << "  ]\n}\n";
}

} // namespace

int main(int argc, char** argv)
{
    leo2_context* context = NULL;
    try
    {
        const Arguments arguments = parse_arguments(argc, argv);
        leo2_context_options options;
        std::memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        const leo2_result create_result = leo2_context_create(&options, &context);
        require(create_result == LEO2_SUCCESS,
            std::string("context create failed: ") +
                leo2_result_string(create_result));

        const BinaryField gf8 = make_gf8();
        const BinaryField gf16 = make_gf16();
        const std::vector<CorrectnessCase> cases = make_correctness_cases();
        uint64_t digest = UINT64_C(14695981039346656037);
        uint64_t direct_checks = 0;
        uint64_t factor_checks = 0;
        std::vector<CorrectnessResult> correctness;
        correctness.reserve(cases.size());
        for (size_t i = 0; i < cases.size(); ++i)
        {
            const BinaryField& field = cases[i].field_bits == 8 ? gf8 : gf16;
            const bool run_direct = cases[i].valid_bytes == 64 ||
                (cases[i].valid_bytes == 1 &&
                 (cases[i].q == 129 || cases[i].q == 257));
            CorrectnessResult result = verify_case(
                field, cases[i], run_direct, &digest);
            direct_checks += result.direct_symbol_checks;
            factor_checks += result.factor_checks;
            correctness.push_back(result);
        }

        std::vector<BenchmarkResult> benchmarks;
        volatile uint64_t sink = digest;
        if (arguments.mode != "correctness")
        {
            const std::vector<BenchmarkSpec> specs = make_benchmark_specs(
                arguments.mode == "backend");
            benchmarks.reserve(specs.size());
            for (size_t i = 0; i < specs.size(); ++i)
            {
                const BinaryField& field = specs[i].field_bits == 8 ? gf8 : gf16;
                benchmarks.push_back(benchmark_case(
                    field, specs[i], arguments.benchmark_scale, &sink));
            }
        }

        std::unique_ptr<std::ofstream> file;
        std::ostream* output = &std::cout;
        if (!arguments.output_path.empty())
        {
            file.reset(new std::ofstream(arguments.output_path.c_str(),
                std::ios::out | std::ios::trunc));
            require(file->good(), "could not open result path");
            output = file.get();
        }
        write_results(*output, arguments,
            leo2_context_backend(context), read_runtime_environment(),
            digest, direct_checks, factor_checks, correctness, benchmarks);
        leo2_context_destroy(context);
        (void)sink;
        return 0;
    }
    catch (const std::exception& error)
    {
        if (context)
            leo2_context_destroy(context);
        std::cerr << "FAIL c5_dyadic_cpp: " << error.what() << std::endl;
        return 1;
    }
}
