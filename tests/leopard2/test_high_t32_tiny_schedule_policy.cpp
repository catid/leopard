/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the conditions in the Leopard-RS
    BSD license are met.  See the repository LICENSE file.
*/

/*
    Production regression for the pass-local GF8/AVX2 native Algorithm 5
    tiny schedule policy.  Reusable AUTO plans retain their exact schedules,
    but qualified tiny multi-loss passes execute the mature regular kernels.
    Exact-byte one-shot plans omit schedules under the same measured policy.
    A forced-tiled codec is the same-executable pruned-schedule control.
*/

#include "Leopard2Direct.h"
#include "Leopard2Dispatch.h"
#include "LeopardFF8.h"
#include "leopard.h"
#include "leopard2.h"

#ifndef LEO2_ENABLE_TEST_HOOKS
#error "The high tiny schedule-policy test requires Leopard2 test hooks"
#endif

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

static const uint8_t kGuardByte = 0xd7;
static const uint8_t kUntouchedByte = 0xa5;
static const size_t kGuardPrefix = 65;
static const size_t kGuardSuffix = 67;

void require(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void require_result(leo2_result result, const std::string& operation)
{
    if (result != LEO2_SUCCESS)
    {
        throw std::runtime_error(operation + ": " +
            leo2_result_string(result));
    }
}

uint32_t mix(uint32_t value)
{
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    return value ^ (value >> 16);
}

uint32_t ceil_pow2(uint32_t value)
{
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

class AlignedBytes
{
public:
    explicit AlignedBytes(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes_ == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes_, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes_) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes_);
        require((reinterpret_cast<uintptr_t>(data_) &
                    (leo2_scratch_alignment() - 1U)) == 0,
            "scratch allocation is not aligned");
    }

    ~AlignedBytes()
    {
#if defined(_MSC_VER)
        _aligned_free(data_);
#else
        free(data_);
#endif
    }

    void* data() { return data_; }
    size_t size() const { return bytes_; }

private:
    AlignedBytes(const AlignedBytes&);
    AlignedBytes& operator=(const AlignedBytes&);

    void* data_;
    size_t bytes_;
};

class GuardedShards
{
public:
    typedef std::vector<std::vector<uint8_t> > Snapshot;

    GuardedShards(uint32_t count, size_t bytes, uint8_t payload)
        : storage_(count)
        , bytes_(bytes)
    {
        for (uint32_t i = 0; i < count; ++i)
        {
            storage_[i].assign(
                kGuardPrefix + bytes_ + kGuardSuffix, kGuardByte);
            memset(data(i), payload, bytes_);
            require((reinterpret_cast<uintptr_t>(data(i)) & 63U) != 0,
                "test shard unexpectedly has 64-byte alignment");
        }
    }

    uint32_t count() const
    {
        return static_cast<uint32_t>(storage_.size());
    }

    uint8_t* data(uint32_t index)
    {
        return &storage_[index][kGuardPrefix];
    }

    const uint8_t* data(uint32_t index) const
    {
        return &storage_[index][kGuardPrefix];
    }

    void fill_originals(uint32_t seed)
    {
        for (uint32_t shard = 0; shard < count(); ++shard)
        {
            for (size_t byte = 0; byte < bytes_; ++byte)
            {
                data(shard)[byte] = static_cast<uint8_t>(mix(
                    seed ^ shard * 0x9e3779b9u ^
                    static_cast<uint32_t>(byte * 257U)));
            }
        }
    }

    std::vector<const void*> const_pointers() const
    {
        std::vector<const void*> result(count());
        for (uint32_t i = 0; i < count(); ++i)
            result[i] = data(i);
        return result;
    }

    std::vector<void*> mutable_pointers()
    {
        std::vector<void*> result(count());
        for (uint32_t i = 0; i < count(); ++i)
            result[i] = data(i);
        return result;
    }

    Snapshot snapshot() const
    {
        Snapshot result(count());
        for (uint32_t i = 0; i < count(); ++i)
            result[i].assign(data(i), data(i) + bytes_);
        return result;
    }

    void require_matches(
        const Snapshot& expected,
        const std::string& operation) const
    {
        require(expected.size() == storage_.size(),
            operation + " snapshot count changed");
        for (uint32_t i = 0; i < count(); ++i)
        {
            require(expected[i].size() == bytes_ &&
                    memcmp(data(i), expected[i].data(), bytes_) == 0,
                operation + " payload changed");
        }
    }

    void require_payload(uint32_t index, uint8_t expected,
        const std::string& operation) const
    {
        for (size_t byte = 0; byte < bytes_; ++byte)
        {
            if (data(index)[byte] != expected)
                throw std::runtime_error(operation + " payload changed");
        }
    }

    void require_guards(const std::string& operation) const
    {
        for (uint32_t shard = 0; shard < count(); ++shard)
        {
            for (size_t i = 0; i < kGuardPrefix; ++i)
            {
                if (storage_[shard][i] != kGuardByte)
                    throw std::runtime_error(operation + " prefix guard changed");
            }
            const size_t suffix = kGuardPrefix + bytes_;
            for (size_t i = suffix; i < storage_[shard].size(); ++i)
            {
                if (storage_[shard][i] != kGuardByte)
                    throw std::runtime_error(operation + " suffix guard changed");
            }
        }
    }

private:
    std::vector<std::vector<uint8_t> > storage_;
    size_t bytes_;
};

class CodecOwner
{
public:
    CodecOwner() : codec_(NULL) {}
    ~CodecOwner() { leo2_codec_destroy(codec_); }
    leo2_codec** output() { return &codec_; }
    leo2_codec* get() const { return codec_; }

private:
    CodecOwner(const CodecOwner&);
    CodecOwner& operator=(const CodecOwner&);
    leo2_codec* codec_;
};

class PlanOwner
{
public:
    PlanOwner() : plan_(NULL) {}
    ~PlanOwner() { leo2_decode_plan_destroy(plan_); }
    leo2_decode_plan** output() { return &plan_; }
    leo2_decode_plan* get() const { return plan_; }

private:
    PlanOwner(const PlanOwner&);
    PlanOwner& operator=(const PlanOwner&);
    leo2_decode_plan* plan_;
};

enum RouteExpectation
{
    RouteRegularOnly,
    RoutePrunedOnly,
    RouteMixed
};

struct ByteCase
{
    size_t bytes;
    RouteExpectation automatic_route;
};

struct ExecutionObservation
{
    leopard::ff8::TestOnlyHighDecodeCounts counts;
    size_t scratch_bytes;
};

std::vector<uint32_t> clustered_missing(uint32_t losses)
{
    std::vector<uint32_t> result(losses);
    for (uint32_t i = 0; i < losses; ++i)
        result[i] = i + 1U;
    return result;
}

std::vector<uint32_t> striped_missing(uint32_t k, uint32_t losses)
{
    require(losses > 1 && losses < k,
        "striped pattern requires one retained original");
    std::vector<uint32_t> result(losses);
    for (uint32_t i = 0; i < losses; ++i)
    {
        result[i] = 1U + static_cast<uint32_t>(
            (static_cast<uint64_t>(i) * (k - 2U)) / (losses - 1U));
        if (i != 0)
            require(result[i - 1] < result[i], "striped indices repeated");
    }
    return result;
}

std::vector<uint32_t> random_missing(
    uint32_t k,
    uint32_t losses,
    uint32_t seed)
{
    std::vector<uint32_t> candidates(k - 1U);
    for (uint32_t i = 0; i + 1U < k; ++i)
        candidates[i] = i + 1U;
    uint32_t state = seed;
    for (size_t remaining = candidates.size(); remaining > 1; --remaining)
    {
        state = mix(state + static_cast<uint32_t>(remaining));
        const size_t other = state % remaining;
        std::swap(candidates[remaining - 1U], candidates[other]);
    }
    candidates.resize(losses);
    std::sort(candidates.begin(), candidates.end());
    return candidates;
}

size_t partial_output_block_count(
    uint32_t k,
    uint32_t r,
    const std::vector<uint32_t>& missing)
{
    const uint32_t t = ceil_pow2(r);
    const uint32_t parent = ceil_pow2(k + t);
    std::vector<uint32_t> requested(parent / t, 0);
    for (size_t i = 0; i < missing.size(); ++i)
    {
        require(missing[i] < k,
            "missing index exceeds the systematic coordinate count");
        if (i != 0)
            require(missing[i - 1] < missing[i],
                "missing pattern is not strictly increasing");
        ++requested[(t + missing[i]) / t];
    }

    size_t blocks = 0;
    for (size_t block = 1; block < requested.size(); ++block)
    {
        if (requested[block] == 0)
            continue;
        require(requested[block] < t,
            "test pattern contains a mature full output block");
        ++blocks;
    }
    return blocks;
}

const char* route_name(RouteExpectation route)
{
    switch (route)
    {
    case RouteRegularOnly: return "regular";
    case RoutePrunedOnly: return "pruned";
    case RouteMixed: return "mixed";
    }
    return "unknown";
}

void create_codec(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    uint32_t flags,
    bool force_native_high,
    CodecOwner& owner,
    const std::string& label)
{
    leo2_codec_options options = {};
    options.struct_size = sizeof(options);
    options.flags = flags;
    options.shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    require_result(leo2_codec_create(context, k, r,
        LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8, &options,
        owner.output()), label + " codec create");
    require(leo2_codec_padded_side(owner.get()) == ceil_pow2(r),
        label + " created the wrong padded redundancy side");
    require(leo2_codec_parent_count(owner.get()) ==
            ceil_pow2(k + ceil_pow2(r)),
        label + " parent count changed");
    if (force_native_high)
    {
        require_result(leo2_test_codec_set_decode_mode(owner.get(),
            LEO2_TEST_DECODE_FORCE_NATIVE_HIGH),
            label + " force native Algorithm 5");
    }
}

void create_plan(
    leo2_codec* codec,
    const std::vector<uint8_t>& original_present,
    const std::vector<uint8_t>& recovery_present,
    size_t output_blocks,
    PlanOwner& owner,
    const std::string& label)
{
    require_result(leo2_decode_plan_create(codec, original_present.data(),
        recovery_present.data(), owner.output()), label + " plan create");
    require(!leo2_test_decode_plan_uses_translated_low(owner.get()),
        label + " unexpectedly selected translated Algorithm 4");

    leopard2_internal::DecodePlanPrunedScheduleInfo info;
    require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(
            owner.get(), &info),
        label + " schedule introspection failed");
    require(info.low_input_plan_count == 0 &&
            info.low_output_plan_count == 0,
        label + " compiled low-profile schedules");
    require(info.high_input_plan_count != 0,
        label + " did not retain an exact input schedule");
    require(info.high_output_plan_count == output_blocks,
        label + " did not retain every partial output-block schedule");
}

void require_route_counts(
    const leopard::ff8::TestOnlyHighDecodeCounts& counts,
    RouteExpectation route,
    size_t expected_output_blocks,
    size_t expected_passes,
    const std::string& label)
{
    require(counts.output_blocks ==
            expected_output_blocks * expected_passes,
        label + " output-block/pass accounting changed");
    require(counts.compatibility_copy_fallbacks == 0,
        label + " entered the whole-block copy fallback");
    require(counts.syndrome_pruned_fallback_blocks == 0,
        label + " rejected a compiled inverse schedule");
    require(counts.pruned_output_blocks + counts.mature_output_blocks ==
            counts.output_blocks,
        label + " output route counters are inconsistent");

    if (route == RouteRegularOnly)
    {
        require(counts.pruned_output_blocks == 0 &&
                counts.mature_output_blocks == counts.output_blocks,
            label + " did not use only mature regular schedules");
        require(counts.syndrome_pruned_accumulated_blocks == 0,
            label + " consumed an exact inverse schedule while bypassed");
    }
    else if (route == RoutePrunedOnly)
    {
        require(counts.pruned_output_blocks == counts.output_blocks &&
                counts.mature_output_blocks == 0,
            label + " did not use only retained exact schedules");
    }
    else
    {
        require(counts.pruned_output_blocks != 0 &&
                counts.mature_output_blocks != 0,
            label + " did not expose both pass-local schedule routes");
    }
}

ExecutionObservation execute_decode(
    uint32_t k,
    size_t bytes,
    const std::string& label,
    leo2_decode_plan* plan,
    leopard2_internal::DecodePath expected_path,
    leopard2_internal::DecodePathRule expected_rule,
    size_t expected_work_slots,
    RouteExpectation expected_route,
    size_t output_blocks,
    const std::vector<uint32_t>& missing,
    const GuardedShards& originals,
    const GuardedShards::Snapshot& original_snapshot,
    const GuardedShards& recovery,
    const GuardedShards::Snapshot& recovery_snapshot)
{
    leopard2_internal::DecodePathInfo path;
    require_result(leopard2_internal::GetDecodePlanPathInfo(
        plan, bytes, false, &path), label + " path query");
    require(path.path == expected_path && path.rule == expected_rule,
        label + " selected an unexpected decode path");
    require(path.required_work_slots == expected_work_slots,
        label + " changed its work-slot geometry");

    std::vector<const void*> original = originals.const_pointers();
    const std::vector<const void*> parity = recovery.const_pointers();
    std::vector<uint8_t> is_missing(k, 0);
    for (size_t i = 0; i < missing.size(); ++i)
    {
        is_missing[missing[i]] = 1;
        original[missing[i]] = NULL;
    }

    GuardedShards restored(k, bytes, kUntouchedByte);
    std::vector<void*> output(k, NULL);
    for (size_t i = 0; i < missing.size(); ++i)
        output[missing[i]] = restored.data(missing[i]);

    size_t scratch_bytes = 0;
    require_result(leo2_decode_plan_scratch_size(
        plan, bytes, &scratch_bytes), label + " scratch query");
    require(scratch_bytes != 0, label + " unexpectedly needs no scratch");
    AlignedBytes scratch(scratch_bytes);

    leopard::ff8::TestOnlyResetHighDecodeCounts();
    require_result(leo2_decode_plan_execute(plan, bytes,
        original.data(), parity.data(), output.data(),
        scratch.data(), scratch.size()), label + " decode");

    ExecutionObservation observation;
    observation.counts = leopard::ff8::TestOnlyGetHighDecodeCounts();
    observation.scratch_bytes = scratch_bytes;
    const size_t aligned_bytes = bytes & ~static_cast<size_t>(63U);
    size_t execution_tile_count = 0;
    size_t maximum_pass_bytes = 0;
    require_result(leopard2_internal::GetDecodePlanExecutionTiles(
        plan, bytes, &execution_tile_count, &maximum_pass_bytes),
        label + " execution-tile query");
    const size_t rounded_bytes = (bytes + 63U) & ~static_cast<size_t>(63U);
    const bool fused_ragged_pass = aligned_bytes != bytes &&
        execution_tile_count == 1 && maximum_pass_bytes == rounded_bytes;
    const size_t expected_passes = execution_tile_count +
        (aligned_bytes != bytes && !fused_ragged_pass ? 1U : 0U);
    require_route_counts(observation.counts, expected_route,
        output_blocks, expected_passes, label);

    for (uint32_t shard = 0; shard < k; ++shard)
    {
        if (is_missing[shard])
        {
            require(memcmp(restored.data(shard), originals.data(shard),
                    bytes) == 0,
                label + " restored original mismatch");
        }
        else
            restored.require_payload(
                shard, kUntouchedByte, label + " unrequested output");
    }
    restored.require_guards(label + " restored");
    originals.require_matches(original_snapshot, label + " original input");
    recovery.require_matches(recovery_snapshot, label + " parity input");
    originals.require_guards(label + " original input");
    recovery.require_guards(label + " parity input");
    return observation;
}

void run_campaign(
    leo2_context* context,
    const char* name,
    uint32_t k,
    uint32_t r,
    const std::vector<uint32_t>& missing,
    const std::vector<ByteCase>& byte_cases,
    bool forced_tiled_control,
    bool force_native_high,
    bool automatic_tiled,
    uint32_t seed)
{
    require(!byte_cases.empty(), std::string(name) + " has no byte cases");
    const size_t maximum_bytes = std::max_element(
        byte_cases.begin(), byte_cases.end(),
        [](const ByteCase& left, const ByteCase& right) {
            return left.bytes < right.bytes;
        })->bytes;
    const uint32_t t = ceil_pow2(r);
    require(t == 16 || t == 32 || t == 64,
        std::string(name) + " escaped the tested transform sides");
    const uint32_t parent = ceil_pow2(k + t);
    const size_t output_blocks =
        partial_output_block_count(k, r, missing);

    CodecOwner automatic_codec;
    create_codec(context, k, r, 0, force_native_high,
        automatic_codec, std::string(name) + " AUTO");
    CodecOwner forced_codec;
    if (forced_tiled_control)
    {
        create_codec(context, k, r, LEO2_CODEC_FORCE_TILED_DECODE,
            force_native_high, forced_codec,
            std::string(name) + " forced-tiled");
    }

    GuardedShards originals(k, maximum_bytes, kUntouchedByte);
    originals.fill_originals(seed);
    GuardedShards recovery(r, maximum_bytes, kUntouchedByte);
    const GuardedShards::Snapshot original_snapshot = originals.snapshot();

    size_t encode_scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(automatic_codec.get(),
        maximum_bytes, &encode_scratch_bytes),
        std::string(name) + " encode scratch query");
    require(encode_scratch_bytes != 0,
        std::string(name) + " unexpectedly needs no encode scratch");
    AlignedBytes encode_scratch(encode_scratch_bytes);
    std::vector<const void*> original_input = originals.const_pointers();
    std::vector<void*> recovery_output = recovery.mutable_pointers();
    require_result(leo2_encode(automatic_codec.get(), maximum_bytes,
        original_input.data(), recovery_output.data(),
        encode_scratch.data(), encode_scratch.size()),
        std::string(name) + " encode");
    const GuardedShards::Snapshot recovery_snapshot = recovery.snapshot();
    originals.require_matches(
        original_snapshot, std::string(name) + " encode input");
    originals.require_guards(std::string(name) + " encode input");
    recovery.require_guards(std::string(name) + " encode output");

    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;

    PlanOwner automatic_plan;
    create_plan(automatic_codec.get(), original_present, recovery_present,
        output_blocks, automatic_plan, std::string(name) + " AUTO");
    PlanOwner forced_plan;
    if (forced_tiled_control)
    {
        create_plan(forced_codec.get(), original_present, recovery_present,
            output_blocks, forced_plan,
            std::string(name) + " forced-tiled");
    }

    const leopard2_internal::DecodePath automatic_path = automatic_tiled
        ? leopard2_internal::kDecodePathTiled
        : leopard2_internal::kDecodePathMaterialized;
    const leopard2_internal::DecodePathRule automatic_rule = automatic_tiled
        ? leopard2_internal::kDecodeRuleWorkspaceTiled
        : leopard2_internal::kDecodeRuleWorkspaceMaterialized;
    const size_t automatic_work_slots = automatic_tiled
        ? static_cast<size_t>(2U * t + missing.size())
        : parent;
    for (size_t i = 0; i < byte_cases.size(); ++i)
    {
        const ByteCase& bytes = byte_cases[i];
        const std::string automatic_label =
            std::string(name) + " AUTO B=" + std::to_string(bytes.bytes);
        const ExecutionObservation automatic = execute_decode(
            k, bytes.bytes, automatic_label, automatic_plan.get(),
            automatic_path, automatic_rule, automatic_work_slots,
            bytes.automatic_route, output_blocks, missing,
            originals, original_snapshot, recovery, recovery_snapshot);

        if (forced_tiled_control)
        {
            const std::string forced_label =
                std::string(name) + " forced-tiled B=" +
                std::to_string(bytes.bytes);
            const ExecutionObservation forced = execute_decode(
                k, bytes.bytes, forced_label, forced_plan.get(),
                leopard2_internal::kDecodePathTiled,
                leopard2_internal::kDecodeRuleForcedTiled,
                static_cast<size_t>(2U * t + missing.size()),
                RoutePrunedOnly, output_blocks, missing,
                originals, original_snapshot, recovery, recovery_snapshot);
            const bool automatic_fuses_ragged_pass =
                automatic_tiled && parent == 128 && t == 32 &&
                bytes.bytes > 64 && bytes.bytes < 512 &&
                (bytes.bytes & 63U) != 0;
            if (automatic_fuses_ragged_pass)
            {
                require(automatic.scratch_bytes > forced.scratch_bytes,
                    std::string(name) +
                        " fused AUTO scratch did not retain rounded rows");
            }
            else
            {
                require(automatic.scratch_bytes == forced.scratch_bytes,
                    std::string(name) +
                        " AUTO and forced-tiled exact scratch sizes differ");
            }
        }

        std::cout << "CASE " << name
                  << " K=" << k << " R=" << r
                  << " B=" << bytes.bytes
                  << " L=" << missing.size()
                  << " auto=" << route_name(bytes.automatic_route)
                  << (forced_tiled_control ? " control=pruned" : "")
                  << " scratch=" << automatic.scratch_bytes << std::endl;
    }
}

std::vector<ByteCase> one_byte_case(
    size_t bytes,
    RouteExpectation route)
{
    std::vector<ByteCase> result(1);
    result[0].bytes = bytes;
    result[0].automatic_route = route;
    return result;
}

size_t transient_high_schedule_count(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    const std::vector<uint32_t>& missing,
    size_t bytes,
    const std::string& label)
{
    CodecOwner codec;
    create_codec(context, k, r, 0, false, codec, label);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;

    PlanOwner plan;
    require_result(
        leopard2_internal::CreateOneShotTransformPlanForDiagnostics(
            codec.get(), bytes, original_present.data(),
            recovery_present.data(), plan.output()),
        label + " transient plan create");
    leopard2_internal::DecodePlanPrunedScheduleInfo info;
    require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(
            plan.get(), &info),
        label + " transient schedule introspection failed");
    require(info.low_input_plan_count == 0 &&
            info.low_output_plan_count == 0,
        label + " transient native-high plan compiled low schedules");
    return info.high_input_plan_count + info.high_output_plan_count;
}

size_t reusable_high_schedule_count(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    const std::vector<uint32_t>& missing,
    const std::string& label)
{
    CodecOwner codec;
    create_codec(context, k, r, 0, false, codec, label);
    std::vector<uint8_t> original_present(k, 1);
    std::vector<uint8_t> recovery_present(r, 1);
    for (size_t i = 0; i < missing.size(); ++i)
        original_present[missing[i]] = 0;

    PlanOwner plan;
    require_result(leo2_decode_plan_create(codec.get(),
            original_present.data(), recovery_present.data(), plan.output()),
        label + " reusable plan create");
    leopard2_internal::DecodePlanPrunedScheduleInfo info;
    require(leopard2_internal::GetDecodePlanPrunedScheduleInfo(
            plan.get(), &info),
        label + " reusable schedule introspection failed");
    require(info.low_input_plan_count == 0 &&
            info.low_output_plan_count == 0,
        label + " reusable native-high plan compiled low schedules");
    return info.high_input_plan_count + info.high_output_plan_count;
}

} // namespace

int main()
{
    try
    {
        require(leo_init() == Leopard_Success, "Leopard initialization");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AVX2;
        options.thread_count = 1;
        leo2_context* avx2_context = NULL;
        const leo2_result context_result =
            leo2_context_create(&options, &avx2_context);
        if (context_result == LEO2_UNSUPPORTED)
        {
            std::cout << "SKIP high_tiny_schedule_policy: "
                         "AVX2 unavailable" << std::endl;
            return 0;
        }
        require_result(context_result, "AVX2 context create");
        require(leo2_context_backend(avx2_context) == LEO2_BACKEND_AVX2,
            "explicit AVX2 context reported another backend");

        std::vector<ByteCase> byte_matrix;
        const size_t regular_sizes[] = {
            1, 63, 64, 65, 511, 512, 513, 575
        };
        for (size_t i = 0;
             i < sizeof(regular_sizes) / sizeof(regular_sizes[0]); ++i)
        {
            ByteCase entry = { regular_sizes[i], RouteRegularOnly };
            byte_matrix.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            byte_matrix.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            byte_matrix.push_back(entry);
        }
        run_campaign(avx2_context, "reused-striped-byte-matrix",
            192, 32, striped_missing(192, 32), byte_matrix,
            true, false, true, 0x54a32000U);

        run_campaign(avx2_context, "lower-k-clustered",
            33, 32, clustered_missing(32),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a32101U);
        run_campaign(avx2_context, "upper-k-random",
            224, 32, random_missing(224, 32, 0xb16483d2U),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a32202U);

        std::vector<ByteCase> route_pair;
        {
            const ByteCase entry = { 64, RouteRegularOnly };
            route_pair.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            route_pair.push_back(entry);
        }
        run_campaign(avx2_context, "clustered-pattern",
            192, 32, clustered_missing(32), route_pair,
            false, false, true, 0x54a32303U);
        run_campaign(avx2_context, "random-pattern",
            192, 32, random_missing(192, 32, 0x73c4e18bU), route_pair,
            false, false, true, 0x54a32404U);

        run_campaign(avx2_context, "k32-neighbor",
            32, 32, striped_missing(32, 31),
            one_byte_case(64, RoutePrunedOnly),
            false, true, false, 0x54a32505U);
        run_campaign(avx2_context, "t32-r29-punctured-byte-matrix",
            57, 29, random_missing(57, 29, 0x6eb496d3U), byte_matrix,
            true, false, true, 0x54a32606U);
        run_campaign(avx2_context, "t32-r17-lower-k",
            33, 17, clustered_missing(17),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a33313U);
        run_campaign(avx2_context, "t32-r31-upper-k",
            224, 31, random_missing(224, 31, 0x9f27a451U),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a33515U);
        run_campaign(avx2_context, "t32-parent-boundary-lower",
            96, 24, striped_missing(96, 24),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a33616U);
        run_campaign(avx2_context, "t32-parent-boundary-upper",
            97, 24, random_missing(97, 24, 0x41bc28d6U),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a33717U);
        run_campaign(avx2_context, "t32-r16-neighbor",
            20, 16, striped_missing(20, 16),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a33414U);
        run_campaign(avx2_context, "t32-punctured-partial-loss-neighbor",
            57, 29, random_missing(57, 28, 0x78a5dd93U),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a33818U);
        run_campaign(avx2_context, "t32-r32-l31",
            192, 32, striped_missing(192, 31),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a32707U);

        std::vector<ByteCase> sparse_partial_thresholds;
        {
            const ByteCase entry = { 64, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 320, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 321, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 384, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 385, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 512, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 513, RouteRegularOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            sparse_partial_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            sparse_partial_thresholds.push_back(entry);
        }
        run_campaign(avx2_context, "t32-sparse-partial-thresholds",
            57, 29, random_missing(57, 2, 0x148c2a91U),
            sparse_partial_thresholds,
            true, false, true, 0x54a33919U);
        run_campaign(avx2_context, "t32-three-loss-b320",
            57, 29, random_missing(57, 3, 0x778bd120U),
            one_byte_case(320, RouteRegularOnly),
            true, false, true, 0x54a33c1cU);
        run_campaign(avx2_context, "t32-four-loss-b320",
            57, 29, random_missing(57, 4, 0x82e4b7a3U),
            one_byte_case(320, RouteRegularOnly),
            true, false, true, 0x54a33d1dU);

        std::vector<ByteCase> half_loss_thresholds;
        {
            const ByteCase entry = { 384, RouteRegularOnly };
            half_loss_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 448, RouteRegularOnly };
            half_loss_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 512, RouteRegularOnly };
            half_loss_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            half_loss_thresholds.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            half_loss_thresholds.push_back(entry);
        }
        run_campaign(avx2_context, "t32-half-loss-thresholds",
            57, 29, random_missing(57, 15, 0x826ba7d4U),
            half_loss_thresholds,
            true, false, true, 0x54a33a1aU);
        run_campaign(avx2_context, "t32-floor-half",
            57, 29, random_missing(57, 14, 0x1f3278c9U),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a33b1bU);

        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
            "select production one-shot setup policy");
        const std::vector<uint32_t> paper_t32_missing =
            random_missing(224, 30, 0xc4b62240U);
        require(transient_high_schedule_count(avx2_context,
                224, 32, paper_t32_missing, 1024,
                "transient R10 T32 B1024") == 0,
            "R10 T32 B1024 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                224, 32, paper_t32_missing, 1025,
                "transient R10 T32 B1025") != 0,
            "R10 T32 B1025 transient plan omitted schedules");
        require(reusable_high_schedule_count(avx2_context,
                224, 32, paper_t32_missing,
                "reusable R10 T32") != 0,
            "R10 T32 reusable plan lost exact schedules");
        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(2),
            "select attribution one-shot setup policy");
        require(transient_high_schedule_count(avx2_context,
                224, 32, paper_t32_missing, 1024,
                "attribution R10 T32 B1024") != 0,
            "attribution R10 T32 B1024 plan omitted schedules");
        require(leopard2_internal::SetOneShotPlanSetupModeForDiagnostics(3),
            "restore production one-shot setup policy");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da19U), 64,
                "transient sparse B64") == 0,
            "eligible sparse B64 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da1aU), 320,
                "transient sparse B320") == 0,
            "eligible sparse B320 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da1bU), 384,
                "transient sparse B384") == 0,
            "eligible one-shot sparse B384 plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 15, 0xc4b6da1cU), 512,
                "transient half-loss B512") == 0,
            "eligible half-loss B512 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 14, 0xc4b6da1dU), 512,
                "transient floor-half B512") == 0,
            "eligible one-shot floor-half B512 plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 15, 0xc4b6da1eU), 513,
                "transient ragged B513") == 0,
            "all-regular ragged B513 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da20U), 575,
                "transient ragged B575") == 0,
            "all-regular ragged B575 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da21U), 576,
                "transient B576") != 0,
            "pruned B576 transient plan omitted schedules");
        require(transient_high_schedule_count(avx2_context,
                57, 29, random_missing(57, 2, 0xc4b6da22U), 577,
                "transient mixed B577") != 0,
            "mixed B577 transient plan omitted schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 25, 0xc4b6da1fU), 64,
                "transient T64 partial B64") == 0,
            "eligible T64 partial transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da23U), 319,
                "transient T64 sparse B319") == 0,
            "eligible T64 sparse B319 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da24U), 320,
                "transient T64 sparse B320") == 0,
            "eligible T64 sparse B320 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da27U), 321,
                "transient T64 sparse B321") == 0,
            "all-regular T64 sparse B321 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da28U), 384,
                "transient T64 sparse B384") == 0,
            "eligible T64 sparse B384 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da29U), 385,
                "transient T64 sparse B385") == 0,
            "all-regular T64 sparse B385 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2aU), 512,
                "transient T64 sparse B512") == 0,
            "one-shot-only T64 sparse B512 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2bU), 575,
                "transient T64 sparse B575") == 0,
            "one-shot-only T64 sparse B575 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2cU), 576,
                "transient T64 sparse B576") == 0,
            "one-shot-only T64 sparse B576 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2dU), 4096,
                "transient T64 sparse B4096") == 0,
            "one-shot-only T64 sparse B4096 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2eU), 4097,
                "transient T64 sparse B4097") == 0,
            "one-shot-only T64 sparse B4097 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 2, 0xc4b6da2fU), 4160,
                "transient T64 sparse B4160") != 0,
            "pruned T64 sparse B4160 transient plan omitted schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 25, 0xc4b6da25U), 575,
                "transient T64 half-loss B575") == 0,
            "eligible T64 half-loss B575 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 25, 0xc4b6da26U), 576,
                "transient T64 half-loss B576") == 0,
            "one-shot-only T64 half-loss B576 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 25, 0xc4b6da30U), 4096,
                "transient T64 half-loss B4096") == 0,
            "one-shot-only T64 half-loss B4096 transient plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                99, 50, random_missing(99, 25, 0xc4b6da31U), 4160,
                "transient T64 half-loss B4160") != 0,
            "pruned T64 half-loss B4160 transient plan omitted schedules");
        require(transient_high_schedule_count(avx2_context,
                192, 62, random_missing(192, 2, 0xc4b6da32U), 448,
                "transient K192 T64 sparse B448") == 0,
            "eligible K192 T64 sparse B448 plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                192, 62, random_missing(192, 2, 0xc4b6da33U), 512,
                "transient K192 T64 sparse B512") == 0,
            "one-shot-only K192 T64 sparse B512 plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                192, 62, random_missing(192, 31, 0xc4b6da34U), 4096,
                "transient K192 T64 half-loss B4096") == 0,
            "one-shot-only K192 T64 B4096 plan retained schedules");
        require(transient_high_schedule_count(avx2_context,
                192, 62, random_missing(192, 31, 0xc4b6da35U), 4160,
                "transient K192 T64 half-loss B4160") != 0,
            "pruned K192 T64 B4160 transient plan omitted schedules");

        std::vector<ByteCase> t64_byte_matrix;
        const size_t t64_regular_sizes[] = {
            1, 63, 64, 65, 127, 128, 129, 511, 512, 513, 575
        };
        for (size_t i = 0;
             i < sizeof(t64_regular_sizes) / sizeof(t64_regular_sizes[0]);
             ++i)
        {
            const ByteCase entry = {
                t64_regular_sizes[i], RouteRegularOnly
            };
            t64_byte_matrix.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            t64_byte_matrix.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            t64_byte_matrix.push_back(entry);
        }
        run_campaign(avx2_context, "t64-striped-byte-matrix",
            99, 50, striped_missing(99, 50), t64_byte_matrix,
            true, false, true, 0x54a32909U);

        std::vector<ByteCase> t64_sparse_partial;
        const size_t t64_sparse_regular_sizes[] = {
            1, 63, 64, 65, 127, 128, 129, 255, 256, 257, 319
        };
        for (size_t i = 0;
             i < sizeof(t64_sparse_regular_sizes) /
                    sizeof(t64_sparse_regular_sizes[0]);
             ++i)
        {
            const ByteCase entry = {
                t64_sparse_regular_sizes[i], RouteRegularOnly
            };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 320, RouteRegularOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 321, RouteRegularOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 383, RouteRegularOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 384, RouteRegularOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 385, RouteRegularOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 512, RoutePrunedOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 513, RouteMixed };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 575, RouteMixed };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            t64_sparse_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            t64_sparse_partial.push_back(entry);
        }
        run_campaign(avx2_context, "t64-sparse-partial-thresholds",
            99, 50, random_missing(99, 2, 0x2377cb51U),
            t64_sparse_partial,
            true, false, true, 0x54a33e1eU);

        std::vector<ByteCase> t64_half_partial;
        const size_t t64_half_regular_sizes[] = {
            320, 384, 448, 512, 513, 575
        };
        for (size_t i = 0;
             i < sizeof(t64_half_regular_sizes) /
                    sizeof(t64_half_regular_sizes[0]);
             ++i)
        {
            const ByteCase entry = {
                t64_half_regular_sizes[i], RouteRegularOnly
            };
            t64_half_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 576, RoutePrunedOnly };
            t64_half_partial.push_back(entry);
        }
        {
            const ByteCase entry = { 577, RouteMixed };
            t64_half_partial.push_back(entry);
        }
        run_campaign(avx2_context, "t64-half-partial-thresholds",
            99, 50, random_missing(99, 25, 0x427f89d3U),
            t64_half_partial,
            true, false, true, 0x54a33f1fU);
        run_campaign(avx2_context, "t64-lower-k-clustered",
            65, 33, clustered_missing(33),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a32a0aU);
        run_campaign(avx2_context, "t64-upper-k-random",
            124, 62, random_missing(124, 62, 0x60c44291U),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a32b0bU);
        run_campaign(avx2_context, "t64-r63-neighbor",
            124, 63, striped_missing(124, 63),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a32c0cU);
        run_campaign(avx2_context, "t64-r64-neighbor",
            124, 64, striped_missing(124, 64),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a32d0dU);
        run_campaign(avx2_context, "t64-r63-partial",
            124, 63, striped_missing(124, 31),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a34020U);
        run_campaign(avx2_context, "t64-r64-partial",
            124, 64, striped_missing(124, 32),
            one_byte_case(512, RouteRegularOnly),
            true, false, true, 0x54a34121U);
        run_campaign(avx2_context, "t64-k64-neighbor",
            64, 50, striped_missing(64, 50),
            one_byte_case(64, RoutePrunedOnly),
            false, true, false, 0x54a32e0eU);
        run_campaign(avx2_context, "t64-upper-parent-regular",
            125, 62, striped_missing(125, 62),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a32f0fU);
        run_campaign(avx2_context, "t64-k191-regular",
            191, 62, random_missing(191, 62, 0x7337ca09U),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a33010U);
        run_campaign(avx2_context, "t64-k192-neighbor",
            192, 62, striped_missing(192, 62),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a33111U);
        run_campaign(avx2_context, "t64-partial-loss-regular",
            99, 50, striped_missing(99, 49),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a33212U);
        run_campaign(avx2_context, "t64-k192-partial-neighbor",
            192, 62, striped_missing(192, 31),
            one_byte_case(64, RouteRegularOnly),
            true, false, true, 0x54a34222U);
        {
            std::vector<ByteCase> k192_partial_boundaries;
            const ByteCase b448 = { 448, RouteRegularOnly };
            const ByteCase b512 = { 512, RouteRegularOnly };
            const ByteCase b4096 = { 4096, RoutePrunedOnly };
            k192_partial_boundaries.push_back(b448);
            k192_partial_boundaries.push_back(b512);
            k192_partial_boundaries.push_back(b4096);
            run_campaign(avx2_context, "t64-k192-half-boundaries",
                192, 62, random_missing(192, 31, 0xa659b17cU),
                k192_partial_boundaries,
                true, false, true, 0x54a34424U);
        }
        {
            std::vector<ByteCase> k192_sparse_boundaries;
            const ByteCase b448 = { 448, RouteRegularOnly };
            const ByteCase b512 = { 512, RoutePrunedOnly };
            k192_sparse_boundaries.push_back(b448);
            k192_sparse_boundaries.push_back(b512);
            run_campaign(avx2_context, "t64-k192-sparse-boundaries",
                192, 62, random_missing(192, 2, 0x9f0cb2d8U),
                k192_sparse_boundaries,
                true, false, true, 0x54a34525U);
        }

        leo2_context_options scalar_options = {};
        scalar_options.struct_size = sizeof(scalar_options);
        scalar_options.backend = LEO2_BACKEND_SCALAR;
        scalar_options.thread_count = 1;
        leo2_context* scalar_context = NULL;
        require_result(leo2_context_create(&scalar_options, &scalar_context),
            "scalar context create");
        run_campaign(scalar_context, "scalar-backend-neighbor",
            192, 32, striped_missing(192, 32),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a32808U);
        run_campaign(scalar_context, "scalar-t64-partial-neighbor",
            99, 50, striped_missing(99, 25),
            one_byte_case(64, RoutePrunedOnly),
            false, false, true, 0x54a34323U);
        leo2_context_destroy(scalar_context);
        leo2_context_destroy(avx2_context);

        std::cout << "PASS high_tiny_schedule_policy" << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL high_tiny_schedule_policy: "
                  << error.what() << std::endl;
        return 1;
    }
}
