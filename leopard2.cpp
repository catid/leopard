/*
    Copyright (c) 2017 Christopher A. Taylor.  All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the name of Leopard-RS nor the names of its contributors may be
      used to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
*/

#include "leopard2.h"

#include "LeopardCommon.h"
#include "Leopard2Backend.h"
#include "LeopardFF8.h"
#include "LeopardFF16.h"
#include "Leopard2Dispatch.h"
#include "Leopard2Direct.h"
#include "leopard.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <new>
#include <string.h>
#include <thread>
#include <utility>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__linux__)
#include <sched.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

class leo2_thread_pool;

struct leo2_context
{
    leo2_backend backend;
    const leopard::backend::Ops* ops;
    uint32_t thread_count;
    leo2_thread_pool* pool;
};

#ifdef LEO2_ENABLE_TEST_HOOKS
namespace leopard { namespace backend {
void TestSetContextOps(leo2_context* context, const Ops* ops)
{
    if (context && ops)
        context->ops = ops;
}
}} // namespace leopard::backend
#endif

struct leo2_codec
{
    leo2_context* context;
    uint32_t original_count;
    uint32_t recovery_count;
    uint32_t parent_count;
    uint32_t padded_side;
    uint32_t parent_dimension;
    leo2_profile profile;
    leo2_field field;
    leo2_shard_layout shard_layout;
    uint32_t flags;
    std::vector<uint8_t> permanent_erased;
    std::vector<uint8_t> permanent_locator8;
    std::vector<uint16_t> permanent_locator16;
    std::vector<uint8_t> fixed_factors8;
    std::vector<uint16_t> fixed_factors16;
    // Barycentric weights for the public systematic coordinates.  These are
    // populated only for the rigorously bounded direct-repair dispatch region.
    std::vector<uint8_t> direct_barycentric8;
    std::vector<uint16_t> direct_barycentric16;
    // Exact row-major generator coefficients, represented as Leopard field
    // logarithms, for the bounded allocation-free direct encoder.
    std::vector<uint8_t> direct_generator_logs8;
    std::vector<uint16_t> direct_generator_logs16;
#ifdef LEO2_ENABLE_TEST_HOOKS
    leo2_test_encode_mode test_encode_mode;
#endif
};

struct leo2_direct_repair_term
{
    // The high bit distinguishes recovery shards from original shards.
    uint32_t tagged_source;
    // Logarithmic fixed multiplier in the selected Leopard field.  Zero is
    // the multiplicative identity and is executed as copy/XOR.
    uint16_t multiplier_log;
};

struct leo2_decode_plan
{
    const leo2_codec* codec;
    std::vector<uint8_t> original_present;
    std::vector<uint8_t> recovery_present;
    std::vector<uint8_t> coordinate_erased;
    std::vector<uint8_t> locator8;
    std::vector<uint16_t> locator16;
    std::vector<uint32_t> missing_originals;
    // Sorted parent coordinates restored by this plan.  The compact transform
    // schedules below are derived once from this immutable list.
    std::vector<uint32_t> requested_coordinates;
    std::vector<uint64_t> generic_output_dependencies;
    std::vector<uint64_t> specialized_output_dependencies;
    std::vector<uint16_t> block_input_counts;
    std::vector<leopard2_internal::DecodeOutputBlock> high_output_blocks;
    // Algorithm 5 schedules are compiled during immutable plan setup.  Input
    // entries exist only for incomplete nonempty T-blocks that actually save
    // work; output entries align one-for-one with high_output_blocks, with an
    // empty plan selecting the mature full transform fallback.
    std::vector<leopard2_internal::PrunedTransformBlock>
        high_pruned_input_blocks;
    std::vector<leopard2_internal::PrunedTransformPlan>
        high_pruned_output_plans;
    std::vector<size_t> direct_term_offsets;
    std::vector<leo2_direct_repair_term> direct_terms;
    uint32_t generic_input_count;
    uint32_t direct_copy_recovery;
    uint32_t missing_original_count;
    bool no_op;
    bool direct_xor;
    bool direct_copy;
    bool direct_repair;
};

static leo2_result DecodePlanExecuteInternal(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes,
    bool multi_item_batch);

namespace {

static const size_t kScratchAlignment = 64;
static const uint32_t kGF8Order = 256;
static const uint32_t kGF16Order = 65536;
static const uint32_t kDirectRecoveryTag = 0x80000000u;
static const uint32_t kDirectMaxOriginals = 16;
static const uint32_t kDirectMaxRecoveries = 16;
static const uint32_t kDirectMaxLosses = 4;
static const uint32_t kDirectMaxParentDimension = 256;
static const uint64_t kDirectSimdTileBytes = 64;
static const uint64_t kDirectMinimumMeasuredBytes = 1024;

typedef leo2_result (*BatchTaskFunction)(void* context, size_t index);

#ifdef LEO2_ENABLE_TEST_HOOKS
static std::atomic<bool> g_test_thread_start_fault(false);
static std::atomic<unsigned> g_test_thread_start_fault_consumptions(0);

static bool TestConsumeThreadStartFault()
{
    if (!g_test_thread_start_fault.exchange(false, std::memory_order_acq_rel))
        return false;
    g_test_thread_start_fault_consumptions.fetch_add(
        1, std::memory_order_relaxed);
    return true;
}
#endif

#if defined(__linux__)
static uint32_t LinuxAllowedThreadCount()
{
    // cpu_set_t covers the ordinary case.  Retry with progressively larger,
    // stack-owned masks so a process restricted to high-numbered CPUs on a
    // very large host does not silently fall back to the online CPU count.
    static const size_t kMaximumAffinityBytes = 16 * 1024;
    unsigned long mask[kMaximumAffinityBytes / sizeof(unsigned long)];
    for (size_t bytes = sizeof(cpu_set_t);
         bytes <= kMaximumAffinityBytes;
         bytes *= 2)
    {
        memset(mask, 0, bytes);
        if (sched_getaffinity(
                0, bytes, reinterpret_cast<cpu_set_t*>(mask)) == 0)
        {
            uint32_t allowed = 0;
            const size_t word_count = bytes / sizeof(unsigned long);
            for (size_t i = 0; i < word_count && allowed < 128; ++i)
            {
                unsigned long word = mask[i];
                while (word != 0 && allowed < 128)
                {
                    word &= word - 1;
                    ++allowed;
                }
            }
            return allowed;
        }
        if (errno != EINVAL || bytes == kMaximumAffinityBytes)
            break;
    }
    return 0;
}
#endif

static uint32_t DefaultThreadCount()
{
    uint32_t allowed = 0;
#if defined(_WIN32)
    DWORD_PTR process_mask = 0;
    DWORD_PTR system_mask = 0;
    if (GetProcessAffinityMask(
            GetCurrentProcess(), &process_mask, &system_mask))
    {
        while (process_mask != 0)
        {
            allowed += static_cast<uint32_t>(process_mask & 1u);
            process_mask >>= 1;
        }
    }
#elif defined(__linux__)
    allowed = LinuxAllowedThreadCount();
#endif
    if (allowed == 0)
        allowed = static_cast<uint32_t>(std::thread::hardware_concurrency());
    if (allowed == 0)
        allowed = 1;
    return std::min<uint32_t>(allowed, 128);
}

} // namespace

class leo2_thread_pool
{
public:
    explicit leo2_thread_pool(uint32_t worker_count)
        : worker_count_(worker_count)
        , start_failed_(false)
        , stopping_(false)
        , generation_(0)
        , function_(NULL)
        , function_context_(NULL)
        , task_count_(0)
        , next_task_(0)
        , failure_(0)
        , completed_workers_(0)
    {}

    ~leo2_thread_pool()
    {
        Stop();
    }

#ifdef LEO2_ENABLE_TEST_HOOKS
    uint32_t TestWorkerCount()
    {
        std::lock_guard<std::mutex> lock(run_mutex_);
        return static_cast<uint32_t>(workers_.size());
    }
#endif

    leo2_result Run(
        size_t task_count,
        BatchTaskFunction function,
        void* function_context)
    {
        if (task_count == 0)
            return LEO2_SUCCESS;
        if (task_count == 1)
        {
            return function(function_context, 0);
        }

        /* One context supports concurrent codecs, but one batch owns its pool
           at a time.  Independent callers remain correct and are serialized at
           this optional scheduling layer rather than racing shared job state. */
        std::lock_guard<std::mutex> run_lock(run_mutex_);
        if (!EnsureStarted(task_count))
            return LEO2_OUT_OF_MEMORY;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            function_ = function;
            function_context_ = function_context;
            task_count_ = task_count;
            next_task_.store(0, std::memory_order_relaxed);
            failure_.store(PackFailure(task_count, LEO2_SUCCESS), std::memory_order_relaxed);
            completed_workers_ = 0;
            ++generation_;
        }
        work_ready_.notify_all();

#ifdef _OPENMP
        const int saved_openmp_threads = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        ExecuteTasks();
#ifdef _OPENMP
        omp_set_num_threads(saved_openmp_threads);
#endif

        {
            std::unique_lock<std::mutex> lock(mutex_);
            work_done_.wait(lock, [this]() {
                return completed_workers_ == workers_.size();
            });
            function_ = NULL;
            function_context_ = NULL;
            task_count_ = 0;
        }
        return UnpackFailureResult(failure_.load(std::memory_order_relaxed));
    }

private:
    bool EnsureStarted(size_t task_count)
    {
        const size_t target_workers = std::min<size_t>(
            worker_count_, task_count - 1);
        if (workers_.size() >= target_workers)
            return true;
        if (start_failed_ || worker_count_ == 0)
            return false;
#ifdef LEO2_ENABLE_TEST_HOOKS
        if (TestConsumeThreadStartFault())
        {
            start_failed_ = true;
            return false;
        }
#endif
        try
        {
            workers_.reserve(worker_count_);
            uint64_t seen_generation = 0;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                seen_generation = generation_;
            }
            while (workers_.size() < target_workers)
                workers_.push_back(std::thread(
                    &leo2_thread_pool::Worker, this, seen_generation));
        }
        catch (...)
        {
            start_failed_ = true;
            Stop();
            return false;
        }
        return true;
    }

    static uint64_t PackFailure(size_t index, leo2_result result)
    {
        const uint64_t bounded_index = index > 0xffffffffu ? 0xffffffffu : index;
        return (bounded_index << 32) | static_cast<uint32_t>(result);
    }

    static leo2_result UnpackFailureResult(uint64_t packed)
    {
        return static_cast<leo2_result>(static_cast<int32_t>(packed));
    }

    void RecordFailure(size_t index, leo2_result result)
    {
        if (result == LEO2_SUCCESS)
            return;
        uint64_t candidate = PackFailure(index, result);
        uint64_t current = failure_.load(std::memory_order_relaxed);
        while ((candidate >> 32) < (current >> 32) &&
               !failure_.compare_exchange_weak(
                   current, candidate, std::memory_order_relaxed, std::memory_order_relaxed))
        {}
    }

    void ExecuteTasks()
    {
        for (;;)
        {
            const size_t index = next_task_.fetch_add(1, std::memory_order_relaxed);
            if (index >= task_count_)
                break;
            RecordFailure(index, function_(function_context_, index));
        }
    }

    void Worker(uint64_t seen_generation)
    {
#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        for (;;)
        {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                work_ready_.wait(lock, [this, seen_generation]() {
                    return stopping_ || generation_ != seen_generation;
                });
                if (stopping_)
                    return;
                seen_generation = generation_;
            }
            ExecuteTasks();
            {
                std::lock_guard<std::mutex> lock(mutex_);
                ++completed_workers_;
                if (completed_workers_ == workers_.size())
                    work_done_.notify_one();
            }
        }
    }

    void Stop()
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        work_ready_.notify_all();
        for (size_t i = 0; i < workers_.size(); ++i)
            if (workers_[i].joinable())
                workers_[i].join();
        workers_.clear();
    }

    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::mutex run_mutex_;
    std::condition_variable work_ready_;
    std::condition_variable work_done_;
    uint32_t worker_count_;
    bool start_failed_;
    bool stopping_;
    uint64_t generation_;
    BatchTaskFunction function_;
    void* function_context_;
    size_t task_count_;
    std::atomic<size_t> next_task_;
    std::atomic<uint64_t> failure_;
    size_t completed_workers_;
};

namespace {

struct AddressRange
{
    uintptr_t begin;
    uintptr_t end;
};

struct ScratchLayout
{
    size_t range_offset;
    size_t pointer_offset;
    size_t data_offset;
    size_t total_bytes;
};

struct DecodeScratchGeometry
{
    ScratchLayout layout;
    size_t rounded_bytes;
    size_t aligned_prefix_bytes;
    size_t tail_bytes;
    size_t work_slot_count;
    size_t work_slot_bytes;
    size_t work_data_offset;
};

struct EncodeScratchGeometry
{
    ScratchLayout layout;
    size_t rounded_bytes;
    size_t aligned_prefix_bytes;
    size_t tail_bytes;
    size_t work_count;
    size_t work_slot_bytes;
    size_t work_data_offset;
};

static bool CheckedAdd(size_t a, size_t b, size_t& result)
{
    if (a > std::numeric_limits<size_t>::max() - b)
        return false;
    result = a + b;
    return true;
}

static bool CheckedMultiply(size_t a, size_t b, size_t& result)
{
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a)
        return false;
    result = a * b;
    return true;
}

static bool AlignUp(size_t value, size_t alignment, size_t& result)
{
    const size_t mask = alignment - 1;
    size_t sum = 0;
    if (!CheckedAdd(value, mask, sum))
        return false;
    result = sum & ~mask;
    return true;
}

static bool RoundShardBytes(uint64_t shard_bytes, size_t& rounded)
{
    if (shard_bytes == 0 || shard_bytes > std::numeric_limits<size_t>::max())
        return false;
    return AlignUp(static_cast<size_t>(shard_bytes), kScratchAlignment, rounded);
}

static leo2_result ResolveWireShardBytes(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    uint64_t& wire_bytes)
{
    if (!codec || payload_bytes == 0)
        return LEO2_INVALID_ARGUMENT;

    if (codec->shard_layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
    {
        if (codec->field != LEO2_FIELD_GF16)
            return LEO2_INTERNAL_ERROR;
        if ((payload_bytes & 1u) == 0)
            return LEO2_UNSUPPORTED;
        if (payload_bytes == std::numeric_limits<uint64_t>::max())
            return LEO2_INVALID_ARGUMENT;
        wire_bytes = payload_bytes + 1;
    }
    else
    {
        if (codec->shard_layout != LEO2_SHARD_LAYOUT_NATIVE_V1)
            return LEO2_INTERNAL_ERROR;
        if (codec->field == LEO2_FIELD_GF16 && (payload_bytes & 1u) != 0)
            return LEO2_UNSUPPORTED;
        wire_bytes = payload_bytes;
    }

    size_t rounded = 0;
    if (!RoundShardBytes(wire_bytes, rounded))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
}

static bool ComputeScratchLayout(
    size_t range_count,
    size_t pointer_count,
    size_t slot_count,
    size_t rounded_bytes,
    ScratchLayout& layout)
{
    size_t ranges_bytes = 0;
    size_t pointers_bytes = 0;
    size_t slots_bytes = 0;
    if (!CheckedMultiply(range_count, sizeof(AddressRange), ranges_bytes) ||
        !CheckedMultiply(pointer_count, sizeof(void*), pointers_bytes) ||
        !CheckedMultiply(slot_count, rounded_bytes, slots_bytes))
        return false;

    layout.range_offset = 0;
    if (!AlignUp(ranges_bytes, alignof(void*), layout.pointer_offset))
        return false;
    size_t end_pointers = 0;
    if (!CheckedAdd(layout.pointer_offset, pointers_bytes, end_pointers) ||
        !AlignUp(end_pointers, kScratchAlignment, layout.data_offset) ||
        !CheckedAdd(layout.data_offset, slots_bytes, layout.total_bytes))
        return false;
    return true;
}

static bool ComputeSplitScratchLayout(
    size_t range_count,
    size_t pointer_count,
    size_t input_slot_count,
    size_t work_slot_count,
    size_t work_slot_bytes,
    ScratchLayout& layout,
    size_t& work_data_offset)
{
    if (!ComputeScratchLayout(
            range_count, pointer_count, 0, 1, layout))
        return false;

    size_t input_bytes = 0;
    size_t work_bytes = 0;
    if (!CheckedMultiply(
            input_slot_count, kScratchAlignment, input_bytes) ||
        !CheckedAdd(layout.data_offset, input_bytes, work_data_offset) ||
        !CheckedMultiply(work_slot_count, work_slot_bytes, work_bytes) ||
        !CheckedAdd(work_data_offset, work_bytes, layout.total_bytes))
        return false;
    return true;
}

static uint32_t CeilPow2(uint64_t value)
{
    if (value == 0 || value > 65536)
        return 0;
    uint32_t result = 1;
    while (result < value)
        result <<= 1;
    return result;
}

static uint32_t CoordinateForOriginal(const leo2_codec* codec, uint32_t index)
{
    return codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? codec->padded_side + index
        : index;
}

static uint32_t CoordinateForRecovery(const leo2_codec* codec, uint32_t index)
{
    return codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? index
        : codec->padded_side + index;
}

static bool PrepareCodecPrunedTransform(
    const leo2_codec* codec,
    uint32_t size,
    uint32_t shift,
    bool inverse,
    const uint8_t* input_mask,
    const uint8_t* output_mask,
    leopard2_internal::PrunedTransformPlan& result)
{
#ifdef LEO_HAS_FF8
    if (codec->field == LEO2_FIELD_GF8)
        return leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, result);
#endif
#ifdef LEO_HAS_FF16
    if (codec->field == LEO2_FIELD_GF16)
        return leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, result);
#endif
    return false;
}

static bool PrunedPlanSavesByteHeavyWork(
    const leopard2_internal::PrunedTransformPlan& plan)
{
    return plan.size >= 2 &&
        (plan.operations.size() < plan.full_butterfly_count ||
         plan.input_zero_specializations != 0 ||
         plan.one_output_butterflies != 0);
}

static bool PreparePlanExecutionMetadata(leo2_decode_plan* plan)
{
    const leo2_codec* codec = plan->codec;
    const bool force_generic =
        (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
    const bool force_specialized =
        (codec->flags & (LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
                         LEO2_CODEC_FORCE_TILED_DECODE |
                         LEO2_CODEC_FORCE_MATERIALIZED_DECODE)) != 0;
    // AUTO owns both immutable paths because dispatch may depend on shard
    // bytes, backend, or a future offline table.  A forced diagnostic path
    // pays only for the metadata it can execute.
    const bool needs_generic = !force_specialized;
    const bool needs_specialized = !force_generic;

    if (plan->requested_coordinates.empty())
        return false;
    for (size_t i = 1; i < plan->requested_coordinates.size(); ++i)
        if (plan->requested_coordinates[i - 1] >= plan->requested_coordinates[i])
            return false;

    plan->generic_input_count = 0;
    std::vector<uint8_t> selected_coordinates;
    if (needs_specialized && codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        selected_coordinates.assign(codec->parent_count, 0);
    if (needs_specialized)
    {
        // A nontrivial parent has at least two blocks, so P/T <= 32768 and a
        // 16-bit prefix stores every possible local count exactly.
        if (codec->padded_side < 2 || codec->padded_side > 32768 ||
            codec->parent_count % codec->padded_side != 0)
            return false;
        plan->block_input_counts.assign(
            codec->parent_count / codec->padded_side, 0);
    }

    const auto observe_selected_coordinate = [&](uint32_t coordinate) {
        if (needs_generic && coordinate + 1 > plan->generic_input_count)
            plan->generic_input_count = coordinate + 1;
        if (needs_specialized)
        {
            const uint32_t block = coordinate / codec->padded_side;
            const uint32_t local_prefix = coordinate % codec->padded_side + 1;
            uint16_t& prefix = plan->block_input_counts[block];
            if (local_prefix > prefix)
                prefix = static_cast<uint16_t>(local_prefix);
            if (!selected_coordinates.empty())
                selected_coordinates[coordinate] = 1;
        }
    };

    // These are the only coordinates staged by execution.  Shortened zeros
    // are absent from the public loops, while punctured and surplus received
    // coordinates are excluded by coordinate_erased.
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        const uint32_t coordinate = CoordinateForOriginal(codec, i);
        if (plan->original_present[i] && !plan->coordinate_erased[coordinate])
            observe_selected_coordinate(coordinate);
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        const uint32_t coordinate = CoordinateForRecovery(codec, i);
        if (plan->recovery_present[i] && !plan->coordinate_erased[coordinate])
            observe_selected_coordinate(coordinate);
    }

    if (needs_generic)
    {
        plan->generic_output_dependencies.resize(
            leopard2_internal::OutputDependencyWordCount(codec->parent_count));
        if (!leopard2_internal::BuildOutputDependencies(
                codec->parent_count, plan->requested_coordinates.data(),
                plan->requested_coordinates.size(),
                plan->generic_output_dependencies.data(),
                plan->generic_output_dependencies.size()))
            return false;
    }

    if (needs_specialized && codec->profile == LEO2_PROFILE_LOW_V1)
    {
        plan->specialized_output_dependencies.resize(
            leopard2_internal::OutputDependencyWordCount(codec->padded_side));
        if (!leopard2_internal::BuildOutputDependencies(
                codec->padded_side, plan->requested_coordinates.data(),
                plan->requested_coordinates.size(),
                plan->specialized_output_dependencies.data(),
                plan->specialized_output_dependencies.size()))
            return false;
    }
    else if (needs_specialized)
    {
        size_t begin = 0;
        while (begin < plan->requested_coordinates.size())
        {
            const uint32_t coordinate = plan->requested_coordinates[begin];
            const uint32_t block = coordinate / codec->padded_side;
            if (block == 0)
                return false;
            size_t end = begin + 1;
            uint32_t requested_prefix = coordinate % codec->padded_side + 1;
            while (end < plan->requested_coordinates.size() &&
                   plan->requested_coordinates[end] / codec->padded_side == block)
            {
                requested_prefix =
                    plan->requested_coordinates[end] % codec->padded_side + 1;
                ++end;
            }
            leopard2_internal::DecodeOutputBlock descriptor = {
                block,
                requested_prefix,
                static_cast<uint32_t>(begin),
                static_cast<uint32_t>(end)
            };
            plan->high_output_blocks.push_back(descriptor);
            begin = end;
        }

        const uint32_t side = codec->padded_side;
        const uint32_t block_count = codec->parent_count / side;
        std::vector<uint8_t> input_mask(side, 0);
        std::vector<uint8_t> output_mask(side, 1);
        for (uint32_t block = 0; block < block_count; ++block)
        {
            const uint32_t offset = block * side;
            uint32_t live_count = 0;
            for (uint32_t i = 0; i < side; ++i)
            {
                input_mask[i] = selected_coordinates[offset + i];
                live_count += input_mask[i];
            }
            if (live_count == 0 || live_count == side)
                continue;
            leopard2_internal::PrunedTransformPlan candidate;
            if (PrepareCodecPrunedTransform(
                    codec, side, offset, true, input_mask.data(),
                    output_mask.data(), candidate) &&
                PrunedPlanSavesByteHeavyWork(candidate))
            {
                leopard2_internal::PrunedTransformBlock entry;
                entry.block = block;
                entry.plan = std::move(candidate);
                plan->high_pruned_input_blocks.push_back(std::move(entry));
            }
        }

        input_mask.assign(side, 1);
        plan->high_pruned_output_plans.reserve(
            plan->high_output_blocks.size());
        for (size_t block_index = 0;
             block_index < plan->high_output_blocks.size();
             ++block_index)
        {
            const leopard2_internal::DecodeOutputBlock& descriptor =
                plan->high_output_blocks[block_index];
            output_mask.assign(side, 0);
            const uint32_t offset = descriptor.block * side;
            for (uint32_t i = descriptor.requested_begin;
                 i < descriptor.requested_end;
                 ++i)
                output_mask[plan->requested_coordinates[i] - offset] = 1;
            leopard2_internal::PrunedTransformPlan candidate;
            if (!PrepareCodecPrunedTransform(
                    codec, side, offset, false, input_mask.data(),
                    output_mask.data(), candidate) ||
                !PrunedPlanSavesByteHeavyWork(candidate))
                candidate = leopard2_internal::PrunedTransformPlan();
            plan->high_pruned_output_plans.push_back(std::move(candidate));
        }
    }
    return true;
}

#ifdef LEO_HAS_FF8
struct DirectField8
{
    typedef leopard::ff8::ffe_t Element;
    static Element Multiply(Element a, Element b)
    {
        return leopard::ff8::MultiplyElements(a, b);
    }
    static Element Inverse(Element value)
    {
        return leopard::ff8::InverseElement(value);
    }
    static Element Log(Element value)
    {
        return leopard::ff8::ElementLog(value);
    }
    static void MultiplyBytes(
        const leopard::backend::Ops& ops,
        void* destination,
        const void* source,
        Element multiplier_log,
        uint64_t byte_count)
    {
        leopard::ff8::MultiplyBytes(
            ops, destination, source, multiplier_log, byte_count);
    }
    static void MultiplyAddBytes(
        const leopard::backend::Ops& ops,
        void* destination,
        const void* source,
        Element multiplier_log,
        uint64_t byte_count)
    {
        leopard::ff8::MultiplyAddBytes(
            ops, destination, source, multiplier_log, byte_count);
    }
};
#endif

#ifdef LEO_HAS_FF16
struct DirectField16
{
    typedef leopard::ff16::ffe_t Element;
    static Element Multiply(Element a, Element b)
    {
        return leopard::ff16::MultiplyElements(a, b);
    }
    static Element Inverse(Element value)
    {
        return leopard::ff16::InverseElement(value);
    }
    static Element Log(Element value)
    {
        return leopard::ff16::ElementLog(value);
    }
    static void MultiplyBytes(
        const leopard::backend::Ops& ops,
        void* destination,
        const void* source,
        Element multiplier_log,
        uint64_t byte_count)
    {
        leopard::ff16::MultiplyBytes(
            ops, destination, source, multiplier_log, byte_count);
    }
    static void MultiplyAddBytes(
        const leopard::backend::Ops& ops,
        void* destination,
        const void* source,
        Element multiplier_log,
        uint64_t byte_count)
    {
        leopard::ff16::MultiplyAddBytes(
            ops, destination, source, multiplier_log, byte_count);
    }
};
#endif

static bool IsDirectEncodeShape(const leo2_codec* codec)
{
    if (!codec || codec->original_count == 0 || codec->recovery_count == 0 ||
        codec->original_count > kDirectMaxOriginals ||
        codec->recovery_count > kDirectMaxRecoveries)
        return false;
    size_t coefficient_count = 0;
    return CheckedMultiply(
        static_cast<size_t>(codec->original_count),
        static_cast<size_t>(codec->recovery_count), coefficient_count) &&
        coefficient_count != 0;
}

static bool CanAutoDirectEncodeCodec(const leo2_codec* codec)
{
    if (!IsDirectEncodeShape(codec) || !codec->context ||
        codec->profile != LEO2_PROFILE_LOW_V1 || codec->original_count < 2)
        return false;
    const leo2_backend backend = codec->context->backend;
    if (backend == LEO2_BACKEND_SCALAR)
        return codec->original_count >= 3;
    return backend == LEO2_BACKEND_SSSE3 || backend == LEO2_BACKEND_AVX2;
}

static bool ShouldPrepareDirectEncode(const leo2_codec* codec)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    // Differential tests exercise both profiles and every bounded fringe.
    return IsDirectEncodeShape(codec);
#else
    return CanAutoDirectEncodeCodec(codec);
#endif
}

static bool CanPrepareDirectRepair(const leo2_codec* codec)
{
    return (codec->flags & (LEO2_CODEC_FORCE_GENERIC_DECODE |
                           LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
                           LEO2_CODEC_FORCE_TILED_DECODE |
                           LEO2_CODEC_FORCE_MATERIALIZED_DECODE)) == 0 &&
        codec->original_count >= 2 &&
        codec->original_count <= kDirectMaxOriginals &&
        codec->parent_dimension <= kDirectMaxParentDimension &&
        codec->padded_side >= 2;
}

template<class Field>
static bool PrepareDirectBarycentricWeights(
    const leo2_codec* codec,
    std::vector<typename Field::Element>& weights)
{
    typedef typename Field::Element Element;
    const uint32_t systematic_begin =
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? codec->padded_side : 0;
    const uint32_t systematic_end = systematic_begin + codec->parent_dimension;
    weights.resize(codec->original_count);
    for (uint32_t original = 0; original < codec->original_count; ++original)
    {
        const uint32_t coordinate = systematic_begin + original;
        Element denominator = 1;
        for (uint32_t other = systematic_begin; other < systematic_end; ++other)
        {
            if (other == coordinate)
                continue;
            denominator = Field::Multiply(
                denominator, static_cast<Element>(coordinate ^ other));
        }
        if (denominator == 0)
        {
            weights.clear();
            return false;
        }
        weights[original] = Field::Inverse(denominator);
    }
    return true;
}

template<class Field>
static bool PrepareDirectGeneratorRow(
    const leo2_codec* codec,
    uint32_t recovery_index,
    const std::vector<typename Field::Element>& barycentric_weights,
    typename Field::Element* row)
{
    typedef typename Field::Element Element;
    if (!codec || recovery_index >= codec->recovery_count || !row ||
        barycentric_weights.size() != codec->original_count)
        return false;

    const uint32_t systematic_begin =
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? codec->padded_side : 0;
    const uint32_t systematic_end = systematic_begin + codec->parent_dimension;
    const uint32_t parity_coordinate =
        CoordinateForRecovery(codec, recovery_index);

    // This is the systematic Lagrange row over the complete parent message
    // set.  Coordinates after the K public originals are shortened zeros, but
    // remain factors in both Z(parity) and each public derivative Z'(x_i).
    Element vanishing_value = 1;
    for (uint32_t coordinate = systematic_begin;
         coordinate < systematic_end;
         ++coordinate)
    {
        vanishing_value = Field::Multiply(vanishing_value,
            static_cast<Element>(parity_coordinate ^ coordinate));
    }
    if (vanishing_value == 0)
        return false;

    for (uint32_t original = 0; original < codec->original_count; ++original)
    {
        const Element difference = static_cast<Element>(
            parity_coordinate ^ (systematic_begin + original));
        if (difference == 0)
            return false;
        row[original] = Field::Multiply(
            Field::Multiply(vanishing_value, Field::Inverse(difference)),
            barycentric_weights[original]);
        if (row[original] == 0)
            return false;
    }
    return true;
}

template<class Field>
static bool PrepareDirectGeneratorLogs(
    const leo2_codec* codec,
    const std::vector<typename Field::Element>& barycentric_weights,
    std::vector<typename Field::Element>& generator_logs)
{
    typedef typename Field::Element Element;
    if (!IsDirectEncodeShape(codec))
        return false;
    size_t coefficient_count = 0;
    if (!CheckedMultiply(
            static_cast<size_t>(codec->original_count),
            static_cast<size_t>(codec->recovery_count), coefficient_count))
        return false;

    generator_logs.resize(coefficient_count);
    Element row[kDirectMaxOriginals];
    for (uint32_t recovery = 0; recovery < codec->recovery_count; ++recovery)
    {
        if (!PrepareDirectGeneratorRow<Field>(
                codec, recovery, barycentric_weights, row))
        {
            generator_logs.clear();
            return false;
        }
        const size_t row_offset =
            static_cast<size_t>(recovery) * codec->original_count;
        for (uint32_t original = 0;
             original < codec->original_count;
             ++original)
        {
            generator_logs[row_offset + original] = Field::Log(row[original]);
        }
    }
    return true;
}

static bool HasDirectGeneratorRows(const leo2_codec* codec)
{
    if (!IsDirectEncodeShape(codec))
        return false;
    size_t coefficient_count = 0;
    if (!CheckedMultiply(
            static_cast<size_t>(codec->original_count),
            static_cast<size_t>(codec->recovery_count), coefficient_count))
        return false;
    return codec->field == LEO2_FIELD_GF8
        ? codec->direct_generator_logs8.size() == coefficient_count
        : codec->field == LEO2_FIELD_GF16 &&
            codec->direct_generator_logs16.size() == coefficient_count;
}

static bool AutoDirectEncodePreferred(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count)
{
    if (!CanAutoDirectEncodeCodec(codec) || !HasDirectGeneratorRows(codec) ||
        requested_recovery_count == 0 ||
        requested_recovery_count > codec->recovery_count ||
        shard_bytes == 0 || shard_bytes > std::numeric_limits<size_t>::max())
        return false;

    if (requested_recovery_count != 1 ||
        shard_bytes < kDirectMinimumMeasuredBytes ||
        shard_bytes % kDirectSimdTileBytes != 0)
        return false;

    // Pinned ABBA measurements promote only the regular SIMD-tile region.
    // Scalar K=2 loses for large shards, and ragged GF8 tails have a sharp
    // sawtooth crossover, so both retain the transform encoder.  Unmeasured
    // backends also stay conservative until they have their own evidence.
    return true;
}

static bool ShouldUseDirectEncode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (codec->test_encode_mode == LEO2_TEST_ENCODE_FORCE_DIRECT)
        return requested_recovery_count != 0 && HasDirectGeneratorRows(codec);
    if (codec->test_encode_mode == LEO2_TEST_ENCODE_FORCE_TRANSFORM)
        return false;
#endif
    return AutoDirectEncodePreferred(
        codec, shard_bytes, requested_recovery_count);
}

template<class Field>
static bool InvertDirectRepairMatrix(
    typename Field::Element matrix[kDirectMaxLosses][kDirectMaxLosses * 2],
    uint32_t size)
{
    typedef typename Field::Element Element;
    for (uint32_t column = 0; column < size; ++column)
    {
        uint32_t pivot = column;
        while (pivot < size && matrix[pivot][column] == 0)
            ++pivot;
        if (pivot == size)
            return false;
        if (pivot != column)
            for (uint32_t j = 0; j < size * 2; ++j)
                std::swap(matrix[pivot][j], matrix[column][j]);

        const Element pivot_value = matrix[column][column];
        if (pivot_value != 1)
        {
            const Element inverse = Field::Inverse(pivot_value);
            for (uint32_t j = 0; j < size * 2; ++j)
                matrix[column][j] = Field::Multiply(matrix[column][j], inverse);
        }
        for (uint32_t row = 0; row < size; ++row)
        {
            if (row == column || matrix[row][column] == 0)
                continue;
            const Element factor = matrix[row][column];
            for (uint32_t j = 0; j < size * 2; ++j)
                matrix[row][j] ^= Field::Multiply(factor, matrix[column][j]);
        }
    }
    return true;
}

template<class Field>
static bool PrepareDirectRepairTerms(
    leo2_decode_plan* plan,
    const std::vector<typename Field::Element>& barycentric_weights)
{
    typedef typename Field::Element Element;
    const leo2_codec* codec = plan->codec;
    const uint32_t losses = plan->missing_original_count;
    if (losses == 0 || losses > kDirectMaxLosses ||
        barycentric_weights.size() != codec->original_count)
        return false;

    uint32_t selected_parities[kDirectMaxLosses] = {};
    uint32_t selected_count = 0;
    for (uint32_t parity = 0;
         parity < codec->recovery_count && selected_count < losses;
         ++parity)
    {
        if (plan->recovery_present[parity])
            selected_parities[selected_count++] = parity;
    }
    if (selected_count != losses)
        return false;

    std::vector<Element> generator_rows(
        static_cast<size_t>(losses) * codec->original_count);
    for (uint32_t equation = 0; equation < losses; ++equation)
    {
        Element* row = &generator_rows[
            static_cast<size_t>(equation) * codec->original_count];
        if (!PrepareDirectGeneratorRow<Field>(
                codec, selected_parities[equation], barycentric_weights, row))
            return false;
    }

    Element augmented[kDirectMaxLosses][kDirectMaxLosses * 2] = {};
    for (uint32_t equation = 0; equation < losses; ++equation)
    {
        for (uint32_t missing = 0; missing < losses; ++missing)
        {
            augmented[equation][missing] = generator_rows[
                static_cast<size_t>(equation) * codec->original_count +
                plan->missing_originals[missing]];
        }
        augmented[equation][losses + equation] = 1;
    }
    if (!InvertDirectRepairMatrix<Field>(augmented, losses))
        return false;

    plan->direct_term_offsets.clear();
    plan->direct_terms.clear();
    plan->direct_term_offsets.push_back(0);
    for (uint32_t output = 0; output < losses; ++output)
    {
        for (uint32_t equation = 0; equation < losses; ++equation)
        {
            const Element coefficient = augmented[output][losses + equation];
            if (coefficient == 0)
                continue;
            leo2_direct_repair_term term = {
                kDirectRecoveryTag | selected_parities[equation],
                static_cast<uint16_t>(Field::Log(coefficient))
            };
            plan->direct_terms.push_back(term);
        }
        for (uint32_t original = 0; original < codec->original_count; ++original)
        {
            if (!plan->original_present[original])
                continue;
            Element coefficient = 0;
            for (uint32_t equation = 0; equation < losses; ++equation)
            {
                coefficient ^= Field::Multiply(
                    augmented[output][losses + equation],
                    generator_rows[static_cast<size_t>(equation) *
                        codec->original_count + original]);
            }
            if (coefficient == 0)
                continue;
            leo2_direct_repair_term term = {
                original, static_cast<uint16_t>(Field::Log(coefficient))
            };
            plan->direct_terms.push_back(term);
        }

        const size_t begin = plan->direct_term_offsets.back();
        const size_t end = plan->direct_terms.size();
        if (begin == end)
            return false;
        for (size_t i = begin + 1; i < end; ++i)
        {
            if (plan->direct_terms[i].multiplier_log == 0)
            {
                std::swap(plan->direct_terms[begin], plan->direct_terms[i]);
                break;
            }
        }
        plan->direct_term_offsets.push_back(end);
    }
    return plan->direct_term_offsets.size() == losses + 1;
}

static leo2_backend RuntimeBackend()
{
    return leopard::backend::ExecutionBackend();
}

static leo2_result EnsureInitialized()
{
    static std::mutex mutex;
    static bool attempted = false;
    static int result = Leopard_CallInitialize;
    std::lock_guard<std::mutex> lock(mutex);
    if (!attempted)
    {
        result = leo_init();
        attempted = true;
    }
    if (result == Leopard_Success)
        return LEO2_SUCCESS;
    if (leopard::backend::StartupQualificationFailure() ==
            leopard::backend::QualificationOutOfMemory)
        return LEO2_OUT_OF_MEMORY;
    return LEO2_INTERNAL_ERROR;
}

static bool MakeRange(const void* pointer, uint64_t bytes, AddressRange& range)
{
    if (!pointer || bytes == 0 || bytes > std::numeric_limits<uintptr_t>::max())
        return false;
    const uintptr_t begin = reinterpret_cast<uintptr_t>(pointer);
    const uintptr_t length = static_cast<uintptr_t>(bytes);
    if (begin > std::numeric_limits<uintptr_t>::max() - length)
        return false;
    range.begin = begin;
    range.end = begin + length;
    return true;
}

static bool RangesOverlap(const AddressRange& a, const AddressRange& b)
{
    return a.begin < b.end && b.begin < a.end;
}

static bool RangeLess(const AddressRange& a, const AddressRange& b)
{
    if (a.begin != b.begin)
        return a.begin < b.begin;
    return a.end < b.end;
}

static size_t MergeRanges(AddressRange* ranges, size_t count)
{
    if (count == 0)
        return 0;
    std::sort(ranges, ranges + count, RangeLess);
    size_t merged = 1;
    for (size_t i = 1; i < count; ++i)
    {
        AddressRange& back = ranges[merged - 1];
        if (ranges[i].begin <= back.end)
        {
            if (ranges[i].end > back.end)
                back.end = ranges[i].end;
        }
        else
            ranges[merged++] = ranges[i];
    }
    return merged;
}

static leo2_result ValidateDisjointRanges(
    AddressRange* input_ranges,
    size_t input_count,
    AddressRange* output_ranges,
    size_t output_count)
{
    if (output_count == 0)
        return LEO2_SUCCESS;

    /*
        Input-input aliasing is permitted.  The overwhelmingly common repair
        and R=1 encode cases have one output, so sorting and merging every
        input range only to compare it with that single range is unnecessary
        setup in the byte-execution path.  Preserve the exact overlap contract
        with a linear scan; the general merge/sweep remains preferable once
        there are multiple outputs.
    */
    if (output_count == 1)
    {
        for (size_t input_i = 0; input_i < input_count; ++input_i)
            if (RangesOverlap(input_ranges[input_i], output_ranges[0]))
                return LEO2_OVERLAP;
        return LEO2_SUCCESS;
    }

    input_count = MergeRanges(input_ranges, input_count);
    std::sort(output_ranges, output_ranges + output_count, RangeLess);
    for (size_t i = 1; i < output_count; ++i)
        if (RangesOverlap(output_ranges[i - 1], output_ranges[i]))
            return LEO2_OVERLAP;

    size_t input_i = 0;
    for (size_t output_i = 0; output_i < output_count; ++output_i)
    {
        while (input_i < input_count && input_ranges[input_i].end <= output_ranges[output_i].begin)
            ++input_i;
        if (input_i < input_count && RangesOverlap(input_ranges[input_i], output_ranges[output_i]))
            return LEO2_OVERLAP;
    }
    return LEO2_SUCCESS;
}

static leo2_result CheckScratch(
    void* scratch,
    size_t scratch_bytes,
    const ScratchLayout& layout,
    AddressRange& scratch_range)
{
    if (scratch_bytes < layout.total_bytes)
        return LEO2_SCRATCH_TOO_SMALL;
    if (layout.total_bytes == 0)
        return LEO2_SUCCESS;
    if (!scratch)
        return LEO2_INVALID_ARGUMENT;
    if ((reinterpret_cast<uintptr_t>(scratch) & (kScratchAlignment - 1)) != 0)
        return LEO2_BAD_ALIGNMENT;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
}

static bool HasValidSystematicPad(
    const leo2_codec* codec,
    const void* shard,
    uint64_t shard_bytes)
{
    if (codec->shard_layout != LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
        return true;
    LEO_DEBUG_ASSERT(shard != NULL && shard_bytes >= 2 && (shard_bytes & 1u) == 0);
    return static_cast<const uint8_t*>(shard)[static_cast<size_t>(shard_bytes - 1)] == 0;
}

static leo2_result ValidateEncodeBuffers(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes)
{
    if (!original || !recovery)
        return LEO2_INVALID_ARGUMENT;
    AddressRange scratch_range;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;

    size_t output_count = 0;
    AddressRange single_output = { 0, 0 };
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (!recovery[i])
            continue;
        AddressRange range;
        if (!MakeRange(recovery[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range))
            return LEO2_OVERLAP;
        if (++output_count == 1)
            single_output = range;
    }

    if (output_count <= 1)
    {
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            AddressRange range;
            if (!MakeRange(original[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range) ||
                (output_count == 1 && RangesOverlap(range, single_output)))
                return LEO2_OVERLAP;
            if (!HasValidSystematicPad(codec, original[i], shard_bytes))
                return LEO2_INVALID_ARGUMENT;
        }
        return LEO2_SUCCESS;
    }

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        AddressRange range;
        if (!MakeRange(original[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range))
            return LEO2_OVERLAP;
        if (!HasValidSystematicPad(codec, original[i], shard_bytes))
            return LEO2_INVALID_ARGUMENT;
    }
    AddressRange* ranges = reinterpret_cast<AddressRange*>(scratch);
    size_t input_count = 0;
    output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        MakeRange(original[i], shard_bytes, ranges[input_count++]);
    AddressRange* outputs = ranges + codec->original_count;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (recovery[i])
            MakeRange(recovery[i], shard_bytes, outputs[output_count++]);
    return ValidateDisjointRanges(ranges, input_count, outputs, output_count);
}

static leo2_result EncodeLayout(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    EncodeScratchGeometry& geometry)
{
    geometry = EncodeScratchGeometry();
    if (!codec || !RoundShardBytes(shard_bytes, geometry.rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    // GF16 shards contain complete two-byte symbols.  A partial final ALTMAP
    // tile is compactly scattered below; an unpaired byte has no GF16 symbol.
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
        return LEO2_UNSUPPORTED;
    geometry.work_count = static_cast<size_t>(codec->padded_side) * 2;
    geometry.tail_bytes = static_cast<size_t>(shard_bytes) &
        (kScratchAlignment - 1);
    geometry.aligned_prefix_bytes =
        static_cast<size_t>(shard_bytes) - geometry.tail_bytes;
    geometry.work_slot_bytes = geometry.tail_bytes == 0
        ? static_cast<size_t>(shard_bytes)
        : std::max(geometry.aligned_prefix_bytes, kScratchAlignment);
    /*
        Transform kernels operate on complete 64-byte tiles.  Aligned shards
        can therefore be read from and written to the caller's disjoint
        buffers directly.  A partial tile is executed separately: each source
        and low-profile parity staging slot is one fixed 64-byte tile, while
        transform work is reused after the aligned prefix.
    */
    const bool stage_partial_tile = geometry.tail_bytes != 0;
    const size_t original_slots = stage_partial_tile ? codec->original_count : 0;
    const size_t parity_pointer_count = codec->profile == LEO2_PROFILE_LOW_V1
        ? codec->recovery_count : 0;
    const size_t parity_slots = stage_partial_tile &&
        codec->profile == LEO2_PROFILE_LOW_V1 ? codec->recovery_count : 0;
    const size_t range_count = static_cast<size_t>(codec->original_count) + codec->recovery_count;
    const size_t pointer_count = static_cast<size_t>(codec->original_count) +
        geometry.work_count + parity_pointer_count;
    size_t staging_slot_count = 0;
    if (!CheckedAdd(original_slots, parity_slots, staging_slot_count) ||
        !ComputeSplitScratchLayout(
            range_count, pointer_count, staging_slot_count,
            geometry.work_count, geometry.work_slot_bytes,
            geometry.layout, geometry.work_data_offset))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static bool UseGenericDecode(
    const leo2_decode_plan* plan,
    size_t rounded_shard_bytes)
{
    const leo2_codec* codec = plan->codec;
    const bool force_generic =
        (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
    const bool force_specialized =
        (codec->flags & (LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
                         LEO2_CODEC_FORCE_TILED_DECODE |
                         LEO2_CODEC_FORCE_MATERIALIZED_DECODE)) != 0;
    const bool measured_balanced_full_recovery =
        leopard2_internal::ShouldUseBalancedGenericDecode(
            codec->profile, codec->field, codec->original_count,
            codec->recovery_count, codec->padded_side, codec->parent_count,
            plan->missing_original_count, rounded_shard_bytes,
            codec->context->backend);
    return force_generic ||
        (!force_specialized && measured_balanced_full_recovery);
}

static bool UseMaterializedDecode(
    const leo2_codec* codec,
    const leo2_decode_plan* plan,
    size_t rounded_shard_bytes)
{
    if ((codec->flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0)
        return false;
    if ((codec->flags & LEO2_CODEC_FORCE_MATERIALIZED_DECODE) != 0)
        return true;
    const uint32_t missing_original_count = plan
        ? plan->missing_original_count
        : 1; // One-shot scratch must cover every AUTO materialized plan.
    return leopard2_internal::ShouldUseMaterializedHighDecode(
        codec->profile, codec->field, codec->original_count,
        codec->recovery_count, codec->padded_side, codec->parent_count,
        missing_original_count, rounded_shard_bytes,
        codec->context->backend);
}

static leo2_result DecodeLayout(
    const leo2_codec* codec,
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    DecodeScratchGeometry& geometry)
{
    geometry = DecodeScratchGeometry();
    if (!codec || !RoundShardBytes(shard_bytes, geometry.rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
        return LEO2_UNSUPPORTED;
    geometry.tail_bytes = static_cast<size_t>(shard_bytes) &
        (kScratchAlignment - 1);
    geometry.aligned_prefix_bytes =
        static_cast<size_t>(shard_bytes) - geometry.tail_bytes;
    geometry.work_slot_bytes = geometry.tail_bytes == 0
        ? static_cast<size_t>(shard_bytes)
        : std::max(geometry.aligned_prefix_bytes, kScratchAlignment);
    const size_t range_count = static_cast<size_t>(codec->original_count) * 2 + codec->recovery_count;
    const bool force_generic =
        (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
    const bool known_generic = plan
        ? UseGenericDecode(plan, geometry.rounded_bytes)
        : force_generic;
    const bool known_materialized = !known_generic &&
        UseMaterializedDecode(codec, plan, geometry.rounded_bytes);
    if (known_generic || known_materialized ||
        (codec->profile != LEO2_PROFILE_LOW_V1 &&
         codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1))
    {
        geometry.work_slot_count = codec->parent_count;
    }
    else
    {
        size_t tiled_slot_count = 0;
        if (!CheckedMultiply(
                static_cast<size_t>(codec->padded_side), 2,
                tiled_slot_count))
            return LEO2_INVALID_COUNTS;
        if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        {
            const size_t output_slots = plan
                ? plan->requested_coordinates.size()
                : static_cast<size_t>(codec->recovery_count);
            if (!CheckedAdd(
                    tiled_slot_count, output_slots, tiled_slot_count))
                return LEO2_INVALID_COUNTS;
        }
        /* Keep the retained materialized specialized kernel when its regular
           N-slot workspace is no larger than the tiled form. */
        const bool force_tiled =
            (codec->flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0;
        geometry.work_slot_count = force_tiled
            ? tiled_slot_count
            : std::min(
                static_cast<size_t>(codec->parent_count), tiled_slot_count);
    }
    size_t pointer_count = 0;
    if (!CheckedAdd(
            static_cast<size_t>(codec->parent_count),
            geometry.work_slot_count,
            pointer_count))
        return LEO2_INVALID_COUNTS;
    /*
        Complete 64-byte tiles already use the transform kernels' native
        layout.  The selected coordinate inputs are immutable, so execution
        can point at the caller shards directly.  Specialized execution uses
        two side-sized transform tiles, plus one retained shard per requested
        high-profile output.  Generic fallback retains N work slots.  A ragged
        final tile is executed separately: K+R public-coordinate staging slots
        are fixed at 64 bytes each, while work slots are sized for the larger
        of the aligned prefix and that one tail tile.
    */
    const size_t input_slot_count = geometry.tail_bytes != 0
        ? static_cast<size_t>(codec->original_count) + codec->recovery_count
        : 0;
    if (!ComputeSplitScratchLayout(
            range_count, pointer_count, input_slot_count,
            geometry.work_slot_count, geometry.work_slot_bytes,
            geometry.layout, geometry.work_data_offset))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static leo2_result DirectDecodeLayout(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    ScratchLayout& layout,
    size_t& rounded_bytes)
{
    if (!codec || !RoundShardBytes(shard_bytes, rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
        return LEO2_UNSUPPORTED;
    // Direct execution writes straight to disjoint restored buffers.  Scratch
    // is needed only for overlap/range validation, never for shard data.
    const size_t range_count = static_cast<size_t>(codec->original_count) * 2 +
        codec->recovery_count;
    if (!ComputeScratchLayout(range_count, 0, 0, rounded_bytes, layout))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static bool UseSingleSideEncodeLayout(const leo2_codec* codec)
{
    if (!codec || codec->padded_side != 1)
        return false;
#ifdef LEO2_ENABLE_TEST_HOOKS
    /* Keep the test-only direct-versus-transform experiment honest.  Both
       production profiles reduce to copy/XOR when their padded side is one,
       but a forced mode must retain the instrumented transform plumbing that
       the diagnostic requested. */
    if (codec->test_encode_mode != LEO2_TEST_ENCODE_AUTO)
        return false;
#endif
    return true;
}

static leo2_result DirectEncodeLayout(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    ScratchLayout& layout,
    size_t& rounded_bytes)
{
    if (!codec || !RoundShardBytes(shard_bytes, rounded_bytes))
        return LEO2_INVALID_ARGUMENT;
    if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
        return LEO2_UNSUPPORTED;
    const size_t range_count =
        static_cast<size_t>(codec->original_count) + codec->recovery_count;
    if (!ComputeScratchLayout(range_count, 0, 0, rounded_bytes, layout))
        return LEO2_INVALID_COUNTS;
    return LEO2_SUCCESS;
}

static void CopyAndPad(void* destination, const void* source, size_t bytes, size_t rounded)
{
    memcpy(destination, source, bytes);
    if (rounded > bytes)
        memset(static_cast<uint8_t*>(destination) + bytes, 0, rounded - bytes);
}

static void ScatterGF16CompactTail(
    void* destination,
    const void* source,
    size_t bytes,
    size_t rounded)
{
    LEO_DEBUG_ASSERT(bytes != 0 && (bytes & 1u) == 0 && rounded >= bytes);
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    if (complete != 0)
        memcpy(output, input, complete);
    const size_t residual = bytes - complete;
    if (residual == 0)
        return;

    LEO_DEBUG_ASSERT(rounded == complete + 64);
    memset(output + complete, 0, rounded - complete);
    const size_t symbols = residual / 2;
    memcpy(output + complete, input + complete, symbols);
    memcpy(output + complete + 32, input + complete + symbols, symbols);
}

static void GatherGF16CompactTail(
    void* destination,
    const void* source,
    size_t bytes)
{
    LEO_DEBUG_ASSERT(bytes != 0 && (bytes & 1u) == 0);
    uint8_t* output = static_cast<uint8_t*>(destination);
    const uint8_t* input = static_cast<const uint8_t*>(source);
    const size_t complete = bytes & ~static_cast<size_t>(63u);
    if (complete != 0)
        memcpy(output, input, complete);
    const size_t residual = bytes - complete;
    if (residual == 0)
        return;

    const size_t symbols = residual / 2;
    memcpy(output + complete, input + complete, symbols);
    memcpy(output + complete + symbols, input + complete + 32, symbols);
}

static void StageShardForKernel(
    const leo2_codec* codec,
    void* destination,
    const void* source,
    size_t bytes,
    size_t rounded)
{
    if (codec->field == LEO2_FIELD_GF16 && rounded != bytes)
        ScatterGF16CompactTail(destination, source, bytes, rounded);
    else
        CopyAndPad(destination, source, bytes, rounded);
}

static void GatherShardFromKernel(
    const leo2_codec* codec,
    void* destination,
    const void* source,
    size_t bytes)
{
    if (codec->field == LEO2_FIELD_GF16 && (bytes & 63u) != 0)
        GatherGF16CompactTail(destination, source, bytes);
    else
        memcpy(destination, source, bytes);
}

static void ExecuteTransformDecodePass(
    const leo2_decode_plan* plan,
    size_t buffer_bytes,
    const void* const* coordinate_input,
    void** work,
    bool use_generic,
    bool use_tiled)
{
    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    const uint32_t* const requested_coordinates =
        plan->requested_coordinates.data();
    const unsigned requested_count =
        static_cast<unsigned>(plan->requested_coordinates.size());
    const leopard2_internal::PrunedTransformBlock* const high_input_plans =
        plan->high_pruned_input_blocks.empty()
            ? NULL : plan->high_pruned_input_blocks.data();
    const unsigned high_input_plan_count = static_cast<unsigned>(
        plan->high_pruned_input_blocks.size());
    const leopard2_internal::PrunedTransformPlan* const high_output_plans =
        plan->high_pruned_output_plans.empty()
            ? NULL : plan->high_pruned_output_plans.data();
    const unsigned high_output_plan_count = static_cast<unsigned>(
        plan->high_pruned_output_plans.size());
    void** const high_requested_output =
        work + static_cast<size_t>(codec->padded_side) * 2;

#ifdef LEO_HAS_FF8
    if (codec->field == LEO2_FIELD_GF8)
    {
        if (use_generic)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->parent_count,
                    plan->generic_output_dependencies.data(),
                    plan->generic_output_dependencies.size());
            leopard::ff8::ReedSolomonDecodePlanned(
                ops, buffer_bytes, codec->parent_count, coordinate_input,
                plan->generic_input_count, requested_coordinates,
                requested_count, dependencies, &plan->locator8[0], work);
        }
        else if (codec->profile == LEO2_PROFILE_LOW_V1 && use_tiled)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            leopard::ff8::ReedSolomonDecodeLowTiledPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, requested_count, dependencies,
                &plan->locator8[0], &codec->fixed_factors8[0], work);
        }
        else if (use_tiled)
            leopard::ff8::ReedSolomonDecodeHighTiledPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator8[0], &codec->fixed_factors8[0],
                high_requested_output, high_input_plans,
                high_input_plan_count, high_output_plans,
                high_output_plan_count, work);
        else if (codec->profile == LEO2_PROFILE_LOW_V1)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            leopard::ff8::ReedSolomonDecodeLowPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, requested_count, dependencies,
                &plan->locator8[0], &codec->fixed_factors8[0], work);
        }
        else
            leopard::ff8::ReedSolomonDecodeHighPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator8[0], &codec->fixed_factors8[0],
                high_input_plans, high_input_plan_count,
                high_output_plans, high_output_plan_count, work);
        return;
    }
#endif
#ifdef LEO_HAS_FF16
    {
        if (use_generic)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->parent_count,
                    plan->generic_output_dependencies.data(),
                    plan->generic_output_dependencies.size());
            leopard::ff16::ReedSolomonDecodePlanned(
                ops, buffer_bytes, codec->parent_count, coordinate_input,
                plan->generic_input_count, requested_coordinates,
                requested_count, dependencies, &plan->locator16[0], work);
        }
        else if (codec->profile == LEO2_PROFILE_LOW_V1 && use_tiled)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            leopard::ff16::ReedSolomonDecodeLowTiledPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, requested_count, dependencies,
                &plan->locator16[0], &codec->fixed_factors16[0], work);
        }
        else if (use_tiled)
            leopard::ff16::ReedSolomonDecodeHighTiledPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator16[0], &codec->fixed_factors16[0],
                high_requested_output, high_input_plans,
                high_input_plan_count, high_output_plans,
                high_output_plan_count, work);
        else if (codec->profile == LEO2_PROFILE_LOW_V1)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            leopard::ff16::ReedSolomonDecodeLowPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, requested_count, dependencies,
                &plan->locator16[0], &codec->fixed_factors16[0], work);
        }
        else
            leopard::ff16::ReedSolomonDecodeHighPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator16[0], &codec->fixed_factors16[0],
                high_input_plans, high_input_plan_count,
                high_output_plans, high_output_plan_count, work);
    }
#endif
}

static void PopulateDecodeCoordinates(
    const leo2_decode_plan* plan,
    const void* const* original,
    const void* const* recovery,
    size_t source_offset,
    size_t pass_bytes,
    uint8_t* staging_slots,
    void** coordinate_data)
{
    const leo2_codec* codec = plan->codec;
    const bool stage_inputs = staging_slots != NULL;
    std::fill(
        coordinate_data, coordinate_data + codec->parent_count,
        static_cast<void*>(NULL));

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        const uint32_t coordinate = CoordinateForOriginal(codec, i);
        if (!original[i] || plan->coordinate_erased[coordinate])
            continue;
        const uint8_t* const source =
            static_cast<const uint8_t*>(original[i]) + source_offset;
        if (stage_inputs)
        {
            uint8_t* const slot = staging_slots +
                static_cast<size_t>(i) * kScratchAlignment;
            StageShardForKernel(
                codec, slot, source, pass_bytes, kScratchAlignment);
            coordinate_data[coordinate] = slot;
        }
        else
            coordinate_data[coordinate] = const_cast<uint8_t*>(source);
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        const uint32_t coordinate = CoordinateForRecovery(codec, i);
        if (!recovery[i] || plan->coordinate_erased[coordinate])
            continue;
        const uint8_t* const source =
            static_cast<const uint8_t*>(recovery[i]) + source_offset;
        if (stage_inputs)
        {
            uint8_t* const slot = staging_slots +
                (static_cast<size_t>(codec->original_count) + i) *
                    kScratchAlignment;
            StageShardForKernel(
                codec, slot, source, pass_bytes, kScratchAlignment);
            coordinate_data[coordinate] = slot;
        }
        else
            coordinate_data[coordinate] = const_cast<uint8_t*>(source);
    }
}

static void GatherTransformDecodePass(
    const leo2_decode_plan* plan,
    void* const* restored_original,
    size_t destination_offset,
    size_t pass_bytes,
    void** work,
    bool use_generic,
    bool use_tiled)
{
    const leo2_codec* codec = plan->codec;
    void** const high_requested_output =
        work + static_cast<size_t>(codec->padded_side) * 2;
    for (size_t i = 0; i < plan->missing_originals.size(); ++i)
    {
        const uint32_t original_index = plan->missing_originals[i];
        const void* const source = use_generic || !use_tiled
            ? work[CoordinateForOriginal(codec, original_index)]
            : codec->profile == LEO2_PROFILE_LOW_V1
                ? work[CoordinateForOriginal(codec, original_index)]
                : high_requested_output[i];
        uint8_t* const destination =
            static_cast<uint8_t*>(restored_original[original_index]) +
            destination_offset;
        GatherShardFromKernel(codec, destination, source, pass_bytes);
    }
}

static void PopulateEncodeInputs(
    const leo2_codec* codec,
    const void* const* original,
    size_t source_offset,
    size_t pass_bytes,
    uint8_t* staging_slots,
    void** pointers)
{
    const bool stage_inputs = staging_slots != NULL;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        const uint8_t* const source =
            static_cast<const uint8_t*>(original[i]) + source_offset;
        if (stage_inputs)
        {
            uint8_t* const slot = staging_slots +
                static_cast<size_t>(i) * kScratchAlignment;
            StageShardForKernel(
                codec, slot, source, pass_bytes, kScratchAlignment);
            pointers[i] = slot;
        }
        else
            pointers[i] = const_cast<uint8_t*>(source);
    }
}

static void ExecuteTransformEncodePass(
    const leo2_codec* codec,
    size_t buffer_bytes,
    uint32_t requested_recovery_prefix,
    const void* const* padded_original,
    void** parity,
    void** work)
{
    const leopard::backend::Ops& ops = *codec->context->ops;
    if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
    {
        if (codec->padded_side == 1)
        {
            uint8_t* const output = static_cast<uint8_t*>(work[0]);
            memcpy(output, padded_original[0], buffer_bytes);
            uint32_t i = 1;
            for (; i + 1 < codec->original_count; i += 2)
                leopard::xor_mem_2to1(
                    ops, output, padded_original[i], padded_original[i + 1],
                    buffer_bytes);
            if (i < codec->original_count)
                leopard::xor_mem(
                    ops, output, padded_original[i], buffer_bytes);
        }
        else if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
            leopard::ff8::ReedSolomonEncode(
                ops, buffer_bytes, codec->original_count,
                requested_recovery_prefix, codec->padded_side,
                padded_original, work);
#else
            return;
#endif
        else
#ifdef LEO_HAS_FF16
            leopard::ff16::ReedSolomonEncode(
                ops, buffer_bytes, codec->original_count,
                requested_recovery_prefix, codec->padded_side,
                padded_original, work);
#else
            return;
#endif
        return;
    }

    if (codec->padded_side == 1)
    {
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (parity[i])
                memcpy(parity[i], padded_original[0], buffer_bytes);
        return;
    }
    if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
        leopard::ff8::ReedSolomonEncodeLow(
            ops, buffer_bytes, codec->original_count,
            codec->recovery_count, codec->padded_side, padded_original,
            parity, work);
#else
        return;
#endif
    else
#ifdef LEO_HAS_FF16
        leopard::ff16::ReedSolomonEncodeLow(
            ops, buffer_bytes, codec->original_count,
            codec->recovery_count, codec->padded_side, padded_original,
            parity, work);
#else
        return;
#endif
}

static void ExecuteSingleSideEncode(
    const leo2_codec* codec,
    size_t shard_bytes,
    const void* const* original,
    void* const* recovery)
{
    const leopard::backend::Ops& ops = *codec->context->ops;
    if (codec->profile == LEO2_PROFILE_LOW_V1)
    {
        LEO_DEBUG_ASSERT(codec->original_count == 1);
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (recovery[i])
                memcpy(recovery[i], original[0], shard_bytes);
        return;
    }

    LEO_DEBUG_ASSERT(codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1);
    LEO_DEBUG_ASSERT(codec->recovery_count == 1);
    if (!recovery[0])
        return;
    uint8_t* const output = static_cast<uint8_t*>(recovery[0]);
    memcpy(output, original[0], shard_bytes);
    uint32_t i = 1;
    for (; i + 1 < codec->original_count; i += 2)
        leopard::xor_mem_2to1(
            ops, output, original[i], original[i + 1], shard_bytes);
    if (i < codec->original_count)
        leopard::xor_mem(ops, output, original[i], shard_bytes);
}

static leo2_result ValidateDecodeBuffers(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    size_t scratch_bytes)
{
    const leo2_codec* codec = plan->codec;
    if (!original || !recovery || !restored)
        return LEO2_INVALID_ARGUMENT;
    AddressRange scratch_range;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;

    /*
        A one-loss plan already identifies its only output.  Validate that
        output first, then compare each received shard against it while the
        shard is otherwise being validated.  This preserves the public alias
        contract without materializing and sorting K+R address ranges on the
        common direct-repair path.
    */
    const bool single_output = plan->missing_original_count == 1;
    AddressRange single_output_range = { 0, 0 };
    uint32_t single_output_index = 0;
    if (single_output)
    {
        if (plan->missing_originals.size() != 1)
            return LEO2_INTERNAL_ERROR;
        single_output_index = plan->missing_originals[0];
        if (single_output_index >= codec->original_count ||
            !restored[single_output_index] ||
            !MakeRange(
                restored[single_output_index], shard_bytes,
                single_output_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(single_output_range, scratch_range))
            return LEO2_OVERLAP;
    }

    size_t input_count = 0;
    size_t output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if ((original[i] != NULL) != (plan->original_present[i] != 0))
            return LEO2_INVALID_ARGUMENT;
        if (original[i])
        {
            AddressRange range;
            if (!MakeRange(original[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range) ||
                (single_output && RangesOverlap(range, single_output_range)))
                return LEO2_OVERLAP;
            if (!HasValidSystematicPad(codec, original[i], shard_bytes))
                return LEO2_INVALID_ARGUMENT;
            ++input_count;
        }
        else
        {
            if (!restored[i])
                return LEO2_INVALID_ARGUMENT;
            if (single_output && i != single_output_index)
                return LEO2_INTERNAL_ERROR;
            AddressRange range;
            if (!MakeRange(restored[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range))
                return LEO2_OVERLAP;
            ++output_count;
        }
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if ((recovery[i] != NULL) != (plan->recovery_present[i] != 0))
            return LEO2_INVALID_ARGUMENT;
        if (recovery[i])
        {
            AddressRange range;
            if (!MakeRange(recovery[i], shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(range, scratch_range) ||
                (single_output && RangesOverlap(range, single_output_range)))
                return LEO2_OVERLAP;
            ++input_count;
        }
    }

    if (single_output)
        return output_count == 1 ? LEO2_SUCCESS : LEO2_INTERNAL_ERROR;

    AddressRange* ranges = reinterpret_cast<AddressRange*>(scratch);
    input_count = 0;
    output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (original[i])
            MakeRange(original[i], shard_bytes, ranges[input_count++]);
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (recovery[i])
            MakeRange(recovery[i], shard_bytes, ranges[input_count++]);
    AddressRange* outputs = ranges + input_count;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (!original[i])
            MakeRange(restored[i], shard_bytes, outputs[output_count++]);
    return ValidateDisjointRanges(ranges, input_count, outputs, output_count);
}

static void XorArbitraryBytes(
    const leopard::backend::Ops& ops,
    void* destination,
    const void* source,
    size_t bytes)
{
    leopard::xor_mem(ops, destination, source, bytes);
}

template<class Field>
static leo2_result ExecuteDirectEncodeRows(
    const leo2_codec* codec,
    size_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    const std::vector<typename Field::Element>& generator_logs)
{
    typedef typename Field::Element Element;
    const leopard::backend::Ops& ops = *codec->context->ops;
    size_t expected_coefficients = 0;
    if (!CheckedMultiply(
            static_cast<size_t>(codec->original_count),
            static_cast<size_t>(codec->recovery_count), expected_coefficients) ||
        generator_logs.size() != expected_coefficients)
        return LEO2_INTERNAL_ERROR;

    for (uint32_t recovery_index = 0;
         recovery_index < codec->recovery_count;
         ++recovery_index)
    {
        void* output = recovery[recovery_index];
        if (!output)
            continue;
        const size_t row_offset =
            static_cast<size_t>(recovery_index) * codec->original_count;
        for (uint32_t original_index = 0;
             original_index < codec->original_count;
             ++original_index)
        {
            const void* source = original[original_index];
            const Element multiplier_log =
                generator_logs[row_offset + original_index];
            if (original_index == 0)
            {
                if (multiplier_log == 0)
                    memcpy(output, source, shard_bytes);
                else
                    Field::MultiplyBytes(
                        ops, output, source, multiplier_log, shard_bytes);
            }
            else if (multiplier_log == 0)
            {
                XorArbitraryBytes(ops, output, source, shard_bytes);
            }
            else
            {
                Field::MultiplyAddBytes(
                    ops, output, source, multiplier_log, shard_bytes);
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ExecuteDirectEncode(
    const leo2_codec* codec,
    size_t shard_bytes,
    const void* const* original,
    void* const* recovery)
{
    if (!HasDirectGeneratorRows(codec))
        return LEO2_INTERNAL_ERROR;
#ifdef LEO_HAS_FF8
    if (codec->field == LEO2_FIELD_GF8)
        return ExecuteDirectEncodeRows<DirectField8>(
            codec, shard_bytes, original, recovery,
            codec->direct_generator_logs8);
#endif
#ifdef LEO_HAS_FF16
    return ExecuteDirectEncodeRows<DirectField16>(
        codec, shard_bytes, original, recovery,
        codec->direct_generator_logs16);
#else
    return LEO2_INTERNAL_ERROR;
#endif
}

static leo2_result ExecuteDirectRepair(
    const leo2_decode_plan* plan,
    size_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original)
{
    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    for (uint32_t output_index = 0;
         output_index < plan->missing_original_count;
         ++output_index)
    {
        const size_t begin = plan->direct_term_offsets[output_index];
        const size_t end = plan->direct_term_offsets[output_index + 1];
        if (begin == end)
            return LEO2_INTERNAL_ERROR;
        void* output = restored_original[plan->missing_originals[output_index]];
        for (size_t term_index = begin; term_index < end; ++term_index)
        {
            const leo2_direct_repair_term& term = plan->direct_terms[term_index];
            const bool parity = (term.tagged_source & kDirectRecoveryTag) != 0;
            const uint32_t source_index = term.tagged_source & ~kDirectRecoveryTag;
            const void* source = parity ? recovery[source_index] : original[source_index];
            if (!source)
                return LEO2_INTERNAL_ERROR;

            if (term_index == begin)
            {
                if (term.multiplier_log == 0)
                    memcpy(output, source, shard_bytes);
                else if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
                    leopard::ff8::MultiplyBytes(ops, output, source,
                        static_cast<leopard::ff8::ffe_t>(term.multiplier_log), shard_bytes);
#else
                    return LEO2_INTERNAL_ERROR;
#endif
                else
#ifdef LEO_HAS_FF16
                    leopard::ff16::MultiplyBytes(ops, output, source,
                        static_cast<leopard::ff16::ffe_t>(term.multiplier_log), shard_bytes);
#else
                    return LEO2_INTERNAL_ERROR;
#endif
            }
            else if (term.multiplier_log == 0)
                XorArbitraryBytes(ops, output, source, shard_bytes);
            else if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
                leopard::ff8::MultiplyAddBytes(ops, output, source,
                    static_cast<leopard::ff8::ffe_t>(term.multiplier_log), shard_bytes);
#else
                return LEO2_INTERNAL_ERROR;
#endif
            else
#ifdef LEO_HAS_FF16
                leopard::ff16::MultiplyAddBytes(ops, output, source,
                    static_cast<leopard::ff16::ffe_t>(term.multiplier_log), shard_bytes);
#else
                return LEO2_INTERNAL_ERROR;
#endif
        }
    }
    return LEO2_SUCCESS;
}

struct EncodeBatchTaskContext
{
    const leo2_codec* codec;
    const leo2_encode_batch_item* items;
};

static leo2_result RunEncodeBatchItem(void* context, size_t index)
{
    const EncodeBatchTaskContext* batch = static_cast<const EncodeBatchTaskContext*>(context);
    const leo2_encode_batch_item& item = batch->items[index];
    return leo2_encode(
        batch->codec, item.shard_bytes, item.original, item.recovery,
        item.scratch, item.scratch_bytes);
}

struct DecodeBatchTaskContext
{
    const leo2_decode_plan* plan;
    const leo2_decode_batch_item* items;
    bool multi_item_batch;
};

static leo2_result RunDecodeBatchItem(void* context, size_t index)
{
    const DecodeBatchTaskContext* batch = static_cast<const DecodeBatchTaskContext*>(context);
    const leo2_decode_batch_item& item = batch->items[index];
    return DecodePlanExecuteInternal(
        batch->plan, item.shard_bytes, item.original, item.recovery,
        item.restored_original, item.scratch, item.scratch_bytes,
        batch->multi_item_batch);
}

} // namespace

extern "C" {

LEO2_EXPORT const char* leo2_result_string(int result)
{
    switch (result)
    {
    case LEO2_SUCCESS: return "Operation succeeded";
    case LEO2_NEED_MORE_DATA: return "Not enough received shards";
    case LEO2_INVALID_ARGUMENT: return "Invalid argument";
    case LEO2_INVALID_COUNTS: return "Invalid or unsupported shard counts";
    case LEO2_UNSUPPORTED: return "Requested profile, field, or backend is unsupported";
    case LEO2_SCRATCH_TOO_SMALL: return "Scratch buffer is too small";
    case LEO2_BAD_ALIGNMENT: return "Scratch buffer has insufficient alignment";
    case LEO2_OVERLAP: return "Unsupported input, output, or scratch overlap";
    case LEO2_OUT_OF_MEMORY: return "Allocation failed during setup";
    case LEO2_INTERNAL_ERROR: return "Internal initialization or execution error";
    }
    return "Unknown Leopard2 result";
}

LEO2_EXPORT uint32_t leo2_context_field_mask(const leo2_context* context)
{
    if (!context)
        return 0;
    uint32_t mask = 0;
#ifdef LEO_HAS_FF8
    mask |= LEO2_FIELD_MASK_GF8;
#endif
#ifdef LEO_HAS_FF16
    mask |= LEO2_FIELD_MASK_GF16;
#endif
    return mask;
}

LEO2_EXPORT leo2_result leo2_context_create(
    const leo2_context_options* options,
    leo2_context** context_out)
{
    if (!context_out)
        return LEO2_INVALID_ARGUMENT;
    *context_out = NULL;
    if (options && (options->struct_size < sizeof(leo2_context_options) ||
                    options->reserved != 0))
        return LEO2_INVALID_ARGUMENT;
    const uint32_t requested_raw = options
        ? options->backend : static_cast<uint32_t>(LEO2_BACKEND_AUTO);
    if (requested_raw > static_cast<uint32_t>(LEO2_BACKEND_NEON))
        return LEO2_INVALID_ARGUMENT;
    uint32_t threads = options ? options->thread_count : 0;
    if (threads > 128)
        return LEO2_INVALID_ARGUMENT;
    if (threads == 0)
        threads = DefaultThreadCount();

    const leo2_result initialized = EnsureInitialized();
    if (initialized != LEO2_SUCCESS)
        return initialized;

    const leo2_backend runtime_backend = RuntimeBackend();
    const leo2_backend requested = static_cast<leo2_backend>(requested_raw);
    const leopard::backend::Ops* ops = NULL;
    leo2_backend effective_backend = runtime_backend;
    if (requested == LEO2_BACKEND_AUTO)
        ops = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AUTO);
    else if (requested == LEO2_BACKEND_NEON)
    {
        // Existing ARM transform kernels are selected outside the fixed-ops
        // table.  Preserve exact NEON selection where it is already active;
        // forcing a lower table on that path requires the native-NEON refactor.
        if (runtime_backend != LEO2_BACKEND_NEON)
            return LEO2_UNSUPPORTED;
        ops = leopard::backend::GetQualifiedOps(LEO2_BACKEND_AUTO);
        effective_backend = LEO2_BACKEND_NEON;
    }
    else
    {
        if (runtime_backend == LEO2_BACKEND_NEON)
            return LEO2_UNSUPPORTED;
#if defined(LEO_TRY_SSSE3) || defined(LEO_TRY_AVX2)
        // A legacy whole-translation-unit ISA build can enter inline kernels
        // before the ops fallback.  AUTO remains compatible, but claiming an
        // exact explicit table would be false.  The ordinary CMake production
        // build isolates those ISAs and does not take this branch.
        return LEO2_UNSUPPORTED;
#else
        leopard::backend::QualificationStatus qualification =
            leopard::backend::QualificationAvailable;
        ops = leopard::backend::GetQualifiedOps(requested, &qualification);
        if (!ops)
        {
            if (qualification == leopard::backend::QualificationOutOfMemory)
                return LEO2_OUT_OF_MEMORY;
            if (qualification ==
                    leopard::backend::QualificationSelfTestFailed)
                return LEO2_INTERNAL_ERROR;
            return LEO2_UNSUPPORTED;
        }
        effective_backend = requested;
#endif
    }
    if (!ops || effective_backend == LEO2_BACKEND_AUTO)
        return LEO2_INTERNAL_ERROR;

    leo2_context* context = new (std::nothrow) leo2_context;
    if (!context)
        return LEO2_OUT_OF_MEMORY;
    context->backend = effective_backend;
    context->ops = ops;
    context->thread_count = threads;
    context->pool = NULL;
    if (threads > 1)
    {
        context->pool = new (std::nothrow) leo2_thread_pool(threads - 1);
        if (!context->pool)
        {
            delete context;
            return LEO2_OUT_OF_MEMORY;
        }
    }
    *context_out = context;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_context_destroy(leo2_context* context)
{
    if (context)
        delete context->pool;
    delete context;
}

LEO2_EXPORT leo2_backend leo2_context_backend(const leo2_context* context)
{
    return context ? context->backend : LEO2_BACKEND_AUTO;
}

LEO2_EXPORT uint32_t leo2_context_thread_count(const leo2_context* context)
{
    return context ? context->thread_count : 0;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
LEO2_EXPORT uint32_t leo2_test_context_worker_count(
    const leo2_context* context)
{
    return context && context->pool
        ? context->pool->TestWorkerCount()
        : 0;
}

LEO2_EXPORT void leo2_test_set_thread_start_fault(int enabled)
{
    g_test_thread_start_fault.store(enabled != 0, std::memory_order_release);
    g_test_thread_start_fault_consumptions.store(0, std::memory_order_release);
}

LEO2_EXPORT unsigned leo2_test_thread_start_fault_consumptions(void)
{
    return g_test_thread_start_fault_consumptions.load(
        std::memory_order_acquire);
}
#endif

LEO2_EXPORT leo2_result leo2_codec_create(
    leo2_context* context,
    uint32_t original_count,
    uint32_t recovery_count,
    leo2_profile profile,
    leo2_field field,
    const leo2_codec_options* options,
    leo2_codec** codec_out)
{
    if (!codec_out)
        return LEO2_INVALID_ARGUMENT;
    *codec_out = NULL;
    if (!context || original_count == 0 || recovery_count == 0)
        return LEO2_INVALID_ARGUMENT;
    leo2_shard_layout shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    if (options)
    {
        const uint32_t supported_flags =
            LEO2_CODEC_FORCE_GENERIC_DECODE |
            LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
            LEO2_CODEC_FORCE_TILED_DECODE |
            LEO2_CODEC_FORCE_MATERIALIZED_DECODE;
        const uint32_t algorithm_flags = options->flags &
            (LEO2_CODEC_FORCE_GENERIC_DECODE |
             LEO2_CODEC_FORCE_SPECIALIZED_DECODE);
        const uint32_t workspace_flags = options->flags &
            (LEO2_CODEC_FORCE_TILED_DECODE |
             LEO2_CODEC_FORCE_MATERIALIZED_DECODE);
        const size_t version1_size = offsetof(leo2_codec_options, shard_layout);
        const size_t layout_field_end = version1_size + sizeof(options->shard_layout);
        if (options->struct_size < version1_size ||
            (options->struct_size > version1_size &&
             options->struct_size < layout_field_end) ||
            options->reserved != 0 ||
            (options->flags & ~supported_flags) != 0 ||
            algorithm_flags == (LEO2_CODEC_FORCE_GENERIC_DECODE |
                                LEO2_CODEC_FORCE_SPECIALIZED_DECODE) ||
            workspace_flags == (LEO2_CODEC_FORCE_TILED_DECODE |
                                LEO2_CODEC_FORCE_MATERIALIZED_DECODE) ||
            ((algorithm_flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0 &&
             workspace_flags != 0))
            return LEO2_INVALID_ARGUMENT;
        if (options->struct_size >= layout_field_end)
        {
            const uint32_t raw_layout = options->shard_layout;
            if (raw_layout != LEO2_SHARD_LAYOUT_NATIVE_V1 &&
                raw_layout != LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
                return LEO2_INVALID_ARGUMENT;
            shard_layout = static_cast<leo2_shard_layout>(raw_layout);
        }
    }
    if (shard_layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1 &&
        field != LEO2_FIELD_GF16)
        return LEO2_INVALID_ARGUMENT;
    if (profile == LEO2_PROFILE_AUTO)
        profile = recovery_count <= original_count
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;
    if (profile != LEO2_PROFILE_LEGACY_HIGH_V1 && profile != LEO2_PROFILE_LOW_V1)
        return profile == LEO2_PROFILE_EXACT_EXPERIMENTAL_V1
            ? LEO2_UNSUPPORTED : LEO2_INVALID_ARGUMENT;

    const uint32_t padded = CeilPow2(
        profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? recovery_count : original_count);
    if (padded == 0)
        return LEO2_INVALID_COUNTS;
    const uint32_t parent = CeilPow2(static_cast<uint64_t>(padded) +
        (profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? original_count : recovery_count));
    if (parent == 0)
        return LEO2_INVALID_COUNTS;

    if (field == LEO2_FIELD_AUTO)
        field = parent <= kGF8Order ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    if (field == LEO2_FIELD_GF8 && parent > kGF8Order)
        return LEO2_INVALID_COUNTS;
    if (field == LEO2_FIELD_GF16 && parent > kGF16Order)
        return LEO2_INVALID_COUNTS;
    if (field != LEO2_FIELD_GF8 && field != LEO2_FIELD_GF16)
        return LEO2_INVALID_ARGUMENT;
#ifndef LEO_HAS_FF8
    if (field == LEO2_FIELD_GF8)
        return LEO2_UNSUPPORTED;
#endif
#ifndef LEO_HAS_FF16
    if (field == LEO2_FIELD_GF16)
        return LEO2_UNSUPPORTED;
#endif

    leo2_codec* codec = new (std::nothrow) leo2_codec;
    if (!codec)
        return LEO2_OUT_OF_MEMORY;
    codec->context = context;
    codec->original_count = original_count;
    codec->recovery_count = recovery_count;
    codec->parent_count = parent;
    codec->padded_side = padded;
    codec->parent_dimension = profile == LEO2_PROFILE_LEGACY_HIGH_V1
        ? parent - padded : padded;
    codec->profile = profile;
    codec->field = field;
    codec->shard_layout = shard_layout;
    codec->flags = options ? options->flags : 0;
#ifdef LEO2_ENABLE_TEST_HOOKS
    codec->test_encode_mode = LEO2_TEST_ENCODE_AUTO;
#endif
    try
    {
        codec->permanent_erased.assign(parent, 0);
        if (profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        {
            for (uint32_t i = recovery_count; i < padded; ++i)
                codec->permanent_erased[i] = 1;
        }
        else
        {
            for (uint32_t i = padded + recovery_count; i < parent; ++i)
                codec->permanent_erased[i] = 1;
        }

        const bool specialized =
            (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) == 0;
        const uint32_t permanent_erasure_count =
            profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? padded - recovery_count
                : parent - padded - recovery_count;

#ifdef LEO_HAS_FF8
        if (field == LEO2_FIELD_GF8)
        {
            const bool prepare_direct_repair = CanPrepareDirectRepair(codec);
            const bool prepare_direct_encode = ShouldPrepareDirectEncode(codec);
            if (prepare_direct_repair || prepare_direct_encode)
            {
                std::vector<DirectField8::Element> weights;
                if (PrepareDirectBarycentricWeights<DirectField8>(codec, weights))
                {
                    if (prepare_direct_encode)
                        PrepareDirectGeneratorLogs<DirectField8>(codec, weights,
                            codec->direct_generator_logs8);
                    if (prepare_direct_repair)
                        codec->direct_barycentric8.swap(weights);
                }
            }
            if (permanent_erasure_count != 0 &&
                leopard::ff8::IsDirectLocatorPreferred(parent, recovery_count))
            {
                codec->permanent_locator8.resize(parent);
                leopard::ff8::PrepareDecode(parent,
                    &codec->permanent_erased[0], &codec->permanent_locator8[0]);
            }
            if (specialized && padded >= 2)
            {
                if (profile == LEO2_PROFILE_LOW_V1)
                {
                    codec->fixed_factors8.resize(parent / padded - 1);
                    leopard::ff8::PrepareLowDecode(
                        parent, padded, &codec->fixed_factors8[0]);
                }
                else
                {
                    codec->fixed_factors8.resize(parent);
                    leopard::ff8::PrepareHighDecode(
                        parent, padded, &codec->fixed_factors8[0]);
                }
            }
        }
#endif
#ifdef LEO_HAS_FF16
        if (field == LEO2_FIELD_GF16)
        {
            const bool prepare_direct_repair = CanPrepareDirectRepair(codec);
            const bool prepare_direct_encode = ShouldPrepareDirectEncode(codec);
            if (prepare_direct_repair || prepare_direct_encode)
            {
                std::vector<DirectField16::Element> weights;
                if (PrepareDirectBarycentricWeights<DirectField16>(codec, weights))
                {
                    if (prepare_direct_encode)
                        PrepareDirectGeneratorLogs<DirectField16>(codec, weights,
                            codec->direct_generator_logs16);
                    if (prepare_direct_repair)
                        codec->direct_barycentric16.swap(weights);
                }
            }
            if (permanent_erasure_count != 0 &&
                leopard::ff16::IsDirectLocatorPreferred(parent, recovery_count))
            {
                codec->permanent_locator16.resize(parent);
                leopard::ff16::PrepareDecode(parent,
                    &codec->permanent_erased[0], &codec->permanent_locator16[0]);
            }
            if (specialized && padded >= 2)
            {
                if (profile == LEO2_PROFILE_LOW_V1)
                {
                    codec->fixed_factors16.resize(parent / padded - 1);
                    leopard::ff16::PrepareLowDecode(
                        parent, padded, &codec->fixed_factors16[0]);
                }
                else
                {
                    codec->fixed_factors16.resize(parent);
                    leopard::ff16::PrepareHighDecode(
                        parent, padded, &codec->fixed_factors16[0]);
                }
            }
        }
#endif
    }
    catch (const std::bad_alloc&)
    {
        delete codec;
        return LEO2_OUT_OF_MEMORY;
    }
    *codec_out = codec;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_codec_destroy(leo2_codec* codec)
{
    delete codec;
}

LEO2_EXPORT uint32_t leo2_codec_original_count(const leo2_codec* codec)
{
    return codec ? codec->original_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_recovery_count(const leo2_codec* codec)
{
    return codec ? codec->recovery_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_parent_count(const leo2_codec* codec)
{
    return codec ? codec->parent_count : 0;
}

LEO2_EXPORT uint32_t leo2_codec_padded_side(const leo2_codec* codec)
{
    return codec ? codec->padded_side : 0;
}

LEO2_EXPORT leo2_profile leo2_codec_profile(const leo2_codec* codec)
{
    return codec ? codec->profile : LEO2_PROFILE_AUTO;
}

LEO2_EXPORT leo2_field leo2_codec_field(const leo2_codec* codec)
{
    return codec ? codec->field : LEO2_FIELD_AUTO;
}

LEO2_EXPORT leo2_shard_layout leo2_codec_shard_layout(const leo2_codec* codec)
{
    return codec ? codec->shard_layout : LEO2_SHARD_LAYOUT_NATIVE_V1;
}

#ifdef LEO2_ENABLE_TEST_HOOKS
LEO2_EXPORT leo2_result leo2_test_codec_set_encode_mode(
    leo2_codec* codec,
    leo2_test_encode_mode mode)
{
    if (!codec || (mode != LEO2_TEST_ENCODE_AUTO &&
                   mode != LEO2_TEST_ENCODE_FORCE_DIRECT &&
                   mode != LEO2_TEST_ENCODE_FORCE_TRANSFORM))
        return LEO2_INVALID_ARGUMENT;
    if (mode == LEO2_TEST_ENCODE_FORCE_DIRECT &&
        !HasDirectGeneratorRows(codec))
        return LEO2_UNSUPPORTED;
    codec->test_encode_mode = mode;
    return LEO2_SUCCESS;
}

LEO2_EXPORT int leo2_test_codec_direct_encode_capable(
    const leo2_codec* codec)
{
    return HasDirectGeneratorRows(codec) ? 1 : 0;
}

LEO2_EXPORT leo2_result leo2_test_codec_encode_path(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count,
    int* direct_out)
{
    if (!direct_out)
        return LEO2_INVALID_ARGUMENT;
    *direct_out = 0;
    if (!codec || requested_recovery_count > codec->recovery_count)
        return LEO2_INVALID_ARGUMENT;
    EncodeScratchGeometry geometry;
    const leo2_result result = EncodeLayout(
        codec, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *direct_out = ShouldUseDirectEncode(
        codec, shard_bytes, requested_recovery_count) ? 1 : 0;
    return LEO2_SUCCESS;
}
#endif

LEO2_EXPORT leo2_result leo2_codec_wire_shard_bytes(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    uint64_t* wire_shard_bytes_out)
{
    if (!wire_shard_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    *wire_shard_bytes_out = 0;
    return ResolveWireShardBytes(codec, payload_bytes, *wire_shard_bytes_out);
}

LEO2_EXPORT leo2_result leo2_pack_systematic_shard(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    const void* payload,
    void* wire_shard,
    uint64_t wire_shard_bytes)
{
    uint64_t expected_wire_bytes = 0;
    const leo2_result result = ResolveWireShardBytes(
        codec, payload_bytes, expected_wire_bytes);
    if (result != LEO2_SUCCESS)
        return result;
    if (!payload || !wire_shard || wire_shard_bytes != expected_wire_bytes)
        return LEO2_INVALID_ARGUMENT;

    memmove(wire_shard, payload, static_cast<size_t>(payload_bytes));
    if (codec->shard_layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
        static_cast<uint8_t*>(wire_shard)[static_cast<size_t>(payload_bytes)] = 0;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_unpack_systematic_shard(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    const void* wire_shard,
    uint64_t wire_shard_bytes,
    void* payload)
{
    uint64_t expected_wire_bytes = 0;
    const leo2_result result = ResolveWireShardBytes(
        codec, payload_bytes, expected_wire_bytes);
    if (result != LEO2_SUCCESS)
        return result;
    if (!wire_shard || !payload || wire_shard_bytes != expected_wire_bytes)
        return LEO2_INVALID_ARGUMENT;
    if (codec->shard_layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1 &&
        static_cast<const uint8_t*>(wire_shard)[static_cast<size_t>(payload_bytes)] != 0)
        return LEO2_INVALID_ARGUMENT;

    memmove(payload, wire_shard, static_cast<size_t>(payload_bytes));
    return LEO2_SUCCESS;
}

LEO2_EXPORT size_t leo2_scratch_alignment(void)
{
    return kScratchAlignment;
}

LEO2_EXPORT leo2_result leo2_encode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!scratch_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    if (UseSingleSideEncodeLayout(codec))
    {
        ScratchLayout direct_layout;
        size_t rounded_bytes = 0;
        const leo2_result result = DirectEncodeLayout(
            codec, shard_bytes, direct_layout, rounded_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        *scratch_bytes_out = direct_layout.total_bytes;
        return LEO2_SUCCESS;
    }
    EncodeScratchGeometry geometry;
    const leo2_result result = EncodeLayout(codec, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = geometry.layout.total_bytes;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_encode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes)
{
    if (UseSingleSideEncodeLayout(codec))
    {
        ScratchLayout direct_layout;
        size_t rounded_bytes = 0;
        leo2_result result = DirectEncodeLayout(
            codec, shard_bytes, direct_layout, rounded_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        AddressRange scratch_range;
        result = CheckScratch(
            scratch, scratch_bytes, direct_layout, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        result = ValidateEncodeBuffers(
            codec, shard_bytes, original, recovery, scratch, scratch_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        ExecuteSingleSideEncode(
            codec, static_cast<size_t>(shard_bytes), original, recovery);
        return LEO2_SUCCESS;
    }

    EncodeScratchGeometry geometry;
    leo2_result result = EncodeLayout(codec, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange scratch_range;
    result = CheckScratch(
        scratch, scratch_bytes, geometry.layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    result = ValidateEncodeBuffers(
        codec, shard_bytes, original, recovery, scratch, scratch_bytes);
    if (result != LEO2_SUCCESS)
        return result;

    uint32_t requested_recovery_count = 0;
    uint32_t requested_recovery_prefix = 0;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (!recovery[i])
            continue;
        ++requested_recovery_count;
        requested_recovery_prefix = i + 1;
    }
    if (requested_recovery_count == 0)
        return LEO2_SUCCESS;
    if (ShouldUseDirectEncode(
            codec, shard_bytes, requested_recovery_count))
    {
        return ExecuteDirectEncode(codec, static_cast<size_t>(shard_bytes),
            original, recovery);
    }

    uint8_t* const base = static_cast<uint8_t*>(scratch);
    void** const pointers = reinterpret_cast<void**>(
        base + geometry.layout.pointer_offset);
    void** const work = pointers + codec->original_count;
    void** const parity = codec->profile == LEO2_PROFILE_LOW_V1
        ? work + geometry.work_count
        : NULL;
    uint8_t* const staging_slots = geometry.tail_bytes != 0
        ? base + geometry.layout.data_offset
        : NULL;
    uint8_t* const parity_staging =
        staging_slots && codec->profile == LEO2_PROFILE_LOW_V1
            ? staging_slots +
                static_cast<size_t>(codec->original_count) *
                    kScratchAlignment
            : NULL;
    uint8_t* const work_storage = base + geometry.work_data_offset;

    if (geometry.aligned_prefix_bytes != 0)
    {
        PopulateEncodeInputs(
            codec, original, 0, geometry.aligned_prefix_bytes,
            NULL, pointers);
        for (size_t i = 0; i < geometry.work_count; ++i)
            work[i] = work_storage + i * geometry.work_slot_bytes;
        if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
        {
            for (uint32_t i = 0; i < requested_recovery_prefix; ++i)
                if (recovery[i])
                    work[i] = static_cast<uint8_t*>(recovery[i]);
        }
        else
        {
            for (uint32_t i = 0; i < codec->recovery_count; ++i)
                parity[i] = recovery[i];
        }
        const void* const* const padded_original =
            const_cast<const void* const*>(pointers);
        ExecuteTransformEncodePass(
            codec, geometry.aligned_prefix_bytes,
            requested_recovery_prefix, padded_original, parity, work);
    }

    if (geometry.tail_bytes != 0)
    {
        PopulateEncodeInputs(
            codec, original, geometry.aligned_prefix_bytes,
            geometry.tail_bytes, staging_slots, pointers);
        for (size_t i = 0; i < geometry.work_count; ++i)
            work[i] = work_storage + i * geometry.work_slot_bytes;
        if (codec->profile == LEO2_PROFILE_LOW_V1)
        {
            for (uint32_t i = 0; i < codec->recovery_count; ++i)
                parity[i] = recovery[i]
                    ? parity_staging +
                        static_cast<size_t>(i) * kScratchAlignment
                    : NULL;
        }
        const void* const* const padded_original =
            const_cast<const void* const*>(pointers);
        ExecuteTransformEncodePass(
            codec, kScratchAlignment, requested_recovery_prefix,
            padded_original, parity, work);

        for (uint32_t i = 0; i < codec->recovery_count; ++i)
        {
            if (!recovery[i])
                continue;
            const void* const source =
                codec->profile == LEO2_PROFILE_LOW_V1
                    ? parity[i]
                    : work[i];
            uint8_t* const destination =
                static_cast<uint8_t*>(recovery[i]) +
                geometry.aligned_prefix_bytes;
            GatherShardFromKernel(
                codec, destination, source, geometry.tail_bytes);
        }
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_encode_batch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count)
{
    if (item_count > 0xffffffffu)
        return LEO2_INVALID_ARGUMENT;
    if (item_count != 0 && !items)
        return LEO2_INVALID_ARGUMENT;
    if (!codec)
        return LEO2_INVALID_ARGUMENT;
    EncodeBatchTaskContext batch = { codec, items };
    if (codec->context->pool)
        return codec->context->pool->Run(item_count, RunEncodeBatchItem, &batch);
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = RunEncodeBatchItem(&batch, i);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_plan_create(
    const leo2_codec* codec,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    leo2_decode_plan** plan_out)
{
    if (!plan_out)
        return LEO2_INVALID_ARGUMENT;
    *plan_out = NULL;
    if (!codec || !original_present || !recovery_present)
        return LEO2_INVALID_ARGUMENT;
    uint32_t present_count = 0;
    uint32_t missing_original_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if (original_present[i] > 1)
            return LEO2_INVALID_ARGUMENT;
        if (original_present[i])
            ++present_count;
        else
            ++missing_original_count;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (recovery_present[i] > 1)
            return LEO2_INVALID_ARGUMENT;
        if (recovery_present[i])
            ++present_count;
    }
    if (present_count < codec->original_count)
        return LEO2_NEED_MORE_DATA;

    leo2_decode_plan* plan = new (std::nothrow) leo2_decode_plan;
    if (!plan)
        return LEO2_OUT_OF_MEMORY;
    try
    {
        plan->codec = codec;
        plan->original_present.assign(original_present, original_present + codec->original_count);
        plan->recovery_present.assign(recovery_present, recovery_present + codec->recovery_count);
        plan->coordinate_erased = codec->permanent_erased;
        plan->generic_input_count = 0;
        plan->direct_copy_recovery = std::numeric_limits<uint32_t>::max();
        plan->missing_original_count = missing_original_count;
        plan->no_op = missing_original_count == 0;
        plan->direct_xor = codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;
        plan->direct_copy = codec->profile == LEO2_PROFILE_LOW_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;
        plan->direct_repair = false;

        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (!original_present[i])
            {
                plan->missing_originals.push_back(i);
                const uint32_t coordinate = CoordinateForOriginal(codec, i);
                plan->coordinate_erased[coordinate] = 1;
                plan->requested_coordinates.push_back(coordinate);
            }
        }
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            if (!recovery_present[i])
                plan->coordinate_erased[CoordinateForRecovery(codec, i)] = 1;

        /*
            Specialized decoders use exactly parent redundancy erasures.  Keep
            every surviving systematic shard, then the lowest-index parity
            shards needed to reach K public survivors; mark surplus received
            parity as deterministic virtual erasures.
        */
        uint32_t public_survivors_needed = codec->original_count;
        for (uint32_t i = 0; i < codec->original_count; ++i)
            if (original_present[i])
                --public_survivors_needed;
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
        {
            if (!recovery_present[i])
                continue;
            if (public_survivors_needed != 0)
                --public_survivors_needed;
            else
                plan->coordinate_erased[CoordinateForRecovery(codec, i)] = 1;
        }
        if (public_survivors_needed != 0)
        {
            delete plan;
            return LEO2_NEED_MORE_DATA;
        }

        if (plan->direct_copy)
        {
            for (uint32_t i = 0; i < codec->recovery_count; ++i)
            {
                const uint32_t coordinate = CoordinateForRecovery(codec, i);
                if (recovery_present[i] && !plan->coordinate_erased[coordinate])
                {
                    plan->direct_copy_recovery = i;
                    break;
                }
            }
            if (plan->direct_copy_recovery ==
                std::numeric_limits<uint32_t>::max())
            {
                delete plan;
                return LEO2_INTERNAL_ERROR;
            }
        }

        if (!plan->no_op && !plan->direct_xor && !plan->direct_copy &&
            missing_original_count <= kDirectMaxLosses)
        {
#ifdef LEO_HAS_FF8
            if (codec->field == LEO2_FIELD_GF8 &&
                !codec->direct_barycentric8.empty())
            {
                plan->direct_repair = PrepareDirectRepairTerms<DirectField8>(
                    plan, codec->direct_barycentric8);
            }
#endif
#ifdef LEO_HAS_FF16
            if (codec->field == LEO2_FIELD_GF16 &&
                     !codec->direct_barycentric16.empty())
            {
                plan->direct_repair = PrepareDirectRepairTerms<DirectField16>(
                    plan, codec->direct_barycentric16);
            }
#endif
        }

        if (!plan->no_op && !plan->direct_xor && !plan->direct_copy &&
            !plan->direct_repair)
        {
            if (!PreparePlanExecutionMetadata(plan))
            {
                delete plan;
                return LEO2_INTERNAL_ERROR;
            }
#ifdef LEO_HAS_FF8
            if (codec->field == LEO2_FIELD_GF8)
            {
                plan->locator8.resize(codec->parent_count);
                if (codec->permanent_locator8.empty())
                    leopard::ff8::PrepareDecode(codec->parent_count,
                        &plan->coordinate_erased[0], &plan->locator8[0]);
                else
                    leopard::ff8::PrepareDecodeWithPermanent(
                        codec->parent_count, &plan->coordinate_erased[0],
                        &codec->permanent_erased[0], &codec->permanent_locator8[0],
                        &plan->locator8[0]);
            }
#endif
#ifdef LEO_HAS_FF16
            if (codec->field == LEO2_FIELD_GF16)
            {
                plan->locator16.resize(codec->parent_count);
                if (codec->permanent_locator16.empty())
                    leopard::ff16::PrepareDecode(codec->parent_count,
                        &plan->coordinate_erased[0], &plan->locator16[0]);
                else
                    leopard::ff16::PrepareDecodeWithPermanent(
                        codec->parent_count, &plan->coordinate_erased[0],
                        &codec->permanent_erased[0], &codec->permanent_locator16[0],
                        &plan->locator16[0]);
            }
#endif
        }
    }
    catch (const std::bad_alloc&)
    {
        delete plan;
        return LEO2_OUT_OF_MEMORY;
    }
    *plan_out = plan;
    return LEO2_SUCCESS;
}

LEO2_EXPORT void leo2_decode_plan_destroy(leo2_decode_plan* plan)
{
    delete plan;
}

LEO2_EXPORT uint32_t leo2_decode_plan_missing_original_count(
    const leo2_decode_plan* plan)
{
    return plan ? plan->missing_original_count : 0;
}

LEO2_EXPORT leo2_result leo2_decode_plan_scratch_size(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!scratch_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    if (!plan || shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
        return LEO2_SUCCESS;
    ScratchLayout direct_layout;
    size_t direct_rounded = 0;
    DecodeScratchGeometry geometry;
    const bool direct_execution =
        plan->direct_repair || plan->direct_xor || plan->direct_copy;
    const leo2_result result = direct_execution
        ? DirectDecodeLayout(
            plan->codec, shard_bytes, direct_layout, direct_rounded)
        : DecodeLayout(plan->codec, plan, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = direct_execution
        ? direct_layout.total_bytes
        : geometry.layout.total_bytes;
    return LEO2_SUCCESS;
}

} // extern "C"

static leo2_result DecodePlanExecuteInternal(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes,
    bool multi_item_batch)
{
    if (!plan || shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
        return LEO2_SUCCESS;

    ScratchLayout direct_layout;
    size_t direct_rounded = 0;
    DecodeScratchGeometry geometry;
    const bool direct_execution =
        plan->direct_repair || plan->direct_xor || plan->direct_copy;
    leo2_result result = direct_execution
        ? DirectDecodeLayout(
            plan->codec, shard_bytes, direct_layout, direct_rounded)
        : DecodeLayout(plan->codec, plan, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    const ScratchLayout& layout = direct_execution
        ? direct_layout
        : geometry.layout;
    AddressRange scratch_range;
    result = CheckScratch(scratch, scratch_bytes, layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    result = ValidateDecodeBuffers(
        plan, shard_bytes, original, recovery, restored_original,
        scratch, scratch_bytes);
    if (result != LEO2_SUCCESS)
        return result;

    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    if (plan->direct_copy)
    {
        const uint32_t recovery_index = plan->direct_copy_recovery;
        if (recovery_index >= codec->recovery_count || !recovery[recovery_index])
            return LEO2_INTERNAL_ERROR;
        memcpy(restored_original[plan->missing_originals[0]],
            recovery[recovery_index], static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }
    if (plan->direct_xor)
    {
        const uint32_t missing = plan->missing_originals[0];
        if (missing >= codec->original_count || !recovery[0])
            return LEO2_INTERNAL_ERROR;
        uint8_t* output = static_cast<uint8_t*>(restored_original[missing]);
        memcpy(output, recovery[0], static_cast<size_t>(shard_bytes));
        const void* waiting = NULL;
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (!original[i])
                continue;
            if (!waiting)
                waiting = original[i];
            else
            {
                leopard::xor_mem_2to1(ops, output, waiting, original[i],
                    static_cast<size_t>(shard_bytes));
                waiting = NULL;
            }
        }
        if (waiting)
            XorArbitraryBytes(ops, output, waiting,
                static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }
    if (plan->direct_repair)
        return ExecuteDirectRepair(plan, static_cast<size_t>(shard_bytes),
            original, recovery, restored_original);

    uint8_t* const base = static_cast<uint8_t*>(scratch);
    void** const coordinate_data =
        reinterpret_cast<void**>(base + geometry.layout.pointer_offset);
    void** const work = coordinate_data + codec->parent_count;
    uint8_t* const staging_slots = geometry.tail_bytes != 0
        ? base + geometry.layout.data_offset
        : NULL;
    uint8_t* const work_storage = base + geometry.work_data_offset;
    for (size_t i = 0; i < geometry.work_slot_count; ++i)
        work[i] = work_storage + i * geometry.work_slot_bytes;

    /*
        At this balanced full-recovery point the generic schedule has lower
        aggregate work despite Algorithm 5's smaller transform side: the
        operation model counts 1280 versus 1856 butterflies.  Two reversed,
        CPU-pinned runs measured the generic decoder 5-32% faster from 256 B
        through 1 MiB on the three production x86 backends.  Keep dispatch
        strictly inside that measured region; neighboring counts, fields,
        backends, and sizes retain the profile-specific decoder.  A ragged
        final tile uses the same path as its aligned prefix.
    */
    const bool use_generic =
        UseGenericDecode(plan, geometry.rounded_bytes);
    const bool force_tiled =
        (codec->flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0;
    const bool force_materialized =
        (codec->flags & LEO2_CODEC_FORCE_MATERIALIZED_DECODE) != 0;
    const bool calibrated_materialized = !force_materialized &&
        UseMaterializedDecode(codec, plan, geometry.rounded_bytes);
    const bool measured_batch_tiled = multi_item_batch &&
        calibrated_materialized &&
        codec->context->backend == LEO2_BACKEND_AVX2;
    const bool use_tiled = !use_generic &&
        (force_tiled || measured_batch_tiled ||
         geometry.work_slot_count < static_cast<size_t>(codec->parent_count));

    if (geometry.aligned_prefix_bytes != 0)
    {
        PopulateDecodeCoordinates(
            plan, original, recovery, 0, geometry.aligned_prefix_bytes,
            NULL, coordinate_data);
        const void* const* const coordinate_input =
            const_cast<const void* const*>(coordinate_data);
        ExecuteTransformDecodePass(
            plan, geometry.aligned_prefix_bytes, coordinate_input, work,
            use_generic, use_tiled);
        GatherTransformDecodePass(
            plan, restored_original, 0, geometry.aligned_prefix_bytes,
            work, use_generic, use_tiled);
    }

    if (geometry.tail_bytes != 0)
    {
        PopulateDecodeCoordinates(
            plan, original, recovery, geometry.aligned_prefix_bytes,
            geometry.tail_bytes, staging_slots, coordinate_data);
        const void* const* const coordinate_input =
            const_cast<const void* const*>(coordinate_data);
        ExecuteTransformDecodePass(
            plan, kScratchAlignment, coordinate_input, work,
            use_generic, use_tiled);
        GatherTransformDecodePass(
            plan, restored_original, geometry.aligned_prefix_bytes,
            geometry.tail_bytes, work, use_generic, use_tiled);
    }
    return LEO2_SUCCESS;
}

extern "C" {

LEO2_EXPORT leo2_result leo2_decode_plan_execute(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes)
{
    return DecodePlanExecuteInternal(
        plan, shard_bytes, original, recovery, restored_original,
        scratch, scratch_bytes, false);
}

LEO2_EXPORT leo2_result leo2_decode_plan_execute_batch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count)
{
    if (item_count > 0xffffffffu)
        return LEO2_INVALID_ARGUMENT;
    if (item_count != 0 && !items)
        return LEO2_INVALID_ARGUMENT;
    if (!plan)
        return LEO2_INVALID_ARGUMENT;
    DecodeBatchTaskContext batch = { plan, items, item_count > 1 };
    if (plan->codec->context->pool)
        return plan->codec->context->pool->Run(item_count, RunDecodeBatchItem, &batch);
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = RunDecodeBatchItem(&batch, i);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_scratch_size(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    size_t* scratch_bytes_out)
{
    if (!scratch_bytes_out)
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    DecodeScratchGeometry geometry;
    const leo2_result result = DecodeLayout(
        codec, NULL, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *scratch_bytes_out = geometry.layout.total_bytes;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const uint8_t* original_present,
    const uint8_t* recovery_present,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes)
{
    leo2_decode_plan* plan = NULL;
    leo2_result result = leo2_decode_plan_create(
        codec, original_present, recovery_present, &plan);
    if (result != LEO2_SUCCESS)
        return result;
    result = leo2_decode_plan_execute(
        plan, shard_bytes, original, recovery, restored_original, scratch, scratch_bytes);
    leo2_decode_plan_destroy(plan);
    return result;
}

} // extern "C"
