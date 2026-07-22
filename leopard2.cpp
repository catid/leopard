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
#include <stdlib.h>
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

/*
    Build-time-only selector used by the bounded AVX2 loss-5-through-8
    experiment.  Zero retains the transform control, one selects the existing
    output-major executor, and two selects the source-major executor.  It is
    absent from the public ABI and remains default-off: exhaustive tiny-tail
    timing found valid byte counts where both direct executors regress.
*/
#ifndef LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE
#define LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE 0
#endif
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE < 0 || \
    LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE > 2
#error "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE must be 0, 1, or 2"
#endif

class leo2_thread_pool;

struct leo2_context
{
    leo2_backend backend;
    const leopard::backend::Ops* ops;
    const leopard::backend::Ops* baseline_ops;
    bool auto_requested;
    bool auto_avx512_encode_host;
    uint32_t thread_count;
    leo2_thread_pool* pool;
    // Lazily populated once per context for the measured GF8/AVX2 K=65
    // multi-loss repair family, then read immutably by every codec/plan.
    std::once_flag direct_repair_rows8_once;
    std::vector<uint8_t> direct_repair_generator_rows8;
};

#ifdef LEO2_ENABLE_TEST_HOOKS
namespace leopard { namespace backend {
void TestSetContextOps(leo2_context* context, const Ops* ops)
{
    if (context && ops)
        context->ops = ops;
}
}} // namespace leopard::backend

#ifndef __has_feature
#define __has_feature(feature) 0
#endif

extern "C" int leo2_test_core_address_sanitizer_compiled()
{
#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)
    return 1;
#else
    return 0;
#endif
}

extern "C" int leo2_test_core_undefined_sanitizer_compiled()
{
#if __has_feature(undefined_behavior_sanitizer)
    return 1;
#else
    return 0;
#endif
}

#if defined(_MSC_VER)
#define LEO2_TEST_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_TEST_NOINLINE __attribute__((noinline))
#else
#define LEO2_TEST_NOINLINE
#endif

extern "C" LEO2_TEST_NOINLINE size_t
leo2_test_core_leak_sanitizer_canary()
{
    // Keep the allocation in the linked codec translation unit: A replay-only
    // leak would prove the wrapper's sanitizer runtime, but not that the core
    // object being exercised participates in the leak check.  Volatile stores
    // make the allocation observable to the compiler without retaining a root
    // that would make LeakSanitizer consider the object reachable.
    static const size_t kCanaryBytes = 12345;
    volatile unsigned char* allocation = static_cast<volatile unsigned char*>(
        malloc(kCanaryBytes));
    if (!allocation)
        return 0;
    allocation[0] = 0x4c;
    allocation[kCanaryBytes - 1] = 0x32;
    return kCanaryBytes;
}

#undef LEO2_TEST_NOINLINE
#endif

struct leo2_codec
{
    leo2_context* context;
    const leopard::backend::Ops* auto_avx512_encode_ops;
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
    uint16_t direct_generator_numerator_log;
    bool direct_generator_numerator_valid;
    // Exact row-major generator coefficients, represented as Leopard field
    // logarithms, for the bounded allocation-free direct encoder.
    std::vector<uint8_t> direct_generator_logs8;
    std::vector<uint16_t> direct_generator_logs16;
    // In the P=T,N=2P legacy-high region, coordinate xor P turns the parent
    // into the Algorithm 4 low-profile view without changing public bytes.
    // Keep that translated view immutable and code-dependent in the codec.
    std::vector<uint8_t> translated_low_permanent_erased;
    std::vector<uint8_t> translated_low_permanent_locator8;
    std::vector<uint16_t> translated_low_permanent_locator16;
    std::vector<uint8_t> translated_low_factors8;
    std::vector<uint16_t> translated_low_factors16;
#ifdef LEO2_ENABLE_TEST_HOOKS
    leo2_test_encode_mode test_encode_mode;
    leo2_test_decode_mode test_decode_mode;
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
    // Algorithm 4 schedules are compiled during immutable plan setup.  Input
    // entries exist only for incomplete nonempty P-blocks that actually save
    // work.  An empty output plan selects the mature dependency-pruned FFT.
    std::vector<leopard2_internal::PrunedTransformBlock>
        low_pruned_input_blocks;
    leopard2_internal::PrunedTransformPlan low_pruned_output_plan;
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
    uint32_t direct_xor_original;
    uint32_t missing_original_count;
    bool no_op;
    bool direct_xor;
    bool direct_copy;
    bool direct_repair;
    /*
        A legacy-high P=T transform plan executes an Algorithm 4 view
        translated by xor P.  The plan captures this decision immutably;
        codec-owned low-profile constants and execution-order plan metadata
        remain private.  Direct repair and forced generic plans leave it false.
    */
    bool translated_low;
};

static bool PlanHasCompactDirectXor(const leo2_decode_plan* plan)
{
    return plan->direct_xor && plan->direct_xor_original !=
        std::numeric_limits<uint32_t>::max();
}

/*
    Compact R=1 plans avoid all pattern-owned vectors, but the mature vector
    representation still produces slightly better large-shard execution for
    K=2..7.  Keep that measured boundary while treating K=1 separately: its
    copy-only validator has no presence vector to scan and compact setup avoids
    three otherwise empty pattern allocations.  The one-shot API has its
    independent allocation-free path for every K.
*/
static const uint32_t kCompactDirectXorMinOriginalCount = 8;

static bool PlanOriginalPresent(
    const leo2_decode_plan* plan,
    uint32_t original)
{
    if (PlanHasCompactDirectXor(plan))
        return original != plan->direct_xor_original;
    return plan->original_present[original] != 0;
}

static bool PlanRecoveryPresent(
    const leo2_decode_plan* plan,
    uint32_t recovery)
{
    if (PlanHasCompactDirectXor(plan))
        return recovery == 0;
    return plan->recovery_present[recovery] != 0;
}

namespace leopard {
int InitializeLibrary(
    backend::QualificationStatus* qualification_failure_out);
}

struct ProtectedMetadataSpan
{
    const void* data;
    size_t bytes;
};

static leo2_result DecodePlanExecuteInternal(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original,
    void* scratch,
    size_t scratch_bytes,
    bool multi_item_batch,
    const ProtectedMetadataSpan* protected_metadata,
    size_t protected_metadata_count,
    bool prevalidated);

static leo2_result EncodeInternal(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes,
    const void* protected_metadata,
    size_t protected_metadata_bytes,
    bool prevalidated);

namespace {

static const size_t kScratchAlignment = 64;
static const uint32_t kGF8Order = 256;
static const uint32_t kGF16Order = 65536;
static const uint32_t kDirectRecoveryTag = 0x80000000u;
static const uint32_t kDirectEncodeMaxOriginals = 16;
static const uint32_t kDirectEncodeMaxRecoveries = 16;
static const uint32_t kDirectLegacyMaxRepairOriginals = 16;
static const uint32_t kDirectMaxRepairLosses = 8;
static const uint32_t kDirectMaxParentDimension = 256;
static const uint64_t kDirectSimdTileBytes = 64;
static const uint64_t kDirectMinimumMeasuredBytes = 1024;
static const size_t kScalableBatchPreflightMinItems = 9;
static const size_t kSparseEncodeScheduleBudget = 65536;

typedef leo2_result (*BatchTaskFunction)(void* context, size_t index);
typedef leo2_result (*BatchPreflightFunction)(void* context);

#ifdef LEO2_ENABLE_TEST_HOOKS
static std::atomic<bool> g_test_thread_start_fault(false);
static std::atomic<unsigned> g_test_thread_start_fault_consumptions(0);
static std::atomic<uint64_t> g_test_generic_direct_reveal_shards(0);
static std::atomic<uint64_t> g_test_low_direct_reveal_shards(0);
static std::atomic<uint64_t> g_test_low_scratch_reveal_shards(0);
static std::atomic<uint64_t> g_test_high_direct_reveal_shards(0);
static std::atomic<uint64_t>
    g_test_high_materialized_direct_reveal_shards(0);
static std::atomic<uint64_t> g_test_high_scratch_reveal_shards(0);

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
        void* function_context,
        BatchPreflightFunction preflight)
    {
        if (task_count == 0)
            return LEO2_SUCCESS;
        if (task_count == 1)
        {
            if (preflight)
            {
                const leo2_result result = preflight(function_context);
                if (result != LEO2_SUCCESS)
                    return result;
            }
            return function(function_context, 0);
        }

        /* One context supports concurrent codecs, but one batch owns its pool
           at a time.  Independent callers remain correct and are serialized at
           this optional scheduling layer rather than racing shared job state. */
        std::lock_guard<std::mutex> run_lock(run_mutex_);
        if (preflight)
        {
            const leo2_result result = preflight(function_context);
            if (result != LEO2_SUCCESS)
                return result;
        }
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

struct BatchRangeRecord
{
    AddressRange range;
    uintptr_t writable;
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
    leopard2_internal::DecodePathInfo selection;
    size_t rounded_bytes;
    size_t aligned_prefix_bytes;
    size_t tail_bytes;
    size_t execution_tile_count;
    size_t execution_tile_bytes;
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
    size_t execution_tile_count;
    size_t execution_tile_bytes;
    size_t work_count;
    size_t work_slot_bytes;
    size_t work_data_offset;
    size_t sparse_bits_offset;
    size_t sparse_workspace_offset;
    size_t sparse_block_count;
    size_t sparse_operation_stride;
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

static bool ComputeBalancedExecutionTiles(
    size_t aligned_bytes,
    size_t requested_tile_bytes,
    size_t& execution_tile_count,
    size_t& maximum_pass_bytes)
{
    execution_tile_count = aligned_bytes == 0 ? 0 : 1;
    maximum_pass_bytes = aligned_bytes;
    if (aligned_bytes == 0)
        return true;
    if (requested_tile_bytes == 0 ||
        (aligned_bytes & (kScratchAlignment - 1U)) != 0 ||
        (requested_tile_bytes & (kScratchAlignment - 1U)) != 0)
        return false;
    if (requested_tile_bytes >= aligned_bytes)
        return true;

    execution_tile_count = 1U +
        (aligned_bytes - 1U) / requested_tile_bytes;
    const size_t total_tiles = aligned_bytes / kScratchAlignment;
    const size_t maximum_pass_tiles =
        total_tiles / execution_tile_count +
        (total_tiles % execution_tile_count != 0);
    return maximum_pass_tiles != 0 &&
        CheckedMultiply(maximum_pass_tiles, kScratchAlignment,
            maximum_pass_bytes);
}

static leo2_result ComputeBatchPreflightScratchBytes(
    size_t item_count,
    size_t ranges_per_item,
    size_t& scratch_bytes)
{
    scratch_bytes = 0;
    if (item_count > 0xffffffffu)
        return LEO2_INVALID_ARGUMENT;
    if (item_count < kScalableBatchPreflightMinItems)
        return LEO2_SUCCESS;

    size_t range_count = 0;
    size_t range_bytes = 0;
    if (!CheckedMultiply(item_count, ranges_per_item, range_count) ||
        !CheckedAdd(range_count, 1, range_count) ||
        !CheckedMultiply(
            range_count, sizeof(BatchRangeRecord), range_bytes) ||
        !AlignUp(range_bytes, kScratchAlignment, scratch_bytes))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
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

static bool CanUseTranslatedLowDecode(const leo2_codec* codec)
{
    if (!codec || codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        codec->padded_side < 2 ||
        codec->parent_count != codec->padded_side * 2 ||
        codec->parent_dimension != codec->padded_side)
        return false;
    return CeilPow2(codec->original_count) == codec->padded_side;
}

static bool CodecMayCreateNativeHighPlan(const leo2_codec* codec)
{
    if (!codec || codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0)
        return false;
    if (!CanUseTranslatedLowDecode(codec))
        return true;
#ifdef LEO2_ENABLE_TEST_HOOKS
    // A hook codec may switch to FORCE_NATIVE_HIGH after construction.
    return true;
#else
    // Every non-direct transform plan captures translated_low=true.  The
    // translated selector then returns before generic/native-high dispatch.
    return false;
#endif
}

static bool CodecMayUseNativeLocator(const leo2_codec* codec)
{
    if (!codec)
        return false;
    if ((codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0)
        return true;
    if (!CanUseTranslatedLowDecode(codec))
        return true;
#ifdef LEO2_ENABLE_TEST_HOOKS
    // The mutable native-high hook still consumes the unshifted locator.
    return true;
#else
    return false;
#endif
}

static bool PlanUsesTranslatedLowDecode(const leo2_decode_plan* plan)
{
    return plan && plan->translated_low;
}

static uint32_t TranslateHalfParentCoordinate(
    const leo2_codec* codec,
    uint32_t coordinate)
{
    LEO_DEBUG_ASSERT(CanUseTranslatedLowDecode(codec));
    LEO_DEBUG_ASSERT(coordinate < codec->parent_count);
    return coordinate ^ codec->padded_side;
}

static const std::vector<uint32_t>& DecodeExecutionRequestedCoordinates(
    const leo2_decode_plan* plan)
{
    if (plan->translated_low)
    {
        // A high-profile systematic coordinate is P+i.  Translation by xor P
        // maps it directly to low-profile systematic coordinate i, which is
        // already stored in missing_originals in increasing order.
        return plan->missing_originals;
    }
    return plan->requested_coordinates;
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

static uint32_t DecodeExecutionCoordinateForOriginal(
    const leo2_decode_plan* plan,
    uint32_t index)
{
    const uint32_t coordinate = CoordinateForOriginal(plan->codec, index);
    return PlanUsesTranslatedLowDecode(plan)
        ? TranslateHalfParentCoordinate(plan->codec, coordinate)
        : coordinate;
}

static uint32_t DecodeExecutionCoordinateForRecovery(
    const leo2_decode_plan* plan,
    uint32_t index)
{
    const uint32_t coordinate = CoordinateForRecovery(plan->codec, index);
    return PlanUsesTranslatedLowDecode(plan)
        ? TranslateHalfParentCoordinate(plan->codec, coordinate)
        : coordinate;
}

static void TranslatePlanErasureMaskForExecution(leo2_decode_plan* plan)
{
    if (!PlanUsesTranslatedLowDecode(plan))
        return;
    const uint32_t side = plan->codec->padded_side;
    LEO_DEBUG_ASSERT(plan->coordinate_erased.size() == side * 2U);
    for (uint32_t i = 0; i < side; ++i)
        std::swap(plan->coordinate_erased[i], plan->coordinate_erased[i + side]);
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
    bool prepared = false;
#ifdef LEO_HAS_FF8
    if (codec->field == LEO2_FIELD_GF8)
        prepared = leopard::ff8::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, result);
#endif
#ifdef LEO_HAS_FF16
    if (codec->field == LEO2_FIELD_GF16)
        prepared = leopard::ff16::PreparePrunedTransformPlan(
            size, shift, inverse, input_mask, output_mask, result);
#endif
    if (!prepared)
        return false;

    // Pinned crossover measurements promote the exact inverse sink for AVX2
    // in both fields and for GF8 SSSE3.  GF16 SSSE3 did not clear the 5%
    // production threshold.  Do not retain unused root metadata in rejected
    // or unmeasured plans: besides wasting setup storage, it can evict the
    // materialized schedule those backends continue to execute.
    const leo2_backend backend = codec->context->ops->kind;
    const bool sink_backend = backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512 ||
        (codec->field == LEO2_FIELD_GF8 &&
         backend == LEO2_BACKEND_SSSE3);
    if (inverse && !sink_backend)
    {
        std::vector<uint8_t>().swap(result.inverse_accumulation_flags);
    }
    return true;
}

static bool PrunedPlanSavesByteHeavyWork(
    const leopard2_internal::PrunedTransformPlan& plan)
{
    return plan.size >= 2 &&
        (plan.operations.size() < plan.full_butterfly_count ||
         plan.input_zero_specializations != 0 ||
         plan.one_output_butterflies != 0);
}

static bool SkipMeasuredDenseGF8AVX2HighPrunedSchedules(
    const leo2_decode_plan* plan)
{
    if (!plan || !plan->codec)
        return false;
    const leo2_codec* codec = plan->codec;

    /*
        Pinned pure-AVX2 measurements found that compiling exact-mask input
        and output schedules for these dense maximum-loss shapes costs more
        than the mature regular transform saves.  Omitting both schedules
        reduced plan setup by 3.0x at T=8 and 5.8--7.9x at T=64, while the
        regular execution path was also 1--5% faster.  T=32 is deliberately
        excluded: its 64-KiB neighbors regressed by more than the production
        tolerance.  AUTO follows the backend selected for decode, while wider
        backends keep their independently measured policy.  This setup choice
        cannot alter the wire profile.
    */
    return codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->field == LEO2_FIELD_GF8 && codec->context &&
        codec->context->ops &&
        codec->context->ops->kind == LEO2_BACKEND_AVX2 &&
        (codec->padded_side == 8 || codec->padded_side == 64) &&
        codec->recovery_count == codec->padded_side &&
        plan->missing_original_count == codec->recovery_count &&
        static_cast<uint64_t>(codec->original_count) >=
            static_cast<uint64_t>(codec->recovery_count) * 2U;
}

static bool SkipMeasuredBoundaryGF8AVX2HighPrunedSchedules(
    const leo2_decode_plan* plan);

static bool PreparePlanExecutionMetadata(leo2_decode_plan* plan)
{
    const leo2_codec* codec = plan->codec;
    const bool translated_low = PlanUsesTranslatedLowDecode(plan);
    const bool force_generic =
        (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0;
    const bool force_specialized =
        (codec->flags & (LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
                         LEO2_CODEC_FORCE_TILED_DECODE |
                         LEO2_CODEC_FORCE_MATERIALIZED_DECODE)) != 0 ||
        translated_low;
    // AUTO owns both immutable paths because dispatch may depend on shard
    // bytes, backend, or a future offline table.  A forced diagnostic path
    // pays only for the metadata it can execute.
    const bool needs_generic = !force_specialized;
    const bool needs_specialized = !force_generic;
    const bool compile_pruned_schedules =
        !SkipMeasuredDenseGF8AVX2HighPrunedSchedules(plan) &&
        !SkipMeasuredBoundaryGF8AVX2HighPrunedSchedules(plan);

    const std::vector<uint32_t>& execution_requested_coordinates =
        DecodeExecutionRequestedCoordinates(plan);
    if (execution_requested_coordinates.empty())
        return false;
    for (size_t i = 1; i < execution_requested_coordinates.size(); ++i)
        if (execution_requested_coordinates[i - 1] >=
            execution_requested_coordinates[i])
            return false;

    plan->generic_input_count = 0;
    std::vector<uint8_t> selected_coordinates;
    if (needs_specialized)
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
            const uint32_t specialized_coordinate = translated_low
                ? TranslateHalfParentCoordinate(codec, coordinate)
                : coordinate;
            const uint32_t block = specialized_coordinate / codec->padded_side;
            const uint32_t local_prefix =
                specialized_coordinate % codec->padded_side + 1;
            uint16_t& prefix = plan->block_input_counts[block];
            if (local_prefix > prefix)
                prefix = static_cast<uint16_t>(local_prefix);
            if (!selected_coordinates.empty())
                selected_coordinates[specialized_coordinate] = 1;
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
                codec->parent_count, execution_requested_coordinates.data(),
                execution_requested_coordinates.size(),
                plan->generic_output_dependencies.data(),
                plan->generic_output_dependencies.size()))
            return false;
    }

    if (needs_specialized &&
        (codec->profile == LEO2_PROFILE_LOW_V1 || translated_low))
    {
        plan->specialized_output_dependencies.resize(
            leopard2_internal::OutputDependencyWordCount(codec->padded_side));
        if (!leopard2_internal::BuildOutputDependencies(
                codec->padded_side, execution_requested_coordinates.data(),
                execution_requested_coordinates.size(),
                plan->specialized_output_dependencies.data(),
                plan->specialized_output_dependencies.size()))
            return false;

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
                plan->low_pruned_input_blocks.push_back(std::move(entry));
            }
        }

        input_mask.assign(side, 1);
        output_mask.assign(side, 0);
        for (size_t i = 0; i < execution_requested_coordinates.size(); ++i)
        {
            const uint32_t coordinate = execution_requested_coordinates[i];
            if (coordinate >= side)
                return false;
            output_mask[coordinate] = 1;
        }
        // A complete requested message block is exactly the mature transform;
        // avoid expanding its full graph merely to discard the candidate.
        if (execution_requested_coordinates.size() < side)
        {
            leopard2_internal::PrunedTransformPlan candidate;
            if (PrepareCodecPrunedTransform(
                    codec, side, 0, false, input_mask.data(),
                    output_mask.data(), candidate) &&
                PrunedPlanSavesByteHeavyWork(candidate))
                plan->low_pruned_output_plan = std::move(candidate);
        }
    }
    else if (needs_specialized)
    {
        size_t begin = 0;
        while (begin < execution_requested_coordinates.size())
        {
            const uint32_t coordinate = execution_requested_coordinates[begin];
            const uint32_t block = coordinate / codec->padded_side;
            if (block == 0)
                return false;
            size_t end = begin + 1;
            uint32_t requested_prefix = coordinate % codec->padded_side + 1;
            while (end < execution_requested_coordinates.size() &&
                   execution_requested_coordinates[end] /
                       codec->padded_side == block)
            {
                requested_prefix =
                    execution_requested_coordinates[end] %
                        codec->padded_side + 1;
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
            if (compile_pruned_schedules && PrepareCodecPrunedTransform(
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
                output_mask[execution_requested_coordinates[i] - offset] = 1;
            leopard2_internal::PrunedTransformPlan candidate;
            if (!compile_pruned_schedules || !PrepareCodecPrunedTransform(
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
    static const uint32_t Modulus = 255;
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
    static Element FromLog(Element value_log)
    {
        return leopard::ff8::MultiplyLogElement(1, value_log);
    }
    static Element AddLogs(Element a, Element b)
    {
        uint32_t sum = static_cast<uint32_t>(a) + b;
        if (sum >= Modulus)
            sum -= Modulus;
        return static_cast<Element>(sum);
    }
    static Element SubtractLogs(Element a, Element b)
    {
        return static_cast<Element>(a >= b
            ? a - b
            : static_cast<uint32_t>(a) + Modulus - b);
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
    static const uint32_t Modulus = 65535;
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
    static Element FromLog(Element value_log)
    {
        return leopard::ff16::MultiplyLogElement(1, value_log);
    }
    static Element AddLogs(Element a, Element b)
    {
        uint32_t sum = static_cast<uint32_t>(a) + b;
        if (sum >= Modulus)
            sum -= Modulus;
        return static_cast<Element>(sum);
    }
    static Element SubtractLogs(Element a, Element b)
    {
        return static_cast<Element>(a >= b
            ? a - b
            : static_cast<uint32_t>(a) + Modulus - b);
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
        codec->original_count > kDirectEncodeMaxOriginals ||
        codec->recovery_count > kDirectEncodeMaxRecoveries)
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
    return backend == LEO2_BACKEND_SSSE3 ||
        backend == LEO2_BACKEND_AVX2 ||
        backend == LEO2_BACKEND_AVX512;
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

static bool IsMeasuredEqualRoundedDirectRepairCodec(const leo2_codec* codec)
{
    /*
        A deterministic 1,950-cell direct-versus-Algorithm-4 screen covered
        every requested K from 17 through 128, all three equal-rounded R
        boundaries, one through sixteen losses, and shard sizes from one byte
        through one MiB.  One-loss repair was the only broad new region that
        won every cell: its weakest execution result was 1.369x and its
        weakest reuse-one result including plan setup was 1.577x.  Keep the
        selector tied to the measured GF8/AVX2 P=T,N=2P construction.
    */
    return codec && codec->context &&
        CanUseTranslatedLowDecode(codec) &&
        codec->field == LEO2_FIELD_GF8 &&
        codec->context->backend == LEO2_BACKEND_AVX2 &&
        codec->original_count > kDirectLegacyMaxRepairOriginals &&
        codec->parent_count <= kGF8Order;
}

static bool IsExpandedDirectRepairCodec(const leo2_codec* codec)
{
    /* K=65 retains its separately measured eight-loss promotion. */
    return IsMeasuredEqualRoundedDirectRepairCodec(codec) &&
        codec->original_count == 65 && codec->padded_side == 128;
}

#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
static bool IsSmallDirectRepairCodec(const leo2_codec* codec)
{
    /*
        A same-source screen covered every K=5..16, R=5..8, and valid loss
        count from five through eight at 64, 65, and 2048 bytes.  Source-major
        direct repair won all 330 targets (weakest directional result 1.166x).
        Five-round head-to-head confirmation against output-major retained a
        1.118x lower confidence bound in its weakest target; the largest
        neighboring four-loss point regression was 1.351 percent.  A later
        1..63-byte screen found non-monotonic regressions in both direct loop
        orders, so retain this region only for explicit diagnostics until a
        byte-size-aware dual transform/direct plan is qualified.
    */
    return codec && codec->context &&
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->field == LEO2_FIELD_GF8 &&
        codec->context->backend == LEO2_BACKEND_AVX2 &&
        codec->original_count >= 5 &&
        codec->original_count <= kDirectLegacyMaxRepairOriginals &&
        codec->recovery_count >= 5 && codec->recovery_count <= 8;
}
#endif

static bool CanPrepareDirectRepair(const leo2_codec* codec)
{
    return (codec->flags & (LEO2_CODEC_FORCE_GENERIC_DECODE |
                           LEO2_CODEC_FORCE_SPECIALIZED_DECODE |
                           LEO2_CODEC_FORCE_TILED_DECODE |
                           LEO2_CODEC_FORCE_MATERIALIZED_DECODE)) == 0 &&
        codec->original_count >= 2 &&
        (codec->original_count <= kDirectLegacyMaxRepairOriginals ||
         IsMeasuredEqualRoundedDirectRepairCodec(codec)) &&
        codec->parent_dimension <= kDirectMaxParentDimension &&
        codec->padded_side >= 2;
}

static uint32_t DirectRepairLossLimit(const leo2_codec* codec)
{
    if (IsExpandedDirectRepairCodec(codec))
        return kDirectMaxRepairLosses;
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE != 0
    if (IsSmallDirectRepairCodec(codec))
        return kDirectMaxRepairLosses;
#endif
    if (IsMeasuredEqualRoundedDirectRepairCodec(codec))
        return 1;
    return 4;
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

    /*
        A complete aligned power-of-two coordinate interval is an additive
        coset in Leopard's Cantor index order.  For every x in that coset,
        xor-by-x maps the other coordinates onto the same nonzero subspace
        elements.  Its vanishing-polynomial derivative, and therefore every
        systematic barycentric denominator, is constant.  Equal-rounded high
        profiles and every low profile have this shape, so compute the product
        once instead of repeating the same O(D) loop for each public original.
    */
    const uint32_t dimension = codec->parent_dimension;
    if (dimension != 0 && (dimension & (dimension - 1)) == 0 &&
        (systematic_begin & (dimension - 1)) == 0)
    {
        Element denominator = 1;
        for (uint32_t difference = 1; difference < dimension; ++difference)
        {
            denominator = Field::Multiply(
                denominator, static_cast<Element>(difference));
        }
        if (denominator == 0)
        {
            weights.clear();
            return false;
        }
        std::fill(weights.begin(), weights.end(), Field::Inverse(denominator));
        return true;
    }

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
    const bool have_numerator_log =
        codec->direct_generator_numerator_valid;
    Element vanishing_log = have_numerator_log
        ? static_cast<Element>(codec->direct_generator_numerator_log)
        : 0;
    if (!have_numerator_log)
    {
        for (uint32_t coordinate = systematic_begin;
             coordinate < systematic_end;
             ++coordinate)
        {
            const Element difference = static_cast<Element>(
                parity_coordinate ^ coordinate);
            if (difference == 0)
                return false;
            vanishing_log =
                Field::AddLogs(vanishing_log, Field::Log(difference));
        }
    }

    for (uint32_t original = 0; original < codec->original_count; ++original)
    {
        const Element difference = static_cast<Element>(
            parity_coordinate ^ (systematic_begin + original));
        if (difference == 0)
            return false;
        const Element weight = barycentric_weights[original];
        if (weight == 0)
            return false;
        Element coefficient_log = Field::SubtractLogs(
            vanishing_log, Field::Log(difference));
        if (!have_numerator_log)
            coefficient_log =
                Field::AddLogs(coefficient_log, Field::Log(weight));
        row[original] = Field::FromLog(coefficient_log);
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
    Element row[kDirectEncodeMaxOriginals];
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

template<class Field>
static bool PrepareDirectRepairGeneratorRows(
    const leo2_codec* codec,
    const std::vector<typename Field::Element>& barycentric_weights,
    std::vector<typename Field::Element>& generator_rows)
{
    typedef typename Field::Element Element;
    size_t coefficient_count = 0;
    if (!codec || !CheckedMultiply(
            static_cast<size_t>(codec->original_count),
            static_cast<size_t>(codec->recovery_count), coefficient_count))
        return false;

    generator_rows.resize(coefficient_count);
    for (uint32_t recovery = 0;
         recovery < codec->recovery_count;
         ++recovery)
    {
        Element* row = generator_rows.data() +
            static_cast<size_t>(recovery) * codec->original_count;
        if (!PrepareDirectGeneratorRow<Field>(
                codec, recovery, barycentric_weights, row))
        {
            generator_rows.clear();
            return false;
        }
    }
    return true;
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
    typename Field::Element matrix
        [kDirectMaxRepairLosses][kDirectMaxRepairLosses * 2],
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
static bool PrepareDirectRepairCauchyInverse(
    const leo2_decode_plan* plan,
    const std::vector<typename Field::Element>& barycentric_weights,
    const std::vector<typename Field::Element>& cached_generator_rows,
    const uint32_t* selected_parities,
    typename Field::Element inverse
        [kDirectMaxRepairLosses][kDirectMaxRepairLosses])
{
    typedef typename Field::Element Element;
    const leo2_codec* codec = plan->codec;
    const uint32_t losses = plan->missing_original_count;
    if (!codec || losses == 0 || losses > kDirectMaxRepairLosses ||
        !selected_parities || barycentric_weights.size() !=
            codec->original_count)
        return false;

    size_t coefficient_count = 0;
    if (!CheckedMultiply(static_cast<size_t>(codec->original_count),
            static_cast<size_t>(codec->recovery_count), coefficient_count) ||
        cached_generator_rows.size() < coefficient_count)
        return false;

    Element parity_coordinates[kDirectMaxRepairLosses] = {};
    Element missing_coordinates[kDirectMaxRepairLosses] = {};
    Element row_factor_logs[kDirectMaxRepairLosses] = {};
    Element column_factor_logs[kDirectMaxRepairLosses] = {};
    Element parity_at_missing_logs[kDirectMaxRepairLosses] = {};
    for (uint32_t equation = 0; equation < losses; ++equation)
    {
        parity_coordinates[equation] = static_cast<Element>(
            CoordinateForRecovery(codec, selected_parities[equation]));
    }
    for (uint32_t missing = 0; missing < losses; ++missing)
    {
        missing_coordinates[missing] = static_cast<Element>(
            CoordinateForOriginal(codec,
                plan->missing_originals[missing]));
    }

    const Element reference_coordinate = static_cast<Element>(
        CoordinateForOriginal(codec, 0));
    const Element reference_weight_log =
        Field::Log(barycentric_weights[0]);
    for (uint32_t equation = 0; equation < losses; ++equation)
    {
        const Element* generator_row = cached_generator_rows.data() +
            static_cast<size_t>(selected_parities[equation]) *
                codec->original_count;
        const Element reference_difference = static_cast<Element>(
            parity_coordinates[equation] ^ reference_coordinate);
        if (reference_difference == 0 || generator_row[0] == 0)
            return false;
        Element row_scale_log = Field::AddLogs(
            Field::Log(generator_row[0]),
            Field::Log(reference_difference));
        row_scale_log =
            Field::SubtractLogs(row_scale_log, reference_weight_log);

        Element message_vanishing_log = 0;
        for (uint32_t missing = 0; missing < losses; ++missing)
        {
            const Element difference = static_cast<Element>(
                parity_coordinates[equation] ^
                missing_coordinates[missing]);
            if (difference == 0)
                return false;
            message_vanishing_log = Field::AddLogs(
                message_vanishing_log, Field::Log(difference));
        }
        Element parity_derivative_log = 0;
        for (uint32_t other = 0; other < losses; ++other)
        {
            if (other == equation)
                continue;
            const Element difference = static_cast<Element>(
                parity_coordinates[equation] ^
                parity_coordinates[other]);
            if (difference == 0)
                return false;
            parity_derivative_log = Field::AddLogs(
                parity_derivative_log, Field::Log(difference));
        }
        row_factor_logs[equation] = Field::SubtractLogs(
            Field::SubtractLogs(message_vanishing_log, row_scale_log),
            parity_derivative_log);
    }

    for (uint32_t missing = 0; missing < losses; ++missing)
    {
        Element parity_vanishing_log = 0;
        for (uint32_t equation = 0; equation < losses; ++equation)
        {
            const Element difference = static_cast<Element>(
                missing_coordinates[missing] ^
                parity_coordinates[equation]);
            if (difference == 0)
                return false;
            parity_vanishing_log = Field::AddLogs(
                parity_vanishing_log, Field::Log(difference));
        }
        parity_at_missing_logs[missing] = parity_vanishing_log;

        Element message_derivative_log = 0;
        for (uint32_t other = 0; other < losses; ++other)
        {
            if (other == missing)
                continue;
            const Element difference = static_cast<Element>(
                missing_coordinates[missing] ^
                missing_coordinates[other]);
            if (difference == 0)
                return false;
            message_derivative_log = Field::AddLogs(
                message_derivative_log, Field::Log(difference));
        }
        column_factor_logs[missing] = Field::SubtractLogs(0,
            Field::AddLogs(Field::Log(barycentric_weights[
                plan->missing_originals[missing]]),
                message_derivative_log));
    }

    for (uint32_t missing = 0; missing < losses; ++missing)
    {
        for (uint32_t equation = 0; equation < losses; ++equation)
        {
            const Element difference = static_cast<Element>(
                missing_coordinates[missing] ^
                parity_coordinates[equation]);
            Element inverse_log = Field::AddLogs(
                row_factor_logs[equation],
                column_factor_logs[missing]);
            inverse_log = Field::AddLogs(
                inverse_log, parity_at_missing_logs[missing]);
            inverse_log =
                Field::SubtractLogs(inverse_log, Field::Log(difference));
            inverse[missing][equation] = Field::FromLog(inverse_log);
        }
    }
    return true;
}

template<class Field>
static bool PrepareDirectRepairTerms(
    leo2_decode_plan* plan,
    const std::vector<typename Field::Element>& barycentric_weights,
    const std::vector<typename Field::Element>* cached_generator_rows)
{
    typedef typename Field::Element Element;
    const leo2_codec* codec = plan->codec;
    const uint32_t losses = plan->missing_original_count;
    if (losses == 0 || losses > DirectRepairLossLimit(codec) ||
        barycentric_weights.size() != codec->original_count)
        return false;

    uint32_t selected_parities[kDirectMaxRepairLosses] = {};
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

    size_t complete_row_count = 0;
    const bool have_cached_rows = CheckedMultiply(
        static_cast<size_t>(codec->original_count),
        static_cast<size_t>(codec->recovery_count), complete_row_count) &&
        cached_generator_rows &&
        cached_generator_rows->size() >= complete_row_count;
    std::vector<Element> generated_rows;
    if (!have_cached_rows)
    {
        generated_rows.resize(
            static_cast<size_t>(losses) * codec->original_count);
        for (uint32_t equation = 0; equation < losses; ++equation)
        {
            Element* row = generated_rows.data() +
                static_cast<size_t>(equation) * codec->original_count;
            if (!PrepareDirectGeneratorRow<Field>(codec,
                    selected_parities[equation], barycentric_weights, row))
                return false;
        }
    }

    const auto generator_row = [&](uint32_t equation) -> const Element* {
        if (have_cached_rows)
            return cached_generator_rows->data() +
                static_cast<size_t>(selected_parities[equation]) *
                    codec->original_count;
        return generated_rows.data() +
            static_cast<size_t>(equation) * codec->original_count;
    };

    Element augmented
        [kDirectMaxRepairLosses][kDirectMaxRepairLosses * 2] = {};
    for (uint32_t equation = 0; equation < losses; ++equation)
    {
        for (uint32_t missing = 0; missing < losses; ++missing)
        {
            augmented[equation][missing] =
                generator_row(equation)[plan->missing_originals[missing]];
        }
        augmented[equation][losses + equation] = 1;
    }
    bool inverted = false;
    if (have_cached_rows && codec->field == LEO2_FIELD_GF8)
    {
        Element cauchy_inverse
            [kDirectMaxRepairLosses][kDirectMaxRepairLosses] = {};
        if (PrepareDirectRepairCauchyInverse<Field>(plan,
                barycentric_weights, *cached_generator_rows,
                selected_parities, cauchy_inverse))
        {
            for (uint32_t output = 0; output < losses; ++output)
            {
                for (uint32_t equation = 0; equation < losses; ++equation)
                {
                    augmented[output][losses + equation] =
                        cauchy_inverse[output][equation];
                }
            }
            inverted = true;
        }
    }
    if (!inverted && !InvertDirectRepairMatrix<Field>(augmented, losses))
        return false;

    plan->direct_term_offsets.clear();
    plan->direct_terms.clear();
    plan->direct_term_offsets.reserve(static_cast<size_t>(losses) + 1U);
    size_t maximum_term_count = 0;
    if (!CheckedMultiply(static_cast<size_t>(losses),
            static_cast<size_t>(codec->original_count), maximum_term_count))
        return false;
    plan->direct_terms.reserve(maximum_term_count);
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

        alignas(kScratchAlignment)
            Element combined_row[kDirectMaxParentDimension] = {};
        bool have_combined_row = false;
        if (have_cached_rows && codec->field == LEO2_FIELD_GF8 &&
            codec->context && codec->context->ops)
        {
            const leopard::backend::Ops& ops = *codec->context->ops;
            for (uint32_t equation = 0; equation < losses; ++equation)
            {
                const Element coefficient =
                    augmented[output][losses + equation];
                if (coefficient == 0)
                    continue;
                const uint16_t multiplier_log =
                    static_cast<uint16_t>(Field::Log(coefficient));
                const Element* row = generator_row(equation);
                if (!have_combined_row)
                {
                    if (multiplier_log == 0)
                        memcpy(combined_row, row, codec->original_count);
                    else
                        ops.ff8_multiply(combined_row, row, multiplier_log,
                            codec->original_count);
                    have_combined_row = true;
                }
                else if (multiplier_log == 0)
                {
                    ops.xor_memory(combined_row, row, codec->original_count);
                }
                else
                {
                    ops.ff8_multiply_add(combined_row, row, multiplier_log,
                        codec->original_count);
                }
            }
        }

        for (uint32_t original = 0;
             original < codec->original_count;
             ++original)
        {
            if (!plan->original_present[original])
                continue;
            Element coefficient = combined_row[original];
            if (!have_combined_row)
            {
                for (uint32_t equation = 0; equation < losses; ++equation)
                {
                    coefficient ^= Field::Multiply(
                        augmented[output][losses + equation],
                        generator_row(equation)[original]);
                }
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
    // Every Leopard2 caller observes the same first process-wide setup result.
    // The detailed qualification status is captured by the shared initializer
    // while it still owns the legacy initialization lock; a concurrent direct
    // leo_init() retry therefore cannot race or rewrite this cached outcome.
    static std::mutex mutex;
    static bool attempted = false;
    static leo2_result cached_result = LEO2_INTERNAL_ERROR;
    std::lock_guard<std::mutex> lock(mutex);
    if (!attempted)
    {
        leopard::backend::QualificationStatus qualification =
            leopard::backend::QualificationAvailable;
        const int result = leopard::InitializeLibrary(&qualification);
        if (result == Leopard_Success)
            cached_result = LEO2_SUCCESS;
        else if (qualification ==
                    leopard::backend::QualificationOutOfMemory)
            cached_result = LEO2_OUT_OF_MEMORY;
        else
            cached_result = LEO2_INTERNAL_ERROR;
        attempted = true;
    }
    return cached_result;
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

static bool MakeArrayRange(
    const void* pointer,
    size_t count,
    size_t element_bytes,
    AddressRange& range)
{
    size_t bytes = 0;
    return pointer && count != 0 &&
        CheckedMultiply(count, element_bytes, bytes) &&
        MakeRange(pointer, static_cast<uint64_t>(bytes), range);
}

template <typename T>
static bool IsRepresentableOutput(T* output)
{
    AddressRange range;
    return MakeArrayRange(output, 1, sizeof(*output), range);
}

static bool RangesOverlap(const AddressRange& a, const AddressRange& b)
{
    return a.begin < b.end && b.begin < a.end;
}

static bool RangeOverlapsAny(
    const AddressRange& range,
    const AddressRange* others,
    size_t other_count)
{
    for (size_t i = 0; i < other_count; ++i)
        if (RangesOverlap(range, others[i]))
            return true;
    return false;
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

static bool RangesSortedByBegin(
    const AddressRange* ranges,
    size_t count)
{
    for (size_t i = 1; i < count; ++i)
        if (ranges[i].begin < ranges[i - 1].begin)
            return false;
    return true;
}

static leo2_result ValidateAlreadySortedDisjointRanges(
    const AddressRange* input_ranges,
    size_t input_count,
    const AddressRange* output_ranges,
    size_t output_count)
{
    /*
        The ordinary application layout stores shards in monotonically
        increasing slabs.  Validate that common case with one linear sweep
        instead of sorting K inputs and R outputs on every encode call.

        Read-only inputs may overlap or alias exactly, so retain the furthest
        end reached by every input whose begin precedes the current output's
        end.  If that frontier crosses the output begin, some input overlaps
        it.  Outputs remain pairwise disjoint and are checked independently.
        Callers reach this helper only after both lists prove monotonic by
        begin; arbitrary pointer orders retain the established sort/merge
        fallback below.
    */
    size_t input_i = 0;
    uintptr_t furthest_input_end = 0;
    bool have_input = false;
    for (size_t output_i = 0; output_i < output_count; ++output_i)
    {
        if (output_i != 0 &&
            RangesOverlap(output_ranges[output_i - 1],
                output_ranges[output_i]))
            return LEO2_OVERLAP;

        const AddressRange& output = output_ranges[output_i];
        while (input_i < input_count &&
            input_ranges[input_i].begin < output.end)
        {
            if (!have_input ||
                input_ranges[input_i].end > furthest_input_end)
                furthest_input_end = input_ranges[input_i].end;
            have_input = true;
            ++input_i;
        }
        if (have_input && furthest_input_end > output.begin)
            return LEO2_OVERLAP;
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateMonotonicEncodePointers(
    const void* const* original,
    size_t original_count,
    void* const* recovery,
    size_t recovery_count,
    uint64_t shard_bytes)
{
    /*
        Every pointer and end address has already passed MakeRange, and both
        caller arrays have already proved monotonic by range begin.  Sweep the
        caller arrays directly so the ordinary packed-shard encode path does
        not stage K+R AddressRange records in scratch and then scan them again.

        Inputs are immutable and may alias.  Preserve the furthest end of all
        inputs beginning before the current output ends; that catches a long
        or exactly aliased input even after the input cursor advances.  Null
        recovery entries are unrequested coordinates and do not participate.
    */
    const uintptr_t bytes = static_cast<uintptr_t>(shard_bytes);
    size_t input_i = 0;
    uintptr_t furthest_input_end = 0;
    uintptr_t furthest_output_end = 0;
    bool have_input = false;
    bool have_output = false;
    for (size_t output_i = 0; output_i < recovery_count; ++output_i)
    {
        if (!recovery[output_i])
            continue;
        const uintptr_t output_begin =
            reinterpret_cast<uintptr_t>(recovery[output_i]);
        const uintptr_t output_end = output_begin + bytes;
        if (have_output && furthest_output_end > output_begin)
            return LEO2_OVERLAP;
        if (!have_output || output_end > furthest_output_end)
            furthest_output_end = output_end;
        have_output = true;

        while (input_i < original_count)
        {
            const uintptr_t input_begin =
                reinterpret_cast<uintptr_t>(original[input_i]);
            if (input_begin >= output_end)
                break;
            const uintptr_t input_end = input_begin + bytes;
            if (!have_input || input_end > furthest_input_end)
                furthest_input_end = input_end;
            have_input = true;
            ++input_i;
        }
        if (have_input && furthest_input_end > output_begin)
            return LEO2_OVERLAP;
    }
    return LEO2_SUCCESS;
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

    if (RangesSortedByBegin(input_ranges, input_count) &&
        RangesSortedByBegin(output_ranges, output_count))
        return ValidateAlreadySortedDisjointRanges(
            input_ranges, input_count, output_ranges, output_count);

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

static leo2_result CheckOptionalScratch(
    void* scratch,
    size_t scratch_bytes,
    AddressRange& scratch_range)
{
    scratch_range.begin = 0;
    scratch_range.end = 0;
    if (scratch_bytes == 0)
        return LEO2_SUCCESS;
    if (!scratch)
        return LEO2_INVALID_ARGUMENT;
    if ((reinterpret_cast<uintptr_t>(scratch) &
            (kScratchAlignment - 1)) != 0)
        return LEO2_BAD_ALIGNMENT;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
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
        return CheckOptionalScratch(
            scratch, scratch_bytes, scratch_range);
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

static leo2_result ValidateK1R1EncodeBuffers(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    const AddressRange& scratch_range,
    const void* protected_metadata,
    size_t protected_metadata_bytes)
{
    LEO_DEBUG_ASSERT(codec && codec->original_count == 1 &&
        codec->recovery_count == 1);
    if (!original || !recovery)
        return LEO2_INVALID_ARGUMENT;

    AddressRange metadata_ranges[3];
    size_t metadata_count = 2;
    if (!MakeArrayRange(original, 1, sizeof(*original), metadata_ranges[0]) ||
        !MakeArrayRange(recovery, 1, sizeof(*recovery), metadata_ranges[1]))
        return LEO2_INVALID_ARGUMENT;
    if (protected_metadata)
    {
        if (protected_metadata_bytes == 0 ||
            !MakeRange(protected_metadata,
                static_cast<uint64_t>(protected_metadata_bytes),
                metadata_ranges[metadata_count]))
            return LEO2_INVALID_ARGUMENT;
        ++metadata_count;
    }
    if (RangeOverlapsAny(scratch_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

    AddressRange output_range = { 0, 0 };
    if (recovery[0])
    {
        if (!MakeRange(recovery[0], shard_bytes, output_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(output_range, scratch_range) ||
            RangeOverlapsAny(output_range, metadata_ranges, metadata_count))
            return LEO2_OVERLAP;
    }

    AddressRange input_range;
    if (!MakeRange(original[0], shard_bytes, input_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangesOverlap(input_range, scratch_range) ||
        (recovery[0] && RangesOverlap(input_range, output_range)))
        return LEO2_OVERLAP;
    if (!HasValidSystematicPad(codec, original[0], shard_bytes))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
}

static leo2_result ValidateEncodeBuffers(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes,
    const void* protected_metadata,
    size_t protected_metadata_bytes)
{
    if (!original || !recovery)
        return LEO2_INVALID_ARGUMENT;

    AddressRange metadata_ranges[3];
    size_t metadata_count = 2;
    if (!MakeArrayRange(original, codec->original_count,
            sizeof(*original), metadata_ranges[0]) ||
        !MakeArrayRange(recovery, codec->recovery_count,
            sizeof(*recovery), metadata_ranges[1]))
        return LEO2_INVALID_ARGUMENT;
    if (protected_metadata)
    {
        if (protected_metadata_bytes == 0 ||
            !MakeRange(protected_metadata,
                static_cast<uint64_t>(protected_metadata_bytes),
                metadata_ranges[metadata_count]))
            return LEO2_INVALID_ARGUMENT;
        ++metadata_count;
    }

    AddressRange scratch_range;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangeOverlapsAny(
            scratch_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

    size_t output_count = 0;
    AddressRange single_output = { 0, 0 };
    bool outputs_monotonic = true;
    bool have_previous_output = false;
    uintptr_t previous_output_begin = 0;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (!recovery[i])
            continue;
        AddressRange range;
        if (!MakeRange(recovery[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range) ||
            RangeOverlapsAny(range, metadata_ranges, metadata_count))
            return LEO2_OVERLAP;
        if (have_previous_output && range.begin < previous_output_begin)
            outputs_monotonic = false;
        previous_output_begin = range.begin;
        have_previous_output = true;
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

    bool inputs_monotonic = true;
    bool have_previous_input = false;
    uintptr_t previous_input_begin = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        AddressRange range;
        if (!MakeRange(original[i], shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(range, scratch_range))
            return LEO2_OVERLAP;
        if (!HasValidSystematicPad(codec, original[i], shard_bytes))
            return LEO2_INVALID_ARGUMENT;
        if (have_previous_input && range.begin < previous_input_begin)
            inputs_monotonic = false;
        previous_input_begin = range.begin;
        have_previous_input = true;
    }
    if (inputs_monotonic && outputs_monotonic)
        return ValidateMonotonicEncodePointers(original,
            codec->original_count, recovery, codec->recovery_count,
            shard_bytes);

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

static bool CodecMayUseAutoAVX512Encode(const leo2_codec* codec);

static bool CodecMayUseMeasuredSparseLowGF16AVX2Encode(
    const leo2_codec* codec,
    uint64_t shard_bytes)
{
#ifdef LEO_HAS_FF16
    /*
        Keep the first exact-output promotion deliberately narrower than the
        arithmetic implementation.  Pinned, setup-inclusive measurements
        established a whole-call win for the K=128,R=896 low-profile AVX2
        shape at and above 1 KiB.  Other fields, backends, count neighbors,
        and the 64-byte cell retain the mature prefix evaluator until they
        pass their own public-API crossover gate.

        This predicate reserves bounded call-local schedule storage.  The
        output bitmap is unavailable to the scratch-size query and is checked
        separately immediately before compilation.
    */
    return codec && codec->profile == LEO2_PROFILE_LOW_V1 &&
        codec->field == LEO2_FIELD_GF16 &&
        codec->original_count == 128 && codec->recovery_count == 896 &&
        codec->padded_side == 128 &&
        codec->context && codec->context->backend == LEO2_BACKEND_AVX2 &&
        shard_bytes >= 1024;
#else
    (void)codec;
    (void)shard_bytes;
    return false;
#endif
}

static bool TestForcesSparseEncodeSchedule(const leo2_codec* codec)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    return codec &&
        codec->test_encode_mode == LEO2_TEST_ENCODE_FORCE_TRANSFORM;
#else
    (void)codec;
    return false;
#endif
}

static bool TestBuildReservesSparseEncodeSchedule(const leo2_codec* codec)
{
#ifdef LEO2_ENABLE_TEST_HOOKS
    // The diagnostic codec may switch modes after a scratch query.  Preserve
    // the historical hook-build upper bound so FORCE_TRANSFORM never grows
    // scratch behind the caller's back.
    return codec != NULL;
#else
    (void)codec;
    return false;
#endif
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
    geometry.execution_tile_count = geometry.aligned_prefix_bytes == 0 ? 0 : 1;
    geometry.execution_tile_bytes = geometry.aligned_prefix_bytes;

    /*
        The complete legacy-high GF16 transform repeatedly sweeps 2T shard
        rows.  The measured K=1000, R=200, 64-KiB cell benefits from executing
        two independent 32-KiB ranges end-to-end so the accumulator and
        temporary half remain cache-resident.  Keep promotion exact to that
        AUTO host cell until neighboring sizes and code shapes have equivalent
        same-source evidence.  The public scratch query is independent of the
        requested parity subset, so partial-output calls reuse this compact
        layout even though their conservative transform selector remains
        AVX2; those calls have their own non-regression measurements.  Explicit
        backends, uncalibrated hosts, ragged tails, and every other shape retain
        their single-pass layouts.
        The balanced split avoids a tiny final pass: every tile in this
        qualified range remains larger than 16 KiB, preserving the established
        aligned-pass source-staging policy.

        This is only an execution/scratch layout choice.  Tiles end on the
        64-byte GF16 ALTMAP boundary, and backend selection still occurs once
        from the complete public shard length below.
    */
    static const size_t kHighGF16EncodeTileBytes = 32U * 1024U;
    static const size_t kHighGF16EncodeCellBytes = 64U * 1024U;
    if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->field == LEO2_FIELD_GF16 &&
        codec->original_count == 1000 && codec->recovery_count == 200 &&
        codec->padded_side == 256 &&
        CodecMayUseAutoAVX512Encode(codec) &&
        codec->auto_avx512_encode_ops != NULL &&
        codec->context->ops == codec->context->baseline_ops &&
        geometry.tail_bytes == 0 &&
        geometry.aligned_prefix_bytes == kHighGF16EncodeCellBytes)
    {
        if (!ComputeBalancedExecutionTiles(
                geometry.aligned_prefix_bytes,
                kHighGF16EncodeTileBytes,
                geometry.execution_tile_count,
                geometry.execution_tile_bytes))
            return LEO2_INVALID_COUNTS;
    }

    // Keep the complete GF8 transform schedule unchanged while bounding the
    // live byte slice of all work shards so large calls remain cache-resident.
    // The side-specific crossover policy is generated offline from pinned
    // exact-backend measurements; smaller transforms retain one full pass.
    size_t high_gf8_tile_bytes = 0;
    size_t high_gf8_min_shard_bytes = 0;
    const bool high_gf8_multiple_message_blocks =
        codec->original_count > codec->padded_side;
    switch (codec->padded_side)
    {
    case 16:
        high_gf8_tile_bytes = 128U * 1024U;
        high_gf8_min_shard_bytes = high_gf8_multiple_message_blocks
            ? 1024U * 1024U : 4U * 1024U * 1024U;
        break;
    case 32:
        high_gf8_tile_bytes = 64U * 1024U;
        high_gf8_min_shard_bytes = high_gf8_multiple_message_blocks
            ? 512U * 1024U : 2U * 1024U * 1024U;
        break;
    case 64:
        high_gf8_tile_bytes = 32U * 1024U;
        high_gf8_min_shard_bytes = high_gf8_multiple_message_blocks
            ? 256U * 1024U : 1024U * 1024U;
        break;
    case 128:
        high_gf8_tile_bytes = 64U * 1024U;
        high_gf8_min_shard_bytes = 512U * 1024U;
        break;
    default:
        break;
    }
    if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->field == LEO2_FIELD_GF8 &&
        codec->context &&
        codec->context->ops &&
        codec->context->ops->kind == LEO2_BACKEND_AVX2 &&
        high_gf8_tile_bytes != 0 &&
        geometry.aligned_prefix_bytes >= high_gf8_min_shard_bytes)
    {
        if (!ComputeBalancedExecutionTiles(
                geometry.aligned_prefix_bytes,
                high_gf8_tile_bytes,
                geometry.execution_tile_count,
                geometry.execution_tile_bytes))
            return LEO2_INVALID_COUNTS;
    }
    geometry.work_slot_bytes = std::max(
        geometry.execution_tile_bytes,
        geometry.tail_bytes == 0 ? size_t(0) : kScratchAlignment);
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

    // A qualifying production codec reserves one bounded call-local exact
    // schedule.  Scratch queries cannot depend on the later output bitmap, so
    // dense and unqualified masks for that codec receive the same conservative
    // upper bound while continuing to execute the mature prefix transform.
    // Test FORCE_TRANSFORM keeps its broader diagnostic coverage.
    if (codec->padded_side >= 2 &&
        (CodecMayUseMeasuredSparseLowGF16AVX2Encode(codec, shard_bytes) ||
         TestBuildReservesSparseEncodeSchedule(codec)))
    {
        geometry.sparse_block_count =
            codec->profile == LEO2_PROFILE_LOW_V1
                ? (static_cast<size_t>(codec->recovery_count) +
                    codec->padded_side - 1U) / codec->padded_side
                : 1U;
        geometry.sparse_operation_stride =
            leopard2_internal::SparseForwardRetainedBytes(
                codec->padded_side);
        if (geometry.sparse_operation_stride == 0)
            return LEO2_INTERNAL_ERROR;
        geometry.sparse_bits_offset = geometry.layout.total_bytes;
        size_t retained_bytes = 0;
        size_t schedule_end = 0;
        if (!CheckedMultiply(geometry.sparse_block_count,
                geometry.sparse_operation_stride, retained_bytes) ||
            !CheckedAdd(geometry.sparse_bits_offset, retained_bytes,
                geometry.sparse_workspace_offset) ||
            !CheckedAdd(geometry.sparse_workspace_offset,
                leopard2_internal::SparseForwardDependencyBytes(
                    codec->padded_side),
                schedule_end))
            return LEO2_INVALID_COUNTS;
        if (schedule_end - geometry.sparse_bits_offset >
            kSparseEncodeScheduleBudget)
        {
            // A call-local exact schedule must stay bounded.  These mid-sized
            // low-rate geometries have many parity cosets; until a compact
            // variable-block representation is measured, retain the mature
            // prefix transform instead of expanding caller scratch.
            geometry.sparse_block_count = 0;
            geometry.sparse_operation_stride = 0;
        }
        else
            geometry.layout.total_bytes = schedule_end;
    }
    return LEO2_SUCCESS;
}

static bool UseFusedGenericRevealScatter(
    const leo2_codec* codec,
    size_t aligned_prefix_bytes)
{
#ifdef LEO_HAS_FF8
    // A counterbalanced screen retained the established scratch reveal below
    // 4 KiB and on scalar/NEON.  At and above 4 KiB the fused
    // SSSE3/AVX2/explicit-AVX512 path removes two full output writes and
    // showed a clear whole-decode win.
    return codec->field == LEO2_FIELD_GF8 &&
        aligned_prefix_bytes >= 4096 &&
        (codec->context->backend == LEO2_BACKEND_SSSE3 ||
         codec->context->backend == LEO2_BACKEND_AVX2 ||
         codec->context->backend == LEO2_BACKEND_AVX512);
#else
    (void)codec;
    (void)aligned_prefix_bytes;
    return false;
#endif
}

static bool UseFusedLowRevealScatter(
    const leo2_decode_plan* plan,
    size_t aligned_prefix_bytes)
{
    // Algorithm 4's final values already occupy the requested low-profile
    // systematic coordinates for both materialized and tiled schedules.
    // Complete tiles can therefore apply the inverse-locator factor directly
    // into the non-overlapping caller destination.  Tails retain the scratch
    // reveal because their padded kernel layout is not the public byte layout.
    return (plan->codec->profile == LEO2_PROFILE_LOW_V1 ||
            PlanUsesTranslatedLowDecode(plan)) &&
        aligned_prefix_bytes != 0;
}

static bool SelectTransformDecodePath(
    const leo2_codec* codec,
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    bool multi_item_batch,
    const DecodeScratchGeometry& geometry,
    leopard2_internal::DecodePathInfo& selection)
{
    if (PlanUsesTranslatedLowDecode(plan))
    {
        selection.path =
            (codec->flags & LEO2_CODEC_FORCE_TILED_DECODE) != 0
                ? leopard2_internal::kDecodePathTiled
                : leopard2_internal::kDecodePathMaterialized;
        selection.rule = leopard2_internal::kDecodeRuleTranslatedLow;
        selection.matching_auto_rules = 0;
        // In this exact duality region N=2P, so both low-profile execution
        // forms own N shard slots.  Forced tiled/materialized controls compare
        // traversal cost, not scratch capacity.
        selection.required_work_slots = codec->parent_count;
        selection.aligned_prefix_bytes = geometry.aligned_prefix_bytes;
        selection.tail_bytes = geometry.tail_bytes;
        selection.rounded_shard_bytes = geometry.rounded_bytes;
        selection.multi_item_batch = multi_item_batch;
        return true;
    }
    leopard2_internal::DecodePathInput input;
    input.profile = codec->profile;
    input.field = codec->field;
    input.backend = codec->context->backend;
    input.original_count = codec->original_count;
    input.recovery_count = codec->recovery_count;
    input.padded_side = codec->padded_side;
    input.parent_count = codec->parent_count;
    input.missing_original_count = plan ? plan->missing_original_count : 0;
    input.requested_output_count = codec->profile ==
            LEO2_PROFILE_LEGACY_HIGH_V1
        ? static_cast<uint32_t>(plan
            ? plan->requested_coordinates.size()
            : codec->recovery_count)
        : 0;
    input.codec_flags = codec->flags;
    input.actual_shard_bytes = shard_bytes;
    input.aligned_prefix_bytes = geometry.aligned_prefix_bytes;
    input.tail_bytes = geometry.tail_bytes;
    input.rounded_shard_bytes = geometry.rounded_bytes;
    input.plan_known = plan != NULL;
    input.multi_item_batch = multi_item_batch;
    return leopard2_internal::SelectDecodePath(input, selection);
}

static size_t GF8AVX2DecodeExecutionTileBytes(
    const leo2_codec* codec,
    const leopard2_internal::DecodePathInfo& selection,
    size_t aligned_prefix_bytes)
{
#ifdef LEO_HAS_FF8
    /*
        Every transform stage is pointwise over the shard payload, so complete
        64-byte payload ranges may execute independently without changing the
        field arithmetic or coordinate schedule.  Bound the live transform
        rows to roughly one cache-resident working set on the measured AVX2
        host.  Counterbalanced, single-core whole-decode measurements promote
        only the legacy-high AVX2 traversal and the four side/byte regions
        below.  At each boundary, shortened, punctured, balanced,
        one-above-side, and high-rate counts were checked with small, half,
        and maximum loss patterns.  The weakest promoted observation was
        still 1.065x, while the immediately lower T=16 and T=64 cells were
        neutral or slower.  A translated Algorithm 4 plan retains the
        legacy-high wire identity and is included; unmeasured low-profile and
        forced regular materialized/generic traversals retain the full-payload
        layout.
    */
    if (!codec || !codec->context ||
        codec->field != LEO2_FIELD_GF8 ||
        codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        codec->context->backend != LEO2_BACKEND_AVX2 ||
        aligned_prefix_bytes < 128U * 1024U ||
        selection.path == leopard2_internal::kDecodePathNoOp ||
        selection.path == leopard2_internal::kDecodePathDirect ||
        (selection.path != leopard2_internal::kDecodePathTiled &&
         selection.rule != leopard2_internal::kDecodeRuleTranslatedLow))
        return aligned_prefix_bytes;

    switch (codec->padded_side)
    {
    case 16:
        return aligned_prefix_bytes >= 768U * 1024U
            ? 128U * 1024U : aligned_prefix_bytes;
    case 32:
        return aligned_prefix_bytes >= 384U * 1024U
            ? 64U * 1024U : aligned_prefix_bytes;
    case 64:
        return aligned_prefix_bytes >= 192U * 1024U
            ? 32U * 1024U : aligned_prefix_bytes;
    case 128: return 32U * 1024U;
    default: return aligned_prefix_bytes;
    }
#else
    (void)codec;
    (void)selection;
    return aligned_prefix_bytes;
#endif
}

static leo2_result DecodeLayout(
    const leo2_codec* codec,
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    bool multi_item_batch,
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
    geometry.execution_tile_count =
        geometry.aligned_prefix_bytes == 0 ? 0 : 1;
    geometry.execution_tile_bytes = geometry.aligned_prefix_bytes;
    const size_t range_count = static_cast<size_t>(codec->original_count) * 2 + codec->recovery_count;
    if (!SelectTransformDecodePath(codec, plan, shard_bytes,
            multi_item_batch, geometry, geometry.selection))
        return LEO2_INTERNAL_ERROR;
    const size_t requested_tile_bytes = GF8AVX2DecodeExecutionTileBytes(
        codec, geometry.selection, geometry.aligned_prefix_bytes);
    if (geometry.aligned_prefix_bytes != 0 &&
        requested_tile_bytes != 0 &&
        requested_tile_bytes < geometry.aligned_prefix_bytes)
    {
        if (!ComputeBalancedExecutionTiles(
                geometry.aligned_prefix_bytes, requested_tile_bytes,
                geometry.execution_tile_count,
                geometry.execution_tile_bytes))
            return LEO2_INVALID_COUNTS;
    }
    geometry.work_slot_bytes = std::max(
        geometry.execution_tile_bytes,
        geometry.tail_bytes == 0 ? size_t(0) : kScratchAlignment);
    geometry.work_slot_count = geometry.selection.required_work_slots;
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
    // is needed only for overlap/range validation, never for shard data.  The
    // K=1 legacy-high copy validator has no range table and needs none.
    const size_t range_count =
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
                codec->original_count == 1 && codec->recovery_count == 1
            ? 0
            : static_cast<size_t>(codec->original_count) * 2 +
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

static bool IsLegacyHighK1R1Codec(const leo2_codec* codec)
{
    return codec && codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->original_count == 1 && codec->recovery_count == 1;
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
    // The K=1 legacy-high copy validator compares its sole input/output
    // directly and therefore does not stage validation ranges.
    const size_t range_count =
        IsLegacyHighK1R1Codec(codec)
            ? 0
            : static_cast<size_t>(codec->original_count) +
                codec->recovery_count;
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
    bool use_tiled,
    bool reveal_outputs_in_place,
    void* const* high_systematic_output)
{
    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    const bool use_low_specialized =
        codec->profile == LEO2_PROFILE_LOW_V1 ||
        PlanUsesTranslatedLowDecode(plan);
    const std::vector<uint32_t>& execution_requested_coordinates =
        DecodeExecutionRequestedCoordinates(plan);
    const uint32_t* const requested_coordinates =
        execution_requested_coordinates.data();
    const unsigned requested_count =
        static_cast<unsigned>(execution_requested_coordinates.size());
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
    const leopard2_internal::PrunedTransformBlock* const low_input_plans =
        plan->low_pruned_input_blocks.empty()
            ? NULL : plan->low_pruned_input_blocks.data();
    const unsigned low_input_plan_count = static_cast<unsigned>(
        plan->low_pruned_input_blocks.size());
    const leopard2_internal::PrunedTransformPlan* const low_output_plan =
        plan->low_pruned_output_plan.size == 0
            ? NULL : &plan->low_pruned_output_plan;
    void** const high_requested_output =
        work + static_cast<size_t>(codec->padded_side) * 2;

#ifdef LEO_HAS_FF8
    if (codec->field == LEO2_FIELD_GF8)
    {
        const uint8_t* low_locator = plan->locator8.data();
        const uint8_t* low_factors = codec->fixed_factors8.data();
        if (plan->translated_low)
            low_factors = codec->translated_low_factors8.data();
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
                requested_count, dependencies, &plan->locator8[0],
                reveal_outputs_in_place, work);
        }
        else if (use_low_specialized && use_tiled)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            if (reveal_outputs_in_place)
                leopard::ff8::ReedSolomonDecodeLowTiledPrunedPlanned(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
            else
                leopard::ff8::
                    ReedSolomonDecodeLowTiledPrunedPlannedUnrevealed(
                        ops, buffer_bytes, codec->parent_count,
                        codec->padded_side, coordinate_input,
                        plan->block_input_counts.data(), requested_coordinates,
                        requested_count, dependencies, low_locator,
                        low_factors, low_input_plans,
                        low_input_plan_count, low_output_plan, work);
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
        else if (use_low_specialized)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            if (reveal_outputs_in_place)
                leopard::ff8::ReedSolomonDecodeLowPrunedPlanned(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
            else
                leopard::ff8::ReedSolomonDecodeLowPrunedPlannedUnrevealed(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
        }
        else
            leopard::ff8::ReedSolomonDecodeHighPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator8[0], &codec->fixed_factors8[0],
                high_input_plans, high_input_plan_count,
                high_output_plans, high_output_plan_count,
                high_systematic_output, work);
        return;
    }
#endif
#ifdef LEO_HAS_FF16
    {
        const uint16_t* low_locator = plan->locator16.data();
        const uint16_t* low_factors = codec->fixed_factors16.data();
        if (plan->translated_low)
            low_factors = codec->translated_low_factors16.data();
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
                requested_count, dependencies, &plan->locator16[0],
                reveal_outputs_in_place, work);
        }
        else if (use_low_specialized && use_tiled)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            if (reveal_outputs_in_place)
                leopard::ff16::ReedSolomonDecodeLowTiledPrunedPlanned(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
            else
                leopard::ff16::
                    ReedSolomonDecodeLowTiledPrunedPlannedUnrevealed(
                        ops, buffer_bytes, codec->parent_count,
                        codec->padded_side, coordinate_input,
                        plan->block_input_counts.data(), requested_coordinates,
                        requested_count, dependencies, low_locator,
                        low_factors, low_input_plans,
                        low_input_plan_count, low_output_plan, work);
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
        else if (use_low_specialized)
        {
            const leopard2_internal::OutputDependencyView dependencies =
                leopard2_internal::MakeOutputDependencyView(
                    codec->padded_side,
                    plan->specialized_output_dependencies.data(),
                    plan->specialized_output_dependencies.size());
            if (reveal_outputs_in_place)
                leopard::ff16::ReedSolomonDecodeLowPrunedPlanned(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
            else
                leopard::ff16::ReedSolomonDecodeLowPrunedPlannedUnrevealed(
                    ops, buffer_bytes, codec->parent_count,
                    codec->padded_side, coordinate_input,
                    plan->block_input_counts.data(), requested_coordinates,
                    requested_count, dependencies, low_locator,
                    low_factors, low_input_plans,
                    low_input_plan_count, low_output_plan, work);
        }
        else
            leopard::ff16::ReedSolomonDecodeHighPrunedPlanned(
                ops, buffer_bytes, codec->parent_count, codec->padded_side,
                coordinate_input, plan->block_input_counts.data(),
                requested_coordinates, plan->high_output_blocks.data(),
                static_cast<unsigned>(plan->high_output_blocks.size()),
                &plan->locator16[0], &codec->fixed_factors16[0],
                high_input_plans, high_input_plan_count,
                high_output_plans, high_output_plan_count,
                high_systematic_output, work);
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
        const uint32_t coordinate =
            DecodeExecutionCoordinateForOriginal(plan, i);
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
        const uint32_t coordinate =
            DecodeExecutionCoordinateForRecovery(plan, i);
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

static LEO_FORCE_INLINE void GatherTransformDecodeOne(
    const leo2_decode_plan* plan,
    void* const* restored_original,
    size_t destination_offset,
    size_t pass_bytes,
    void** work,
    bool use_generic,
    bool use_tiled,
    bool outputs_revealed_in_place,
    const leopard::backend::Ops& ops,
    void** high_requested_output,
    size_t output_index)
{
    const leo2_codec* codec = plan->codec;
    const uint32_t original_index = plan->missing_originals[output_index];
    const bool use_low_specialized = !use_generic &&
        (codec->profile == LEO2_PROFILE_LOW_V1 ||
         PlanUsesTranslatedLowDecode(plan));
    const uint32_t execution_coordinate =
        DecodeExecutionCoordinateForOriginal(plan, original_index);
    const void* const source = use_generic || !use_tiled
        ? work[execution_coordinate]
        : use_low_specialized
            ? work[execution_coordinate]
            : high_requested_output[output_index];
    uint8_t* const destination =
        static_cast<uint8_t*>(restored_original[original_index]) +
        destination_offset;
    // Complete kernel tiles can reveal a generic or Algorithm 4 output directly
    // into its caller-owned destination.  This removes both the alias-safe
    // in-place multiply's temporary copy and the following scatter copy.
    // Ragged tails still reveal in scratch first because their 64-byte
    // kernel layout can be larger (and, for GF16, differently arranged)
    // than the public destination.
    if ((use_generic || use_low_specialized) &&
        !outputs_revealed_in_place)
    {
        LEO_DEBUG_ASSERT((pass_bytes & (kScratchAlignment - 1)) == 0);
#ifdef LEO_HAS_FF8
        if (codec->field == LEO2_FIELD_GF8)
        {
            const uint8_t locator_log =
                plan->locator8[execution_coordinate];
            LEO_DEBUG_ASSERT(ops.ff8_multiply != NULL);
            ops.ff8_multiply(destination, source,
                static_cast<uint16_t>(255U - locator_log),
                pass_bytes);
#ifdef LEO2_ENABLE_TEST_HOOKS
            if (use_low_specialized)
                g_test_low_direct_reveal_shards.fetch_add(
                    1, std::memory_order_relaxed);
            else
                g_test_generic_direct_reveal_shards.fetch_add(
                    1, std::memory_order_relaxed);
#endif
            return;
        }
#endif
#ifdef LEO_HAS_FF16
        LEO_DEBUG_ASSERT(codec->field == LEO2_FIELD_GF16 &&
            ops.ff16_multiply != NULL);
        const uint16_t locator_log = plan->locator16[execution_coordinate];
        ops.ff16_multiply(destination, source,
            static_cast<uint16_t>(65535U - locator_log),
            pass_bytes);
#ifdef LEO2_ENABLE_TEST_HOOKS
        if (use_low_specialized)
            g_test_low_direct_reveal_shards.fetch_add(
                1, std::memory_order_relaxed);
        else
            g_test_generic_direct_reveal_shards.fetch_add(
                1, std::memory_order_relaxed);
#endif
        return;
#endif
    }
#ifdef LEO2_ENABLE_TEST_HOOKS
    if (use_low_specialized)
        g_test_low_scratch_reveal_shards.fetch_add(
            1, std::memory_order_relaxed);
#endif
    GatherShardFromKernel(codec, destination, source, pass_bytes);
}

static void GatherTransformDecodePass(
    const leo2_decode_plan* plan,
    void* const* restored_original,
    size_t destination_offset,
    size_t pass_bytes,
    void** work,
    bool use_generic,
    bool use_tiled,
    bool outputs_revealed_in_place)
{
    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    void** const high_requested_output =
        work + static_cast<size_t>(codec->padded_side) * 2;

#if defined(_OPENMP) && defined(LEO_HAS_FF16)
    const bool direct_ff16_reveal = !outputs_revealed_in_place &&
        (use_generic || codec->profile == LEO2_PROFILE_LOW_V1 ||
         PlanUsesTranslatedLowDecode(plan)) &&
        codec->field == LEO2_FIELD_GF16;
    if (direct_ff16_reveal)
    {
#pragma omp parallel for
        for (int i = 0;
             i < static_cast<int>(plan->missing_originals.size()); ++i)
            GatherTransformDecodeOne(
                plan, restored_original, destination_offset, pass_bytes, work,
                use_generic, use_tiled, outputs_revealed_in_place, ops,
                high_requested_output, static_cast<size_t>(i));
        return;
    }
#endif
    for (size_t i = 0; i < plan->missing_originals.size(); ++i)
        GatherTransformDecodeOne(
            plan, restored_original, destination_offset, pass_bytes, work,
            use_generic, use_tiled, outputs_revealed_in_place, ops,
            high_requested_output, i);
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

static bool MatchesMeasuredSparseLowGF16OutputMask(
    const leo2_codec* codec,
    void* const* recovery,
    uint32_t requested_recovery_count)
{
    if (!codec || !recovery || codec->recovery_count != 896)
        return false;

    static const uint16_t kEdgeSparse[] = {
        0, 127, 128, 383, 895
    };
    static const uint16_t kScatteredSparse[] = {
        7, 63, 135, 255, 519, 895
    };
    const uint16_t* expected = NULL;
    size_t expected_count = 0;
    if (requested_recovery_count ==
        sizeof(kEdgeSparse) / sizeof(kEdgeSparse[0]))
    {
        expected = kEdgeSparse;
        expected_count = sizeof(kEdgeSparse) / sizeof(kEdgeSparse[0]);
    }
    else if (requested_recovery_count ==
        sizeof(kScatteredSparse) / sizeof(kScatteredSparse[0]))
    {
        expected = kScatteredSparse;
        expected_count =
            sizeof(kScatteredSparse) / sizeof(kScatteredSparse[0]);
    }
    else
        return false;

    // EncodeInternal already counted every non-null output.  Matching all Q
    // expected entries therefore proves there cannot be an unexpected extra
    // output and avoids another O(R) pointer-array pass in the hot call.
    for (size_t i = 0; i < expected_count; ++i)
        if (!recovery[expected[i]])
            return false;
    return true;
}

static bool PrepareSparseEncodePlans(
    const leo2_codec* codec,
    const EncodeScratchGeometry& geometry,
    void* const* recovery,
    uint32_t requested_recovery_count,
    uint8_t* scratch_base,
    leopard2_internal::SparseForwardPlanBatchView& view)
{
    view.operation_masks = NULL;
    view.operation_stride = 0;
    view.block_count = 0;
    const bool forced = TestForcesSparseEncodeSchedule(codec);
    const bool measured =
        CodecMayUseMeasuredSparseLowGF16AVX2Encode(codec,
            geometry.aligned_prefix_bytes + geometry.tail_bytes) &&
        MatchesMeasuredSparseLowGF16OutputMask(
            codec, recovery, requested_recovery_count);
    if (!forced && !measured)
        return true;
    // A one-symbol side has no butterflies to schedule.  Keep the empty view
    // valid so the existing XOR/copy direct transform path remains available.
    if (codec->padded_side < 2)
        return true;
    // EncodeLayout deliberately publishes an empty pair when the exact masks
    // would exceed the bounded call-local budget.  That is a valid prefix
    // fallback, not an execution error.
    if (geometry.sparse_block_count == 0 &&
        geometry.sparse_operation_stride == 0)
        return true;
    if (geometry.sparse_block_count == 0 ||
        geometry.sparse_operation_stride == 0 ||
        geometry.sparse_block_count > 0xffffffffu)
        return false;

    uint8_t* const operation_masks =
        scratch_base + geometry.sparse_bits_offset;
    uint8_t* const workspace =
        scratch_base + geometry.sparse_workspace_offset;

    const uint32_t side = codec->padded_side;
    const size_t dependency_bytes =
        leopard2_internal::SparseForwardDependencyBytes(side);
    if (dependency_bytes == 0)
        return false;
    for (size_t block = 0; block < geometry.sparse_block_count; ++block)
    {
        const uint32_t recovery_offset = codec->profile == LEO2_PROFILE_LOW_V1
            ? static_cast<uint32_t>(block) * side
            : 0;
        const uint32_t block_count = std::min<uint32_t>(
            side, codec->recovery_count - recovery_offset);
        uint32_t requested_count = 0;
        memset(workspace, 0, dependency_bytes);
        for (uint32_t i = 0; i < block_count; ++i)
        {
            if (!recovery[recovery_offset + i])
                continue;
            workspace[i >> 3] |= static_cast<uint8_t>(1U << (i & 7U));
            ++requested_count;
        }
        if (requested_count == 0)
            continue;

        leopard2_internal::SparseForwardPlanStats stats;
        uint8_t* const block_bits = operation_masks +
            block * geometry.sparse_operation_stride;
        const uint32_t shift = codec->profile == LEO2_PROFILE_LOW_V1
            ? side + recovery_offset
            : 0;
        bool compiled = false;
#ifdef LEO_HAS_FF8
        if (codec->field == LEO2_FIELD_GF8)
        {
            compiled = leopard::ff8::PrepareSparseForwardPlan(
                side, shift, workspace, dependency_bytes, block_bits,
                geometry.sparse_operation_stride, stats);
        }
#endif
#ifdef LEO_HAS_FF16
        if (codec->field == LEO2_FIELD_GF16)
        {
            compiled = leopard::ff16::PrepareSparseForwardPlan(
                side, shift, workspace, dependency_bytes, block_bits,
                geometry.sparse_operation_stride, stats);
        }
#endif
        if (!compiled)
            return false;

    }

    view.operation_masks = operation_masks;
    view.operation_stride = geometry.sparse_operation_stride;
    view.block_count = static_cast<uint32_t>(geometry.sparse_block_count);
    return true;
}

static void ExecuteTransformEncodePass(
    const leo2_codec* codec,
    const leopard::backend::Ops& ops,
    size_t buffer_bytes,
    size_t policy_buffer_bytes,
    uint32_t requested_recovery_count,
    uint32_t requested_recovery_prefix,
    const void* const* padded_original,
    void** parity,
    void** work,
    const leopard2_internal::SparseForwardPlanBatchView* sparse_plans)
{
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
                requested_recovery_prefix, requested_recovery_count,
                codec->padded_side,
                padded_original, work, sparse_plans);
#else
            return;
#endif
        else
#ifdef LEO_HAS_FF16
            leopard::ff16::ReedSolomonEncodeWithSourcePolicy(
                ops, buffer_bytes, policy_buffer_bytes,
                codec->original_count,
                requested_recovery_prefix, requested_recovery_count,
                codec->padded_side,
                padded_original, work, sparse_plans);
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
            parity, work, sparse_plans);
#else
        return;
#endif
    else
#ifdef LEO_HAS_FF16
        leopard::ff16::ReedSolomonEncodeLow(
            ops, buffer_bytes, codec->original_count,
            codec->recovery_count, codec->padded_side, padded_original,
            parity, work, sparse_plans);
#else
        return;
#endif
}

static bool CodecMayUseAutoAVX512Encode(const leo2_codec* codec)
{
    if (!codec || !codec->context ||
        !codec->context->auto_requested ||
        !codec->context->auto_avx512_encode_host ||
        codec->context->backend != LEO2_BACKEND_AVX2 ||
        codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1)
        return false;

    if (codec->field == LEO2_FIELD_GF16)
    {
        if (codec->original_count < 8 ||
            codec->recovery_count > 4096 ||
            codec->parent_count < 16)
            return false;

        // T=1 uses the direct copy/XOR path and never reaches this transform
        // selector.  The remaining bounds exactly match the retained Zen 5
        // calibration screen: K >= 8, R <= 4096, and N >= 16.
        return codec->padded_side >= 2;
    }

    if (codec->field != LEO2_FIELD_GF8)
        return false;

    const uint32_t side = codec->padded_side;
    if (side == 8)
    {
        // The register-resident T=8 kernel clears the exact Leopard1 gate for
        // the balanced code.  Keep shortened and punctured neighbors out of
        // AUTO: their explicit-callback correctness is covered, but their
        // small-shard exact-main results did not clear the promotion margin.
        return codec->original_count == side &&
            codec->recovery_count == side;
    }
    if (side == 64)
    {
        // The same-source AVX2 gate qualified all four combinations of an
        // exact or one-coordinate-shortened message block and an exact or
        // one-coordinate-punctured parity block.  UseAutoAVX512Encode applies
        // the stricter exact-Leopard-main byte gate to comparable R <= K
        // neighbors.  The K=63,R=64 shortened shape has no legal Leopard1
        // comparison because its recovery count exceeds its original count.
        return (codec->original_count == side ||
                codec->original_count + 1U == side) &&
            (codec->recovery_count == side ||
             codec->recovery_count + 1U == side);
    }

    // Retain the earlier exact balanced T=16 and T=32 gates.  Their ragged
    // neighbors were not part of the isolated promotion screen.
    return (side == 16 || side == 32) &&
        codec->original_count == side &&
        codec->recovery_count == side;
}

static bool UseAutoAVX512Encode(
    const leo2_codec* codec,
    size_t buffer_bytes,
    uint32_t requested_recovery_count,
    uint32_t requested_recovery_prefix)
{
    if (!CodecMayUseAutoAVX512Encode(codec) ||
        !codec->auto_avx512_encode_ops ||
        codec->context->ops != codec->context->baseline_ops ||
        requested_recovery_count != codec->recovery_count ||
        requested_recovery_prefix != codec->recovery_count ||
        (buffer_bytes & (kScratchAlignment - 1U)) != 0)
        return false;

    if (codec->field == LEO2_FIELD_GF8)
    {
        // The whole-transform AVX-512VL calibration screen is positive in
        // these cache-resident regions.  Exact T=8 clears the exact-main gate
        // at 2 KiB and from 8 through 64 KiB, but its negative 4-KiB cell
        // requires a disjoint predicate.  Exact T=32 and the previously
        // qualified exact T=64 shape clear it at 2 KiB; T=16 needs 4 KiB.  A
        // one-row-punctured T=64 code clears the exact-main threshold only at
        // 64 KiB.  The shortened K=63,R=64 shape has no legal Leopard1
        // counterpart and retains the positive same-source AVX2 interval.
        // Above the upper bound the backends converge as shard traffic becomes
        // the bottleneck.
        if (codec->padded_side == 64 &&
            codec->recovery_count + 1U == codec->padded_side)
            return buffer_bytes == 64U * 1024U;
        if (codec->padded_side == 8)
            return buffer_bytes == 2U * 1024U ||
                (buffer_bytes >= 8U * 1024U &&
                 buffer_bytes <= 64U * 1024U);
        const size_t minimum_bytes = codec->padded_side == 16
            ? 4U * 1024U : 2U * 1024U;
        return buffer_bytes >= minimum_bytes &&
            buffer_bytes <= 64U * 1024U;
    }

    if (buffer_bytes > 4U * 1024U * 1024U)
        return false;

    // The AVX-512VL translation unit deliberately retains 256-bit data
    // operations and uses the expanded register file.  Counterbalanced
    // screens across both cache complexes were positive at every complete
    // tile from 64 B through 4 MiB inclusive.  Larger shards and ragged tails
    // retain AVX2 until they have equally broad evidence.
    return buffer_bytes >= kScratchAlignment;
}

static const leopard::backend::Ops& SelectTransformEncodeOps(
    const leo2_codec* codec,
    size_t buffer_bytes,
    uint32_t requested_recovery_count,
    uint32_t requested_recovery_prefix)
{
    if (UseAutoAVX512Encode(codec, buffer_bytes,
            requested_recovery_count, requested_recovery_prefix))
        return *codec->auto_avx512_encode_ops;
    return *codec->context->ops;
}

static bool UseAVX2GF8PairwiseLargeR1Xor(
    uint32_t original_count,
    size_t shard_bytes);

static bool UseCoarseR1Xor(
    const leo2_codec* codec,
    size_t shard_bytes)
{
    if (shard_bytes < 1024)
        return false;
    uint32_t minimum_originals = 0;
    switch (codec->context->backend)
    {
    case LEO2_BACKEND_SCALAR:
    case LEO2_BACKEND_SSSE3:
        minimum_originals = 5;
        break;
    case LEO2_BACKEND_AVX2:
        minimum_originals = 7;
        break;
    case LEO2_BACKEND_AVX512:
        /* The AVX-512 translation unit currently reuses the AVX2 kernel, but
           that is not enough evidence to promote this end-to-end codec path:
           frequency effects and the distinct execution environment still
           need a retained, reproducible benchmark matrix.  Keep the mature
           paired reduction until that gate has passed. */
        return false;
    default:
        /* Native NEON and future backends retain the mature paired loop until
           a coarse implementation has backend-specific measurements. */
        return false;
    }
    if (codec->original_count < minimum_originals)
        return false;
    if (codec->context->backend != LEO2_BACKEND_AVX2 ||
        codec->original_count <= 7 || codec->field != LEO2_FIELD_GF8 ||
        shard_bytes < 1024U * 1024U)
        return true;
    return !UseAVX2GF8PairwiseLargeR1Xor(
        codec->original_count, shard_bytes);
}

#if defined(_MSC_VER)
#define LEO2_R1_FUSED_POLICY_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__ELF__)
#define LEO2_R1_FUSED_POLICY_NOINLINE \
    __attribute__((noinline, section(".text.leo2_r1_fused_policy")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_R1_FUSED_POLICY_NOINLINE __attribute__((noinline))
#else
#define LEO2_R1_FUSED_POLICY_NOINLINE
#endif

static LEO2_R1_FUSED_POLICY_NOINLINE bool UseFusedFinalR1Xor(
    uint32_t original_count,
    size_t shard_bytes)
{
    /* The remainder-specific callback folds the final 3..6 sources into the
       coarse reduction instead of returning to the generic loop.  Retained
       Zen 5 AVX2 screens promoted the exact K groups below.  K=7 crosses back
       to the mature coarse loop above 4 MiB; K=20/21 missed the 5% threshold;
       K=8 Group7, K=31, explicit prefetch, and byte-loop unrolling all
       regressed and remain excluded. */
    /* This layout-isolated policy is reached only after the inline eligibility
       gate has proved the field, backend, callback, lower byte bound, and K
       set.  Keep those preconditions explicit while avoiding duplicate loads
       and branches at both encode and decode call sites. */
    LEO_DEBUG_ASSERT(shard_bytes >= 4096);
    LEO_DEBUG_ASSERT(original_count == 7 ||
        (original_count >= 12 && original_count <= 15) ||
        (original_count >= 22 && original_count <= 23));
    if (original_count == 7)
        return shard_bytes <= 4U * 1024U * 1024U;
    return shard_bytes <= 1024U * 1024U;
}

#undef LEO2_R1_FUSED_POLICY_NOINLINE

static LEO_FORCE_INLINE bool IsFusedFinalR1XorEligible(
    const leo2_codec* codec,
    const leopard::backend::Ops& ops,
    size_t shard_bytes)
{
#if defined(LEO2_R1_FUSED_FINAL_CONTROL)
    /* Diagnostic same-source benchmark control.  This definition is never
       enabled by the default build: it preserves every implementation object
       while disabling only the new dispatch edge so retained ABBA evidence
       can attribute the difference to the fused-final selector. */
    (void)codec;
    (void)ops;
    (void)shard_bytes;
    return false;
#else
    if (shard_bytes < 4096 || codec->field != LEO2_FIELD_GF8 ||
        codec->context->backend != LEO2_BACKEND_AVX2 ||
        !ops.xor_memory_sources_fused_final)
        return false;

    const uint32_t k = codec->original_count;
    return k == 7 || (k >= 12 && k <= 15) || (k >= 22 && k <= 23);
#endif
}

static LEO_FORCE_INLINE bool UseDenseR1Xor(
    const leo2_codec* codec,
    size_t shard_bytes)
{
    const bool eligible =
        (codec->original_count == 2 && shard_bytes >= 4096) ||
        (codec->original_count == 4 && shard_bytes >= 2048);
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_expect(!eligible, 1))
#else
    if (!eligible)
#endif
        return false;
    return codec->context->ops->xor_memory_dense != NULL;
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
    if (UseCoarseR1Xor(codec, shard_bytes))
    {
        const leopard::backend::XorMemorySources reduce =
            IsFusedFinalR1XorEligible(codec, ops, shard_bytes) &&
                    UseFusedFinalR1Xor(
                        codec->original_count, shard_bytes) ?
                ops.xor_memory_sources_fused_final : ops.xor_memory_sources;
        reduce(recovery[0], original[0], original + 1,
            codec->original_count - 1, shard_bytes);
        return;
    }
    if (UseDenseR1Xor(codec, shard_bytes))
    {
        ops.xor_memory_dense(recovery[0], original,
            codec->original_count, shard_bytes);
        return;
    }
    uint8_t* const output = static_cast<uint8_t*>(recovery[0]);
    memcpy(output, original[0], shard_bytes);
    uint32_t i = 1;
    for (; i + 1 < codec->original_count; i += 2)
        leopard::xor_mem_2to1(
            ops, output, original[i], original[i + 1], shard_bytes);
    if (i < codec->original_count)
        leopard::xor_mem(ops, output, original[i], shard_bytes);
}

static LEO_FORCE_INLINE void ExecuteDirectXorDecode(
    const leo2_codec* codec,
    uint32_t missing,
    size_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original)
{
    const leopard::backend::Ops& ops = *codec->context->ops;
    if (UseCoarseR1Xor(codec, shard_bytes))
    {
        const leopard::backend::XorMemorySources reduce =
            IsFusedFinalR1XorEligible(codec, ops, shard_bytes) &&
                    UseFusedFinalR1Xor(
                        codec->original_count, shard_bytes) ?
                ops.xor_memory_sources_fused_final : ops.xor_memory_sources;
        reduce(restored_original[missing], recovery[0],
            original, codec->original_count, shard_bytes);
        return;
    }
    if (UseDenseR1Xor(codec, shard_bytes))
    {
        const void* sources[4];
        uint32_t source_count = 0;
        sources[source_count++] = recovery[0];
        for (uint32_t i = 0; i < codec->original_count; ++i)
            if (original[i])
                sources[source_count++] = original[i];
        LEO_DEBUG_ASSERT(source_count == codec->original_count);
        ops.xor_memory_dense(restored_original[missing], sources,
            source_count, shard_bytes);
        return;
    }
    uint8_t* output =
        static_cast<uint8_t*>(restored_original[missing]);
    memcpy(output, recovery[0], shard_bytes);
    const void* waiting = NULL;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if (!original[i])
            continue;
        if (!waiting)
            waiting = original[i];
        else
        {
            leopard::xor_mem_2to1(
                ops, output, waiting, original[i], shard_bytes);
            waiting = NULL;
        }
    }
    if (waiting)
        leopard::xor_mem(ops, output, waiting, shard_bytes);
}

static leo2_result ValidateK1R1DecodeBuffers(
    uint64_t shard_bytes,
    const AddressRange& scratch_range,
    const AddressRange* protected_ranges,
    size_t protected_count,
    const void* const* original,
    const void* const* recovery,
    void* const* restored)
{
    if (!original || !recovery || !restored)
        return LEO2_INVALID_ARGUMENT;
    if (protected_count > 2 || (protected_count != 0 && !protected_ranges))
        return LEO2_INVALID_ARGUMENT;

    AddressRange metadata_ranges[5];
    if (!MakeArrayRange(original, 1, sizeof(*original), metadata_ranges[0]) ||
        !MakeArrayRange(recovery, 1, sizeof(*recovery), metadata_ranges[1]) ||
        !MakeArrayRange(restored, 1, sizeof(*restored), metadata_ranges[2]))
        return LEO2_INVALID_ARGUMENT;
    for (size_t i = 0; i < protected_count; ++i)
        metadata_ranges[3 + i] = protected_ranges[i];
    const size_t metadata_count = 3 + protected_count;
    if (RangeOverlapsAny(scratch_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

    AddressRange output_range;
    if (!restored[0] ||
        !MakeRange(restored[0], shard_bytes, output_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangesOverlap(output_range, scratch_range) ||
        RangeOverlapsAny(output_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;
    if (original[0])
        return LEO2_INVALID_ARGUMENT;
    if (!recovery[0])
        return LEO2_INVALID_ARGUMENT;

    AddressRange recovery_range;
    if (!MakeRange(recovery[0], shard_bytes, recovery_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangesOverlap(recovery_range, scratch_range) ||
        RangesOverlap(recovery_range, output_range))
        return LEO2_OVERLAP;

    return LEO2_SUCCESS;
}

static leo2_result DecodeDirectXorValidated(
    const leo2_codec* codec,
    uint32_t missing,
    uint64_t shard_bytes,
    const AddressRange* protected_ranges,
    size_t protected_count,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    size_t scratch_bytes)
{
    if (shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (IsLegacyHighK1R1Codec(codec) && missing == 0)
    {
        size_t rounded_bytes = 0;
        if (!RoundShardBytes(shard_bytes, rounded_bytes))
            return LEO2_INVALID_ARGUMENT;
        if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
            return LEO2_UNSUPPORTED;
        AddressRange scratch_range;
        leo2_result result = CheckOptionalScratch(
            scratch, scratch_bytes, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        result = ValidateK1R1DecodeBuffers(
            shard_bytes, scratch_range, protected_ranges, protected_count,
            original, recovery, restored);
        if (result != LEO2_SUCCESS)
            return result;
        memcpy(restored[0], recovery[0], static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }
    ScratchLayout layout;
    size_t rounded_bytes = 0;
    leo2_result result = DirectDecodeLayout(
        codec, shard_bytes, layout, rounded_bytes);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange scratch_range = { 0, 0 };
    result = CheckScratch(
        scratch, scratch_bytes, layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    if (codec->original_count == 1 && codec->recovery_count == 1 &&
        missing == 0)
    {
        result = ValidateK1R1DecodeBuffers(
            shard_bytes, scratch_range, protected_ranges, protected_count,
            original, recovery, restored);
        if (result != LEO2_SUCCESS)
            return result;
        memcpy(restored[0], recovery[0], static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }
    if (!original || !recovery || !restored ||
        missing >= codec->original_count)
        return LEO2_INVALID_ARGUMENT;

    if (protected_count > 2 || (protected_count != 0 && !protected_ranges))
        return LEO2_INVALID_ARGUMENT;
    AddressRange metadata_ranges[5];
    if (!MakeArrayRange(original, codec->original_count,
            sizeof(*original), metadata_ranges[0]) ||
        !MakeArrayRange(recovery, codec->recovery_count,
            sizeof(*recovery), metadata_ranges[1]) ||
        !MakeArrayRange(restored, codec->original_count,
            sizeof(*restored), metadata_ranges[2]))
        return LEO2_INVALID_ARGUMENT;
    for (size_t i = 0; i < protected_count; ++i)
        metadata_ranges[3 + i] = protected_ranges[i];
    const size_t metadata_count = 3 + protected_count;
    if (RangeOverlapsAny(
            scratch_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

    AddressRange output_range;
    if (!restored[missing] ||
        !MakeRange(restored[missing], shard_bytes, output_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangesOverlap(output_range, scratch_range) ||
        RangeOverlapsAny(output_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        const bool expected = i != missing;
        if ((original[i] != NULL) != expected)
            return LEO2_INVALID_ARGUMENT;
        if (!expected)
            continue;
        AddressRange input_range;
        if (!MakeRange(original[i], shard_bytes, input_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(input_range, scratch_range) ||
            RangesOverlap(input_range, output_range))
            return LEO2_OVERLAP;
        if (!HasValidSystematicPad(codec, original[i], shard_bytes))
            return LEO2_INVALID_ARGUMENT;
    }

    if (!recovery[0])
        return LEO2_INVALID_ARGUMENT;
    AddressRange recovery_range;
    if (!MakeRange(recovery[0], shard_bytes, recovery_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangesOverlap(recovery_range, scratch_range) ||
        RangesOverlap(recovery_range, output_range))
        return LEO2_OVERLAP;

    ExecuteDirectXorDecode(codec, missing, static_cast<size_t>(shard_bytes),
        original, recovery, restored);
    return LEO2_SUCCESS;
}

static leo2_result ValidateDecodeBuffers(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored,
    void* scratch,
    size_t scratch_bytes,
    const ProtectedMetadataSpan* protected_metadata,
    size_t protected_metadata_count,
    void** coordinate_data)
{
    const leo2_codec* codec = plan->codec;
    if (!original || !recovery || !restored)
        return LEO2_INVALID_ARGUMENT;

    /* Three pointer arrays plus the one-shot presence arrays (or one batch
       item-array span).  All are immutable metadata that writable scratch
       and restored shard ranges must not overlap. */
    AddressRange metadata_ranges[5];
    size_t metadata_count = 3;
    if (!MakeArrayRange(original, codec->original_count,
            sizeof(*original), metadata_ranges[0]) ||
        !MakeArrayRange(recovery, codec->recovery_count,
            sizeof(*recovery), metadata_ranges[1]) ||
        !MakeArrayRange(restored, codec->original_count,
            sizeof(*restored), metadata_ranges[2]))
        return LEO2_INVALID_ARGUMENT;
    if (protected_metadata_count > 2 ||
        (protected_metadata_count != 0 && !protected_metadata))
        return LEO2_INVALID_ARGUMENT;
    for (size_t i = 0; i < protected_metadata_count; ++i)
    {
        if (!protected_metadata[i].data || protected_metadata[i].bytes == 0 ||
            !MakeRange(protected_metadata[i].data,
                static_cast<uint64_t>(protected_metadata[i].bytes),
                metadata_ranges[metadata_count]))
            return LEO2_INVALID_ARGUMENT;
        ++metadata_count;
    }

    AddressRange scratch_range;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
    if (RangeOverlapsAny(
            scratch_range, metadata_ranges, metadata_count))
        return LEO2_OVERLAP;

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
        /* Compact R=1 plans deliberately omit every pattern vector.  Batch
           preflight still shares this validator, so recover the sole output
           from the compact scalar instead of requiring missing_originals. */
        if (PlanHasCompactDirectXor(plan))
            single_output_index = plan->direct_xor_original;
        else
        {
            if (plan->missing_originals.size() != 1)
                return LEO2_INTERNAL_ERROR;
            single_output_index = plan->missing_originals[0];
        }
        if (single_output_index >= codec->original_count ||
            !restored[single_output_index] ||
            !MakeRange(
                restored[single_output_index], shard_bytes,
                single_output_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(single_output_range, scratch_range))
            return LEO2_OVERLAP;
        if (RangeOverlapsAny(
                single_output_range, metadata_ranges, metadata_count))
            return LEO2_OVERLAP;
    }

    size_t input_count = 0;
    size_t output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if ((original[i] != NULL) != PlanOriginalPresent(plan, i))
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
            if (RangeOverlapsAny(range, metadata_ranges, metadata_count))
                return LEO2_OVERLAP;
            ++output_count;
        }
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if ((recovery[i] != NULL) != PlanRecoveryPresent(plan, i))
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
    if (coordinate_data)
        std::fill(coordinate_data,
            coordinate_data + codec->parent_count,
            static_cast<void*>(NULL));
    input_count = 0;
    output_count = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if (original[i])
        {
            MakeRange(original[i], shard_bytes, ranges[input_count++]);
            if (coordinate_data)
            {
                const uint32_t coordinate =
                    DecodeExecutionCoordinateForOriginal(plan, i);
                if (!plan->coordinate_erased[coordinate])
                    coordinate_data[coordinate] =
                        const_cast<void*>(original[i]);
            }
        }
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (recovery[i])
        {
            MakeRange(recovery[i], shard_bytes, ranges[input_count++]);
            if (coordinate_data)
            {
                const uint32_t coordinate =
                    DecodeExecutionCoordinateForRecovery(plan, i);
                if (!plan->coordinate_erased[coordinate])
                    coordinate_data[coordinate] =
                        const_cast<void*>(recovery[i]);
            }
        }
    }
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

#ifdef LEO_HAS_FF8
#if defined(_MSC_VER)
#define LEO2_SOURCE_MAJOR_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_SOURCE_MAJOR_NOINLINE __attribute__((noinline))
#else
#define LEO2_SOURCE_MAJOR_NOINLINE
#endif

static LEO2_SOURCE_MAJOR_NOINLINE leo2_result
ExecuteGF8SourceMajorDirectRepair(
    const leo2_decode_plan* plan,
    size_t shard_bytes,
    const void* const* original,
    const void* const* recovery,
    void* const* restored_original)
{
    const leo2_codec* codec = plan->codec;
    const leopard::backend::Ops& ops = *codec->context->ops;
    const uint32_t output_count = plan->missing_original_count;
    if (output_count < 2 || output_count > kDirectMaxRepairLosses ||
        codec->original_count > kDirectMaxParentDimension ||
        codec->recovery_count > kDirectMaxParentDimension ||
        plan->missing_originals.size() != output_count ||
        plan->direct_term_offsets.size() !=
            static_cast<size_t>(output_count) + 1U ||
        plan->direct_term_offsets.front() != 0 ||
        plan->direct_term_offsets.back() != plan->direct_terms.size())
        return LEO2_INTERNAL_ERROR;

    // Preserve the compact immutable output-major plan because it is faster
    // for short shards and has no extra cold-plan allocation.  For a large
    // execution, transpose its at-most K*L coefficients into a bounded stack
    // schedule.  Direct repair uses L selected parities and K-L surviving
    // originals, so the union contains at most K rows.
    uint32_t sources[kDirectMaxParentDimension] = {};
    uint16_t original_source_rows[kDirectMaxParentDimension];
    uint16_t recovery_source_rows[kDirectMaxParentDimension];
    uint16_t source_logs[
        kDirectMaxParentDimension * kDirectMaxRepairLosses];
    std::fill(original_source_rows,
        original_source_rows + codec->original_count, UINT16_MAX);
    std::fill(recovery_source_rows,
        recovery_source_rows + codec->recovery_count, UINT16_MAX);
    const size_t maximum_source_logs =
        static_cast<size_t>(codec->original_count) * output_count;
    std::fill(source_logs, source_logs + maximum_source_logs, UINT16_MAX);
    size_t source_count = 0;
    for (uint32_t output_index = 0;
         output_index < output_count; ++output_index)
    {
        const size_t begin = plan->direct_term_offsets[output_index];
        const size_t end = plan->direct_term_offsets[output_index + 1];
        if (begin >= end || end > plan->direct_terms.size())
            return LEO2_INTERNAL_ERROR;

        // The output-major plan places a unit coefficient first when one
        // exists.  Preserve that efficient copy/multiply initializer and
        // transpose only the remaining accumulation.
        for (size_t term_index = begin + 1;
             term_index < end; ++term_index)
        {
            const leo2_direct_repair_term& term =
                plan->direct_terms[term_index];
            const bool parity =
                (term.tagged_source & kDirectRecoveryTag) != 0;
            const uint32_t source_index =
                term.tagged_source & ~kDirectRecoveryTag;
            if ((parity && source_index >= codec->recovery_count) ||
                (!parity && source_index >= codec->original_count))
                return LEO2_INTERNAL_ERROR;
            uint16_t* source_row = parity
                ? &recovery_source_rows[source_index]
                : &original_source_rows[source_index];
            if (*source_row == UINT16_MAX)
            {
                if (source_count >= codec->original_count)
                    return LEO2_INTERNAL_ERROR;
                *source_row = static_cast<uint16_t>(source_count);
                sources[source_count++] = term.tagged_source;
            }
            source_logs[
                static_cast<size_t>(*source_row) * output_count +
                    output_index] = term.multiplier_log;
        }
    }
    if (source_count == 0)
        return LEO2_INTERNAL_ERROR;

    void* outputs[kDirectMaxRepairLosses] = {};
    for (uint32_t output_index = 0;
         output_index < output_count; ++output_index)
    {
        const uint32_t missing_original =
            plan->missing_originals[output_index];
        if (missing_original >= codec->original_count)
            return LEO2_INTERNAL_ERROR;
        outputs[output_index] = restored_original[missing_original];
        if (!outputs[output_index])
            return LEO2_INTERNAL_ERROR;
        const size_t begin = plan->direct_term_offsets[output_index];
        const leo2_direct_repair_term& initial = plan->direct_terms[begin];
        const bool parity =
            (initial.tagged_source & kDirectRecoveryTag) != 0;
        const uint32_t source_index =
            initial.tagged_source & ~kDirectRecoveryTag;
        if ((parity && source_index >= codec->recovery_count) ||
            (!parity && source_index >= codec->original_count))
            return LEO2_INTERNAL_ERROR;
        const void* source = parity
            ? recovery[source_index] : original[source_index];
        if (!source)
            return LEO2_INTERNAL_ERROR;
        if (initial.multiplier_log == 0)
            memcpy(outputs[output_index], source, shard_bytes);
        else
            leopard::ff8::MultiplyBytes(ops, outputs[output_index], source,
                static_cast<leopard::ff8::ffe_t>(initial.multiplier_log),
                shard_bytes);
    }

    for (size_t source_row = 0;
         source_row < source_count; ++source_row)
    {
        const uint32_t tagged_source = sources[source_row];
        const bool parity = (tagged_source & kDirectRecoveryTag) != 0;
        const uint32_t source_index = tagged_source & ~kDirectRecoveryTag;
        const void* source = parity
            ? recovery[source_index] : original[source_index];
        if (!source)
            return LEO2_INTERNAL_ERROR;
        ops.ff8_multiply_add_outputs(outputs, source,
            source_logs + source_row * output_count,
            output_count, shard_bytes);
    }
    return LEO2_SUCCESS;
}

#undef LEO2_SOURCE_MAJOR_NOINLINE
#endif

static leopard2_internal::DirectRepairExecutor SelectDirectRepairExecutor(
    const leo2_decode_plan* plan,
    size_t shard_bytes)
{
    if (!plan || !plan->direct_repair || !plan->codec)
        return leopard2_internal::kDirectRepairExecutorNone;

#ifdef LEO_HAS_FF8
    const leo2_codec* codec = plan->codec;
    const uint32_t output_count = plan->missing_original_count;
    const leopard::backend::Ops& ops = *codec->context->ops;
    const bool representable =
        output_count >= 2 && output_count <= kDirectMaxRepairLosses &&
        codec->field == LEO2_FIELD_GF8 &&
        ops.ff8_multiply_add_outputs &&
        codec->original_count <= kDirectMaxParentDimension &&
        codec->recovery_count <= kDirectMaxParentDimension;
    if (representable)
    {
        const bool expanded_k65 = shard_bytes >= 2048 &&
            IsExpandedDirectRepairCodec(codec);
#if LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE == 2
        const bool measured_small = output_count >= 5 &&
            IsSmallDirectRepairCodec(codec);
#else
        const bool measured_small = false;
#endif
        if (expanded_k65 || measured_small)
            return leopard2_internal::kDirectRepairExecutorSourceMajor;
    }
#else
    (void)shard_bytes;
#endif
    return leopard2_internal::kDirectRepairExecutorOutputMajor;
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
#ifdef LEO_HAS_FF8
    if (SelectDirectRepairExecutor(plan, shard_bytes) ==
            leopard2_internal::kDirectRepairExecutorSourceMajor)
        return ExecuteGF8SourceMajorDirectRepair(plan, shard_bytes,
            original, recovery, restored_original);
#endif
#ifdef LEO_HAS_FF8
    const bool tiled_gf8_tail =
        IsMeasuredEqualRoundedDirectRepairCodec(codec) &&
        (shard_bytes & (kScratchAlignment - 1U)) != 0;
    const size_t aligned_bytes = tiled_gf8_tail
        ? shard_bytes & ~(kScratchAlignment - 1U)
        : shard_bytes;
#else
    const size_t aligned_bytes = shard_bytes;
#endif
    for (uint32_t output_index = 0;
         output_index < plan->missing_original_count;
         ++output_index)
    {
        const size_t begin = plan->direct_term_offsets[output_index];
        const size_t end = plan->direct_term_offsets[output_index + 1];
        if (begin == end)
            return LEO2_INTERNAL_ERROR;
        void* output = restored_original[plan->missing_originals[output_index]];
        for (size_t term_index = begin;
             aligned_bytes != 0 && term_index < end;
             ++term_index)
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
                    memcpy(output, source, aligned_bytes);
                else if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
                    leopard::ff8::MultiplyBytes(ops, output, source,
                        static_cast<leopard::ff8::ffe_t>(term.multiplier_log), aligned_bytes);
#else
                    return LEO2_INTERNAL_ERROR;
#endif
                else
#ifdef LEO_HAS_FF16
                    leopard::ff16::MultiplyBytes(ops, output, source,
                        static_cast<leopard::ff16::ffe_t>(term.multiplier_log), aligned_bytes);
#else
                    return LEO2_INTERNAL_ERROR;
#endif
            }
            else if (term.multiplier_log == 0)
                XorArbitraryBytes(ops, output, source, aligned_bytes);
            else if (codec->field == LEO2_FIELD_GF8)
#ifdef LEO_HAS_FF8
                leopard::ff8::MultiplyAddBytes(ops, output, source,
                    static_cast<leopard::ff8::ffe_t>(term.multiplier_log), aligned_bytes);
#else
                return LEO2_INTERNAL_ERROR;
#endif
            else
#ifdef LEO_HAS_FF16
                leopard::ff16::MultiplyAddBytes(ops, output, source,
                    static_cast<leopard::ff16::ffe_t>(term.multiplier_log), aligned_bytes);
#else
                return LEO2_INTERNAL_ERROR;
#endif
        }
#ifdef LEO_HAS_FF8
        if (tiled_gf8_tail)
        {
            const size_t tail_bytes = shard_bytes - aligned_bytes;
            alignas(kScratchAlignment) uint8_t source_tile[kScratchAlignment];
            alignas(kScratchAlignment) uint8_t output_tile[kScratchAlignment];
            for (size_t term_index = begin; term_index < end; ++term_index)
            {
                const leo2_direct_repair_term& term =
                    plan->direct_terms[term_index];
                const bool parity =
                    (term.tagged_source & kDirectRecoveryTag) != 0;
                const uint32_t source_index =
                    term.tagged_source & ~kDirectRecoveryTag;
                const uint8_t* source = static_cast<const uint8_t*>(
                    parity ? recovery[source_index] : original[source_index]);
                if (!source)
                    return LEO2_INTERNAL_ERROR;
                memcpy(source_tile, source + aligned_bytes, tail_bytes);
                memset(source_tile + tail_bytes, 0,
                    sizeof(source_tile) - tail_bytes);
                if (term_index == begin)
                {
                    if (term.multiplier_log == 0)
                        memcpy(output_tile, source_tile, sizeof(output_tile));
                    else
                        leopard::ff8::MultiplyBytes(ops,
                            output_tile, source_tile,
                            static_cast<leopard::ff8::ffe_t>(
                                term.multiplier_log),
                            sizeof(output_tile));
                }
                else if (term.multiplier_log == 0)
                {
                    XorArbitraryBytes(
                        ops, output_tile, source_tile, sizeof(output_tile));
                }
                else
                {
                    leopard::ff8::MultiplyAddBytes(ops,
                        output_tile, source_tile,
                        static_cast<leopard::ff8::ffe_t>(
                            term.multiplier_log),
                        sizeof(output_tile));
                }
            }
            memcpy(static_cast<uint8_t*>(output) + aligned_bytes,
                output_tile, tail_bytes);
        }
#endif
    }
    return LEO2_SUCCESS;
}

/*
    Batch execution must not start until every writable range is disjoint from
    every other item's readable and writable ranges.  Otherwise a serial batch
    can consume bytes produced by an earlier item while a pooled batch races on
    the same bytes.  These helpers deliberately use each item's already
    required scratch range table as temporary preflight storage: batch setup
    remains allocation-free and does not impose a batch-size-dependent stack or
    heap bound.  For B items the pairwise sorted-list sweeps cost
    O(B^2 * (K + R)); this is deliberate setup work that prevents data races
    without silently lowering the documented UINT32_MAX item-count limit.

    Scratch is written only after the first pass proves that every supplied
    scratch span is disjoint from all batch metadata and shard spans.  Immutable
    input-input and metadata-input aliases remain supported.
*/
static bool SortedRangeListsOverlap(
    const AddressRange* a,
    size_t a_count,
    const AddressRange* b,
    size_t b_count)
{
    size_t a_i = 0;
    size_t b_i = 0;
    while (a_i < a_count && b_i < b_count)
    {
        if (RangesOverlap(a[a_i], b[b_i]))
            return true;
        if (a[a_i].end <= b[b_i].begin)
            ++a_i;
        else
            ++b_i;
    }
    return false;
}

static bool SortedRangesOverlapEachOther(
    const AddressRange* ranges,
    size_t count)
{
    for (size_t i = 1; i < count; ++i)
        if (RangesOverlap(ranges[i - 1], ranges[i]))
            return true;
    return false;
}

static bool BatchRangeLess(
    const BatchRangeRecord& a,
    const BatchRangeRecord& b)
{
    if (a.range.begin != b.range.begin)
        return a.range.begin < b.range.begin;
    if (a.range.end != b.range.end)
        return a.range.end < b.range.end;
    return a.writable < b.writable;
}

static void SiftBatchRangeHeap(
    BatchRangeRecord* ranges,
    size_t root,
    size_t count)
{
    for (;;)
    {
        if (root > (count - 1) / 2)
            return;
        const size_t left = root * 2 + 1;
        if (left >= count)
            return;
        size_t child = left;
        const size_t right = left + 1;
        if (right < count && BatchRangeLess(ranges[child], ranges[right]))
            child = right;
        if (!BatchRangeLess(ranges[root], ranges[child]))
            return;
        std::swap(ranges[root], ranges[child]);
        root = child;
    }
}

static void SortBatchRangesWithoutAllocation(
    BatchRangeRecord* ranges,
    size_t count)
{
    if (count < 2)
        return;
    for (size_t root = count / 2; root != 0; --root)
        SiftBatchRangeHeap(ranges, root - 1, count);
    for (size_t end = count; end > 1; --end)
    {
        std::swap(ranges[0], ranges[end - 1]);
        SiftBatchRangeHeap(ranges, 0, end - 1);
    }
}

static leo2_result ValidateSortedBatchRanges(
    const BatchRangeRecord* ranges,
    size_t count)
{
    uintptr_t furthest_end = 0;
    uintptr_t furthest_writable_end = 0;
    bool have_range = false;
    bool have_writable = false;
    for (size_t i = 0; i < count; ++i)
    {
        const BatchRangeRecord& current = ranges[i];
        if ((current.writable && have_range &&
                furthest_end > current.range.begin) ||
            (!current.writable && have_writable &&
                furthest_writable_end > current.range.begin))
            return LEO2_OVERLAP;
        if (!have_range || current.range.end > furthest_end)
            furthest_end = current.range.end;
        have_range = true;
        if (current.writable &&
            (!have_writable || current.range.end > furthest_writable_end))
        {
            furthest_writable_end = current.range.end;
            have_writable = true;
        }
    }
    return LEO2_SUCCESS;
}

static bool AppendBatchRange(
    BatchRangeRecord* ranges,
    size_t& count,
    size_t capacity,
    const AddressRange& range,
    bool writable)
{
    if (count >= capacity)
        return false;
    ranges[count].range = range;
    ranges[count].writable = writable ? 1 : 0;
    ++count;
    return true;
}

static leo2_result PrepareBatchPreflightScratch(
    void* scratch,
    size_t scratch_bytes,
    size_t required_bytes,
    AddressRange& scratch_range)
{
    scratch_range.begin = scratch_range.end = 0;
    if (scratch_bytes < required_bytes)
        return LEO2_SCRATCH_TOO_SMALL;
    if (required_bytes == 0)
        return LEO2_SUCCESS;
    if (!scratch)
        return LEO2_INVALID_ARGUMENT;
    if ((reinterpret_cast<uintptr_t>(scratch) &
            (kScratchAlignment - 1)) != 0)
        return LEO2_BAD_ALIGNMENT;
    if (!MakeRange(scratch, scratch_bytes, scratch_range))
        return LEO2_INVALID_ARGUMENT;
    return LEO2_SUCCESS;
}

static leo2_result EncodeBatchPreflightScratchBytes(
    const leo2_codec* codec,
    size_t item_count,
    size_t& scratch_bytes)
{
    scratch_bytes = 0;
    if (!codec)
        return LEO2_INVALID_ARGUMENT;
    if (item_count > 0xffffffffu)
        return LEO2_INVALID_ARGUMENT;
    size_t ranges_per_item = 0;
    if (!CheckedAdd(codec->original_count, codec->recovery_count,
            ranges_per_item) ||
        !CheckedAdd(ranges_per_item, 3, ranges_per_item))
        return LEO2_INVALID_ARGUMENT;
    return ComputeBatchPreflightScratchBytes(
        item_count, ranges_per_item, scratch_bytes);
}

static leo2_result DecodeBatchPreflightScratchBytes(
    const leo2_decode_plan* plan,
    size_t item_count,
    size_t& scratch_bytes)
{
    scratch_bytes = 0;
    if (!plan)
        return LEO2_INVALID_ARGUMENT;
    if (item_count > 0xffffffffu)
        return LEO2_INVALID_ARGUMENT;
    if (item_count < kScalableBatchPreflightMinItems || plan->no_op)
        return LEO2_SUCCESS;

    size_t recovery_input_count = 0;
    for (uint32_t i = 0; i < plan->codec->recovery_count; ++i)
        recovery_input_count += PlanRecoveryPresent(plan, i);
    size_t ranges_per_item = 0;
    if (!CheckedAdd(plan->codec->original_count,
            recovery_input_count, ranges_per_item) ||
        !CheckedAdd(ranges_per_item, 4, ranges_per_item))
        return LEO2_INVALID_ARGUMENT;
    return ComputeBatchPreflightScratchBytes(
        item_count, ranges_per_item, scratch_bytes);
}

static leo2_result ValidateEncodeBatchItemAddressability(
    const leo2_codec* codec,
    const leo2_encode_batch_item& item)
{
    ScratchLayout direct_layout;
    size_t rounded_bytes = 0;
    EncodeScratchGeometry geometry;
    const bool direct = UseSingleSideEncodeLayout(codec);
    leo2_result result = direct
        ? DirectEncodeLayout(
            codec, item.shard_bytes, direct_layout, rounded_bytes)
        : EncodeLayout(codec, item.shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;

    const ScratchLayout& layout = direct ? direct_layout : geometry.layout;
    /* DirectEncodeLayout and EncodeLayout both reserve exactly K+R range
       entries at range_offset zero.  Batch preflight uses no more than K
       inputs plus the requested subset of R outputs. */
    LEO_DEBUG_ASSERT(layout.range_offset == 0);
    AddressRange scratch_range;
    result = CheckScratch(item.scratch, item.scratch_bytes,
        layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    if (!item.original || !item.recovery)
        return LEO2_INVALID_ARGUMENT;

    AddressRange metadata_range;
    if (!MakeArrayRange(item.original, codec->original_count,
            sizeof(*item.original), metadata_range) ||
        !MakeArrayRange(item.recovery, codec->recovery_count,
            sizeof(*item.recovery), metadata_range))
        return LEO2_INVALID_ARGUMENT;

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        AddressRange range;
        if (!MakeRange(item.original[i], item.shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
        if (!HasValidSystematicPad(
                codec, item.original[i], item.shard_bytes))
            return LEO2_INVALID_ARGUMENT;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if (!item.recovery[i])
            continue;
        AddressRange range;
        if (!MakeRange(item.recovery[i], item.shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateDecodeBatchItemAddressability(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item& item,
    bool multi_item_batch)
{
    const leo2_codec* codec = plan->codec;
    ScratchLayout direct_layout;
    size_t rounded_bytes = 0;
    DecodeScratchGeometry geometry;
    const bool direct =
        plan->direct_repair || plan->direct_xor || plan->direct_copy;
    leo2_result result = direct
        ? DirectDecodeLayout(
            codec, item.shard_bytes, direct_layout, rounded_bytes)
        : DecodeLayout(codec, plan, item.shard_bytes,
            multi_item_batch, geometry);
    if (result != LEO2_SUCCESS)
        return result;

    const ScratchLayout& layout = direct ? direct_layout : geometry.layout;
    /* DirectDecodeLayout and DecodeLayout reserve 2K+R range entries at
       range_offset zero.  At most K+R received inputs plus K restored outputs
       are staged in that table. */
    LEO_DEBUG_ASSERT(layout.range_offset == 0);
    AddressRange scratch_range;
    result = CheckScratch(item.scratch, item.scratch_bytes,
        layout, scratch_range);
    if (result != LEO2_SUCCESS)
        return result;
    if (!item.original || !item.recovery || !item.restored_original)
        return LEO2_INVALID_ARGUMENT;

    AddressRange metadata_range;
    if (!MakeArrayRange(item.original, codec->original_count,
            sizeof(*item.original), metadata_range) ||
        !MakeArrayRange(item.recovery, codec->recovery_count,
            sizeof(*item.recovery), metadata_range) ||
        !MakeArrayRange(item.restored_original, codec->original_count,
            sizeof(*item.restored_original), metadata_range))
        return LEO2_INVALID_ARGUMENT;

    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if ((item.original[i] != NULL) !=
            PlanOriginalPresent(plan, i))
            return LEO2_INVALID_ARGUMENT;
        AddressRange range;
        if (item.original[i])
        {
            if (!MakeRange(item.original[i], item.shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
            if (!HasValidSystematicPad(
                    codec, item.original[i], item.shard_bytes))
                return LEO2_INVALID_ARGUMENT;
        }
        else if (!item.restored_original[i] ||
                 !MakeRange(item.restored_original[i],
                     item.shard_bytes, range))
            return LEO2_INVALID_ARGUMENT;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
    {
        if ((item.recovery[i] != NULL) !=
            PlanRecoveryPresent(plan, i))
            return LEO2_INVALID_ARGUMENT;
        if (item.recovery[i])
        {
            AddressRange range;
            if (!MakeRange(item.recovery[i], item.shard_bytes, range))
                return LEO2_INVALID_ARGUMENT;
        }
    }
    return LEO2_SUCCESS;
}

static bool EncodeScratchOverlapsItem(
    const leo2_codec* codec,
    const AddressRange& scratch_range,
    const leo2_encode_batch_item& item,
    bool compare_scratch)
{
    AddressRange range;
    if (MakeArrayRange(item.original, codec->original_count,
            sizeof(*item.original), range) &&
        RangesOverlap(scratch_range, range))
        return true;
    if (MakeArrayRange(item.recovery, codec->recovery_count,
            sizeof(*item.recovery), range) &&
        RangesOverlap(scratch_range, range))
        return true;
    if (compare_scratch &&
        MakeRange(item.scratch, item.scratch_bytes, range) &&
        RangesOverlap(scratch_range, range))
        return true;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (MakeRange(item.original[i], item.shard_bytes, range) &&
            RangesOverlap(scratch_range, range))
            return true;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (item.recovery[i] &&
            MakeRange(item.recovery[i], item.shard_bytes, range) &&
            RangesOverlap(scratch_range, range))
            return true;
    return false;
}

static bool DecodeScratchOverlapsItem(
    const leo2_decode_plan* plan,
    const AddressRange& scratch_range,
    const leo2_decode_batch_item& item,
    bool compare_scratch)
{
    const leo2_codec* codec = plan->codec;
    AddressRange range;
    if (MakeArrayRange(item.original, codec->original_count,
            sizeof(*item.original), range) &&
        RangesOverlap(scratch_range, range))
        return true;
    if (MakeArrayRange(item.recovery, codec->recovery_count,
            sizeof(*item.recovery), range) &&
        RangesOverlap(scratch_range, range))
        return true;
    if (MakeArrayRange(item.restored_original, codec->original_count,
            sizeof(*item.restored_original), range) &&
        RangesOverlap(scratch_range, range))
        return true;
    if (compare_scratch &&
        MakeRange(item.scratch, item.scratch_bytes, range) &&
        RangesOverlap(scratch_range, range))
        return true;
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        const void* const data = item.original[i]
            ? item.original[i] : item.restored_original[i];
        if (MakeRange(data, item.shard_bytes, range) &&
            RangesOverlap(scratch_range, range))
            return true;
    }
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (item.recovery[i] &&
            MakeRange(item.recovery[i], item.shard_bytes, range) &&
            RangesOverlap(scratch_range, range))
            return true;
    return false;
}

static size_t PrepareEncodeBatchRanges(
    const leo2_codec* codec,
    const leo2_encode_batch_item& item)
{
    AddressRange* const inputs =
        static_cast<AddressRange*>(item.scratch);
    for (uint32_t i = 0; i < codec->original_count; ++i)
        MakeRange(item.original[i], item.shard_bytes, inputs[i]);
    std::sort(inputs, inputs + codec->original_count, RangeLess);

    AddressRange* const outputs = inputs + codec->original_count;
    size_t output_count = 0;
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (item.recovery[i])
            MakeRange(item.recovery[i], item.shard_bytes,
                outputs[output_count++]);
    LEO_DEBUG_ASSERT(output_count <= codec->recovery_count);
    std::sort(outputs, outputs + output_count, RangeLess);
    return output_count;
}

static size_t DecodeBatchInputCount(const leo2_decode_plan* plan)
{
    size_t input_count = 0;
    for (uint32_t i = 0; i < plan->codec->original_count; ++i)
        input_count += PlanOriginalPresent(plan, i);
    for (uint32_t i = 0; i < plan->codec->recovery_count; ++i)
        input_count += PlanRecoveryPresent(plan, i);
    return input_count;
}

static void PrepareDecodeBatchRanges(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item& item,
    size_t input_count)
{
    const leo2_codec* codec = plan->codec;
    AddressRange* const inputs =
        static_cast<AddressRange*>(item.scratch);
    size_t input_i = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (item.original[i])
            MakeRange(item.original[i], item.shard_bytes, inputs[input_i++]);
    for (uint32_t i = 0; i < codec->recovery_count; ++i)
        if (item.recovery[i])
            MakeRange(item.recovery[i], item.shard_bytes, inputs[input_i++]);
    LEO_DEBUG_ASSERT(input_i == input_count);
    LEO_DEBUG_ASSERT(input_count <=
        static_cast<size_t>(codec->original_count) + codec->recovery_count);
    std::sort(inputs, inputs + input_count, RangeLess);

    AddressRange* const outputs = inputs + input_count;
    size_t output_i = 0;
    for (uint32_t i = 0; i < codec->original_count; ++i)
        if (!item.original[i])
            MakeRange(item.restored_original[i], item.shard_bytes,
                outputs[output_i++]);
    LEO_DEBUG_ASSERT(output_i == plan->missing_original_count);
    LEO_DEBUG_ASSERT(input_count + output_i <=
        static_cast<size_t>(codec->original_count) * 2 +
            codec->recovery_count);
    std::sort(outputs, outputs + output_i, RangeLess);
}

static leo2_result ValidateEncodeBatchOutputMetadataAliases(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range)
{
    /* This pass must remain read-only.  ValidateEncodeBuffers may use scratch
       as a range table, so reject every cross-item output/metadata alias before
       invoking even the first ordinary per-item validator. */
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        for (uint32_t output_i = 0;
             output_i < codec->recovery_count;
             ++output_i)
        {
            if (!items[item_i].recovery[output_i])
                continue;
            AddressRange output_range;
            if (!MakeRange(items[item_i].recovery[output_i],
                    items[item_i].shard_bytes, output_range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(output_range, item_range))
                return LEO2_OVERLAP;
            for (size_t other_i = 0; other_i < item_count; ++other_i)
            {
                AddressRange metadata_range;
                if ((MakeArrayRange(items[other_i].original,
                         codec->original_count,
                         sizeof(*items[other_i].original), metadata_range) &&
                     RangesOverlap(output_range, metadata_range)) ||
                    (MakeArrayRange(items[other_i].recovery,
                         codec->recovery_count,
                         sizeof(*items[other_i].recovery), metadata_range) &&
                     RangesOverlap(output_range, metadata_range)))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateDecodeBatchOutputMetadataAliases(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range)
{
    const leo2_codec* codec = plan->codec;
    /* As above, preserve the public guarantee that metadata-overlap rejection
       happens before any scratch byte is modified. */
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        for (uint32_t output_i = 0;
             output_i < codec->original_count;
             ++output_i)
        {
            if (items[item_i].original[output_i])
                continue;
            AddressRange output_range;
            if (!MakeRange(items[item_i].restored_original[output_i],
                    items[item_i].shard_bytes, output_range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(output_range, item_range))
                return LEO2_OVERLAP;
            for (size_t other_i = 0; other_i < item_count; ++other_i)
            {
                AddressRange metadata_range;
                if ((MakeArrayRange(items[other_i].original,
                         codec->original_count,
                         sizeof(*items[other_i].original), metadata_range) &&
                     RangesOverlap(output_range, metadata_range)) ||
                    (MakeArrayRange(items[other_i].recovery,
                         codec->recovery_count,
                         sizeof(*items[other_i].recovery), metadata_range) &&
                     RangesOverlap(output_range, metadata_range)) ||
                    (MakeArrayRange(items[other_i].restored_original,
                         codec->original_count,
                         sizeof(*items[other_i].restored_original),
                         metadata_range) &&
                     RangesOverlap(output_range, metadata_range)))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateK1R1EncodeBatch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    size_t item_bytes)
{
    LEO_DEBUG_ASSERT(IsLegacyHighK1R1Codec(codec));
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_result result =
            ValidateEncodeBatchItemAddressability(codec, items[item_i]);
        if (result != LEO2_SUCCESS)
            return result;
    }

    leo2_result result = ValidateEncodeBatchOutputMetadataAliases(
        codec, items, item_count, item_range);
    if (result != LEO2_SUCCESS)
        return result;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange scratch_range;
        result = CheckOptionalScratch(items[item_i].scratch,
            items[item_i].scratch_bytes, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        result = ValidateK1R1EncodeBuffers(codec,
            items[item_i].shard_bytes, items[item_i].original,
            items[item_i].recovery, scratch_range, items, item_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        if (items[item_i].scratch_bytes != 0)
        {
            if (RangesOverlap(scratch_range, item_range))
                return LEO2_OVERLAP;
            for (size_t other_i = 0; other_i < item_count; ++other_i)
                if (EncodeScratchOverlapsItem(codec, scratch_range,
                        items[other_i], item_i != other_i))
                    return LEO2_OVERLAP;
        }
    }
    for (size_t output_i = 0; output_i < item_count; ++output_i)
    {
        if (!items[output_i].recovery[0])
            continue;
        AddressRange output_range;
        if (!MakeRange(items[output_i].recovery[0],
                items[output_i].shard_bytes, output_range))
            return LEO2_INVALID_ARGUMENT;
        for (size_t other_i = 0; other_i < item_count; ++other_i)
        {
            AddressRange input_range;
            if (!MakeRange(items[other_i].original[0],
                    items[other_i].shard_bytes, input_range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(output_range, input_range))
                return LEO2_OVERLAP;
            if (other_i > output_i && items[other_i].recovery[0])
            {
                AddressRange other_output_range;
                if (!MakeRange(items[other_i].recovery[0],
                        items[other_i].shard_bytes, other_output_range))
                    return LEO2_INVALID_ARGUMENT;
                if (RangesOverlap(output_range, other_output_range))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateK1R1DecodeBatch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range)
{
    LEO_DEBUG_ASSERT(plan && IsLegacyHighK1R1Codec(plan->codec) &&
        PlanHasCompactDirectXor(plan) && plan->direct_xor_original == 0);
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_result result = ValidateDecodeBatchItemAddressability(
            plan, items[item_i], true);
        if (result != LEO2_SUCCESS)
            return result;
    }

    leo2_result result = ValidateDecodeBatchOutputMetadataAliases(
        plan, items, item_count, item_range);
    if (result != LEO2_SUCCESS)
        return result;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange scratch_range;
        result = CheckOptionalScratch(items[item_i].scratch,
            items[item_i].scratch_bytes, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        result = ValidateK1R1DecodeBuffers(items[item_i].shard_bytes,
            scratch_range, &item_range, 1, items[item_i].original,
            items[item_i].recovery, items[item_i].restored_original);
        if (result != LEO2_SUCCESS)
            return result;
        if (items[item_i].scratch_bytes != 0)
        {
            if (RangesOverlap(scratch_range, item_range))
                return LEO2_OVERLAP;
            for (size_t other_i = 0; other_i < item_count; ++other_i)
                if (DecodeScratchOverlapsItem(plan, scratch_range,
                        items[other_i], item_i != other_i))
                    return LEO2_OVERLAP;
        }
    }
    for (size_t output_i = 0; output_i < item_count; ++output_i)
    {
        AddressRange output_range;
        if (!MakeRange(items[output_i].restored_original[0],
                items[output_i].shard_bytes, output_range))
            return LEO2_INVALID_ARGUMENT;
        for (size_t other_i = 0; other_i < item_count; ++other_i)
        {
            AddressRange input_range;
            if (!MakeRange(items[other_i].recovery[0],
                    items[other_i].shard_bytes, input_range))
                return LEO2_INVALID_ARGUMENT;
            if (RangesOverlap(output_range, input_range))
                return LEO2_OVERLAP;
            if (other_i > output_i)
            {
                AddressRange other_output_range;
                if (!MakeRange(items[other_i].restored_original[0],
                        items[other_i].shard_bytes, other_output_range))
                    return LEO2_INVALID_ARGUMENT;
                if (RangesOverlap(output_range, other_output_range))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateEncodeBatchAliases(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    size_t item_bytes)
{
    leo2_result result = ValidateEncodeBatchOutputMetadataAliases(
        codec, items, item_count, item_range);
    if (result != LEO2_SUCCESS)
        return result;

    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange scratch_range = { 0, 0 };
        if (!MakeRange(items[item_i].scratch,
                items[item_i].scratch_bytes, scratch_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(scratch_range, item_range))
            return LEO2_OVERLAP;
        for (size_t other_i = 0; other_i < item_count; ++other_i)
            if (EncodeScratchOverlapsItem(codec, scratch_range,
                    items[other_i], item_i != other_i))
                return LEO2_OVERLAP;
    }

    /* Scratch is globally safe now, so the ordinary per-item validator may
       use it for range sorting without corrupting another item. */
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_encode_batch_item& item = items[item_i];
        result = ValidateEncodeBuffers(codec,
            item.shard_bytes, item.original, item.recovery, item.scratch,
            item.scratch_bytes, items, item_bytes);
        if (result != LEO2_SUCCESS)
            return result;
    }

    for (size_t item_i = 0; item_i < item_count; ++item_i)
        PrepareEncodeBatchRanges(codec, items[item_i]);

    const size_t input_count = codec->original_count;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange* const item_inputs =
            static_cast<AddressRange*>(items[item_i].scratch);
        AddressRange* const item_outputs = item_inputs + input_count;
        size_t item_output_count = 0;
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
            item_output_count += items[item_i].recovery[i] != NULL;
        if (SortedRangesOverlapEachOther(
                item_outputs, item_output_count))
            return LEO2_OVERLAP;

        for (size_t other_i = 0; other_i < item_count; ++other_i)
        {
            AddressRange* const other_inputs =
                static_cast<AddressRange*>(items[other_i].scratch);
            if (SortedRangeListsOverlap(item_outputs, item_output_count,
                    other_inputs, input_count))
                return LEO2_OVERLAP;
            if (other_i > item_i)
            {
                AddressRange* const other_outputs =
                    other_inputs + input_count;
                size_t other_output_count = 0;
                for (uint32_t i = 0; i < codec->recovery_count; ++i)
                    other_output_count +=
                        items[other_i].recovery[i] != NULL;
                if (SortedRangeListsOverlap(item_outputs,
                        item_output_count, other_outputs,
                        other_output_count))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result ValidateDecodeBatchAliases(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    size_t item_bytes)
{
    leo2_result result = ValidateDecodeBatchOutputMetadataAliases(
        plan, items, item_count, item_range);
    if (result != LEO2_SUCCESS)
        return result;

    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange scratch_range = { 0, 0 };
        if (!MakeRange(items[item_i].scratch,
                items[item_i].scratch_bytes, scratch_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(scratch_range, item_range))
            return LEO2_OVERLAP;
        for (size_t other_i = 0; other_i < item_count; ++other_i)
            if (DecodeScratchOverlapsItem(plan, scratch_range,
                    items[other_i], item_i != other_i))
                return LEO2_OVERLAP;
    }

    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_decode_batch_item& item = items[item_i];
        const ProtectedMetadataSpan item_metadata = { items, item_bytes };
        result = ValidateDecodeBuffers(plan,
            item.shard_bytes, item.original, item.recovery,
            item.restored_original, item.scratch, item.scratch_bytes,
            &item_metadata, 1, NULL);
        if (result != LEO2_SUCCESS)
            return result;
    }

    const size_t input_count = DecodeBatchInputCount(plan);
    const size_t output_count = plan->missing_original_count;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
        PrepareDecodeBatchRanges(
            plan, items[item_i], input_count);

    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        AddressRange* const item_inputs =
            static_cast<AddressRange*>(items[item_i].scratch);
        AddressRange* const item_outputs = item_inputs + input_count;
        if (SortedRangesOverlapEachOther(item_outputs, output_count))
            return LEO2_OVERLAP;

        for (size_t other_i = 0; other_i < item_count; ++other_i)
        {
            AddressRange* const other_inputs =
                static_cast<AddressRange*>(items[other_i].scratch);
            if (SortedRangeListsOverlap(item_outputs, output_count,
                    other_inputs, input_count))
                return LEO2_OVERLAP;
            if (other_i > item_i)
            {
                AddressRange* const other_outputs =
                    other_inputs + input_count;
                if (SortedRangeListsOverlap(item_outputs, output_count,
                        other_outputs, output_count))
                    return LEO2_OVERLAP;
            }
        }
    }
    return LEO2_SUCCESS;
}

static leo2_result PreflightEncodeBatch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    size_t item_bytes)
{
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result =
            ValidateEncodeBatchItemAddressability(codec, items[i]);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return ValidateEncodeBatchAliases(
        codec, items, item_count, item_range, item_bytes);
}

static leo2_result PreflightDecodeBatch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    size_t item_bytes)
{
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result =
            ValidateDecodeBatchItemAddressability(
                plan, items[i], item_count > 1);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return ValidateDecodeBatchAliases(
        plan, items, item_count, item_range, item_bytes);
}

static leo2_result PreflightEncodeBatchScalable(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    void* preflight_scratch,
    size_t preflight_scratch_bytes)
{
    size_t required_bytes = 0;
    leo2_result result = EncodeBatchPreflightScratchBytes(
        codec, item_count, required_bytes);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange preflight_range;
    result = PrepareBatchPreflightScratch(preflight_scratch,
        preflight_scratch_bytes, required_bytes, preflight_range);
    if (result != LEO2_SUCCESS)
        return result;

    /* Finish every potentially failing read-only check before writing the
       caller's batch-wide workspace. */
    for (size_t i = 0; i < item_count; ++i)
    {
        result = ValidateEncodeBatchItemAddressability(codec, items[i]);
        if (result != LEO2_SUCCESS)
            return result;
    }
    if (RangesOverlap(preflight_range, item_range))
        return LEO2_OVERLAP;
    for (size_t i = 0; i < item_count; ++i)
        if (EncodeScratchOverlapsItem(
                codec, preflight_range, items[i], true))
            return LEO2_OVERLAP;

    BatchRangeRecord* const ranges =
        static_cast<BatchRangeRecord*>(preflight_scratch);
    const size_t capacity = required_bytes / sizeof(*ranges);
    size_t count = 0;
    const auto append = [&](const AddressRange& range, bool writable) {
        return AppendBatchRange(
            ranges, count, capacity, range, writable);
    };
    if (!append(item_range, false))
        return LEO2_INTERNAL_ERROR;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_encode_batch_item& item = items[item_i];
        AddressRange range;
        if (!MakeArrayRange(item.original, codec->original_count,
                sizeof(*item.original), range))
            return LEO2_INTERNAL_ERROR;
        if (!append(range, false))
            return LEO2_INTERNAL_ERROR;
        if (!MakeArrayRange(item.recovery, codec->recovery_count,
                sizeof(*item.recovery), range))
            return LEO2_INTERNAL_ERROR;
        if (!append(range, false))
            return LEO2_INTERNAL_ERROR;
        if (item.scratch_bytes != 0)
        {
            if (!MakeRange(item.scratch, item.scratch_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, true))
                return LEO2_INTERNAL_ERROR;
        }
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (!MakeRange(item.original[i], item.shard_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, false))
                return LEO2_INTERNAL_ERROR;
        }
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
        {
            if (!item.recovery[i])
                continue;
            if (!MakeRange(item.recovery[i], item.shard_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, true))
                return LEO2_INTERNAL_ERROR;
        }
    }
    if (count > capacity)
        return LEO2_INTERNAL_ERROR;
    SortBatchRangesWithoutAllocation(ranges, count);
    return ValidateSortedBatchRanges(ranges, count);
}

static leo2_result PreflightDecodeBatchScalable(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    const AddressRange& item_range,
    void* preflight_scratch,
    size_t preflight_scratch_bytes)
{
    size_t required_bytes = 0;
    leo2_result result = DecodeBatchPreflightScratchBytes(
        plan, item_count, required_bytes);
    if (result != LEO2_SUCCESS)
        return result;
    AddressRange preflight_range;
    result = PrepareBatchPreflightScratch(preflight_scratch,
        preflight_scratch_bytes, required_bytes, preflight_range);
    if (result != LEO2_SUCCESS)
        return result;

    for (size_t i = 0; i < item_count; ++i)
    {
        result = ValidateDecodeBatchItemAddressability(
            plan, items[i], item_count > 1);
        if (result != LEO2_SUCCESS)
            return result;
    }
    if (RangesOverlap(preflight_range, item_range))
        return LEO2_OVERLAP;
    for (size_t i = 0; i < item_count; ++i)
        if (DecodeScratchOverlapsItem(
                plan, preflight_range, items[i], true))
            return LEO2_OVERLAP;

    const leo2_codec* const codec = plan->codec;
    BatchRangeRecord* const ranges =
        static_cast<BatchRangeRecord*>(preflight_scratch);
    const size_t capacity = required_bytes / sizeof(*ranges);
    size_t count = 0;
    const auto append = [&](const AddressRange& range, bool writable) {
        return AppendBatchRange(
            ranges, count, capacity, range, writable);
    };
    if (!append(item_range, false))
        return LEO2_INTERNAL_ERROR;
    for (size_t item_i = 0; item_i < item_count; ++item_i)
    {
        const leo2_decode_batch_item& item = items[item_i];
        AddressRange range;
        if (!MakeArrayRange(item.original, codec->original_count,
                sizeof(*item.original), range))
            return LEO2_INTERNAL_ERROR;
        if (!append(range, false))
            return LEO2_INTERNAL_ERROR;
        if (!MakeArrayRange(item.recovery, codec->recovery_count,
                sizeof(*item.recovery), range))
            return LEO2_INTERNAL_ERROR;
        if (!append(range, false))
            return LEO2_INTERNAL_ERROR;
        if (!MakeArrayRange(item.restored_original, codec->original_count,
                sizeof(*item.restored_original), range))
            return LEO2_INTERNAL_ERROR;
        if (!append(range, false))
            return LEO2_INTERNAL_ERROR;
        if (item.scratch_bytes != 0)
        {
            if (!MakeRange(item.scratch, item.scratch_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, true))
                return LEO2_INTERNAL_ERROR;
        }
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            const bool writable = item.original[i] == NULL;
            const void* const shard = writable
                ? item.restored_original[i] : item.original[i];
            if (!MakeRange(shard, item.shard_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, writable))
                return LEO2_INTERNAL_ERROR;
        }
        for (uint32_t i = 0; i < codec->recovery_count; ++i)
        {
            if (!item.recovery[i])
                continue;
            if (!MakeRange(item.recovery[i], item.shard_bytes, range))
                return LEO2_INTERNAL_ERROR;
            if (!append(range, false))
                return LEO2_INTERNAL_ERROR;
        }
    }
    if (count > capacity)
        return LEO2_INTERNAL_ERROR;
    SortBatchRangesWithoutAllocation(ranges, count);
    return ValidateSortedBatchRanges(ranges, count);
}

struct EncodeBatchTaskContext
{
    const leo2_codec* codec;
    const leo2_encode_batch_item* items;
    size_t item_bytes;
    AddressRange item_range;
};

static leo2_result PreflightEncodeBatchTask(void* context)
{
    const EncodeBatchTaskContext* batch =
        static_cast<const EncodeBatchTaskContext*>(context);
    return PreflightEncodeBatch(batch->codec, batch->items,
        batch->item_bytes / sizeof(*batch->items), batch->item_range,
        batch->item_bytes);
}

static leo2_result PreflightK1R1EncodeBatchTask(void* context)
{
    const EncodeBatchTaskContext* batch =
        static_cast<const EncodeBatchTaskContext*>(context);
    return ValidateK1R1EncodeBatch(batch->codec, batch->items,
        batch->item_bytes / sizeof(*batch->items), batch->item_range,
        batch->item_bytes);
}

static leo2_result RunEncodeBatchItem(void* context, size_t index);

struct ScalableEncodeBatchPreflightContext
{
    EncodeBatchTaskContext batch;
    void* scratch;
    size_t scratch_bytes;
};

static leo2_result PreflightEncodeBatchScalableTask(void* context)
{
    const ScalableEncodeBatchPreflightContext* scalable =
        static_cast<const ScalableEncodeBatchPreflightContext*>(context);
    const EncodeBatchTaskContext& batch = scalable->batch;
    return PreflightEncodeBatchScalable(batch.codec, batch.items,
        batch.item_bytes / sizeof(*batch.items), batch.item_range,
        scalable->scratch, scalable->scratch_bytes);
}

static leo2_result RunScalableEncodeBatchItem(void* context, size_t index)
{
    ScalableEncodeBatchPreflightContext* scalable =
        static_cast<ScalableEncodeBatchPreflightContext*>(context);
    return RunEncodeBatchItem(&scalable->batch, index);
}

static leo2_result RunEncodeBatchItem(void* context, size_t index)
{
    const EncodeBatchTaskContext* batch = static_cast<const EncodeBatchTaskContext*>(context);
    const leo2_encode_batch_item& item = batch->items[index];
    return EncodeInternal(
        batch->codec, item.shard_bytes, item.original, item.recovery,
        item.scratch, item.scratch_bytes, batch->items, batch->item_bytes,
        true);
}

struct DecodeBatchTaskContext
{
    const leo2_decode_plan* plan;
    const leo2_decode_batch_item* items;
    size_t item_bytes;
    AddressRange item_range;
    bool multi_item_batch;
};

static leo2_result PreflightDecodeBatchTask(void* context)
{
    const DecodeBatchTaskContext* batch =
        static_cast<const DecodeBatchTaskContext*>(context);
    return PreflightDecodeBatch(batch->plan, batch->items,
        batch->item_bytes / sizeof(*batch->items), batch->item_range,
        batch->item_bytes);
}

static leo2_result PreflightK1R1DecodeBatchTask(void* context)
{
    const DecodeBatchTaskContext* batch =
        static_cast<const DecodeBatchTaskContext*>(context);
    return ValidateK1R1DecodeBatch(batch->plan, batch->items,
        batch->item_bytes / sizeof(*batch->items), batch->item_range);
}

static leo2_result RunDecodeBatchItem(void* context, size_t index);

struct ScalableDecodeBatchPreflightContext
{
    DecodeBatchTaskContext batch;
    void* scratch;
    size_t scratch_bytes;
};

static leo2_result PreflightDecodeBatchScalableTask(void* context)
{
    const ScalableDecodeBatchPreflightContext* scalable =
        static_cast<const ScalableDecodeBatchPreflightContext*>(context);
    const DecodeBatchTaskContext& batch = scalable->batch;
    return PreflightDecodeBatchScalable(batch.plan, batch.items,
        batch.item_bytes / sizeof(*batch.items), batch.item_range,
        scalable->scratch, scalable->scratch_bytes);
}

static leo2_result RunScalableDecodeBatchItem(void* context, size_t index)
{
    ScalableDecodeBatchPreflightContext* scalable =
        static_cast<ScalableDecodeBatchPreflightContext*>(context);
    return RunDecodeBatchItem(&scalable->batch, index);
}

static leo2_result RunDecodeBatchItem(void* context, size_t index)
{
    const DecodeBatchTaskContext* batch = static_cast<const DecodeBatchTaskContext*>(context);
    const leo2_decode_batch_item& item = batch->items[index];
    return DecodePlanExecuteInternal(
        batch->plan, item.shard_bytes, item.original, item.recovery,
        item.restored_original, item.scratch, item.scratch_bytes,
        batch->multi_item_batch, NULL, 0, true);
}

} // namespace

namespace leopard2_internal {

static void FillTerminalDecodePathInfo(
    DecodePath path,
    DecodePathRule rule,
    uint64_t shard_bytes,
    bool multi_item_batch,
    DecodePathInfo& info)
{
    info.path = path;
    info.rule = rule;
    info.direct_executor = leopard2_internal::kDirectRepairExecutorNone;
    info.matching_auto_rules = 0;
    info.required_work_slots = 0;
    info.aligned_prefix_bytes = 0;
    info.tail_bytes = 0;
    info.rounded_shard_bytes = 0;
    info.multi_item_batch = multi_item_batch;
    size_t rounded = 0;
    if (shard_bytes != 0 && RoundShardBytes(shard_bytes, rounded))
    {
        info.tail_bytes = static_cast<size_t>(shard_bytes) & 63u;
        info.aligned_prefix_bytes =
            static_cast<size_t>(shard_bytes) - info.tail_bytes;
        info.rounded_shard_bytes = rounded;
    }
}

leo2_result GetDecodePlanPathInfo(
    const leo2_decode_plan* plan,
    uint64_t shard_bytes,
    bool multi_item_batch,
    DecodePathInfo* info_out)
{
    if (!IsRepresentableOutput(info_out))
        return LEO2_INVALID_ARGUMENT;
    FillTerminalDecodePathInfo(kDecodePathNoOp, kDecodeRuleNoOp,
        0, multi_item_batch, *info_out);
    if (!plan)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
    {
        FillTerminalDecodePathInfo(kDecodePathNoOp, kDecodeRuleNoOp,
            shard_bytes, multi_item_batch, *info_out);
        return LEO2_SUCCESS;
    }
    if (shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (plan->direct_repair || plan->direct_xor || plan->direct_copy)
    {
        ScratchLayout layout;
        size_t rounded = 0;
        const leo2_result result = DirectDecodeLayout(
            plan->codec, shard_bytes, layout, rounded);
        if (result != LEO2_SUCCESS)
            return result;
        FillTerminalDecodePathInfo(kDecodePathDirect, kDecodeRuleDirect,
            shard_bytes, multi_item_batch, *info_out);
        if (plan->direct_repair)
            info_out->direct_executor =
                SelectDirectRepairExecutor(
                    plan, static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }
    DecodeScratchGeometry geometry;
    const leo2_result result = DecodeLayout(
        plan->codec, plan, shard_bytes, multi_item_batch, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *info_out = geometry.selection;
    return LEO2_SUCCESS;
}

leo2_result GetDecodeCodecScratchPathInfo(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    DecodePathInfo* info_out)
{
    if (!IsRepresentableOutput(info_out))
        return LEO2_INVALID_ARGUMENT;
    FillTerminalDecodePathInfo(kDecodePathNoOp, kDecodeRuleNoOp,
        0, false, *info_out);
    if (IsLegacyHighK1R1Codec(codec))
    {
        ScratchLayout layout;
        size_t rounded_bytes = 0;
        const leo2_result result = DirectDecodeLayout(
            codec, shard_bytes, layout, rounded_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        FillTerminalDecodePathInfo(kDecodePathDirect, kDecodeRuleDirect,
            shard_bytes, false, *info_out);
        return LEO2_SUCCESS;
    }
    DecodeScratchGeometry geometry;
    const leo2_result result = DecodeLayout(
        codec, NULL, shard_bytes, false, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    *info_out = geometry.selection;
    return LEO2_SUCCESS;
}

} // namespace leopard2_internal

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

    AddressRange output_range;
    if (!MakeArrayRange(
            context_out, 1, sizeof(*context_out), output_range))
        return LEO2_INVALID_ARGUMENT;
    size_t options_size = 0;
    if (options)
    {
        AddressRange options_prefix_range;
        if (!MakeRange(options, sizeof(options->struct_size),
                options_prefix_range))
            return LEO2_INVALID_ARGUMENT;
        options_size = options->struct_size;
        const size_t protected_size = std::max(
            options_size, sizeof(options->struct_size));
        AddressRange options_range;
        if (!MakeRange(options, static_cast<uint64_t>(protected_size),
                options_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(options_range, output_range))
            return LEO2_OVERLAP;
    }
    *context_out = NULL;
    if (options && (options_size < sizeof(leo2_context_options) ||
                    options->reserved != 0))
        return LEO2_INVALID_ARGUMENT;
    const uint32_t requested_raw = options
        ? options->backend : static_cast<uint32_t>(LEO2_BACKEND_AUTO);
    if (requested_raw > static_cast<uint32_t>(LEO2_BACKEND_AVX512))
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
    context->baseline_ops = ops;
    context->auto_requested = requested == LEO2_BACKEND_AUTO;
    context->auto_avx512_encode_host = context->auto_requested &&
        leopard::backend::IsCalibratedAutoAVX512EncodeHost();
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
LEO2_EXPORT leo2_result leo2_test_balanced_execution_tiles(
    uint64_t aligned_bytes,
    uint64_t requested_tile_bytes,
    size_t* execution_tile_count_out,
    size_t* maximum_pass_bytes_out)
{
    if (!execution_tile_count_out || !maximum_pass_bytes_out ||
        aligned_bytes > std::numeric_limits<size_t>::max() ||
        requested_tile_bytes > std::numeric_limits<size_t>::max())
        return LEO2_INVALID_ARGUMENT;
    size_t execution_tile_count = 0;
    size_t maximum_pass_bytes = 0;
    if (!ComputeBalancedExecutionTiles(
            static_cast<size_t>(aligned_bytes),
            static_cast<size_t>(requested_tile_bytes),
            execution_tile_count, maximum_pass_bytes))
        return LEO2_INVALID_ARGUMENT;
    *execution_tile_count_out = execution_tile_count;
    *maximum_pass_bytes_out = maximum_pass_bytes;
    return LEO2_SUCCESS;
}

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

LEO2_EXPORT void leo2_test_reset_generic_reveal_counts(void)
{
    g_test_generic_direct_reveal_shards.store(0, std::memory_order_release);
}

LEO2_EXPORT uint64_t leo2_test_generic_direct_reveal_shards(void)
{
    return g_test_generic_direct_reveal_shards.load(std::memory_order_acquire);
}

LEO2_EXPORT void leo2_test_reset_low_reveal_counts(void)
{
    g_test_low_direct_reveal_shards.store(0, std::memory_order_release);
    g_test_low_scratch_reveal_shards.store(0, std::memory_order_release);
}

LEO2_EXPORT uint64_t leo2_test_low_direct_reveal_shards(void)
{
    return g_test_low_direct_reveal_shards.load(std::memory_order_acquire);
}

LEO2_EXPORT uint64_t leo2_test_low_scratch_reveal_shards(void)
{
    return g_test_low_scratch_reveal_shards.load(std::memory_order_acquire);
}

LEO2_EXPORT void leo2_test_reset_high_reveal_counts(void)
{
    g_test_high_direct_reveal_shards.store(0, std::memory_order_release);
    g_test_high_materialized_direct_reveal_shards.store(
        0, std::memory_order_release);
    g_test_high_scratch_reveal_shards.store(0, std::memory_order_release);
}

LEO2_EXPORT uint64_t leo2_test_high_direct_reveal_shards(void)
{
    return g_test_high_direct_reveal_shards.load(std::memory_order_acquire);
}

LEO2_EXPORT uint64_t
leo2_test_high_materialized_direct_reveal_shards(void)
{
    return g_test_high_materialized_direct_reveal_shards.load(
        std::memory_order_acquire);
}

LEO2_EXPORT uint64_t leo2_test_high_scratch_reveal_shards(void)
{
    return g_test_high_scratch_reveal_shards.load(std::memory_order_acquire);
}
#endif

LEO2_EXPORT leo2_result leo2_codec_create(
    leo2_context* context,
    uint32_t original_count,
    uint32_t recovery_count,
    int profile_value,
    int field_value,
    const leo2_codec_options* options,
    leo2_codec** codec_out)
{
    if (!codec_out)
        return LEO2_INVALID_ARGUMENT;

    AddressRange output_range;
    if (!MakeArrayRange(codec_out, 1, sizeof(*codec_out), output_range))
        return LEO2_INVALID_ARGUMENT;
    size_t options_size = 0;
    if (options)
    {
        AddressRange options_prefix_range;
        if (!MakeRange(options, sizeof(options->struct_size),
                options_prefix_range))
            return LEO2_INVALID_ARGUMENT;
        options_size = options->struct_size;
        const size_t protected_size = std::max(
            options_size, sizeof(options->struct_size));
        AddressRange options_range;
        if (!MakeRange(options, static_cast<uint64_t>(protected_size),
                options_range))
            return LEO2_INVALID_ARGUMENT;
        if (RangesOverlap(options_range, output_range))
            return LEO2_OVERLAP;
    }
    *codec_out = NULL;
    if (!context || original_count == 0 || recovery_count == 0)
        return LEO2_INVALID_ARGUMENT;
    leo2_shard_layout shard_layout = LEO2_SHARD_LAYOUT_NATIVE_V1;
    if (options)
    {
        const size_t version1_size = offsetof(leo2_codec_options, shard_layout);
        const size_t layout_field_end =
            version1_size + sizeof(options->shard_layout);
        // Validate the readable v1 prefix before accessing flags or reserved.
        // This permits callers to pass an intentionally short versioned
        // prefix and receive INVALID_ARGUMENT without an out-of-range read.
        if (options_size < version1_size ||
            (options_size > version1_size &&
             options_size < layout_field_end))
            return LEO2_INVALID_ARGUMENT;

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
        if (options->reserved != 0 ||
            (options->flags & ~supported_flags) != 0 ||
            algorithm_flags == (LEO2_CODEC_FORCE_GENERIC_DECODE |
                                LEO2_CODEC_FORCE_SPECIALIZED_DECODE) ||
            workspace_flags == (LEO2_CODEC_FORCE_TILED_DECODE |
                                LEO2_CODEC_FORCE_MATERIALIZED_DECODE) ||
            ((algorithm_flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0 &&
             workspace_flags != 0))
            return LEO2_INVALID_ARGUMENT;
        if (options_size >= layout_field_end)
        {
            const uint32_t raw_layout = options->shard_layout;
            if (raw_layout != LEO2_SHARD_LAYOUT_NATIVE_V1 &&
                raw_layout != LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
                return LEO2_INVALID_ARGUMENT;
            shard_layout = static_cast<leo2_shard_layout>(raw_layout);
        }
    }
    if (shard_layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1 &&
        field_value != LEO2_FIELD_GF16)
        return LEO2_INVALID_ARGUMENT;
    if (profile_value == LEO2_PROFILE_AUTO)
        profile_value = recovery_count <= original_count
            ? LEO2_PROFILE_LEGACY_HIGH_V1 : LEO2_PROFILE_LOW_V1;
    if (profile_value != LEO2_PROFILE_LEGACY_HIGH_V1 &&
        profile_value != LEO2_PROFILE_LOW_V1)
        return profile_value == LEO2_PROFILE_EXACT_EXPERIMENTAL_V1
            ? LEO2_UNSUPPORTED : LEO2_INVALID_ARGUMENT;
    const leo2_profile profile = static_cast<leo2_profile>(profile_value);

    const uint32_t padded = CeilPow2(
        profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? recovery_count : original_count);
    if (padded == 0)
        return LEO2_INVALID_COUNTS;
    const uint32_t parent = CeilPow2(static_cast<uint64_t>(padded) +
        (profile == LEO2_PROFILE_LEGACY_HIGH_V1 ? original_count : recovery_count));
    if (parent == 0)
        return LEO2_INVALID_COUNTS;

    if (field_value == LEO2_FIELD_AUTO)
        field_value = parent <= kGF8Order ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    if (field_value == LEO2_FIELD_GF8 && parent > kGF8Order)
        return LEO2_INVALID_COUNTS;
    if (field_value == LEO2_FIELD_GF16 && parent > kGF16Order)
        return LEO2_INVALID_COUNTS;
    if (field_value != LEO2_FIELD_GF8 && field_value != LEO2_FIELD_GF16)
        return LEO2_INVALID_ARGUMENT;
    const leo2_field field = static_cast<leo2_field>(field_value);
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
    codec->auto_avx512_encode_ops = NULL;
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
    codec->direct_generator_numerator_log = 0;
    codec->direct_generator_numerator_valid = false;
#ifdef LEO2_ENABLE_TEST_HOOKS
    codec->test_encode_mode = LEO2_TEST_ENCODE_AUTO;
    codec->test_decode_mode = LEO2_TEST_DECODE_AUTO;
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
        const bool prepare_translated_low =
            specialized && CanUseTranslatedLowDecode(codec);
        const bool prepare_native_locator =
            CodecMayUseNativeLocator(codec);
        const bool prepare_native_high =
            CodecMayCreateNativeHighPlan(codec);
        if (prepare_translated_low)
        {
            codec->translated_low_permanent_erased.resize(parent);
            for (uint32_t i = 0; i < parent; ++i)
            {
                codec->translated_low_permanent_erased[i] =
                    codec->permanent_erased[
                        TranslateHalfParentCoordinate(codec, i)];
            }
        }
        const uint32_t permanent_erasure_count =
            profile == LEO2_PROFILE_LEGACY_HIGH_V1
                ? padded - recovery_count
                : parent - padded - recovery_count;
        /*
            Every transform plan later retains exactly K of K+R transmitted
            coordinates, so its non-permanent missing plus virtual count is
            exactly recovery_count.  Cache the permanent coordinate-domain
            locator precisely when that dynamic contribution will use direct
            construction.  Dense Walsh construction from the explicit union
            avoids the extra N-entry modular-add pass required by either a
            coordinate- or transform-domain permanent cache.
        */

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
                    if (IsMeasuredEqualRoundedDirectRepairCodec(codec) &&
                        !weights.empty() &&
                        std::all_of(weights.begin(), weights.end(),
                            [&](DirectField8::Element value) {
                                return value == weights[0];
                            }))
                    {
                        DirectField8::Element numerator_log =
                            DirectField8::Log(weights[0]);
                        for (uint32_t coordinate = codec->padded_side;
                             coordinate < codec->parent_count;
                             ++coordinate)
                        {
                            numerator_log = DirectField8::AddLogs(
                                numerator_log, DirectField8::Log(
                                    static_cast<DirectField8::Element>(
                                        coordinate)));
                        }
                        codec->direct_generator_numerator_log = numerator_log;
                        codec->direct_generator_numerator_valid = true;
                    }
                    if (prepare_direct_encode)
                        PrepareDirectGeneratorLogs<DirectField8>(codec, weights,
                            codec->direct_generator_logs8);
                    if (prepare_direct_repair)
                    {
                        if (IsExpandedDirectRepairCodec(codec))
                        {
                            std::call_once(
                                codec->context->direct_repair_rows8_once,
                                [&]() {
                                    leo2_codec canonical = {};
                                    canonical.original_count = 65;
                                    canonical.recovery_count = 128;
                                    canonical.parent_count = 256;
                                    canonical.padded_side = 128;
                                    canonical.parent_dimension = 128;
                                    canonical.profile =
                                        LEO2_PROFILE_LEGACY_HIGH_V1;
                                    canonical.field = LEO2_FIELD_GF8;
                                    canonical.direct_generator_numerator_log =
                                        codec->direct_generator_numerator_log;
                                    canonical.direct_generator_numerator_valid =
                                        codec->
                                            direct_generator_numerator_valid;
                                    PrepareDirectRepairGeneratorRows<
                                        DirectField8>(&canonical, weights,
                                        codec->context->
                                            direct_repair_generator_rows8);
                                });
                        }
                        codec->direct_barycentric8.swap(weights);
                    }
                }
            }
            if (permanent_erasure_count != 0 &&
                leopard::ff8::IsDirectLocatorPreferred(parent, recovery_count))
            {
                if (prepare_native_locator)
                {
                    codec->permanent_locator8.resize(parent);
                    leopard::ff8::PrepareDecode(parent,
                        &codec->permanent_erased[0],
                        &codec->permanent_locator8[0]);
                }
                if (!codec->translated_low_permanent_erased.empty())
                {
                    codec->translated_low_permanent_locator8.resize(parent);
                    leopard::ff8::PrepareDecode(parent,
                        codec->translated_low_permanent_erased.data(),
                        codec->translated_low_permanent_locator8.data());
                }
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
                    if (prepare_native_high)
                    {
                        codec->fixed_factors8.resize(parent);
                        leopard::ff8::PrepareHighDecode(
                            parent, padded, &codec->fixed_factors8[0]);
                    }
                    if (prepare_translated_low)
                    {
                        codec->translated_low_factors8.resize(
                            parent / padded - 1);
                        leopard::ff8::PrepareLowDecode(parent, padded,
                            codec->translated_low_factors8.data());
                    }
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
                if (prepare_native_locator)
                {
                    codec->permanent_locator16.resize(parent);
                    leopard::ff16::PrepareDecode(parent,
                        &codec->permanent_erased[0],
                        &codec->permanent_locator16[0]);
                }
                if (!codec->translated_low_permanent_erased.empty())
                {
                    codec->translated_low_permanent_locator16.resize(parent);
                    leopard::ff16::PrepareDecode(parent,
                        codec->translated_low_permanent_erased.data(),
                        codec->translated_low_permanent_locator16.data());
                }
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
                    if (prepare_native_high)
                    {
                        codec->fixed_factors16.resize(parent);
                        leopard::ff16::PrepareHighDecode(
                            parent, padded, &codec->fixed_factors16[0]);
                    }
                    if (prepare_translated_low)
                    {
                        codec->translated_low_factors16.resize(
                            parent / padded - 1);
                        leopard::ff16::PrepareLowDecode(parent, padded,
                            codec->translated_low_factors16.data());
                    }
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
    if (CodecMayUseAutoAVX512Encode(codec))
    {
        // Qualification is setup-only, serialized, cached, and includes the
        // complete backend known-answer test.  AUTO safely retains its
        // baseline table when AVX-512 is absent or fails qualification.
        codec->auto_avx512_encode_ops =
            leopard::backend::GetQualifiedOps(LEO2_BACKEND_AVX512);
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
    int mode_value)
{
    if (!codec || (mode_value != LEO2_TEST_ENCODE_AUTO &&
                   mode_value != LEO2_TEST_ENCODE_FORCE_DIRECT &&
                   mode_value != LEO2_TEST_ENCODE_FORCE_TRANSFORM))
        return LEO2_INVALID_ARGUMENT;
    const leo2_test_encode_mode mode =
        static_cast<leo2_test_encode_mode>(mode_value);
    if (mode == LEO2_TEST_ENCODE_FORCE_DIRECT &&
        !HasDirectGeneratorRows(codec))
        return LEO2_UNSUPPORTED;
    codec->test_encode_mode = mode;
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_test_codec_set_decode_mode(
    leo2_codec* codec,
    int mode_value)
{
    if (!codec || (mode_value != LEO2_TEST_DECODE_AUTO &&
                   mode_value != LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW &&
                   mode_value != LEO2_TEST_DECODE_FORCE_NATIVE_HIGH))
        return LEO2_INVALID_ARGUMENT;
    const leo2_test_decode_mode mode =
        static_cast<leo2_test_decode_mode>(mode_value);
    if (mode != LEO2_TEST_DECODE_AUTO)
    {
        if (!CanUseTranslatedLowDecode(codec))
            return LEO2_UNSUPPORTED;
        if ((codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) != 0)
            return LEO2_INVALID_ARGUMENT;
    }
    codec->test_decode_mode = mode;
    return LEO2_SUCCESS;
}

LEO2_EXPORT int leo2_test_codec_translated_low_capable(
    const leo2_codec* codec)
{
    return CanUseTranslatedLowDecode(codec) ? 1 : 0;
}

LEO2_EXPORT int leo2_test_decode_plan_uses_translated_low(
    const leo2_decode_plan* plan)
{
    return PlanUsesTranslatedLowDecode(plan) ? 1 : 0;
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
    if (!IsRepresentableOutput(direct_out))
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

LEO2_EXPORT leo2_result leo2_test_codec_transform_encode_backend(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    uint32_t requested_recovery_count,
    uint32_t requested_recovery_prefix,
    leo2_backend* backend_out)
{
    if (!IsRepresentableOutput(backend_out))
        return LEO2_INVALID_ARGUMENT;
    *backend_out = LEO2_BACKEND_AUTO;
    if (!codec || shard_bytes == 0 ||
        shard_bytes > static_cast<uint64_t>(SIZE_MAX) ||
        requested_recovery_count > codec->recovery_count ||
        requested_recovery_prefix > codec->recovery_count ||
        requested_recovery_count > requested_recovery_prefix)
        return LEO2_INVALID_ARGUMENT;
    *backend_out = SelectTransformEncodeOps(codec,
        static_cast<size_t>(shard_bytes), requested_recovery_count,
        requested_recovery_prefix).kind;
    return LEO2_SUCCESS;
}
#endif

LEO2_EXPORT leo2_result leo2_codec_wire_shard_bytes(
    const leo2_codec* codec,
    uint64_t payload_bytes,
    uint64_t* wire_shard_bytes_out)
{
    if (!IsRepresentableOutput(wire_shard_bytes_out))
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
    AddressRange payload_range;
    AddressRange wire_range;
    if (!MakeRange(payload, payload_bytes, payload_range) ||
        !MakeRange(wire_shard, expected_wire_bytes, wire_range))
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
    AddressRange wire_range;
    AddressRange payload_range;
    if (!MakeRange(wire_shard, expected_wire_bytes, wire_range) ||
        !MakeRange(payload, payload_bytes, payload_range))
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
    if (!IsRepresentableOutput(scratch_bytes_out))
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

} // extern "C"

namespace leopard2_internal {

bool GetCodecDecodeMetadataInfo(
    const leo2_codec* codec,
    CodecDecodeMetadataInfo* info_out)
{
    if (!codec || !info_out)
        return false;
    CodecDecodeMetadataInfo info;
    info.permanent_erased_bytes = codec->permanent_erased.size();
    info.native_locator_bytes = codec->permanent_locator8.size() +
        codec->permanent_locator16.size() * sizeof(uint16_t);
    info.native_factor_bytes = codec->fixed_factors8.size() +
        codec->fixed_factors16.size() * sizeof(uint16_t);
    info.translated_permanent_erased_bytes =
        codec->translated_low_permanent_erased.size();
    info.translated_locator_bytes =
        codec->translated_low_permanent_locator8.size() +
        codec->translated_low_permanent_locator16.size() * sizeof(uint16_t);
    info.translated_factor_bytes = codec->translated_low_factors8.size() +
        codec->translated_low_factors16.size() * sizeof(uint16_t);
    *info_out = info;
    return true;
}

bool GetDecodePlanPrunedScheduleInfo(
    const leo2_decode_plan* plan,
    DecodePlanPrunedScheduleInfo* info_out)
{
    if (!plan || !info_out)
        return false;
    DecodePlanPrunedScheduleInfo info;
    info.low_input_plan_count = plan->low_pruned_input_blocks.size();
    info.low_output_plan_count =
        plan->low_pruned_output_plan.size == 0 ? 0 : 1;
    info.high_input_plan_count = plan->high_pruned_input_blocks.size();
    info.high_output_plan_count = 0;
    for (size_t i = 0; i < plan->high_pruned_output_plans.size(); ++i)
    {
        if (plan->high_pruned_output_plans[i].size != 0)
            ++info.high_output_plan_count;
    }
    *info_out = info;
    return true;
}

} // namespace leopard2_internal

static leo2_result EncodeInternal(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes,
    const void* protected_metadata,
    size_t protected_metadata_bytes,
    bool prevalidated)
{
    const bool single_side = UseSingleSideEncodeLayout(codec);
    if (single_side && IsLegacyHighK1R1Codec(codec))
    {
        size_t rounded_bytes = 0;
        if (!RoundShardBytes(shard_bytes, rounded_bytes))
            return LEO2_INVALID_ARGUMENT;
        if (codec->field == LEO2_FIELD_GF16 && (shard_bytes & 1u) != 0)
            return LEO2_UNSUPPORTED;
        if (!prevalidated)
        {
            AddressRange scratch_range;
            leo2_result result = CheckOptionalScratch(
                scratch, scratch_bytes, scratch_range);
            if (result != LEO2_SUCCESS)
                return result;
            result = ValidateK1R1EncodeBuffers(
                codec, shard_bytes, original, recovery, scratch_range,
                protected_metadata, protected_metadata_bytes);
            if (result != LEO2_SUCCESS)
                return result;
        }
        if (recovery[0])
            memcpy(recovery[0], original[0], static_cast<size_t>(shard_bytes));
        return LEO2_SUCCESS;
    }

    if (single_side)
    {
        ScratchLayout direct_layout;
        size_t rounded_bytes = 0;
        leo2_result result = DirectEncodeLayout(
            codec, shard_bytes, direct_layout, rounded_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        if (!prevalidated)
        {
            AddressRange scratch_range;
            result = CheckScratch(
                scratch, scratch_bytes, direct_layout, scratch_range);
            if (result != LEO2_SUCCESS)
                return result;
            result = ValidateEncodeBuffers(
                codec, shard_bytes, original, recovery, scratch,
                scratch_bytes, protected_metadata,
                protected_metadata_bytes);
            if (result != LEO2_SUCCESS)
                return result;
        }
        ExecuteSingleSideEncode(
            codec, static_cast<size_t>(shard_bytes), original, recovery);
        return LEO2_SUCCESS;
    }

    EncodeScratchGeometry geometry;
    leo2_result result = EncodeLayout(codec, shard_bytes, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    if (!prevalidated)
    {
        AddressRange scratch_range;
        result = CheckScratch(
            scratch, scratch_bytes, geometry.layout, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        result = ValidateEncodeBuffers(
            codec, shard_bytes, original, recovery, scratch, scratch_bytes,
            protected_metadata, protected_metadata_bytes);
        if (result != LEO2_SUCCESS)
            return result;
    }

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

    // Choose once for the complete call so aligned and ragged passes cannot
    // observe different arithmetic tables.  The choice is deterministic and
    // depends only on immutable codec state, output shape, and shard length.
    const leopard::backend::Ops& transform_ops = SelectTransformEncodeOps(
        codec, static_cast<size_t>(shard_bytes), requested_recovery_count,
        requested_recovery_prefix);

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
    leopard2_internal::SparseForwardPlanBatchView sparse_plans;
    if (!PrepareSparseEncodePlans(
            codec, geometry, recovery, requested_recovery_count,
            base, sparse_plans))
        return LEO2_INTERNAL_ERROR;

    if (geometry.aligned_prefix_bytes != 0)
    {
        if (geometry.execution_tile_count == 1)
        {
            if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
            {
                for (size_t i = 0; i < geometry.work_count; ++i)
                {
                    work[i] = i < requested_recovery_prefix && recovery[i]
                        ? recovery[i]
                        : work_storage + i * geometry.work_slot_bytes;
                }
            }
            else
            {
                for (size_t i = 0; i < geometry.work_count; ++i)
                    work[i] = work_storage + i * geometry.work_slot_bytes;
                for (uint32_t i = 0; i < codec->recovery_count; ++i)
                    parity[i] = recovery[i];
            }
            ExecuteTransformEncodePass(
                codec, transform_ops, geometry.aligned_prefix_bytes,
                geometry.aligned_prefix_bytes,
                requested_recovery_count, requested_recovery_prefix,
                original, parity, work, &sparse_plans);
        }
        else
        {
            const size_t total_tiles =
                geometry.aligned_prefix_bytes / kScratchAlignment;
            size_t first_tile = 0;
            for (size_t pass_index = 0;
                pass_index < geometry.execution_tile_count; ++pass_index)
            {
                const size_t passes_left =
                    geometry.execution_tile_count - pass_index;
                const size_t tiles_left = total_tiles - first_tile;
                const size_t pass_tiles =
                    tiles_left / passes_left +
                    (tiles_left % passes_left != 0);
                const size_t pass_bytes = pass_tiles * kScratchAlignment;
                LEO_DEBUG_ASSERT(pass_bytes <= geometry.work_slot_bytes);
                const size_t source_offset = first_tile * kScratchAlignment;
                PopulateEncodeInputs(
                    codec, original, source_offset, pass_bytes, NULL, pointers);

                if (codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1)
                {
                    for (size_t i = 0; i < geometry.work_count; ++i)
                    {
                        work[i] = i < requested_recovery_prefix && recovery[i]
                            ? static_cast<uint8_t*>(recovery[i]) + source_offset
                            : work_storage + i * geometry.work_slot_bytes;
                    }
                }
                else
                {
                    for (size_t i = 0; i < geometry.work_count; ++i)
                        work[i] = work_storage + i * geometry.work_slot_bytes;
                    for (uint32_t i = 0; i < codec->recovery_count; ++i)
                        parity[i] = recovery[i]
                            ? static_cast<uint8_t*>(recovery[i]) + source_offset
                            : NULL;
                }
                const void* const* const pass_original =
                    const_cast<const void* const*>(pointers);
                ExecuteTransformEncodePass(
                    codec, transform_ops, pass_bytes,
                    geometry.aligned_prefix_bytes,
                    requested_recovery_count, requested_recovery_prefix,
                    pass_original, parity, work, &sparse_plans);
                first_tile += pass_tiles;
            }
            LEO_DEBUG_ASSERT(first_tile == total_tiles);
        }
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
            codec, transform_ops, kScratchAlignment, kScratchAlignment,
            requested_recovery_count,
            requested_recovery_prefix,
            padded_original, parity, work, &sparse_plans);

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

extern "C" {

LEO2_EXPORT leo2_result leo2_encode(
    const leo2_codec* codec,
    uint64_t shard_bytes,
    const void* const* original,
    void* const* recovery,
    void* scratch,
    size_t scratch_bytes)
{
    return EncodeInternal(codec, shard_bytes, original, recovery,
        scratch, scratch_bytes, NULL, 0, false);
}

static leo2_result EncodeBatchCompatibilityInternal(
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
    size_t item_bytes = 0;
    AddressRange item_range = { 0, 0 };
    if (item_count != 0 &&
        (!CheckedMultiply(item_count, sizeof(*items), item_bytes) ||
         !MakeRange(items, static_cast<uint64_t>(item_bytes), item_range)))
        return LEO2_INVALID_ARGUMENT;
    /*
        A one-item batch has no cross-item aliasing to validate.  Route it
        through the ordinary one-shot validator while protecting the batch
        descriptor itself.  This preserves the batch metadata-overlap
        contract and avoids repeating the K/R addressability scans, overlap
        scans, and range sorts performed by the general batch preflight.
    */
    if (item_count == 1)
    {
        return EncodeInternal(codec, items[0].shard_bytes,
            items[0].original, items[0].recovery, items[0].scratch,
            items[0].scratch_bytes, items, item_bytes, false);
    }
    EncodeBatchTaskContext batch = {
        codec, items, item_bytes, item_range
    };
    leo2_result (*const preflight)(void*) =
        UseSingleSideEncodeLayout(codec) && IsLegacyHighK1R1Codec(codec)
            ? PreflightK1R1EncodeBatchTask
            : PreflightEncodeBatchTask;
    if (codec->context->pool)
        return codec->context->pool->Run(item_count,
            RunEncodeBatchItem, &batch, preflight);
    if (item_count != 0)
    {
        const leo2_result result = preflight(&batch);
        if (result != LEO2_SUCCESS)
            return result;
    }
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = RunEncodeBatchItem(&batch, i);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_encode_batch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count)
{
    return EncodeBatchCompatibilityInternal(codec, items, item_count);
}

LEO2_EXPORT leo2_result leo2_encode_batch_preflight_scratch_size(
    const leo2_codec* codec,
    size_t item_count,
    size_t* scratch_bytes_out)
{
    if (!IsRepresentableOutput(scratch_bytes_out))
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    return EncodeBatchPreflightScratchBytes(
        codec, item_count, *scratch_bytes_out);
}

LEO2_EXPORT leo2_result leo2_encode_batch_with_preflight_scratch(
    const leo2_codec* codec,
    const leo2_encode_batch_item* items,
    size_t item_count,
    void* preflight_scratch,
    size_t preflight_scratch_bytes)
{
    /* Preserve the established empty/single-item validation, metadata
       protection, and lazy-worker behavior without imposing a new workspace. */
    if (item_count < kScalableBatchPreflightMinItems)
        return EncodeBatchCompatibilityInternal(codec, items, item_count);
    if (item_count > 0xffffffffu || !items || !codec)
        return LEO2_INVALID_ARGUMENT;
    size_t item_bytes = 0;
    AddressRange item_range = { 0, 0 };
    if (!CheckedMultiply(item_count, sizeof(*items), item_bytes) ||
        !MakeRange(items, static_cast<uint64_t>(item_bytes), item_range))
        return LEO2_INVALID_ARGUMENT;

    ScalableEncodeBatchPreflightContext scalable = {
        { codec, items, item_bytes, item_range },
        preflight_scratch,
        preflight_scratch_bytes
    };
    if (codec->context->pool)
        return codec->context->pool->Run(item_count,
            RunScalableEncodeBatchItem, &scalable,
            PreflightEncodeBatchScalableTask);
    const leo2_result preflight =
        PreflightEncodeBatchScalableTask(&scalable);
    if (preflight != LEO2_SUCCESS)
        return preflight;
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result =
            RunScalableEncodeBatchItem(&scalable, i);
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
    AddressRange output_range;
    if (!MakeArrayRange(plan_out, 1, sizeof(*plan_out), output_range))
        return LEO2_INVALID_ARGUMENT;
    if (!codec)
    {
        *plan_out = NULL;
        return LEO2_INVALID_ARGUMENT;
    }

    AddressRange original_range = { 0, 0 };
    AddressRange recovery_range = { 0, 0 };
    if ((original_present &&
         !MakeArrayRange(original_present, codec->original_count,
             sizeof(*original_present), original_range)) ||
        (recovery_present &&
         !MakeArrayRange(recovery_present, codec->recovery_count,
             sizeof(*recovery_present), recovery_range)))
        return LEO2_INVALID_ARGUMENT;
    if ((original_present && RangesOverlap(output_range, original_range)) ||
        (recovery_present && RangesOverlap(output_range, recovery_range)))
        return LEO2_OVERLAP;
    *plan_out = NULL;
    if (!original_present || !recovery_present)
        return LEO2_INVALID_ARGUMENT;
    uint32_t present_count = 0;
    uint32_t missing_original_count = 0;
    uint32_t single_missing_original = std::numeric_limits<uint32_t>::max();
    for (uint32_t i = 0; i < codec->original_count; ++i)
    {
        if (original_present[i] > 1)
            return LEO2_INVALID_ARGUMENT;
        if (original_present[i])
            ++present_count;
        else
        {
            ++missing_original_count;
            single_missing_original = i;
        }
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
        plan->translated_low = false;
#ifdef LEO2_ENABLE_TEST_HOOKS
        if (codec->test_decode_mode != LEO2_TEST_DECODE_AUTO &&
            !CanUseTranslatedLowDecode(codec))
        {
            delete plan;
            return LEO2_INTERNAL_ERROR;
        }
#endif
        plan->generic_input_count = 0;
        plan->direct_copy_recovery = std::numeric_limits<uint32_t>::max();
        plan->direct_xor_original = std::numeric_limits<uint32_t>::max();
        plan->missing_original_count = missing_original_count;
        plan->no_op = missing_original_count == 0;
        plan->direct_xor = codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;
        plan->direct_copy = codec->profile == LEO2_PROFILE_LOW_V1 &&
            codec->padded_side == 1 && missing_original_count == 1;
        plan->direct_repair = false;

        if (plan->direct_xor &&
            (codec->original_count == 1 ||
             codec->original_count >= kCompactDirectXorMinOriginalCount))
        {
            if (single_missing_original >= codec->original_count)
            {
                delete plan;
                return LEO2_INTERNAL_ERROR;
            }
            plan->direct_xor_original = single_missing_original;
            *plan_out = plan;
            return LEO2_SUCCESS;
        }

        plan->original_present.assign(
            original_present, original_present + codec->original_count);
        plan->recovery_present.assign(
            recovery_present, recovery_present + codec->recovery_count);
        plan->missing_originals.reserve(missing_original_count);
        for (uint32_t i = 0; i < codec->original_count; ++i)
            if (!original_present[i])
                plan->missing_originals.push_back(i);

        bool force_test_transform = false;
#ifdef LEO2_ENABLE_TEST_HOOKS
        force_test_transform =
            codec->test_decode_mode != LEO2_TEST_DECODE_AUTO;
#endif
        if (!plan->no_op && !plan->direct_xor && !plan->direct_copy &&
            !force_test_transform &&
            missing_original_count <= DirectRepairLossLimit(codec))
        {
#ifdef LEO_HAS_FF8
            if (codec->field == LEO2_FIELD_GF8 &&
                !codec->direct_barycentric8.empty())
            {
                plan->direct_repair = PrepareDirectRepairTerms<DirectField8>(
                    plan, codec->direct_barycentric8,
                    IsExpandedDirectRepairCodec(codec)
                        ? &codec->context->direct_repair_generator_rows8
                        : NULL);
            }
#endif
#ifdef LEO_HAS_FF16
            if (codec->field == LEO2_FIELD_GF16 &&
                !codec->direct_barycentric16.empty())
            {
                plan->direct_repair = PrepareDirectRepairTerms<DirectField16>(
                    plan, codec->direct_barycentric16, NULL);
            }
#endif
        }

        if (!plan->direct_repair)
        {
            plan->coordinate_erased = codec->permanent_erased;
            plan->requested_coordinates.reserve(missing_original_count);
            for (uint32_t i = 0; i < missing_original_count; ++i)
            {
                const uint32_t original = plan->missing_originals[i];
                const uint32_t coordinate =
                    CoordinateForOriginal(codec, original);
                plan->coordinate_erased[coordinate] = 1;
                plan->requested_coordinates.push_back(coordinate);
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

#ifdef LEO_DEBUG
        /*
            The transmitted code has K + R coordinates and the deterministic
            selection above retains exactly K of them.  Consequently the
            pattern-specific missing plus virtual coordinates outside the
            codec's permanent shortening/puncturing mask must always total R.
            Permanent-locator cache policy relies on this invariant when it
            predicts the direct/Walsh crossover from recovery_count at codec
            setup rather than rescanning each future plan.
        */
        uint32_t dynamic_erasure_count = 0;
        for (uint32_t i = 0; i < codec->parent_count; ++i)
        {
            if (plan->coordinate_erased[i] && !codec->permanent_erased[i])
                ++dynamic_erasure_count;
        }
        LEO_DEBUG_ASSERT(dynamic_erasure_count == codec->recovery_count);
#endif

            if (plan->direct_copy)
            {
                for (uint32_t i = 0; i < codec->recovery_count; ++i)
                {
                    const uint32_t coordinate = CoordinateForRecovery(codec, i);
                    if (recovery_present[i] &&
                        !plan->coordinate_erased[coordinate])
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
        }

        if (!plan->no_op && !plan->direct_xor && !plan->direct_copy &&
            !plan->direct_repair)
        {
            plan->translated_low =
                (codec->flags & LEO2_CODEC_FORCE_GENERIC_DECODE) == 0 &&
                CanUseTranslatedLowDecode(codec);
#ifdef LEO2_ENABLE_TEST_HOOKS
            if (codec->test_decode_mode == LEO2_TEST_DECODE_FORCE_NATIVE_HIGH)
                plan->translated_low = false;
            else if (codec->test_decode_mode ==
                     LEO2_TEST_DECODE_FORCE_TRANSLATED_LOW)
                plan->translated_low = true;
#endif
            if (!PreparePlanExecutionMetadata(plan))
            {
                delete plan;
                return LEO2_INTERNAL_ERROR;
            }
            // From this point onward coordinate_erased follows the execution
            // view.  A translated plan has no generic fallback metadata, so
            // constructing its locator directly in low-profile order avoids
            // retaining and touching a second N-entry locator vector.
            TranslatePlanErasureMaskForExecution(plan);
#ifdef LEO_HAS_FF8
            if (codec->field == LEO2_FIELD_GF8)
            {
                plan->locator8.resize(codec->parent_count);
                const bool translated = PlanUsesTranslatedLowDecode(plan);
                const std::vector<uint8_t>& permanent_locator = translated
                    ? codec->translated_low_permanent_locator8
                    : codec->permanent_locator8;
                const std::vector<uint8_t>& permanent_erased = translated
                    ? codec->translated_low_permanent_erased
                    : codec->permanent_erased;
                if (permanent_locator.empty())
                    leopard::ff8::PrepareDecode(codec->parent_count,
                        &plan->coordinate_erased[0], &plan->locator8[0]);
                else
                    leopard::ff8::PrepareDecodeWithPermanent(
                        codec->parent_count, &plan->coordinate_erased[0],
                        permanent_erased.data(), permanent_locator.data(),
                        &plan->locator8[0]);
            }
#endif
#ifdef LEO_HAS_FF16
            if (codec->field == LEO2_FIELD_GF16)
            {
                plan->locator16.resize(codec->parent_count);
                const bool translated = PlanUsesTranslatedLowDecode(plan);
                const std::vector<uint16_t>& permanent_locator = translated
                    ? codec->translated_low_permanent_locator16
                    : codec->permanent_locator16;
                const std::vector<uint8_t>& permanent_erased = translated
                    ? codec->translated_low_permanent_erased
                    : codec->permanent_erased;
                if (permanent_locator.empty())
                    leopard::ff16::PrepareDecode(codec->parent_count,
                        &plan->coordinate_erased[0], &plan->locator16[0]);
                else
                    leopard::ff16::PrepareDecodeWithPermanent(
                        codec->parent_count, &plan->coordinate_erased[0],
                        permanent_erased.data(), permanent_locator.data(),
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
    if (!IsRepresentableOutput(scratch_bytes_out))
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
        : DecodeLayout(plan->codec, plan, shard_bytes, false, geometry);
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
    bool multi_item_batch,
    const ProtectedMetadataSpan* protected_metadata,
    size_t protected_metadata_count,
    bool prevalidated)
{
    if (!plan)
        return LEO2_INVALID_ARGUMENT;
    if (plan->no_op)
        return LEO2_SUCCESS;
    if (shard_bytes == 0)
        return LEO2_INVALID_ARGUMENT;
    if (PlanHasCompactDirectXor(plan))
    {
        if (!prevalidated)
        {
            if (protected_metadata_count > 2 ||
                (protected_metadata_count != 0 && !protected_metadata))
                return LEO2_INVALID_ARGUMENT;
            AddressRange protected_ranges[2];
            for (size_t i = 0; i < protected_metadata_count; ++i)
            {
                if (!protected_metadata[i].data ||
                    protected_metadata[i].bytes == 0 ||
                    !MakeRange(protected_metadata[i].data,
                        static_cast<uint64_t>(protected_metadata[i].bytes),
                        protected_ranges[i]))
                    return LEO2_INVALID_ARGUMENT;
            }
            return DecodeDirectXorValidated(plan->codec,
                plan->direct_xor_original, shard_bytes, protected_ranges,
                protected_metadata_count, original, recovery,
                restored_original, scratch, scratch_bytes);
        }
        if (!recovery || !recovery[0])
            return LEO2_INTERNAL_ERROR;
        if (IsLegacyHighK1R1Codec(plan->codec) &&
            plan->direct_xor_original == 0)
        {
            memcpy(restored_original[0], recovery[0],
                static_cast<size_t>(shard_bytes));
            return LEO2_SUCCESS;
        }
        ExecuteDirectXorDecode(plan->codec, plan->direct_xor_original,
            static_cast<size_t>(shard_bytes), original, recovery,
            restored_original);
        return LEO2_SUCCESS;
    }

    ScratchLayout direct_layout;
    size_t direct_rounded = 0;
    DecodeScratchGeometry geometry;
    const bool direct_execution =
        plan->direct_repair || plan->direct_xor || plan->direct_copy;
    leo2_result result = direct_execution
        ? DirectDecodeLayout(
            plan->codec, shard_bytes, direct_layout, direct_rounded)
        : DecodeLayout(plan->codec, plan, shard_bytes,
            multi_item_batch, geometry);
    if (result != LEO2_SUCCESS)
        return result;
    const ScratchLayout& layout = direct_execution
        ? direct_layout
        : geometry.layout;
    bool coordinates_prepopulated = false;
    if (!prevalidated)
    {
        AddressRange scratch_range;
        result = CheckScratch(scratch, scratch_bytes, layout, scratch_range);
        if (result != LEO2_SUCCESS)
            return result;
        void** validation_coordinate_data = NULL;
        if (!direct_execution && geometry.aligned_prefix_bytes != 0 &&
            plan->missing_original_count > 1)
        {
            validation_coordinate_data = reinterpret_cast<void**>(
                static_cast<uint8_t*>(scratch) + layout.pointer_offset);
        }
        result = ValidateDecodeBuffers(
            plan, shard_bytes, original, recovery, restored_original,
            scratch, scratch_bytes, protected_metadata,
            protected_metadata_count, validation_coordinate_data);
        if (result != LEO2_SUCCESS)
            return result;
        coordinates_prepopulated = validation_coordinate_data != NULL;
    }

    const leo2_codec* codec = plan->codec;
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
        ExecuteDirectXorDecode(codec, missing,
            static_cast<size_t>(shard_bytes), original, recovery,
            restored_original);
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
        through 1 MiB on scalar, SSSE3, and AVX2.  The explicit AVX-512VL
        table recompiles the AVX2 traversal with the expanded register file
        and follows the same rule.  Keep dispatch strictly inside that region;
        neighboring counts, fields, backends, and sizes retain the
        profile-specific decoder.  A ragged final tile uses the same path as
        its aligned prefix.
    */
    const bool use_generic = geometry.selection.path ==
        leopard2_internal::kDecodePathGeneric;
    const bool fuse_generic_reveal_scatter = use_generic &&
        UseFusedGenericRevealScatter(codec, geometry.aligned_prefix_bytes);
    const bool fuse_low_reveal_scatter = !use_generic &&
        UseFusedLowRevealScatter(plan, geometry.aligned_prefix_bytes);
    const bool reveal_aligned_outputs_in_place =
        !(fuse_generic_reveal_scatter || fuse_low_reveal_scatter);
    const bool use_tiled = geometry.selection.path ==
        leopard2_internal::kDecodePathTiled;
    const bool fuse_high_reveal_scatter = !use_generic &&
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        !PlanUsesTranslatedLowDecode(plan) &&
        geometry.aligned_prefix_bytes != 0;
    void** const high_requested_output =
        work + static_cast<size_t>(codec->padded_side) * 2;

    if (geometry.aligned_prefix_bytes != 0)
    {
        if (fuse_high_reveal_scatter)
        {
            /*
                Complete kernel tiles have the public byte layout, and
                validation has proved that every restored output is disjoint
                from all received shards, scratch, metadata, and other
                outputs.  Tiled Algorithm 5 accepts compact requested-output
                pointers.  A single materialized Algorithm 5 pass indexes the
                public systematic array by coordinate-T.  In both cases the
                fixed reveal multiply writes the final destination instead of
                a scratch shard that GatherTransformDecodePass would copy.
            */
            if (use_tiled && geometry.execution_tile_count == 1)
            {
                for (size_t i = 0; i < plan->missing_originals.size(); ++i)
                {
                    const uint32_t original_index =
                        plan->missing_originals[i];
                    high_requested_output[i] = static_cast<uint8_t*>(
                        restored_original[original_index]);
                }
            }
#ifdef LEO2_ENABLE_TEST_HOOKS
            g_test_high_direct_reveal_shards.fetch_add(
                plan->missing_originals.size(), std::memory_order_relaxed);
            if (!use_tiled)
                g_test_high_materialized_direct_reveal_shards.fetch_add(
                    plan->missing_originals.size(),
                    std::memory_order_relaxed);
#endif
        }
        if (geometry.execution_tile_count == 1)
        {
            if (!coordinates_prepopulated)
            {
                PopulateDecodeCoordinates(
                    plan, original, recovery, 0,
                    geometry.aligned_prefix_bytes, NULL, coordinate_data);
            }
            const void* const* const coordinate_input =
                const_cast<const void* const*>(coordinate_data);
            ExecuteTransformDecodePass(
                plan, geometry.aligned_prefix_bytes, coordinate_input, work,
                use_generic, use_tiled, reveal_aligned_outputs_in_place,
                fuse_high_reveal_scatter && !use_tiled
                    ? restored_original : NULL);
            if (!fuse_high_reveal_scatter)
                GatherTransformDecodePass(
                    plan, restored_original, 0,
                    geometry.aligned_prefix_bytes, work, use_generic,
                    use_tiled, reveal_aligned_outputs_in_place);
        }
        else
        {
            /*
                GF8AVX2DecodeExecutionTileBytes admits a multi-pass native
                high decoder only when dispatch selected Algorithm 5's tiled
                traversal.  The only admitted non-tiled rule is TranslatedLow,
                for which fuse_high_reveal_scatter is false.  Keep this
                invariant explicit so payload tiling never needs a second K
                pointer array for a hypothetical materialized-high path.
            */
            LEO_DEBUG_ASSERT(!fuse_high_reveal_scatter || use_tiled);
            const size_t total_tiles =
                geometry.aligned_prefix_bytes / kScratchAlignment;
            size_t first_tile = 0;
            for (size_t pass_index = 0;
                pass_index < geometry.execution_tile_count; ++pass_index)
            {
                const size_t passes_left =
                    geometry.execution_tile_count - pass_index;
                const size_t tiles_left = total_tiles - first_tile;
                const size_t pass_tiles =
                    tiles_left / passes_left +
                    (tiles_left % passes_left != 0);
                const size_t pass_bytes =
                    pass_tiles * kScratchAlignment;
                LEO_DEBUG_ASSERT(pass_bytes <= geometry.work_slot_bytes);
                const size_t source_offset =
                    first_tile * kScratchAlignment;

                PopulateDecodeCoordinates(
                    plan, original, recovery, source_offset, pass_bytes,
                    NULL, coordinate_data);
                if (fuse_high_reveal_scatter)
                {
                    for (size_t i = 0;
                        i < plan->missing_originals.size(); ++i)
                    {
                        const uint32_t original_index =
                            plan->missing_originals[i];
                        high_requested_output[i] =
                            static_cast<uint8_t*>(
                                restored_original[original_index]) +
                            source_offset;
                    }
                }
                const void* const* const coordinate_input =
                    const_cast<const void* const*>(coordinate_data);
                ExecuteTransformDecodePass(
                    plan, pass_bytes, coordinate_input, work,
                    use_generic, use_tiled, reveal_aligned_outputs_in_place,
                    NULL);
                if (!fuse_high_reveal_scatter)
                    GatherTransformDecodePass(
                        plan, restored_original, source_offset, pass_bytes,
                        work, use_generic, use_tiled,
                        reveal_aligned_outputs_in_place);
                first_tile += pass_tiles;
            }
            LEO_DEBUG_ASSERT(first_tile == total_tiles);
        }
    }

    if (geometry.tail_bytes != 0)
    {
        if (fuse_high_reveal_scatter && use_tiled)
        {
            /* The compact tail still needs its retained kernel-layout shards. */
            const size_t high_output_begin =
                static_cast<size_t>(codec->padded_side) * 2;
            for (size_t i = 0; i < plan->missing_originals.size(); ++i)
            {
                const size_t slot = high_output_begin + i;
                high_requested_output[i] =
                    work_storage + slot * geometry.work_slot_bytes;
            }
        }
#ifdef LEO2_ENABLE_TEST_HOOKS
        if (!use_generic &&
            codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
            !PlanUsesTranslatedLowDecode(plan))
        {
            g_test_high_scratch_reveal_shards.fetch_add(
                plan->missing_originals.size(), std::memory_order_relaxed);
        }
#endif
        PopulateDecodeCoordinates(
            plan, original, recovery, geometry.aligned_prefix_bytes,
            geometry.tail_bytes, staging_slots, coordinate_data);
        const void* const* const coordinate_input =
            const_cast<const void* const*>(coordinate_data);
        ExecuteTransformDecodePass(
            plan, kScratchAlignment, coordinate_input, work,
            use_generic, use_tiled, true, NULL);
        GatherTransformDecodePass(
            plan, restored_original, geometry.aligned_prefix_bytes,
            geometry.tail_bytes, work, use_generic, use_tiled, true);
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
        scratch, scratch_bytes, false, NULL, 0, false);
}

static leo2_result DecodeBatchCompatibilityInternal(
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
    size_t item_bytes = 0;
    AddressRange item_range = { 0, 0 };
    if (item_count != 0 &&
        (!CheckedMultiply(item_count, sizeof(*items), item_bytes) ||
         !MakeRange(items, static_cast<uint64_t>(item_bytes), item_range)))
        return LEO2_INVALID_ARGUMENT;
    /* A no-loss plan is a true no-op: after validating only the top-level
       plan/item-array address arithmetic, do not inspect per-item pointers,
       scratch, sizes, or aliases and do not start the worker pool. */
    if (item_count == 0 || plan->no_op)
        return LEO2_SUCCESS;
    /* A one-item batch has no cross-item aliases.  Match the encoder's
       one-item fast path: run the ordinary validator once while protecting
       the batch descriptor, rather than repeating addressability scans and
       constructing/sorting K+R range tables in the general batch preflight. */
    if (item_count == 1)
    {
        const ProtectedMetadataSpan item_metadata = { items, item_bytes };
        return DecodePlanExecuteInternal(plan, items[0].shard_bytes,
            items[0].original, items[0].recovery,
            items[0].restored_original, items[0].scratch,
            items[0].scratch_bytes, false, &item_metadata, 1, false);
    }
    DecodeBatchTaskContext batch = {
        plan, items, item_bytes, item_range, item_count > 1
    };
    leo2_result (*const preflight)(void*) =
        PlanHasCompactDirectXor(plan) &&
                IsLegacyHighK1R1Codec(plan->codec) &&
                plan->direct_xor_original == 0
            ? PreflightK1R1DecodeBatchTask
            : PreflightDecodeBatchTask;
    if (plan->codec->context->pool)
        return plan->codec->context->pool->Run(item_count,
            RunDecodeBatchItem, &batch, preflight);
    const leo2_result preflight_result = preflight(&batch);
    if (preflight_result != LEO2_SUCCESS)
        return preflight_result;
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result = RunDecodeBatchItem(&batch, i);
        if (result != LEO2_SUCCESS)
            return result;
    }
    return LEO2_SUCCESS;
}

LEO2_EXPORT leo2_result leo2_decode_plan_execute_batch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count)
{
    return DecodeBatchCompatibilityInternal(plan, items, item_count);
}

LEO2_EXPORT leo2_result leo2_decode_plan_batch_preflight_scratch_size(
    const leo2_decode_plan* plan,
    size_t item_count,
    size_t* scratch_bytes_out)
{
    if (!IsRepresentableOutput(scratch_bytes_out))
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    return DecodeBatchPreflightScratchBytes(
        plan, item_count, *scratch_bytes_out);
}

LEO2_EXPORT leo2_result
leo2_decode_plan_execute_batch_with_preflight_scratch(
    const leo2_decode_plan* plan,
    const leo2_decode_batch_item* items,
    size_t item_count,
    void* preflight_scratch,
    size_t preflight_scratch_bytes)
{
    /* The compatibility entry point owns the exact true-no-op and one-item
       semantics, including validation of only top-level item-array arithmetic
       for a no-loss plan. */
    if (item_count < kScalableBatchPreflightMinItems ||
        (plan && plan->no_op))
        return DecodeBatchCompatibilityInternal(plan, items, item_count);
    if (item_count > 0xffffffffu || !items || !plan)
        return LEO2_INVALID_ARGUMENT;
    size_t item_bytes = 0;
    AddressRange item_range = { 0, 0 };
    if (!CheckedMultiply(item_count, sizeof(*items), item_bytes) ||
        !MakeRange(items, static_cast<uint64_t>(item_bytes), item_range))
        return LEO2_INVALID_ARGUMENT;

    ScalableDecodeBatchPreflightContext scalable = {
        { plan, items, item_bytes, item_range, true },
        preflight_scratch,
        preflight_scratch_bytes
    };
    if (plan->codec->context->pool)
        return plan->codec->context->pool->Run(item_count,
            RunScalableDecodeBatchItem, &scalable,
            PreflightDecodeBatchScalableTask);
    const leo2_result preflight =
        PreflightDecodeBatchScalableTask(&scalable);
    if (preflight != LEO2_SUCCESS)
        return preflight;
    for (size_t i = 0; i < item_count; ++i)
    {
        const leo2_result result =
            RunScalableDecodeBatchItem(&scalable, i);
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
    if (!IsRepresentableOutput(scratch_bytes_out))
        return LEO2_INVALID_ARGUMENT;
    *scratch_bytes_out = 0;
    if (IsLegacyHighK1R1Codec(codec))
    {
        ScratchLayout layout;
        size_t rounded_bytes = 0;
        const leo2_result result = DirectDecodeLayout(
            codec, shard_bytes, layout, rounded_bytes);
        if (result != LEO2_SUCCESS)
            return result;
        *scratch_bytes_out = layout.total_bytes;
        return LEO2_SUCCESS;
    }
    DecodeScratchGeometry geometry;
    const leo2_result result = DecodeLayout(
        codec, NULL, shard_bytes, false, geometry);
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
    /*
        In a legacy-high R=1 codeword, exactly one missing original implies
        that the parity and every other original are present.  Validate and
        execute that terminal operation directly in the one-shot API instead
        of allocating and immediately destroying a plan.  Unsafe or
        non-terminal presence metadata falls through to the canonical
        constructor so established error and no-op semantics remain there.
    */
    AddressRange original_presence_range = { 0, 0 };
    AddressRange recovery_presence_range = { 0, 0 };
    bool direct_xor_one_shot = codec &&
        codec->profile == LEO2_PROFILE_LEGACY_HIGH_V1 &&
        codec->padded_side == 1 && original_present && recovery_present &&
        MakeArrayRange(original_present, codec->original_count,
            sizeof(*original_present), original_presence_range) &&
        MakeArrayRange(recovery_present, codec->recovery_count,
            sizeof(*recovery_present), recovery_presence_range) &&
        recovery_present[0] == 1;
#ifdef LEO2_ENABLE_TEST_HOOKS
    direct_xor_one_shot = direct_xor_one_shot &&
        codec->test_decode_mode == LEO2_TEST_DECODE_AUTO;
#endif
    if (direct_xor_one_shot)
    {
        uint32_t missing = std::numeric_limits<uint32_t>::max();
        uint32_t missing_count = 0;
        bool valid = true;
        for (uint32_t i = 0; i < codec->original_count; ++i)
        {
            if (original_present[i] > 1)
            {
                valid = false;
                break;
            }
            if (!original_present[i])
            {
                missing = i;
                ++missing_count;
            }
        }
        if (valid && missing_count == 1)
        {
            const AddressRange protected_ranges[2] = {
                original_presence_range, recovery_presence_range
            };
            return DecodeDirectXorValidated(codec, missing, shard_bytes,
                protected_ranges, 2, original, recovery, restored_original,
                scratch, scratch_bytes);
        }
    }

    leo2_decode_plan* plan = NULL;
    leo2_result result = leo2_decode_plan_create(
        codec, original_present, recovery_present, &plan);
    if (result != LEO2_SUCCESS)
        return result;
    const ProtectedMetadataSpan protected_metadata[2] = {
        { original_present,
          codec->original_count * sizeof(*original_present) },
        { recovery_present,
          codec->recovery_count * sizeof(*recovery_present) }
    };
    result = DecodePlanExecuteInternal(
        plan, shard_bytes, original, recovery, restored_original, scratch,
        scratch_bytes, false, protected_metadata, 2, false);
    leo2_decode_plan_destroy(plan);
    return result;
}

} // extern "C"

namespace {

#if defined(_MSC_VER)
#define LEO2_R1_PAIR_POLICY_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__ELF__)
#define LEO2_R1_PAIR_POLICY_NOINLINE \
    __attribute__((noinline, section(".text.leo2_r1_pair_policy")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_R1_PAIR_POLICY_NOINLINE __attribute__((noinline))
#else
#define LEO2_R1_PAIR_POLICY_NOINLINE
#endif

static LEO2_R1_PAIR_POLICY_NOINLINE bool
UseAVX2GF8PairwiseLargeR1Xor(
    uint32_t original_count,
    size_t shard_bytes)
{
    return (original_count >= 51 && shard_bytes >= 1024U * 1024U) ||
        (original_count >= 31 && shard_bytes >= 2U * 1024U * 1024U) ||
        (original_count >= 10 && shard_bytes >= 4U * 1024U * 1024U) ||
        (original_count >= 8 && shard_bytes >= 8U * 1024U * 1024U);
}

#undef LEO2_R1_PAIR_POLICY_NOINLINE

#if defined(_MSC_VER)
#define LEO2_PLAN_NOINLINE __declspec(noinline)
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__ELF__)
#define LEO2_PLAN_NOINLINE \
    __attribute__((noinline, section(".text.leo2_plan_policy")))
#elif defined(__GNUC__) || defined(__clang__)
#define LEO2_PLAN_NOINLINE __attribute__((noinline))
#else
#define LEO2_PLAN_NOINLINE
#endif

static LEO2_PLAN_NOINLINE bool
SkipMeasuredBoundaryGF8AVX2HighPrunedSchedules(
    const leo2_decode_plan* plan)
{
    if (!plan || !plan->codec)
        return false;
    const leo2_codec* codec = plan->codec;
    if (codec->profile != LEO2_PROFILE_LEGACY_HIGH_V1 ||
        codec->field != LEO2_FIELD_GF8 || !codec->context ||
        !codec->context->ops ||
        codec->context->ops->kind != LEO2_BACKEND_AVX2 ||
        plan->missing_original_count != codec->recovery_count)
        return false;

    /*
        This measured extension covers the one-coordinate K boundary and one
        punctured parity at T=8, T=16, and T=64.  Omitting both schedules cut
        setup by 3.0--7.3x; the 48-cell pattern screen's worst execution change
        was a 0.7% loss.  Keeping the extension out of line and, on ELF,
        isolated in its own section preserves the existing hot-function layout.
    */
    const uint64_t side = codec->padded_side;
    const uint64_t recovery = codec->recovery_count;
    const uint64_t originals = codec->original_count;
    return (side == 8 || side == 16 || side == 64) &&
        recovery + 1U >= side && originals + 1U >= recovery * 2U;
}

#undef LEO2_PLAN_NOINLINE

} // namespace
