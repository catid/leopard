// leopard-79h.38.5.4.19.1.1: paired timer adapter; derived from 8f854b1.
// Qualification uses synthetic/abort clocks or clock-free exercise ONLY.
#ifdef LEO_PAIRED_NATIVE
#include "leopard.h"
#include "LeopardCommon.h"
#else
#include "leopard2.h"
#include "Leopard2Direct.h"
#endif
#include <chrono>
#include "PairedGroupTiming.h"
extern "C" const char* LeoPairedClockKind();
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>
#if defined(__SANITIZE_ADDRESS__)
#include <sanitizer/asan_interface.h>
#endif
#ifndef LEO_PAIRED_CODEC
#error Pin the unchanged codec archive identity.
#endif

namespace {
void Require(bool ok, const char* what) { if (!ok) throw std::runtime_error(what); }
struct Buffer
{
    uint8_t* raw;
    uint8_t* data;
    const size_t bytes;
    explicit Buffer(size_t size) : raw(NULL), data(NULL), bytes(size)
    {
        void* ptr = NULL;
        Require(posix_memalign(&ptr, 64, bytes + 128) == 0, "allocation");
        raw = static_cast<uint8_t*>(ptr); data = raw + 64;
        std::memset(raw, 0xd3, bytes + 128);
        std::memset(data, 0, bytes);
        Poison(true);
    }
    void Poison(bool enabled)
    {
#if defined(__SANITIZE_ADDRESS__)
        if (enabled) {
            __asan_poison_memory_region(raw, 64);
            __asan_poison_memory_region(data + bytes, 64);
        } else __asan_unpoison_memory_region(raw, bytes + 128);
#else
        (void)enabled;
#endif
    }
    void Check()
    {
        Poison(false);
        bool ok = true;
        for (unsigned i = 0; i < 64; ++i)
            ok &= raw[i] == 0xd3 && data[bytes + i] == 0xd3;
        Poison(true);
        Require(ok, "allocation guard changed");
    }
    ~Buffer() { Poison(false); std::free(raw); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
};
uint64_t Hash(const uint8_t* data, size_t bytes)
{
    uint64_t h = UINT64_C(14695981039346656037);
    for (size_t i = 0; i < bytes; ++i) h = (h ^ data[i]) * UINT64_C(1099511628211);
    return h;
}
struct Cell { unsigned k, r; size_t bytes; };
const Cell kCells[] = {
    {1000,199,32768}, {1000,199,32768}, {1000,200,32768},
    {1000,199,65536}, {1000,200,65536}, {1000,198,32768},
    {1000,199,32768}, {4096,512,4096}, {17,7,64}
};
}

int main(int argc, char** argv)
{
    try {
        Require((argc == 5 || argc == 6) &&
            (!std::strcmp(argv[1], "--check") || !std::strcmp(argv[1], "--exercise") ||
             !std::strcmp(argv[1], "--clock-guard") || !std::strcmp(argv[1], "--clock-exercise") ||
             !std::strcmp(argv[1], "--measure")),
            "usage: --check|--exercise|--clock-exercise|--clock-guard|--measure cell schedule 1|256 [new_parity_file]");
        Require(std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '8', "cell");
        const unsigned index = static_cast<unsigned>(argv[2][0] - '0');
        const Cell cell = kCells[index];
        const char* const schedule = argv[3];
#ifdef LEO_PAIRED_NATIVE
        Require(!std::strcmp(schedule, "NNNN"), "native schedule must be NNNN");
#else
        Require(!std::strcmp(schedule, "0110") || !std::strcmp(schedule, "1001") ||
                !std::strcmp(schedule, "0000") || !std::strcmp(schedule, "1111"), "schedule");
#endif
        Require(!std::strcmp(argv[4], "1") || (index == 8 && !std::strcmp(argv[4], "256")), "group");
        const unsigned group = !std::strcmp(argv[4], "1") ? 1 : 256;
        const bool exercise = std::strcmp(argv[1], "--check") != 0;
        const bool clock_guard = !std::strcmp(argv[1], "--clock-guard");
        const bool measured = !std::strcmp(argv[1], "--measure");
        const bool synthetic = !std::strcmp(argv[1], "--clock-exercise");
        const bool clocks_enabled = measured || synthetic || clock_guard;
        const char* const clock_kind = LeoPairedClockKind();
        Require(!synthetic || !std::strcmp(clock_kind, "synthetic"), "synthetic clock required");
        Require(!measured || !std::strcmp(clock_kind, "steady"), "steady clock required");
        Require(!clock_guard || !std::strcmp(clock_kind, "abort"), "abort clock required");
        Require(!(clock_guard || measured) || argc == 5, "no measured or clock-guard parity dump");
        const size_t input_bytes = static_cast<size_t>(cell.k) * cell.bytes;
        const size_t output_bytes = static_cast<size_t>(cell.r) * cell.bytes;
        Buffer source(input_bytes), reference(output_bytes);
        uint32_t random = 20260906;
        for (size_t i = 0; i < input_bytes; ++i) {
            random ^= random << 13; random ^= random >> 17; random ^= random << 5;
            source.data[i] = static_cast<uint8_t>(random);
        }
        std::vector<const void*> inputs(cell.k);
        for (unsigned i = 0; i < cell.k; ++i) inputs[i] = source.data + i * cell.bytes;
        const uint64_t input_hash = Hash(source.data, input_bytes);
        unsigned calls = 0, selections = 0, probes[4] = {};
        unsigned per_slot[4] = {};
        std::vector<leopard_paired::Sample> samples;
        samples.reserve(84);
        const auto now = []() -> int64_t {
            static_assert(std::ratio_equal<std::chrono::steady_clock::period, std::nano>::value,
                "This pinned Linux frontend requires nanosecond steady-clock ticks");
            return std::chrono::steady_clock::now().time_since_epoch().count();
        };
#ifdef LEO_PAIRED_NATIVE
#ifndef LEO_TRY_AVX2
#error Original native comparator must include AVX2.
#endif
        Require(leo_init() == 0 && leopard::CpuHasAVX2, "native initialization");
        const char* const api = "leo_encode";
        const unsigned count = leo_encode_work_count(cell.k, cell.r);
        Require(count >= cell.r, "native scratch count");
        const size_t scratch_bytes = static_cast<size_t>(count) * cell.bytes;
        Buffer scratch(scratch_bytes);
        uint8_t* const parity = scratch.data;
        std::vector<void*> outputs(count);
        for (unsigned i = 0; i < count; ++i) outputs[i] = scratch.data + i * cell.bytes;
        const auto select = [&](unsigned) { ++selections; };
        const auto encode = [&]() {
            Require(leo_encode(cell.bytes, cell.k, cell.r, count, inputs.data(), outputs.data())
                == Leopard_Success, "native encode");
            ++calls;
        };
#else
        namespace diag = leopard2_internal;
        Require(!diag::AutoGF16GFNIR19932EnabledForDiagnostics(), "candidate must default OFF");
        leo2_context_options options = {};
        options.struct_size = sizeof(options);
        options.backend = index == 6 ? LEO2_BACKEND_AVX2 : LEO2_BACKEND_AUTO;
        options.thread_count = 1;
        leo2_context* raw_context = NULL;
        Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "context");
        std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)> context(raw_context, leo2_context_destroy);
        Require(leo2_context_backend(context.get()) == LEO2_BACKEND_AVX2 &&
            leo2_context_field_mask(context.get()) == (LEO2_FIELD_MASK_GF8 | LEO2_FIELD_MASK_GF16),
            "AVX2 context and both fields required");
        leo2_codec* raw_codec = NULL;
        Require(leo2_codec_create(context.get(), cell.k, cell.r, LEO2_PROFILE_LEGACY_HIGH_V1,
            index == 8 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16, NULL, &raw_codec) == LEO2_SUCCESS, "codec");
        std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)> codec(raw_codec, leo2_codec_destroy);
        size_t scratch_bytes = 0;
        Require(leo2_encode_scratch_size(codec.get(), cell.bytes, &scratch_bytes) == LEO2_SUCCESS, "scratch query");
        Buffer scratch(scratch_bytes), output(output_bytes);
        uint8_t* const parity = output.data;
        std::vector<void*> outputs(cell.r);
        for (unsigned i = 0; i < cell.r; ++i) outputs[i] = parity + i * cell.bytes;
        const char* const api = index == 1 ? "leo2_encode_batch_one_item" : "leo2_encode";
        leo2_encode_batch_item item = {};
        item.shard_bytes = cell.bytes; item.original = inputs.data(); item.recovery = outputs.data();
        item.scratch = scratch.data; item.scratch_bytes = scratch_bytes;
        const auto selected_gfni = [&](unsigned slot) {
            return (index >= 2 && index <= 4) || (index < 2 && schedule[slot] == '1');
        };
        const auto select = [&](unsigned slot) {
            // No execution or inspection is in flight in this single-threaded driver.
            Require(diag::SetAutoGF16GFNIR19932EnabledForDiagnostics(schedule[slot] == '1'), "select state");
            Require(diag::AutoGF16GFNIEncodeSelectedForDiagnostics(codec.get(), cell.bytes) == selected_gfni(slot),
                "selected route");
            ++selections;
        };
        const auto encode = [&]() {
            Require((index == 1 ? leo2_encode_batch(codec.get(), &item, 1) :
                leo2_encode(codec.get(), cell.bytes, inputs.data(), outputs.data(), scratch.data,
                    scratch_bytes)) == LEO2_SUCCESS, "Leopard2 encode");
            ++calls;
        };
#endif
        const auto check_buffers = [&]() {
            Require(std::memcmp(parity, reference.data, output_bytes) == 0, "paired full parity differs");
            source.Check(); reference.Check(); scratch.Check();
#ifndef LEO_PAIRED_NATIVE
            output.Check();
#endif
        };
        for (unsigned slot = 0; slot < 4; ++slot) {
            select(slot);
#ifndef LEO_PAIRED_NATIVE
            Require(diag::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true), "arm route probe");
#endif
            encode(); ++per_slot[slot];
            if (slot == 0) std::memcpy(reference.data, parity, output_bytes);
            check_buffers();
#ifndef LEO_PAIRED_NATIVE
            probes[slot] = diag::AutoGF16GFNIEncodeCallCountForDiagnostics();
            Require(probes[slot] == static_cast<unsigned>(selected_gfni(slot)), "actual probe route count");
            Require(diag::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics() &&
                diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1, "normalize route probe");
#endif
        }
        // Sample-major order: each pass executes the entire four-state schedule.
        // A group means N complete public API calls, NOT one multi-item batch.
        for (unsigned pass = 0; pass < (exercise ? 25U : 0U); ++pass) {
            for (unsigned slot = 0; slot < 4; ++slot) {
                select(slot);
                // State selection, output comparison and bookkeeping are outside
                // the span. Public result checks, per-call counters and the inner
                // repetition loop are inside, matching the public-call cost model.
                const bool sampled = clocks_enabled && pass >= 4;
                const leopard_paired::Sample sample =
                    leopard_paired::RunGroup(group, sampled, encode, now);
                if (sampled) samples.push_back(sample);
                per_slot[slot] += group;
                check_buffers();
            }
        }
        const unsigned expected = exercise ? 1 + 25 * group : 1;
        for (unsigned slot = 0; slot < 4; ++slot) Require(per_slot[slot] == expected, "slot calls");
        Require(calls == 4 * expected && selections == (exercise ? 104U : 4U), "total accounting");
        Require(samples.size() == (clocks_enabled ? 84U : 0U), "sample accounting");
        Require(Hash(source.data, input_bytes) == input_hash, "input changed");
#ifndef LEO_PAIRED_NATIVE
        Require(diag::AutoGF16GFNIEncodeCallCountForDiagnostics() == probes[3] &&
            diag::AutoGF16GFNIEncodeModeForDiagnostics() == 1, "probe leaked into exercise");
        Require(diag::SetAutoGF16GFNIR19932EnabledForDiagnostics(false) &&
            !diag::AutoGF16GFNIR19932EnabledForDiagnostics(), "restore default OFF");
#endif
        if (argc == 6) {
            FILE* file = std::fopen(argv[5], "wbx");
            Require(file != NULL, "parity file exists or cannot be created");
            const bool written = std::fwrite(parity, 1, output_bytes, file) == output_bytes;
            const int closed = std::fclose(file);
            Require(written && closed == 0, "parity write");
        }
        std::printf("{\"schema\":\"leopard-paired-timer-r19932/v1\",\"codec\":\"%s\",\"cell\":%u,"
            "\"k\":%u,\"r\":%u,\"bytes\":%zu,\"api\":\"%s\",\"schedule\":\"%s\","
            "\"group\":%u,\"warmup_passes\":%u,\"exercise_passes\":%u,\"encode_calls\":%u,"
            "\"selections\":%u,\"per_slot_calls\":[%u,%u,%u,%u],\"probes\":[%u,%u,%u,%u],"
            "\"scratch_bytes\":%zu,\"input_hash\":\"%016llx\",\"output_hash\":\"%016llx\","
            "\"clock_source\":\"%s\",\"samples\":[",
            LEO_PAIRED_CODEC, index, cell.k, cell.r, cell.bytes, api, schedule, group,
            exercise ? 4U : 0U, exercise ? 21U : 0U, calls, selections,
            per_slot[0], per_slot[1], per_slot[2], per_slot[3], probes[0], probes[1], probes[2], probes[3],
            scratch_bytes, static_cast<unsigned long long>(input_hash),
            static_cast<unsigned long long>(Hash(parity, output_bytes)), clock_kind);
        for (size_t i = 0; i < samples.size(); ++i)
            std::printf("%s[%llu,%.8f]", i ? "," : "",
                static_cast<unsigned long long>(samples[i].elapsed_ns), samples[i].ns_per_call);
        std::printf("],\"timed\":%s,\"default_enabled\":false}\n", measured ? "true" : "false");
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
