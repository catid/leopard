// Untimed public-API copy witness for leopard-79h.38.5.4.10.1.
// Reuses the exact .10 workload without changing either codec archive.
// GNU ld --wrap observes external memcpy/memset and fortified references. NOT
// count inline copies, SIMD reads/writes, cache traffic, or execution time.
#define main CurrentRouteCheckMain
#include "current_route_screen.cpp"
#undef main
#include <omp.h>

namespace {
struct CopyCounts
{
    uint64_t input_calls, input_bytes, other_calls, other_bytes;
    uint64_t output_calls, output_bytes, zero_calls, zero_bytes;
    uint64_t other_set_calls, other_set_bytes, encode_calls;
    uint64_t checked_copy_calls, checked_set_calls;
};
CopyCounts counts = {};
volatile bool recording = false;
uintptr_t input_begin = 0, input_end = 0, output_begin = 0, output_end = 0;
unsigned selected_cell = 0;

bool Inside(const void* pointer, size_t bytes, uintptr_t begin, uintptr_t end)
{
    const uintptr_t value = reinterpret_cast<uintptr_t>(pointer);
    return value >= begin && value <= end && bytes <= end - value;
}

void Begin(const void* first_input, void* first_output)
{
    Require(!recording && counts.encode_calls == 0, "one untimed encode only");
    const Cell& cell = kCells[selected_cell];
    input_begin = reinterpret_cast<uintptr_t>(first_input);
    output_begin = reinterpret_cast<uintptr_t>(first_output);
    const size_t inputs = static_cast<size_t>(cell.k) * cell.bytes;
    const size_t outputs = static_cast<size_t>(cell.r) * cell.bytes;
    Require(input_begin <= UINTPTR_MAX - inputs &&
            output_begin <= UINTPTR_MAX - outputs, "slab address overflow");
    input_end = input_begin + inputs;
    output_end = output_begin + outputs;
    ++counts.encode_calls;
    recording = true;
}
}

extern "C" void* __real_memcpy(void*, const void*, size_t);
extern "C" void* __real_memset(void*, int, size_t);
extern "C" void* __real___memcpy_chk(void*, const void*, size_t, size_t);
extern "C" void* __real___memset_chk(void*, int, size_t, size_t);
extern "C" void* CopyChecked(void*, const void*, size_t, size_t) asm("__memcpy_chk");
extern "C" void* SetChecked(void*, int, size_t, size_t) asm("__memset_chk");

static void RecordCopy(void* destination, const void* source, size_t bytes)
{
    if (recording)
    {
        if (Inside(source, bytes, input_begin, input_end))
        {
            ++counts.input_calls;
            counts.input_bytes += bytes;
        }
        else
        {
            ++counts.other_calls;
            counts.other_bytes += bytes;
        }
        if (Inside(destination, bytes, output_begin, output_end))
        {
            ++counts.output_calls;
            counts.output_bytes += bytes;
        }
    }
}

static void RecordSet(int value, size_t bytes)
{
    if (recording)
    {
        if (static_cast<unsigned char>(value) == 0)
        {
            ++counts.zero_calls;
            counts.zero_bytes += bytes;
        }
        else
        {
            ++counts.other_set_calls;
            counts.other_set_bytes += bytes;
        }
    }
}

extern "C" void* __wrap_memcpy(void* destination, const void* source, size_t bytes)
{
    RecordCopy(destination, source, bytes);
    return __real_memcpy(destination, source, bytes);
}

extern "C" void* __wrap_memset(void* destination, int value, size_t bytes)
{
    RecordSet(value, bytes);
    return __real_memset(destination, value, bytes);
}

extern "C" void* __wrap___memcpy_chk(void* destination, const void* source,
    size_t bytes, size_t capacity)
{
    RecordCopy(destination, source, bytes);
    if (recording) ++counts.checked_copy_calls;
    return __real___memcpy_chk(destination, source, bytes, capacity);
}

extern "C" void* __wrap___memset_chk(void* destination, int value,
    size_t bytes, size_t capacity)
{
    RecordSet(value, bytes);
    if (recording) ++counts.checked_set_calls;
    return __real___memset_chk(destination, value, bytes, capacity);
}

#ifdef LEO_CURRENT_SCREEN_MAIN
extern "C" LeopardResult __real_leo_encode(uint64_t, unsigned, unsigned,
    unsigned, const void* const*, void**);
extern "C" LeopardResult __wrap_leo_encode(uint64_t bytes, unsigned k,
    unsigned r, unsigned work_count, const void* const* original, void** work)
{
    Begin(original[0], work[0]);
    const LeopardResult result = __real_leo_encode(bytes, k, r, work_count,
        original, work);
    recording = false;
    return result;
}
#else
extern "C" leo2_result __real_leo2_encode(const leo2_codec*, uint64_t,
    const void* const*, void* const*, void*, size_t);
extern "C" leo2_result __wrap_leo2_encode(const leo2_codec* codec, uint64_t bytes,
    const void* const* original, void* const* recovery, void* scratch,
    size_t scratch_bytes)
{
    Begin(original[0], recovery[0]);
    const leo2_result result = __real_leo2_encode(codec, bytes, original,
        recovery, scratch, scratch_bytes);
    recording = false;
    return result;
}
#endif

int main(int argc, char** argv)
{
    try
    {
        Require((argc == 3 || argc == 4) &&
                std::strcmp(argv[1], "--check") == 0 &&
                std::strlen(argv[2]) == 1 && argv[2][0] >= '0' && argv[2][0] <= '5',
                "copy probe accepts --check cell[0..5] [parity_file] only");
        // Counters are intentionally process-local and single-threaded.
        // Do not rely on an inherited OpenMP environment for this invariant.
        omp_set_dynamic(0);
        omp_set_num_threads(1);
        // Volatile indirect calls also prevent fortified libc inline wrappers
        // from erasing the very symbol calls this self-check must observe.
        // Compile only this driver with -fno-builtin-memcpy/-memset.
        void* (*volatile copy_call)(void*, const void*, size_t) = &std::memcpy;
        void* (*volatile set_call)(void*, int, size_t) = &std::memset;
        void* (*volatile checked_copy_call)(void*, const void*, size_t, size_t)
            = &CopyChecked;
        void* (*volatile checked_set_call)(void*, int, size_t, size_t)
            = &SetChecked;
        uint8_t source[64] = {}, destination[64] = {}, other[64] = {};
        input_begin = reinterpret_cast<uintptr_t>(source);
        input_end = input_begin + sizeof(source);
        output_begin = reinterpret_cast<uintptr_t>(destination);
        output_end = output_begin + sizeof(destination);
        recording = true;
        copy_call(destination, source, sizeof(source));
        copy_call(other, destination, sizeof(other));
        set_call(other, 0, sizeof(other));
        checked_copy_call(destination, source, sizeof(source), sizeof(destination));
        checked_set_call(other, 0, sizeof(other), sizeof(other));
        set_call(other, 257, sizeof(other));
        recording = false;
        Require(counts.input_calls == 2 && counts.input_bytes == 128 &&
                counts.other_calls == 1 && counts.other_bytes == 64 &&
                counts.output_calls == 2 && counts.output_bytes == 128 &&
                counts.zero_calls == 2 && counts.zero_bytes == 128 &&
                counts.other_set_calls == 1 && counts.other_set_bytes == 64 &&
                counts.checked_copy_calls == 1 && counts.checked_set_calls == 1 &&
                counts.encode_calls == 0 && other[63] == 1,
                "linker copy/set wrapper self-check");
        counts = CopyCounts{};
        selected_cell = static_cast<unsigned>(argv[2][0] - '0');
        const int result = CurrentRouteCheckMain(argc, argv);
        Require(result == 0 && !recording && counts.encode_calls == 1,
                "untimed workload or public wrapper failed");
        std::fprintf(stderr,
            "{\"cell\":%u,\"encode_calls\":1,\"input_copy_calls\":%llu,"
            "\"input_copy_bytes\":%llu,\"other_copy_calls\":%llu,"
            "\"other_copy_bytes\":%llu,\"output_copy_calls\":%llu,"
            "\"output_copy_bytes\":%llu,\"zero_calls\":%llu,\"zero_bytes\":%llu,"
            "\"other_set_calls\":%llu,\"other_set_bytes\":%llu,"
            "\"checked_copy_calls\":%llu,\"checked_set_calls\":%llu,"
            "\"external_calls_only\":true,\"timed\":false}\n",
            selected_cell,
            static_cast<unsigned long long>(counts.input_calls),
            static_cast<unsigned long long>(counts.input_bytes),
            static_cast<unsigned long long>(counts.other_calls),
            static_cast<unsigned long long>(counts.other_bytes),
            static_cast<unsigned long long>(counts.output_calls),
            static_cast<unsigned long long>(counts.output_bytes),
            static_cast<unsigned long long>(counts.zero_calls),
            static_cast<unsigned long long>(counts.zero_bytes),
            static_cast<unsigned long long>(counts.other_set_calls),
            static_cast<unsigned long long>(counts.other_set_bytes),
            static_cast<unsigned long long>(counts.checked_copy_calls),
            static_cast<unsigned long long>(counts.checked_set_calls));
        return 0;
    }
    catch (const std::exception& error)
    {
        recording = false;
        std::fprintf(stderr, "copy probe: %s\n", error.what());
        return 1;
    }
}
