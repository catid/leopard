// Untimed endpoint observations only. leopard-79h.38.5.4.19.1.4.1.
// Included after the pinned driver's Buffer definition. No per-call observers.
#ifndef LEOPARD_PAIRED_RUNTIME_METADATA_H
#define LEOPARD_PAIRED_RUNTIME_METADATA_H
#include <dlfcn.h>
#include <link.h>
#include <unistd.h>
#include <limits>

int main(int, char**);
#ifdef LEO_PAIRED_NATIVE
extern "C" LeopardResult __real_leo_encode(uint64_t, unsigned, unsigned, unsigned, const void* const*, void**);
#else
extern "C" leo2_result __real_leo2_encode(const leo2_codec*, uint64_t, const void* const*, void* const*, void*, size_t);
extern "C" leo2_result __real_leo2_encode_batch(const leo2_codec*, const leo2_encode_batch_item*, size_t);
#endif

namespace paired_metadata {
static_assert(sizeof(uintptr_t) == 8 && sizeof(void*) == 8, "qualified ELF64 host required");
const unsigned kInputs = 4096, kOutputs = 1024, kSelections = 104, kSegments = 16;
struct Span { uintptr_t address; size_t bytes; };
struct Allocation { uintptr_t raw, data; size_t bytes; };
struct Segment { uintptr_t address; uint64_t bytes; unsigned flags; };
struct Snapshot {
    Allocation source, reference, scratch, output;
    Span parity, input_array, output_array;
    unsigned input_count, output_count;
    uintptr_t inputs[kInputs], outputs[kOutputs];
    uintptr_t main_address, encode_address, batch_address, anchor_address;
    uintptr_t dladdr_image_base, load_bias;
    unsigned segment_count;
    Segment segments[kSegments];
    long page_bytes;
};
struct Selection { unsigned index, slot; int requested, context, state, gfni; };
struct Store {
    Snapshot snapshots[2];
    Selection selections[kSelections];
    unsigned snapshot_count, selection_count;
};
static Store records = {};

inline uintptr_t End(uintptr_t address, uint64_t bytes)
{
    Require(bytes <= std::numeric_limits<uintptr_t>::max() - address, "metadata span overflow");
    return address + bytes;
}
inline Allocation Describe(const Buffer& buffer)
{
    const uintptr_t raw = reinterpret_cast<uintptr_t>(buffer.raw);
    const uintptr_t data = reinterpret_cast<uintptr_t>(buffer.data);
    Require(raw && data == End(raw, 64) && raw % 64 == 0, "metadata allocation geometry");
    End(End(data, buffer.bytes), 64);
    return Allocation{raw, data, buffer.bytes};
}
// Deliberately named as an anchor, NOT as the main/group loop address.
__attribute__((noinline)) static void Anchor() { asm volatile("" ::: "memory"); }

struct ImageQuery { Snapshot* snapshot; bool found; };
static int FindImage(dl_phdr_info* info, size_t, void* opaque)
{
    ImageQuery& query = *static_cast<ImageQuery*>(opaque);
    Snapshot& s = *query.snapshot;
    const uintptr_t bias = static_cast<uintptr_t>(info->dlpi_addr);
    bool contains = false;
    for (unsigned i = 0; i < info->dlpi_phnum; ++i) {
        const ElfW(Phdr)& p = info->dlpi_phdr[i];
        if (p.p_type != PT_LOAD) continue;
        const uintptr_t begin = End(bias, p.p_vaddr);
        contains |= (p.p_flags & PF_X) && s.main_address >= begin && s.main_address < End(begin, p.p_memsz);
    }
    if (!contains) return 0;
    Require(!query.found, "metadata duplicate main image");
    query.found = true; s.load_bias = bias;
    for (unsigned i = 0; i < info->dlpi_phnum; ++i) {
        const ElfW(Phdr)& p = info->dlpi_phdr[i];
        if (p.p_type != PT_LOAD) continue;
        Require(s.segment_count < kSegments, "metadata segment capacity");
        s.segments[s.segment_count++] = Segment{End(bias, p.p_vaddr), p.p_memsz, p.p_flags};
    }
    return 1;
}
inline void Select(unsigned index, unsigned slot, int requested, int context, int state, int gfni)
{
    Require(index == records.selection_count && index < kSelections && slot == index % 4,
        "metadata selection order/capacity");
    records.selections[records.selection_count++] = Selection{index, slot, requested, context, state, gfni};
}
inline bool Same(const Allocation& a, const Allocation& b)
{
    return a.raw == b.raw && a.data == b.data && a.bytes == b.bytes;
}
inline bool Same(const Span& a, const Span& b) { return a.address == b.address && a.bytes == b.bytes; }
inline bool Same(const Snapshot& a, const Snapshot& b)
{
    if (!Same(a.source,b.source) || !Same(a.reference,b.reference) || !Same(a.scratch,b.scratch) ||
        !Same(a.output,b.output) || !Same(a.parity,b.parity) || !Same(a.input_array,b.input_array) ||
        !Same(a.output_array,b.output_array) || a.input_count != b.input_count || a.output_count != b.output_count ||
        a.main_address != b.main_address || a.encode_address != b.encode_address || a.batch_address != b.batch_address ||
        a.anchor_address != b.anchor_address || a.dladdr_image_base != b.dladdr_image_base ||
        a.load_bias != b.load_bias || a.page_bytes != b.page_bytes || a.segment_count != b.segment_count) return false;
    for (unsigned i = 0; i < a.input_count; ++i) if (a.inputs[i] != b.inputs[i]) return false;
    for (unsigned i = 0; i < a.output_count; ++i) if (a.outputs[i] != b.outputs[i]) return false;
    for (unsigned i = 0; i < a.segment_count; ++i)
        if (a.segments[i].address != b.segments[i].address || a.segments[i].bytes != b.segments[i].bytes ||
            a.segments[i].flags != b.segments[i].flags) return false;
    return true;
}

inline void Capture(unsigned endpoint, const Buffer& source, const Buffer& reference,
    const Buffer& scratch, const Buffer* output, const uint8_t* parity, size_t parity_bytes,
    const std::vector<const void*>& inputs, const std::vector<void*>& outputs, size_t shard_bytes)
{
    Require(endpoint == records.snapshot_count && endpoint < 2, "metadata endpoint order/capacity");
    Require(!inputs.empty() && inputs.size() <= kInputs && !outputs.empty() && outputs.size() <= kOutputs,
        "metadata pointer capacity");
    Snapshot& s = records.snapshots[endpoint];
    s.source = Describe(source); s.reference = Describe(reference); s.scratch = Describe(scratch);
    if (output) s.output = Describe(*output);
    s.parity = Span{reinterpret_cast<uintptr_t>(parity), parity_bytes};
    End(s.parity.address, s.parity.bytes);
    s.input_array = Span{reinterpret_cast<uintptr_t>(inputs.data()), inputs.size() * sizeof(void*)};
    s.output_array = Span{reinterpret_cast<uintptr_t>(outputs.data()), outputs.size() * sizeof(void*)};
    End(s.input_array.address, s.input_array.bytes); End(s.output_array.address, s.output_array.bytes);
    s.input_count = inputs.size(); s.output_count = outputs.size();
    for (unsigned i = 0; i < s.input_count; ++i) {
        s.inputs[i] = reinterpret_cast<uintptr_t>(inputs[i]);
        Require(s.inputs[i] == End(s.source.data, static_cast<uint64_t>(i) * shard_bytes), "metadata input rows");
    }
    for (unsigned i = 0; i < s.output_count; ++i) {
        s.outputs[i] = reinterpret_cast<uintptr_t>(outputs[i]);
        Require(s.outputs[i] == End(s.parity.address, static_cast<uint64_t>(i) * shard_bytes), "metadata output rows");
    }
    s.main_address = reinterpret_cast<uintptr_t>(&main);
    s.anchor_address = reinterpret_cast<uintptr_t>(&Anchor);
#ifdef LEO_PAIRED_NATIVE
    s.encode_address = reinterpret_cast<uintptr_t>(&__real_leo_encode);
#else
    s.encode_address = reinterpret_cast<uintptr_t>(&__real_leo2_encode);
    s.batch_address = reinterpret_cast<uintptr_t>(&__real_leo2_encode_batch);
#endif
    Dl_info image = {};
    Require(dladdr(reinterpret_cast<const void*>(s.main_address), &image) != 0, "metadata dladdr");
    s.dladdr_image_base = reinterpret_cast<uintptr_t>(image.dli_fbase);
    ImageQuery query = {&s, false};
    Require(dl_iterate_phdr(FindImage, &query) == 1 && query.found, "metadata executable image");
    s.page_bytes = sysconf(_SC_PAGESIZE);
    Require(s.page_bytes > 0, "metadata page size");
    ++records.snapshot_count;
    // Compare fields, never padding introduced by aggregate assignments.
    if (endpoint == 1) Require(Same(records.snapshots[0], s), "metadata endpoints differ");
}
inline unsigned long long Number(uintptr_t value) { return static_cast<unsigned long long>(value); }
inline void PrintAllocation(const char* name, const Allocation& a)
{
    std::printf("\"%s\":{\"raw\":%llu,\"data\":%llu,\"bytes\":%zu}", name, Number(a.raw), Number(a.data), a.bytes);
}
inline void PrintSpan(const char* name, const Span& s)
{
    std::printf("\"%s\":{\"address\":%llu,\"bytes\":%zu}", name, Number(s.address), s.bytes);
}
inline void Print(unsigned selections, const unsigned* probes)
{
    Require(records.snapshot_count == 2 && records.selection_count == selections, "metadata inventory");
    std::printf("{\"schema\":\"paired-runtime-metadata/v1\",\"bead\":\"leopard-79h.38.5.4.19.1.4.1\","
        "\"timed\":false,\"observation\":\"new_frontend_endpoints\",\"pointer_bytes\":8,"
        "\"preflight_gfni_counts\":[%u,%u,%u,%u],\"selections\":[", probes[0],probes[1],probes[2],probes[3]);
    for (unsigned i = 0; i < records.selection_count; ++i) {
        const Selection& s = records.selections[i];
        std::printf("%s{\"index\":%u,\"phase\":\"%s\",\"pass\":%d,\"slot\":%u,"
            "\"requested_backend\":%d,\"context_backend\":%d,\"candidate_state\":%d,\"operation_gfni\":%d}",
            i ? "," : "", i, i < 4 ? "preflight" : "exercise", i < 4 ? -1 : static_cast<int>((i-4)/4),
            s.slot, s.requested, s.context, s.state, s.gfni);
    }
    std::printf("],\"native_route_label\":\"%s\",\"snapshots\":[",
#ifdef LEO_PAIRED_NATIVE
        "original_native_compiler_policy"
#else
        "not_native"
#endif
    );
    for (unsigned i = 0; i < 2; ++i) {
        const Snapshot& s = records.snapshots[i];
        std::printf("%s{\"endpoint\":%u,\"allocations\":{", i ? "," : "", i);
        PrintAllocation("source", s.source); std::printf(","); PrintAllocation("reference", s.reference);
        std::printf(","); PrintAllocation("scratch", s.scratch); std::printf(","); PrintAllocation("output", s.output);
        std::printf("},\"spans\":{"); PrintSpan("parity", s.parity); std::printf(",");
        PrintSpan("input_array", s.input_array); std::printf(","); PrintSpan("output_array", s.output_array);
        std::printf("},\"inputs\":[");
        for (unsigned n = 0; n < s.input_count; ++n) std::printf("%s%llu",n ? "," : "",Number(s.inputs[n]));
        std::printf("],\"outputs\":[");
        for (unsigned n = 0; n < s.output_count; ++n) std::printf("%s%llu",n ? "," : "",Number(s.outputs[n]));
        std::printf("],\"functions\":{\"main\":%llu,\"encode\":%llu,\"batch\":%llu,\"metadata_anchor\":%llu},"
            "\"dladdr_image_base\":%llu,\"load_bias\":%llu,\"page_bytes\":%ld,\"segments\":[",
            Number(s.main_address),Number(s.encode_address),Number(s.batch_address),Number(s.anchor_address),
            Number(s.dladdr_image_base),Number(s.load_bias),s.page_bytes);
        for (unsigned n = 0; n < s.segment_count; ++n)
            std::printf("%s{\"address\":%llu,\"bytes\":%llu,\"flags\":%u}", n ? "," : "",
                Number(s.segments[n].address),static_cast<unsigned long long>(s.segments[n].bytes),s.segments[n].flags);
        std::printf("]}");
    }
    std::puts("]}");
}
}
#endif
