// Untimed public integration qualification; leopard-79h.38.5.4.18.4.
// Reuse allocation, canary and codec helpers, but never execute old main.
#define main PreviousBoundaryMain
#include "test_gfni_boundary.cpp"
#undef main
#include "tower_encoder.h"
#include <atomic>
#include <thread>
#include <omp.h>

#if defined(LEO_TOWER_ORIGINAL_CONTROL)
// Link the exact original archive without any experimental module. These
// stubs only resolve unused candidate-check functions; original_shape below
// never calls them. The original encoder has neither overlay nor tower Ops.
namespace tower_encoder {
void SetEnabled(bool) {}
bool TraceAvailable() { return false; }
void ResetCounts() {}
Counts GetCounts() { return {}; }
unsigned InitializationCount() { return 0; }
}
#endif

namespace {
struct Shape { unsigned k, r; size_t bytes, offset; leo2_field field; leo2_backend backend; };
const Shape shapes[] = {
    {1000,200,32768,0,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,199,65536,0,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,200,65536,0,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,199,32768,0,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,200,32770,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,199,65598,17,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,200,16384,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,200,16448,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,200,128,17,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,200,130,17,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,128,32768,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {300,129,32768,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {4096,512,4096,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {1000,200,32768,0,LEO2_FIELD_GF16,LEO2_BACKEND_AUTO},
    {1000,200,32768,0,LEO2_FIELD_GF16,LEO2_BACKEND_GFNI},
    {300,200,16448,1,LEO2_FIELD_GF16,LEO2_BACKEND_SCALAR},
    {17,7,65,1,LEO2_FIELD_GF8,LEO2_BACKEND_AVX2},
    {256,129,16448,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {129,129,16448,17,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {129,129,16450,17,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
    {513,257,16448,1,LEO2_FIELD_GF16,LEO2_BACKEND_AVX2},
};

struct Config {
    Codec::Context context;
    Codec::Handle codec;
    const unsigned k, r;
    Config(unsigned data, unsigned parity, leo2_field field,
           leo2_profile profile = LEO2_PROFILE_LEGACY_HIGH_V1,
           leo2_shard_layout layout = LEO2_SHARD_LAYOUT_NATIVE_V1, unsigned threads = 1)
        : context(nullptr, leo2_context_destroy), codec(nullptr, leo2_codec_destroy), k(data), r(parity)
    {
        leo2_context_options options{};
        options.struct_size = sizeof(options); options.backend = LEO2_BACKEND_AVX2; options.thread_count = threads;
        leo2_context* c = nullptr;
        Require(leo2_context_create(&options, &c) == LEO2_SUCCESS, "special context");
        context.reset(c);
        leo2_codec_options settings{};
        settings.struct_size = sizeof(settings); settings.shard_layout = layout;
        leo2_codec* code = nullptr;
        Require(leo2_codec_create(c, k, r, profile, field, &settings, &code) == LEO2_SUCCESS, "special codec");
        codec.reset(code);
    }
};

struct Buffers {
    Config& config;
    const size_t bytes;
    std::vector<std::unique_ptr<Guard>> source, output;
    std::vector<const void*> inputs;
    std::vector<void*> outputs;
    std::vector<uint64_t> hashes;
    std::unique_ptr<Guard> scratch;
    std::vector<uint8_t> reference;
    Buffers(Config& c, size_t n, unsigned seed, bool padded = false)
        : config(c), bytes(n), inputs(c.k), outputs(c.r), hashes(c.k), reference(size_t(c.r)*n)
    {
        uint32_t random = seed;
        for (unsigned i = 0; i < c.k; ++i) {
            source.emplace_back(new Guard(n, 1));
            Fill(source.back()->data, n, random);
            if (padded) source.back()->data[n-1] = 0;
            inputs[i] = source.back()->data;
            hashes[i] = Hash(source.back()->data, n);
        }
        for (unsigned i = 0; i < c.r; ++i) {
            output.emplace_back(new Guard(n, 17)); outputs[i] = output.back()->data;
        }
        size_t need = 0;
        Require(leo2_encode_scratch_size(c.codec.get(), n, &need) == LEO2_SUCCESS && need, "special scratch");
        scratch.reset(new Guard(need));
    }
    void encode() {
        Require(leo2_encode(config.codec.get(), bytes, inputs.data(), outputs.data(), scratch->data,
                            scratch->bytes) == LEO2_SUCCESS, "special encode");
    }
    void remember() {
        for (unsigned i = 0; i < config.r; ++i) std::memcpy(reference.data()+size_t(i)*bytes, output[i]->data, bytes);
    }
    void poison_outputs() {
        for (auto& row : output) std::memset(row->data,0x5a,bytes);
    }
    void check(bool parity = true) {
        for (unsigned i = 0; i < config.k; ++i) {
            source[i]->Check(); Require(Hash(source[i]->data, bytes) == hashes[i], "special source changed");
        }
        for (unsigned i = 0; i < config.r; ++i) {
            output[i]->Check();
            if (parity && outputs[i]) Require(!std::memcmp(reference.data()+size_t(i)*bytes, output[i]->data, bytes), "special parity");
            if (parity && !outputs[i]) for (size_t j = 0; j < bytes; ++j)
                Require(output[i]->data[j] == 0x5a, "special unrequested output changed");
        }
        scratch->Check();
    }
    void decode(unsigned losses) {
        std::vector<uint8_t> present(config.k, 1), recovery_present(config.r, 1);
        std::vector<const void*> received = inputs, recovery(outputs.begin(), outputs.end());
        std::vector<void*> restored(config.k, nullptr);
        std::vector<std::unique_ptr<Guard>> recovered;
        const unsigned missing[] = {config.k-1, 0, config.k/2};
        for (unsigned i = 0; i < losses; ++i) {
            present[missing[i]] = 0; received[missing[i]] = nullptr;
            recovered.emplace_back(new Guard(bytes, 1)); restored[missing[i]] = recovered.back()->data;
        }
        size_t need = 0;
        Require(leo2_decode_scratch_size(config.codec.get(), bytes, &need) == LEO2_SUCCESS, "decode scratch");
        Guard work(need);
        tower_encoder::ResetCounts();
        Require(leo2_decode(config.codec.get(), bytes, present.data(), recovery_present.data(),
                            received.data(), recovery.data(), restored.data(), work.data, need) == LEO2_SUCCESS,
                "canonical decode roundtrip");
        Require(tower_encoder::GetCounts().selected_passes == 0, "decode entered tower encoder");
        for (unsigned i = 0; i < losses; ++i) {
            recovered[i]->Check();
            Require(!std::memcmp(recovered[i]->data, source[missing[i]]->data, bytes), "restored source differs");
        }
        work.Check(); check();
    }
};

void special(unsigned index)
{
    Require(index < 6, "special index");
    const bool gf8 = index == 3, low = index == 1, padded = index == 0;
    Config c(gf8 ? 17 : 300, gf8 ? 7 : 200, gf8 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16,
             low ? LEO2_PROFILE_LOW_V1 : LEO2_PROFILE_LEGACY_HIGH_V1,
             padded ? LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1 : LEO2_SHARD_LAYOUT_NATIVE_V1, index == 5 ? 2 : 1);
    Buffers first(c, gf8 ? 65 : padded ? 16450 : 16448, 20260910, padded);
    tower_encoder::SetEnabled(false);
    first.encode(); first.remember(); first.check();
    Require(tower_encoder::InitializationCount() == 0, "special cold OFF");
    if (index == 5) {
        Buffers second(c,16450,20260911);
        second.encode(); second.remember(); second.check();
        for (unsigned i = 0; i < c.r; ++i) {
            second.outputs[i] = i == c.r/2 || i == c.r-1 ? second.output[i]->data : nullptr;
            std::memset(second.output[i]->data,0x5a,second.bytes);
        }
        leo2_encode_batch_item items[] = {
            {first.bytes,first.inputs.data(),first.outputs.data(),first.scratch->data,first.scratch->bytes},
            {second.bytes,second.inputs.data(),second.outputs.data(),second.scratch->data,second.scratch->bytes}};
        tower_encoder::SetEnabled(true);
        first.poison_outputs(); second.poison_outputs();
        Require(leo2_encode_batch(c.codec.get(),items,2) == LEO2_SUCCESS,"two-item ordinary batch");
        first.check(); second.check();
        size_t need = 0;
        Require(leo2_encode_batch_preflight_scratch_size(c.codec.get(),2,&need) == LEO2_SUCCESS,"preflight scratch");
        Guard preflight(need);
        first.poison_outputs(); second.poison_outputs();
        Require(leo2_encode_batch_with_preflight_scratch(c.codec.get(),items,2,preflight.data,need) == LEO2_SUCCESS,
                "two-item preflight batch");
        first.check(); second.check(); preflight.Check();
        leo2_encode_batch_binding* raw = nullptr;
        first.poison_outputs(); second.poison_outputs();
        Require(leo2_encode_batch_binding_create(c.codec.get(),items,2,&raw) == LEO2_SUCCESS,"batch binding create");
        std::unique_ptr<leo2_encode_batch_binding,decltype(&leo2_encode_batch_binding_destroy)>
            binding(raw,leo2_encode_batch_binding_destroy);
        for (const auto& row : first.output) for (size_t j = 0; j < first.bytes; ++j)
            Require(row->data[j] == 0x5a,"binding creation encoded first item");
        for (const auto& row : second.output) for (size_t j = 0; j < second.bytes; ++j)
            Require(row->data[j] == 0x5a,"binding creation encoded second item");
        Require(leo2_encode_batch_binding_execute(raw) == LEO2_SUCCESS,"batch binding execute");
        first.check(); second.check();
        // Distinct poison makes premature first-item execution observable,
        // even if a repeated encode would otherwise produce identical bytes.
        for (auto& row : first.output) std::memset(row->data,0xa9,first.bytes);
        for (auto& row : second.output) std::memset(row->data,0xb7,second.bytes);
        std::memset(first.scratch->data,0xc6,first.scratch->bytes);
        std::memset(second.scratch->data,0xd8,second.scratch->bytes);
        const auto a = Hash(first.scratch->data,first.scratch->bytes);
        const auto b = Hash(second.scratch->data,second.scratch->bytes);
        items[1].scratch_bytes -= 1;
        tower_encoder::ResetCounts();
        Require(leo2_encode_batch(c.codec.get(),items,2) == LEO2_SCRATCH_TOO_SMALL,"invalid second batch item");
        Require(Hash(first.scratch->data,first.scratch->bytes) == a &&
                Hash(second.scratch->data,second.scratch->bytes) == b &&
                tower_encoder::GetCounts().selected_passes == 0,"batch atomic rejection");
        first.check(false); second.check(false);
        for (const auto& row : first.output) for (size_t j = 0; j < first.bytes; ++j)
            Require(row->data[j] == 0xa9,"invalid batch changed first output");
        for (const auto& row : second.output) for (size_t j = 0; j < second.bytes; ++j)
            Require(row->data[j] == 0xb7,"invalid batch changed second output");
        Require(tower_encoder::InitializationCount() == 1,"batch cache initialized once");
    } else if (index == 4) {
        Buffers second(c, 16448, 20260911);
        second.encode(); second.remember(); second.check();
        std::atomic<unsigned> ready{0}; std::atomic<bool> go{false};
        std::exception_ptr failures[2]; tower_encoder::Counts worker_counts[2]{};
        tower_encoder::SetEnabled(true);
        auto worker = [&](unsigned id, Buffers& b) {
            omp_set_num_threads(1); omp_set_dynamic(0);
            ready.fetch_add(1);
            while (!go.load()) std::this_thread::yield();
            try {
                tower_encoder::ResetCounts();
                for (unsigned i = 0; i < 3; ++i) { b.poison_outputs(); b.encode(); b.check(); }
                worker_counts[id] = tower_encoder::GetCounts();
            } catch (...) { failures[id] = std::current_exception(); }
        };
        std::thread a(worker, 0, std::ref(first)), b(worker, 1, std::ref(second));
        while (ready.load() != 2) std::this_thread::yield();
        go.store(true); a.join(); b.join();
        for (unsigned i = 0; i < 2; ++i) {
            if (failures[i]) std::rethrow_exception(failures[i]);
            if (tower_encoder::TraceAvailable()) Require(worker_counts[i].selected_passes == 3 &&
                worker_counts[i].source_bytes == UINT64_C(3)*300*16448 &&
                worker_counts[i].output_bytes == UINT64_C(3)*200*16448, "worker-local counts");
        }
        Require(tower_encoder::InitializationCount() == 1, "concurrent initialization exactly once");
        // Warm-cache shared-codec encode/decode overlap. The switch is not
        // changed while either call is active; buffers remain disjoint.
        ready.store(0); go.store(false);
        auto mixed = [&](unsigned id) {
            omp_set_num_threads(1); omp_set_dynamic(0);
            ready.fetch_add(1);
            while (!go.load()) std::this_thread::yield();
            try {
                if (id == 0) { for (unsigned i = 0; i < 3; ++i) { first.poison_outputs(); first.encode(); first.check(); } }
                else { second.decode(1); second.decode(3); }
            } catch (...) { failures[id] = std::current_exception(); }
        };
        std::thread x(mixed, 0), y(mixed, 1);
        while (ready.load() != 2) std::this_thread::yield();
        go.store(true); x.join(); y.join();
        for (const auto& error : failures) if (error) std::rethrow_exception(error);
        first.check(); second.check();
    } else {
        tower_encoder::SetEnabled(true); tower_encoder::ResetCounts();
        first.poison_outputs();
        first.encode(); first.check();
        Require(tower_encoder::InitializationCount() == unsigned(!gf8 && !low), "special selection");
        if (tower_encoder::TraceAvailable()) {
            const auto counts = tower_encoder::GetCounts();
            Require(counts.selected_passes == unsigned(!gf8 && !low), "special pass count");
            Require(counts.source_bytes == (!gf8 && !low ? UINT64_C(300)*16448 : 0), "special source boundary");
        }
        first.decode(1); first.decode(3);
        if (padded) {
            first.source.back()->data[first.bytes-1] = 1;
            const auto source_bad = Hash(first.source.back()->data, first.bytes);
            const auto work_hash = Hash(first.scratch->data, first.scratch->bytes);
            tower_encoder::ResetCounts();
            Require(leo2_encode(c.codec.get(), first.bytes, first.inputs.data(), first.outputs.data(),
                                first.scratch->data, first.scratch->bytes) == LEO2_INVALID_ARGUMENT,
                    "nonzero systematic pad rejection");
            Require(Hash(first.source.back()->data, first.bytes) == source_bad &&
                    Hash(first.scratch->data, first.scratch->bytes) == work_hash &&
                    tower_encoder::GetCounts().selected_passes == 0, "bad pad failure atomicity");
            first.source.back()->data[first.bytes-1] = 0; first.check();
        }
    }
    std::printf("{\"tracker\":\"leopard-79h.38.5.4.18.4\",\"special\":%u,\"initializations\":%u,"
                "\"trace\":%s,\"timed\":false}\n", index, tower_encoder::InitializationCount(),
                tower_encoder::TraceAvailable() ? "true" : "false");
}

void dump(const char* filename, const uint8_t* bytes, size_t count)
{
    FILE* file = std::fopen(filename, "wbx");
    Require(file != nullptr, "exclusive parity file");
    const bool written = std::fwrite(bytes, 1, count, file) == count;
    const int closed = std::fclose(file);
    Require(written && closed == 0, "parity file write");
}

#if defined(LEO_TOWER_ORIGINAL_CONTROL)
void original_shape(unsigned index, const char* filename)
{
    Require(index < sizeof(shapes)/sizeof(shapes[0]) && filename, "original shape arguments");
    const Shape s = shapes[index];
    // First run the original guarded six-mask and invalid-call check. Its
    // separate record is retained, not confused with the parity record below.
    CheckShape(index,s.k,s.r,s.bytes,s.offset,s.field,s.backend,s.backend);
    Codec codec(s.k,s.r,s.field,s.backend);
    Guard source(size_t(s.k)*s.bytes,s.offset), output(size_t(s.r)*s.bytes,s.offset), scratch(codec.Scratch(s.bytes));
    std::vector<const void*> inputs(s.k); std::vector<void*> outputs(s.r);
    uint32_t random = 20260906; Fill(source.data,source.bytes,random);
    for (unsigned i = 0; i < s.k; ++i) inputs[i] = source.data+size_t(i)*s.bytes;
    for (unsigned i = 0; i < s.r; ++i) outputs[i] = output.data+size_t(i)*s.bytes;
    const auto before = Hash(source.data,source.bytes);
    Require(leo2_encode(codec.codec.get(),s.bytes,inputs.data(),outputs.data(),scratch.data,scratch.bytes) == LEO2_SUCCESS,
            "original archive encode");
    Require(Hash(source.data,source.bytes) == before,"original source unchanged");
    source.Check(); output.Check(); scratch.Check(); dump(filename,output.data,output.bytes);
    std::printf("{\"tracker\":\"leopard-79h.38.5.4.18.4\",\"shape\":%u,\"k\":%u,\"r\":%u,\"bytes\":%zu,"
                "\"scratch_bytes\":%zu,\"original\":true,\"trace\":false,\"timed\":false,\"parity_hash\":\"%016llx\"}\n",
                index,s.k,s.r,s.bytes,scratch.bytes,(unsigned long long)Hash(output.data,output.bytes));
}
#endif

void public_shape(unsigned index, const char* parity_file)
{
    Require(index < sizeof(shapes)/sizeof(shapes[0]), "shape index");
    const Shape s = shapes[index];
    Codec codec(s.k, s.r, s.field, s.backend);
    const size_t scratch_bytes = codec.Scratch(s.bytes);
    Guard scratch(scratch_bytes);
    std::vector<std::unique_ptr<Guard>> source, output;
    std::vector<const void*> inputs(s.k);
    std::vector<void*> outputs(s.r);
    std::vector<uint64_t> hashes(s.k);
    uint32_t random = 20260906;
    for (unsigned i = 0; i < s.k; ++i) {
        source.emplace_back(new Guard(s.bytes, s.offset));
        Fill(source.back()->data, s.bytes, random);
        inputs[i] = source.back()->data;
        hashes[i] = Hash(source.back()->data, s.bytes);
    }
    for (unsigned i = 0; i < s.r; ++i) {
        output.emplace_back(new Guard(s.bytes, s.offset));
        outputs[i] = output.back()->data;
    }
    auto encode = [&] {
        Require(leo2_encode(codec.codec.get(), s.bytes, inputs.data(), outputs.data(),
                            scratch.data, scratch_bytes) == LEO2_SUCCESS, "public encode");
    };
    tower_encoder::SetEnabled(false);
    Require(tower_encoder::InitializationCount() == 0, "fresh process cache");
    encode();
    Require(tower_encoder::InitializationCount() == 0, "OFF initialized tower cache");
    Aligned reference(size_t(s.r)*s.bytes);
    for (unsigned i = 0; i < s.r; ++i)
        std::memcpy(reference.data+i*s.bytes, output[i]->data, s.bytes);
    tower_encoder::SetEnabled(true);
    tower_encoder::Counts full{};
    tower_encoder::Counts mask_counts[6]{};
    bool any_selected = false;
    const size_t converted = s.bytes-s.bytes%64;
    // Full, one-prefix, shortened-prefix, sparse, alternating, none. Every
    // selected call must leave its public outputs in canonical coordinates.
    for (unsigned mask = 0; mask < 6; ++mask) {
        unsigned prefix = 0;
        for (unsigned i = 0; i < s.r; ++i) {
            const bool selected = mask == 0 || (mask == 1 && i == 0) ||
                (mask == 2 && i < s.r-1) || (mask == 3 && (i == s.r/2 || i == s.r-1)) ||
                (mask == 4 && i%2 == 0);
            outputs[i] = selected ? output[i]->data : nullptr;
            if (selected) prefix = i+1;
            std::memset(output[i]->data, 0x5a, s.bytes);
        }
        std::memset(scratch.data, 0xa6, scratch_bytes);
        tower_encoder::ResetCounts();
        encode();
        const auto counts = tower_encoder::GetCounts();
        mask_counts[mask] = counts;
        if (mask == 0) full = counts;
        // These high-profile production shapes do not use sparse-plan
        // execution. AUTO's qualified full-output GFNI cell falls back to
        // the context AVX2 backend for every nonempty partial mask.
        const bool selected = prefix != 0 && s.field == LEO2_FIELD_GF16 &&
            (s.backend == LEO2_BACKEND_AVX2 || (s.backend == LEO2_BACKEND_AUTO && mask != 0)) &&
            s.r > 128 && converted > 16384;
        any_selected |= selected;
        Require(tower_encoder::InitializationCount() == unsigned(any_selected), "per-mask cache initialization");
        if (tower_encoder::TraceAvailable()) {
            Require((counts.selected_passes != 0) == selected, "per-mask route eligibility");
            Require(counts.source_bytes == (selected ? converted*s.k : 0) &&
                    counts.output_bytes == (selected ? converted*prefix : 0) &&
                    counts.source_rows == counts.selected_passes*s.k &&
                    counts.output_rows == counts.selected_passes*prefix, "per-mask conversion boundaries");
        }
        for (unsigned i = 0; i < s.r; ++i) {
            if (outputs[i]) Require(!std::memcmp(output[i]->data, reference.data+i*s.bytes, s.bytes),
                                    "canonical output parity differs");
            else for (size_t j = 0; j < s.bytes; ++j)
                Require(output[i]->data[j] == 0x5a, "unrequested output changed");
            output[i]->Check();
        }
        for (unsigned i = 0; i < s.k; ++i) {
            source[i]->Check();
            Require(Hash(source[i]->data, s.bytes) == hashes[i], "source changed");
        }
        scratch.Check();
        if (mask == 5 && tower_encoder::TraceAvailable())
            Require(counts.selected_passes == 0, "zero outputs entered tower path");
    }
    // Restore output pointers before testing failure atomicity.
    for (unsigned i = 0; i < s.r; ++i) outputs[i] = output[i]->data;
    std::memset(scratch.data, 0xa6, scratch_bytes);
    const uint64_t scratch_hash = Hash(scratch.data, scratch_bytes);
    tower_encoder::ResetCounts();
    Require(leo2_encode(codec.codec.get(), s.bytes, inputs.data(), outputs.data(), scratch.data,
                        scratch_bytes-1) == LEO2_SCRATCH_TOO_SMALL, "short scratch not rejected");
    if (s.field == LEO2_FIELD_GF16)
        Require(leo2_encode(codec.codec.get(), s.bytes-1, inputs.data(), outputs.data(), scratch.data,
                            scratch_bytes) == LEO2_UNSUPPORTED, "odd GF16 not rejected");
    Require(Hash(scratch.data, scratch_bytes) == scratch_hash, "size rejection mutated scratch");
    // Payload-overlap rejection can use the caller's scratch range tables.
    // Only metadata-overlap rejection promises untouched scratch. Preserve
    // public source/output bytes and guards, not temporary scratch contents.
    outputs[0] = source[0]->data;
    Require(leo2_encode(codec.codec.get(),s.bytes,inputs.data(),outputs.data(),scratch.data,scratch_bytes) == LEO2_OVERLAP,
            "source/output overlap not rejected");
    outputs[0] = output[0]->data;
    outputs[1] = output[0]->data;
    Require(leo2_encode(codec.codec.get(),s.bytes,inputs.data(),outputs.data(),scratch.data,scratch_bytes) == LEO2_OVERLAP,
            "duplicate output not rejected");
    outputs[1] = output[1]->data;
    inputs[0] = scratch.data;
    Require(leo2_encode(codec.codec.get(),s.bytes,inputs.data(),outputs.data(),scratch.data,scratch_bytes) == LEO2_OVERLAP,
            "source/scratch overlap not rejected");
    inputs[0] = source[0]->data;
    Require(tower_encoder::GetCounts().selected_passes == 0, "invalid call entered tower");
    for (unsigned i = 0; i < s.k; ++i) {
        source[i]->Check();
        Require(Hash(source[i]->data, s.bytes) == hashes[i], "invalid call mutated source");
    }
    for (unsigned i = 0; i < s.r; ++i) {
        for (size_t j = 0; j < s.bytes; ++j)
            Require(output[i]->data[j] == 0x5a, "invalid call mutated output");
        output[i]->Check();
    }
    scratch.Check();
    // Explicit expected full-route region, independent of the implementation
    // Select helper. The AUTO test is a previously promoted GFNI cell.
    const bool eligible = s.field == LEO2_FIELD_GF16 && s.backend == LEO2_BACKEND_AVX2 &&
        s.r > 128 && s.bytes-s.bytes%64 > 16384;
    if (tower_encoder::TraceAvailable()) {
        Require((full.selected_passes != 0) == eligible, "full route eligibility");
        if (eligible) {
            Require(full.source_bytes == converted*s.k && full.output_bytes == converted*s.r &&
                    full.source_rows == full.selected_passes*s.k &&
                    full.output_rows == full.selected_passes*s.r, "public conversion boundary counts");
            unsigned side = 1;
            while (side < s.r) side *= 2;
            Require(full.inverse_pairs && full.forward_pairs &&
                    bool(full.accumulating_pairs) == (s.k > side),
                    "full transform phases not observed");
        } else Require(full.source_bytes == 0 && full.output_bytes == 0, "excluded conversion");
    }
    Require(tower_encoder::InitializationCount() == unsigned(any_selected), "cache initialization count");
    // The output file is the successful OFF canonical oracle; a fresh ON
    // full encode must match it byte-for-byte before that file is emitted.
    encode();
    for (unsigned i = 0; i < s.r; ++i)
        Require(!std::memcmp(output[i]->data, reference.data+i*s.bytes, s.bytes), "final ON canonical parity");
    // Actual public one-item batch, full and a hole-containing output mask.
    for (unsigned partial = 0; partial < 2; ++partial) {
        for (unsigned i = 0; i < s.r; ++i) {
            outputs[i] = !partial || i == s.r/2 || i == s.r-1 ? output[i]->data : nullptr;
            std::memset(output[i]->data, 0x5a, s.bytes);
        }
        leo2_encode_batch_item item{s.bytes, inputs.data(), outputs.data(), scratch.data, scratch_bytes};
        Require(leo2_encode_batch(codec.codec.get(), &item, 1) == LEO2_SUCCESS, "one-item batch");
        for (unsigned i = 0; i < s.r; ++i) {
            if (outputs[i]) Require(!std::memcmp(output[i]->data, reference.data+i*s.bytes, s.bytes), "batch parity");
            else for (size_t j = 0; j < s.bytes; ++j) Require(output[i]->data[j] == 0x5a, "batch unrequested output");
            output[i]->Check();
        }
        for (unsigned i = 0; i < s.k; ++i) {
            source[i]->Check();
            Require(Hash(source[i]->data, s.bytes) == hashes[i], "batch source changed");
        }
        scratch.Check();
    }
    if (parity_file) dump(parity_file, reference.data, size_t(s.r)*s.bytes);
    std::printf("{\"tracker\":\"leopard-79h.38.5.4.18.4\",\"shape\":%u,\"k\":%u,\"r\":%u,"
                "\"bytes\":%zu,\"offset\":%zu,\"field\":%u,\"backend\":%u,\"masks\":6,"
                "\"scratch_bytes\":%zu,\"selected_passes\":%llu,\"source_rows\":%llu,"
                "\"source_bytes\":%llu,\"output_rows\":%llu,\"output_bytes\":%llu,"
                "\"inverse_pairs\":%llu,\"forward_pairs\":%llu,\"accumulating_pairs\":%llu,"
                "\"initializations\":%u,\"trace\":%s,\"parity_hash\":\"%016llx\",\"batch_calls\":2,\"timed\":false,\"mask_counts\":[",
                index,s.k,s.r,s.bytes,s.offset,unsigned(s.field),unsigned(s.backend),scratch_bytes,
                (unsigned long long)full.selected_passes,(unsigned long long)full.source_rows,
                (unsigned long long)full.source_bytes,(unsigned long long)full.output_rows,
                (unsigned long long)full.output_bytes,(unsigned long long)full.inverse_pairs,
                (unsigned long long)full.forward_pairs,(unsigned long long)full.accumulating_pairs,
                tower_encoder::InitializationCount(),tower_encoder::TraceAvailable()?"true":"false",
                (unsigned long long)Hash(reference.data,size_t(s.r)*s.bytes));
    for (unsigned mask = 0; mask < 6; ++mask) {
        const auto c = mask_counts[mask];
        std::printf("%s{\"passes\":%llu,\"source_bytes\":%llu,\"output_bytes\":%llu}", mask ? "," : "",
                    (unsigned long long)c.selected_passes, (unsigned long long)c.source_bytes,
                    (unsigned long long)c.output_bytes);
    }
    std::puts("]}");
}
} // namespace

int main(int argc, char** argv)
{
    try {
        Require((argc == 3 || argc == 4) && (!std::strcmp(argv[1], "--shape") ||
                (argc == 3 && !std::strcmp(argv[1], "--special"))),
                "usage: --shape INDEX [exclusive_parity_file] or --special INDEX; no timing supported");
        unsigned index = 0;
        Require(*argv[2] && std::strlen(argv[2]) <= 2, "shape syntax");
        for (const char* p = argv[2]; *p; ++p) {
            Require(*p >= '0' && *p <= '9', "shape syntax");
            index = index*10 + unsigned(*p-'0');
        }
#if defined(LEO_TOWER_ORIGINAL_CONTROL)
        Require(argc == 4 && !std::strcmp(argv[1], "--shape"), "original shape-only control");
        original_shape(index,argv[3]);
#else
        if (!std::strcmp(argv[1], "--special")) special(index);
        else public_shape(index, argc == 4 ? argv[3] : nullptr);
#endif
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what()); return 1;
    }
}
