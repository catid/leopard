// Bounded correctness checks for the driver-only .38.5.4.11 contrast.
#include "gfni_source_stage_probe.h"
#include "leopard2.h"
#include "leopard.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>
#include <omp.h>

namespace {
void Require(bool value, const char* message)
{
    if (!value) throw std::runtime_error(message);
}

struct Rows
{
    size_t prefix, bytes;
    std::vector<std::vector<uint8_t> > storage;
    std::vector<void*> pointers;
    Rows(unsigned count, size_t length, size_t offset, uint32_t seed)
        : prefix(offset ? offset : (length ? 0 : 17)), bytes(length),
          storage(count), pointers(count)
    {
        for (unsigned row = 0; row < count; ++row)
        {
            storage[row] = std::vector<uint8_t>(prefix + bytes, 0xa5);
            Require(storage[row].capacity() == storage[row].size(),
                    "exact allocation-end fixture");
            for (size_t i = prefix; i < storage[row].size(); ++i)
            {
                seed ^= seed << 13; seed ^= seed >> 17; seed ^= seed << 5;
                storage[row][i] = static_cast<uint8_t>(seed);
            }
            pointers[row] = storage[row].data() + prefix;
        }
    }
    void CheckPrefix() const
    {
        for (const auto& row : storage)
            for (size_t i = 0; i < prefix; ++i)
                Require(row[i] == 0xa5, "row prefix modified");
    }
    uint64_t Hash() const
    {
        uint64_t hash = UINT64_C(14695981039346656037);
        for (const auto& row : storage)
            for (uint8_t value : row)
                hash = (hash ^ value) * UINT64_C(1099511628211);
        return hash;
    }
};

void CheckPredicate()
{
    using namespace gfni_source_stage_probe;
    const Call base = {LEO2_BACKEND_GFNI, 1000, 200, 200, 256, 0,
        32768, 65536, 65536};
    Require(Matches(base), "exact predicate");
    unsigned count = 0;
    unsigned Call::* fields[] = {&Call::kind, &Call::k, &Call::r,
        &Call::requested, &Call::side, &Call::sparse_blocks};
    for (auto member : fields)
    {
        Call bad = base; ++(bad.*member);
        Require(!Matches(bad), "upper predicate neighbor"); ++count;
        bad = base; --(bad.*member);
        Require(!Matches(bad), "lower predicate neighbor"); ++count;
    }
    uint64_t Call::* sizes[] = {&Call::bytes, &Call::source_policy};
    for (auto member : sizes)
    {
        Call bad = base; ++(bad.*member);
        Require(!Matches(bad), "upper size neighbor"); ++count;
        bad = base; --(bad.*member);
        Require(!Matches(bad), "lower size neighbor"); ++count;
    }
    std::printf("predicate: exact match and %u negative neighbors passed\n", count);
}

void CheckKernels()
{
    Require(leo_init() == Leopard_Success, "initialization");
    const auto* gfni = leopard::backend::GetQualifiedOps(LEO2_BACKEND_GFNI);
    const auto* scalar = leopard::backend::GetQualifiedOps(LEO2_BACKEND_SCALAR);
    Require(gfni && scalar, "qualified GFNI and scalar required, no skip");
    const size_t lengths[] = {0, 2, 62, 64, 66, 128, 8190, 8192,
        16384, 16386, 32768, 32770, 65536, 65598};
    unsigned count = 0;
    for (size_t length : lengths)
    for (size_t offset : {size_t(0), size_t(17)})
    for (unsigned zero_mask = 0; zero_mask < 8; ++zero_mask)
    {
        Rows source(4, length, offset, 19 + zero_mask);
        Rows actual(4, length, offset, 31), reference(4, length, offset, 31);
        const uint64_t source_hash = source.Hash();
        uint16_t logs[] = {static_cast<uint16_t>((zero_mask & 1) ? 65535 : 0),
            static_cast<uint16_t>((zero_mask & 2) ? 65535 : 32768),
            static_cast<uint16_t>((zero_mask & 4) ? 65535 : 65534)};
        for (unsigned i = 0; i < 4; ++i)
            std::memcpy(reference.pointers[i], source.pointers[i], length);
        scalar->ff16_ifft_butterfly4(reference.pointers[0], reference.pointers[1],
            reference.pointers[2], reference.pointers[3], logs[0], logs[1], logs[2], length);
        gfni->ff16_ifft_butterfly4_out(source.pointers[0], source.pointers[1],
            source.pointers[2], source.pointers[3], actual.pointers[0],
            actual.pointers[1], actual.pointers[2], actual.pointers[3],
            logs[0], logs[1], logs[2], length);
        Require(actual.storage == reference.storage, "GFNI source kernel vs scalar");
        Require(source.Hash() == source_hash, "kernel source modified");
        source.CheckPrefix(); actual.CheckPrefix(); reference.CheckPrefix();
        ++count;
    }
    std::printf("GFNI first-stage kernel: %u exact-end/unaligned/scalar cases passed\n", count);
}

void CheckPublic(unsigned index)
{
    using namespace gfni_source_stage_probe;
    struct Shape { unsigned k, r; size_t bytes; leo2_backend backend; bool partial; unsigned changes; };
    const Shape shapes[] = {
        {1000, 200, 65536, LEO2_BACKEND_AUTO, false, 2},
        {1000, 200, 65535, LEO2_BACKEND_AUTO, false, 0},
        {1000, 200, 65537, LEO2_BACKEND_AUTO, false, 0},
        {1000, 200, 65536, LEO2_BACKEND_GFNI, false, 2},
        {1000, 200, 65538, LEO2_BACKEND_GFNI, false, 2},
        {1000, 200, 65536, LEO2_BACKEND_AUTO, true, 0},
        {1000, 200, 65536, LEO2_BACKEND_GFNI, true, 0},
        {1000, 199, 65536, LEO2_BACKEND_AUTO, false, 0},
        {999, 200, 65536, LEO2_BACKEND_GFNI, false, 0},
        {1000, 201, 65536, LEO2_BACKEND_GFNI, false, 0},
        {1000, 200, 32768, LEO2_BACKEND_GFNI, false, 0},
        {1000, 200, 65536, LEO2_BACKEND_AVX512, false, 0},
        {1000, 200, 65534, LEO2_BACKEND_AUTO, false, 0},
        {1000, 200, 65538, LEO2_BACKEND_AUTO, false, 0}
    };
    Require(index < sizeof(shapes) / sizeof(shapes[0]), "public case");
    const Shape& shape = shapes[index];
    leo2_context_options options = {};
    options.struct_size = sizeof(options); options.backend = shape.backend; options.thread_count = 1;
    leo2_context* raw_context = NULL;
    Require(leo2_context_create(&options, &raw_context) == LEO2_SUCCESS, "context");
    std::unique_ptr<leo2_context, decltype(&leo2_context_destroy)> context(raw_context, leo2_context_destroy);
    leo2_codec* raw_codec = NULL;
    Require(leo2_codec_create(context.get(), shape.k, shape.r, LEO2_PROFILE_LEGACY_HIGH_V1,
        LEO2_FIELD_GF16, NULL, &raw_codec) == LEO2_SUCCESS, "codec");
    std::unique_ptr<leo2_codec, decltype(&leo2_codec_destroy)> codec(raw_codec, leo2_codec_destroy);
    if ((shape.bytes & 1U) != 0)
    {
        // Native GF16 represents complete symbols. Odd application payloads
        // require the separate padded-odd framing API, not this native route.
        for (bool enabled : {false, true})
        {
            Reset(enabled);
            size_t rejected_scratch = 77;
            Require(leo2_encode_scratch_size(codec.get(), shape.bytes, &rejected_scratch)
                == LEO2_UNSUPPORTED && rejected_scratch == 0, "native odd query rejection");
            Require(leo2_encode(codec.get(), shape.bytes, NULL, NULL, NULL, 0)
                == LEO2_UNSUPPORTED && Get().calls == 0, "native odd execution rejection");
        }
        std::printf("public case %u: native odd B%zu rejected in both modes before transform\n",
            index, shape.bytes);
        return;
    }
    // Read-only inputs may repeat. This bounds ASan allocator-class overhead
    // while every distinct input and each output still has an exact row end.
    Rows source(17, shape.bytes, 17, 20260906);
    Rows control(shape.r, shape.bytes, 17, 29), candidate(shape.r, shape.bytes, 17, 29);
    std::vector<const void*> inputs(shape.k);
    for (unsigned i = 0; i < shape.k; ++i) inputs[i] = source.pointers[i % 17];
    if (shape.partial) control.pointers[17] = candidate.pointers[17] = NULL;
    const uint64_t source_hash = source.Hash();
    const auto omitted = candidate.storage[17];
    size_t scratch_bytes = 0;
    Require(leo2_encode_scratch_size(codec.get(), shape.bytes, &scratch_bytes) == LEO2_SUCCESS,
            "scratch query");
    void* raw_scratch = NULL;
    Require(posix_memalign(&raw_scratch, 64, scratch_bytes) == 0, "scratch allocation");
    std::unique_ptr<void, decltype(&std::free)> scratch(raw_scratch, std::free);
    if (index == 0)
    {
        const uint64_t before = candidate.Hash();
        Reset(true);
        Require(leo2_encode(codec.get(), shape.bytes, inputs.data(), candidate.pointers.data(),
            scratch.get(), scratch_bytes - 1) == LEO2_SCRATCH_TOO_SMALL && Get().calls == 0,
            "short scratch rejected before transform");
        void* saved = candidate.pointers[0]; candidate.pointers[0] = source.pointers[0];
        Require(leo2_encode(codec.get(), shape.bytes, inputs.data(), candidate.pointers.data(),
            scratch.get(), scratch_bytes) == LEO2_OVERLAP && Get().calls == 0,
            "input/output overlap rejected before transform");
        candidate.pointers[0] = saved;
        Require(candidate.Hash() == before && source.Hash() == source_hash,
                "rejected call modified bytes");
    }
    Reset(false);
    std::memset(scratch.get(), 0x5a, scratch_bytes);
    Require(leo2_encode(codec.get(), shape.bytes, inputs.data(), control.pointers.data(),
        scratch.get(), scratch_bytes) == LEO2_SUCCESS, "control encode");
    const State control_state = Get();
    Require(control_state.calls != 0 && control_state.changed == 0 &&
            control_state.matches == shape.changes, "control interception");
    Reset(true);
    std::memset(scratch.get(), 0xa7, scratch_bytes);
    Require(leo2_encode(codec.get(), shape.bytes, inputs.data(), candidate.pointers.data(),
        scratch.get(), scratch_bytes) == LEO2_SUCCESS, "candidate encode");
    Require(Get().calls == control_state.calls && Get().matches == shape.changes &&
            Get().changed == shape.changes, "candidate exact predicate/pass count");
    Require(control.storage == candidate.storage, "public parity/canaries");
    if (shape.partial) Require(candidate.storage[17] == omitted, "omitted output modified");
    Require(source.Hash() == source_hash, "public source modified");
    source.CheckPrefix(); control.CheckPrefix(); candidate.CheckPrefix();
    std::printf("public case %u: K%u/R%u/B%zu backend%d partial%d calls%u changed%u parity%016llx\n",
        index, shape.k, shape.r, shape.bytes, shape.backend, shape.partial,
        Get().calls, Get().changed, static_cast<unsigned long long>(candidate.Hash()));
}
}

int main(int argc, char** argv)
{
    try
    {
        omp_set_dynamic(0); omp_set_num_threads(1);
        Require(argc == 2 || argc == 3, "--kernel or --public [0..13]");
        if (argc == 2 && std::strcmp(argv[1], "--kernel") == 0)
        {
            CheckPredicate(); CheckKernels();
        }
        else
        {
            Require(argc == 3 && std::strcmp(argv[1], "--public") == 0, "public arguments");
            char* end = NULL;
            const unsigned long index = std::strtoul(argv[2], &end, 10);
            Require(argv[2][0] && *end == 0 && index <= 13, "public index");
            CheckPublic(static_cast<unsigned>(index));
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "source-stage test: %s\n", error.what());
        return 1;
    }
}
