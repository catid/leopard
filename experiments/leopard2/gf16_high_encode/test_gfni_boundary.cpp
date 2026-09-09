// Focused public-API bounds and differential checks; never reads a clock.
#define LEO_BOUNDARY_NO_MAIN 1
#include "gfni_boundary_screen.cpp"
#if defined(__SANITIZE_ADDRESS__)
#include <sanitizer/asan_interface.h>
#endif

namespace {
struct Guard
{
    Aligned storage;
    size_t bytes, offset;
    uint8_t* data;
    Guard(size_t size, size_t misalignment = 0)
        : storage(size + 192), bytes(size), offset(64 + misalignment),
          data(storage.data + offset)
    {
        std::memset(storage.data, 0xd3, bytes + 192);
        Poison(true);
    }
    ~Guard() { Poison(false); }
    void Poison(bool enable)
    {
#if defined(__SANITIZE_ADDRESS__)
        if (enable)
        {
            __asan_poison_memory_region(storage.data, offset);
            __asan_poison_memory_region(data + bytes, 192 - offset);
        }
        else __asan_unpoison_memory_region(storage.data, bytes + 192);
#else
        (void)enable;
#endif
    }
    void Check()
    {
        Poison(false);
        bool valid = true;
        for (size_t i = 0; i < offset; ++i) valid &= storage.data[i] == 0xd3;
        for (size_t i = offset + bytes; i < bytes + 192; ++i)
            valid &= storage.data[i] == 0xd3;
        Poison(true);
        Require(valid, "guard overwritten");
    }
};

void Check(unsigned cell, leo2_backend candidate_backend = LEO2_BACKEND_GFNI)
{
    Require(cell < 8, "guarded cell");
    const unsigned k = cell < 6 ? 1000 : 17;
    const unsigned r = cell < 6 ? (cell % 2 ? 199 : 200) : 7;
    const size_t bytes = cell < 6 ? (cell % 2 ? 65536 : 32768) +
        (cell >= 4 ? 2 : 0) : (cell == 6 ? 65 : 66);
    const size_t offset = cell < 2 ? 0 : cell < 4 ? 1 : cell < 6 ? 2 : 1;
    const leo2_field field = cell == 6 ? LEO2_FIELD_GF8 : LEO2_FIELD_GF16;
    Codec baseline(k, r, field, LEO2_BACKEND_AVX2);
    Codec candidate(k, r, field, candidate_backend);
    Require(leo2_context_field_mask(candidate.context.get()) ==
        (LEO2_FIELD_MASK_GF8 | LEO2_FIELD_MASK_GF16), "both fields required");
    const size_t scratch_bytes = candidate.Scratch(bytes);
    Require(scratch_bytes == baseline.Scratch(bytes) && scratch_bytes > 0,
        "scratch geometry changed");
    Guard scratch(scratch_bytes);
    std::vector<std::unique_ptr<Guard> > source, output;
    std::vector<const void*> inputs(k);
    std::vector<void*> outputs(r);
    std::vector<uint64_t> source_hashes(k);
    uint32_t random = 20260906;
    for (unsigned i = 0; i < k; ++i)
    {
        source.emplace_back(new Guard(bytes, offset));
        Fill(source.back()->data, bytes, random);
        inputs[i] = source.back()->data;
        source_hashes[i] = Hash(source.back()->data, bytes);
    }
    for (unsigned i = 0; i < r; ++i)
    {
        output.emplace_back(new Guard(bytes, offset));
        outputs[i] = output.back()->data;
    }
    Require(leo2_encode(baseline.codec.get(), bytes, inputs.data(), outputs.data(),
        scratch.data, scratch_bytes) == LEO2_SUCCESS, "baseline full encode");
    Aligned reference(static_cast<size_t>(r) * bytes);
    for (unsigned i = 0; i < r; ++i)
        std::memcpy(reference.data + i * bytes, output[i]->data, bytes);
    // Full, prefix-1, prefix-(R-1), sparse-middle/last, alternating, no outputs.
    for (unsigned mask = 0; mask < 6; ++mask)
    {
        for (unsigned i = 0; i < r; ++i)
        {
            const bool selected = mask == 0 || (mask == 1 && i == 0) ||
                (mask == 2 && i < r - 1) || (mask == 3 && (i == r / 2 || i == r - 1)) ||
                (mask == 4 && i % 2 == 0);
            outputs[i] = selected ? output[i]->data : NULL;
            std::memset(output[i]->data, 0x5a, bytes);
        }
        Require(leo2_encode(candidate.codec.get(), bytes, inputs.data(), outputs.data(),
            scratch.data, scratch_bytes) == LEO2_SUCCESS, "GFNI full/subset encode");
        for (unsigned i = 0; i < r; ++i)
        {
            if (outputs[i]) Require(std::memcmp(output[i]->data,
                reference.data + i * bytes, bytes) == 0, "GFNI parity differs");
            else for (size_t j = 0; j < bytes; ++j)
                Require(output[i]->data[j] == 0x5a, "unrequested output changed");
            output[i]->Check();
        }
        for (unsigned i = 0; i < k; ++i)
        {
            source[i]->Check();
            Require(Hash(source[i]->data, bytes) == source_hashes[i], "input changed");
        }
        scratch.Check();
    }
    for (unsigned i = 0; i < r; ++i) outputs[i] = output[i]->data;
    std::memset(scratch.data, 0xa6, scratch_bytes);
    const uint64_t scratch_hash = Hash(scratch.data, scratch_bytes);
    Require(leo2_encode(candidate.codec.get(), bytes, inputs.data(), outputs.data(),
        scratch.data, scratch_bytes - 1) == LEO2_SCRATCH_TOO_SMALL,
        "short scratch not rejected");
    if (field == LEO2_FIELD_GF16)
        Require(leo2_encode(candidate.codec.get(), bytes - 1, inputs.data(), outputs.data(),
            scratch.data, scratch_bytes) == LEO2_UNSUPPORTED, "odd GF16 accepted");
    Require(Hash(scratch.data, scratch_bytes) == scratch_hash, "failed call changed scratch");
    for (unsigned i = 0; i < r; ++i)
    {
        for (size_t j = 0; j < bytes; ++j)
            Require(output[i]->data[j] == 0x5a, "failed call changed output");
        output[i]->Check();
    }
    scratch.Check();
    std::printf("{\"schema\":\"leopard-gfni-boundary-guards/v1\",\"cell\":%u,"
        "\"k\":%u,\"r\":%u,\"bytes\":%zu,\"misalignment\":%zu,\"field\":%u,"
        "\"subset_masks\":6,\"scratch_bytes\":%zu,\"timed\":false}\n",
        cell, k, r, bytes, offset, static_cast<unsigned>(field), scratch_bytes);
}
}

#ifndef LEO_BOUNDARY_GUARD_NO_MAIN
int main(int argc, char** argv)
{
    try
    {
        Require(argc == 2 && std::strlen(argv[1]) == 1 && argv[1][0] >= '0' &&
            argv[1][0] <= '7', "usage: test_gfni_boundary cell[0..7]");
        Check(static_cast<unsigned>(argv[1][0] - '0'));
        return 0;
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "%s\n", error.what());
        return 1;
    }
}
#endif
