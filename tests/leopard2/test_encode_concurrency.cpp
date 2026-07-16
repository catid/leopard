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

#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if defined(_MSC_VER)
#include <malloc.h>
#endif

namespace {

class AlignedBuffer
{
public:
    explicit AlignedBuffer(size_t bytes)
        : data_(NULL)
        , bytes_(bytes)
    {
        if (bytes == 0)
            return;
#if defined(_MSC_VER)
        data_ = _aligned_malloc(bytes, leo2_scratch_alignment());
#else
        if (posix_memalign(&data_, leo2_scratch_alignment(), bytes) != 0)
            data_ = NULL;
#endif
        if (!data_)
            throw std::bad_alloc();
        memset(data_, 0, bytes);
    }

    ~AlignedBuffer()
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
    AlignedBuffer(const AlignedBuffer&);
    AlignedBuffer& operator=(const AlignedBuffer&);

    void* data_;
    size_t bytes_;
};

void require(bool condition, const char* operation)
{
    if (!condition)
        throw std::runtime_error(operation);
}

void require_result(leo2_result result, const char* operation)
{
    if (result != LEO2_SUCCESS)
    {
        throw std::runtime_error(
            std::string(operation) + ": " + leo2_result_string(result));
    }
}

typedef std::vector<std::vector<uint8_t> > Storage;

void run_profile(
    leo2_context* context,
    uint32_t k,
    uint32_t r,
    leo2_profile profile,
    leo2_field field,
    leo2_shard_layout layout,
    size_t shard_bytes,
    uint64_t& compared_bytes,
    uint64_t& executions)
{
    const unsigned thread_count = 8;
    const unsigned repetitions = 16;
    const unsigned batch_count = 4;

    leo2_codec* codec = NULL;
    leo2_codec_options codec_options;
    memset(&codec_options, 0, sizeof(codec_options));
    codec_options.struct_size = sizeof(codec_options);
    codec_options.shard_layout = layout;
    require_result(leo2_codec_create(context, k, r, profile, field,
        &codec_options, &codec), "concurrent codec create");

    Storage original_storage(k, std::vector<uint8_t>(shard_bytes + 2, 0));
    std::vector<const void*> original(k);
    for (uint32_t shard = 0; shard < k; ++shard)
    {
        uint8_t* data = &original_storage[shard][1];
        for (size_t byte_i = 0; byte_i < shard_bytes; ++byte_i)
        {
            data[byte_i] = static_cast<uint8_t>(
                shard * 31u + byte_i * 17u + (byte_i >> 3));
        }
        if (layout == LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1)
            data[shard_bytes - 1] = 0;
        original[shard] = data;
    }

    size_t scratch_bytes = 0;
    require_result(leo2_encode_scratch_size(codec, shard_bytes, &scratch_bytes),
        "concurrent scratch query");
    AlignedBuffer reference_scratch(scratch_bytes);
    Storage expected(r, std::vector<uint8_t>(shard_bytes + 2, 0));
    std::vector<void*> expected_output(r);
    for (uint32_t shard = 0; shard < r; ++shard)
        expected_output[shard] = &expected[shard][1];
    require_result(leo2_encode(codec, shard_bytes, &original[0],
        &expected_output[0], reference_scratch.data(), reference_scratch.size()),
        "concurrent reference encode");
    for (uint32_t shard = 0; shard < r; ++shard)
        require(expected[shard].front() == 0 && expected[shard].back() == 0,
            "serial encode modified an output guard");

    std::vector<Storage> batch_output(batch_count,
        Storage(r, std::vector<uint8_t>(shard_bytes + 2, 0x5a)));
    std::vector<std::vector<void*> > batch_output_pointers(batch_count,
        std::vector<void*>(r, NULL));
    std::vector<std::unique_ptr<AlignedBuffer> > batch_scratch(batch_count);
    std::vector<leo2_encode_batch_item> batch_items(batch_count);
    for (unsigned item_i = 0; item_i < batch_count; ++item_i)
    {
        for (uint32_t shard = 0; shard < r; ++shard)
            batch_output_pointers[item_i][shard] =
                &batch_output[item_i][shard][1];
        batch_scratch[item_i].reset(new AlignedBuffer(scratch_bytes));
        batch_items[item_i].shard_bytes = shard_bytes;
        batch_items[item_i].original = &original[0];
        batch_items[item_i].recovery = &batch_output_pointers[item_i][0];
        batch_items[item_i].scratch = batch_scratch[item_i]->data();
        batch_items[item_i].scratch_bytes = batch_scratch[item_i]->size();
    }
    require_result(leo2_encode_batch(codec, &batch_items[0], batch_items.size()),
        "concurrent-profile batch encode");
    for (unsigned item_i = 0; item_i < batch_count; ++item_i)
        for (uint32_t shard = 0; shard < r; ++shard)
        {
            require(memcmp(&batch_output[item_i][shard][1],
                &expected[shard][1], shard_bytes) == 0,
                "batch encode differed from serial reference");
            require(batch_output[item_i][shard].front() == 0x5a &&
                batch_output[item_i][shard].back() == 0x5a,
                "batch encode modified an output guard");
        }

    std::atomic<unsigned> ready(0);
    std::atomic<bool> go(false);
    std::atomic<bool> failed(false);
    std::vector<std::thread> threads;
    threads.reserve(thread_count);
    for (unsigned thread_i = 0; thread_i < thread_count; ++thread_i)
    {
        threads.push_back(std::thread([&, thread_i]() {
            ready.fetch_add(1, std::memory_order_release);
            while (!go.load(std::memory_order_acquire))
                std::this_thread::yield();
            try
            {
                Storage output(r,
                    std::vector<uint8_t>(shard_bytes + 2, 0xa5));
                std::vector<void*> output_ptrs(r);
                for (uint32_t shard = 0; shard < r; ++shard)
                    output_ptrs[shard] = &output[shard][1];
                AlignedBuffer scratch(scratch_bytes);

                for (unsigned repetition = 0; repetition < repetitions;
                     ++repetition)
                {
                    const leo2_result result = leo2_encode(codec, shard_bytes,
                        &original[0], &output_ptrs[0], scratch.data(),
                        scratch.size());
                    if (result != LEO2_SUCCESS)
                    {
                        failed.store(true, std::memory_order_relaxed);
                        return;
                    }
                    for (uint32_t shard = 0; shard < r; ++shard)
                    {
                        if (memcmp(&output[shard][1], &expected[shard][1],
                                shard_bytes) != 0 ||
                            output[shard].front() != 0xa5 ||
                            output[shard].back() != 0xa5)
                        {
                            failed.store(true, std::memory_order_relaxed);
                            return;
                        }
                    }
                }
                (void)thread_i;
            }
            catch (...)
            {
                failed.store(true, std::memory_order_relaxed);
            }
        }));
    }
    while (ready.load(std::memory_order_acquire) != thread_count)
        std::this_thread::yield();
    go.store(true, std::memory_order_release);
    for (size_t i = 0; i < threads.size(); ++i)
        threads[i].join();
    require(!failed.load(std::memory_order_relaxed),
        "concurrent encode differed from serial reference");
    for (uint32_t shard = 0; shard < k; ++shard)
        require(original_storage[shard].front() == 0 &&
            original_storage[shard].back() == 0,
            "encode modified an input guard");

    executions += static_cast<uint64_t>(thread_count) * repetitions;
    executions += batch_count;
    compared_bytes += static_cast<uint64_t>(thread_count) * repetitions * r *
        shard_bytes;
    compared_bytes += static_cast<uint64_t>(batch_count) * r * shard_bytes;
    leo2_codec_destroy(codec);
}

} // namespace

int main()
{
    try
    {
        leo2_context_options options;
        memset(&options, 0, sizeof(options));
        options.struct_size = sizeof(options);
        options.backend = LEO2_BACKEND_AUTO;
        options.thread_count = 4;
        leo2_context* context = NULL;
        require_result(leo2_context_create(&options, &context),
            "concurrent context create");

        uint64_t compared_bytes = 0;
        uint64_t executions = 0;
        run_profile(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1, 257,
            compared_bytes, executions);
        run_profile(context, 5, 11, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF8, LEO2_SHARD_LAYOUT_NATIVE_V1, 257,
            compared_bytes, executions);
        run_profile(context, 9, 7, LEO2_PROFILE_LEGACY_HIGH_V1,
            LEO2_FIELD_GF16, LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 66,
            compared_bytes, executions);
        run_profile(context, 5, 11, LEO2_PROFILE_LOW_V1,
            LEO2_FIELD_GF16, LEO2_SHARD_LAYOUT_GF16_PADDED_ODD_V1, 66,
            compared_bytes, executions);

        leo2_context_destroy(context);
        std::cout << "leopard2 concurrent encode passed: profiles=4"
                  << " executions=" << executions
                  << " compared_bytes=" << compared_bytes << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "leopard2 concurrent encode failed: "
                  << error.what() << std::endl;
        return 1;
    }
}
