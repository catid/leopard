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

#include "leopard.h"
#include "leopard2.h"

#include <dirent.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

namespace {

static bool LinuxTaskCount(unsigned& count)
{
    count = 0;
    DIR* directory = opendir("/proc/self/task");
    if (!directory)
        return false;

    errno = 0;
    for (;;)
    {
        const dirent* entry = readdir(directory);
        if (!entry)
            break;
        if (entry->d_name[0] != '.')
            ++count;
    }
    const int read_error = errno;
    const int close_error = closedir(directory);
    return read_error == 0 && close_error == 0;
}

static bool RequireTaskCount(unsigned expected, const char* operation)
{
    unsigned actual = 0;
    if (!LinuxTaskCount(actual))
    {
        fprintf(stderr, "%s: cannot read /proc/self/task\n", operation);
        return false;
    }
    if (actual != expected)
    {
        fprintf(stderr, "%s: task count changed from %u to %u\n",
            operation, expected, actual);
        return false;
    }
    return true;
}

static int TestLegacyInitialization(unsigned initial_tasks)
{
    if (leo_init() != Leopard_Success)
    {
        fprintf(stderr, "legacy leo_init failed\n");
        return 10;
    }
    return RequireTaskCount(initial_tasks, "legacy leo_init") ? 0 : 11;
}

static int TestContextInitialization(bool explicit_single_thread,
                                     unsigned initial_tasks)
{
    leo2_context_options options;
    memset(&options, 0, sizeof(options));
    options.struct_size = sizeof(options);
    options.thread_count = explicit_single_thread ? 1 : 0;

    leo2_context* context = NULL;
    const leo2_result create_result = leo2_context_create(
        explicit_single_thread ? &options : NULL, &context);
    if (create_result != LEO2_SUCCESS)
    {
        fprintf(stderr, "context creation failed: %s\n",
            leo2_result_string(create_result));
        return 20;
    }
    if (!RequireTaskCount(initial_tasks, "context creation"))
    {
        leo2_context_destroy(context);
        return 21;
    }

    const uint32_t effective_threads = leo2_context_thread_count(context);
    if ((explicit_single_thread && effective_threads != 1) ||
        (!explicit_single_thread &&
         (effective_threads == 0 || effective_threads > 128)))
    {
        fprintf(stderr, "unexpected effective thread count: %u\n",
            effective_threads);
        leo2_context_destroy(context);
        return 22;
    }

    leo2_codec* codec = NULL;
    const leo2_result codec_result = leo2_codec_create(
        context, 3, 1, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF8,
        NULL, &codec);
    if (codec_result != LEO2_SUCCESS)
    {
        fprintf(stderr, "codec creation failed: %s\n",
            leo2_result_string(codec_result));
        leo2_context_destroy(context);
        return 23;
    }
    if (!RequireTaskCount(initial_tasks, "codec creation"))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 24;
    }

    if (leo2_encode_batch(codec, NULL, 0) != LEO2_SUCCESS ||
        !RequireTaskCount(initial_tasks, "empty batch"))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 25;
    }

    static const size_t kShardBytes = 64;
    uint8_t original_storage[3][kShardBytes];
    uint8_t recovery_storage[kShardBytes];
    const void* original[3] = {
        original_storage[0], original_storage[1], original_storage[2]
    };
    void* recovery[1] = { recovery_storage };
    for (size_t byte_i = 0; byte_i < kShardBytes; ++byte_i)
    {
        original_storage[0][byte_i] = static_cast<uint8_t>(byte_i * 3u + 1u);
        original_storage[1][byte_i] = static_cast<uint8_t>(byte_i * 5u + 7u);
        original_storage[2][byte_i] = static_cast<uint8_t>(byte_i * 11u + 9u);
        recovery_storage[byte_i] = 0;
    }

    size_t scratch_bytes = 0;
    if (leo2_encode_scratch_size(codec, kShardBytes, &scratch_bytes) !=
            LEO2_SUCCESS)
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 26;
    }
    void* scratch = NULL;
    const size_t allocation_bytes = scratch_bytes == 0 ? 64 : scratch_bytes;
    if (posix_memalign(&scratch, leo2_scratch_alignment(), allocation_bytes) != 0)
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 27;
    }

    leo2_encode_batch_item item;
    item.shard_bytes = kShardBytes;
    item.original = original;
    item.recovery = recovery;
    item.scratch = scratch;
    item.scratch_bytes = scratch_bytes;
    const leo2_result batch_result = leo2_encode_batch(codec, &item, 1);
    bool parity_matches = batch_result == LEO2_SUCCESS;
    for (size_t byte_i = 0; byte_i < kShardBytes; ++byte_i)
    {
        const uint8_t expected = static_cast<uint8_t>(
            original_storage[0][byte_i] ^ original_storage[1][byte_i] ^
            original_storage[2][byte_i]);
        if (recovery_storage[byte_i] != expected)
            parity_matches = false;
    }
    free(scratch);
    if (!parity_matches ||
        !RequireTaskCount(initial_tasks, "single-item batch"))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 28;
    }

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return RequireTaskCount(initial_tasks, "context destruction") ? 0 : 29;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "usage: %s legacy|explicit|default\n", argv[0]);
        return 2;
    }

    // Re-exec once after removing thread-limit variables.  This ensures the
    // OpenMP runtime is loaded in a fresh process with its default behavior,
    // even if it parsed the inherited environment before main().
    static const char* const kCleanEnvironment =
        "LEO2_INITIALIZATION_THREADS_CLEAN_ENV";
    if (!getenv(kCleanEnvironment))
    {
        if (setenv(kCleanEnvironment, "1", 1) != 0 ||
            unsetenv("OMP_NUM_THREADS") != 0 ||
            unsetenv("OMP_THREAD_LIMIT") != 0)
        {
            fprintf(stderr, "cannot prepare clean OpenMP environment\n");
            return 5;
        }
        execv("/proc/self/exe", argv);
        fprintf(stderr, "cannot re-exec test process: %s\n", strerror(errno));
        return 6;
    }

    unsigned initial_tasks = 0;
    if (!LinuxTaskCount(initial_tasks))
    {
        fprintf(stderr, "cannot read initial /proc/self/task count\n");
        return 3;
    }

    if (strcmp(argv[1], "legacy") == 0)
        return TestLegacyInitialization(initial_tasks);
    if (strcmp(argv[1], "explicit") == 0)
        return TestContextInitialization(true, initial_tasks);
    if (strcmp(argv[1], "default") == 0)
        return TestContextInitialization(false, initial_tasks);

    fprintf(stderr, "unknown mode: %s\n", argv[1]);
    return 4;
}
