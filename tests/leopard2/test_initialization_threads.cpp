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
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <vector>

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

struct ReexecToken
{
    uint64_t magic;
    uint64_t parent_pid;
    uint64_t child_pid;
};

static const uint64_t kReexecMagic = UINT64_C(0x4c454f324f4d5054);

static bool WriteAll(int fd, const void* data, size_t bytes)
{
    const uint8_t* next = static_cast<const uint8_t*>(data);
    while (bytes != 0)
    {
        const ssize_t written = write(fd, next, bytes);
        if (written < 0)
        {
            if (errno == EINTR)
                continue;
            return false;
        }
        if (written == 0)
            return false;
        next += written;
        bytes -= static_cast<size_t>(written);
    }
    return true;
}

static bool ReadAll(int fd, void* data, size_t bytes)
{
    uint8_t* next = static_cast<uint8_t*>(data);
    while (bytes != 0)
    {
        const ssize_t received = read(fd, next, bytes);
        if (received < 0)
        {
            if (errno == EINTR)
                continue;
            return false;
        }
        if (received == 0)
            return false;
        next += received;
        bytes -= static_cast<size_t>(received);
    }
    return true;
}

static int LaunchWithCleanOpenMPEnvironment(const char* mode)
{
    int token_pipe[2];
    if (pipe(token_pipe) != 0)
    {
        fprintf(stderr, "cannot create re-exec token pipe: %s\n",
            strerror(errno));
        return 5;
    }

    const pid_t parent_pid = getpid();
    const pid_t child_pid = fork();
    if (child_pid < 0)
    {
        fprintf(stderr, "cannot fork clean test process: %s\n",
            strerror(errno));
        close(token_pipe[0]);
        close(token_pipe[1]);
        return 6;
    }

    if (child_pid == 0)
    {
        close(token_pipe[1]);
        if (unsetenv("OMP_NUM_THREADS") != 0 ||
            unsetenv("OMP_THREAD_LIMIT") != 0 ||
            unsetenv("OMP_DYNAMIC") != 0)
        {
            fprintf(stderr, "cannot prepare clean OpenMP environment\n");
            _exit(125);
        }

        char fd_string[32];
        snprintf(fd_string, sizeof(fd_string), "%d", token_pipe[0]);
        char* const child_argv[] = {
            const_cast<char*>("leopard2_initialization_threads_test"),
            const_cast<char*>("--clean-child"),
            fd_string,
            const_cast<char*>(mode),
            NULL
        };
        execv("/proc/self/exe", child_argv);
        fprintf(stderr, "cannot re-exec clean test process: %s\n",
            strerror(errno));
        _exit(126);
    }

    close(token_pipe[0]);
    const ReexecToken token = {
        kReexecMagic,
        static_cast<uint64_t>(parent_pid),
        static_cast<uint64_t>(child_pid)
    };
    const bool wrote_token = WriteAll(token_pipe[1], &token, sizeof(token));
    const int close_error = close(token_pipe[1]);
    const bool token_sent = wrote_token && close_error == 0;
    if (!token_sent)
        fprintf(stderr, "cannot send re-exec token: %s\n", strerror(errno));

    int status = 0;
    while (waitpid(child_pid, &status, 0) < 0)
    {
        if (errno == EINTR)
            continue;
        fprintf(stderr, "cannot wait for clean test process: %s\n",
            strerror(errno));
        return 8;
    }
    if (!token_sent)
        return 7;
    if (WIFEXITED(status))
        return WEXITSTATUS(status);
    if (WIFSIGNALED(status))
    {
        fprintf(stderr, "clean test process terminated by signal %d\n",
            WTERMSIG(status));
        return 128 + WTERMSIG(status);
    }
    return 9;
}

static bool ValidateCleanChildToken(const char* fd_text)
{
    char* end = NULL;
    errno = 0;
    const long parsed_fd = strtol(fd_text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || parsed_fd < 0 ||
        parsed_fd > INT_MAX)
        return false;

    const int fd = static_cast<int>(parsed_fd);
    struct stat descriptor_status;
    if (fstat(fd, &descriptor_status) != 0 ||
        !S_ISFIFO(descriptor_status.st_mode))
    {
        close(fd);
        return false;
    }

    // Unlike the former environment marker, entering child mode requires a
    // live kernel pipe inherited from our immediately preceding fork plus the
    // exact PID pair supplied through that one-shot channel.
    ReexecToken token = {};
    const bool received = ReadAll(fd, &token, sizeof(token));
    const int close_error = close(fd);
    return received && close_error == 0 && token.magic == kReexecMagic &&
        token.parent_pid == static_cast<uint64_t>(getppid()) &&
        token.child_pid == static_cast<uint64_t>(getpid());
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

static bool CreateDefaultContext(leo2_context** context_out,
                                 unsigned initial_tasks,
                                 int& error_code)
{
    *context_out = NULL;
    const leo2_result result = leo2_context_create(NULL, context_out);
    if (result != LEO2_SUCCESS)
    {
        fprintf(stderr, "context creation failed: %s\n",
            leo2_result_string(result));
        error_code = 30;
        return false;
    }
    if (!RequireTaskCount(initial_tasks, "GF16 context creation"))
    {
        leo2_context_destroy(*context_out);
        *context_out = NULL;
        error_code = 31;
        return false;
    }
    return true;
}

static int TestGF16HighCodecCreation(unsigned initial_tasks)
{
    leo2_context* context = NULL;
    int error_code = 0;
    if (!CreateDefaultContext(&context, initial_tasks, error_code))
        return error_code;

    leo2_codec* codec = NULL;
    const leo2_result result = leo2_codec_create(
        context, 300, 100, LEO2_PROFILE_LEGACY_HIGH_V1, LEO2_FIELD_GF16,
        NULL, &codec);
    if (result != LEO2_SUCCESS)
    {
        fprintf(stderr, "GF16 high codec creation failed: %s\n",
            leo2_result_string(result));
        leo2_context_destroy(context);
        return 32;
    }
    if (leo2_codec_parent_count(codec) != 512 ||
        leo2_codec_padded_side(codec) != 128 ||
        !RequireTaskCount(initial_tasks, "GF16 high codec creation"))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 33;
    }

    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return RequireTaskCount(initial_tasks, "GF16 high codec destruction") ?
        0 : 34;
}

static int TestGF16DensePlan(bool low_profile, unsigned initial_tasks)
{
    const uint32_t original_count = low_profile ? 200 : 300;
    const uint32_t recovery_count = low_profile ? 300 : 100;
    const uint32_t missing_count = low_profile ? 150 : 100;
    const leo2_profile profile = low_profile ? LEO2_PROFILE_LOW_V1 :
        LEO2_PROFILE_LEGACY_HIGH_V1;

    leo2_context* context = NULL;
    int error_code = 0;
    if (!CreateDefaultContext(&context, initial_tasks, error_code))
        return error_code;

    leo2_codec* codec = NULL;
    const leo2_result codec_result = leo2_codec_create(
        context, original_count, recovery_count, profile, LEO2_FIELD_GF16,
        NULL, &codec);
    if (codec_result != LEO2_SUCCESS)
    {
        fprintf(stderr, "GF16 %s codec creation failed: %s\n",
            low_profile ? "low" : "high", leo2_result_string(codec_result));
        leo2_context_destroy(context);
        return 35;
    }
    if (!RequireTaskCount(initial_tasks, "GF16 dense-plan codec creation"))
    {
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 36;
    }

    // More than the GF16 Walsh/direct crossover is deliberately missing, so
    // this reaches dense active-parent Walsh locator construction rather than
    // the bounded direct-repair or sparse-direct setup paths.
    std::vector<uint8_t> original_present(original_count, 1);
    std::vector<uint8_t> recovery_present(recovery_count, 0);
    for (uint32_t i = 0; i < missing_count; ++i)
    {
        original_present[i] = 0;
        recovery_present[i] = 1;
    }

    leo2_decode_plan* plan = NULL;
    const leo2_result plan_result = leo2_decode_plan_create(codec,
        &original_present[0], &recovery_present[0], &plan);
    if (plan_result != LEO2_SUCCESS)
    {
        fprintf(stderr, "GF16 dense %s plan creation failed: %s\n",
            low_profile ? "low" : "high", leo2_result_string(plan_result));
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 37;
    }
    if (leo2_decode_plan_missing_original_count(plan) != missing_count ||
        !RequireTaskCount(initial_tasks,
            low_profile ? "GF16 dense low plan creation" :
                          "GF16 dense high plan creation"))
    {
        leo2_decode_plan_destroy(plan);
        leo2_codec_destroy(codec);
        leo2_context_destroy(context);
        return 38;
    }

    leo2_decode_plan_destroy(plan);
    leo2_codec_destroy(codec);
    leo2_context_destroy(context);
    return RequireTaskCount(initial_tasks, "GF16 dense plan destruction") ?
        0 : 39;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2)
        return LaunchWithCleanOpenMPEnvironment(argv[1]);

    if (argc != 4 || strcmp(argv[1], "--clean-child") != 0 ||
        !ValidateCleanChildToken(argv[2]))
    {
        fprintf(stderr,
            "usage: %s legacy|explicit|default|gf16-high-codec|"
            "gf16-high-plan|gf16-low-plan\n", argv[0]);
        return 2;
    }
    if (getenv("OMP_NUM_THREADS") || getenv("OMP_THREAD_LIMIT") ||
        getenv("OMP_DYNAMIC"))
    {
        fprintf(stderr, "clean child retained an OpenMP control variable\n");
        return 10;
    }

    unsigned initial_tasks = 0;
    if (!LinuxTaskCount(initial_tasks))
    {
        fprintf(stderr, "cannot read initial /proc/self/task count\n");
        return 3;
    }

    const char* const mode = argv[3];
    if (strcmp(mode, "legacy") == 0)
        return TestLegacyInitialization(initial_tasks);
    if (strcmp(mode, "explicit") == 0)
        return TestContextInitialization(true, initial_tasks);
    if (strcmp(mode, "default") == 0)
        return TestContextInitialization(false, initial_tasks);
    if (strcmp(mode, "gf16-high-codec") == 0)
        return TestGF16HighCodecCreation(initial_tasks);
    if (strcmp(mode, "gf16-high-plan") == 0)
        return TestGF16DensePlan(false, initial_tasks);
    if (strcmp(mode, "gf16-low-plan") == 0)
        return TestGF16DensePlan(true, initial_tasks);

    fprintf(stderr, "unknown mode: %s\n", mode);
    return 4;
}
