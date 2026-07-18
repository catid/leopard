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

#pragma once

/*
 * Clang's ThreadSanitizer runtime supplies strong replacement new/delete
 * definitions.  An executable cannot also install the allocation-audit
 * replacements used by a few Leopard2 tests.  Detect instrumentation from the
 * compiler rather than from CMake flag strings: this covers global,
 * configuration-specific, environment-provided, and target-local flags.
 *
 * ASan does not have this restriction, so allocation auditing remains enabled
 * there as well as in ordinary Release builds.
 */
#if defined(__has_feature)
# if __has_feature(thread_sanitizer)
#  define LEO2_TEST_THREAD_SANITIZED 1
# endif
#endif

#if defined(__SANITIZE_THREAD__) && !defined(LEO2_TEST_THREAD_SANITIZED)
# define LEO2_TEST_THREAD_SANITIZED 1
#endif

#ifndef LEO2_TEST_THREAD_SANITIZED
# define LEO2_TEST_THREAD_SANITIZED 0
#endif

#if LEO2_TEST_THREAD_SANITIZED
# define LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE 0
#else
# define LEO2_TEST_ALLOCATION_AUDIT_AVAILABLE 1
#endif
