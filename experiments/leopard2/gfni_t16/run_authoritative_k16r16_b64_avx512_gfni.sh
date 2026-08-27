#!/usr/bin/bash
set -Eeuo pipefail
umask 077

if [[ ${1:-} != --leopard-t16-clean-env-internal ]]; then
    authoritative_script="$(/usr/bin/readlink -f "$0")"
    exec /usr/bin/env -i \
        PATH=/usr/bin:/bin \
        LANG=C \
        LC_ALL=C \
        TZ=UTC \
        TMPDIR=/tmp \
        OMP_DYNAMIC=FALSE \
        OMP_NUM_THREADS=1 \
        OMP_THREAD_LIMIT=1 \
        PYTHONNOUSERSITE=1 \
        PYTHONDONTWRITEBYTECODE=1 \
        CC=/usr/bin/cc \
        CXX=/usr/bin/c++ \
        CFLAGS= \
        CXXFLAGS= \
        CPPFLAGS= \
        LDFLAGS= \
        /usr/bin/bash "$authoritative_script" \
        --leopard-t16-clean-env-internal "$@"
fi
shift

for required_pair in \
    'PATH=/usr/bin:/bin' \
    'LANG=C' \
    'LC_ALL=C' \
    'TZ=UTC' \
    'TMPDIR=/tmp' \
    'OMP_DYNAMIC=FALSE' \
    'OMP_NUM_THREADS=1' \
    'OMP_THREAD_LIMIT=1' \
    'PYTHONNOUSERSITE=1' \
    'PYTHONDONTWRITEBYTECODE=1' \
    'CC=/usr/bin/cc' \
    'CXX=/usr/bin/c++' \
    'CFLAGS=' \
    'CXXFLAGS=' \
    'CPPFLAGS=' \
    'LDFLAGS='; do
    required_name=${required_pair%%=*}
    required_value=${required_pair#*=}
    if ! [[ -v $required_name ]]; then
        /usr/bin/printf 'missing authoritative environment key: %s\n' \
            "$required_name" >&2
        exit 2
    fi
    test "${!required_name}" = "$required_value"
done
while IFS= read -r -d '' environment_entry; do
    environment_name=${environment_entry%%=*}
    case "$environment_name" in
        PATH|LANG|LC_ALL|TZ|TMPDIR|OMP_DYNAMIC|OMP_NUM_THREADS|OMP_THREAD_LIMIT|PYTHONNOUSERSITE|PYTHONDONTWRITEBYTECODE|CC|CXX|CFLAGS|CXXFLAGS|CPPFLAGS|LDFLAGS|PWD|SHLVL|_)
            ;;
        *)
            /usr/bin/printf 'unexpected authoritative environment key: %s\n' \
                "$environment_name" >&2
            exit 2
            ;;
    esac
done < <(/usr/bin/env -0)

repo=/home/catid/leopard
relative_wrapper=experiments/leopard2/gfni_t16/run_authoritative_k16r16_b64_avx512_gfni.sh
relative_runner=experiments/leopard2/gfni_t16/run_k16r16_b64_avx512_gfni_abba.py
relative_auditor=experiments/leopard2/gfni_t16/audit_k16r16_b64_avx512_gfni_abba.py
relative_validator=tools/leopard2_benchmark_json_test.py
lock=/tmp/leopard-gf8-authoritative.lock

require_empty_output()
{
    local observed_output
    observed_output="$("$@")"
    test -z "$observed_output"
}

verify_envelope()
{
    verified_envelope=$1
    verified_core="$verified_envelope/core"
    test -d "$verified_envelope"
    test ! -L "$verified_envelope"
    test -d "$verified_core"
    test ! -L "$verified_core"
    test -f "$verified_envelope/COMPLETED.json"
    test -f "$verified_envelope/SHA256SUMS"
    test -f "$verified_core/SHA256SUMS"
    test -f "$verified_core/campaign.json"
    test -f "$verified_core/campaign.jsonl"
    test -f "$verified_core/audit.json"
    test -f "$verified_envelope/postseal-audit.json"
    require_empty_output /usr/bin/find "$verified_envelope" \
        -type l -print -quit
    require_empty_output /usr/bin/find "$verified_envelope" \
        -type f -links +1 -print -quit
    require_empty_output /usr/bin/find "$verified_envelope" \
        -type f -perm /222 -print -quit
    require_empty_output /usr/bin/find "$verified_envelope" \
        -type d -perm /222 -print -quit
    (
        cd "$verified_core"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null
    (
        cd "$verified_envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null
    /usr/bin/cmp "$verified_core/audit.json" \
        "$verified_envelope/postseal-audit.json"
    completed_status="$(/usr/bin/jq -er '.status | strings' \
        "$verified_envelope/COMPLETED.json")"
    completed_schema="$(/usr/bin/jq -er '.schema | strings' \
        "$verified_envelope/COMPLETED.json")"
    completed_promotion="$(/usr/bin/jq -er '.promotion_passed | booleans' \
        "$verified_envelope/COMPLETED.json")"
    completed_exit="$(/usr/bin/jq -er '.campaign_exit_status | numbers' \
        "$verified_envelope/COMPLETED.json")"
    completed_claim="$(/usr/bin/jq -er '.claim_passed | booleans' \
        "$verified_envelope/COMPLETED.json")"
    test "$completed_schema" = \
        leopard2-k16r16-b64-avx512-gfni-completion-envelope/v1
    test "$completed_status" = complete
    test "$completed_promotion" = true
    test "$completed_exit" = 0
    test "$completed_claim" = true
    test "$(/usr/bin/jq -r '.promotion_requires_completion_envelope' \
        "$verified_core/manifest.json")" = true
    test "$(/usr/bin/jq -r '.promotion_passed' \
        "$verified_core/manifest.json")" = false
    test "$(/usr/bin/sha256sum "$verified_core/manifest.json" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.core_manifest_sha256' \
            "$verified_envelope/COMPLETED.json")"
    test "$(/usr/bin/sha256sum "$verified_core/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.core_sha256sums_sha256' \
            "$verified_envelope/COMPLETED.json")"
    test "$(/usr/bin/sha256sum "$verified_core/audit.json" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.audit_sha256' \
            "$verified_envelope/COMPLETED.json")"
    test "$(/usr/bin/sha256sum "$verified_core/run-authoritative.sh" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.wrapper_sha256' \
            "$verified_core/manifest.json")"
    test "$(/usr/bin/sha256sum \
        "$verified_core/audit_k16r16_b64_avx512_gfni_abba.py" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.auditor_sha256' \
            "$verified_core/manifest.json")"
    /usr/bin/cmp "$verified_core/run-authoritative.sh" \
        "$verified_core/build-closure/committed-wrapper.sh"
    /usr/bin/cmp \
        "$verified_core/audit_k16r16_b64_avx512_gfni_abba.py" \
        "$verified_core/build-closure/committed-auditor.py"
    test "$(/usr/bin/sha256sum "$verified_envelope/postseal-audit.json" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -r '.postseal_audit_sha256' \
            "$verified_envelope/COMPLETED.json")"
    verifier_tmp="$(/usr/bin/mktemp -d /tmp/leopard-t16-envelope-verify.XXXXXX)"
    /usr/bin/python3 -I -S -B \
        "$verified_core/audit_k16r16_b64_avx512_gfni_abba.py" \
        --archive-only-source-closure \
        --report "$verified_core/campaign.json" \
        --journal "$verified_core/campaign.jsonl" \
        --output "$verifier_tmp/audit.json" \
        > "$verifier_tmp/audit-summary.json" \
        2> "$verifier_tmp/audit-stderr.log"
    /usr/bin/cmp "$verified_core/audit.json" "$verifier_tmp/audit.json"
    /usr/bin/printf 'authoritative envelope verified: %s\n' \
        "$verified_envelope"
}

if [[ $# -eq 2 && $1 == --verify && $2 == /* ]]; then
    verify_envelope "$2"
    exit 0
fi

if [[ $# -ne 1 || $1 != /* ]]; then
    /usr/bin/printf 'usage: %s /absolute/repository/.research/envelope\n' "$0" >&2
    /usr/bin/printf '       %s --verify /absolute/repository/.research/envelope\n' "$0" >&2
    exit 2
fi
envelope=$1
case "$envelope" in
    "$repo"/.research/*) ;;
    *)
        /usr/bin/printf 'artifact lane must be below %s/.research\n' "$repo" >&2
        exit 2
        ;;
esac

export GIT_CONFIG_GLOBAL=/dev/null
export GIT_CONFIG_SYSTEM=/dev/null
export GIT_CONFIG_NOSYSTEM=1
export GIT_NO_REPLACE_OBJECTS=1
export GIT_OPTIONAL_LOCKS=0
export LANG=C
export LC_ALL=C
export TZ=UTC
export OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

stage=initialization

next_stage()
{
    stage=$1
    /usr/bin/printf 'AUTHORITATIVE_STAGE %s\n' "$stage" | \
        /usr/bin/tee -a "$lane/wrapper-stage.log"
}

failure_record()
{
    status=$?
    set +e
    if [[ -d "$lane" && -w "$lane" ]]; then
        /usr/bin/jq -n \
            --arg stage "$stage" \
            --argjson exit_status "$status" \
            --arg commit "${commit:-unknown}" \
            --arg tree "${tree:-unknown}" \
            '{schema:"leopard2-t16-authoritative-preflight-failure/v1",stage:$stage,exit_status:$exit_status,source_commit:$commit,source_tree:$tree}' \
            > "$lane/preflight-failure.json" || true
        if [[ ! -e "$lane/SHA256SUMS" ]]; then
            (
                cd "$lane" || exit
                /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
                    /usr/bin/sort -z | \
                    /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
            ) || true
        fi
        /usr/bin/find "$lane" -type f -perm /222 \
            -exec /usr/bin/chmod a-w {} + || true
        /usr/bin/find "$lane" -type d -perm /222 \
            -exec /usr/bin/chmod a-w {} + || true
    fi
    /usr/bin/printf 'authoritative wrapper failed in stage %s with status %s\n' \
        "$stage" "$status" >&2
    exit "$status"
}

if [[ -e "$envelope" ]]; then
    /usr/bin/printf 'refusing to reuse artifact envelope: %s\n' "$envelope" >&2
    exit 2
fi
/usr/bin/mkdir -p "$(/usr/bin/dirname "$envelope")"
/usr/bin/mkdir -m 0700 "$envelope"
lane="$envelope/core"
/usr/bin/mkdir -m 0700 "$lane"
trap failure_record ERR

next_stage canonical_lock
exec 9>> "$lock"
/usr/bin/flock -n 9
test -f "$lock"
test ! -L "$lock"
test "$(/usr/bin/stat -c %h "$lock")" = 1

next_stage pre_tool_closure
/usr/bin/sha256sum \
    /usr/bin/bash /usr/bin/env /usr/bin/c++ /usr/bin/cc /usr/bin/ld \
    /usr/bin/ar /usr/bin/ranlib /usr/bin/cmake /usr/bin/ninja /usr/bin/git \
    /usr/bin/objdump /usr/bin/readelf /usr/bin/nm /usr/bin/python3 \
    /usr/bin/taskset /usr/bin/flock /usr/bin/jq /usr/bin/sha256sum \
    /usr/bin/find /usr/bin/grep /usr/bin/cmp /usr/bin/cp /usr/bin/stat \
    /usr/bin/chmod /usr/bin/sort /usr/bin/xargs /usr/bin/cut /usr/bin/tee \
    /usr/bin/ctest /usr/bin/cat /usr/bin/dirname /usr/bin/lscpu \
    /usr/bin/mkdir /usr/bin/mktemp /usr/bin/printf /usr/bin/readlink \
    /usr/bin/uname \
    > "$lane/pre-tool-SHA256SUMS"
/usr/bin/python3 -I -S -B -c \
    'import json,sys; assert sys.flags.isolated == 1 and sys.flags.no_site == 1 and sys.flags.ignore_environment == 1; print(json.dumps({"executable":sys.executable,"flags":{"ignore_environment":sys.flags.ignore_environment,"isolated":sys.flags.isolated,"no_site":sys.flags.no_site,"no_user_site":sys.flags.no_user_site}},sort_keys=True,separators=(",",":")))' \
    > "$lane/python-isolation.json"

next_stage clean_source_preflight
commit="$(/usr/bin/git -C "$repo" rev-parse HEAD)"
tree="$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')"
test "${#commit}" = 40
test "${#tree}" = 40
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/git -C "$repo" rev-parse --show-toplevel)" = "$repo"
test "$(/usr/bin/readlink -f "$repo/$relative_wrapper")" = \
    "$(/usr/bin/readlink -f "$0")"
/usr/bin/git -C "$repo" ls-tree HEAD sse2neon > "$lane/live-sse2neon-gitlink.txt"
/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal \
    > "$lane/live-git-status.txt"
/usr/bin/printf '%s\n' "$commit" > "$lane/source-commit.txt"
/usr/bin/printf '%s\n' "$tree" > "$lane/source-tree.txt"
/usr/bin/printf '%s\n' "$(/usr/bin/cat /sys/devices/system/cpu/cpu10/topology/thread_siblings_list)" \
    > "$lane/cpu10-thread-siblings.txt"
test "$(/usr/bin/cat "$lane/cpu10-thread-siblings.txt")" = 10,74

next_stage fresh_build_roots
live_root="$(/usr/bin/mktemp -d /tmp/leopard-t16-live.XXXXXX)"
live_build="$live_root/build"
replay_root="$(/usr/bin/mktemp -d /tmp/leopard-t16-replay.XXXXXX)"
replay_source="$replay_root/source"
replay_build="$replay_root/build"
/usr/bin/jq -n \
    --arg live_root "$live_root" \
    --arg live_build "$live_build" \
    --arg replay_root "$replay_root" \
    --arg replay_source "$replay_source" \
    --arg replay_build "$replay_build" \
    '{live_root:$live_root,live_build:$live_build,replay_root:$replay_root,replay_source:$replay_source,replay_build:$replay_build}' \
    > "$lane/build-paths.json"

next_stage replay_clone
/usr/bin/git clone --no-hardlinks --no-checkout "$repo" "$replay_source" \
    > "$lane/replay-clone.log" 2>&1
/usr/bin/git -C "$replay_source" checkout --detach "$commit" \
    > "$lane/replay-checkout.log" 2>&1
test "$(/usr/bin/git -C "$replay_source" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$replay_source" rev-parse 'HEAD^{tree}')" = "$tree"
require_empty_output /usr/bin/git -C "$replay_source" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$replay_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
/usr/bin/git -C "$replay_source" ls-tree HEAD sse2neon > "$lane/replay-sse2neon-gitlink.txt"
/usr/bin/cmp "$lane/live-sse2neon-gitlink.txt" "$lane/replay-sse2neon-gitlink.txt"

configure_common=(
    -G Ninja
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    -DCMAKE_C_COMPILER=/usr/bin/cc
    -DCMAKE_CXX_COMPILER=/usr/bin/c++
    -DCMAKE_AR=/usr/bin/ar
    -DCMAKE_RANLIB=/usr/bin/ranlib
    -DPython3_EXECUTABLE=/usr/bin/python3
    -DCMAKE_C_FLAGS=
    -DCMAKE_CXX_FLAGS=
    -DCMAKE_EXE_LINKER_FLAGS=
    -DCMAKE_SHARED_LINKER_FLAGS=
    -DLEO2_BUILD_TESTS=ON
    -DLEO2_BUILD_BENCHMARKS=ON
    -DLEO2_ENABLE_CUDA=OFF
    -DLEO2_BACKEND_VARIANT=auto
    -DLEOPARD_ENABLE_GF8=ON
    -DLEOPARD_ENABLE_GF16=ON
    -DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=ON
)
build_targets=(
    bench_leopard2
    bench_leopard2_prevalidated_batch
    leopard2_balanced_b64_terminal_test
    leopard2_balanced_b64_terminal_production_test
    leopard2_backend_ops_test
    leopard2_avx512_gfni_t16_prototype_test
)
focused_regex='^leopard2_(portable_isa|portable_isa_registration|balanced_b64_terminal|balanced_b64_terminal_production|backend_ops|avx512_gfni_t16_prototype|benchmark_json_regression)$'

assert_pinned_cache()
{
    cache=$1/CMakeCache.txt
    /usr/bin/grep -Ex \
        'CMAKE_C_COMPILER:(FILEPATH|STRING)=/usr/bin/cc' "$cache"
    /usr/bin/grep -Ex \
        'CMAKE_CXX_COMPILER:(FILEPATH|STRING)=/usr/bin/c\+\+' "$cache"
    /usr/bin/grep -Fx 'CMAKE_AR:FILEPATH=/usr/bin/ar' "$cache"
    /usr/bin/grep -Fx 'CMAKE_RANLIB:FILEPATH=/usr/bin/ranlib' "$cache"
    /usr/bin/grep -Ex \
        'Python3_EXECUTABLE:(UNINITIALIZED|FILEPATH)=/usr/bin/python3' "$cache"
    /usr/bin/grep -Fx 'CMAKE_C_FLAGS:STRING=' "$cache"
    /usr/bin/grep -Fx 'CMAKE_CXX_FLAGS:STRING=' "$cache"
    /usr/bin/grep -Fx 'CMAKE_EXE_LINKER_FLAGS:STRING=' "$cache"
    /usr/bin/grep -Fx 'CMAKE_SHARED_LINKER_FLAGS:STRING=' "$cache"
}

assert_effective_compiler_configuration()
{
    configuration=$1/generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt
    /usr/bin/grep -Fx 'CMAKE_CXX_COMPILER=/usr/bin/c++' "$configuration"
}

assert_focused_test_census()
{
    build=$1
    output=$2
    /usr/bin/ctest --test-dir "$build" -N -R "$focused_regex" \
        --show-only=json-v1 > "$output"
    /usr/bin/jq -e '
        [.tests[].name] == [
            "leopard2_portable_isa",
            "leopard2_portable_isa_registration",
            "leopard2_balanced_b64_terminal",
            "leopard2_balanced_b64_terminal_production",
            "leopard2_backend_ops",
            "leopard2_avx512_gfni_t16_prototype",
            "leopard2_benchmark_json_regression"
        ]
    ' "$output" >/dev/null
}

next_stage live_configure
/usr/bin/cmake -S "$repo" -B "$live_build" "${configure_common[@]}" \
    > "$lane/live-configure.log" 2>&1
assert_pinned_cache "$live_build" >> "$lane/live-configure.log"
assert_effective_compiler_configuration "$live_build" \
    >> "$lane/live-configure.log"

next_stage live_build
/usr/bin/cmake --build "$live_build" --parallel 2 --target "${build_targets[@]}" \
    > "$lane/live-build.log" 2>&1

next_stage live_focused_test_census
assert_focused_test_census "$live_build" \
    "$lane/live-focused-test-census.json"

next_stage live_focused_tests
/usr/bin/ctest --test-dir "$live_build" -j1 --output-on-failure \
    --no-tests=error -R "$focused_regex" \
    > "$lane/live-focused-tests.log" 2>&1

next_stage replay_configure
/usr/bin/cmake -S "$replay_source" -B "$replay_build" "${configure_common[@]}" \
    > "$lane/replay-configure.log" 2>&1
assert_pinned_cache "$replay_build" >> "$lane/replay-configure.log"
assert_effective_compiler_configuration "$replay_build" \
    >> "$lane/replay-configure.log"

next_stage replay_build
/usr/bin/cmake --build "$replay_build" --parallel 2 --target "${build_targets[@]}" \
    > "$lane/replay-build.log" 2>&1

next_stage replay_focused_test_census
assert_focused_test_census "$replay_build" \
    "$lane/replay-focused-test-census.json"

next_stage replay_focused_tests
/usr/bin/ctest --test-dir "$replay_build" -j1 --output-on-failure \
    --no-tests=error -R "$focused_regex" \
    > "$lane/replay-focused-tests.log" 2>&1

live_t16_object="$live_build/CMakeFiles/leopard2_backend_avx512_gfni_t16.dir/Leopard2BackendAVX512GFNIT16.cpp.o"
replay_t16_object="$replay_build/CMakeFiles/leopard2_backend_avx512_gfni_t16.dir/Leopard2BackendAVX512GFNIT16.cpp.o"
live_benchmark_object="$live_build/CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o"
replay_benchmark_object="$replay_build/CMakeFiles/bench_leopard2.dir/bench/leopard2/benchmark.cpp.o"
live_attestation="$live_build/generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h"
replay_attestation="$replay_build/generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h"
live_configuration="$live_build/generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt"
replay_configuration="$replay_build/generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt"

next_stage byte_replay
/usr/bin/cmp "$live_build/bench_leopard2" "$replay_build/bench_leopard2"
/usr/bin/cmp "$live_build/bench_leopard2_prevalidated_batch" \
    "$replay_build/bench_leopard2_prevalidated_batch"
/usr/bin/cmp "$live_build/libleopard.a" "$replay_build/libleopard.a"
/usr/bin/cmp "$live_t16_object" "$replay_t16_object"
/usr/bin/cmp "$live_benchmark_object" "$replay_benchmark_object"
/usr/bin/cmp "$live_attestation" "$replay_attestation"
/usr/bin/cmp "$live_configuration" "$replay_configuration"

next_stage canonical_source_archive
/usr/bin/git -C "$repo" archive --format=tar --prefix=source/ \
    --output="$lane/source.tar" "$commit"
/usr/bin/git -C "$replay_source" archive --format=tar --prefix=source/ \
    --output="$replay_root/replay-source.tar" "$commit"
/usr/bin/cmp "$lane/source.tar" "$replay_root/replay-source.tar"

next_stage freeze_core
/usr/bin/cp --reflink=never "$replay_build/bench_leopard2" "$lane/bench_leopard2"
/usr/bin/cp --reflink=never "$repo/$relative_runner" \
    "$lane/run_k16r16_b64_avx512_gfni_abba.py"
/usr/bin/cp --reflink=never "$repo/$relative_auditor" \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py"
/usr/bin/cp --reflink=never "$repo/$relative_validator" \
    "$lane/leopard2_benchmark_json_test.py"
/usr/bin/cp --reflink=never "$repo/$relative_wrapper" \
    "$lane/run-authoritative.sh"
/usr/bin/chmod 0555 \
    "$lane/bench_leopard2" \
    "$lane/run_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/run-authoritative.sh"
/usr/bin/chmod 0444 "$lane/leopard2_benchmark_json_test.py" "$lane/source.tar"
for frozen in \
    "$lane/bench_leopard2" \
    "$lane/run_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/leopard2_benchmark_json_test.py" \
    "$lane/run-authoritative.sh" \
    "$lane/source.tar"; do
    test ! -L "$frozen"
    test "$(/usr/bin/stat -c %h "$frozen")" = 1
done
test "$(/usr/bin/stat -c %i "$lane/bench_leopard2")" != \
    "$(/usr/bin/stat -c %i "$replay_build/bench_leopard2")"

next_stage retain_build_closure
/usr/bin/mkdir -m 0700 "$lane/build-closure"
/usr/bin/mkdir -m 0700 "$lane/build-closure/live" "$lane/build-closure/replay"

copy_build_closure()
{
    build=$1
    destination=$2
    t16_object=$3
    benchmark_object=$4
    attestation=$5
    configuration=$6
    /usr/bin/cp --reflink=never "$build/CMakeCache.txt" "$destination/CMakeCache.txt"
    /usr/bin/cp --reflink=never "$build/compile_commands.json" "$destination/compile_commands.json"
    /usr/bin/cp --reflink=never "$build/build.ninja" "$destination/build.ninja"
    /usr/bin/cp --reflink=never "$build/libleopard.a" "$destination/libleopard.a"
    /usr/bin/cp --reflink=never "$build/bench_leopard2" "$destination/bench_leopard2"
    /usr/bin/cp --reflink=never "$build/bench_leopard2_prevalidated_batch" \
        "$destination/bench_leopard2_prevalidated_batch"
    /usr/bin/cp --reflink=never "$t16_object" \
        "$destination/Leopard2BackendAVX512GFNIT16.cpp.o"
    /usr/bin/cp --reflink=never "$benchmark_object" "$destination/benchmark.cpp.o"
    /usr/bin/cp --reflink=never "$attestation" \
        "$destination/leopard2_benchmark_source_attestation.h"
    /usr/bin/cp --reflink=never "$configuration" \
        "$destination/leopard2_benchmark_build_configuration.txt"
}

copy_build_closure "$live_build" "$lane/build-closure/live" \
    "$live_t16_object" "$live_benchmark_object" \
    "$live_attestation" "$live_configuration"
copy_build_closure "$replay_build" "$lane/build-closure/replay" \
    "$replay_t16_object" "$replay_benchmark_object" \
    "$replay_attestation" "$replay_configuration"

validate_retained_build_closure()
{
    /usr/bin/cmp "$lane/build-closure/live/CMakeCache.txt" \
        "$live_build/CMakeCache.txt"
    /usr/bin/cmp "$lane/build-closure/live/compile_commands.json" \
        "$live_build/compile_commands.json"
    /usr/bin/cmp "$lane/build-closure/live/build.ninja" \
        "$live_build/build.ninja"
    /usr/bin/cmp "$lane/build-closure/live/libleopard.a" \
        "$live_build/libleopard.a"
    /usr/bin/cmp "$lane/build-closure/live/bench_leopard2" \
        "$live_build/bench_leopard2"
    /usr/bin/cmp "$lane/build-closure/live/bench_leopard2_prevalidated_batch" \
        "$live_build/bench_leopard2_prevalidated_batch"
    /usr/bin/cmp "$lane/build-closure/live/Leopard2BackendAVX512GFNIT16.cpp.o" \
        "$live_t16_object"
    /usr/bin/cmp "$lane/build-closure/live/benchmark.cpp.o" \
        "$live_benchmark_object"
    /usr/bin/cmp "$lane/build-closure/live/leopard2_benchmark_source_attestation.h" \
        "$live_attestation"
    /usr/bin/cmp "$lane/build-closure/live/leopard2_benchmark_build_configuration.txt" \
        "$live_configuration"
    /usr/bin/cmp "$lane/build-closure/replay/CMakeCache.txt" \
        "$replay_build/CMakeCache.txt"
    /usr/bin/cmp "$lane/build-closure/replay/compile_commands.json" \
        "$replay_build/compile_commands.json"
    /usr/bin/cmp "$lane/build-closure/replay/build.ninja" \
        "$replay_build/build.ninja"
    /usr/bin/cmp "$lane/build-closure/replay/libleopard.a" \
        "$replay_build/libleopard.a"
    /usr/bin/cmp "$lane/build-closure/replay/bench_leopard2" \
        "$replay_build/bench_leopard2"
    /usr/bin/cmp "$lane/build-closure/replay/bench_leopard2_prevalidated_batch" \
        "$replay_build/bench_leopard2_prevalidated_batch"
    /usr/bin/cmp "$lane/build-closure/replay/Leopard2BackendAVX512GFNIT16.cpp.o" \
        "$replay_t16_object"
    /usr/bin/cmp "$lane/build-closure/replay/benchmark.cpp.o" \
        "$replay_benchmark_object"
    /usr/bin/cmp "$lane/build-closure/replay/leopard2_benchmark_source_attestation.h" \
        "$replay_attestation"
    /usr/bin/cmp "$lane/build-closure/replay/leopard2_benchmark_build_configuration.txt" \
        "$replay_configuration"
    /usr/bin/cmp "$lane/build-closure/live/libleopard.a" \
        "$lane/build-closure/replay/libleopard.a"
    /usr/bin/cmp "$lane/build-closure/live/bench_leopard2" \
        "$lane/build-closure/replay/bench_leopard2"
    /usr/bin/cmp "$lane/build-closure/live/bench_leopard2_prevalidated_batch" \
        "$lane/build-closure/replay/bench_leopard2_prevalidated_batch"
    /usr/bin/cmp "$lane/build-closure/live/Leopard2BackendAVX512GFNIT16.cpp.o" \
        "$lane/build-closure/replay/Leopard2BackendAVX512GFNIT16.cpp.o"
    /usr/bin/cmp "$lane/build-closure/live/benchmark.cpp.o" \
        "$lane/build-closure/replay/benchmark.cpp.o"
    /usr/bin/cmp "$lane/build-closure/live/leopard2_benchmark_source_attestation.h" \
        "$lane/build-closure/replay/leopard2_benchmark_source_attestation.h"
    /usr/bin/cmp "$lane/build-closure/live/leopard2_benchmark_build_configuration.txt" \
        "$lane/build-closure/replay/leopard2_benchmark_build_configuration.txt"
}

validate_retained_build_closure

/usr/bin/ninja -C "$live_build" -t commands bench_leopard2 \
    > "$lane/build-closure/live-ninja-commands.txt"
/usr/bin/ninja -C "$replay_build" -t commands bench_leopard2 \
    > "$lane/build-closure/replay-ninja-commands.txt"
/usr/bin/objdump -drwC -Mintel "$replay_t16_object" \
    > "$lane/build-closure/t16-object-disassembly.txt"
/usr/bin/objdump -drwC -Mintel -j .text.leo2_avx512_gfni_t16 \
    "$replay_t16_object" > "$lane/build-closure/t16-hot-section-disassembly.txt"
/usr/bin/objdump -h "$replay_t16_object" \
    > "$lane/build-closure/t16-object-sections.txt"
/usr/bin/readelf -SW "$replay_t16_object" \
    > "$lane/build-closure/t16-object-readelf-sections.txt"
/usr/bin/readelf -Ws "$replay_t16_object" \
    > "$lane/build-closure/t16-object-readelf-symbols.txt"
/usr/bin/nm -S --size-sort "$replay_t16_object" \
    > "$lane/build-closure/t16-object-nm-sizes.txt"
/usr/bin/sha256sum \
    /usr/bin/bash /usr/bin/env /usr/bin/c++ /usr/bin/cc /usr/bin/ld \
    /usr/bin/ar /usr/bin/ranlib /usr/bin/cmake /usr/bin/ninja /usr/bin/git \
    /usr/bin/objdump /usr/bin/readelf /usr/bin/nm /usr/bin/python3 \
    /usr/bin/taskset /usr/bin/flock /usr/bin/jq /usr/bin/sha256sum \
    /usr/bin/find /usr/bin/grep /usr/bin/cmp /usr/bin/cp /usr/bin/stat \
    /usr/bin/chmod /usr/bin/sort /usr/bin/xargs /usr/bin/cut /usr/bin/tee \
    /usr/bin/ctest /usr/bin/cat /usr/bin/dirname /usr/bin/lscpu \
    /usr/bin/mkdir /usr/bin/mktemp /usr/bin/printf /usr/bin/readlink \
    /usr/bin/uname \
    > "$lane/build-closure/tool-SHA256SUMS"
/usr/bin/c++ --version > "$lane/build-closure/compiler-version.txt"
/usr/bin/cc --version > "$lane/build-closure/c-compiler-version.txt"
/usr/bin/cmake --version > "$lane/build-closure/cmake-version.txt"
/usr/bin/ninja --version > "$lane/build-closure/ninja-version.txt"
/usr/bin/git --version > "$lane/build-closure/git-version.txt"
/usr/bin/python3 --version > "$lane/build-closure/python-version.txt" 2>&1
/usr/bin/lscpu --json > "$lane/build-closure/lscpu.json"
/usr/bin/uname -a > "$lane/build-closure/uname.txt"

/usr/bin/git -C "$repo" show "$commit:$relative_runner" \
    > "$lane/build-closure/committed-runner.py"
/usr/bin/git -C "$repo" show "$commit:$relative_auditor" \
    > "$lane/build-closure/committed-auditor.py"
/usr/bin/git -C "$repo" show "$commit:$relative_validator" \
    > "$lane/build-closure/committed-validator.py"
/usr/bin/git -C "$repo" show "$commit:$relative_wrapper" \
    > "$lane/build-closure/committed-wrapper.sh"
/usr/bin/cmp "$lane/run_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/build-closure/committed-runner.py"
/usr/bin/cmp "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/build-closure/committed-auditor.py"
/usr/bin/cmp "$lane/leopard2_benchmark_json_test.py" \
    "$lane/build-closure/committed-validator.py"
/usr/bin/cmp "$lane/run-authoritative.sh" \
    "$lane/build-closure/committed-wrapper.sh"

binary_hash="$(/usr/bin/sha256sum "$lane/bench_leopard2" | /usr/bin/cut -d' ' -f1)"
runner_hash="$(/usr/bin/sha256sum "$lane/run_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)"
auditor_hash="$(/usr/bin/sha256sum "$lane/audit_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)"
validator_hash="$(/usr/bin/sha256sum "$lane/leopard2_benchmark_json_test.py" | /usr/bin/cut -d' ' -f1)"
wrapper_hash="$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | /usr/bin/cut -d' ' -f1)"
archive_hash="$(/usr/bin/sha256sum "$lane/source.tar" | /usr/bin/cut -d' ' -f1)"
library_hash="$(/usr/bin/sha256sum "$lane/build-closure/replay/libleopard.a" | /usr/bin/cut -d' ' -f1)"
t16_object_hash="$(/usr/bin/sha256sum "$lane/build-closure/replay/Leopard2BackendAVX512GFNIT16.cpp.o" | /usr/bin/cut -d' ' -f1)"
benchmark_object_hash="$(/usr/bin/sha256sum "$lane/build-closure/replay/benchmark.cpp.o" | /usr/bin/cut -d' ' -f1)"
attestation_hash="$(/usr/bin/sha256sum "$lane/build-closure/replay/leopard2_benchmark_source_attestation.h" | /usr/bin/cut -d' ' -f1)"
configuration_hash="$(/usr/bin/sha256sum "$lane/build-closure/replay/leopard2_benchmark_build_configuration.txt" | /usr/bin/cut -d' ' -f1)"

/usr/bin/jq -n \
    --arg repository "$repo" \
    --arg live_root "$live_root" \
    --arg live_build "$live_build" \
    --arg replay_root "$replay_root" \
    --arg replay_source "$replay_source" \
    --arg replay_build "$replay_build" \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg binary_sha256 "$binary_hash" \
    --arg library_sha256 "$library_hash" \
    --arg t16_object_sha256 "$t16_object_hash" \
    --arg benchmark_object_sha256 "$benchmark_object_hash" \
    --arg attestation_sha256 "$attestation_hash" \
    --arg configuration_sha256 "$configuration_hash" \
    --arg runner_sha256 "$runner_hash" \
    --arg auditor_sha256 "$auditor_hash" \
    --arg validator_sha256 "$validator_hash" \
    --arg wrapper_sha256 "$wrapper_hash" \
    --arg archive_sha256 "$archive_hash" \
    '{schema:"leopard2-k16r16-b64-avx512-gfni-build-closure/v1",repository:$repository,live_root:$live_root,live_build:$live_build,replay_root:$replay_root,replay_source:$replay_source,replay_build:$replay_build,source_commit:$commit,source_tree:$tree,byte_identical:{live_replay_binary:true,live_replay_prevalidated_binary:true,live_replay_library:true,live_replay_t16_object:true,live_replay_benchmark_object:true,live_replay_attestation:true,live_replay_configuration:true,live_replay_archive:true},frozen:{binary_sha256:$binary_sha256,library_sha256:$library_sha256,t16_object_sha256:$t16_object_sha256,benchmark_object_sha256:$benchmark_object_sha256,attestation_sha256:$attestation_sha256,configuration_sha256:$configuration_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256},configure:["/usr/bin/cmake","-G","Ninja","-DCMAKE_BUILD_TYPE=Release","-DCMAKE_EXPORT_COMPILE_COMMANDS=ON","-DCMAKE_C_COMPILER=/usr/bin/cc","-DCMAKE_CXX_COMPILER=/usr/bin/c++","-DCMAKE_AR=/usr/bin/ar","-DCMAKE_RANLIB=/usr/bin/ranlib","-DPython3_EXECUTABLE=/usr/bin/python3","-DCMAKE_C_FLAGS=","-DCMAKE_CXX_FLAGS=","-DCMAKE_EXE_LINKER_FLAGS=","-DCMAKE_SHARED_LINKER_FLAGS=","-DLEO2_BUILD_TESTS=ON","-DLEO2_BUILD_BENCHMARKS=ON","-DLEO2_ENABLE_CUDA=OFF","-DLEO2_BACKEND_VARIANT=auto","-DLEOPARD_ENABLE_GF8=ON","-DLEOPARD_ENABLE_GF16=ON","-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=ON"],build:["/usr/bin/cmake","--build","<build>","--parallel","2","--target","bench_leopard2","bench_leopard2_prevalidated_batch","leopard2_balanced_b64_terminal_test","leopard2_balanced_b64_terminal_production_test","leopard2_backend_ops_test","leopard2_avx512_gfni_t16_prototype_test"],focused_test_regex:"^leopard2_(portable_isa|portable_isa_registration|balanced_b64_terminal|balanced_b64_terminal_production|backend_ops|avx512_gfni_t16_prototype|benchmark_json_regression)$",focused_tests:["leopard2_portable_isa","leopard2_portable_isa_registration","leopard2_balanced_b64_terminal","leopard2_balanced_b64_terminal_production","leopard2_backend_ops","leopard2_avx512_gfni_t16_prototype","leopard2_benchmark_json_regression"],python_controller:["/usr/bin/python3","-I","-S","-B"],environment_policy:"env -i allowlist with pinned compiler paths, empty compile/link flags, no Python user site, and pre/post tool hashes"}' \
    > "$lane/build-closure.json"
/usr/bin/cmp "$lane/pre-tool-SHA256SUMS" \
    "$lane/build-closure/tool-SHA256SUMS"
/usr/bin/find "$lane/build-closure" -type f -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$lane/build-closure" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/chmod a-w "$lane/build-closure.json"

next_stage campaign
campaign_command=(
    /usr/bin/taskset -c 10
    /usr/bin/python3 -I -S -B
    "$lane/run_k16r16_b64_avx512_gfni_abba.py"
    --binary "$lane/bench_leopard2"
    --build-binary "$replay_build/bench_leopard2"
    --binary-sha256 "$binary_hash"
    --runner-sha256 "$runner_hash"
    --validator "$lane/leopard2_benchmark_json_test.py"
    --validator-sha256 "$validator_hash"
    --source-archive "$lane/source.tar"
    --source-archive-sha256 "$archive_hash"
    --source-commit "$commit"
    --source-tree "$tree"
    --repository "$repo"
    --cpu 10
    --sibling 74
    --output "$lane/campaign.json"
    --journal "$lane/campaign.jsonl"
    --invocations "$lane/invocations"
    --lock-fd 9
)
/usr/bin/jq -n --args '$ARGS.positional' -- "${campaign_command[@]}" \
    > "$lane/campaign-command.json"
trap - ERR
set +e
"${campaign_command[@]}" > "$lane/campaign-summary.json" \
    2> "$lane/campaign-stderr.log"
campaign_status=$?
set -e
trap failure_record ERR

next_stage post_campaign_closure
/usr/bin/cmp "$live_build/bench_leopard2" "$replay_build/bench_leopard2"
/usr/bin/cmp "$lane/bench_leopard2" "$replay_build/bench_leopard2"
/usr/bin/cmp "$live_build/libleopard.a" "$replay_build/libleopard.a"
/usr/bin/cmp "$live_t16_object" "$replay_t16_object"
/usr/bin/cmp "$live_benchmark_object" "$replay_benchmark_object"
/usr/bin/cmp "$live_attestation" "$replay_attestation"
/usr/bin/cmp "$live_configuration" "$replay_configuration"
validate_retained_build_closure
require_empty_output /usr/bin/find "$lane/build-closure" \
    -type f -perm /222 -print -quit
require_empty_output /usr/bin/find "$lane/build-closure" \
    -type d -perm /222 -print -quit
test ! -w "$lane/build-closure.json"
/usr/bin/sha256sum \
    /usr/bin/bash /usr/bin/env /usr/bin/c++ /usr/bin/cc /usr/bin/ld \
    /usr/bin/ar /usr/bin/ranlib /usr/bin/cmake /usr/bin/ninja /usr/bin/git \
    /usr/bin/objdump /usr/bin/readelf /usr/bin/nm /usr/bin/python3 \
    /usr/bin/taskset /usr/bin/flock /usr/bin/jq /usr/bin/sha256sum \
    /usr/bin/find /usr/bin/grep /usr/bin/cmp /usr/bin/cp /usr/bin/stat \
    /usr/bin/chmod /usr/bin/sort /usr/bin/xargs /usr/bin/cut /usr/bin/tee \
    /usr/bin/ctest /usr/bin/cat /usr/bin/dirname /usr/bin/lscpu \
    /usr/bin/mkdir /usr/bin/mktemp /usr/bin/printf /usr/bin/readlink \
    /usr/bin/uname \
    > "$lane/post-tool-SHA256SUMS"
/usr/bin/cmp "$lane/pre-tool-SHA256SUMS" "$lane/post-tool-SHA256SUMS"
test "$(/usr/bin/sha256sum "$lane/bench_leopard2" | /usr/bin/cut -d' ' -f1)" = "$binary_hash"
test "$(/usr/bin/sha256sum "$lane/run_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)" = "$runner_hash"
test "$(/usr/bin/sha256sum "$lane/audit_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)" = "$auditor_hash"
test "$(/usr/bin/sha256sum "$lane/leopard2_benchmark_json_test.py" | /usr/bin/cut -d' ' -f1)" = "$validator_hash"
test "$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | /usr/bin/cut -d' ' -f1)" = "$wrapper_hash"
/usr/bin/cmp "$lane/run-authoritative.sh" \
    "$lane/build-closure/committed-wrapper.sh"
test "$(/usr/bin/sha256sum "$lane/source.tar" | /usr/bin/cut -d' ' -f1)" = "$archive_hash"
test "$(/usr/bin/git -C "$repo" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" = "$tree"
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
/usr/bin/git -C "$repo" archive --format=tar --prefix=source/ \
    --output="$replay_root/post-source.tar" "$commit"
/usr/bin/cmp "$lane/source.tar" "$replay_root/post-source.tar"
test -f "$lane/campaign.json"
report_status="$(/usr/bin/jq -er '.status | strings' "$lane/campaign.json")"
if [[ "$report_status" == failed ]]; then
    next_stage retain_failed_campaign
    failed_report_sha="$(/usr/bin/sha256sum "$lane/campaign.json" | /usr/bin/cut -d' ' -f1)"
    failed_journal_sha=null
    if [[ -f "$lane/campaign.jsonl" ]]; then
        failed_journal_sha="$(/usr/bin/sha256sum "$lane/campaign.jsonl" | /usr/bin/cut -d' ' -f1)"
    fi
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg binary_sha256 "$binary_hash" \
        --arg runner_sha256 "$runner_hash" \
        --arg auditor_sha256 "$auditor_hash" \
        --arg validator_sha256 "$validator_hash" \
        --arg wrapper_sha256 "$wrapper_hash" \
        --arg archive_sha256 "$archive_hash" \
        --arg report_sha256 "$failed_report_sha" \
        --arg journal_sha256 "$failed_journal_sha" \
        '{schema:"leopard2-k16r16-b64-avx512-gfni-failed-manifest/v1",campaign_exit_status:$campaign_exit_status,report_status:"failed",claim_passed:false,audit_performed:false,source_commit:$commit,source_tree:$tree,binary_sha256_pre_post:$binary_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256,campaign_sha256:$report_sha256,journal_sha256:($journal_sha256 | if . == "null" then null else . end),canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:10,sibling:74}' \
        > "$lane/manifest.json"
    next_stage seal_failed_campaign
    trap - ERR
    (
        cd "$lane"
        /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
            /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum \
            > SHA256SUMS
    )
    /usr/bin/find "$lane" -type f -perm /222 -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$lane" -type d -perm /222 -exec /usr/bin/chmod a-w {} +
    failed_verification="$(/usr/bin/mktemp -d /tmp/leopard-t16-failed-seal.XXXXXX)"
    (
        cd "$lane"
        /usr/bin/sha256sum -c SHA256SUMS
    ) > "$failed_verification/sha256-check.txt"
    require_empty_output /usr/bin/find "$lane" -type l -print -quit
    require_empty_output /usr/bin/find "$lane" \
        -type f -links +1 -print -quit
    require_empty_output /usr/bin/find "$lane" \
        -type f -perm /222 -print -quit
    require_empty_output /usr/bin/find "$lane" \
        -type d -perm /222 -print -quit
    failed_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg core_sha256sums_sha256 "$failed_core_sha" \
        '{schema:"leopard2-k16r16-b64-avx512-gfni-failed-envelope/v1",status:"failed",promotion_passed:false,campaign_exit_status:$campaign_exit_status,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256}' \
        > "$envelope/FAILED.json"
    (
        cd "$envelope"
        /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
            /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum \
            > SHA256SUMS
    )
    /usr/bin/find "$envelope" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$envelope" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    (
        cd "$envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) > "$failed_verification/envelope-sha256-check.txt"
    /usr/bin/printf 'sealed_failed_envelope=%s\nverification=%s\n' \
        "$envelope" "$failed_verification"
    if [[ "$campaign_status" -eq 0 ]]; then
        exit 1
    fi
    exit "$campaign_status"
fi
test "$report_status" = complete
reported_claim="$(/usr/bin/jq -er \
    '.claim_passed as $value | if ($value | type) == "boolean" then ($value | tostring) else error("claim_passed is not boolean") end' \
    "$lane/campaign.json")"
if [[ "$reported_claim" == true ]]; then
    test "$campaign_status" -eq 0
else
    test "$campaign_status" -ne 0
fi

next_stage independent_preseal_audit
/usr/bin/python3 -I -S -B \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    --archive-only-source-closure \
    --report "$lane/campaign.json" \
    --journal "$lane/campaign.jsonl" \
    --output "$lane/audit.json" \
    > "$lane/audit-summary.json" 2> "$lane/audit-stderr.log"
test "$(/usr/bin/jq -r '.audit_passed' "$lane/audit.json")" = true
test "$(/usr/bin/jq -r '.claim_passed' "$lane/audit.json")" = \
    "$(/usr/bin/jq -r '.claim_passed' "$lane/campaign.json")"

next_stage final_manifest
claim_passed="$(/usr/bin/jq -r '.claim_passed' "$lane/campaign.json")"
campaign_sha="$(/usr/bin/sha256sum "$lane/campaign.json" | /usr/bin/cut -d' ' -f1)"
journal_sha="$(/usr/bin/sha256sum "$lane/campaign.jsonl" | /usr/bin/cut -d' ' -f1)"
audit_sha="$(/usr/bin/sha256sum "$lane/audit.json" | /usr/bin/cut -d' ' -f1)"
/usr/bin/jq -n \
    --argjson campaign_exit_status "$campaign_status" \
    --argjson claim_passed "$claim_passed" \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg binary_sha256 "$binary_hash" \
    --arg runner_sha256 "$runner_hash" \
    --arg auditor_sha256 "$auditor_hash" \
    --arg validator_sha256 "$validator_hash" \
    --arg wrapper_sha256 "$wrapper_hash" \
    --arg archive_sha256 "$archive_hash" \
    --arg campaign_sha256 "$campaign_sha" \
    --arg journal_sha256 "$journal_sha" \
    --arg audit_sha256 "$audit_sha" \
    --arg replay_root "$replay_root" \
    '{schema:"leopard2-k16r16-b64-avx512-gfni-core-manifest/v1",campaign_exit_status:$campaign_exit_status,claim_passed:$claim_passed,preseal_audit_passed:true,promotion_passed:false,promotion_requires_completion_envelope:true,source_commit:$commit,source_tree:$tree,binary_sha256_pre_post:$binary_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256,campaign_sha256:$campaign_sha256,journal_sha256:$journal_sha256,audit_sha256:$audit_sha256,replay_root:$replay_root,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:10,sibling:74,postseal_policy:"promotion requires the enclosing COMPLETED.json written only after a byte-identical postseal independent audit, core SHA verification, clean-source recheck, and zero controller exit"}' \
    > "$lane/manifest.json"

next_stage seal_core
trap - ERR
(
    cd "$lane"
    /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
        /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
)
/usr/bin/find "$lane" -type f -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$lane" -type d -perm /222 -exec /usr/bin/chmod a-w {} +

/usr/bin/printf 'AUTHORITATIVE_STAGE independent_postseal_audit\n'
/usr/bin/python3 -I -S -B \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    --archive-only-source-closure \
    --report "$lane/campaign.json" \
    --journal "$lane/campaign.jsonl" \
    --output "$envelope/postseal-audit.json" \
    > "$envelope/postseal-audit-summary.json" \
    2> "$envelope/postseal-audit-stderr.log"
/usr/bin/cmp "$lane/audit.json" "$envelope/postseal-audit.json"
(
    cd "$lane"
    /usr/bin/sha256sum -c SHA256SUMS
) > "$envelope/core-sha256-check.txt"
require_empty_output /usr/bin/find "$lane" -type l -print -quit
require_empty_output /usr/bin/find "$lane" \
    -type f -links +1 -print -quit
require_empty_output /usr/bin/find "$lane" \
    -type f -perm /222 -print -quit
require_empty_output /usr/bin/find "$lane" \
    -type d -perm /222 -print -quit
test "$(/usr/bin/git -C "$repo" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" = "$tree"
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/sha256sum \
    "$lane/audit_k16r16_b64_avx512_gfni_abba.py" | \
    /usr/bin/cut -d' ' -f1)" = "$auditor_hash"
/usr/bin/cmp "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    "$lane/build-closure/committed-auditor.py"

if [[ "$claim_passed" != true ]]; then
    /usr/bin/printf 'AUTHORITATIVE_STAGE publish_not_promoted_envelope\n'
    negative_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | /usr/bin/cut -d' ' -f1)"
    negative_audit_sha="$(/usr/bin/sha256sum "$envelope/postseal-audit.json" | /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg core_sha256sums_sha256 "$negative_core_sha" \
        --arg postseal_audit_sha256 "$negative_audit_sha" \
        '{schema:"leopard2-k16r16-b64-avx512-gfni-not-promoted-envelope/v1",status:"complete",promotion_passed:false,claim_passed:false,campaign_exit_status:$campaign_exit_status,postseal_audit_passed:true,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256,postseal_audit_sha256:$postseal_audit_sha256}' \
        > "$envelope/NOT_PROMOTED.json"
    (
        cd "$envelope"
        /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
            /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum \
            > SHA256SUMS
    )
    /usr/bin/find "$envelope" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$envelope" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    negative_verification="$(/usr/bin/mktemp -d /tmp/leopard-t16-negative.XXXXXX)"
    (
        cd "$envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) > "$negative_verification/sha256-check.txt"
    /usr/bin/printf 'sealed_not_promoted_envelope=%s\nverification=%s\n' \
        "$envelope" "$negative_verification"
    exit "$campaign_status"
fi

/usr/bin/printf 'AUTHORITATIVE_STAGE publish_completion_envelope\n'
core_manifest_sha="$(/usr/bin/sha256sum "$lane/manifest.json" | /usr/bin/cut -d' ' -f1)"
core_sha256sums_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | /usr/bin/cut -d' ' -f1)"
postseal_audit_sha="$(/usr/bin/sha256sum "$envelope/postseal-audit.json" | /usr/bin/cut -d' ' -f1)"
test "$campaign_status" -eq 0
test "$claim_passed" = true
completion_json="$(/usr/bin/jq -c -n \
    --argjson campaign_exit_status "$campaign_status" \
    --argjson claim_passed "$claim_passed" \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg binary_sha256 "$binary_hash" \
    --arg campaign_sha256 "$campaign_sha" \
    --arg audit_sha256 "$audit_sha" \
    --arg postseal_audit_sha256 "$postseal_audit_sha" \
    --arg core_manifest_sha256 "$core_manifest_sha" \
    --arg core_sha256sums_sha256 "$core_sha256sums_sha" \
    '{schema:"leopard2-k16r16-b64-avx512-gfni-completion-envelope/v1",status:"complete",promotion_passed:true,campaign_exit_status:$campaign_exit_status,claim_passed:$claim_passed,preseal_audit_passed:true,postseal_audit_passed:true,postseal_audit_byte_identical:true,core_sha256sums_verified:true,source_commit:$commit,source_tree:$tree,binary_sha256:$binary_sha256,campaign_sha256:$campaign_sha256,audit_sha256:$audit_sha256,postseal_audit_sha256:$postseal_audit_sha256,core_manifest_sha256:$core_manifest_sha256,core_sha256sums_sha256:$core_sha256sums_sha256,verification_command:["core/run-authoritative.sh","--verify","<absolute-envelope-path>"]}')"
completion_hash="$(/usr/bin/printf '%s\n' "$completion_json" | \
    /usr/bin/sha256sum | /usr/bin/cut -d' ' -f1)"
exec 8> "$envelope/COMPLETED.json"
(
    cd "$envelope"
    /usr/bin/find . -type f ! -path './SHA256SUMS' \
        ! -path './COMPLETED.json' -print0 | \
        /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum \
        > SHA256SUMS
    /usr/bin/printf '%s  ./COMPLETED.json\n' "$completion_hash" \
        >> SHA256SUMS
)
/usr/bin/find "$envelope" -type f -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$envelope" -type d -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/python3 -I -S -B -c \
    'import os, sys
data = (sys.argv[1] + "\n").encode()
view = memoryview(data)
while view:
    written = os.write(8, view)
    if written <= 0:
        raise RuntimeError("short completion write")
    view = view[written:]
os.fsync(8)
directory = os.open(sys.argv[2], os.O_RDONLY | os.O_DIRECTORY)
try:
    os.fsync(directory)
finally:
    os.close(directory)' \
    "$completion_json" "$envelope"
exec 8>&-

/usr/bin/printf 'AUTHORITATIVE_STAGE verify_completion_envelope\n'
verification_root="$(/usr/bin/mktemp -d /tmp/leopard-t16-envelope-verification.XXXXXX)"
"$lane/run-authoritative.sh" --verify "$envelope" \
    > "$verification_root/verification.txt" \
    2> "$verification_root/verification-stderr.log"

/usr/bin/printf 'AUTHORITATIVE_STAGE complete\n'
/usr/bin/jq '{claim_passed,gate_results,cells:[.cells[]|select(.role=="target")|{id,encode_execution,one_shot_encode}]}' \
    "$lane/campaign.json"
/usr/bin/printf 'sealed_envelope=%s\nexternal_verification=%s\n' \
    "$envelope" "$verification_root"
exit "$campaign_status"
