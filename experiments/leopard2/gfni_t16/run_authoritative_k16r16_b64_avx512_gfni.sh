#!/usr/bin/env bash
set -Eeuo pipefail
umask 077

repo=/home/catid/leopard
relative_wrapper=experiments/leopard2/gfni_t16/run_authoritative_k16r16_b64_avx512_gfni.sh
relative_runner=experiments/leopard2/gfni_t16/run_k16r16_b64_avx512_gfni_abba.py
relative_auditor=experiments/leopard2/gfni_t16/audit_k16r16_b64_avx512_gfni_abba.py
relative_validator=tools/leopard2_benchmark_json_test.py
lock=/tmp/leopard-gf8-authoritative.lock

if [[ $# -ne 1 || $1 != /* ]]; then
    /usr/bin/printf 'usage: %s /absolute/repository/.research/lane\n' "$0" >&2
    exit 2
fi
lane=$1
case "$lane" in
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

if [[ -e "$lane" ]]; then
    /usr/bin/printf 'refusing to reuse artifact lane: %s\n' "$lane" >&2
    exit 2
fi
/usr/bin/mkdir -p "$(/usr/bin/dirname "$lane")"
/usr/bin/mkdir -m 0700 "$lane"
trap failure_record ERR

next_stage canonical_lock
exec 9>> "$lock"
/usr/bin/flock -n 9
test -f "$lock"
test ! -L "$lock"
test "$(/usr/bin/stat -c %h "$lock")" = 1

next_stage clean_source_preflight
commit="$(/usr/bin/git -C "$repo" rev-parse HEAD)"
tree="$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')"
test "${#commit}" = 40
test "${#tree}" = 40
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal)"
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal --ignore-submodules=none)"
test "$(/usr/bin/git -C "$repo" rev-parse --show-toplevel)" = "$repo"
test "$(/usr/bin/readlink -f "$repo/$relative_wrapper")" = \
    "$(/usr/bin/readlink -f "$0")"
/usr/bin/git -C "$repo" ls-tree HEAD sse2neon > "$lane/live-sse2neon-gitlink.txt"
/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal \
    > "$lane/live-git-status.txt"
/usr/bin/printf '%s\n' "$commit" > "$lane/source-commit.txt"
/usr/bin/printf '%s\n' "$tree" > "$lane/source-tree.txt"
/usr/bin/printf '%s\n' "$(/usr/bin/cat /sys/devices/system/cpu/cpu13/topology/thread_siblings_list)" \
    > "$lane/cpu13-thread-siblings.txt"
test "$(/usr/bin/cat "$lane/cpu13-thread-siblings.txt")" = 13,77

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
test -z "$(/usr/bin/git -C "$replay_source" status --porcelain=v1 --untracked-files=normal)"
test -z "$(/usr/bin/git -C "$replay_source" status --porcelain=v1 --untracked-files=normal --ignore-submodules=none)"
/usr/bin/git -C "$replay_source" ls-tree HEAD sse2neon > "$lane/replay-sse2neon-gitlink.txt"
/usr/bin/cmp "$lane/live-sse2neon-gitlink.txt" "$lane/replay-sse2neon-gitlink.txt"

configure_common=(
    -G Ninja
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
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
focused_regex='^leopard2_(portable_isa|portable_isa_registration|portable_isa_checker_self_test|balanced_b64_terminal|balanced_b64_terminal_production|backend_ops|avx512_gfni_t16_prototype|benchmark_json_regression)$'

next_stage live_configure
/usr/bin/cmake -S "$repo" -B "$live_build" "${configure_common[@]}" \
    > "$lane/live-configure.log" 2>&1

next_stage live_build
/usr/bin/cmake --build "$live_build" --parallel 2 --target "${build_targets[@]}" \
    > "$lane/live-build.log" 2>&1

next_stage live_focused_tests
/usr/bin/ctest --test-dir "$live_build" -j1 --output-on-failure \
    --no-tests=error -R "$focused_regex" \
    > "$lane/live-focused-tests.log" 2>&1

next_stage replay_configure
/usr/bin/cmake -S "$replay_source" -B "$replay_build" "${configure_common[@]}" \
    > "$lane/replay-configure.log" 2>&1

next_stage replay_build
/usr/bin/cmake --build "$replay_build" --parallel 2 --target "${build_targets[@]}" \
    > "$lane/replay-build.log" 2>&1

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
    /usr/bin/c++ /usr/bin/cc /usr/bin/ld /usr/bin/ar /usr/bin/cmake \
    /usr/bin/ninja /usr/bin/git /usr/bin/objdump /usr/bin/readelf \
    /usr/bin/nm /usr/bin/python3 /usr/bin/taskset /usr/bin/flock \
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
library_hash="$(/usr/bin/sha256sum "$replay_build/libleopard.a" | /usr/bin/cut -d' ' -f1)"
t16_object_hash="$(/usr/bin/sha256sum "$replay_t16_object" | /usr/bin/cut -d' ' -f1)"
benchmark_object_hash="$(/usr/bin/sha256sum "$replay_benchmark_object" | /usr/bin/cut -d' ' -f1)"
attestation_hash="$(/usr/bin/sha256sum "$replay_attestation" | /usr/bin/cut -d' ' -f1)"
configuration_hash="$(/usr/bin/sha256sum "$replay_configuration" | /usr/bin/cut -d' ' -f1)"

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
    '{schema:"leopard2-k16r16-b64-avx512-gfni-build-closure/v1",repository:$repository,live_root:$live_root,live_build:$live_build,replay_root:$replay_root,replay_source:$replay_source,replay_build:$replay_build,source_commit:$commit,source_tree:$tree,byte_identical:{live_replay_binary:true,live_replay_prevalidated_binary:true,live_replay_library:true,live_replay_t16_object:true,live_replay_benchmark_object:true,live_replay_attestation:true,live_replay_configuration:true,live_replay_archive:true},frozen:{binary_sha256:$binary_sha256,library_sha256:$library_sha256,t16_object_sha256:$t16_object_sha256,benchmark_object_sha256:$benchmark_object_sha256,attestation_sha256:$attestation_sha256,configuration_sha256:$configuration_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256},configure:["/usr/bin/cmake","-G","Ninja","-DCMAKE_BUILD_TYPE=Release","-DCMAKE_EXPORT_COMPILE_COMMANDS=ON","-DLEO2_BUILD_TESTS=ON","-DLEO2_BUILD_BENCHMARKS=ON","-DLEO2_ENABLE_CUDA=OFF","-DLEO2_BACKEND_VARIANT=auto","-DLEOPARD_ENABLE_GF8=ON","-DLEOPARD_ENABLE_GF16=ON","-DLEO2_EXPERIMENT_HIGH_T16_B64_GENERATED=ON"],build:["/usr/bin/cmake","--build","<build>","--parallel","2","--target","bench_leopard2","bench_leopard2_prevalidated_batch","leopard2_balanced_b64_terminal_test","leopard2_balanced_b64_terminal_production_test","leopard2_backend_ops_test","leopard2_avx512_gfni_t16_prototype_test"],focused_test_regex:"^leopard2_(portable_isa|portable_isa_registration|portable_isa_checker_self_test|balanced_b64_terminal|balanced_b64_terminal_production|backend_ops|avx512_gfni_t16_prototype|benchmark_json_regression)$"}' \
    > "$lane/build-closure.json"

next_stage campaign
campaign_command=(
    /usr/bin/taskset -c 13
    /usr/bin/python3 "$lane/run_k16r16_b64_avx512_gfni_abba.py"
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
    --cpu 13
    --sibling 77
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
test "$(/usr/bin/sha256sum "$lane/bench_leopard2" | /usr/bin/cut -d' ' -f1)" = "$binary_hash"
test "$(/usr/bin/sha256sum "$lane/run_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)" = "$runner_hash"
test "$(/usr/bin/sha256sum "$lane/audit_k16r16_b64_avx512_gfni_abba.py" | /usr/bin/cut -d' ' -f1)" = "$auditor_hash"
test "$(/usr/bin/sha256sum "$lane/leopard2_benchmark_json_test.py" | /usr/bin/cut -d' ' -f1)" = "$validator_hash"
test "$(/usr/bin/sha256sum "$lane/source.tar" | /usr/bin/cut -d' ' -f1)" = "$archive_hash"
test "$(/usr/bin/git -C "$repo" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" = "$tree"
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal)"
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal --ignore-submodules=none)"
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
        '{schema:"leopard2-k16r16-b64-avx512-gfni-failed-manifest/v1",campaign_exit_status:$campaign_exit_status,report_status:"failed",claim_passed:false,audit_performed:false,source_commit:$commit,source_tree:$tree,binary_sha256_pre_post:$binary_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256,campaign_sha256:$report_sha256,journal_sha256:($journal_sha256 | if . == "null" then null else . end),canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:13,sibling:77}' \
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
    test -z "$(/usr/bin/find "$lane" -type l -print -quit)"
    test -z "$(/usr/bin/find "$lane" -type f -links +1 -print -quit)"
    test -z "$(/usr/bin/find "$lane" -type f -perm /222 -print -quit)"
    test -z "$(/usr/bin/find "$lane" -type d -perm /222 -print -quit)"
    /usr/bin/printf 'sealed_failed_lane=%s\nverification=%s\n' \
        "$lane" "$failed_verification"
    if [[ "$campaign_status" -eq 0 ]]; then
        exit 1
    fi
    exit "$campaign_status"
fi
test "$report_status" = complete

next_stage independent_preseal_audit
/usr/bin/python3 "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
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
    '{schema:"leopard2-k16r16-b64-avx512-gfni-final-manifest/v1",campaign_exit_status:$campaign_exit_status,claim_passed:$claim_passed,audit_passed:true,source_commit:$commit,source_tree:$tree,binary_sha256_pre_post:$binary_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,validator_sha256:$validator_sha256,wrapper_sha256:$wrapper_sha256,source_archive_sha256:$archive_sha256,campaign_sha256:$campaign_sha256,journal_sha256:$journal_sha256,audit_sha256:$audit_sha256,replay_root:$replay_root,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:13,sibling:77,postseal_policy:"rerun the independent auditor outside the sealed lane, require byte-identical audit output, and verify SHA256SUMS while the live repository remains exact and clean"}' \
    > "$lane/manifest.json"

next_stage seal
trap - ERR
(
    cd "$lane"
    /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
        /usr/bin/sort -z | /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
)
/usr/bin/find "$lane" -type f -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$lane" -type d -perm /222 -exec /usr/bin/chmod a-w {} +

/usr/bin/printf 'AUTHORITATIVE_STAGE independent_postseal_audit\n'
postseal_root="$(/usr/bin/mktemp -d /tmp/leopard-t16-postseal.XXXXXX)"
/usr/bin/python3 "$lane/audit_k16r16_b64_avx512_gfni_abba.py" \
    --report "$lane/campaign.json" \
    --journal "$lane/campaign.jsonl" \
    --output "$postseal_root/audit.json" \
    > "$postseal_root/audit-summary.json" 2> "$postseal_root/audit-stderr.log"
/usr/bin/cmp "$lane/audit.json" "$postseal_root/audit.json"
(
    cd "$lane"
    /usr/bin/sha256sum -c SHA256SUMS
) > "$postseal_root/sha256-check.txt"
test -z "$(/usr/bin/find "$lane" -type l -print -quit)"
test -z "$(/usr/bin/find "$lane" -type f -links +1 -print -quit)"
test -z "$(/usr/bin/find "$lane" -type f -perm /222 -print -quit)"
test -z "$(/usr/bin/find "$lane" -type d -perm /222 -print -quit)"
test "$(/usr/bin/git -C "$repo" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" = "$tree"
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal)"
test -z "$(/usr/bin/git -C "$repo" status --porcelain=v1 --untracked-files=normal --ignore-submodules=none)"

/usr/bin/printf 'AUTHORITATIVE_STAGE complete\n'
/usr/bin/jq '{claim_passed,gate_results,cells:[.cells[]|select(.role=="target")|{id,encode_execution,one_shot_encode}]}' \
    "$lane/campaign.json"
/usr/bin/printf 'sealed_lane=%s\npostseal_verification=%s\n' \
    "$lane" "$postseal_root"
exit "$campaign_status"
