#!/usr/bin/bash
set -Eeuo pipefail
umask 077

if [[ ${1:-} != --leopard-v17-gfni-main-clean-env-internal ]]; then
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
        SOURCE_DATE_EPOCH=0 \
        CC=/usr/bin/cc \
        CXX=/usr/bin/c++ \
        CFLAGS= \
        CXXFLAGS= \
        CPPFLAGS= \
        LDFLAGS= \
        /usr/bin/bash "$authoritative_script" \
        --leopard-v17-gfni-main-clean-env-internal "$@"
fi
shift

passive_mode=false
attempt=
attempt_budget=
if [[ ${1:-} == --passive-shared-host ]]; then
    passive_mode=true
    shift
    if [[ ${1:-} != --attempt || ${3:-} != --attempt-budget ]]; then
        /usr/bin/printf '%s requires --attempt N --attempt-budget 3\n' \
            --passive-shared-host >&2
        exit 2
    fi
    attempt=${2:-}
    attempt_budget=${4:-}
    shift 4
fi

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
    'SOURCE_DATE_EPOCH=0' \
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
        PATH|LANG|LC_ALL|TZ|TMPDIR|OMP_DYNAMIC|OMP_NUM_THREADS|OMP_THREAD_LIMIT|PYTHONNOUSERSITE|PYTHONDONTWRITEBYTECODE|SOURCE_DATE_EPOCH|CC|CXX|CFLAGS|CXXFLAGS|CPPFLAGS|LDFLAGS|PWD|SHLVL|_)
            ;;
        *)
            /usr/bin/printf 'unexpected authoritative environment key: %s\n' \
                "$environment_name" >&2
            exit 2
            ;;
    esac
done < <(/usr/bin/env -0)

repo=/home/catid/leopard
main_commit=6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198
relative_wrapper=experiments/leopard2/main_compare/run_authoritative_v17_gfni_main_compare.sh
relative_runner=experiments/leopard2/main_compare/run_abba.py
relative_git_capture=experiments/leopard2/main_compare/git_capture.py
relative_helper=experiments/leopard2/decoder_dispatch/balanced_evidence_common.py
relative_build_provenance=tools/leopard2_build_provenance.py
relative_auditor=experiments/leopard2/main_compare/audit_v17_gfni_main_compare.py
relative_census=experiments/leopard2/main_compare/passive_environment_census.py
relative_supervisor=tools/leopard2_affinity_supervisor.py
lock=/tmp/leopard-gf8-authoritative.lock
passive_housekeeping_jq='[.allowed_cpus[] | select(. != 52 and . != 116)]
    | map(tostring) | join(",") | select(length > 0)'
passive_timeout_argument_jq='.campaign.timeout_seconds
    | select(type == "number" and . == 600)
    | "600"'

validate_attempt_contract()
{
    local observed_attempt=$1
    local observed_budget=$2
    [[ "$observed_attempt" =~ ^[1-3]$ ]] || return 1
    test "$observed_budget" = 3
}

require_path_absent()
{
    test ! -e "$1" || return 1
    test ! -L "$1" || return 1
    return 0
}

v18_envelope_has_observational_output()
{
    local observed_envelope=$1
    local observed_core="$observed_envelope/core"
    local observed_path=
    for observed_path in \
        campaign/raw.json \
        campaign/manifest.json \
        campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374 \
        campaign/.leopard2-evidence-pending-ffa5b716b5a57837f7929dfcca4b4dfdeb97210a7fd5a12d2f1978846d6f1743 \
        audit.json \
        result-summary.json \
        candidate-test-temporal-closure.json \
        passive-environment-policy.json; do
        if [[ -e "$observed_core/$observed_path" ||
              -L "$observed_core/$observed_path" ]]; then
            return 0
        fi
    done
    if [[ -e "$observed_envelope/postseal-audit.json" ||
          -L "$observed_envelope/postseal-audit.json" ]]; then
        return 0
    fi
    if [[ -f "$observed_core/manifest.json" ]] &&
       /usr/bin/jq -e -s '
           length == 1 and
           .[0].schema ==
               "leopard2-v18-gfni-main-passive-core-manifest/v1" and
           .[0].status == "complete" and .[0].evidence_valid == true
       ' "$observed_core/manifest.json" >/dev/null 2>&1; then
        return 0
    fi
    return 1
}

validate_single_json_value()
{
    local json_path=$1
    local expected_kind=$2
    /usr/bin/python3 -I -S -B -c '
import json
import math
import sys

def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON object key")
        result[key] = value
    return result

def reject_constant(value):
    raise ValueError(f"non-finite JSON number: {value}")

def finite_float(token):
    value = float(token)
    if not math.isfinite(value):
        raise ValueError(f"non-finite JSON number: {token}")
    return value

limit = 16 * 1024 * 1024
chunks = []
size = 0
with open(sys.argv[1], "rb", buffering=0) as stream:
    while True:
        chunk = stream.read(min(65536, limit + 1 - size))
        if not chunk:
            break
        chunks.append(chunk)
        size += len(chunk)
        if size > limit:
            raise ValueError("wrapper-owned JSON exceeds the validation bound")
data = b"".join(chunks)
value = json.loads(
    data.decode("utf-8"), object_pairs_hook=unique_object,
    parse_constant=reject_constant, parse_float=finite_float)
expected = {"object": dict, "array": list}.get(sys.argv[2])
if expected is None or type(value) is not expected:
    raise ValueError(f"wrapper-owned JSON is not one {sys.argv[2]}")
' "$json_path" "$expected_kind" >/dev/null 2>&1 || return 1
    return 0
}

validate_single_json_object()
{
    validate_single_json_value "$1" object
}

validate_single_json_array()
{
    validate_single_json_value "$1" array
}

validate_exact_json_schema()
{
    local json_path=$1
    local expected_schema=$2
    /usr/bin/jq -e --arg expected_schema "$expected_schema" \
        '.schema == $expected_schema' "$json_path" >/dev/null
}

validate_v18_terminal_static()
{
    local status_path=$1
    local generation=$2
    validate_single_json_object "$status_path" || return 1
    if [[ "$generation" == passive-v2 ]]; then
        /usr/bin/jq -e '
            .schema ==
                "leopard2-v18-gfni-main-passive-not-promoted-envelope/v1" and
            .status == "complete" and
            .acquisition_generation == "passive-v2" and
            (.attempt | type) == "number" and
            .attempt == (.attempt | floor) and
            .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
            (.attempt_lineage_sha256 | type) == "string" and
            (.attempt_lineage_sha256 | length) == 64 and
            (.attempt_lineage_sha256 | test("^[0-9a-f]{64}$")) and
            .promotion_passed == false and .campaign_exit_status == 0 and
            (.source_commit | type) == "string" and
            (.source_commit | length) == 40 and
            (.source_commit | test("^[0-9a-f]{40}$")) and
            (.source_tree | type) == "string" and
            (.source_tree | length) == 40 and
            (.source_tree | test("^[0-9a-f]{40}$")) and
            .baseline_commit ==
                "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
            all([.campaign_manifest_sha256,.audit_sha256,
                .postseal_audit_sha256,.core_manifest_sha256,
                .core_sha256sums_sha256][];
                type == "string" and length == 64 and
                test("^[0-9a-f]{64}$"))
        ' "$status_path" >/dev/null || return 1
    elif [[ "$generation" == failed-v2 ]]; then
        /usr/bin/jq -e '
            .schema == "leopard2-v18-gfni-main-failed-envelope/v1" and
            .status == "failed" and .acquisition_generation == "passive-v2" and
            (.attempt | type) == "number" and
            .attempt == (.attempt | floor) and
            .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
            (.attempt_lineage_sha256 | type) == "string" and
            (.attempt_lineage_sha256 | length) == 64 and
            (.attempt_lineage_sha256 | test("^[0-9a-f]{64}$")) and
            .promotion_passed == false and
            (.campaign_exit_status | type) == "number" and
            .campaign_exit_status == (.campaign_exit_status | floor) and
            .campaign_exit_status >= 1 and .campaign_exit_status <= 255 and
            (.source_commit | type) == "string" and
            (.source_commit | length) == 40 and
            (.source_commit | test("^[0-9a-f]{40}$")) and
            (.source_tree | type) == "string" and
            (.source_tree | length) == 40 and
            (.source_tree | test("^[0-9a-f]{40}$")) and
            (.core_sha256sums_sha256 | type) == "string" and
            (.core_sha256sums_sha256 | length) == 64 and
            (.core_sha256sums_sha256 | test("^[0-9a-f]{64}$"))
        ' "$status_path" >/dev/null || return 1
    else
        return 1
    fi
    return 0
}

validate_v18_wrapper_failure_binding()
{
    local retained_failure=$1
    local outer_failure=$2
    validate_single_json_object "$retained_failure" || return 1
    validate_single_json_object "$outer_failure" || return 1
    /usr/bin/jq -e --slurpfile status "$outer_failure" '
        ($status | length) == 1 and
        keys == (["acquisition_generation","attempt","attempt_budget",
            "attempt_lineage_sha256","exit_status","promotion_passed",
            "schema","source_commit","source_tree","stage","status"] | sort) and
        .schema ==
            "leopard2-v18-gfni-main-passive-wrapper-failure/v1" and
        .acquisition_generation == "passive-v2" and
        (.attempt | type) == "number" and .attempt == (.attempt | floor) and
        .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
        (.attempt_lineage_sha256 | type) == "string" and
        (.attempt_lineage_sha256 | length) == 64 and
        (.attempt_lineage_sha256 | test("^[0-9a-f]{64}$")) and
        .status == "failed" and .promotion_passed == false and
        (.exit_status | type) == "number" and
        .exit_status == (.exit_status | floor) and
        .exit_status >= 1 and .exit_status <= 255 and
        (.stage | type) == "string" and (.stage | length) > 0 and
        (.source_commit | type) == "string" and
        (.source_commit | length) == 40 and
        (.source_commit | test("^[0-9a-f]{40}$")) and
        (.source_tree | type) == "string" and
        (.source_tree | length) == 40 and
        (.source_tree | test("^[0-9a-f]{40}$")) and
        .attempt == $status[0].attempt and
        .attempt_budget == $status[0].attempt_budget and
        .attempt_lineage_sha256 == $status[0].attempt_lineage_sha256 and
        .exit_status == $status[0].campaign_exit_status and
        .stage == $status[0].stage and
        .source_commit == $status[0].source_commit and
        .source_tree == $status[0].source_tree
    ' "$retained_failure" >/dev/null
}

validate_v18_complete_core_static()
{
    validate_single_json_object "$1" || return 1
    /usr/bin/jq -e '
        keys == ([
            "acquisition_generation", "active_affinity_supervisor_executed",
            "attempt", "attempt_budget", "attempt_lineage_sha256",
            "audit_sha256", "auditor_sha256",
            "baseline_archive_sha256_pre_post",
            "baseline_binary_sha256_pre_post", "baseline_commit",
            "build_order_normalizer_sha256", "build_reproduction",
            "campaign_exit_status", "campaign_manifest_sha256",
            "campaign_raw_sha256", "candidate_archive_sha256_pre_post",
            "candidate_binary_sha256_pre_post",
            "candidate_test_normalization_report_sha256",
            "candidate_timing_normalization_report_sha256",
            "canonical_lock", "causal_performance_claim_eligible",
            "controller_affinity_sha256", "cpu", "cpu_pair_exclusive",
            "environment_census_post_sha256",
            "environment_census_pre_sha256", "evidence_class",
            "evidence_valid", "independent_auditor_scope",
            "independent_auditor_supervision_mode",
            "independent_preseal_audit_passed", "isolation_claim",
            "passive_census_sha256", "passive_environment_policy_sha256",
            "performance_gate_passed", "postseal_policy",
            "producer_verification_mode", "producer_verification_passed",
            "promotion_eligible", "promotion_passed", "ratio_policy",
            "runner_sha256", "schema", "sibling", "source_commit",
            "source_tree", "sse2neon_commit",
            "sse2neon_source_archive_sha256", "status",
            "supervisor_role", "supervisor_sha256",
            "out_of_window_benchmark_cpu_nonidle_jiffies",
            "out_of_window_reserved_sibling_nonidle_jiffies",
            "windowed_benchmark_cpu_nonidle_excess_jiffies",
            "windowed_reserved_sibling_nonidle_jiffies", "wrapper_sha256"
        ] | sort) and
        .schema == "leopard2-v18-gfni-main-passive-core-manifest/v1" and
        .status == "complete" and .acquisition_generation == "passive-v2" and
        (.attempt | type) == "number" and .attempt == (.attempt | floor) and
        .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
        (.attempt_lineage_sha256 | type) == "string" and
        (.attempt_lineage_sha256 | length) == 64 and
        (.attempt_lineage_sha256 | test("^[0-9a-f]{64}$")) and
        .campaign_exit_status == 0 and .evidence_valid == true and
        .evidence_class == "passive-windowed-shared-host-observation/v1" and
        (.performance_gate_passed | type) == "boolean" and
        .promotion_eligible == false and .promotion_passed == false and
        .causal_performance_claim_eligible == false and
        .cpu_pair_exclusive == false and
        .windowed_benchmark_cpu_nonidle_excess_jiffies == 0 and
        .windowed_reserved_sibling_nonidle_jiffies == 0 and
        (.out_of_window_benchmark_cpu_nonidle_jiffies | type) == "number" and
        .out_of_window_benchmark_cpu_nonidle_jiffies ==
            (.out_of_window_benchmark_cpu_nonidle_jiffies | floor) and
        .out_of_window_benchmark_cpu_nonidle_jiffies >= 0 and
        (.out_of_window_reserved_sibling_nonidle_jiffies | type) == "number" and
        .out_of_window_reserved_sibling_nonidle_jiffies ==
            (.out_of_window_reserved_sibling_nonidle_jiffies | floor) and
        .out_of_window_reserved_sibling_nonidle_jiffies >= 0 and
        (.source_commit | type) == "string" and
        (.source_commit | length) == 40 and
        (.source_commit | test("^[0-9a-f]{40}$")) and
        (.source_tree | type) == "string" and
        (.source_tree | length) == 40 and
        (.source_tree | test("^[0-9a-f]{40}$")) and
        .baseline_commit == "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
        (.sse2neon_commit | type) == "string" and
        (.sse2neon_commit | length) == 40 and
        (.sse2neon_commit | test("^[0-9a-f]{40}$")) and
        all([
            .sse2neon_source_archive_sha256,
            .candidate_binary_sha256_pre_post,
            .candidate_archive_sha256_pre_post,
            .baseline_binary_sha256_pre_post,
            .baseline_archive_sha256_pre_post,
            .runner_sha256, .auditor_sha256, .supervisor_sha256,
            .passive_census_sha256, .wrapper_sha256,
            .build_order_normalizer_sha256,
            .candidate_timing_normalization_report_sha256,
            .candidate_test_normalization_report_sha256,
            .campaign_manifest_sha256, .campaign_raw_sha256, .audit_sha256,
            .controller_affinity_sha256, .environment_census_pre_sha256,
            .environment_census_post_sha256,
            .passive_environment_policy_sha256
        ][]; type == "string" and length == 64 and
            test("^[0-9a-f]{64}$")) and
        .isolation_claim == {
            schema:"leopard2-v18-passive-windowed-isolation-claim/v1",
            mechanism:(
                "owned-process taskset pinning, reservation file, pair " +
                "lease, and exact per-invocation CPU52/CPU116 jiffy " +
                "rejection screens"
            ),
            campaign_supervision:null,
            foreign_process_affinity_mutation_claimed:false,
            foreign_process_signalling_claimed:false,
            benchmark_cpu_exclusive_ownership_claimed:false,
            same_uid_pair_exclusion_certificate:false,
            cgroup_or_os_exclusive_certificate:false,
            every_retained_window_benchmark_cpu_excess_zero:true,
            every_retained_window_reserved_sibling_nonidle_zero:true,
            whole_campaign_and_out_of_window_activity_gated:false,
            windowed_screen:"per retained benchmark invocation",
            out_of_window_activity_gated:false,
            reserved_sibling_zero_nonidle_jiffies_in_every_retained_window:true,
            counter_resolution:(
                "frozen x86-64 Linux /proc/stat USER_HZ=100 contract"
            ),
            interval_complete_task_observation:false,
            benchmark_cpu_foreign_work_attributable:false,
            causal_performance_claim_eligible:false,
            promotion_eligible:false,
            promotion_passed:false,
            unmitigated_confounders:[
                "foreign CPU52 work hidden below the owned-wall ceiling",
                "package boost and thermal residency",
                "shared LLC and memory bandwidth",
                "root-owned, other-UID, and kernel work",
                "sub-jiffy transient work"
            ],
            scope:(
                "host, compiler, API, and workload specific; noncausal, " +
                "nonexclusive, not generalizable, and never promotion " +
                "evidence"
            )
        } and
        .build_reproduction == {
            candidate_timing_qualified_semantic_equal:true,
            candidate_timing_raw_tree_file_bytes_identical:
                .build_reproduction.candidate_timing_raw_tree_file_bytes_identical,
            candidate_tests_qualified_semantic_equal:true,
            candidate_tests_raw_tree_file_bytes_identical:
                .build_reproduction.candidate_tests_raw_tree_file_bytes_identical,
            baseline_raw_byte_identical:true
        } and
        (.build_reproduction.candidate_timing_raw_tree_file_bytes_identical |
            type) == "boolean" and
        (.build_reproduction.candidate_tests_raw_tree_file_bytes_identical |
            type) == "boolean" and
        .producer_verification_passed == true and
        .producer_verification_mode == "manifest-without-affinity-binding" and
        .independent_preseal_audit_passed == true and
        .independent_auditor_supervision_mode == "windowed" and
        .independent_auditor_scope ==
            "campaign semantics only; build-order qualification is separately recomputed by the frozen wrapper normalizer" and
        .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
        .cpu == 52 and .sibling == 116 and
        .supervisor_role == "retained-active-v1-verifier-only" and
        .active_affinity_supervisor_executed == false and
        .ratio_policy == {
            ordinary_and_one_shot_are_separate:true,
            ratios_are_separate_correlated_and_must_not_be_multiplied:true,
            combined_or_stacked_ratio_emitted:false,
            same_binary_ratio_is_another_campaign:true,
            clustered_t_interval_reported_as_nominal_under_shared_host_load:true
        } and
        .postseal_policy ==
            "passive-v2 shared-host observations always publish NOT_PROMOTED.json, independent of the observed speed gate; at most three preregistered attempts and every outcome is retained"
    ' "$1" >/dev/null
}

validate_v18_failed_core_static()
{
    validate_single_json_object "$1" || return 1
    /usr/bin/jq -e '
        keys == (["acquisition_generation","attempt","attempt_budget",
            "attempt_lineage_sha256","baseline_binary_sha256",
            "baseline_commit","campaign_exit_status",
            "candidate_binary_sha256","canonical_lock","cpu",
            "failure_sha256","failure_verified","failure_verify_status",
            "promotion_passed","schema","sibling","source_commit",
            "source_tree","status"] | sort) and
        .schema ==
            "leopard2-v18-gfni-main-passive-failed-core-manifest/v1" and
        .status == "failed" and .acquisition_generation == "passive-v2" and
        (.attempt | type) == "number" and .attempt == (.attempt | floor) and
        .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
        (.attempt_lineage_sha256 | type) == "string" and
        (.attempt_lineage_sha256 | length) == 64 and
        (.attempt_lineage_sha256 | test("^[0-9a-f]{64}$")) and
        .promotion_passed == false and
        (.campaign_exit_status | type) == "number" and
        .campaign_exit_status == (.campaign_exit_status | floor) and
        .campaign_exit_status >= 1 and .campaign_exit_status <= 255 and
        (.failure_verify_status | type) == "number" and
        .failure_verify_status == (.failure_verify_status | floor) and
        .failure_verify_status >= 0 and .failure_verify_status <= 255 and
        (.failure_verified | type) == "boolean" and
        .failure_verified == (.failure_verify_status == 0) and
        (.source_commit | type) == "string" and
        (.source_commit | length) == 40 and
        (.source_commit | test("^[0-9a-f]{40}$")) and
        (.source_tree | type) == "string" and
        (.source_tree | length) == 40 and
        (.source_tree | test("^[0-9a-f]{40}$")) and
        .baseline_commit == "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
        all([.candidate_binary_sha256,.baseline_binary_sha256][];
            type == "string" and length == 64 and
            test("^[0-9a-f]{64}$")) and
        ((.failure_sha256 == null and .failure_verify_status == 1 and
             .failure_verified == false) or
            ((.failure_sha256 | type) == "string" and
             (.failure_sha256 | length) == 64 and
             (.failure_sha256 | test("^[0-9a-f]{64}$")))) and
        .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
        .cpu == 52 and .sibling == 116
    ' "$1" >/dev/null
}

v18_failed_core_static_self_test()
(
    local self_test_root=
    local self_test_manifest=
    local command_status=0
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-failed-core-static-self-test.XXXXXX)" || exit 1
    self_test_manifest="$self_test_root/manifest.json"
    cleanup_failed_core_static_self_test()
    {
        command_status=$?
        trap - EXIT
        if ! /usr/bin/find "$self_test_root" -depth -delete; then
            command_status=1
        fi
        exit "$command_status"
    }
    trap cleanup_failed_core_static_self_test EXIT
    /usr/bin/jq -n \
        '{schema:"leopard2-v18-gfni-main-passive-failed-core-manifest/v1",status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,failure_verify_status:1,failure_verified:false,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",candidate_binary_sha256:("c" * 64),baseline_binary_sha256:("e" * 64),failure_sha256:null,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$self_test_manifest" || exit 1
    validate_v18_failed_core_static "$self_test_manifest" || exit 1
    /usr/bin/jq -n \
        '{schema:"leopard2-v18-gfni-main-passive-failed-core-manifest/v1",status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,failure_verify_status:0,failure_verified:true,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",candidate_binary_sha256:("c" * 64),baseline_binary_sha256:("e" * 64),failure_sha256:("f" * 64),canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$self_test_manifest" || exit 1
    validate_v18_failed_core_static "$self_test_manifest" || exit 1
    /usr/bin/jq -n \
        '{schema:"leopard2-v18-gfni-main-passive-failed-core-manifest/v1",status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,failure_verify_status:0,failure_verified:true,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",candidate_binary_sha256:("c" * 64),baseline_binary_sha256:("e" * 64),failure_sha256:null,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$self_test_manifest" || exit 1
    if validate_v18_failed_core_static "$self_test_manifest"; then
        exit 1
    fi
    /usr/bin/jq -n \
        '{schema:"leopard2-v18-gfni-main-passive-failed-core-manifest/v1",status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,failure_verify_status:255,failure_verified:false,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",candidate_binary_sha256:("c" * 64),baseline_binary_sha256:("e" * 64),failure_sha256:null,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$self_test_manifest" || exit 1
    if validate_v18_failed_core_static "$self_test_manifest"; then
        exit 1
    fi
    exit 0
)

v18_observational_output_self_test()
(
    local self_test_root=
    local command_status=0
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-observational-output-self-test.XXXXXX)" || exit 1
    cleanup_observational_output_self_test()
    {
        command_status=$?
        trap - EXIT
        if ! /usr/bin/find "$self_test_root" -depth -delete; then
            command_status=1
        fi
        exit "$command_status"
    }
    trap cleanup_observational_output_self_test EXIT
    /usr/bin/mkdir -p "$self_test_root/core/campaign" || exit 1
    if v18_envelope_has_observational_output "$self_test_root"; then
        exit 1
    fi
    /usr/bin/printf '%s\n' '{"diagnostic":"raw"}' \
        > "$self_test_root/core/campaign/raw.json" || exit 1
    v18_envelope_has_observational_output "$self_test_root" || exit 1
    /usr/bin/rm "$self_test_root/core/campaign/raw.json" || exit 1
    /usr/bin/printf '%s\n' '{"diagnostic":"pending raw"}' \
        > "$self_test_root/core/campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374" || \
        exit 1
    v18_envelope_has_observational_output "$self_test_root" || exit 1
    /usr/bin/rm \
        "$self_test_root/core/campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374" || \
        exit 1
    /usr/bin/printf '%s\n' \
        '{"schema":"leopard2-v18-gfni-main-passive-core-manifest/v1","status":"complete","evidence_valid":true}' \
        > "$self_test_root/core/manifest.json" || exit 1
    v18_envelope_has_observational_output "$self_test_root" || exit 1
    exit 0
)

verify_core_file_at_source_commit()
{
    local core_path=$1
    local source_commit=$2
    local relative_path=$3
    local expected_hash=
    local observed_hash=
    [[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || return 1
    test -f "$core_path" || return 1
    test ! -L "$core_path" || return 1
    expected_hash="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        /usr/bin/git -C "$repo" cat-file blob \
        "$source_commit:$relative_path" | /usr/bin/sha256sum | \
        /usr/bin/cut -d' ' -f1)" || return 1
    observed_hash="$(/usr/bin/sha256sum "$core_path" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    test "$observed_hash" = "$expected_hash" || return 1
    return 0
}

verify_source_tree_at_commit()
{
    local source_commit=$1
    local source_tree=$2
    local observed_commit=
    local observed_tree=
    [[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || return 1
    [[ "$source_tree" =~ ^[0-9a-f]{40}$ ]] || return 1
    observed_commit="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        /usr/bin/git -C "$repo" rev-parse --verify \
        "$source_commit^{commit}")" || return 1
    test "$observed_commit" = "$source_commit" || return 1
    observed_tree="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        /usr/bin/git -C "$repo" rev-parse --verify \
        "$source_commit^{tree}")" || return 1
    test "$observed_tree" = "$source_tree" || return 1
    return 0
}

canonical_git_archive_stream()
(
    local source_repo=$1
    local source_commit=$2
    local archive_prefix=$3
    local archive_root=
    local source_common_dir=
    local source_objects=
    local observed_commit=
    local command_status=0
    [[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || exit 1
    case "$archive_prefix" in
        candidate-source|leopard1-source|sse2neon-source) ;;
        *) exit 1 ;;
    esac
    test -d "$source_repo" || exit 1
    archive_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-canonical-git-archive.XXXXXX)" || exit 1
    cleanup_canonical_archive()
    {
        command_status=$?
        trap - EXIT
        if ! /usr/bin/find "$archive_root" -depth -delete; then
            command_status=1
        fi
        exit "$command_status"
    }
    trap cleanup_canonical_archive EXIT
    /usr/bin/mkdir "$archive_root/empty-template" || exit 1
    source_common_dir="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git -C "$source_repo" rev-parse \
        --path-format=absolute --git-common-dir)" || exit 1
    source_objects="$(/usr/bin/readlink -f \
        "$source_common_dir/objects")" || exit 1
    test -d "$source_objects" || exit 1
    case "$source_objects" in
        *:*) exit 1 ;;
    esac
    /usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git init --bare --quiet \
        --template="$archive_root/empty-template" \
        "$archive_root/repository.git" || exit 1
    require_path_absent \
        "$archive_root/repository.git/info/attributes" || exit 1
    observed_commit="$(/usr/bin/env \
        GIT_DIR="$archive_root/repository.git" \
        GIT_ALTERNATE_OBJECT_DIRECTORIES="$source_objects" \
        GIT_ATTR_SOURCE="$source_commit" \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git -c core.attributesFile=/dev/null \
        -c tar.umask=0002 rev-parse --verify \
        "$source_commit^{commit}")" || exit 1
    test "$observed_commit" = "$source_commit" || exit 1
    /usr/bin/env \
        GIT_DIR="$archive_root/repository.git" \
        GIT_ALTERNATE_OBJECT_DIRECTORIES="$source_objects" \
        GIT_ATTR_SOURCE="$source_commit" \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git -c core.attributesFile=/dev/null \
        -c tar.umask=0002 archive --format=tar \
        --prefix="$archive_prefix/" "$source_commit"
)

write_canonical_git_archive()
{
    local source_repo=$1
    local source_commit=$2
    local archive_prefix=$3
    local archive_path=$4
    local archive_parent=
    local archive_root=
    local staged_archive=
    [[ "$archive_path" == /* ]] || return 1
    archive_parent=${archive_path%/*}
    test -d "$archive_parent" || return 1
    test ! -L "$archive_parent" || return 1
    require_path_absent "$archive_path" || return 1
    archive_root="$(/usr/bin/mktemp -d \
        "$archive_parent/.leopard-canonical-git-archive-output.XXXXXX")" || \
        return 1
    staged_archive="$archive_root/archive.tar"
    if ! canonical_git_archive_stream \
            "$source_repo" "$source_commit" "$archive_prefix" \
            > "$staged_archive"; then
        /usr/bin/find "$archive_root" -depth -delete || true
        return 1
    fi
    require_path_absent "$archive_path" || {
        /usr/bin/find "$archive_root" -depth -delete || true
        return 1
    }
    if ! /usr/bin/mv -T -- "$staged_archive" "$archive_path"; then
        /usr/bin/find "$archive_root" -depth -delete || true
        return 1
    fi
    /usr/bin/find "$archive_root" -depth -delete || return 1
    return 0
}

verify_source_archive_at_commit()
{
    local source_repo=$1
    local source_commit=$2
    local archive_prefix=$3
    local archive_path=$4
    local expected_hash=
    local observed_hash=
    test -f "$archive_path" || return 1
    test ! -L "$archive_path" || return 1
    expected_hash="$(canonical_git_archive_stream \
        "$source_repo" "$source_commit" "$archive_prefix" | \
        /usr/bin/sha256sum | /usr/bin/cut -d' ' -f1)" || return 1
    observed_hash="$(/usr/bin/sha256sum "$archive_path" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    test "$observed_hash" = "$expected_hash" || return 1
    return 0
}

canonical_git_archive_self_test()
(
    local source_commit=$1
    local self_test_root=
    local source_common_dir=
    local source_objects=
    local clean_hash=
    local poisoned_hash=
    local isolated_hash=
    local command_status=0
    [[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || exit 1
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-canonical-git-archive-self-test.XXXXXX)" || exit 1
    cleanup_canonical_archive_self_test()
    {
        command_status=$?
        trap - EXIT
        if ! /usr/bin/find "$self_test_root" -depth -delete; then
            command_status=1
        fi
        exit "$command_status"
    }
    trap cleanup_canonical_archive_self_test EXIT
    /usr/bin/mkdir "$self_test_root/empty-template" || exit 1
    source_common_dir="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git -C "$repo" rev-parse \
        --path-format=absolute --git-common-dir)" || exit 1
    source_objects="$(/usr/bin/readlink -f \
        "$source_common_dir/objects")" || exit 1
    test -d "$source_objects" || exit 1
    /usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git init --bare --quiet \
        --template="$self_test_root/empty-template" \
        "$self_test_root/poisoned.git" || exit 1
    /usr/bin/mkdir -p "$self_test_root/poisoned.git/info" \
        "$self_test_root/poisoned.git/objects/info" || exit 1
    /usr/bin/printf '%s\n' "$source_objects" \
        > "$self_test_root/poisoned.git/objects/info/alternates" || exit 1
    /usr/bin/git --git-dir="$self_test_root/poisoned.git" \
        config tar.umask 0077 || exit 1
    /usr/bin/printf '%s\n' '* export-ignore' \
        > "$self_test_root/poisoned.git/info/attributes" || exit 1
    clean_hash="$(canonical_git_archive_stream \
        "$repo" "$source_commit" candidate-source | \
        /usr/bin/sha256sum | /usr/bin/cut -d' ' -f1)" || exit 1
    poisoned_hash="$(/usr/bin/env \
        GIT_DIR="$self_test_root/poisoned.git" \
        GIT_ATTR_SOURCE="$source_commit" \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        GIT_ATTR_NOSYSTEM=1 \
        /usr/bin/git -c core.attributesFile=/dev/null \
        archive --format=tar --prefix=candidate-source/ \
        "$source_commit" | /usr/bin/sha256sum | \
        /usr/bin/cut -d' ' -f1)" || exit 1
    test "$poisoned_hash" != "$clean_hash" || exit 1
    isolated_hash="$(canonical_git_archive_stream \
        "$self_test_root/poisoned.git" "$source_commit" candidate-source | \
        /usr/bin/sha256sum | /usr/bin/cut -d' ' -f1)" || exit 1
    test "$isolated_hash" = "$clean_hash" || exit 1
    /usr/bin/mkdir "$self_test_root/output" || exit 1
    write_canonical_git_archive "$self_test_root/poisoned.git" \
        "$source_commit" candidate-source \
        "$self_test_root/output/archive.tar" || exit 1
    test "$(/usr/bin/sha256sum "$self_test_root/output/archive.tar" | \
        /usr/bin/cut -d' ' -f1)" = "$clean_hash" || exit 1
    test -z "$(/usr/bin/find "$self_test_root/output" -maxdepth 1 \
        -name '.leopard-canonical-git-archive-output.*' -print)" || exit 1
    if canonical_git_archive_stream \
            "$repo" "$source_commit" bad-prefix >/dev/null 2>&1; then
        exit 1
    fi
    exit 0
)

verify_submodule_gitlink_at_commit()
{
    local source_commit=$1
    local submodule_path=$2
    local submodule_commit=$3
    local observed_record=
    [[ "$source_commit" =~ ^[0-9a-f]{40}$ ]] || return 1
    [[ "$submodule_commit" =~ ^[0-9a-f]{40}$ ]] || return 1
    observed_record="$(/usr/bin/env \
        GIT_CONFIG_GLOBAL=/dev/null \
        GIT_CONFIG_SYSTEM=/dev/null \
        GIT_CONFIG_NOSYSTEM=1 \
        GIT_NO_REPLACE_OBJECTS=1 \
        GIT_OPTIONAL_LOCKS=0 \
        /usr/bin/git -C "$repo" ls-tree "$source_commit" -- \
        "$submodule_path")" || return 1
    test "$observed_record" = \
        "160000 commit $submodule_commit"$'\t'"$submodule_path" || return 1
    return 0
}

verify_v18_complete_core_claim_bindings_preflight()
{
    local verified_core=$1
    local hash_binding=
    local hash_field=
    local hash_path=
    local expected_hash=
    local observed_hash=
    local performance_gate=
    local source_commit=
    local source_tree=
    local sse2neon_commit=
    local candidate_source_archive_hash=
    local baseline_source_archive_hash=
    validate_v18_complete_core_static \
        "$verified_core/manifest.json" || return 1
    for hash_path in \
        build-closure.json \
        build-order-normalizer-reconstruction.json \
        candidate-test-temporal-closure.json \
        result-summary.json \
        wrapper-launch-affinity.json; do
        validate_single_json_object \
            "$verified_core/$hash_path" || return 1
    done
    validate_single_json_array \
        "$verified_core/campaign-command.json" || return 1
    source_commit="$(/usr/bin/jq -er '.source_commit | strings' \
        "$verified_core/manifest.json")" || return 1
    source_tree="$(/usr/bin/jq -er '.source_tree | strings' \
        "$verified_core/manifest.json")" || return 1
    sse2neon_commit="$(/usr/bin/jq -er '.sse2neon_commit | strings' \
        "$verified_core/manifest.json")" || return 1
    verify_source_tree_at_commit "$source_commit" "$source_tree" || return 1
    verify_submodule_gitlink_at_commit \
        "$source_commit" sse2neon "$sse2neon_commit" || return 1
    verify_submodule_gitlink_at_commit \
        "$main_commit" sse2neon "$sse2neon_commit" || return 1
    verify_source_archive_at_commit "$repo" "$source_commit" \
        candidate-source "$verified_core/candidate-source.tar" || return 1
    verify_source_archive_at_commit "$repo" "$main_commit" \
        leopard1-source "$verified_core/leopard1-source.tar" || return 1
    verify_source_archive_at_commit "$repo/sse2neon" "$sse2neon_commit" \
        sse2neon-source "$verified_core/sse2neon-source.tar" || return 1
    candidate_source_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/candidate-source.tar" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    baseline_source_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/leopard1-source.tar" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    /usr/bin/diff -qr \
        "$verified_core/build-closure/live-baseline" \
        "$verified_core/build-closure/replay-baseline" \
        >/dev/null || return 1
    for hash_binding in \
        run_abba.py:"$relative_runner" \
        audit_v17_gfni_main_compare.py:"$relative_auditor" \
        passive_environment_census.py:"$relative_census" \
        leopard2_affinity_supervisor.py:"$relative_supervisor" \
        run-authoritative.sh:"$relative_wrapper" \
        git_capture.py:"$relative_git_capture" \
        balanced_evidence_common.py:"$relative_helper" \
        leopard2_build_provenance.py:"$relative_build_provenance"; do
        verify_core_file_at_source_commit \
            "$verified_core/${hash_binding%%:*}" "$source_commit" \
            "${hash_binding#*:}" || return 1
    done
    for hash_binding in \
        attempt_lineage_sha256:attempt-lineage.json \
        sse2neon_source_archive_sha256:sse2neon-source.tar \
        candidate_binary_sha256_pre_post:build-closure/replay-candidate/bench_leopard2 \
        candidate_archive_sha256_pre_post:build-closure/replay-candidate/libleopard.a \
        baseline_binary_sha256_pre_post:build-closure/replay-baseline/leopard_main_benchmark \
        baseline_archive_sha256_pre_post:build-closure/replay-baseline/libleopard_main_exact.a \
        runner_sha256:run_abba.py \
        auditor_sha256:audit_v17_gfni_main_compare.py \
        supervisor_sha256:leopard2_affinity_supervisor.py \
        passive_census_sha256:passive_environment_census.py \
        wrapper_sha256:run-authoritative.sh \
        build_order_normalizer_sha256:build-order-normalizer.py \
        candidate_timing_normalization_report_sha256:build-closure/candidate-timing-order-normalization/normalized/report.json \
        candidate_test_normalization_report_sha256:build-closure/candidate-test-order-normalization/normalized/report.json \
        campaign_manifest_sha256:campaign/manifest.json \
        campaign_raw_sha256:campaign/raw.json \
        audit_sha256:audit.json \
        controller_affinity_sha256:controller-affinity.json \
        environment_census_pre_sha256:environment-census-pre.json \
        environment_census_post_sha256:environment-census-post.json \
        passive_environment_policy_sha256:passive-environment-policy.json; do
        hash_field=${hash_binding%%:*}
        hash_path=${hash_binding#*:}
        test -f "$verified_core/$hash_path" || return 1
        expected_hash="$(/usr/bin/jq -er --arg field "$hash_field" \
            '.[$field] | strings |
                select(length == 64 and test("^[0-9a-f]{64}$"))' \
            "$verified_core/manifest.json")" || return 1
        observed_hash="$(/usr/bin/sha256sum "$verified_core/$hash_path" | \
            /usr/bin/cut -d' ' -f1)" || return 1
        test "$observed_hash" = "$expected_hash" || return 1
    done
    performance_gate="$(/usr/bin/jq -er '
        [
            .analysis["gf16-high-full"].encode
                .promotion_lower_bound_at_least_1_05,
            .analysis["gf16-high-full"].one_shot_encode
                .promotion_lower_bound_at_least_1_05
        ] |
        if all(.[]; type == "boolean") then all
        else error("invalid gates") end | tostring
    ' "$verified_core/campaign/manifest.json")" || return 1
    /usr/bin/jq -e \
        --slurpfile closure "$verified_core/build-closure.json" \
        --slurpfile campaign "$verified_core/campaign/manifest.json" \
        --slurpfile audit "$verified_core/audit.json" \
        --slurpfile policy "$verified_core/passive-environment-policy.json" \
        --slurpfile temporal \
            "$verified_core/candidate-test-temporal-closure.json" \
        --arg main_commit "$main_commit" --argjson gate "$performance_gate" \
        --arg candidate_source_archive_hash \
            "$candidate_source_archive_hash" \
        --arg baseline_source_archive_hash \
            "$baseline_source_archive_hash" '
        ($closure | length) == 1 and ($campaign | length) == 1 and
        ($audit | length) == 1 and ($policy | length) == 1 and
        ($temporal | length) == 1 and
        $closure[0].schema ==
            "leopard2-v17-gfni-main-build-closure/v3" and
        $campaign[0].schema == "leopard2-main-compare-manifest/v18" and
        $audit[0].schema ==
            "leopard2-main-compare-v18-passive-independent-audit/v1" and
        $policy[0].schema == "leopard2-passive-shared-host-policy/v2" and
        $temporal[0].schema ==
            "leopard2-v17-gfni-main-candidate-test-temporal-closure/v2" and
        .source_commit == $closure[0].candidate.commit and
        .source_commit == $campaign[0].identities.candidate_source.head and
        .source_tree == $closure[0].candidate.tree and
        .source_tree == $campaign[0].identities.candidate_source.tree and
        .baseline_commit == $main_commit and
        .baseline_commit == $closure[0].baseline.commit and
        .baseline_commit == $campaign[0].identities.baseline_source.head and
        .baseline_commit == $audit[0].contract.baseline_main_commit and
        .sse2neon_commit == $closure[0].sse2neon.commit and
        .sse2neon_source_archive_sha256 ==
            $closure[0].sse2neon.source_archive_sha256 and
        $candidate_source_archive_hash ==
            $closure[0].candidate.source_archive_sha256 and
        $baseline_source_archive_hash ==
            $closure[0].baseline.source_archive_sha256 and
        .candidate_binary_sha256_pre_post ==
            $closure[0].candidate.binary_sha256 and
        .candidate_binary_sha256_pre_post ==
            $campaign[0].identities.candidate_executable.sha256 and
        .candidate_archive_sha256_pre_post ==
            $closure[0].candidate.archive_sha256 and
        .candidate_archive_sha256_pre_post ==
            $campaign[0].identities.candidate_archive.sha256 and
        .baseline_binary_sha256_pre_post ==
            $closure[0].baseline.binary_sha256 and
        .baseline_binary_sha256_pre_post ==
            $campaign[0].identities.baseline_executable.sha256 and
        .baseline_archive_sha256_pre_post ==
            $closure[0].baseline.archive_sha256 and
        .baseline_archive_sha256_pre_post ==
            $campaign[0].identities.baseline_archive.sha256 and
        .campaign_raw_sha256 == $campaign[0].raw.sha256 and
        .campaign_raw_sha256 == $audit[0].raw.sha256 and
        .campaign_manifest_sha256 == $audit[0].manifest.sha256 and
        .auditor_sha256 == $audit[0].auditor.sha256 and
        .runner_sha256 == $closure[0].controllers.runner_sha256 and
        .auditor_sha256 == $closure[0].controllers.auditor_sha256 and
        .supervisor_sha256 == $closure[0].controllers.supervisor_sha256 and
        .wrapper_sha256 == $closure[0].controllers.wrapper_sha256 and
        .build_order_normalizer_sha256 ==
            $closure[0].controllers.build_order_normalizer_sha256 and
        .candidate_timing_normalization_report_sha256 ==
            $closure[0].candidate.reproduction.normalization_report_sha256 and
        .candidate_test_normalization_report_sha256 ==
            $closure[0].candidate_tests.reproduction.normalization_report_sha256 and
        .candidate_test_normalization_report_sha256 ==
            $temporal[0].normalization_report_sha256 and
        .performance_gate_passed == $gate and
        .build_reproduction == {
            candidate_timing_qualified_semantic_equal:
                $closure[0].build_reproduction.candidate_timing_qualified_semantic_equal,
            candidate_timing_raw_tree_file_bytes_identical:
                $closure[0].build_reproduction.candidate_timing_raw_tree_file_bytes_identical,
            candidate_tests_qualified_semantic_equal:
                $closure[0].build_reproduction.candidate_tests_qualified_semantic_equal,
            candidate_tests_raw_tree_file_bytes_identical:
                $closure[0].build_reproduction.candidate_tests_raw_tree_file_bytes_identical,
            baseline_raw_byte_identical:
                $closure[0].build_reproduction.baseline_raw_byte_identical
        } and
        .isolation_claim == $audit[0].isolation_claim and
        .windowed_benchmark_cpu_nonidle_excess_jiffies ==
            $audit[0].windowed_observation.windowed
                .benchmark_cpu_nonidle_excess_jiffies and
        .windowed_benchmark_cpu_nonidle_excess_jiffies ==
            $policy[0].windowed_contamination.windowed
                .benchmark_cpu_nonidle_excess_jiffies and
        .windowed_reserved_sibling_nonidle_jiffies ==
            $audit[0].windowed_observation.windowed
                .reserved_sibling_nonidle_jiffies and
        .windowed_reserved_sibling_nonidle_jiffies ==
            $policy[0].windowed_contamination.windowed
                .reserved_sibling_nonidle_jiffies and
        .out_of_window_benchmark_cpu_nonidle_jiffies ==
            $audit[0].windowed_observation.out_of_window
                .benchmark_cpu_nonidle_jiffies and
        .out_of_window_benchmark_cpu_nonidle_jiffies ==
            $policy[0].outer_disclosure.isolation_out_of_window
                .benchmark_cpu_nonidle_jiffies and
        .out_of_window_reserved_sibling_nonidle_jiffies ==
            $audit[0].windowed_observation.out_of_window
                .reserved_sibling_nonidle_jiffies and
        .out_of_window_reserved_sibling_nonidle_jiffies ==
            $policy[0].outer_disclosure.isolation_out_of_window
                .reserved_sibling_nonidle_jiffies
    ' "$verified_core/manifest.json" >/dev/null || return 1
    return 0
}

verify_v18_complete_core_claim_bindings()
{
    local verified_core=$1
    local semantic_root=
    local semantic_envelope=
    local semantic_core=
    local campaign_manifest_hash=
    local audit_hash=
    local postseal_audit_hash=
    local core_manifest_hash=
    local core_sha256sums_hash=
    verify_v18_complete_core_claim_bindings_preflight \
        "$verified_core" || return 1

    # A checksum-valid complete core reached only through a later post-seal
    # failure must satisfy exactly the same semantic verifier as an ordinary
    # NOT_PROMOTED result.  Synthesize only the missing success terminal and
    # post-seal audit views, then delegate against the original canonical core
    # path in a fresh shell.  This avoids maintaining a weaker second
    # implementation of the runner, auditor, census, and build-normalization
    # replay contract.
    semantic_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-complete-core-replay.XXXXXX)" || return 1
    semantic_envelope=${verified_core%/core}
    semantic_core=$verified_core
    test "$semantic_envelope/core" = "$verified_core" || return 1
    /usr/bin/cp --reflink=never "$semantic_core/audit.json" \
        "$semantic_root/postseal-audit.json" || return 1

    campaign_manifest_hash="$(/usr/bin/sha256sum \
        "$semantic_core/campaign/manifest.json" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    audit_hash="$(/usr/bin/sha256sum "$semantic_core/audit.json" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    postseal_audit_hash="$(/usr/bin/sha256sum \
        "$semantic_root/postseal-audit.json" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    core_manifest_hash="$(/usr/bin/sha256sum \
        "$semantic_core/manifest.json" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    core_sha256sums_hash="$(/usr/bin/sha256sum \
        "$semantic_core/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    test "$audit_hash" = "$postseal_audit_hash" || return 1

    /usr/bin/jq -n \
        --slurpfile core "$semantic_core/manifest.json" \
        --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
        --arg audit_sha256 "$audit_hash" \
        --arg postseal_audit_sha256 "$postseal_audit_hash" \
        --arg core_manifest_sha256 "$core_manifest_hash" \
        --arg core_sha256sums_sha256 "$core_sha256sums_hash" '
        ($core | length) == 1 or error("complete core count differs") |
        $core[0] as $core |
        {
            schema:
                "leopard2-v18-gfni-main-passive-not-promoted-envelope/v1",
            status:"complete",
            acquisition_generation:"passive-v2",
            attempt:$core.attempt,
            attempt_budget:$core.attempt_budget,
            attempt_lineage_sha256:$core.attempt_lineage_sha256,
            evidence_class:$core.evidence_class,
            promotion_eligible:false,
            promotion_passed:false,
            causal_performance_claim_eligible:false,
            cpu_pair_exclusive:false,
            windowed_benchmark_cpu_nonidle_excess_jiffies:
                $core.windowed_benchmark_cpu_nonidle_excess_jiffies,
            windowed_reserved_sibling_nonidle_jiffies:
                $core.windowed_reserved_sibling_nonidle_jiffies,
            out_of_window_benchmark_cpu_nonidle_jiffies:
                $core.out_of_window_benchmark_cpu_nonidle_jiffies,
            out_of_window_reserved_sibling_nonidle_jiffies:
                $core.out_of_window_reserved_sibling_nonidle_jiffies,
            performance_gate_passed:$core.performance_gate_passed,
            evidence_valid:true,
            campaign_exit_status:0,
            preseal_audit_passed:true,
            postseal_audit_passed:true,
            postseal_audit_byte_identical:true,
            source_commit:$core.source_commit,
            source_tree:$core.source_tree,
            baseline_commit:$core.baseline_commit,
            campaign_manifest_sha256:$campaign_manifest_sha256,
            audit_sha256:$audit_sha256,
            postseal_audit_sha256:$postseal_audit_sha256,
            core_manifest_sha256:$core_manifest_sha256,
            core_sha256sums_sha256:$core_sha256sums_sha256,
            isolation_claim:$core.isolation_claim
        }
    ' > "$semantic_root/NOT_PROMOTED.json" || return 1
    /usr/bin/chmod 0400 "$semantic_root/NOT_PROMOTED.json" \
        "$semantic_root/postseal-audit.json" || return 1
    /usr/bin/bash "$0" --verify-v18-complete-core-semantics \
        "$semantic_envelope" "$semantic_root/NOT_PROMOTED.json" \
        "$semantic_root/postseal-audit.json" \
        > "$semantic_root/verification.log" \
        2> "$semantic_root/verification-stderr.log" || return 1
    return 0
}

verify_v18_failed_core_claim_bindings()
(
    local verified_core=$1
    local controller_binding=
    local controller_field=
    local controller_path=
    local controller_hash=
    local expected_controller_hash=
    local candidate_hash=
    local candidate_archive_hash=
    local candidate_source_archive_hash=
    local baseline_hash=
    local baseline_archive_hash=
    local baseline_source_archive_hash=
    local runner_hash=
    local evidence_helper_hash=
    local expected_candidate_hash=
    local expected_baseline_hash=
    local expected_lineage_hash=
    local observed_lineage_hash=
    local retained_failure_hash=
    local expected_failure_hash=
    local replay_root=
    local replay_controller_root=
    local replay_runner=
    local replay_status=0
    local reconstructed_normalizer_root=
    local reconstructed_normalizer=
    local recorded_verify_status=
    local source_commit=
    local source_tree=
    local sse2neon_commit=
    local sse2neon_source_archive_hash=
    local cleanup_status=0
    cleanup_v18_failed_core_claim_bindings()
    {
        cleanup_status=$?
        trap - EXIT
        for controller_path in "$replay_root" \
                "$reconstructed_normalizer_root"; do
            if [[ -n "$controller_path" && -e "$controller_path" ]] &&
               ! /usr/bin/find "$controller_path" -depth -delete; then
                cleanup_status=1
            fi
        done
        exit "$cleanup_status"
    }
    trap cleanup_v18_failed_core_claim_bindings EXIT
    validate_v18_failed_core_static \
        "$verified_core/manifest.json" || return 1
    for controller_path in \
        audit.json \
        result-summary.json \
        candidate-test-temporal-closure.json \
        passive-environment-policy.json; do
        require_path_absent "$verified_core/$controller_path" || return 1
    done
    for controller_path in \
        campaign/raw.json \
        campaign/manifest.json \
        campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374 \
        campaign/.leopard2-evidence-pending-ffa5b716b5a57837f7929dfcca4b4dfdeb97210a7fd5a12d2f1978846d6f1743; do
        if [[ -e "$verified_core/$controller_path" ||
              -L "$verified_core/$controller_path" ]]; then
            test -f "$verified_core/$controller_path" || return 1
            test ! -L "$verified_core/$controller_path" || return 1
        fi
    done
    if [[ -e "$verified_core/campaign/manifest.json" ]]; then
        test -f "$verified_core/campaign/raw.json" || return 1
        test ! -L "$verified_core/campaign/raw.json" || return 1
        require_path_absent \
            "$verified_core/campaign/.leopard2-evidence-pending-ffa5b716b5a57837f7929dfcca4b4dfdeb97210a7fd5a12d2f1978846d6f1743" || \
            return 1
    fi
    if [[ -e "$verified_core/campaign/raw.json" ]]; then
        require_path_absent \
            "$verified_core/campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374" || \
            return 1
    fi
    if [[ -e "$verified_core/campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374" ]]; then
        require_path_absent \
            "$verified_core/campaign/raw.json" || return 1
        require_path_absent \
            "$verified_core/campaign/manifest.json" || return 1
        require_path_absent \
            "$verified_core/campaign/.leopard2-evidence-pending-ffa5b716b5a57837f7929dfcca4b4dfdeb97210a7fd5a12d2f1978846d6f1743" || \
            return 1
    fi
    if [[ -e "$verified_core/campaign/.leopard2-evidence-pending-ffa5b716b5a57837f7929dfcca4b4dfdeb97210a7fd5a12d2f1978846d6f1743" ]]; then
        test -f "$verified_core/campaign/raw.json" || return 1
        test ! -L "$verified_core/campaign/raw.json" || return 1
        require_path_absent \
            "$verified_core/campaign/manifest.json" || return 1
        require_path_absent \
            "$verified_core/campaign/.leopard2-evidence-pending-4f44a53c068755081419b7c7bc70e241854a1cc34ff9340529ccad707497f374" || \
            return 1
    fi
    validate_single_json_object \
        "$verified_core/build-closure.json" || return 1
    /usr/bin/jq -e -s '
        length == 1 and
        .[0].schema == "leopard2-v17-gfni-main-build-closure/v3"
    ' "$verified_core/build-closure.json" >/dev/null || return 1
    source_commit="$(/usr/bin/jq -er '.source_commit | strings' \
        "$verified_core/manifest.json")" || return 1
    source_tree="$(/usr/bin/jq -er '.source_tree | strings' \
        "$verified_core/manifest.json")" || return 1
    sse2neon_commit="$(/usr/bin/jq -er \
        '.sse2neon.commit | strings |
            select(length == 40 and test("^[0-9a-f]{40}$"))' \
        "$verified_core/build-closure.json")" || return 1
    verify_source_tree_at_commit "$source_commit" "$source_tree" || return 1
    verify_submodule_gitlink_at_commit \
        "$source_commit" sse2neon "$sse2neon_commit" || return 1
    verify_submodule_gitlink_at_commit \
        "$main_commit" sse2neon "$sse2neon_commit" || return 1
    verify_source_archive_at_commit "$repo" "$source_commit" \
        candidate-source "$verified_core/candidate-source.tar" || return 1
    verify_source_archive_at_commit "$repo" "$main_commit" \
        leopard1-source "$verified_core/leopard1-source.tar" || return 1
    verify_source_archive_at_commit "$repo/sse2neon" "$sse2neon_commit" \
        sse2neon-source "$verified_core/sse2neon-source.tar" || return 1
    for controller_binding in \
        run_abba.py:"$relative_runner" \
        audit_v17_gfni_main_compare.py:"$relative_auditor" \
        passive_environment_census.py:"$relative_census" \
        leopard2_affinity_supervisor.py:"$relative_supervisor" \
        git_capture.py:"$relative_git_capture" \
        balanced_evidence_common.py:"$relative_helper" \
        leopard2_build_provenance.py:"$relative_build_provenance" \
        run-authoritative.sh:"$relative_wrapper"; do
        verify_core_file_at_source_commit \
            "$verified_core/${controller_binding%%:*}" "$source_commit" \
            "${controller_binding#*:}" || return 1
    done
    reconstructed_normalizer_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-failed-normalizer.XXXXXX)" || return 1
    reconstructed_normalizer="$reconstructed_normalizer_root/normalizer.py"
    install_build_order_normalizer "$reconstructed_normalizer" || return 1
    /usr/bin/cmp "$verified_core/build-order-normalizer.py" \
        "$reconstructed_normalizer" || return 1
    for controller_binding in \
        runner_sha256:run_abba.py \
        auditor_sha256:audit_v17_gfni_main_compare.py \
        supervisor_sha256:leopard2_affinity_supervisor.py \
        wrapper_sha256:run-authoritative.sh \
        build_order_normalizer_sha256:build-order-normalizer.py; do
        controller_field=${controller_binding%%:*}
        controller_path=${controller_binding#*:}
        controller_hash="$(/usr/bin/sha256sum \
            "$verified_core/$controller_path" | \
            /usr/bin/cut -d' ' -f1)" || return 1
        expected_controller_hash="$(/usr/bin/jq -er \
            --arg field "$controller_field" \
            '.controllers[$field] | strings |
                select(length == 64 and test("^[0-9a-f]{64}$"))' \
            "$verified_core/build-closure.json")" || return 1
        test "$controller_hash" = "$expected_controller_hash" || return 1
    done
    candidate_hash="$(/usr/bin/sha256sum \
        "$verified_core/build-closure/replay-candidate/bench_leopard2" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    baseline_hash="$(/usr/bin/sha256sum \
        "$verified_core/build-closure/replay-baseline/leopard_main_benchmark" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    candidate_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/build-closure/replay-candidate/libleopard.a" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    candidate_source_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/candidate-source.tar" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    baseline_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/build-closure/replay-baseline/libleopard_main_exact.a" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    baseline_source_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/leopard1-source.tar" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    sse2neon_source_archive_hash="$(/usr/bin/sha256sum \
        "$verified_core/sse2neon-source.tar" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    runner_hash="$(/usr/bin/sha256sum "$verified_core/run_abba.py" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    evidence_helper_hash="$(/usr/bin/sha256sum \
        "$verified_core/balanced_evidence_common.py" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    expected_candidate_hash="$(/usr/bin/jq -er \
        '.candidate_binary_sha256 | strings' \
        "$verified_core/manifest.json")" || return 1
    expected_baseline_hash="$(/usr/bin/jq -er \
        '.baseline_binary_sha256 | strings' \
        "$verified_core/manifest.json")" || return 1
    test "$candidate_hash" = "$expected_candidate_hash" || return 1
    test "$baseline_hash" = "$expected_baseline_hash" || return 1
    observed_lineage_hash="$(/usr/bin/sha256sum \
        "$verified_core/attempt-lineage.json" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    expected_lineage_hash="$(/usr/bin/jq -er \
        '.attempt_lineage_sha256 | strings' \
        "$verified_core/manifest.json")" || return 1
    test "$observed_lineage_hash" = "$expected_lineage_hash" || return 1
    /usr/bin/jq -e --slurpfile closure "$verified_core/build-closure.json" \
        --arg main_commit "$main_commit" \
        --arg candidate_archive_hash "$candidate_archive_hash" \
        --arg candidate_source_archive_hash \
            "$candidate_source_archive_hash" \
        --arg baseline_archive_hash "$baseline_archive_hash" \
        --arg baseline_source_archive_hash \
            "$baseline_source_archive_hash" \
        --arg sse2neon_commit "$sse2neon_commit" \
        --arg sse2neon_source_archive_hash \
            "$sse2neon_source_archive_hash" \
        --arg runner_hash "$runner_hash" '
        ($closure | length) == 1 and
        $closure[0].schema ==
            "leopard2-v17-gfni-main-build-closure/v3" and
        .source_commit == $closure[0].candidate.commit and
        .source_tree == $closure[0].candidate.tree and
        .baseline_commit == $main_commit and
        .baseline_commit == $closure[0].baseline.commit and
        .candidate_binary_sha256 == $closure[0].candidate.binary_sha256 and
        .baseline_binary_sha256 == $closure[0].baseline.binary_sha256 and
        $candidate_archive_hash == $closure[0].candidate.archive_sha256 and
        $candidate_source_archive_hash ==
            $closure[0].candidate.source_archive_sha256 and
        $baseline_archive_hash == $closure[0].baseline.archive_sha256 and
        $baseline_source_archive_hash ==
            $closure[0].baseline.source_archive_sha256 and
        $sse2neon_commit == $closure[0].sse2neon.commit and
        $closure[0].sse2neon.archive_prefix == "sse2neon-source/" and
        $closure[0].sse2neon.reproduced_from_candidate_and_baseline_clones ==
            true and
        $sse2neon_source_archive_hash ==
            $closure[0].sse2neon.source_archive_sha256 and
        $runner_hash == $closure[0].controllers.runner_sha256
    ' "$verified_core/manifest.json" >/dev/null || return 1
    if [[ -f "$verified_core/campaign/failure.json" ||
          -f "$verified_core/campaign/raw.json" ||
          -f "$verified_core/campaign/manifest.json" ]]; then
        replay_root="$(/usr/bin/mktemp -d \
            /tmp/leopard-v18-failed-claim-replay.XXXXXX)" || return 1
        reconstruct_owner_only_campaign_tree "$verified_core/campaign" \
            "$replay_root/campaign" || return 1
        replay_controller_root="$replay_root/controller"
        /usr/bin/mkdir -p \
            "$replay_controller_root/experiments/leopard2/main_compare" \
            "$replay_controller_root/experiments/leopard2/decoder_dispatch" \
            "$replay_controller_root/tools" || return 1
        /usr/bin/cp --reflink=never "$verified_core/run_abba.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/run_abba.py" || \
            return 1
        /usr/bin/cp --reflink=never "$verified_core/git_capture.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/git_capture.py" || \
            return 1
        /usr/bin/cp --reflink=never \
            "$verified_core/balanced_evidence_common.py" \
            "$replay_controller_root/experiments/leopard2/decoder_dispatch/balanced_evidence_common.py" || \
            return 1
        /usr/bin/cp --reflink=never \
            "$verified_core/leopard2_build_provenance.py" \
            "$replay_controller_root/tools/leopard2_build_provenance.py" || \
            return 1
        /usr/bin/find "$replay_controller_root" -type d \
            -exec /usr/bin/chmod 0700 {} + || return 1
        /usr/bin/find "$replay_controller_root" -type f \
            -exec /usr/bin/chmod 0600 {} + || return 1
        replay_runner=
        replay_runner+="$replay_controller_root/experiments/leopard2/"
        replay_runner+=main_compare/run_abba.py
    fi
    if [[ -f "$verified_core/campaign/failure.json" ]]; then
        retained_failure_hash="$(/usr/bin/sha256sum \
            "$verified_core/campaign/failure.json" | \
            /usr/bin/cut -d' ' -f1)" || return 1
        expected_failure_hash="$(/usr/bin/jq -er \
            '.failure_sha256 | strings |
                select(length == 64 and test("^[0-9a-f]{64}$"))' \
            "$verified_core/manifest.json")" || return 1
        test "$retained_failure_hash" = "$expected_failure_hash" || return 1
        recorded_verify_status="$(/usr/bin/jq -er \
            '.failure_verify_status | numbers' \
            "$verified_core/manifest.json")" || return 1
        if [[ "$recorded_verify_status" -eq 0 ]]; then
            /usr/bin/jq -e --slurpfile failure \
                "$verified_core/campaign/failure.json" \
                --arg runner_hash "$runner_hash" \
                --arg evidence_helper_hash "$evidence_helper_hash" \
                --arg candidate_archive_hash "$candidate_archive_hash" \
                --arg baseline_archive_hash "$baseline_archive_hash" '
            ($failure | length) == 1 and
            $failure[0].schema ==
                "leopard2-main-compare-failure/v18" and
            $failure[0].evidence_contract ==
                "leopard2-main-compare-failure-evidence-contract/v18" and
            .source_commit ==
                $failure[0].input_specification.candidate_commit and
            .cpu == $failure[0].campaign.benchmark_cpu and
            .sibling == $failure[0].campaign.reserved_sibling and
            (if ($failure[0].identities_initial | type) == "object" then
                .source_commit ==
                    $failure[0].identities_initial.candidate_source.head and
                .source_tree ==
                    $failure[0].identities_initial.candidate_source.tree and
                .baseline_commit ==
                    $failure[0].identities_initial.baseline_source.head and
                .candidate_binary_sha256 ==
                    $failure[0].identities_initial.candidate_executable.sha256 and
                .baseline_binary_sha256 ==
                    $failure[0].identities_initial.baseline_executable.sha256 and
                $candidate_archive_hash ==
                    $failure[0].identities_initial.candidate_archive.sha256 and
                $baseline_archive_hash ==
                    $failure[0].identities_initial.baseline_archive.sha256 and
                $runner_hash ==
                    $failure[0].identities_initial.runner.sha256 and
                $evidence_helper_hash ==
                    $failure[0].identities_initial.evidence_helper.sha256
            else
                $failure[0].identities_initial == null
            end)
            ' "$verified_core/manifest.json" >/dev/null || return 1
        fi
        if /usr/bin/python3 -I -S -B \
                "$replay_runner" \
                verify-failure \
                --failure "$replay_root/campaign/failure.json" \
                >/dev/null 2>&1; then
            replay_status=0
        else
            replay_status=$?
        fi
        test "$replay_status" -eq "$recorded_verify_status" || return 1
    else
        require_path_absent \
            "$verified_core/campaign/failure.json" || return 1
        test "$(/usr/bin/jq -er \
            '.failure_sha256 == null and .failure_verify_status == 1 and
             .failure_verified == false' \
            "$verified_core/manifest.json")" = true || return 1
    fi
    if [[ -f "$verified_core/campaign/raw.json" ]]; then
        /usr/bin/python3 -I -S -B -c '
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1])
import run_abba

raw_path = Path(sys.argv[2])
directory = run_abba.EvidenceDirectory.open_existing(raw_path.parent)
try:
    _, raw_bytes = directory.snapshot(
        raw_path.name, run_abba.MAX_IDENTITY_FILE_BYTES)
    raw = run_abba.strict_json_loads(raw_bytes, "retained v18 raw diagnostic")
    run_abba.require(
        isinstance(raw, dict) and raw.get("schema") == run_abba.RAW_SCHEMA_V18,
        "failed v18 campaign retained a non-v18 raw diagnostic")
    run_abba.validate_raw(
        raw, raw_path.parent, check_files=True,
        check_current_inputs=False, evidence_directory=directory)
finally:
    directory.close()
' "${replay_runner%/*}" "$replay_root/campaign/raw.json" \
            >/dev/null 2>&1 || return 1
    fi
    if [[ -f "$verified_core/campaign/manifest.json" ]]; then
        /usr/bin/python3 -I -S -B "$replay_runner" verify \
            --manifest "$replay_root/campaign/manifest.json" \
            --no-current-input-check >/dev/null 2>&1 || return 1
    fi
    return 0
)

write_v18_attempt_lineage()
{
    local lineage_path="$lane/attempt-lineage.json"
    test -n "$attempt_lineage_json" || return 1
    if [[ ! -e "$lineage_path" ]]; then
        /usr/bin/printf '%s\n' "$attempt_lineage_json" > "$lineage_path" || \
            return 1
    fi
    test -f "$lineage_path" || return 1
    test ! -L "$lineage_path" || return 1
    test "$(/usr/bin/cat "$lineage_path")" = "$attempt_lineage_json" || \
        return 1
    attempt_lineage_sha256="$(/usr/bin/sha256sum "$lineage_path" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    [[ "$attempt_lineage_sha256" =~ ^[0-9a-f]{64}$ ]] || return 1
    return 0
}

validate_v18_attempt_lineage_shape()
{
    local lineage=$1
    local current_status=$2
    validate_single_json_object "$lineage" || return 1
    validate_single_json_object "$current_status" || return 1
    /usr/bin/jq -e --arg repo "$repo" --slurpfile status "$current_status" '
        . as $lineage |
        keys == (["acquisition_generation","attempt","attempt_budget",
            "prior_attempts","schema","source_commit","source_tree"] | sort) and
        .schema == "leopard2-v18-gfni-main-attempt-lineage/v1" and
        .acquisition_generation == "passive-v2" and
        (.attempt | type) == "number" and .attempt == (.attempt | floor) and
        .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
        (.source_commit | type) == "string" and
        (.source_commit | length) == 40 and
        (.source_commit | test("^[0-9a-f]{40}$")) and
        (.source_tree | type) == "string" and
        (.source_tree | length) == 40 and
        (.source_tree | test("^[0-9a-f]{40}$")) and
        .attempt == $status[0].attempt and
        .attempt_budget == $status[0].attempt_budget and
        .source_commit == $status[0].source_commit and
        .source_tree == $status[0].source_tree and
        (.prior_attempts | type) == "array" and
        (.prior_attempts | length) == (.attempt - 1) and
        [.prior_attempts[].attempt] == [range(1; $lineage.attempt)] and
        all(.prior_attempts[];
            keys == (["attempt","envelope","envelope_sha256sums_sha256",
                "terminal","terminal_schema"] | sort) and
            (.attempt | type) == "number" and .attempt == (.attempt | floor) and
            (.envelope | type) == "string" and
            .terminal == "FAILED.json" and
            .terminal_schema == "leopard2-v18-gfni-main-failed-envelope/v1" and
            (.envelope_sha256sums_sha256 | type) == "string" and
            (.envelope_sha256sums_sha256 | length) == 64 and
            (.envelope_sha256sums_sha256 | test("^[0-9a-f]{64}$"))) and
        all(range(0; $lineage.attempt - 1);
            . as $index |
            $lineage.prior_attempts[$index].envelope ==
                ($repo + "/.research/leopard-79h/" +
                 $lineage.source_commit[0:7] +
                 "-v18-passive-main-a" + (($index + 1) | tostring)))
    ' "$lineage" >/dev/null
}

validate_v18_attempt_lineage()
{
    local current_envelope=$1
    local current_core=$2
    local current_status=$3
    local lineage="$current_core/attempt-lineage.json"
    local current_attempt=
    local current_commit=
    local current_tree=
    local prior_index=
    local prior_path=
    local expected_prior_path=
    local prior_sha_before=
    local prior_sha_after=
    local recorded_prior_sha=
    local observed_lineage_sha=
    local recorded_lineage_sha=
    test -f "$lineage" || return 1
    test ! -L "$lineage" || return 1
    observed_lineage_sha="$(/usr/bin/sha256sum "$lineage" | \
        /usr/bin/cut -d' ' -f1)" || return 1
    recorded_lineage_sha="$(/usr/bin/jq -er \
        '.attempt_lineage_sha256 | strings |
            select(length == 64 and test("^[0-9a-f]{64}$"))' \
        "$current_status")" || return 1
    test "$observed_lineage_sha" = "$recorded_lineage_sha" || return 1
    validate_v18_attempt_lineage_shape "$lineage" "$current_status" || return 1
    current_attempt="$(/usr/bin/jq -er '.attempt | numbers' "$lineage")" || \
        return 1
    current_commit="$(/usr/bin/jq -er '.source_commit | strings' "$lineage")" || \
        return 1
    current_tree="$(/usr/bin/jq -er '.source_tree | strings' "$lineage")" || \
        return 1
    for ((prior_index = 1; prior_index < current_attempt; ++prior_index)); do
        prior_path="$(/usr/bin/jq -er --argjson index "$prior_index" \
            '.prior_attempts[$index - 1].envelope | strings' "$lineage")" || \
            return 1
        expected_prior_path="$repo/.research/leopard-79h/${current_commit:0:7}-v18-passive-main-a${prior_index}"
        test "$prior_path" = "$expected_prior_path" || return 1
        test "$prior_path" != "$current_envelope" || return 1
        test -f "$prior_path/FAILED.json" || return 1
        test ! -e "$prior_path/COMPLETED.json" || return 1
        test ! -e "$prior_path/NOT_PROMOTED.json" || return 1
        recorded_prior_sha="$(/usr/bin/jq -er --argjson index "$prior_index" \
            '.prior_attempts[$index - 1].envelope_sha256sums_sha256 | strings' \
            "$lineage")" || return 1
        prior_sha_before="$(/usr/bin/sha256sum "$prior_path/SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)" || return 1
        test "$prior_sha_before" = "$recorded_prior_sha" || return 1
        /usr/bin/cmp "$prior_path/core/run-authoritative.sh" \
            "$current_core/run-authoritative.sh" || return 1
        /usr/bin/cmp "$current_core/run-authoritative.sh" "$0" || return 1
        /usr/bin/jq -e \
            --argjson attempt "$prior_index" \
            --arg commit "$current_commit" --arg tree "$current_tree" '
            .schema == "leopard2-v18-gfni-main-failed-envelope/v1" and
            .acquisition_generation == "passive-v2" and
            .attempt == $attempt and .attempt_budget == 3 and
            .source_commit == $commit and .source_tree == $tree
        ' "$prior_path/FAILED.json" >/dev/null || return 1
        /usr/bin/bash "$0" \
            --verify "$prior_path" >/dev/null || return 1
        if v18_envelope_has_observational_output "$prior_path"; then
            return 1
        fi
        prior_sha_after="$(/usr/bin/sha256sum "$prior_path/SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)" || return 1
        test "$prior_sha_after" = "$prior_sha_before" || return 1
    done
    return 0
}

v18_attempt_lineage_self_test()
{
    local failure_schema=$1
    local self_test_root=
    local self_test_core=
    local lineage_hash=
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-attempt-lineage-self-test.XXXXXX)"
    self_test_core="$self_test_root/core"
    /usr/bin/mkdir -m 0700 "$self_test_core"
    /usr/bin/cp --reflink=never "$repo/$relative_wrapper" \
        "$self_test_core/run-authoritative.sh"
    /usr/bin/jq -cS -n \
        '{schema:"leopard2-v18-gfni-main-attempt-lineage/v1",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",prior_attempts:[]}' \
        > "$self_test_core/attempt-lineage.json"
    lineage_hash="$(/usr/bin/sha256sum \
        "$self_test_core/attempt-lineage.json" | /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n --arg schema "$failure_schema" \
        --arg lineage_hash "$lineage_hash" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:$lineage_hash,promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$self_test_root/FAILED.json"
    validate_v18_attempt_lineage \
        "$self_test_root" "$self_test_core" "$self_test_root/FAILED.json"
    /usr/bin/jq -n --arg schema "$failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("e" * 64),promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$self_test_root/BAD-FAILED.json"
    if validate_v18_attempt_lineage \
            "$self_test_root" "$self_test_core" \
            "$self_test_root/BAD-FAILED.json"; then
        return 1
    fi
    /usr/bin/jq -cS -n \
        '{schema:"leopard2-v18-gfni-main-attempt-lineage/v1",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,source_commit:(("a" * 40) + "\n"),source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",prior_attempts:[]}' \
        > "$self_test_core/attempt-lineage.json"
    lineage_hash="$(/usr/bin/sha256sum \
        "$self_test_core/attempt-lineage.json" | /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n --arg schema "$failure_schema" \
        --arg lineage_hash "$lineage_hash" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:$lineage_hash,promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:(("a" * 40) + "\n"),source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$self_test_root/BAD-NEWLINE-FAILED.json"
    if validate_v18_attempt_lineage \
            "$self_test_root" "$self_test_core" \
            "$self_test_root/BAD-NEWLINE-FAILED.json"; then
        return 1
    fi
    /usr/bin/jq -cS -n \
        --arg envelope \
            "$repo/.research/leopard-79h/aaaaaaa-v18-passive-main-a1" \
        '{schema:"leopard2-v18-gfni-main-attempt-lineage/v1",acquisition_generation:"passive-v2",attempt:2,attempt_budget:3,source_commit:("a" * 40),source_tree:("b" * 40),prior_attempts:[{attempt:1,envelope:($envelope + "\n"),terminal:"FAILED.json",terminal_schema:"leopard2-v18-gfni-main-failed-envelope/v1",envelope_sha256sums_sha256:("c" * 64)}]}' \
        > "$self_test_core/BAD-PATH-LINEAGE.json"
    /usr/bin/jq -n --arg schema "$failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:2,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:("a" * 40),source_tree:("b" * 40),core_sha256sums_sha256:("c" * 64)}' \
        > "$self_test_root/BAD-PATH-FAILED.json"
    if validate_v18_attempt_lineage_shape \
            "$self_test_core/BAD-PATH-LINEAGE.json" \
            "$self_test_root/BAD-PATH-FAILED.json"; then
        return 1
    fi
    /usr/bin/rm -rf -- "$self_test_root"
}

tools_to_hash=(
    /usr/bin/bash /usr/bin/env /usr/bin/c++ /usr/bin/cc /usr/bin/ld
    /usr/bin/ar /usr/bin/ranlib /usr/bin/cmake /usr/bin/gmake /usr/bin/git
    /usr/bin/objdump /usr/bin/readelf /usr/bin/nm /usr/bin/python3
    /usr/bin/taskset /usr/bin/ldd /usr/bin/flock /usr/bin/jq
    /usr/bin/sha256sum /usr/bin/find /usr/bin/grep /usr/bin/cmp /usr/bin/cp
    /usr/bin/stat /usr/bin/chmod /usr/bin/sort /usr/bin/xargs /usr/bin/cut
    /usr/bin/tee /usr/bin/ctest /usr/bin/cat /usr/bin/dirname /usr/bin/lscpu
    /usr/bin/mkdir /usr/bin/mktemp /usr/bin/printf /usr/bin/readlink
    /usr/bin/uname /usr/bin/mv /usr/bin/diff /usr/bin/rm
)

require_empty_output()
{
    local observed_output
    observed_output="$("$@")" || return 1
    test -z "$observed_output"
}

passive_contract_self_test()
{
    local observed_housekeeping=
    local observed_timeout_argument=
    local expected_core_schema=leopard2-v18-gfni-main-
    local expected_terminal_schema=leopard2-v18-gfni-main-
    local expected_result_schema=leopard2-v18-gfni-main-
    local expected_failure_schema=leopard2-v18-gfni-main-
    local expected_wrapper_failure_schema=leopard2-v18-gfni-main-
    local absence_self_test_root=
    local validation_self_test_root=
    local self_test_commit=
    local self_test_tree=
    expected_core_schema+=passive-core-manifest/v1
    expected_terminal_schema+=passive-not-promoted-envelope/v1
    expected_result_schema+=passive-result-summary/v1
    expected_failure_schema+=failed-envelope/v1
    expected_wrapper_failure_schema+=passive-wrapper-failure/v1
    self_test_commit="$(/usr/bin/git -C "$repo" rev-parse HEAD)" || return 1
    self_test_tree="$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" || \
        return 1
    verify_source_tree_at_commit "$self_test_commit" "$self_test_tree" || \
        return 1
    canonical_git_archive_self_test "$self_test_commit" || return 1
    if verify_source_tree_at_commit "$self_test_tree" "$self_test_tree" \
            >/dev/null 2>&1; then
        return 1
    fi
    if verify_source_tree_at_commit "$self_test_commit" \
            0000000000000000000000000000000000000000 \
            >/dev/null 2>&1; then
        return 1
    fi
    observed_housekeeping="$(/usr/bin/printf '%s\n' \
        '{"allowed_cpus":[0,1,52,116]}' | \
        /usr/bin/jq -er "$passive_housekeeping_jq")" || return 1
    test "$observed_housekeeping" = 0,1
    if /usr/bin/printf '%s\n' '{"allowed_cpus":[52,116]}' | \
            /usr/bin/jq -er "$passive_housekeeping_jq" >/dev/null 2>&1; then
        return 1
    fi
    observed_timeout_argument="$(/usr/bin/printf '%s\n' \
        '{"campaign":{"timeout_seconds":600.0}}' | \
        /usr/bin/jq -er "$passive_timeout_argument_jq")" || return 1
    test "$observed_timeout_argument" = 600
    if /usr/bin/printf '%s\n' \
            '{"campaign":{"timeout_seconds":600.5}}' | \
            /usr/bin/jq -er "$passive_timeout_argument_jq" \
                >/dev/null 2>&1; then
        return 1
    fi
    validate_attempt_contract 1 3
    validate_attempt_contract 2 3
    validate_attempt_contract 3 3
    local bad_attempt=
    local bad_budget=
    for bad_attempt in 0 4 01 -1 true ''; do
        if validate_attempt_contract "$bad_attempt" 3; then
            return 1
        fi
    done
    for bad_budget in 0 2 4 03 true ''; do
        if validate_attempt_contract 1 "$bad_budget"; then
            return 1
        fi
    done
    absence_self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-absence-self-test.XXXXXX)" || return 1
    /usr/bin/python3 -I -S -B -c \
        'import os, sys; os.symlink("missing", sys.argv[1])' \
        "$absence_self_test_root/dangling" || return 1
    if require_path_absent "$absence_self_test_root/dangling"; then
        return 1
    fi
    /usr/bin/find "$absence_self_test_root" -depth -delete || return 1
    validate_single_json_object \
        <(/usr/bin/printf '%s\n' '{"schema":"fixture"}') || return 1
    if validate_single_json_object \
            <(/usr/bin/printf '%s\n%s\n' '{}' \
                '{"schema":"fixture"}'); then
        return 1
    fi
    if validate_single_json_object \
            <(/usr/bin/printf '%s\n' \
                '{"schema":"first","schema":"second"}'); then
        return 1
    fi
    if validate_single_json_object \
            <(/usr/bin/printf '%s\n' '{"value":NaN}'); then
        return 1
    fi
    if validate_single_json_object \
            <(/usr/bin/printf '%s\n' '{"value":1e309}'); then
        return 1
    fi
    validate_single_json_array \
        <(/usr/bin/printf '%s\n' '[]') || return 1
    if validate_single_json_array \
            <(/usr/bin/printf '%s\n' '{}'); then
        return 1
    fi
    validate_exact_json_schema \
        <(/usr/bin/jq -n --arg schema "$expected_terminal_schema" \
            '{schema:$schema}') "$expected_terminal_schema"
    if validate_exact_json_schema \
            <(/usr/bin/jq -n \
                --arg schema "${expected_terminal_schema}"$'\n' \
                '{schema:$schema}') "$expected_terminal_schema"; then
        return 1
    fi
    validation_self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v18-validation-self-test.XXXXXX)" || return 1
    /usr/bin/jq -n --arg schema "$expected_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,source_commit:("a" * 40),source_tree:("b" * 40),core_sha256sums_sha256:("c" * 64)}' \
        > "$validation_self_test_root/terminal.json" || return 1
    validate_v18_terminal_static \
        "$validation_self_test_root/terminal.json" failed-v2 || return 1
    /usr/bin/jq -n --arg schema "$expected_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,source_commit:("a" * 40),source_tree:("b" * 40),core_sha256sums_sha256:(("c" * 64) + "\n")}' \
        > "$validation_self_test_root/terminal.json" || return 1
    if validate_v18_terminal_static \
            "$validation_self_test_root/terminal.json" failed-v2; then
        return 1
    fi
    /usr/bin/jq -n --arg schema "$expected_terminal_schema" \
        '{schema:$schema,status:"complete",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:0,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",campaign_manifest_sha256:("c" * 64),audit_sha256:("e" * 64),postseal_audit_sha256:("f" * 64),core_manifest_sha256:("1" * 64),core_sha256sums_sha256:("2" * 64)}' \
        > "$validation_self_test_root/terminal.json" || return 1
    validate_v18_terminal_static \
        "$validation_self_test_root/terminal.json" passive-v2 || return 1
    /usr/bin/jq -n --arg schema "$expected_terminal_schema" \
        '{schema:$schema,status:"complete",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:0,source_commit:("a" * 40),source_tree:("b" * 40),baseline_commit:"6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",campaign_manifest_sha256:(("c" * 64) + "\n"),audit_sha256:("e" * 64),postseal_audit_sha256:("f" * 64),core_manifest_sha256:("1" * 64),core_sha256sums_sha256:("2" * 64)}' \
        > "$validation_self_test_root/terminal.json" || return 1
    if validate_v18_terminal_static \
            "$validation_self_test_root/terminal.json" passive-v2; then
        return 1
    fi
    /usr/bin/jq -n --arg schema "$expected_wrapper_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,stage:"fixture",exit_status:7,source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"}' \
        > "$validation_self_test_root/wrapper.json" || return 1
    /usr/bin/jq -n --arg schema "$expected_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$validation_self_test_root/status.json" || return 1
    validate_v18_wrapper_failure_binding \
        "$validation_self_test_root/wrapper.json" \
        "$validation_self_test_root/status.json" || return 1
    /usr/bin/jq -n --arg schema "$expected_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:2,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$validation_self_test_root/status.json" || return 1
    if validate_v18_wrapper_failure_binding \
            "$validation_self_test_root/wrapper.json" \
            "$validation_self_test_root/status.json"; then
        return 1
    fi
    /usr/bin/jq -n --arg schema "$expected_wrapper_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,stage:"fixture",exit_status:7,source_commit:(("a" * 40) + "\n"),source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"}' \
        > "$validation_self_test_root/wrapper.json" || return 1
    /usr/bin/jq -n --arg schema "$expected_failure_schema" \
        '{schema:$schema,status:"failed",acquisition_generation:"passive-v2",attempt:1,attempt_budget:3,attempt_lineage_sha256:("d" * 64),promotion_passed:false,campaign_exit_status:7,stage:"fixture",source_commit:(("a" * 40) + "\n"),source_tree:"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",core_sha256sums_sha256:("c" * 64)}' \
        > "$validation_self_test_root/status.json" || return 1
    if validate_v18_wrapper_failure_binding \
            "$validation_self_test_root/wrapper.json" \
            "$validation_self_test_root/status.json"; then
        return 1
    fi
    /usr/bin/find "$validation_self_test_root" -depth -delete || return 1
    v18_attempt_lineage_self_test "$expected_failure_schema"
    v18_failed_core_static_self_test
    v18_observational_output_self_test
    /usr/bin/grep -Fq "$expected_core_schema" \
        "$repo/$relative_wrapper" || return 1
    /usr/bin/grep -Fq "$expected_terminal_schema" \
        "$repo/$relative_wrapper" || return 1
    /usr/bin/grep -Fq "$expected_result_schema" \
        "$repo/$relative_wrapper" || return 1
    /usr/bin/grep -Fq "$expected_failure_schema" \
        "$repo/$relative_wrapper" || return 1
    /usr/bin/printf 'v18 passive wrapper contract self-test passed\n'
}

install_build_order_normalizer()
{
    local normalizer_output=$1
    test ! -e "$normalizer_output" || return 1
    /usr/bin/tee "$normalizer_output" >/dev/null <<'PY'
#!/usr/bin/python3
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import stat
import tempfile

SCHEMA = "leopard2-v17-candidate-build-order-normalization/v3"
EXCEPTIONS = ("compile_commands.json", "CMakeFiles/Makefile2")
COMPILE_KEYS = {"command", "directory", "file", "output"}
TEMPLATE_BUILD = b"${BUILD}"
TEMPLATE_SOURCE = b"${SOURCE}"
TEMPLATE_COMMIT = b"${COMMIT}"
MAKE_RULE = re.compile(
    r"^(CMakeFiles/[A-Za-z0-9_.+@-]+\.dir/all): "
    r"(CMakeFiles/[A-Za-z0-9_.+@-]+\.dir/all)\n$")
MAKE_RULE_LIKE = re.compile(
    r"^(CMakeFiles/[A-Za-z0-9_.+@-]+\.dir/all)[ \t]*:+[ \t]*"
    r"(CMakeFiles/[A-Za-z0-9_.+@-]+\.dir/all)[ \t]*\n$")
ORDER_SENSITIVE_AUTOMATIC = re.compile(
    r"\$(?:[<^+?|]|\([<^+?|](?:[DF])?\)|\{[<^+?|](?:[DF])?\})")
COMMON_MAKE_BLOCKS = (
    (
        "CMakeFiles/leopard.dir/all",
        (
            "CMakeFiles/leopard2_backend_avx2.dir/all",
            "CMakeFiles/leopard2_backend_avx2_t16_b64.dir/all",
            "CMakeFiles/leopard2_backend_avx2_t16_k66.dir/all",
            "CMakeFiles/leopard2_backend_avx2_t2_k4.dir/all",
            "CMakeFiles/leopard2_backend_avx2_t32_b256.dir/all",
            "CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/all",
            "CMakeFiles/leopard2_backend_gfni.dir/all",
            "CMakeFiles/leopard2_backend_ssse3.dir/all",
            "CMakeFiles/leopard2_low_p32_b64_avx2.dir/all",
        ),
    ),
    (
        "CMakeFiles/bench_leopard2.dir/all",
        (
            "CMakeFiles/leopard.dir/all",
            "CMakeFiles/leopard2_benchmark_source_attestation_refresh.dir/all",
        ),
    ),
    (
        "CMakeFiles/bench_leopard2_prevalidated_batch.dir/all",
        (
            "CMakeFiles/leopard.dir/all",
            "CMakeFiles/leopard2_benchmark_source_attestation_refresh.dir/all",
        ),
    ),
    (
        "CMakeFiles/bench_leopard2_sparse_encode.dir/all",
        (
            "CMakeFiles/leopard.dir/all",
            "CMakeFiles/leopard2_sparse_encode_benchmark_object.dir/all",
            "CMakeFiles/leopard2_sparse_encode_oracle_object.dir/all",
        ),
    ),
)
TEST_HOOKS_BLOCK = (
    "CMakeFiles/leopard_test_hooks.dir/all",
    (
        "CMakeFiles/leopard2_backend_avx2_t16_b64.dir/all",
        "CMakeFiles/leopard2_backend_avx2_t16_k66.dir/all",
        "CMakeFiles/leopard2_backend_avx2_t2_k4.dir/all",
        "CMakeFiles/leopard2_backend_avx2_t8_k8_b1024.dir/all",
        "CMakeFiles/leopard2_backend_avx2_test_hooks.dir/all",
        "CMakeFiles/leopard2_backend_gfni_test_hooks.dir/all",
        "CMakeFiles/leopard2_backend_ssse3_test_hooks.dir/all",
        "CMakeFiles/leopard2_low_p32_b64_avx2.dir/all",
    ),
)
DIRECT_ENCODE_BLOCK = (
    "CMakeFiles/bench_leopard2_direct_encode.dir/all",
    (
        "CMakeFiles/leopard2_benchmark_source_attestation_refresh.dir/all",
        "CMakeFiles/leopard_test_hooks.dir/all",
    ),
)
PROFILES = {
    "candidate-timing": {
        "compile_entry_count": 30,
        "compiler_counts": {"/usr/bin/c++": 30},
        "compile_output_list_sha256":
            "a39c7e87cc0a506148faaa4023b11f96c2e26b22c103584a62b62ae91f2387b7",
        "compile_template_sha256":
            "6ad3a97ffe31f4eaf5c8f2ff5459dd310635e6e430181c460d9803f11cc5fb28",
        "make_blocks": COMMON_MAKE_BLOCKS,
        "make_normalized_targets": frozenset(
            target for target, _prerequisites in COMMON_MAKE_BLOCKS),
        "make_template_sha256":
            "5ab77166d931546b1cb998cfdf4eb5c7dd1850994b2eef65d8c5e627a3899a62",
        "node_census_count": 131,
        "node_census_sha256":
            "d3194ea1135fbdc0067cd3af72693cbf27d1e28a4734804b52bf5be1a3f37aca",
        "profile_sha256":
            "d2a07561d2e5919051c2b7d7464f6aa844752e42e92099190feeea004f2716b3",
        "single_c_record": None,
        "template_origin":
            "rederived from sealed v1 candidate timing closure",
    },
    "candidate-tests": {
        "compile_entry_count": 170,
        "compiler_counts": {"/usr/bin/c++": 169, "/usr/bin/cc": 1},
        "compile_output_list_sha256":
            "8aeea7120bbbe5e3f5e71059b29007b4fb3e01b153c55f0aa874c45de9e48574",
        "compile_template_sha256":
            "dda169bcd458db8b17c09c1e8873c9342731a4c5639874849ad87d3c1003badb",
        "make_blocks":
            COMMON_MAKE_BLOCKS + (TEST_HOOKS_BLOCK, DIRECT_ENCODE_BLOCK),
        "make_normalized_targets": frozenset(
            target for target, _prerequisites in
            COMMON_MAKE_BLOCKS + (TEST_HOOKS_BLOCK, DIRECT_ENCODE_BLOCK)),
        "make_template_sha256":
            "c56eb276850e7acbc5242e91c26b8fecf993040feeeb88b56596a189924a8479",
        "node_census_count": 556,
        "node_census_sha256":
            "2fa17fb9bb1458ea2bb7c361de077323e6f0a55359ba137eb445d231dd13d24b",
        "profile_sha256":
            "bead85b797c0e9966d4f09ba5b18b48739a350f550f831d1afa5bd2ae68dfa69",
        "single_c_record": {
            "output":
                "CMakeFiles/leopard2_codec_options_abi_test.dir/"
                "tests/leopard2/test_codec_options_abi.c.o",
            "source": "tests/leopard2/test_codec_options_abi.c",
        },
        "template_origin":
            "rederived from sealed v1/v2 candidate correctness closures",
    },
}


class ContractError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise ContractError(message)


def sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def canonical_json(value, final_lf=True):
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return encoded + (b"\n" if final_lf else b"")


def profile_document(profile_name, profile):
    return {
        "compile_commands": {
            "compiler_counts": profile["compiler_counts"],
            "entry_count": profile["compile_entry_count"],
            "output_list_sha256": profile["compile_output_list_sha256"],
            "single_c_record": profile["single_c_record"],
            "whole_file_template_sha256":
                profile["compile_template_sha256"],
        },
        "makefile2": {
            "normalized_targets":
                sorted(profile["make_normalized_targets"]),
            "validated_blocks": [
                {"prerequisites": list(prerequisites), "target": target}
                for target, prerequisites in profile["make_blocks"]
            ],
            "whole_file_template_sha256": profile["make_template_sha256"],
        },
        "name": profile_name,
        "node_census": {
            "count": profile["node_census_count"],
            "relative_path_type_mode_sha256":
                profile["node_census_sha256"],
        },
        "schema": "leopard2-v17-build-order-profile/v2",
        "template_path_tokens": {
            "build": TEMPLATE_BUILD.decode("ascii"),
            "commit": TEMPLATE_COMMIT.decode("ascii"),
            "source": TEMPLATE_SOURCE.decode("ascii"),
        },
    }


def profile_contract(profile_name):
    require(profile_name in PROFILES, "unknown normalization profile")
    profile = PROFILES[profile_name]
    require(profile["compile_entry_count"] ==
            sum(profile["compiler_counts"].values()),
            "normalization profile compiler census is inconsistent")
    require(set(profile["make_normalized_targets"]).issubset(
            {target for target, _dependencies in profile["make_blocks"]}),
            "normalization profile target map is inconsistent")
    require(isinstance(profile["node_census_count"], int) and
            profile["node_census_count"] > 0,
            "normalization profile node census is invalid")
    for key in (
            "compile_output_list_sha256", "compile_template_sha256",
            "make_template_sha256", "node_census_sha256", "profile_sha256"):
        require(re.fullmatch(r"[0-9a-f]{64}", profile[key]) is not None,
                "normalization profile has an invalid template digest")
    require(sha256_bytes(canonical_json(
        profile_document(profile_name, profile))) == profile["profile_sha256"],
        "normalization profile digest differs")
    return profile


def exact_commit(value, label):
    require(isinstance(value, str) and
            re.fullmatch(r"[0-9a-f]{40}", value) is not None,
            label + " is not an exact lowercase 40-hex commit")
    encoded = value.encode("ascii")
    require(all(token not in encoded for token in (
                TEMPLATE_BUILD, TEMPLATE_SOURCE, TEMPLATE_COMMIT)),
            label + " collides with a reserved template token")
    return value


def canonical_template(
        value, expected_build, expected_source, expected_commit, label,
        commit_required=False):
    build = expected_build.encode("utf-8")
    source = expected_source.encode("utf-8")
    commit = exact_commit(expected_commit, "expected commit").encode("ascii")
    require(build != source and build not in source and source not in build,
            label + " has overlapping expected path prefixes")
    require(commit not in build and commit not in source,
            label + " has an expected commit embedded in a path prefix")
    require(all(token not in value for token in (
                TEMPLATE_BUILD, TEMPLATE_SOURCE, TEMPLATE_COMMIT)),
            label + " already contains a reserved template token")
    require(build in value and source in value,
            label + " omits an expected path prefix")
    if commit_required:
        require(commit in value, label + " omits the expected commit")
    templated = value.replace(build, TEMPLATE_BUILD)
    templated = templated.replace(source, TEMPLATE_SOURCE)
    templated = templated.replace(commit, TEMPLATE_COMMIT)
    if commit_required:
        require(TEMPLATE_COMMIT in templated,
                label + " did not produce the expected commit token")
    require(build not in templated and source not in templated and
            commit not in templated,
            label + " retained an expected input after templating")
    return templated


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def exact_absolute(value, label):
    require(isinstance(value, str) and value, label + " is not a string")
    require("\x00" not in value and "\n" not in value and "\r" not in value,
            label + " contains a control separator")
    require(os.path.isabs(value) and os.path.normpath(value) == value,
            label + " is not an exact normalized absolute path")
    return value


def render_cmake_compile_commands(records):
    pieces = ["[\n"]
    ordered_keys = ("directory", "command", "file", "output")
    for index, record in enumerate(records):
        pieces.append("{\n")
        for key_index, key in enumerate(ordered_keys):
            pieces.append("  " + json.dumps(key) + ": " +
                          json.dumps(record[key], ensure_ascii=False))
            pieces.append(",\n" if key_index + 1 < len(ordered_keys) else "\n")
        pieces.append("},\n" if index + 1 < len(records) else "}\n")
    pieces.append("]")
    return "".join(pieces).encode("utf-8")


def stable_bytes(path, label):
    path = Path(path)
    metadata = os.lstat(path)
    require(stat.S_ISREG(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode),
            label + " is not a regular file")
    require(metadata.st_uid == os.geteuid() and metadata.st_gid == os.getegid(),
            label + " has unexpected ownership")
    require(metadata.st_nlink == 1, label + " is hard-linked")
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        identity = (
            metadata.st_dev, metadata.st_ino, metadata.st_mode,
            metadata.st_uid, metadata.st_gid, metadata.st_nlink,
            metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns,
        )
        require(identity == (
            opened.st_dev, opened.st_ino, opened.st_mode, opened.st_uid,
            opened.st_gid, opened.st_nlink, opened.st_size,
            opened.st_mtime_ns, opened.st_ctime_ns,
        ), label + " changed before read")
        remaining = metadata.st_size
        chunks = []
        while remaining:
            block = os.read(descriptor, min(1024 * 1024, remaining))
            require(block, label + " had a short read")
            chunks.append(block)
            remaining -= len(block)
        require(not os.read(descriptor, 1), label + " grew during read")
        final = os.fstat(descriptor)
        require(identity == (
            final.st_dev, final.st_ino, final.st_mode, final.st_uid,
            final.st_gid, final.st_nlink, final.st_size,
            final.st_mtime_ns, final.st_ctime_ns,
        ), label + " changed during read")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def snapshot_tree(root, label):
    root = Path(root)
    require(root.is_absolute() and os.path.normpath(str(root)) == str(root),
            label + " root is not an exact absolute path")
    metadata = os.lstat(root)
    require(stat.S_ISDIR(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode),
            label + " root is not a directory")
    require(metadata.st_uid == os.geteuid() and metadata.st_gid == os.getegid(),
            label + " root has unexpected ownership")
    require(not stat.S_IMODE(metadata.st_mode) & 0o222,
            label + " root remains writable")

    def node_record(relative, kind, item):
        return {
            "gid": item.st_gid,
            "mode": format(stat.S_IMODE(item.st_mode), "04o"),
            "nlink": item.st_nlink,
            "path": relative,
            "type": kind,
            "uid": item.st_uid,
        }

    nodes = {".": node_record(".", "directory", metadata)}
    file_bytes = {}

    def visit(directory, relative_parent):
        with os.scandir(directory) as stream:
            entries = sorted(stream, key=lambda entry: entry.name)
        for entry in entries:
            relative = (relative_parent + "/" + entry.name
                        if relative_parent else entry.name)
            path = Path(directory) / entry.name
            item = os.lstat(path)
            require(item.st_uid == os.geteuid() and item.st_gid == os.getegid(),
                    label + " node has unexpected ownership: " + relative)
            require(not stat.S_IMODE(item.st_mode) & 0o222,
                    label + " node remains writable: " + relative)
            if stat.S_ISDIR(item.st_mode) and not stat.S_ISLNK(item.st_mode):
                nodes[relative] = node_record(relative, "directory", item)
                visit(path, relative)
            elif stat.S_ISREG(item.st_mode) and not stat.S_ISLNK(item.st_mode):
                require(item.st_nlink == 1,
                        label + " file is hard-linked: " + relative)
                if entry.name in {"compile_commands.json", "Makefile2"}:
                    require(relative in EXCEPTIONS,
                            label + " has a non-allowlisted ordering file: " +
                            relative)
                nodes[relative] = node_record(relative, "file", item)
                file_bytes[relative] = stable_bytes(path, label + " " + relative)
            else:
                raise ContractError(label + " has a special node: " + relative)

    visit(root, "")
    require(all(path in file_bytes for path in EXCEPTIONS),
            label + " omits an exact ordering exception")
    return nodes, file_bytes


def relative_node_census(nodes):
    return [
        {
            "mode": nodes[path]["mode"],
            "path": path,
            "type": nodes[path]["type"],
        }
        for path in sorted(nodes)
    ]


def normalize_compile_commands(
        raw, expected_build, expected_source, expected_commit, profile_name):
    profile = profile_contract(profile_name)
    require(raw and not raw.endswith(b"\n") and b"\r" not in raw,
            "compile_commands.json source framing changed")
    try:
        value = json.loads(
            raw.decode("utf-8"), object_pairs_hook=strict_object,
            parse_constant=lambda token: (_ for _ in ()).throw(
                ContractError("non-finite JSON token: " + token)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError("compile_commands.json is not strict UTF-8 JSON") from error
    require(isinstance(value, list) and
            len(value) == profile["compile_entry_count"],
            "compile_commands.json entry count differs for profile " +
            profile_name)
    require(render_cmake_compile_commands(value) == raw,
            "compile_commands.json differs from the exact CMake formatter")
    outputs = set()
    compiler_counts = Counter()
    c_records = []
    normalized = []
    source_prefix = expected_source + os.sep
    for index, record in enumerate(value):
        label = "compile command {}".format(index)
        require(isinstance(record, dict) and set(record) == COMPILE_KEYS,
                label + " does not have the exact four-key schema")
        require(all(isinstance(record[key], str) and record[key]
                    for key in COMPILE_KEYS),
                label + " has a non-string or empty member")
        directory = exact_absolute(record["directory"], label + " directory")
        source_file = exact_absolute(record["file"], label + " file")
        require(directory == expected_build,
                label + " directory differs from the canonical build path")
        require(source_file.startswith(source_prefix),
                label + " file escapes the canonical source root")
        output = record["output"]
        require("\x00" not in output and "\n" not in output and
                "\r" not in output and not os.path.isabs(output) and
                os.path.normpath(output) == output and
                output not in (".", "..") and
                not output.startswith(".." + os.sep),
                label + " output is not an exact safe relative path")
        require(output not in outputs,
                "compile_commands.json has a duplicate output: " + output)
        outputs.add(output)
        command = record["command"]
        require("\x00" not in command and "\n" not in command and
                "\r" not in command,
                label + " command contains a control separator")
        try:
            tokens = shlex.split(command, posix=True)
        except ValueError as error:
            raise ContractError(label + " command is not shell-tokenizable") from error
        require(tokens and tokens[0] in profile["compiler_counts"] and
                tokens.count("-o") == 1 and tokens.count("-c") == 1,
                label + " command has an unexpected compiler/output/source form")
        compiler_counts[tokens[0]] += 1
        output_index = tokens.index("-o")
        source_index = tokens.index("-c")
        require(output_index + 1 < len(tokens) and
                tokens[output_index + 1] == output and
                source_index + 1 < len(tokens) and
                tokens[source_index + 1] == source_file,
                label + " command/output/source fields disagree")
        if tokens[0] == "/usr/bin/cc":
            c_records.append({"output": output, "source": source_file})
        normalized.append(record)
    require(dict(compiler_counts) == profile["compiler_counts"],
            "compile_commands.json compiler census differs for profile " +
            profile_name)
    expected_c_record = profile["single_c_record"]
    if expected_c_record is None:
        require(not c_records,
                "compile_commands.json unexpectedly contains a C compiler record")
    else:
        require(c_records == [{
            "output": expected_c_record["output"],
            "source": expected_source + os.sep + expected_c_record["source"],
        }], "compile_commands.json C compiler record differs")
    normalized.sort(key=lambda record: record["output"])
    canonical = render_cmake_compile_commands(normalized)
    templated = canonical_template(
        canonical, expected_build, expected_source, expected_commit,
        "compile_commands.json canonical template", commit_required=True)
    output_list = canonical_json(sorted(outputs))
    template_sha256 = sha256_bytes(templated)
    output_list_sha256 = sha256_bytes(output_list)
    require(template_sha256 == profile["compile_template_sha256"],
            "compile_commands.json closed-world template differs for profile " +
            profile_name)
    require(output_list_sha256 == profile["compile_output_list_sha256"],
            "compile_commands.json output census differs for profile " +
            profile_name)
    return canonical, templated, output_list, {
        "compiler_counts": dict(sorted(compiler_counts.items())),
        "entry_count": len(normalized),
        "exact_cmake_rendering": True,
        "expected_output_list_sha256":
            profile["compile_output_list_sha256"],
        "expected_template_sha256": profile["compile_template_sha256"],
        "output_list_sha256": output_list_sha256,
        "source_final_lf": False,
        "template_algorithm":
            "canonicalize then replace exact expected path bytes with "
            "${BUILD} and ${SOURCE}, then replace the exact expected commit "
            "with ${COMMIT}",
        "template_sha256": template_sha256,
        "template_size": len(templated),
        "unique_output_count": len(outputs),
        "whole_file_template": True,
    }


def normalize_makefile2(
        raw, expected_build, expected_source, expected_commit, profile_name):
    profile = profile_contract(profile_name)
    require(raw.endswith(b"\n") and b"\r" not in raw and b"\x00" not in raw,
            "Makefile2 source framing changed")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ContractError("Makefile2 is not UTF-8") from error
    lines = text.splitlines(keepends=True)
    require(lines and all(line.endswith("\n") for line in lines),
            "Makefile2 has a non-LF-terminated line")
    require(lines[0] == "# CMAKE generated file: DO NOT EDIT!\n",
            "Makefile2 header differs")
    require(lines.count("CMAKE_SOURCE_DIR = " + expected_source + "\n") == 1,
            "Makefile2 source path binding differs")
    require(lines.count("CMAKE_BINARY_DIR = " + expected_build + "\n") == 1,
            "Makefile2 build path binding differs")
    require(ORDER_SENSITIVE_AUTOMATIC.search(text) is None,
            "Makefile2 contains an order-sensitive automatic variable")
    expected_blocks = {target: tuple(prerequisites)
                       for target, prerequisites in profile["make_blocks"]}
    all_positions = defaultdict(list)
    for index, line in enumerate(lines):
        rule_like = MAKE_RULE_LIKE.fullmatch(line)
        if rule_like is not None:
            require(MAKE_RULE.fullmatch(line) is not None,
                    "Makefile2 has a noncanonical .dir/all dependency rule: " +
                    rule_like.group(1))
        match = MAKE_RULE.fullmatch(line)
        if match:
            all_positions[match.group(1)].append((index, match.group(2)))
    multi_prerequisite_targets = {
        target for target, observed in all_positions.items()
        if len(observed) >= 2
    }
    require(multi_prerequisite_targets == set(expected_blocks),
            "Makefile2 multi-prerequisite target allowlist is incomplete")
    positions = {target: all_positions[target] for target in expected_blocks}
    output = list(lines)
    block_report = []
    for target, prerequisites in profile["make_blocks"]:
        observed = positions[target]
        require(len(observed) == len(prerequisites),
                "Makefile2 allowlisted block length differs: " + target)
        indexes = [index for index, _dependency in observed]
        require(indexes == list(range(indexes[0], indexes[0] + len(indexes))),
                "Makefile2 allowlisted block is not consecutive: " + target)
        observed_dependencies = [dependency for _index, dependency in observed]
        require(len(set(observed_dependencies)) == len(observed_dependencies) and
                set(observed_dependencies) == set(prerequisites),
                "Makefile2 has a moved/nonallowlisted prerequisite: " + target)
        start = indexes[0]
        end = start + len(indexes)
        phony = ".PHONY : " + target + "\n"
        phony_indexes = [index for index, line in enumerate(lines)
                         if line == phony]
        require(len(phony_indexes) == 1 and phony_indexes[0] > end,
                "Makefile2 allowlisted target has no unique later PHONY marker: " +
                target)
        recipe = lines[end:phony_indexes[0]]
        require(recipe and all(line.startswith("\t") for line in recipe),
                "Makefile2 allowlisted target recipe context differs: " + target)
        require(ORDER_SENSITIVE_AUTOMATIC.search("".join(recipe)) is None,
                "Makefile2 allowlisted recipe is order-sensitive: " + target)
        normalized = target in profile["make_normalized_targets"]
        if normalized:
            output[start:end] = sorted(output[start:end])
        block_report.append({
            "line_count": len(indexes),
            "normalized": normalized,
            "phony_marker_count": len(phony_indexes),
            "recipe_line_count": len(recipe),
            "target": target,
        })
    normalized_count = sum(block["normalized"] for block in block_report)
    require(len(block_report) == len(profile["make_blocks"]),
            "Makefile2 validated block census differs for profile " +
            profile_name)
    require(normalized_count == len(profile["make_normalized_targets"]),
            "Makefile2 normalized block census differs for profile " +
            profile_name)
    canonical = "".join(output).encode("utf-8")
    templated = canonical_template(
        canonical, expected_build, expected_source, expected_commit,
        "Makefile2 canonical template")
    template_sha256 = sha256_bytes(templated)
    require(template_sha256 == profile["make_template_sha256"],
            "Makefile2 closed-world template differs for profile " +
            profile_name)
    return canonical, templated, {
        "all_allowlisted_lhs_colon_forms_validated": True,
        "all_multi_prerequisite_targets_allowlisted": True,
        "expected_template_sha256": profile["make_template_sha256"],
        "multi_prerequisite_target_count":
            len(multi_prerequisite_targets),
        "multi_prerequisite_targets":
            sorted(multi_prerequisite_targets),
        "normalized_block_count": normalized_count,
        "source_final_lf": True,
        "template_algorithm":
            "canonicalize then replace exact expected path bytes with "
            "${BUILD} and ${SOURCE}, then replace the exact expected commit "
            "with ${COMMIT}",
        "template_sha256": template_sha256,
        "template_size": len(templated),
        "validated_block_count": len(block_report),
        "validated_blocks": block_report,
        "whole_file_template": True,
    }


def write_exclusive(path, value):
    path = Path(path)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        view = memoryview(value)
        while view:
            written = os.write(descriptor, view)
            require(written > 0, "short canonical output write")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def compare_pair(
        left, right, expected_build, expected_source, expected_commit,
        output_dir, profile_name):
    profile = profile_contract(profile_name)
    expected_build = exact_absolute(expected_build, "expected build path")
    expected_source = exact_absolute(expected_source, "expected source path")
    expected_commit = exact_commit(expected_commit, "expected commit")
    left = Path(exact_absolute(str(left), "left closure path"))
    right = Path(exact_absolute(str(right), "right closure path"))
    output_dir = Path(exact_absolute(str(output_dir), "normalization output path"))
    require(not output_dir.exists(), "normalization output already exists")
    output_dir.mkdir(mode=0o700)
    left_nodes, left_files = snapshot_tree(left, "left closure")
    right_nodes, right_files = snapshot_tree(right, "right closure")
    require(left_nodes == right_nodes,
            "closure path/type/mode/ownership/link census differs")
    node_census = relative_node_census(left_nodes)
    node_census_sha256 = sha256_bytes(canonical_json(node_census))
    require(len(node_census) == profile["node_census_count"] and
            node_census_sha256 == profile["node_census_sha256"],
            "closed-world relative node/type/mode census differs for profile " +
            profile_name)
    for relative in sorted(left_files):
        if relative not in EXCEPTIONS:
            require(left_files[relative] == right_files[relative],
                    "nonexception closure bytes differ: " + relative)

    left_compile, left_compile_template, left_compile_outputs, \
        left_compile_report = \
        normalize_compile_commands(
            left_files[EXCEPTIONS[0]], expected_build, expected_source,
            expected_commit, profile_name)
    right_compile, right_compile_template, right_compile_outputs, \
        right_compile_report = \
        normalize_compile_commands(
            right_files[EXCEPTIONS[0]], expected_build, expected_source,
            expected_commit, profile_name)
    left_make, left_make_template, left_make_report = normalize_makefile2(
        left_files[EXCEPTIONS[1]], expected_build, expected_source,
        expected_commit, profile_name)
    right_make, right_make_template, right_make_report = normalize_makefile2(
        right_files[EXCEPTIONS[1]], expected_build, expected_source,
        expected_commit, profile_name)

    outputs = {
        "left-compile_commands.canonical.json": left_compile,
        "left-compile_outputs.canonical.json": left_compile_outputs,
        "left-compile_commands.template.json": left_compile_template,
        "right-compile_commands.canonical.json": right_compile,
        "right-compile_outputs.canonical.json": right_compile_outputs,
        "right-compile_commands.template.json": right_compile_template,
        "left-Makefile2.canonical": left_make,
        "left-Makefile2.template": left_make_template,
        "right-Makefile2.canonical": right_make,
        "right-Makefile2.template": right_make_template,
    }
    for name, value in outputs.items():
        write_exclusive(output_dir / name, value)
    require(left_compile == right_compile,
            "compile_commands canonical semantics differ")
    require(left_compile_outputs == right_compile_outputs,
            "compile_commands output censuses differ")
    require(left_compile_template == right_compile_template,
            "compile_commands closed-world templates differ")
    require(left_make == right_make, "Makefile2 canonical semantics differ")
    require(left_make_template == right_make_template,
            "Makefile2 closed-world templates differ")

    raw_equal = all(left_files[path] == right_files[path]
                    for path in left_files)
    inventory = [left_nodes[path] for path in sorted(left_nodes)]
    exact_byte_files = [
        {
            "path": path,
            "sha256": sha256_bytes(left_files[path]),
            "size": len(left_files[path]),
        }
        for path in sorted(left_files)
        if path not in EXCEPTIONS
    ]
    report = {
        "compile_commands": {
            "canonical_sha256": sha256_bytes(left_compile),
            "left": left_compile_report,
            "left_raw_sha256": sha256_bytes(left_files[EXCEPTIONS[0]]),
            "left_raw_size": len(left_files[EXCEPTIONS[0]]),
            "raw_byte_identical":
                left_files[EXCEPTIONS[0]] == right_files[EXCEPTIONS[0]],
            "right": right_compile_report,
            "right_raw_sha256": sha256_bytes(right_files[EXCEPTIONS[0]]),
            "right_raw_size": len(right_files[EXCEPTIONS[0]]),
        },
        "contract": {
            "all_other_relative_file_bytes_identical": True,
            "closed_world_templates_verified": True,
            "compile_commands_order_only":
                "sort exactly {} unique exact four-key records by output".format(
                    profile["compile_entry_count"]),
            "exception_paths": list(EXCEPTIONS),
            "makefile2_order_only":
                "sort only exact profile-allowlisted consecutive dependency blocks",
            "operation": "compare",
            "profile": profile_name,
            "raw_tree_file_bytes_identical": raw_equal,
            "template_path_tokens": {
                "build": TEMPLATE_BUILD.decode("ascii"),
                "commit": TEMPLATE_COMMIT.decode("ascii"),
                "source": TEMPLATE_SOURCE.decode("ascii"),
            },
        },
        "exact_byte_files": {
            "count": len(exact_byte_files),
            "manifest": exact_byte_files,
            "manifest_sha256": sha256_bytes(canonical_json(exact_byte_files)),
        },
        "expected_build": expected_build,
        "expected_commit": expected_commit,
        "expected_source": expected_source,
        "makefile2": {
            "canonical_sha256": sha256_bytes(left_make),
            "left": left_make_report,
            "left_raw_sha256": sha256_bytes(left_files[EXCEPTIONS[1]]),
            "left_raw_size": len(left_files[EXCEPTIONS[1]]),
            "raw_byte_identical":
                left_files[EXCEPTIONS[1]] == right_files[EXCEPTIONS[1]],
            "right": right_make_report,
            "right_raw_sha256": sha256_bytes(right_files[EXCEPTIONS[1]]),
            "right_raw_size": len(right_files[EXCEPTIONS[1]]),
        },
        "node_inventory": {
            "count": len(inventory),
            "manifest": inventory,
            "manifest_sha256": sha256_bytes(canonical_json(inventory)),
            "profile_relative_path_type_mode_count": len(node_census),
            "profile_relative_path_type_mode_manifest": node_census,
            "profile_relative_path_type_mode_sha256": node_census_sha256,
        },
        "profile_contract": {
            "compile_entry_count": profile["compile_entry_count"],
            "compile_output_list_sha256":
                profile["compile_output_list_sha256"],
            "compile_template_sha256": profile["compile_template_sha256"],
            "compiler_counts": profile["compiler_counts"],
            "make_normalized_targets":
                sorted(profile["make_normalized_targets"]),
            "make_template_sha256": profile["make_template_sha256"],
            "make_validated_targets": [
                target for target, _dependencies in profile["make_blocks"]],
            "node_census_count": profile["node_census_count"],
            "node_census_sha256": profile["node_census_sha256"],
            "profile_sha256": profile["profile_sha256"],
            "profile_document": profile_document(profile_name, profile),
            "single_c_record": profile["single_c_record"],
            "template_origin": profile["template_origin"],
        },
        "schema": SCHEMA,
        "semantic_equal": True,
    }
    write_exclusive(output_dir / "report.json", canonical_json(report))
    return report


def fixture_compile(
        build, source, expected_commit, base_profile_name, reverse=False):
    profile = profile_contract(base_profile_name)
    records = []
    c_record = profile["single_c_record"]
    c_index = profile["compile_entry_count"] - 1 if c_record else None
    for index in range(profile["compile_entry_count"]):
        if index == c_index:
            compiler = "/usr/bin/cc"
            source_file = source + os.sep + c_record["source"]
            output = c_record["output"]
        else:
            compiler = "/usr/bin/c++"
            source_file = source + "/fixture/f{:03d}.cpp".format(index)
            output = (
                "CMakeFiles/fixture_{:03d}.dir/fixture/f{:03d}.cpp.o".format(
                    index, index))
        records.append({
            "command":
                "{} -DLEO2_FIXTURE_SOURCE_GIT_SHA={} -DVALUE={} "
                "-o {} -c {}".format(
                    compiler, expected_commit, index, output, source_file),
            "directory": build,
            "file": source_file,
            "output": output,
        })
    if reverse:
        records.reverse()
    return render_cmake_compile_commands(records)


def fixture_make(build, source, base_profile_name, reverse=False):
    profile = profile_contract(base_profile_name)
    lines = [
        "# CMAKE generated file: DO NOT EDIT!\n",
        "# Generated by fixture\n",
        "CMAKE_SOURCE_DIR = " + source + "\n",
        "CMAKE_BINARY_DIR = " + build + "\n",
        "\n",
    ]
    for target, prerequisites in profile["make_blocks"]:
        ordered = list(prerequisites)
        if reverse and target in profile["make_normalized_targets"]:
            ordered.reverse()
        lines.extend(target + ": " + dependency + "\n"
                     for dependency in ordered)
        lines.extend((
            "\t@echo stable recipe\n",
            ".PHONY : " + target + "\n",
            "\n",
        ))
    return "".join(lines).encode("utf-8")


def install_fixture(
        root, build, source, expected_commit, base_profile_name,
        reverse=False):
    root.mkdir()
    (root / "CMakeFiles").mkdir()
    (root / "compile_commands.json").write_bytes(
        fixture_compile(
            build, source, expected_commit, base_profile_name, reverse))
    (root / "CMakeFiles/Makefile2").write_bytes(
        fixture_make(build, source, base_profile_name, reverse))
    (root / "selected-object.o").write_bytes(b"exact-object-bytes\x00\xff")


def install_fixture_profile(
        profile_name, base_profile_name, build, source, expected_commit,
        fixture_root):
    require(profile_name not in PROFILES, "fixture profile already exists")
    base_profile = profile_contract(base_profile_name)
    profile = dict(base_profile)
    records = json.loads(fixture_compile(
        build, source, expected_commit, base_profile_name, False))
    compile_canonical = render_cmake_compile_commands(
        sorted(records, key=lambda record: record["output"]))
    make_canonical = fixture_make(
        build, source, base_profile_name, False)
    profile["compile_output_list_sha256"] = sha256_bytes(
        canonical_json(sorted(record["output"] for record in records)))
    profile["compile_template_sha256"] = sha256_bytes(canonical_template(
        compile_canonical, build, source, expected_commit,
        "fixture compile template", commit_required=True))
    profile["make_template_sha256"] = sha256_bytes(canonical_template(
        make_canonical, build, source, expected_commit,
        "fixture Makefile2 template"))
    nodes, _files = snapshot_tree(fixture_root, "fixture profile closure")
    node_census = relative_node_census(nodes)
    profile["node_census_count"] = len(node_census)
    profile["node_census_sha256"] = sha256_bytes(canonical_json(node_census))
    profile["template_origin"] = "embedded adversarial self-test fixture"
    profile["profile_sha256"] = sha256_bytes(canonical_json(
        profile_document(profile_name, profile)))
    PROFILES[profile_name] = profile


def set_fixture_writable(root):
    for directory, directories, files in os.walk(root, topdown=True):
        os.chmod(directory, stat.S_IMODE(os.lstat(directory).st_mode) | 0o200)
        for name in files:
            path = Path(directory) / name
            if not path.is_symlink():
                os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) | 0o200)


def seal_fixture(root):
    for directory, directories, files in os.walk(root, topdown=False):
        for name in files:
            path = Path(directory) / name
            if not path.is_symlink():
                os.chmod(path, stat.S_IMODE(os.lstat(path).st_mode) & ~0o222)
        os.chmod(directory, stat.S_IMODE(os.lstat(directory).st_mode) & ~0o222)


def self_test():
    with tempfile.TemporaryDirectory(prefix="leopard-v17-order-self-test.") as tmp:
        base = Path(tmp)
        build = "/tmp/leopard-v17-order-self-test/build"
        source = "/tmp/leopard-v17-order-self-test/source"
        expected_commit = "0123456789abcdef0123456789abcdef01234567"
        alternate_commit = "89abcdef0123456789abcdef0123456789abcdef"
        timing_profile = "self-test-candidate-timing"
        tests_profile = "self-test-candidate-tests"

        closures = {}
        for label, base_profile_name, fixture_profile_name in (
                ("timing", "candidate-timing", timing_profile),
                ("tests", "candidate-tests", tests_profile)):
            left = base / (label + "-left")
            reordered = base / (label + "-reordered")
            install_fixture(
                left, build, source, expected_commit, base_profile_name, False)
            install_fixture(
                reordered, build, source, expected_commit, base_profile_name,
                True)
            seal_fixture(left)
            seal_fixture(reordered)
            install_fixture_profile(
                fixture_profile_name, base_profile_name, build, source,
                expected_commit, left)
            report = compare_pair(
                left, reordered, build, source, expected_commit,
                base / (label + "-accepted"), fixture_profile_name)
            expected_validated = 4 if label == "timing" else 6
            expected_normalized = 4 if label == "timing" else 6
            require(report["semantic_equal"] is True and
                    report["contract"]["profile"] == fixture_profile_name and
                    report["contract"]["operation"] == "compare" and
                    report["contract"]["raw_tree_file_bytes_identical"] is
                    False and
                    report["compile_commands"]["left"]["entry_count"] ==
                    PROFILES[fixture_profile_name]["compile_entry_count"] and
                    report["makefile2"]["left"]["validated_block_count"] ==
                    expected_validated and
                    report["makefile2"]["left"]["normalized_block_count"] ==
                    expected_normalized and
                    report["makefile2"]["left"]
                        ["all_multi_prerequisite_targets_allowlisted"] is
                    True and
                    report["makefile2"]["left"]
                        ["multi_prerequisite_target_count"] ==
                    expected_validated and
                    report["contract"]["closed_world_templates_verified"] is
                    True and
                    report["contract"]["template_path_tokens"] == {
                        "build": "${BUILD}", "commit": "${COMMIT}",
                        "source": "${SOURCE}"} and
                    report["expected_commit"] == expected_commit,
                    label + " reorder-only fixture was not qualified honestly")
            alternate_raw = fixture_compile(
                build, source, alternate_commit, base_profile_name, False)
            _alternate_canonical, alternate_template, \
                _alternate_outputs, _alternate_report = \
                normalize_compile_commands(
                    alternate_raw, build, source, alternate_commit,
                    fixture_profile_name)
            require(alternate_template == (base / (
                        label + "-accepted/left-compile_commands.template.json"
                    )).read_bytes(),
                    label + " compile template is not commit-invariant")
            closures[label] = (left, reordered, fixture_profile_name)

        def rejected(
                name, label="timing", mutate_left=None, mutate_right=None,
                profile_name=None, expected_commit_override=None):
            clean_left, clean_right, default_profile = closures[label]
            left = base / (name + "-left")
            right = base / (name + "-right")
            shutil.copytree(clean_left, left, symlinks=True)
            shutil.copytree(clean_right, right, symlinks=True)
            set_fixture_writable(left)
            set_fixture_writable(right)
            if mutate_left is not None:
                mutate_left(left)
            if mutate_right is not None:
                mutate_right(right)
            seal_fixture(left)
            seal_fixture(right)
            try:
                compare_pair(
                    left, right, build, source,
                    expected_commit_override or expected_commit,
                    base / (name + "-out"),
                    profile_name or default_profile)
            except ContractError:
                return
            raise ContractError("adversarial fixture was accepted: " + name)

        def duplicate_key(root):
            path = root / "compile_commands.json"
            raw = path.read_bytes()
            raw = raw.replace(
                b'"directory":', b'"directory":"duplicate","directory":', 1)
            path.write_bytes(raw)

        def duplicate_output(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            records[1]["output"] = records[0]["output"]
            path.write_bytes(render_cmake_compile_commands(records))

        def count_mutation(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            records.pop()
            path.write_bytes(render_cmake_compile_commands(records))

        def compiler_action_mutation(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            records[0]["command"] = records[0]["command"].replace(
                " -c ", " -E ", 1)
            path.write_bytes(render_cmake_compile_commands(records))

        def token_mutation(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            record = min(records, key=lambda item: item["output"])
            record["command"] = record["command"].replace(
                "/usr/bin/c++ ", "/usr/bin/c++ -DSEMANTIC_MUTATION=1 ", 1)
            path.write_bytes(render_cmake_compile_commands(records))

        def commit_token_collision(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            record = min(records, key=lambda item: item["output"])
            record["command"] = record["command"].replace(
                "/usr/bin/c++ ", "/usr/bin/c++ -DRESERVED=${COMMIT} ", 1)
            path.write_bytes(render_cmake_compile_commands(records))

        def c_record_mutation(root):
            path = root / "compile_commands.json"
            records = json.loads(path.read_bytes())
            record = next(record for record in records
                          if record["command"].startswith("/usr/bin/cc "))
            changed = record["file"] + ".changed"
            record["command"] = record["command"].replace(
                record["file"], changed)
            record["file"] = changed
            path.write_bytes(render_cmake_compile_commands(records))

        def nonallowlisted_prerequisite(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes()
            raw = raw.replace(
                b"CMakeFiles/leopard2_backend_avx2.dir/all",
                b"CMakeFiles/not_allowlisted.dir/all", 1)
            path.write_bytes(raw)

        def moved_prerequisite(root):
            path = root / "CMakeFiles/Makefile2"
            lines = path.read_text().splitlines(keepends=True)
            needle = ("CMakeFiles/leopard.dir/all: "
                      "CMakeFiles/leopard2_backend_avx2.dir/all\n")
            lines.remove(needle)
            lines.append(needle)
            path.write_text("".join(lines))

        def alternate_colon_form(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes().replace(
                b"CMakeFiles/leopard.dir/all: ",
                b"CMakeFiles/leopard.dir/all : ", 1)
            path.write_bytes(raw)

        def recipe_mutation(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes().replace(
                b"\t@echo stable recipe\n", b"\t@echo changed recipe\n", 1)
            path.write_bytes(raw)

        def automatic_variable(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes().replace(
                b"\t@echo stable recipe\n", b"\t@echo $^\n", 1)
            path.write_bytes(raw)

        def automatic_variable_derived(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes().replace(
                b"\t@echo stable recipe\n", b"\t@echo $(^D)\n", 1)
            path.write_bytes(raw)

        def duplicate_phony(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes()
            marker = b".PHONY : CMakeFiles/leopard.dir/all\n"
            raw = raw.replace(marker, marker + marker, 1)
            path.write_bytes(raw)

        def mode_drift(root):
            os.chmod(root / "selected-object.o", 0o440)

        def extra_node(root):
            (root / "unexpected.bin").write_bytes(b"unexpected")

        def unallowlisted_multi_prerequisite(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes()
            raw += (
                b"\nCMakeFiles/unallowlisted.dir/all: "
                b"CMakeFiles/unexpected_dep_a.dir/all\n"
                b"CMakeFiles/unallowlisted.dir/all: "
                b"CMakeFiles/unexpected_dep_b.dir/all\n"
                b"\t@echo stable recipe\n"
                b".PHONY : CMakeFiles/unallowlisted.dir/all\n"
            )
            path.write_bytes(raw)

        def unallowlisted_multi_prerequisite_alternate_colon(root):
            path = root / "CMakeFiles/Makefile2"
            raw = path.read_bytes()
            raw += (
                b"\nCMakeFiles/unallowlisted.dir/all : "
                b"CMakeFiles/unexpected_dep_a.dir/all\n"
                b"CMakeFiles/unallowlisted.dir/all : "
                b"CMakeFiles/unexpected_dep_b.dir/all\n"
                b"\t@echo stable recipe\n"
                b".PHONY : CMakeFiles/unallowlisted.dir/all\n"
            )
            path.write_bytes(raw)

        def special_node(root):
            os.symlink("selected-object.o", root / "unexpected-link")

        def hardlink_node(root):
            os.link(root / "selected-object.o", root / "unexpected-hardlink")

        def misplaced_exception(root):
            (root / "nested").mkdir()
            (root / "nested/compile_commands.json").write_bytes(b"[]")

        for name, mutation in (
            ("duplicate-key", duplicate_key),
            ("duplicate-output", duplicate_output),
            ("count-mutation", count_mutation),
            ("compiler-action-mutation", compiler_action_mutation),
            ("one-sided-semantic-mutation", token_mutation),
            ("commit-token-collision", commit_token_collision),
            ("nonallowlisted-prerequisite", nonallowlisted_prerequisite),
            ("moved-prerequisite", moved_prerequisite),
            ("alternate-colon-form", alternate_colon_form),
            ("recipe-mutation", recipe_mutation),
            ("automatic-variable", automatic_variable),
            ("automatic-variable-derived", automatic_variable_derived),
            ("duplicate-phony", duplicate_phony),
            ("mode-drift", mode_drift),
            ("extra-node", extra_node),
            ("special-node", special_node),
            ("hardlink-node", hardlink_node),
            ("misplaced-exception", misplaced_exception),
        ):
            rejected(name, mutate_right=mutation)

        for label in ("timing", "tests"):
            rejected(
                label + "-identical-compile-template-mutation", label=label,
                mutate_left=token_mutation, mutate_right=token_mutation)
            rejected(
                label + "-identical-make-template-mutation", label=label,
                mutate_left=recipe_mutation, mutate_right=recipe_mutation)
            rejected(
                label + "-identical-extra-node", label=label,
                mutate_left=extra_node, mutate_right=extra_node)
            rejected(
                label + "-identical-mode-drift", label=label,
                mutate_left=mode_drift, mutate_right=mode_drift)
            rejected(
                label + "-identical-unallowlisted-multi-prerequisite",
                label=label,
                mutate_left=unallowlisted_multi_prerequisite,
                mutate_right=unallowlisted_multi_prerequisite)
            rejected(
                label + "-identical-unallowlisted-multi-alternate-colon",
                label=label,
                mutate_left=unallowlisted_multi_prerequisite_alternate_colon,
                mutate_right=unallowlisted_multi_prerequisite_alternate_colon)
        rejected(
            "candidate-tests-c-record-mutation", label="tests",
            mutate_right=c_record_mutation)
        rejected(
            "timing-closure-under-tests-profile", label="timing",
            profile_name=tests_profile)
        rejected(
            "tests-closure-under-timing-profile", label="tests",
            profile_name=timing_profile)
        rejected(
            "uppercase-expected-commit",
            expected_commit_override=expected_commit.upper())
        rejected(
            "short-expected-commit",
            expected_commit_override=expected_commit[:-1])
        rejected(
            "wrong-expected-commit",
            expected_commit_override=alternate_commit)
        try:
            profile_contract("candidate-unknown")
        except ContractError:
            pass
        else:
            raise ContractError("unknown profile was accepted")
    print("v17 candidate build order normalizer self-test passed")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="operation", required=True)
    compare_parser = subparsers.add_parser("compare")
    compare_parser.add_argument("--left", required=True, type=Path)
    compare_parser.add_argument("--right", required=True, type=Path)
    compare_parser.add_argument("--expected-build", required=True)
    compare_parser.add_argument("--expected-commit", required=True)
    compare_parser.add_argument("--expected-source", required=True)
    compare_parser.add_argument("--output-dir", required=True, type=Path)
    compare_parser.add_argument(
        "--profile", required=True, choices=tuple(sorted(PROFILES)))
    subparsers.add_parser("self-test")
    options = parser.parse_args()
    try:
        if options.operation == "self-test":
            self_test()
        else:
            compare_pair(
                options.left, options.right, options.expected_build,
                options.expected_source, options.expected_commit,
                options.output_dir, options.profile)
    except ContractError as error:
        print("build order normalization rejected: " + str(error), file=os.sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY
    test "${PIPESTATUS[0]}" -eq 0 || return 1
    /usr/bin/chmod 0400 "$normalizer_output" || return 1
    return 0
}

capture_allowed_diff()
{
    local output=$1
    shift
    local diff_status=0
    if /usr/bin/diff "$@" > "$output"; then
        diff_status=0
    else
        diff_status=$?
        test "$diff_status" -eq 1
    fi
    /usr/bin/printf '%s\n' "$diff_status" > "$output.status"
}

normalized_evidence_files=(
    left-Makefile2.canonical
    left-Makefile2.template
    left-compile_commands.canonical.json
    left-compile_commands.template.json
    left-compile_outputs.canonical.json
    report.json
    right-Makefile2.canonical
    right-Makefile2.template
    right-compile_commands.canonical.json
    right-compile_commands.template.json
    right-compile_outputs.canonical.json
)

verify_normalized_evidence_directory()
{
    local normalized=$1
    local expected_census=
    local observed_census=
    test -d "$normalized" || return 1
    test ! -L "$normalized" || return 1
    require_empty_output /usr/bin/find "$normalized" \
        -mindepth 1 -maxdepth 1 ! -type f -print -quit || return 1
    require_empty_output /usr/bin/find "$normalized" \
        -mindepth 2 -print -quit || return 1
    require_empty_output /usr/bin/find "$normalized" \
        -type f -links +1 -print -quit || return 1
    expected_census="$(/usr/bin/printf '%s\n' \
        "${normalized_evidence_files[@]}" | /usr/bin/sort)" || return 1
    observed_census="$(/usr/bin/find "$normalized" \
        -mindepth 1 -maxdepth 1 -type f -printf '%f\n' | /usr/bin/sort)" \
        || return 1
    test "$observed_census" = "$expected_census" || return 1
}

compare_normalized_evidence_closures()
{
    local retained=$1
    local recomputed=$2
    verify_normalized_evidence_directory "$retained/normalized" || return 1
    verify_normalized_evidence_directory "$recomputed/normalized" || return 1
    /usr/bin/diff -qr "$retained/normalized" \
        "$recomputed/normalized" >/dev/null || return 1
    /usr/bin/cmp "$retained/canonical-SHA256SUMS" \
        "$recomputed/canonical-SHA256SUMS" || return 1
    (
        cd "$retained/normalized"
        /usr/bin/sha256sum -c ../canonical-SHA256SUMS
    ) >/dev/null || return 1
    (
        cd "$recomputed/normalized"
        /usr/bin/sha256sum -c ../canonical-SHA256SUMS
    ) >/dev/null || return 1
}

normalized_evidence_closure_self_test()
{
    local self_test_root=
    local left=
    local right=
    local relative=
    local failed_census=
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v17-normalized-evidence-self-test.XXXXXX)"
    if require_empty_output /usr/bin/bash -c 'exit 1'; then
        /usr/bin/printf '%s\n' \
            'empty output from a failed command was accepted' >&2
        return 1
    fi
    if failed_census="$(/usr/bin/find \
            "$self_test_root/missing-normalized-evidence" \
            -mindepth 1 -maxdepth 1 -type f -printf '%f\n' \
            2>/dev/null | /usr/bin/sort)"; then
        /usr/bin/printf '%s\n' \
            'empty output from a failed normalized census was accepted' >&2
        return 1
    fi
    test -z "$failed_census" || return 1
    left="$self_test_root/left"
    right="$self_test_root/right"
    /usr/bin/mkdir -m 0700 "$left" "$right"
    /usr/bin/mkdir -m 0700 "$left/normalized" "$right/normalized"
    for relative in "${normalized_evidence_files[@]}"; do
        /usr/bin/printf 'canonical fixture: %s\n' "$relative" \
            > "$left/normalized/$relative"
        /usr/bin/cp --reflink=never "$left/normalized/$relative" \
            "$right/normalized/$relative"
    done
    (
        cd "$left/normalized"
        /usr/bin/sha256sum "${normalized_evidence_files[@]}" \
            > ../canonical-SHA256SUMS
    )
    /usr/bin/cp --reflink=never "$left/canonical-SHA256SUMS" \
        "$right/canonical-SHA256SUMS"
    compare_normalized_evidence_closures "$left" "$right"
    /usr/bin/printf 'adversarial replacement\n' \
        >> "$right/normalized/left-compile_commands.template.json"
    (
        cd "$right/normalized"
        /usr/bin/sha256sum "${normalized_evidence_files[@]}" \
            > ../canonical-SHA256SUMS.replaced
    )
    /usr/bin/mv "$right/canonical-SHA256SUMS.replaced" \
        "$right/canonical-SHA256SUMS"
    if compare_normalized_evidence_closures "$left" "$right" \
            >/dev/null 2>&1; then
        /usr/bin/printf '%s\n' \
            'mutated normalized evidence and regenerated checksums were accepted' \
            >&2
        return 1
    fi
    /usr/bin/printf '%s\n' \
        'normalized evidence closure adversarial self-test passed'
}

verify_candidate_commit_binding()
{
    local build_closure_json=$1
    local core_manifest_json=$2
    local status_json=$3
    local bound_commit=
    local manifest_commit=
    local status_commit=
    bound_commit="$(/usr/bin/jq -er '.candidate.commit | strings' \
        "$build_closure_json")" || return 1
    [[ "$bound_commit" =~ ^[0-9a-f]{40}$ ]] || return 1
    manifest_commit="$(/usr/bin/jq -er '.source_commit | strings' \
        "$core_manifest_json")" || return 1
    status_commit="$(/usr/bin/jq -er '.source_commit | strings' \
        "$status_json")" || return 1
    test "$manifest_commit" = "$bound_commit" || return 1
    test "$status_commit" = "$bound_commit" || return 1
    /usr/bin/printf '%s\n' "$bound_commit"
}

candidate_commit_binding_self_test()
{
    local self_test_root=
    local expected_commit=0123456789abcdef0123456789abcdef01234567
    local tampered_commit=89abcdef0123456789abcdef0123456789abcdef
    local observed_commit=
    self_test_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v17-commit-binding-self-test.XXXXXX)" || return 1
    /usr/bin/jq -n --arg commit "$expected_commit" \
        '{candidate:{commit:$commit}}' \
        > "$self_test_root/build-closure.json" || return 1
    /usr/bin/jq -n --arg commit "$expected_commit" \
        '{source_commit:$commit}' \
        > "$self_test_root/manifest.json" || return 1
    /usr/bin/jq -n --arg commit "$expected_commit" \
        '{source_commit:$commit}' \
        > "$self_test_root/status.json" || return 1
    observed_commit="$(verify_candidate_commit_binding \
        "$self_test_root/build-closure.json" \
        "$self_test_root/manifest.json" \
        "$self_test_root/status.json")" || return 1
    test "$observed_commit" = "$expected_commit" || return 1
    /usr/bin/jq --arg commit "$tampered_commit" \
        '.source_commit = $commit' "$self_test_root/manifest.json" \
        > "$self_test_root/manifest-tampered.json" || return 1
    if verify_candidate_commit_binding \
            "$self_test_root/build-closure.json" \
            "$self_test_root/manifest-tampered.json" \
            "$self_test_root/status.json" >/dev/null 2>&1; then
        /usr/bin/printf '%s\n' \
            'mismatched candidate/core/status commits were accepted' >&2
        return 1
    fi
    /usr/bin/jq --arg commit "$tampered_commit" \
        '.source_commit = $commit' "$self_test_root/status.json" \
        > "$self_test_root/status-tampered.json" || return 1
    if verify_candidate_commit_binding \
            "$self_test_root/build-closure.json" \
            "$self_test_root/manifest.json" \
            "$self_test_root/status-tampered.json" >/dev/null 2>&1; then
        /usr/bin/printf '%s\n' \
            'mismatched candidate/core/status commits were accepted' >&2
        return 1
    fi
    /usr/bin/jq --arg commit "${expected_commit^^}" \
        '.candidate.commit = $commit' "$self_test_root/build-closure.json" \
        > "$self_test_root/build-closure-uppercase.json" || return 1
    if verify_candidate_commit_binding \
            "$self_test_root/build-closure-uppercase.json" \
            "$self_test_root/manifest.json" \
            "$self_test_root/status.json" >/dev/null 2>&1; then
        /usr/bin/printf '%s\n' \
            'non-lowercase candidate commit was accepted' >&2
        return 1
    fi
    /usr/bin/cp --reflink=never "$self_test_root/manifest.json" \
        "$self_test_root/manifest-trailing-garbage.json" || return 1
    /usr/bin/printf '%s\n' 'not-json' \
        >> "$self_test_root/manifest-trailing-garbage.json" || return 1
    if verify_candidate_commit_binding \
            "$self_test_root/build-closure.json" \
            "$self_test_root/manifest-trailing-garbage.json" \
            "$self_test_root/status.json" >/dev/null 2>&1; then
        /usr/bin/printf '%s\n' \
            'failed commit-binding jq parse with valid prefix was accepted' >&2
        return 1
    fi
    /usr/bin/printf '%s\n' \
        'candidate commit binding adversarial self-test passed'
}

compare_candidate_build_closures()
{
    local normalizer=$1
    local profile=$2
    local left=$3
    local right=$4
    local expected_build=$5
    local expected_source=$6
    local expected_commit=$7
    local evidence=$8
    test ! -e "$evidence"
    /usr/bin/mkdir -m 0700 "$evidence"
    capture_allowed_diff "$evidence/raw-tree.diff" -qr "$left" "$right"
    capture_allowed_diff "$evidence/raw-compile_commands.diff" -u \
        "$left/compile_commands.json" "$right/compile_commands.json"
    capture_allowed_diff "$evidence/raw-Makefile2.diff" -u \
        "$left/CMakeFiles/Makefile2" "$right/CMakeFiles/Makefile2"
    /usr/bin/diff -qr --exclude=compile_commands.json --exclude=Makefile2 \
        "$left" "$right" > "$evidence/required-byte-identity.diff"
    /usr/bin/python3 -I -S -B "$normalizer" compare \
        --profile "$profile" \
        --left "$left" \
        --right "$right" \
        --expected-build "$expected_build" \
        --expected-source "$expected_source" \
        --expected-commit "$expected_commit" \
        --output-dir "$evidence/normalized" \
        > "$evidence/normalizer-summary.txt" \
        2> "$evidence/normalizer-stderr.log"
    /usr/bin/diff -u \
        "$evidence/normalized/left-compile_commands.canonical.json" \
        "$evidence/normalized/right-compile_commands.canonical.json" \
        > "$evidence/canonical-compile_commands.diff"
    /usr/bin/diff -u \
        "$evidence/normalized/left-Makefile2.canonical" \
        "$evidence/normalized/right-Makefile2.canonical" \
        > "$evidence/canonical-Makefile2.diff"
    /usr/bin/sha256sum \
        "$left/compile_commands.json" \
        "$right/compile_commands.json" \
        "$left/CMakeFiles/Makefile2" \
        "$right/CMakeFiles/Makefile2" \
        > "$evidence/raw-ordering-file-SHA256SUMS"
    verify_normalized_evidence_directory "$evidence/normalized"
    (
        cd "$evidence/normalized"
        /usr/bin/sha256sum "${normalized_evidence_files[@]}" \
            > ../canonical-SHA256SUMS
        /usr/bin/sha256sum -c ../canonical-SHA256SUMS
    ) > "$evidence/canonical-sha256-check.txt"
}

verify_sealed_tree()
{
    local verified=$1
    test -d "$verified"
    test ! -L "$verified"
    require_empty_output /usr/bin/find "$verified" -type l -print -quit
    require_empty_output /usr/bin/find "$verified" \
        -type f -links +1 -print -quit
    require_empty_output /usr/bin/find "$verified" \
        -type f -perm /222 -print -quit
    require_empty_output /usr/bin/find "$verified" \
        -type d -perm /222 -print -quit
    require_empty_output /usr/bin/find "$verified" \
        ! -type d ! -type f -print -quit
}

write_tree_metadata_manifest()
{
    local metadata_root=$1
    local metadata_output="$metadata_root/TREE-METADATA.json"
    test -d "$metadata_root"
    test ! -L "$metadata_root"
    test -f "$metadata_root/SHA256SUMS"
    test ! -e "$metadata_output"
    /usr/bin/python3 -I -S -B - "$metadata_root" "$metadata_output" <<'PY'
import json
import os
import stat
import sys

root = os.path.abspath(sys.argv[1])
output = os.path.abspath(sys.argv[2])
if output != os.path.join(root, "TREE-METADATA.json"):
    raise SystemExit("tree metadata output path differs")
root_stat = os.lstat(root)
if (not stat.S_ISDIR(root_stat.st_mode) or stat.S_ISLNK(root_stat.st_mode) or
        root_stat.st_uid != os.geteuid() or root_stat.st_gid != os.getegid()):
    raise SystemExit("tree metadata root ownership/type differs")

records = []

def visit(path, relative):
    metadata = os.lstat(path)
    if stat.S_ISLNK(metadata.st_mode):
        raise SystemExit("symbolic link in canonical tree: " + relative)
    if stat.S_ISDIR(metadata.st_mode):
        kind = "directory"
    elif stat.S_ISREG(metadata.st_mode):
        kind = "file"
        if metadata.st_nlink != 1:
            raise SystemExit("non-canonical file link count: " + relative)
    else:
        raise SystemExit("special node in canonical tree: " + relative)
    if metadata.st_uid != root_stat.st_uid or metadata.st_gid != root_stat.st_gid:
        raise SystemExit("non-canonical node ownership: " + relative)
    if relative != "TREE-METADATA.json":
        records.append({
            "gid": metadata.st_gid,
            "mode": format(stat.S_IMODE(metadata.st_mode) & ~0o222, "04o"),
            "nlink": metadata.st_nlink,
            "path": relative,
            "type": kind,
            "uid": metadata.st_uid,
        })
    if kind == "directory":
        with os.scandir(path) as children:
            entries = sorted(children, key=lambda entry: entry.name)
        for entry in entries:
            child_relative = entry.name if relative == "." else \
                relative + "/" + entry.name
            visit(entry.path, child_relative)

visit(root, ".")
records.sort(key=lambda record: record["path"])
document = {
    "entries": records,
    "excluded_paths": ["TREE-METADATA.json"],
    "final_mode_policy": "observed mode with all write bits removed",
    "root": ".",
    "schema": "leopard2-authoritative-tree-metadata/v1",
    "self_policy": {
        "gid": root_stat.st_gid,
        "mode": "0400",
        "nlink": 1,
        "sha256_binding": "exactly one ./TREE-METADATA.json checksum entry",
        "type": "file",
        "uid": root_stat.st_uid,
    },
    "uid_gid_policy": {
        "gid": root_stat.st_gid,
        "rule": "every retained node has the invoking effective uid and gid",
        "uid": root_stat.st_uid,
    },
}
payload = json.dumps(
    document, sort_keys=True, separators=(",", ":"), ensure_ascii=True
).encode("ascii") + b"\n"
flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
if hasattr(os, "O_NOFOLLOW"):
    flags |= os.O_NOFOLLOW
descriptor = os.open(output, flags, 0o600)
try:
    os.fchmod(descriptor, 0o600)
    view = memoryview(payload)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise RuntimeError("short tree metadata write")
        view = view[written:]
    os.fsync(descriptor)
finally:
    os.close(descriptor)
directory = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
try:
    os.fsync(directory)
finally:
    os.close(directory)
PY
}

verify_tree_metadata_manifest()
{
    local metadata_root=$1
    /usr/bin/python3 -I -S -B - "$metadata_root" <<'PY'
import hashlib
import json
import os
import stat
import sys

root = os.path.abspath(sys.argv[1])
manifest_path = os.path.join(root, "TREE-METADATA.json")

def strict_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result

with open(manifest_path, "rb") as stream:
    raw = stream.read()
document = json.loads(raw, object_pairs_hook=strict_object)
canonical = json.dumps(
    document, sort_keys=True, separators=(",", ":"), ensure_ascii=True
).encode("ascii") + b"\n"
if raw != canonical:
    raise SystemExit("tree metadata is not canonical JSON")
if set(document) != {
        "entries", "excluded_paths", "final_mode_policy", "root", "schema",
        "self_policy", "uid_gid_policy"}:
    raise SystemExit("tree metadata keys differ")
if (document["schema"] != "leopard2-authoritative-tree-metadata/v1" or
        document["root"] != "." or
        document["excluded_paths"] != ["TREE-METADATA.json"] or
        document["final_mode_policy"] !=
        "observed mode with all write bits removed"):
    raise SystemExit("tree metadata contract differs")

root_stat = os.lstat(root)
policy = document["uid_gid_policy"]
if (policy != {
        "gid": os.getegid(),
        "rule": "every retained node has the invoking effective uid and gid",
        "uid": os.geteuid(),
        } or root_stat.st_uid != policy["uid"] or
        root_stat.st_gid != policy["gid"]):
    raise SystemExit("tree metadata uid/gid policy differs")
self_policy = document["self_policy"]
manifest_stat = os.lstat(manifest_path)
if (self_policy != {
        "gid": policy["gid"],
        "mode": "0400",
        "nlink": 1,
        "sha256_binding": "exactly one ./TREE-METADATA.json checksum entry",
        "type": "file",
        "uid": policy["uid"],
        } or not stat.S_ISREG(manifest_stat.st_mode) or
        stat.S_ISLNK(manifest_stat.st_mode) or
        manifest_stat.st_uid != self_policy["uid"] or
        manifest_stat.st_gid != self_policy["gid"] or
        format(stat.S_IMODE(manifest_stat.st_mode), "04o") !=
        self_policy["mode"] or manifest_stat.st_nlink != self_policy["nlink"]):
    raise SystemExit("tree metadata self policy differs")

actual = []

def visit(path, relative):
    metadata = os.lstat(path)
    if stat.S_ISLNK(metadata.st_mode):
        raise SystemExit("symbolic link in sealed tree: " + relative)
    if stat.S_ISDIR(metadata.st_mode):
        kind = "directory"
    elif stat.S_ISREG(metadata.st_mode):
        kind = "file"
        if metadata.st_nlink != 1:
            raise SystemExit("non-canonical file link count: " + relative)
    else:
        raise SystemExit("special node in sealed tree: " + relative)
    if metadata.st_uid != policy["uid"] or metadata.st_gid != policy["gid"]:
        raise SystemExit("sealed node ownership differs: " + relative)
    if stat.S_IMODE(metadata.st_mode) & 0o222:
        raise SystemExit("sealed node remains writable: " + relative)
    if relative != "TREE-METADATA.json":
        actual.append({
            "gid": metadata.st_gid,
            "mode": format(stat.S_IMODE(metadata.st_mode), "04o"),
            "nlink": metadata.st_nlink,
            "path": relative,
            "type": kind,
            "uid": metadata.st_uid,
        })
    if kind == "directory":
        with os.scandir(path) as children:
            entries = sorted(children, key=lambda entry: entry.name)
        for entry in entries:
            child_relative = entry.name if relative == "." else \
                relative + "/" + entry.name
            visit(entry.path, child_relative)

visit(root, ".")
actual.sort(key=lambda record: record["path"])
if type(document["entries"]) is not list or document["entries"] != actual:
    raise SystemExit("sealed tree has extra, missing, or changed nodes")

manifest_hash = hashlib.sha256(raw).hexdigest()
with open(os.path.join(root, "SHA256SUMS"), "rb") as stream:
    checksum_lines = stream.read().decode("ascii").splitlines()
expected = manifest_hash + "  ./TREE-METADATA.json"
if checksum_lines.count(expected) != 1:
    raise SystemExit("tree metadata is not uniquely SHA-bound")
PY
}

reconstruct_owner_only_campaign_tree()
{
    local replay_source=$1
    local replay_destination=$2
    /usr/bin/python3 -I -S -B - \
        "$replay_source" "$replay_destination" <<'PY'
import os
import stat
import sys

source_root = os.path.abspath(sys.argv[1])
destination_root = os.path.abspath(sys.argv[2])
if os.path.exists(destination_root):
    raise SystemExit("retained replay destination already exists")
source_root_stat = os.lstat(source_root)
if (not stat.S_ISDIR(source_root_stat.st_mode) or
        stat.S_ISLNK(source_root_stat.st_mode) or
        source_root_stat.st_uid != os.geteuid()):
    raise SystemExit("retained campaign source root is unsafe")
os.mkdir(destination_root, 0o700)
os.chmod(destination_root, 0o700, follow_symlinks=False)

def copy_directory(source, destination):
    before = os.lstat(source)
    if (not stat.S_ISDIR(before.st_mode) or stat.S_ISLNK(before.st_mode) or
            before.st_uid != os.geteuid() or
            stat.S_IMODE(before.st_mode) & 0o222):
        raise SystemExit("retained campaign directory is unsafe: " + source)
    with os.scandir(source) as stream:
        entries = sorted(stream, key=lambda entry: entry.name)
    for entry in entries:
        source_path = os.path.join(source, entry.name)
        destination_path = os.path.join(destination, entry.name)
        metadata = os.lstat(source_path)
        if stat.S_ISDIR(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
            os.mkdir(destination_path, 0o700)
            os.chmod(destination_path, 0o700, follow_symlinks=False)
            copy_directory(source_path, destination_path)
        elif stat.S_ISREG(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
            if (metadata.st_uid != os.geteuid() or metadata.st_nlink != 1 or
                    stat.S_IMODE(metadata.st_mode) & 0o222):
                raise SystemExit("retained campaign file is unsafe: " + source_path)
            source_flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
            destination_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
            if hasattr(os, "O_NOFOLLOW"):
                source_flags |= os.O_NOFOLLOW
                destination_flags |= os.O_NOFOLLOW
            source_descriptor = os.open(source_path, source_flags)
            destination_descriptor = os.open(
                destination_path, destination_flags, 0o600)
            try:
                opened = os.fstat(source_descriptor)
                identity = (
                    metadata.st_dev, metadata.st_ino, metadata.st_mode,
                    metadata.st_nlink, metadata.st_size, metadata.st_mtime_ns,
                    metadata.st_ctime_ns,
                )
                if identity != (
                        opened.st_dev, opened.st_ino, opened.st_mode,
                        opened.st_nlink, opened.st_size, opened.st_mtime_ns,
                        opened.st_ctime_ns):
                    raise RuntimeError("retained campaign file changed before copy")
                remaining = metadata.st_size
                while remaining:
                    block = os.read(source_descriptor, min(1024 * 1024, remaining))
                    if not block:
                        raise RuntimeError("short retained campaign read")
                    view = memoryview(block)
                    while view:
                        written = os.write(destination_descriptor, view)
                        if written <= 0:
                            raise RuntimeError("short retained campaign write")
                        view = view[written:]
                    remaining -= len(block)
                if os.read(source_descriptor, 1):
                    raise RuntimeError("retained campaign file grew during copy")
                final = os.fstat(source_descriptor)
                if identity != (
                        final.st_dev, final.st_ino, final.st_mode,
                        final.st_nlink, final.st_size, final.st_mtime_ns,
                        final.st_ctime_ns):
                    raise RuntimeError("retained campaign file changed during copy")
                os.fchmod(destination_descriptor, 0o600)
                os.fsync(destination_descriptor)
            finally:
                os.close(destination_descriptor)
                os.close(source_descriptor)
        else:
            raise SystemExit("special node in retained campaign: " + source_path)
    after = os.lstat(source)
    with os.scandir(source) as stream:
        after_names = sorted(entry.name for entry in stream)
    if (before.st_dev, before.st_ino, before.st_mode, before.st_nlink,
            before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_dev, after.st_ino, after.st_mode, after.st_nlink,
            after.st_mtime_ns, after.st_ctime_ns) or \
            [entry.name for entry in entries] != after_names:
        raise RuntimeError("retained campaign directory changed during copy")
    descriptor = os.open(destination, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)

copy_directory(source_root, destination_root)
PY
}

verify_envelope()
{
    local verified_envelope=$1
    local semantic_complete_status=${2:-}
    local semantic_postseal_audit=${3:-}
    local verified_core="$verified_envelope/core"
    local verified_postseal_audit="$verified_envelope/postseal-audit.json"
    local status_file=
    local status_schema=
    local evidence_generation=
    local core_schema=
    local audit_schema=
    local status_value=
    local promotion_value=
    local campaign_exit_value=
    local campaign_manifest_sha=
    local candidate_raw_tree_value=
    local candidate_test_normalization_report_sha256=
    local candidate_test_raw_tree_value=
    local verified_candidate_commit=
    local replay_campaign=
    local replay_controller_root=
    local verifier_tmp=
    local passive_file=
    local passive_field=
    local passive_hash_binding=
    local artifact_hash_binding=
    local artifact_hash_field=
    local artifact_hash_path=
    local artifact_hash_expected=
    local verified_performance_gate_passed=
    local active_required_artifact=
    local passive_only_artifact=
    local active_only_artifact=
    local failure_path=
    local retained_failure_sha=
    local failed_verify_root=
    local replay_failure_status=
    local failed_core_schema=
    local verified_timeout_argument=
    local controller_schema=
    local controller_generation=
    local observed_core_sha256sums_verified=false
    local -a audit_replay_args=()
    local -a census_verify_args=()

    test -d "$verified_envelope"
    test ! -L "$verified_envelope"
    test "$(/usr/bin/readlink -f "$verified_envelope")" = \
        "$verified_envelope"
    test -d "$verified_core"
    test ! -L "$verified_core"
    test "$(/usr/bin/readlink -f "$verified_core")" = "$verified_core"
    test -f "$verified_envelope/SHA256SUMS"
    test -f "$verified_envelope/TREE-METADATA.json"
    test -f "$verified_core/SHA256SUMS"
    test -f "$verified_core/run-authoritative.sh"
    (
        cd "$verified_envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null
    verify_sealed_tree "$verified_envelope"
    verify_tree_metadata_manifest "$verified_envelope"

    if [[ -n "$semantic_complete_status" || \
          -n "$semantic_postseal_audit" ]]; then
        test -n "$semantic_complete_status"
        test -n "$semantic_postseal_audit"
        test "${semantic_complete_status##*/}" = NOT_PROMOTED.json
        test "$semantic_complete_status" = \
            "$(/usr/bin/readlink -f "$semantic_complete_status")"
        test "$semantic_postseal_audit" = \
            "$(/usr/bin/readlink -f "$semantic_postseal_audit")"
        test -f "$semantic_complete_status"
        test ! -L "$semantic_complete_status"
        test -f "$semantic_postseal_audit"
        test ! -L "$semantic_postseal_audit"
        test -f "$verified_envelope/FAILED.json"
        test ! -e "$verified_envelope/COMPLETED.json"
        test ! -e "$verified_envelope/NOT_PROMOTED.json"
        validate_v18_terminal_static \
            "$verified_envelope/FAILED.json" failed-v2
        /usr/bin/jq -e '
            keys == (["acquisition_generation","attempt","attempt_budget",
                "attempt_lineage_sha256","baseline_commit",
                "campaign_exit_status","core_sha256sums_sha256",
                "core_sha256sums_verified","promotion_passed","schema",
                "source_commit","source_tree","stage","status"] | sort) and
            .schema == "leopard2-v18-gfni-main-failed-envelope/v1" and
            .status == "failed" and
            .acquisition_generation == "passive-v2" and
            .core_sha256sums_verified == true and
            (.stage as $stage |
                ["seal_core", "independent_postseal_audit",
                 "publish_not_promoted_envelope"] |
                index($stage) != null)
        ' "$verified_envelope/FAILED.json" >/dev/null
        /usr/bin/jq -e \
            --slurpfile failed "$verified_envelope/FAILED.json" '
            ($failed | length) == 1 and
            .attempt == $failed[0].attempt and
            .attempt_budget == $failed[0].attempt_budget and
            .attempt_lineage_sha256 ==
                $failed[0].attempt_lineage_sha256 and
            .source_commit == $failed[0].source_commit and
            .source_tree == $failed[0].source_tree and
            .baseline_commit == $failed[0].baseline_commit and
            .core_sha256sums_sha256 ==
                $failed[0].core_sha256sums_sha256
        ' "$semantic_complete_status" >/dev/null
        status_file="$semantic_complete_status"
        status_value=complete
        promotion_value=false
        verified_postseal_audit="$semantic_postseal_audit"
    elif [[ -f "$verified_envelope/COMPLETED.json" ]]; then
        test ! -e "$verified_envelope/NOT_PROMOTED.json"
        test ! -e "$verified_envelope/FAILED.json"
        status_file="$verified_envelope/COMPLETED.json"
        status_value=complete
        promotion_value=true
    elif [[ -f "$verified_envelope/NOT_PROMOTED.json" ]]; then
        test ! -e "$verified_envelope/COMPLETED.json"
        test ! -e "$verified_envelope/FAILED.json"
        status_file="$verified_envelope/NOT_PROMOTED.json"
        status_value=complete
        promotion_value=false
    else
        test -f "$verified_envelope/FAILED.json"
        test ! -e "$verified_envelope/COMPLETED.json"
        test ! -e "$verified_envelope/NOT_PROMOTED.json"
        status_file="$verified_envelope/FAILED.json"
        status_value=failed
        promotion_value=false
    fi
    status_schema="$(/usr/bin/jq -er '.schema | strings' "$status_file")"
    case "${status_file##*/}:$status_schema" in
        COMPLETED.json:leopard2-v17-gfni-main-completion-envelope/v1)
            evidence_generation=active-v1
            core_schema=leopard2-v17-gfni-main-core-manifest/v3
            audit_schema=leopard2-main-compare-v17-independent-audit/v1
            ;;
        NOT_PROMOTED.json:leopard2-v17-gfni-main-not-promoted-envelope/v1)
            evidence_generation=active-v1
            core_schema=leopard2-v17-gfni-main-core-manifest/v3
            audit_schema=leopard2-main-compare-v17-independent-audit/v1
            ;;
        NOT_PROMOTED.json:leopard2-v17-gfni-main-passive-not-promoted-envelope/v2)
            evidence_generation=passive-v1
            core_schema=leopard2-v17-gfni-main-passive-core-manifest/v4
            audit_schema=leopard2-main-compare-v17-passive-independent-audit/v1
            ;;
        NOT_PROMOTED.json:leopard2-v18-gfni-main-passive-not-promoted-envelope/v1)
            evidence_generation=passive-v2
            core_schema=leopard2-v18-gfni-main-passive-core-manifest/v1
            audit_schema=leopard2-main-compare-v18-passive-independent-audit/v1
            ;;
        FAILED.json:leopard2-v17-gfni-main-failed-envelope/v1)
            evidence_generation=failed-v1
            ;;
        FAILED.json:leopard2-v18-gfni-main-failed-envelope/v1)
            evidence_generation=failed-v2
            ;;
        *)
            return 1
            ;;
    esac
    if [[ "$evidence_generation" == passive-v2 ||
          "$evidence_generation" == failed-v2 ]]; then
        validate_v18_terminal_static "$status_file" "$evidence_generation"
    fi
    if (
        cd "$verified_core"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null 2>&1; then
        observed_core_sha256sums_verified=true
    fi
    if [[ "$evidence_generation" == failed-* ]] && \
       /usr/bin/jq -e 'has("core_sha256sums_verified")' \
           "$status_file" >/dev/null; then
        test "$observed_core_sha256sums_verified" = \
            "$(/usr/bin/jq -er \
                '.core_sha256sums_verified | booleans | tostring' \
                "$status_file")"
    else
        test "$observed_core_sha256sums_verified" = true
    fi
    test "$(/usr/bin/jq -er '.status | strings' "$status_file")" = \
        "$status_value"
    test "$(/usr/bin/jq -er \
        '.promotion_passed | booleans | tostring' "$status_file")" = \
        "$promotion_value"
    if [[ "$evidence_generation" == passive-v1 ]]; then
        /usr/bin/jq -e '
            keys == ([
                "acquisition_generation", "audit_sha256", "baseline_commit",
                "benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies",
                "campaign_exit_status", "campaign_manifest_sha256",
                "causal_performance_claim_eligible", "core_manifest_sha256",
                "core_sha256sums_sha256", "cpu_pair_exclusive",
                "evidence_class", "evidence_valid", "isolation_claim",
                "performance_gate_passed", "postseal_audit_byte_identical",
                "postseal_audit_passed", "postseal_audit_sha256",
                "preseal_audit_passed", "promotion_eligible",
                "promotion_passed", "schema", "source_commit", "source_tree",
                "status"
            ] | sort) and
            .acquisition_generation == "passive-v1" and
            .evidence_class == "passive-shared-host-observation/v1" and
            .promotion_eligible == false and .promotion_passed == false and
            .causal_performance_claim_eligible == false and
            .cpu_pair_exclusive == false and .evidence_valid == true and
            .preseal_audit_passed == true and
            .postseal_audit_passed == true and
            .postseal_audit_byte_identical == true and
            (.performance_gate_passed | type == "boolean") and
            (.benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies |
                type) == "number" and
            (.benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies |
                floor) ==
                .benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies and
            .benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies >= 0
        ' "$status_file" >/dev/null
    elif [[ "$evidence_generation" == passive-v2 ]]; then
        /usr/bin/jq -e '
            keys == ([
                "acquisition_generation", "attempt", "attempt_budget",
                "attempt_lineage_sha256",
                "audit_sha256", "baseline_commit", "campaign_exit_status",
                "campaign_manifest_sha256",
                "causal_performance_claim_eligible", "core_manifest_sha256",
                "core_sha256sums_sha256", "cpu_pair_exclusive",
                "evidence_class", "evidence_valid", "isolation_claim",
                "performance_gate_passed", "postseal_audit_byte_identical",
                "postseal_audit_passed", "postseal_audit_sha256",
                "preseal_audit_passed", "promotion_eligible",
                "promotion_passed", "schema", "source_commit", "source_tree",
                "out_of_window_benchmark_cpu_nonidle_jiffies",
                "out_of_window_reserved_sibling_nonidle_jiffies", "status",
                "windowed_benchmark_cpu_nonidle_excess_jiffies",
                "windowed_reserved_sibling_nonidle_jiffies"
            ] | sort) and
            .acquisition_generation == "passive-v2" and
            (.attempt | type) == "number" and .attempt == (.attempt | floor) and
            .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
            .evidence_class == "passive-windowed-shared-host-observation/v1" and
            .promotion_eligible == false and .promotion_passed == false and
            .causal_performance_claim_eligible == false and
            .cpu_pair_exclusive == false and .evidence_valid == true and
            .preseal_audit_passed == true and
            .postseal_audit_passed == true and
            .postseal_audit_byte_identical == true and
            (.performance_gate_passed | type) == "boolean" and
            .windowed_benchmark_cpu_nonidle_excess_jiffies == 0 and
            .windowed_reserved_sibling_nonidle_jiffies == 0 and
            (.out_of_window_benchmark_cpu_nonidle_jiffies | type) == "number" and
            .out_of_window_benchmark_cpu_nonidle_jiffies ==
                (.out_of_window_benchmark_cpu_nonidle_jiffies | floor) and
            .out_of_window_benchmark_cpu_nonidle_jiffies >= 0 and
            (.out_of_window_reserved_sibling_nonidle_jiffies | type) == "number" and
            .out_of_window_reserved_sibling_nonidle_jiffies ==
                (.out_of_window_reserved_sibling_nonidle_jiffies | floor) and
            .out_of_window_reserved_sibling_nonidle_jiffies >= 0
        ' "$status_file" >/dev/null
        test "$verified_envelope" = \
            "$repo/.research/leopard-79h/$(/usr/bin/jq -er '.source_commit | strings' "$status_file" | /usr/bin/cut -c1-7)-v18-passive-main-a$(/usr/bin/jq -er '.attempt | numbers' "$status_file")"
    fi
    if [[ "$evidence_generation" == failed-v2 ]]; then
        /usr/bin/jq -e '
            .acquisition_generation == "passive-v2" and
            (.attempt | type) == "number" and
            .attempt == (.attempt | floor) and
            .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3
        ' "$status_file" >/dev/null
        test "$verified_envelope" = \
            "$repo/.research/leopard-79h/$(/usr/bin/jq -er '.source_commit | strings' "$status_file" | /usr/bin/cut -c1-7)-v18-passive-main-a$(/usr/bin/jq -er '.attempt | numbers' "$status_file")"
    fi
    if [[ "$evidence_generation" == passive-v2 ||
          "$evidence_generation" == failed-v2 ]]; then
        validate_v18_attempt_lineage \
            "$verified_envelope" "$verified_core" "$status_file"
    fi
    campaign_exit_value="$(/usr/bin/jq -er \
        '.campaign_exit_status | numbers' "$status_file")"
    test "$(/usr/bin/sha256sum "$verified_core/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -er '.core_sha256sums_sha256 | strings' \
            "$status_file")"

    if [[ "$status_value" == complete ]]; then
        test "$campaign_exit_value" -eq 0
        verifier_tmp="$(/usr/bin/mktemp -d \
            /tmp/leopard-v17-gfni-main-verify.XXXXXX)"
        test -f "$verified_core/manifest.json"
        test -f "$verified_core/campaign/manifest.json"
        test -f "$verified_core/audit.json"
        test -f "$verified_core/result-summary.json"
        if [[ "$evidence_generation" == active-v1 ]]; then
            for active_required_artifact in \
                affinity-report.json affinity-binding.json \
                affinity-report.json.accepted.json \
                affinity-report-verification.log \
                affinity-binding-create.log \
                affinity-binding-verification.log supervisor-self-test.log; do
                test -f "$verified_core/$active_required_artifact"
            done
            test ! -e "$verified_core/affinity-report.json.ambiguous"
            for passive_only_artifact in \
                passive_environment_census.py \
                environment-census-pre.json environment-census-post.json \
                passive-environment-policy.json controller-affinity.json \
                passive-census-self-test.log passive-auditor-self-test.log \
                passive-controller-taskset.log \
                attempt-lineage.json \
                build-closure/committed-passive-census.py; do
                test ! -e "$verified_core/$passive_only_artifact"
            done
        else
            test "$evidence_generation" = passive-v1 ||
                test "$evidence_generation" = passive-v2
            for active_only_artifact in \
                affinity-report.json affinity-binding.json \
                affinity-report.json.accepted.json \
                affinity-report.json.ambiguous \
                affinity-report-verification.log \
                affinity-binding-create.log \
                affinity-binding-verification.log supervisor-self-test.log \
                passive-auditor-self-test.log \
                passive-controller-taskset.log; do
                test ! -e "$verified_core/$active_only_artifact"
            done
            test -f "$verified_core/passive_environment_census.py"
            test -f "$verified_core/environment-census-pre.json"
            test -f "$verified_core/environment-census-post.json"
            test -f "$verified_core/passive-environment-policy.json"
            test -f "$verified_core/controller-affinity.json"
            test -f "$verified_core/passive-census-self-test.log"
            test -f "$verified_core/campaign-command.json"
            test -f "$verified_core/wrapper-launch-affinity.json"
            if [[ "$evidence_generation" == passive-v2 ]]; then
                test -f "$verified_core/attempt-lineage.json"
            else
                test ! -e "$verified_core/attempt-lineage.json"
            fi
            test -f \
                "$verified_core/build-closure/committed-passive-census.py"
            /usr/bin/cmp "$verified_core/passive_environment_census.py" \
                "$verified_core/build-closure/committed-passive-census.py"
        fi
        test -f "$verified_core/run_abba.py"
        test -f "$verified_core/git_capture.py"
        test -f "$verified_core/balanced_evidence_common.py"
        test -f "$verified_core/leopard2_build_provenance.py"
        test -f "$verified_core/build-order-normalizer.py"
        test -f "$verified_core/build-order-normalizer-reconstruction.json"
        test -f "$verified_core/leopard2_affinity_supervisor.py"
        test -f "$verified_core/sse2neon-source.tar"
        test -f "$verified_core/build-closure.json"
        test -f "$verified_core/candidate-test-temporal-closure.json"
        test -f \
            "$verified_core/build-closure/candidate-test-selected-SHA256SUMS"
        test -f \
            "$verified_core/build-closure/candidate-test-order-normalization/normalized/report.json"
        test -f \
            "$verified_core/build-closure/candidate-test-order-normalization/canonical-SHA256SUMS"
        test -f \
            "$verified_core/build-closure/candidate-timing-order-normalization/canonical-SHA256SUMS"
        test -f \
            "$verified_core/build-closure/live-candidate-tests/leopard2_backend_failures_test"
        test -f \
            "$verified_core/build-closure/replay-candidate-tests/leopard2_backend_failures_test"
        test -f \
            "$verified_core/build-closure/live-candidate-tests/bench_leopard2_prevalidated_batch"
        test -f \
            "$verified_core/build-closure/replay-candidate-tests/bench_leopard2_prevalidated_batch"
        test -f "$verified_postseal_audit"
        test -f "$verified_core/audit_v17_gfni_main_compare.py"
        /usr/bin/cmp "$verified_core/audit.json" \
            "$verified_postseal_audit"
        /usr/bin/cmp "$verified_core/run-authoritative.sh" \
            "$verified_core/build-closure/committed-wrapper.sh"
        /usr/bin/cmp "$verified_core/audit_v17_gfni_main_compare.py" \
            "$verified_core/build-closure/committed-auditor.py"
        /usr/bin/cmp "$verified_core/leopard2_build_provenance.py" \
            "$verified_core/build-closure/committed-build-provenance.py"
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/build-closure.json")" = \
            leopard2-v17-gfni-main-build-closure/v3
        candidate_commit_binding_self_test \
            > "$verifier_tmp/candidate-commit-binding-self-test.log" 2>&1
        verified_candidate_commit="$(verify_candidate_commit_binding \
            "$verified_core/build-closure.json" \
            "$verified_core/manifest.json" "$status_file")" || return 1
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-order-normalizer.py" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.controllers.build_order_normalizer_sha256 | strings' \
                "$verified_core/build-closure.json")"
        install_build_order_normalizer \
            "$verifier_tmp/reconstructed-build-order-normalizer.py"
        /usr/bin/cmp "$verified_core/build-order-normalizer.py" \
            "$verifier_tmp/reconstructed-build-order-normalizer.py"
        /usr/bin/jq -e '
            .schema ==
                "leopard2-v17-build-order-normalizer-reconstruction/v3" and
            .embedded_wrapper_reconstruction_byte_identical == true and
            .self_test_passed == true and
            .normalized_evidence_closure_self_test_passed == true and
            .candidate_commit_binding_self_test_passed == true and
            .timing_performed == false and
            .profiles == ["candidate-tests", "candidate-timing"] and
            .profile_sha256s == {
                "candidate-tests":"bead85b797c0e9966d4f09ba5b18b48739a350f550f831d1afa5bd2ae68dfa69",
                "candidate-timing":"d2a07561d2e5919051c2b7d7464f6aa844752e42e92099190feeea004f2716b3"
            } and
            .closed_world_profile_templates == true and
            .normalized_evidence_exact_file_count == 11 and
            .template_path_tokens ==
                {build:"${BUILD}",commit:"${COMMIT}",source:"${SOURCE}"}
        ' "$verified_core/build-order-normalizer-reconstruction.json" >/dev/null
        /usr/bin/python3 -I -S -B \
            "$verified_core/build-order-normalizer.py" self-test \
            > "$verifier_tmp/build-order-normalizer-self-test.log" \
            2> "$verifier_tmp/build-order-normalizer-self-test-stderr.log"
        normalized_evidence_closure_self_test \
            > "$verifier_tmp/normalized-evidence-closure-self-test.log" 2>&1
        compare_candidate_build_closures \
            "$verified_core/build-order-normalizer.py" candidate-timing \
            "$verified_core/build-closure/live-candidate" \
            "$verified_core/build-closure/replay-candidate" \
            "$(/usr/bin/jq -er '.candidate.build | strings' \
                "$verified_core/build-closure.json")" \
            "$(/usr/bin/jq -er '.candidate.source | strings' \
                "$verified_core/build-closure.json")" \
            "$verified_candidate_commit" \
            "$verifier_tmp/candidate-timing-order-normalization"
        compare_normalized_evidence_closures \
            "$verified_core/build-closure/candidate-timing-order-normalization" \
            "$verifier_tmp/candidate-timing-order-normalization"
        /usr/bin/cmp \
            "$verified_core/build-closure/candidate-timing-order-normalization/normalized/report.json" \
            "$verifier_tmp/candidate-timing-order-normalization/normalized/report.json"
        /usr/bin/jq -e '
            .schema ==
                "leopard2-v17-candidate-build-order-normalization/v3" and
            .contract.operation == "compare" and
            .contract.profile == "candidate-timing" and
            .contract.template_path_tokens ==
                {build:"${BUILD}",commit:"${COMMIT}",source:"${SOURCE}"} and
            .makefile2.left.all_multi_prerequisite_targets_allowlisted == true and
            .makefile2.right.all_multi_prerequisite_targets_allowlisted == true and
            .makefile2.left.multi_prerequisite_target_count == 4 and
            .makefile2.right.multi_prerequisite_target_count == 4 and
            .semantic_equal == true
        ' "$verifier_tmp/candidate-timing-order-normalization/normalized/report.json" \
            >/dev/null
        candidate_raw_tree_value="$(/usr/bin/jq -er \
            '.contract.raw_tree_file_bytes_identical | booleans | tostring' \
            "$verifier_tmp/candidate-timing-order-normalization/normalized/report.json")"
        test "$(/usr/bin/jq -er '.expected_commit | strings' \
            "$verifier_tmp/candidate-timing-order-normalization/normalized/report.json")" = \
            "$verified_candidate_commit"
        test "$candidate_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.candidate.reproduction.raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/build-closure.json")"
        test "$candidate_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.build_reproduction.candidate_timing_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/build-closure.json")"
        test "$candidate_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.build_reproduction.candidate_timing_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/manifest.json")"
        if [[ "$candidate_raw_tree_value" == true ]]; then
            test "$(/usr/bin/cat \
                "$verifier_tmp/candidate-timing-order-normalization/raw-tree.diff.status")" = 0
        else
            test "$(/usr/bin/cat \
                "$verifier_tmp/candidate-timing-order-normalization/raw-tree.diff.status")" = 1
        fi
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-closure/candidate-timing-order-normalization/normalized/report.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.candidate.reproduction.normalization_report_sha256 | strings' \
                "$verified_core/build-closure.json")"
        compare_candidate_build_closures \
            "$verified_core/build-order-normalizer.py" candidate-tests \
            "$verified_core/build-closure/live-candidate-tests" \
            "$verified_core/build-closure/replay-candidate-tests" \
            "$(/usr/bin/jq -er '.candidate_tests.build | strings' \
                "$verified_core/build-closure.json")" \
            "$(/usr/bin/jq -er '.candidate.source | strings' \
                "$verified_core/build-closure.json")" \
            "$verified_candidate_commit" \
            "$verifier_tmp/candidate-test-order-normalization"
        compare_normalized_evidence_closures \
            "$verified_core/build-closure/candidate-test-order-normalization" \
            "$verifier_tmp/candidate-test-order-normalization"
        /usr/bin/cmp \
            "$verified_core/build-closure/candidate-test-order-normalization/normalized/report.json" \
            "$verifier_tmp/candidate-test-order-normalization/normalized/report.json"
        /usr/bin/jq -e '
            .schema ==
                "leopard2-v17-candidate-build-order-normalization/v3" and
            .contract.operation == "compare" and
            .contract.profile == "candidate-tests" and
            .contract.template_path_tokens ==
                {build:"${BUILD}",commit:"${COMMIT}",source:"${SOURCE}"} and
            .makefile2.left.all_multi_prerequisite_targets_allowlisted == true and
            .makefile2.right.all_multi_prerequisite_targets_allowlisted == true and
            .makefile2.left.multi_prerequisite_target_count == 6 and
            .makefile2.right.multi_prerequisite_target_count == 6 and
            .semantic_equal == true
        ' "$verifier_tmp/candidate-test-order-normalization/normalized/report.json" \
            >/dev/null
        candidate_test_raw_tree_value="$(/usr/bin/jq -er \
            '.contract.raw_tree_file_bytes_identical | booleans | tostring' \
            "$verifier_tmp/candidate-test-order-normalization/normalized/report.json")"
        test "$(/usr/bin/jq -er '.expected_commit | strings' \
            "$verifier_tmp/candidate-test-order-normalization/normalized/report.json")" = \
            "$verified_candidate_commit"
        test "$candidate_test_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.reproduction.raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/build-closure.json")"
        test "$candidate_test_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.two_clean_builds_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/build-closure.json")"
        test "$candidate_test_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.build_reproduction.candidate_tests_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/build-closure.json")"
        test "$candidate_test_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.two_clean_builds_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/candidate-test-temporal-closure.json")"
        test "$candidate_test_raw_tree_value" = \
            "$(/usr/bin/jq -er \
                '.build_reproduction.candidate_tests_raw_tree_file_bytes_identical | booleans | tostring' \
                "$verified_core/manifest.json")"
        if [[ "$candidate_test_raw_tree_value" == true ]]; then
            test "$(/usr/bin/cat \
                "$verifier_tmp/candidate-test-order-normalization/raw-tree.diff.status")" = 0
        else
            test "$(/usr/bin/cat \
                "$verifier_tmp/candidate-test-order-normalization/raw-tree.diff.status")" = 1
        fi
        candidate_test_normalization_report_sha256="$(/usr/bin/sha256sum \
            "$verified_core/build-closure/candidate-test-order-normalization/normalized/report.json" | \
            /usr/bin/cut -d' ' -f1)"
        test "$candidate_test_normalization_report_sha256" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.reproduction.normalization_report_sha256 | strings' \
                "$verified_core/build-closure.json")"
        test "$candidate_test_normalization_report_sha256" = \
            "$(/usr/bin/jq -er '.normalization_report_sha256 | strings' \
                "$verified_core/candidate-test-temporal-closure.json")"
        test "$candidate_test_normalization_report_sha256" = \
            "$(/usr/bin/jq -er \
                '.candidate_test_normalization_report_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/jq -cer '[
                .profile_contract.profile_sha256,
                .profile_contract.node_census_count,
                .profile_contract.node_census_sha256,
                .profile_contract.compile_template_sha256,
                .profile_contract.compile_output_list_sha256,
                .profile_contract.make_template_sha256
            ]' \
            "$verified_core/build-closure/candidate-timing-order-normalization/normalized/report.json")" = \
            "$(/usr/bin/jq -cer '[
                .candidate.reproduction.profile_sha256,
                .candidate.reproduction.node_census_count,
                .candidate.reproduction.node_census_sha256,
                .candidate.reproduction.closed_world_template_sha256s.compile_commands,
                .candidate.reproduction.closed_world_template_sha256s.compile_output_list,
                .candidate.reproduction.closed_world_template_sha256s.makefile2
            ]' "$verified_core/build-closure.json")"
        test "$(/usr/bin/jq -cer '[
                .profile_contract.profile_sha256,
                .profile_contract.node_census_count,
                .profile_contract.node_census_sha256,
                .profile_contract.compile_template_sha256,
                .profile_contract.compile_output_list_sha256,
                .profile_contract.make_template_sha256
            ]' \
            "$verified_core/build-closure/candidate-test-order-normalization/normalized/report.json")" = \
            "$(/usr/bin/jq -cer '[
                .candidate_tests.reproduction.profile_sha256,
                .candidate_tests.reproduction.node_census_count,
                .candidate_tests.reproduction.node_census_sha256,
                .candidate_tests.reproduction.closed_world_template_sha256s.compile_commands,
                .candidate_tests.reproduction.closed_world_template_sha256s.compile_output_list,
                .candidate_tests.reproduction.closed_world_template_sha256s.makefile2
            ]' "$verified_core/build-closure.json")"
        /usr/bin/jq -e '
            .candidate.reproduction.qualified_semantic_equal == true and
            .candidate.reproduction.normalizer_profile == "candidate-timing" and
            .candidate.reproduction.all_nonexception_file_bytes_identical == true and
            .candidate.reproduction.order_normalized_exception_paths ==
                ["compile_commands.json", "CMakeFiles/Makefile2"] and
            .candidate.reproduction.compile_commands_exact_entry_count == 30 and
            .candidate.reproduction.compile_commands_compiler_counts ==
                {"/usr/bin/c++":30} and
            .candidate.reproduction.makefile2_exact_validated_block_count == 4 and
            .candidate.reproduction.makefile2_exact_normalized_block_count == 4 and
            .candidate.reproduction.closed_world_template_sha256s == {
                compile_commands:"6ad3a97ffe31f4eaf5c8f2ff5459dd310635e6e430181c460d9803f11cc5fb28",
                compile_output_list:"a39c7e87cc0a506148faaa4023b11f96c2e26b22c103584a62b62ae91f2387b7",
                makefile2:"5ab77166d931546b1cb998cfdf4eb5c7dd1850994b2eef65d8c5e627a3899a62"
            } and
            .candidate.reproduction.node_census_count == 131 and
            .candidate.reproduction.node_census_sha256 ==
                "d3194ea1135fbdc0067cd3af72693cbf27d1e28a4734804b52bf5be1a3f37aca" and
            .candidate.reproduction.profile_sha256 ==
                "d2a07561d2e5919051c2b7d7464f6aa844752e42e92099190feeea004f2716b3" and
            .candidate.reproduction.normalized_evidence_exact_file_count == 11 and
            .candidate.reproduction.normalized_evidence_relative_canonical_sha256sums == true and
            .candidate_tests.reproduction.qualified_semantic_equal == true and
            .candidate_tests.reproduction.normalizer_profile == "candidate-tests" and
            .candidate_tests.reproduction.all_nonexception_file_bytes_identical == true and
            .candidate_tests.reproduction.order_normalized_exception_paths ==
                ["compile_commands.json", "CMakeFiles/Makefile2"] and
            .candidate_tests.reproduction.compile_commands_exact_entry_count == 170 and
            .candidate_tests.reproduction.compile_commands_compiler_counts ==
                {"/usr/bin/c++":169,"/usr/bin/cc":1} and
            .candidate_tests.reproduction.single_c_record == {
                source:"tests/leopard2/test_codec_options_abi.c",
                output:"CMakeFiles/leopard2_codec_options_abi_test.dir/tests/leopard2/test_codec_options_abi.c.o"
            } and
            .candidate_tests.reproduction.makefile2_exact_validated_block_count == 6 and
            .candidate_tests.reproduction.makefile2_exact_normalized_block_count == 6 and
            .candidate_tests.reproduction.closed_world_template_sha256s == {
                compile_commands:"dda169bcd458db8b17c09c1e8873c9342731a4c5639874849ad87d3c1003badb",
                compile_output_list:"8aeea7120bbbe5e3f5e71059b29007b4fb3e01b153c55f0aa874c45de9e48574",
                makefile2:"c56eb276850e7acbc5242e91c26b8fecf993040feeeb88b56596a189924a8479"
            } and
            .candidate_tests.reproduction.node_census_count == 556 and
            .candidate_tests.reproduction.node_census_sha256 ==
                "2fa17fb9bb1458ea2bb7c361de077323e6f0a55359ba137eb445d231dd13d24b" and
            .candidate_tests.reproduction.profile_sha256 ==
                "bead85b797c0e9966d4f09ba5b18b48739a350f550f831d1afa5bd2ae68dfa69" and
            .candidate_tests.reproduction.normalized_evidence_exact_file_count == 11 and
            .candidate_tests.reproduction.normalized_evidence_relative_canonical_sha256sums == true and
            .candidate_tests.two_clean_builds_qualified_semantic_equal == true and
            .build_reproduction.candidate_timing_qualified_semantic_equal == true and
            .build_reproduction.candidate_tests_qualified_semantic_equal == true and
            .build_reproduction.baseline_raw_byte_identical == true and
            .baseline.two_clean_builds_raw_byte_identical == true and
            .build_reproduction.identical_absolute_source_and_build_paths == true and
            .build_reproduction.objects_archives_binaries_link_recipes_cache_and_nonordering_generated_inputs_raw_byte_identical == true and
            .controllers.build_order_normalizer_origin ==
                "embedded deterministic helper from committed wrapper" and
            .audit_boundary ==
                "independent campaign auditor covers retained campaign semantics; build-order qualification is separately recomputed by the frozen wrapper normalizer during durable verification"
        ' "$verified_core/build-closure.json" >/dev/null
        /usr/bin/cmp \
            "$verified_core/build-closure/live-candidate-tests/leopard2_backend_failures_test" \
            "$verified_core/build-closure/replay-candidate-tests/leopard2_backend_failures_test"
        /usr/bin/cmp \
            "$verified_core/build-closure/live-candidate-tests/bench_leopard2_prevalidated_batch" \
            "$verified_core/build-closure/replay-candidate-tests/bench_leopard2_prevalidated_batch"
        /usr/bin/diff -qr \
            --exclude=compile_commands.json --exclude=Makefile2 \
            "$verified_core/build-closure/live-candidate-tests" \
            "$verified_core/build-closure/replay-candidate-tests" >/dev/null
        (
            cd "$verified_core/build-closure/live-candidate-tests"
            /usr/bin/sha256sum -c \
                "$verified_core/build-closure/candidate-test-selected-SHA256SUMS"
        ) >/dev/null
        (
            cd "$verified_core/build-closure/replay-candidate-tests"
            /usr/bin/sha256sum -c \
                "$verified_core/build-closure/candidate-test-selected-SHA256SUMS"
        ) >/dev/null
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-closure/candidate-test-selected-SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.selected_sha256sums_sha256 | strings' \
                "$verified_core/build-closure.json")"
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-closure/replay-candidate-tests/leopard2_backend_failures_test" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.backend_failures_test_sha256 | strings' \
                "$verified_core/build-closure.json")"
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-closure/replay-candidate-tests/bench_leopard2_prevalidated_batch" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.candidate_tests.prevalidated_batch_test_sha256 | strings' \
                "$verified_core/build-closure.json")"
        test "$(/usr/bin/jq -er \
            '.candidate_tests.two_clean_builds_qualified_semantic_equal | booleans | tostring' \
            "$verified_core/build-closure.json")" = true
        test "$(/usr/bin/jq -er \
            '.candidate_tests.complete_object_link_cache_closure | booleans | tostring' \
            "$verified_core/build-closure.json")" = true
        /usr/bin/jq -e '
            .candidate_tests.selected_files == [
                "bench_leopard2",
                "bench_leopard2_prevalidated_batch",
                "leopard2_auto_gf16_gfni_production_test",
                "leopard2_backend_failures_test",
                "leopard2_legacy_golden_test",
                "libleopard.a",
                "libleopard_test_hooks.a",
                "generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h",
                "generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt"
            ] and
            .candidate_tests.posttest_raw_byte_identical == true and
            .candidate_tests.postcampaign_revalidation_required == true
        ' "$verified_core/build-closure.json" >/dev/null
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/candidate-test-temporal-closure.json")" = \
            leopard2-v17-gfni-main-candidate-test-temporal-closure/v2
        test "$(/usr/bin/jq -er \
            '.selected_sha256sums_sha256 | strings' \
            "$verified_core/candidate-test-temporal-closure.json")" = \
            "$(/usr/bin/jq -er \
            '.candidate_tests.selected_sha256sums_sha256 | strings' \
            "$verified_core/build-closure.json")"
        /usr/bin/jq -e '
            .two_clean_builds_qualified_semantic_equal == true and
            (.two_clean_builds_raw_tree_file_bytes_identical |
                type == "boolean") and
            .posttest_byte_identical == true and
            .postcampaign_byte_identical == true and
            .canonical_test_build_frozen_during_campaign == true
        ' "$verified_core/candidate-test-temporal-closure.json" >/dev/null
        test "$(/usr/bin/jq -er '.evidence_valid | booleans | tostring' \
            "$status_file")" = true
        if [[ "$evidence_generation" == active-v1 ]]; then
            test "$(/usr/bin/jq -er \
                '.performance_gate_passed | booleans | tostring' \
                "$status_file")" = "$promotion_value"
        else
            /usr/bin/jq -e \
                '.performance_gate_passed | type == "boolean"' \
                "$status_file" >/dev/null
        fi
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/manifest.json")" = "$core_schema"
        test "$(/usr/bin/jq -er '.campaign_exit_status | numbers' \
            "$verified_core/manifest.json")" -eq 0
        test "$(/usr/bin/jq -er '.evidence_valid | booleans | tostring' \
            "$verified_core/manifest.json")" = true
        if [[ "$evidence_generation" == active-v1 ]]; then
            test "$(/usr/bin/jq -er \
                '.promotion_requires_completion_envelope | booleans | tostring' \
                "$verified_core/manifest.json")" = true
        elif [[ "$evidence_generation" == passive-v1 ]]; then
            /usr/bin/jq -e '
                keys == ([
                    "acquisition_generation", "active_affinity_supervisor_executed",
                    "audit_sha256", "auditor_sha256", "baseline_archive_sha256_pre_post",
                    "baseline_binary_sha256_pre_post", "baseline_commit",
                    "benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies",
                    "build_order_normalizer_sha256", "build_reproduction",
                    "campaign_exit_status", "campaign_manifest_sha256",
                    "campaign_raw_sha256", "candidate_archive_sha256_pre_post",
                    "candidate_binary_sha256_pre_post",
                    "candidate_test_normalization_report_sha256",
                    "candidate_timing_normalization_report_sha256",
                    "canonical_lock", "causal_performance_claim_eligible",
                    "controller_affinity_sha256", "cpu", "cpu_pair_exclusive",
                    "environment_census_post_sha256",
                    "environment_census_pre_sha256", "evidence_class",
                    "evidence_valid", "independent_auditor_scope",
                    "independent_auditor_supervision_mode",
                    "independent_preseal_audit_passed", "isolation_claim",
                    "passive_census_sha256", "passive_environment_policy_sha256",
                    "performance_gate_passed", "postseal_policy",
                    "producer_verification_mode", "producer_verification_passed",
                    "promotion_eligible", "promotion_passed", "ratio_policy",
                    "runner_sha256", "schema", "sibling", "source_commit",
                    "source_tree", "sse2neon_commit",
                    "sse2neon_source_archive_sha256", "status",
                    "supervisor_role", "supervisor_sha256", "wrapper_sha256"
                ] | sort) and
                .acquisition_generation == "passive-v1" and
                .evidence_class == "passive-shared-host-observation/v1" and
                .promotion_eligible == false and .promotion_passed == false and
                .causal_performance_claim_eligible == false and
                .cpu_pair_exclusive == false and .status == "complete" and
                .campaign_exit_status == 0 and .evidence_valid == true and
                .active_affinity_supervisor_executed == false and
                .supervisor_role == "retained-active-v1-verifier-only" and
                .producer_verification_passed == true and
                .producer_verification_mode == "manifest-without-affinity-binding" and
                .independent_preseal_audit_passed == true and
                .independent_auditor_supervision_mode == "absent" and
                .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
                .cpu == 52 and .sibling == 116 and
                .ratio_policy == {
                    ordinary_and_one_shot_are_separate:true,
                    ratios_are_separate_correlated_and_must_not_be_multiplied:true,
                    combined_or_stacked_ratio_emitted:false,
                    same_binary_ratio_is_another_campaign:true,
                    clustered_t_interval_reported_as_nominal_under_shared_host_load:true
                } and
                .postseal_policy ==
                    "passive shared-host observations always publish NOT_PROMOTED.json, independent of the observed speed gate"
            ' "$verified_core/manifest.json" >/dev/null
        else
            test "$evidence_generation" = passive-v2
            validate_v18_complete_core_static "$verified_core/manifest.json"
            if [[ -z "$semantic_complete_status" ]]; then
                verify_v18_complete_core_claim_bindings_preflight \
                    "$verified_core"
            fi
            /usr/bin/jq -e --slurpfile status "$status_file" '
                keys == ([
                    "acquisition_generation", "active_affinity_supervisor_executed",
                    "attempt", "attempt_budget", "attempt_lineage_sha256",
                    "audit_sha256", "auditor_sha256",
                    "baseline_archive_sha256_pre_post",
                    "baseline_binary_sha256_pre_post", "baseline_commit",
                    "build_order_normalizer_sha256", "build_reproduction",
                    "campaign_exit_status", "campaign_manifest_sha256",
                    "campaign_raw_sha256", "candidate_archive_sha256_pre_post",
                    "candidate_binary_sha256_pre_post",
                    "candidate_test_normalization_report_sha256",
                    "candidate_timing_normalization_report_sha256",
                    "canonical_lock", "causal_performance_claim_eligible",
                    "controller_affinity_sha256", "cpu", "cpu_pair_exclusive",
                    "environment_census_post_sha256",
                    "environment_census_pre_sha256", "evidence_class",
                    "evidence_valid", "independent_auditor_scope",
                    "independent_auditor_supervision_mode",
                    "independent_preseal_audit_passed", "isolation_claim",
                    "passive_census_sha256", "passive_environment_policy_sha256",
                    "performance_gate_passed", "postseal_policy",
                    "producer_verification_mode", "producer_verification_passed",
                    "promotion_eligible", "promotion_passed", "ratio_policy",
                    "runner_sha256", "schema", "sibling", "source_commit",
                    "source_tree", "sse2neon_commit",
                    "sse2neon_source_archive_sha256", "status",
                    "supervisor_role", "supervisor_sha256",
                    "out_of_window_benchmark_cpu_nonidle_jiffies",
                    "out_of_window_reserved_sibling_nonidle_jiffies",
                    "windowed_benchmark_cpu_nonidle_excess_jiffies",
                    "windowed_reserved_sibling_nonidle_jiffies", "wrapper_sha256"
                ] | sort) and
                .acquisition_generation == "passive-v2" and
                .attempt == $status[0].attempt and
                .attempt_budget == $status[0].attempt_budget and
                .attempt_lineage_sha256 ==
                    $status[0].attempt_lineage_sha256 and
                .attempt_budget == 3 and
                .evidence_class ==
                    "passive-windowed-shared-host-observation/v1" and
                .promotion_eligible == false and .promotion_passed == false and
                .causal_performance_claim_eligible == false and
                .cpu_pair_exclusive == false and .status == "complete" and
                .campaign_exit_status == 0 and .evidence_valid == true and
                .active_affinity_supervisor_executed == false and
                .supervisor_role == "retained-active-v1-verifier-only" and
                .producer_verification_passed == true and
                .producer_verification_mode == "manifest-without-affinity-binding" and
                .independent_preseal_audit_passed == true and
                .independent_auditor_supervision_mode == "windowed" and
                .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
                .cpu == 52 and .sibling == 116 and
                .windowed_benchmark_cpu_nonidle_excess_jiffies == 0 and
                .windowed_benchmark_cpu_nonidle_excess_jiffies ==
                    $status[0].windowed_benchmark_cpu_nonidle_excess_jiffies and
                .windowed_reserved_sibling_nonidle_jiffies == 0 and
                .windowed_reserved_sibling_nonidle_jiffies ==
                    $status[0].windowed_reserved_sibling_nonidle_jiffies and
                (.out_of_window_benchmark_cpu_nonidle_jiffies | type) == "number" and
                .out_of_window_benchmark_cpu_nonidle_jiffies ==
                    (.out_of_window_benchmark_cpu_nonidle_jiffies | floor) and
                .out_of_window_benchmark_cpu_nonidle_jiffies >= 0 and
                .out_of_window_benchmark_cpu_nonidle_jiffies ==
                    $status[0].out_of_window_benchmark_cpu_nonidle_jiffies and
                (.out_of_window_reserved_sibling_nonidle_jiffies | type) == "number" and
                .out_of_window_reserved_sibling_nonidle_jiffies ==
                    (.out_of_window_reserved_sibling_nonidle_jiffies | floor) and
                .out_of_window_reserved_sibling_nonidle_jiffies >= 0 and
                .out_of_window_reserved_sibling_nonidle_jiffies ==
                    $status[0].out_of_window_reserved_sibling_nonidle_jiffies and
                .ratio_policy == {
                    ordinary_and_one_shot_are_separate:true,
                    ratios_are_separate_correlated_and_must_not_be_multiplied:true,
                    combined_or_stacked_ratio_emitted:false,
                    same_binary_ratio_is_another_campaign:true,
                    clustered_t_interval_reported_as_nominal_under_shared_host_load:true
                } and
                .postseal_policy ==
                    "passive-v2 shared-host observations always publish NOT_PROMOTED.json, independent of the observed speed gate; at most three preregistered attempts and every outcome is retained"
            ' "$verified_core/manifest.json" >/dev/null
        fi
        for artifact_hash_binding in \
            candidate_binary_sha256_pre_post:build-closure/replay-candidate/bench_leopard2 \
            candidate_archive_sha256_pre_post:build-closure/replay-candidate/libleopard.a \
            baseline_binary_sha256_pre_post:build-closure/replay-baseline/leopard_main_benchmark \
            baseline_archive_sha256_pre_post:build-closure/replay-baseline/libleopard_main_exact.a \
            campaign_raw_sha256:campaign/raw.json; do
            artifact_hash_field=${artifact_hash_binding%%:*}
            artifact_hash_path=${artifact_hash_binding#*:}
            artifact_hash_expected="$(/usr/bin/jq -er \
                --arg field "$artifact_hash_field" \
                '.[$field] | strings |
                    select(length == 64 and test("^[0-9a-f]{64}$"))' \
                "$verified_core/manifest.json")" || return 1
            test "$(/usr/bin/sha256sum \
                "$verified_core/$artifact_hash_path" | \
                /usr/bin/cut -d' ' -f1)" = "$artifact_hash_expected"
        done
        verified_performance_gate_passed="$(/usr/bin/jq -er '
            [
                .analysis["gf16-high-full"].encode
                    .promotion_lower_bound_at_least_1_05,
                .analysis["gf16-high-full"].one_shot_encode
                    .promotion_lower_bound_at_least_1_05
            ] |
            if all(.[]; type == "boolean") then all
            else error("invalid gates") end |
            tostring
        ' "$verified_core/campaign/manifest.json")" || return 1
        /usr/bin/jq -e \
            --slurpfile status "$status_file" \
            --slurpfile closure "$verified_core/build-closure.json" \
            --slurpfile campaign "$verified_core/campaign/manifest.json" \
            --slurpfile audit "$verified_core/audit.json" \
            --arg main_commit "$main_commit" \
            --argjson gate "$verified_performance_gate_passed" '
            .campaign_raw_sha256 == $campaign[0].raw.sha256 and
            .campaign_raw_sha256 == $audit[0].raw.sha256 and
            .candidate_binary_sha256_pre_post ==
                $closure[0].candidate.binary_sha256 and
            .candidate_binary_sha256_pre_post ==
                $campaign[0].identities.candidate_executable.sha256 and
            .candidate_archive_sha256_pre_post ==
                $closure[0].candidate.archive_sha256 and
            .candidate_archive_sha256_pre_post ==
                $campaign[0].identities.candidate_archive.sha256 and
            .baseline_binary_sha256_pre_post ==
                $closure[0].baseline.binary_sha256 and
            .baseline_binary_sha256_pre_post ==
                $campaign[0].identities.baseline_executable.sha256 and
            .baseline_archive_sha256_pre_post ==
                $closure[0].baseline.archive_sha256 and
            .baseline_archive_sha256_pre_post ==
                $campaign[0].identities.baseline_archive.sha256 and
            .performance_gate_passed == $gate and
            $status[0].performance_gate_passed == $gate and
            (.source_tree | type) == "string" and
            (.source_tree | length) == 40 and
            (.source_tree | test("^[0-9a-f]{40}$")) and
            .source_tree == $closure[0].candidate.tree and
            .source_tree == $campaign[0].identities.candidate_source.tree and
            .source_tree == $status[0].source_tree and
            .source_commit == $closure[0].candidate.commit and
            .source_commit == $campaign[0].identities.candidate_source.head and
            .source_commit == $status[0].source_commit and
            .baseline_commit == $main_commit and
            .baseline_commit == $closure[0].baseline.commit and
            .baseline_commit == $campaign[0].identities.baseline_source.head and
            .baseline_commit == $audit[0].contract.baseline_main_commit and
            .baseline_commit == $status[0].baseline_commit
        ' "$verified_core/manifest.json" >/dev/null
        if [[ "${status_file##*/}" == COMPLETED.json ]]; then
            test "$(/usr/bin/jq -er '.candidate_binary_sha256 | strings' \
                "$status_file")" = \
                "$(/usr/bin/jq -er \
                    '.candidate_binary_sha256_pre_post | strings' \
                    "$verified_core/manifest.json")"
            test "$(/usr/bin/jq -er '.baseline_binary_sha256 | strings' \
                "$status_file")" = \
                "$(/usr/bin/jq -er \
                    '.baseline_binary_sha256_pre_post | strings' \
                    "$verified_core/manifest.json")"
        fi
        test "$(/usr/bin/jq -er '.promotion_passed | booleans | tostring' \
            "$verified_core/manifest.json")" = false
        test "$(/usr/bin/jq -er \
            '.performance_gate_passed | booleans | tostring' \
            "$verified_core/manifest.json")" = \
            "$(/usr/bin/jq -er \
                '.performance_gate_passed | booleans | tostring' \
                "$status_file")"
        test "$(/usr/bin/jq -er '.source_commit | strings' \
            "$verified_core/manifest.json")" = \
            "$(/usr/bin/jq -er '.source_commit | strings' "$status_file")"
        test "$(/usr/bin/jq -er '.source_tree | strings' \
            "$verified_core/manifest.json")" = \
            "$(/usr/bin/jq -er '.source_tree | strings' "$status_file")"
        test "$(/usr/bin/jq -er '.baseline_commit | strings' \
            "$verified_core/manifest.json")" = \
            "$(/usr/bin/jq -er '.baseline_commit | strings' "$status_file")"
        test "$(/usr/bin/sha256sum "$verified_core/run-authoritative.sh" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.wrapper_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/sha256sum \
            "$verified_core/audit_v17_gfni_main_compare.py" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.auditor_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/sha256sum "$verified_core/run_abba.py" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.runner_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/sha256sum \
            "$verified_core/leopard2_affinity_supervisor.py" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.supervisor_sha256 | strings' \
                "$verified_core/manifest.json")"
        if [[ "$evidence_generation" == passive-v1 ||
              "$evidence_generation" == passive-v2 ]]; then
            test "$(/usr/bin/sha256sum \
                "$verified_core/passive_environment_census.py" | \
                /usr/bin/cut -d' ' -f1)" = \
                "$(/usr/bin/jq -er '.passive_census_sha256 | strings' \
                    "$verified_core/manifest.json")"
            for passive_hash_binding in \
                controller-affinity.json:controller_affinity_sha256 \
                environment-census-pre.json:environment_census_pre_sha256 \
                environment-census-post.json:environment_census_post_sha256 \
                passive-environment-policy.json:passive_environment_policy_sha256; do
                passive_file=${passive_hash_binding%%:*}
                passive_field=${passive_hash_binding#*:}
                test "$(/usr/bin/sha256sum \
                    "$verified_core/$passive_file" | /usr/bin/cut -d' ' -f1)" = \
                    "$(/usr/bin/jq -er --arg field "$passive_field" \
                        '.[$field] | strings' "$verified_core/manifest.json")"
            done
            if [[ "$evidence_generation" == passive-v1 ]]; then
                /usr/bin/jq -e --slurpfile audit "$verified_core/audit.json" \
                --slurpfile status "$status_file" \
                --slurpfile census_pre \
                    "$verified_core/environment-census-pre.json" \
                --slurpfile census_post \
                    "$verified_core/environment-census-post.json" '
                .acquisition_generation == "passive-v1" and
                .acquisition_generation == $audit[0].acquisition_generation and
                .acquisition_generation == $status[0].acquisition_generation and
                .evidence_class == "passive-shared-host-observation/v1" and
                .evidence_class == $audit[0].evidence_class and
                .evidence_class == $status[0].evidence_class and
                .active_affinity_supervisor_executed == false and
                .supervisor_role == "retained-active-v1-verifier-only" and
                .promotion_eligible == false and
                .promotion_eligible == $audit[0].promotion_eligible and
                .promotion_eligible == $status[0].promotion_eligible and
                .promotion_passed == false and
                .promotion_passed == $audit[0].promotion_passed and
                .promotion_passed == $status[0].promotion_passed and
                .causal_performance_claim_eligible == false and
                .causal_performance_claim_eligible ==
                    $audit[0].causal_performance_claim_eligible and
                .causal_performance_claim_eligible ==
                    $status[0].causal_performance_claim_eligible and
                .cpu_pair_exclusive == false and
                .cpu_pair_exclusive == $audit[0].cpu_pair_exclusive and
                .cpu_pair_exclusive == $status[0].cpu_pair_exclusive and
                $audit[0].contamination.clock_ticks_per_second == 100 and
                $audit[0].contamination.clock_ticks_per_second ==
                    $census_pre[0].clock_ticks_per_second and
                $audit[0].contamination.clock_ticks_per_second ==
                    $census_post[0].clock_ticks_per_second and
                .benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies ==
                    $audit[0].contamination.benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies and
                .benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies ==
                    $status[0].benchmark_cpu_nonidle_excess_over_child_wall_ceiling_jiffies and
                .isolation_claim == $audit[0].isolation_claim and
                .isolation_claim == $status[0].isolation_claim
                ' "$verified_core/manifest.json" >/dev/null
            else
                /usr/bin/jq -e \
                    --slurpfile audit "$verified_core/audit.json" \
                    --slurpfile status "$status_file" \
                    --slurpfile policy \
                        "$verified_core/passive-environment-policy.json" '
                    .acquisition_generation == "passive-v2" and
                    .acquisition_generation == $audit[0].acquisition_generation and
                    .acquisition_generation == $status[0].acquisition_generation and
                    .acquisition_generation == $policy[0].acquisition_generation and
                    .attempt == $status[0].attempt and
                    .attempt_budget == $status[0].attempt_budget and
                    .evidence_class ==
                        "passive-windowed-shared-host-observation/v1" and
                    .evidence_class == $audit[0].evidence_class and
                    .evidence_class == $status[0].evidence_class and
                    .evidence_class == $policy[0].evidence_class and
                    .active_affinity_supervisor_executed == false and
                    .promotion_eligible == false and
                    .promotion_eligible == $audit[0].promotion_eligible and
                    .promotion_eligible == $status[0].promotion_eligible and
                    .promotion_eligible == $policy[0].promotion_eligible and
                    .promotion_passed == false and
                    .causal_performance_claim_eligible == false and
                    .cpu_pair_exclusive == false and
                    .windowed_benchmark_cpu_nonidle_excess_jiffies ==
                        $audit[0].windowed_observation.windowed.benchmark_cpu_nonidle_excess_jiffies and
                    .windowed_benchmark_cpu_nonidle_excess_jiffies ==
                        $policy[0].windowed_contamination.windowed.benchmark_cpu_nonidle_excess_jiffies and
                    .windowed_benchmark_cpu_nonidle_excess_jiffies ==
                        $status[0].windowed_benchmark_cpu_nonidle_excess_jiffies and
                    .windowed_benchmark_cpu_nonidle_excess_jiffies == 0 and
                    .windowed_reserved_sibling_nonidle_jiffies ==
                        $audit[0].windowed_observation.windowed.reserved_sibling_nonidle_jiffies and
                    .windowed_reserved_sibling_nonidle_jiffies ==
                        $policy[0].windowed_contamination.windowed.reserved_sibling_nonidle_jiffies and
                    .windowed_reserved_sibling_nonidle_jiffies ==
                        $status[0].windowed_reserved_sibling_nonidle_jiffies and
                    .windowed_reserved_sibling_nonidle_jiffies == 0 and
                    .out_of_window_benchmark_cpu_nonidle_jiffies ==
                        $audit[0].windowed_observation.out_of_window.benchmark_cpu_nonidle_jiffies and
                    .out_of_window_benchmark_cpu_nonidle_jiffies ==
                        $policy[0].outer_disclosure.isolation_out_of_window.benchmark_cpu_nonidle_jiffies and
                    .out_of_window_benchmark_cpu_nonidle_jiffies ==
                        $status[0].out_of_window_benchmark_cpu_nonidle_jiffies and
                    .out_of_window_reserved_sibling_nonidle_jiffies ==
                        $audit[0].windowed_observation.out_of_window.reserved_sibling_nonidle_jiffies and
                    .out_of_window_reserved_sibling_nonidle_jiffies ==
                        $policy[0].outer_disclosure.isolation_out_of_window.reserved_sibling_nonidle_jiffies and
                    .out_of_window_reserved_sibling_nonidle_jiffies ==
                        $status[0].out_of_window_reserved_sibling_nonidle_jiffies and
                    .isolation_claim == $audit[0].isolation_claim and
                    .isolation_claim == $status[0].isolation_claim
                ' "$verified_core/manifest.json" >/dev/null
            fi
        fi
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-order-normalizer.py" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.build_order_normalizer_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/sha256sum \
            "$verified_core/build-closure/candidate-timing-order-normalization/normalized/report.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.candidate_timing_normalization_report_sha256 | strings' \
                "$verified_core/manifest.json")"
        /usr/bin/jq -e '
            .build_reproduction.candidate_timing_qualified_semantic_equal == true and
            (.build_reproduction.candidate_timing_raw_tree_file_bytes_identical |
                type == "boolean") and
            .build_reproduction.candidate_tests_qualified_semantic_equal == true and
            (.build_reproduction.candidate_tests_raw_tree_file_bytes_identical |
                type == "boolean") and
            .build_reproduction.baseline_raw_byte_identical == true and
            .independent_auditor_scope ==
                "campaign semantics only; build-order qualification is separately recomputed by the frozen wrapper normalizer"
        ' "$verified_core/manifest.json" >/dev/null
        test "$(/usr/bin/sha256sum "$verified_core/manifest.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.core_manifest_sha256 | strings' \
                "$status_file")"
        test "$(/usr/bin/sha256sum "$verified_core/campaign/manifest.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.campaign_manifest_sha256 | strings' \
                "$status_file")"
        campaign_manifest_sha="$(/usr/bin/sha256sum \
            "$verified_core/campaign/manifest.json" | \
            /usr/bin/cut -d' ' -f1)"
        test "$campaign_manifest_sha" = \
            "$(/usr/bin/jq -er '.campaign_manifest_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/sha256sum "$verified_core/sse2neon-source.tar" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er \
                '.sse2neon_source_archive_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/jq -er '.sse2neon.commit | strings' \
            "$verified_core/build-closure.json")" = \
            "$(/usr/bin/jq -er '.sse2neon_commit | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/jq -er \
            '.sse2neon.source_archive_sha256 | strings' \
            "$verified_core/build-closure.json")" = \
            "$(/usr/bin/jq -er \
                '.sse2neon_source_archive_sha256 | strings' \
                "$verified_core/manifest.json")"
        test "$(/usr/bin/jq -er '.sse2neon.archive_prefix | strings' \
            "$verified_core/build-closure.json")" = sse2neon-source/
        test "$(/usr/bin/jq -er \
            '.sse2neon.reproduced_from_candidate_and_baseline_clones | booleans | tostring' \
            "$verified_core/build-closure.json")" = true
        test "$(/usr/bin/sha256sum "$verified_core/audit.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.audit_sha256 | strings' \
                "$status_file")"
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/audit.json")" = "$audit_schema"
        /usr/bin/jq -e --slurpfile campaign \
            "$verified_core/campaign/manifest.json" '
            .ordinary_one_item_batch ==
                $campaign[0].analysis["gf16-high-full"].encode and
            .one_shot_encode ==
                $campaign[0].analysis["gf16-high-full"].one_shot_encode and
            .ratios_are_separate_correlated_and_must_not_be_multiplied == true
        ' "$verified_core/result-summary.json" >/dev/null
        if [[ "$evidence_generation" == passive-v1 ]]; then
            test "$(/usr/bin/jq -er '.schema | strings' \
                "$verified_core/result-summary.json")" = \
                leopard2-v17-gfni-main-passive-result-summary/v2
            /usr/bin/jq -e --slurpfile audit "$verified_core/audit.json" \
                --slurpfile status "$status_file" \
                --slurpfile policy \
                    "$verified_core/passive-environment-policy.json" '
                keys == ([
                    "acquisition_generation", "active_affinity_supervisor_executed",
                    "causal_performance_claim_eligible",
                    "clustered_t_interval_reported_as_nominal_under_shared_host_load",
                    "contamination", "cpu_pair_exclusive", "evidence_class",
                    "isolation_claim", "one_shot_encode", "ordinary_one_item_batch",
                    "outer_contamination",
                    "promotion_eligible", "promotion_passed", "ratio_orientation",
                    "ratios_are_separate_correlated_and_must_not_be_multiplied",
                    "same_binary_gfni_over_avx2_is_a_separate_campaign", "schema"
                ] | sort) and
                .acquisition_generation == "passive-v1" and
                .active_affinity_supervisor_executed == false and
                .acquisition_generation == $audit[0].acquisition_generation and
                .acquisition_generation == $status[0].acquisition_generation and
                .evidence_class == "passive-shared-host-observation/v1" and
                .evidence_class == $audit[0].evidence_class and
                .evidence_class == $status[0].evidence_class and
                .promotion_eligible == false and
                .promotion_passed == false and
                .causal_performance_claim_eligible == false and
                .cpu_pair_exclusive == false and
                .ratio_orientation ==
                    "exact_leopard1_native_time_over_leopard2_candidate_time" and
                .ratios_are_separate_correlated_and_must_not_be_multiplied == true and
                .clustered_t_interval_reported_as_nominal_under_shared_host_load == true and
                .same_binary_gfni_over_avx2_is_a_separate_campaign == true and
                .contamination == $audit[0].contamination and
                .outer_contamination == $policy[0].outer_contamination and
                .isolation_claim == $audit[0].isolation_claim
            ' "$verified_core/result-summary.json" >/dev/null
        elif [[ "$evidence_generation" == passive-v2 ]]; then
            test "$(/usr/bin/jq -er '.schema | strings' \
                "$verified_core/result-summary.json")" = \
                leopard2-v18-gfni-main-passive-result-summary/v1
            /usr/bin/jq -e --slurpfile audit "$verified_core/audit.json" \
                --slurpfile core "$verified_core/manifest.json" \
                --slurpfile status "$status_file" \
                --slurpfile lineage "$verified_core/attempt-lineage.json" \
                --slurpfile policy \
                    "$verified_core/passive-environment-policy.json" '
                keys == ([
                    "acquisition_generation", "active_affinity_supervisor_executed",
                    "attempt", "attempt_budget", "attempt_lineage_sha256",
                    "attempt_statement",
                    "causal_performance_claim_eligible",
                    "clustered_t_interval_reported_as_nominal_under_shared_host_load",
                    "cpu_pair_exclusive", "evidence_class", "isolation_claim",
                    "one_shot_encode", "ordinary_one_item_batch",
                    "outer_disclosure", "promotion_eligible", "promotion_passed",
                    "ratio_orientation",
                    "ratios_are_separate_correlated_and_must_not_be_multiplied",
                    "same_binary_gfni_over_avx2_is_a_separate_campaign", "schema",
                    "sealed_attempt_envelopes", "shared_host_exposure",
                    "windowed_observation"
                ] | sort) and
                .schema ==
                    "leopard2-v18-gfni-main-passive-result-summary/v1" and
                .acquisition_generation == "passive-v2" and
                .attempt == $core[0].attempt and .attempt == $status[0].attempt and
                .attempt_budget == 3 and
                .attempt_budget == $core[0].attempt_budget and
                .attempt_budget == $status[0].attempt_budget and
                .attempt_lineage_sha256 ==
                    $core[0].attempt_lineage_sha256 and
                .attempt_lineage_sha256 ==
                    $status[0].attempt_lineage_sha256 and
                .attempt_statement ==
                    ("attempt \(.attempt) of at most \(.attempt_budget)") and
                .sealed_attempt_envelopes ==
                    ([$lineage[0].prior_attempts[].envelope] + [
                        "/home/catid/leopard/.research/leopard-79h/\($status[0].source_commit[0:7])-v18-passive-main-a\(.attempt)"
                    ]) and
                .active_affinity_supervisor_executed == false and
                .evidence_class ==
                    "passive-windowed-shared-host-observation/v1" and
                .promotion_eligible == false and .promotion_passed == false and
                .causal_performance_claim_eligible == false and
                .cpu_pair_exclusive == false and
                .ratio_orientation ==
                    "exact_leopard1_native_time_over_leopard2_candidate_time" and
                .ratios_are_separate_correlated_and_must_not_be_multiplied == true and
                .clustered_t_interval_reported_as_nominal_under_shared_host_load == true and
                .same_binary_gfni_over_avx2_is_a_separate_campaign == true and
                .windowed_observation == $audit[0].windowed_observation and
                .outer_disclosure == $policy[0].outer_disclosure and
                .shared_host_exposure == $policy[0].shared_host_exposure and
                .isolation_claim == $audit[0].isolation_claim
            ' "$verified_core/result-summary.json" >/dev/null
        else
            test "$(/usr/bin/jq -er '.schema | strings' \
                "$verified_core/result-summary.json")" = \
                leopard2-v17-gfni-main-result-summary/v1
        fi
        test "$(/usr/bin/sha256sum \
            "$verified_postseal_audit" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.postseal_audit_sha256 | strings' \
                "$status_file")"
        if [[ "$evidence_generation" == active-v1 ]]; then
            test "$(/usr/bin/sha256sum \
                "$verified_core/affinity-report.json" | \
                /usr/bin/cut -d' ' -f1)" = \
                "$(/usr/bin/jq -er '.affinity_report_sha256 | strings' \
                    "$verified_core/manifest.json")"
            test "$(/usr/bin/sha256sum \
                "$verified_core/affinity-binding.json" | \
                /usr/bin/cut -d' ' -f1)" = \
                "$(/usr/bin/jq -er '.affinity_binding_sha256 | strings' \
                    "$verified_core/manifest.json")"
            /usr/bin/python3 -I -S -B \
                "$verified_core/leopard2_affinity_supervisor.py" verify-report \
                --report "$verified_core/affinity-report.json" \
                > "$verifier_tmp/supervisor-report-verification.log" \
                2> "$verifier_tmp/supervisor-report-verification-stderr.log"
            /usr/bin/python3 -I -S -B \
                "$verified_core/leopard2_affinity_supervisor.py" verify-binding \
                --binding "$verified_core/affinity-binding.json" \
                --manifest "$verified_core/campaign/manifest.json" \
                --manifest-sha256 "$campaign_manifest_sha" \
                > "$verifier_tmp/supervisor-binding-verification.log" \
                2> "$verifier_tmp/supervisor-binding-verification-stderr.log"
        else
            test "$(/usr/bin/jq -er '.supervision == null' \
                "$verified_core/campaign/manifest.json")" = true
            test "$(/usr/bin/jq -er '.supervision == null' \
                "$verified_core/campaign/raw.json")" = true
            census_verify_args=()
            controller_schema=leopard2-v17-passive-controller-affinity/v1
            controller_generation=passive-v1
            if [[ "$evidence_generation" == passive-v2 ]]; then
                census_verify_args=(--generation passive-v2)
                controller_schema=leopard2-v18-passive-controller-affinity/v1
                controller_generation=passive-v2
            fi
            /usr/bin/python3 -I -S -B \
                "$verified_core/passive_environment_census.py" verify \
                --input "$verified_core/environment-census-pre.json" --phase pre \
                "${census_verify_args[@]}"
            /usr/bin/python3 -I -S -B \
                "$verified_core/passive_environment_census.py" verify \
                --input "$verified_core/environment-census-post.json" --phase post \
                "${census_verify_args[@]}"
            /usr/bin/python3 -I -S -B \
                "$verified_core/passive_environment_census.py" compare \
                --pre "$verified_core/environment-census-pre.json" \
                --post "$verified_core/environment-census-post.json" \
                --raw "$verified_core/campaign/raw.json" \
                --controller "$verified_core/controller-affinity.json" \
                --output "$verifier_tmp/passive-environment-policy.json"
            /usr/bin/cmp "$verified_core/passive-environment-policy.json" \
                "$verifier_tmp/passive-environment-policy.json"
            if [[ "$evidence_generation" == passive-v2 ]]; then
                /usr/bin/jq -e '
                    .schema == "leopard2-passive-shared-host-policy/v2" and
                    .status == "complete" and
                    .acquisition_generation == "passive-v2" and
                    .evidence_class ==
                        "passive-windowed-shared-host-observation/v1" and
                    .policy_evaluation_complete == true and
                    .promotion_eligible == false and .promotion_passed == false and
                    .causal_performance_claim_eligible == false and
                    .cpu_pair_exclusive == false and
                    .windowed_contamination.gated == true and
                    .windowed_contamination.retained_window_count == 12 and
                    .windowed_contamination.all_benchmark_cpu_excess_zero == true and
                    .windowed_contamination.all_reserved_sibling_nonidle_zero == true and
                    .outer_disclosure.gated == false
                ' "$verified_core/passive-environment-policy.json" >/dev/null
            fi
            test "$(/usr/bin/jq -er \
                '.promotion_eligible | booleans | tostring' \
                "$verified_core/passive-environment-policy.json")" = false
            test "$(/usr/bin/jq -er \
                '.cpu_pair_exclusive | booleans | tostring' \
                "$verified_core/passive-environment-policy.json")" = false
            /usr/bin/jq -e \
                --arg controller_schema "$controller_schema" \
                --arg controller_generation "$controller_generation" '
                keys == (["acquisition_generation","active_affinity_supervisor_executed","affinity_mutation_scope","after_allowed_cpus","before_allowed_cpus","benchmark_cpu","runner_launch_allowed_cpus","schema","reserved_sibling","wrapper_pid"] | sort) and
                .schema == $controller_schema and
                .acquisition_generation == $controller_generation and
                .active_affinity_supervisor_executed == false and
                .affinity_mutation_scope ==
                    "wrapper-process-and-owned-descendants-only" and
                .benchmark_cpu == 52 and .reserved_sibling == 116 and
                .wrapper_pid > 0 and
                .runner_launch_allowed_cpus == .before_allowed_cpus and
                (.before_allowed_cpus | index(52)) != null and
                (.before_allowed_cpus | index(116)) != null and
                (.after_allowed_cpus | length) > 0 and
                (.after_allowed_cpus | index(52)) == null and
                (.after_allowed_cpus | index(116)) == null and
                (.after_allowed_cpus ==
                    [.before_allowed_cpus[] | select(. != 52 and . != 116)])
            ' "$verified_core/controller-affinity.json" >/dev/null
            /usr/bin/jq -e \
                --slurpfile controller "$verified_core/controller-affinity.json" \
                --slurpfile raw "$verified_core/campaign/raw.json" \
                --slurpfile campaign_manifest \
                    "$verified_core/campaign/manifest.json" '
                keys == ["allowed_cpus"] and
                .allowed_cpus == $controller[0].before_allowed_cpus and
                $controller[0].runner_launch_allowed_cpus ==
                    $raw[0].campaign.allowed_cpu_set_at_launch and
                $controller[0].runner_launch_allowed_cpus ==
                    $campaign_manifest[0].campaign.allowed_cpu_set_at_launch
            ' "$verified_core/wrapper-launch-affinity.json" >/dev/null
            verified_timeout_argument="$(/usr/bin/jq -er \
                "$passive_timeout_argument_jq" \
                "$verified_core/campaign/raw.json")" || return 1
            /usr/bin/jq -e \
                --arg campaign_dir "$verified_core/campaign" \
                --arg timeout_argument "$verified_timeout_argument" \
                --slurpfile controller "$verified_core/controller-affinity.json" \
                --slurpfile raw "$verified_core/campaign/raw.json" '
                . == [
                    "/usr/bin/taskset", "-c",
                    ($controller[0].runner_launch_allowed_cpus |
                        map(tostring) | join(",")),
                    "/usr/bin/python3", "-I", "-S", "-B",
                    $raw[0].input_specification.runner, "run",
                    "--baseline", $raw[0].input_specification.baseline_executable,
                    "--candidate", $raw[0].input_specification.candidate_executable,
                    "--baseline-archive",
                        $raw[0].input_specification.baseline_archive,
                    "--candidate-archive",
                        $raw[0].input_specification.candidate_archive,
                    "--baseline-build-dir",
                        $raw[0].input_specification.baseline_build_dir,
                    "--candidate-build-dir",
                        $raw[0].input_specification.candidate_build_dir,
                    "--baseline-source-root",
                        $raw[0].input_specification.baseline_source_root,
                    "--candidate-source-root",
                        $raw[0].input_specification.candidate_source_root,
                    "--candidate-commit",
                        $raw[0].input_specification.candidate_commit,
                    "--baseline-native", "--candidate-mode",
                        $raw[0].campaign.candidate_mode,
                    "--reservation-file", $raw[0].reservation.path,
                    "--output", $campaign_dir,
                    "--cpu", ($raw[0].campaign.benchmark_cpu | tostring),
                    "--reserved-sibling",
                        ($raw[0].campaign.reserved_sibling | tostring),
                    "--taskset", $raw[0].input_specification.taskset,
                    "--ldd", $raw[0].input_specification.ldd,
                    "--preset", "v17-gfni-encode",
                    "--reuse", ($raw[0].campaign.reuse | tostring),
                    "--iterations", ($raw[0].campaign.iterations | tostring),
                    "--warmup", ($raw[0].campaign.warmup | tostring),
                    "--timeout", $timeout_argument
                ] and $raw[0].input_specification.baseline_pure_avx2 == false
            ' "$verified_core/campaign-command.json" >/dev/null
        fi
        replay_campaign="$verifier_tmp/campaign"
        reconstruct_owner_only_campaign_tree \
            "$verified_core/campaign" "$replay_campaign"
        /usr/bin/diff -qr "$verified_core/campaign" "$replay_campaign" \
            > "$verifier_tmp/reconstructed-campaign-byte-diff.txt"
        replay_controller_root="$verifier_tmp/controller"
        /usr/bin/mkdir -p \
            "$replay_controller_root/experiments/leopard2/main_compare" \
            "$replay_controller_root/experiments/leopard2/decoder_dispatch" \
            "$replay_controller_root/tools"
        /usr/bin/find "$replay_controller_root" -type d \
            -exec /usr/bin/chmod 0700 {} +
        /usr/bin/cp --reflink=never "$verified_core/run_abba.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/run_abba.py"
        /usr/bin/cp --reflink=never "$verified_core/git_capture.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/git_capture.py"
        /usr/bin/cp --reflink=never \
            "$verified_core/balanced_evidence_common.py" \
            "$replay_controller_root/experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
        /usr/bin/cp --reflink=never \
            "$verified_core/leopard2_affinity_supervisor.py" \
            "$replay_controller_root/tools/leopard2_affinity_supervisor.py"
        /usr/bin/cp --reflink=never \
            "$verified_core/leopard2_build_provenance.py" \
            "$replay_controller_root/tools/leopard2_build_provenance.py"
        /usr/bin/find "$replay_controller_root" -type f \
            -exec /usr/bin/chmod 0400 {} +
        /usr/bin/cmp "$verified_core/run_abba.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/run_abba.py"
        /usr/bin/cmp "$verified_core/git_capture.py" \
            "$replay_controller_root/experiments/leopard2/main_compare/git_capture.py"
        /usr/bin/cmp "$verified_core/balanced_evidence_common.py" \
            "$replay_controller_root/experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
        # The frozen supervisor above verifies the v8 binding against the
        # original sealed absolute manifest path and exact manifest hash.  The
        # producer requires an owner-only 0700/0600 evidence tree, so replay it
        # only against this byte-identical reconstruction.  A relocated copy
        # cannot honestly reuse the absolute-path-bound v8 binding.
        /usr/bin/python3 -I -S -B \
            "$replay_controller_root/experiments/leopard2/main_compare/run_abba.py" \
            verify --manifest "$replay_campaign/manifest.json" \
            --no-current-input-check \
            > "$verifier_tmp/producer-retained-verification.log" \
            2> "$verifier_tmp/producer-retained-verification-stderr.log"
        audit_replay_args=()
        if [[ "$evidence_generation" == passive-v1 ]]; then
            audit_replay_args=(--supervision absent)
        elif [[ "$evidence_generation" == passive-v2 ]]; then
            audit_replay_args=(--supervision windowed)
        fi
        /usr/bin/python3 -I -S -B \
            "$verified_core/audit_v17_gfni_main_compare.py" \
            --manifest "$verified_core/campaign/manifest.json" \
            "${audit_replay_args[@]}" \
            --output "$verifier_tmp/audit.json" \
            > "$verifier_tmp/audit-summary.json" \
            2> "$verifier_tmp/audit-stderr.log"
        /usr/bin/cmp "$verified_core/audit.json" "$verifier_tmp/audit.json"
        (
            cd "$verified_core"
            /usr/bin/sha256sum -c SHA256SUMS
        ) > "$verifier_tmp/post-replay-core-sha256-check.txt"
        (
            cd "$verified_envelope"
            /usr/bin/sha256sum -c SHA256SUMS
        ) > "$verifier_tmp/post-replay-envelope-sha256-check.txt"
        verify_sealed_tree "$verified_core"
        verify_tree_metadata_manifest "$verified_envelope"
    else
        test "$campaign_exit_value" -ne 0
        if [[ "$evidence_generation" == failed-v1 ]]; then
            test ! -e "$verified_core/attempt-lineage.json"
        fi
        if /usr/bin/jq -e \
                'has("stage") and (has("core_sha256sums_verified") | not)' \
                "$status_file" >/dev/null; then
            if [[ "$evidence_generation" == failed-v1 ]]; then
                /usr/bin/jq -e '
                    keys == (["campaign_exit_status","core_sha256sums_sha256",
                        "promotion_passed","schema","source_commit","source_tree",
                        "stage","status"] | sort) and
                    .schema == "leopard2-v17-gfni-main-failed-envelope/v1" and
                    .status == "failed" and .promotion_passed == false and
                    (.campaign_exit_status | type == "number") and
                    .campaign_exit_status != 0 and
                    (.stage | type == "string" and length > 0)
                ' "$status_file" >/dev/null
            else
                /usr/bin/jq -e '
                    keys == (["acquisition_generation","attempt","attempt_budget",
                        "attempt_lineage_sha256",
                        "campaign_exit_status","core_sha256sums_sha256",
                        "promotion_passed","schema","source_commit","source_tree",
                        "stage","status"] | sort) and
                    .schema ==
                        "leopard2-v18-gfni-main-failed-envelope/v1" and
                    .acquisition_generation == "passive-v2" and
                    (.attempt | type) == "number" and
                    .attempt == (.attempt | floor) and
                    .attempt >= 1 and .attempt <= 3 and .attempt_budget == 3 and
                    .status == "failed" and .promotion_passed == false and
                    (.campaign_exit_status | type) == "number" and
                    .campaign_exit_status == (.campaign_exit_status | floor) and
                    .campaign_exit_status >= 1 and
                    .campaign_exit_status <= 255 and
                    (.stage | type) == "string" and (.stage | length) > 0
                ' "$status_file" >/dev/null
                test -f "$verified_core/WRAPPER_FAILURE.json"
                validate_v18_wrapper_failure_binding \
                    "$verified_core/WRAPPER_FAILURE.json" "$status_file"
            fi
        else
            if /usr/bin/jq -e 'has("failure_verified")' \
                    "$status_file" >/dev/null; then
                test -f "$verified_core/manifest.json"
                if [[ "$evidence_generation" == failed-v1 ]]; then
                    test "$(/usr/bin/jq -er '.schema | strings' \
                        "$verified_core/manifest.json")" = \
                        leopard2-v17-gfni-main-failed-core-manifest/v1
                    /usr/bin/jq -e --slurpfile status "$status_file" '
                    keys == (["baseline_binary_sha256","baseline_commit",
                        "campaign_exit_status","candidate_binary_sha256",
                        "canonical_lock","cpu","failure_sha256",
                        "failure_verified","failure_verify_status",
                        "promotion_passed","schema","sibling","source_commit",
                        "source_tree","status"] | sort) and
                    .status == "failed" and .promotion_passed == false and
                    .campaign_exit_status != 0 and
                    (.failure_verify_status | type == "number") and
                    (.failure_verified | type == "boolean") and
                    .failure_verified == $status[0].failure_verified and
                    .source_commit == $status[0].source_commit and
                    .source_tree == $status[0].source_tree and
                    .campaign_exit_status == $status[0].campaign_exit_status and
                    .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
                    .baseline_commit ==
                        "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
                    .cpu == 52 and .sibling == 116
                    ' "$verified_core/manifest.json" >/dev/null
                    /usr/bin/jq -e '
                    keys == (["campaign_exit_status","core_sha256sums_sha256",
                        "failure_verified","promotion_passed","schema",
                        "source_commit","source_tree","status"] | sort) and
                    .status == "failed" and .promotion_passed == false and
                    .campaign_exit_status != 0 and
                    (.failure_verified | type == "boolean")
                    ' "$status_file" >/dev/null
                else
                    test "$(/usr/bin/jq -er '.schema | strings' \
                        "$verified_core/manifest.json")" = \
                        leopard2-v18-gfni-main-passive-failed-core-manifest/v1
                    validate_v18_failed_core_static \
                        "$verified_core/manifest.json"
                    /usr/bin/jq -e --slurpfile status "$status_file" '
                        keys == (["acquisition_generation","attempt","attempt_budget",
                            "attempt_lineage_sha256",
                            "baseline_binary_sha256","baseline_commit",
                            "campaign_exit_status","candidate_binary_sha256",
                            "canonical_lock","cpu","failure_sha256",
                            "failure_verified","failure_verify_status",
                            "promotion_passed","schema","sibling","source_commit",
                            "source_tree","status"] | sort) and
                        .acquisition_generation == "passive-v2" and
                        .attempt == $status[0].attempt and
                        .attempt_budget == $status[0].attempt_budget and
                        .attempt_lineage_sha256 ==
                            $status[0].attempt_lineage_sha256 and
                        .attempt_budget == 3 and .status == "failed" and
                        .promotion_passed == false and
                        (.campaign_exit_status | type) == "number" and
                        .campaign_exit_status ==
                            (.campaign_exit_status | floor) and
                        .campaign_exit_status >= 1 and
                        .campaign_exit_status <= 255 and
                        (.failure_verify_status | type) == "number" and
                        .failure_verify_status ==
                            (.failure_verify_status | floor) and
                        .failure_verify_status >= 0 and
                        .failure_verify_status <= 255 and
                        (.failure_verified | type) == "boolean" and
                        .failure_verified == (.failure_verify_status == 0) and
                        .failure_verified == $status[0].failure_verified and
                        .source_commit == $status[0].source_commit and
                        .source_tree == $status[0].source_tree and
                        .campaign_exit_status == $status[0].campaign_exit_status and
                        .canonical_lock == "/tmp/leopard-gf8-authoritative.lock" and
                        .baseline_commit ==
                            "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
                        .cpu == 52 and .sibling == 116
                    ' "$verified_core/manifest.json" >/dev/null
                    /usr/bin/jq -e '
                        keys == (["acquisition_generation","attempt","attempt_budget",
                            "attempt_lineage_sha256",
                            "campaign_exit_status","core_sha256sums_sha256",
                            "failure_verified","promotion_passed","schema",
                            "source_commit","source_tree","status"] | sort) and
                        .acquisition_generation == "passive-v2" and
                        .attempt_budget == 3 and .status == "failed" and
                        .promotion_passed == false and
                        (.campaign_exit_status | type) == "number" and
                        .campaign_exit_status ==
                            (.campaign_exit_status | floor) and
                        .campaign_exit_status >= 1 and
                        .campaign_exit_status <= 255 and
                        (.failure_verified | type) == "boolean"
                    ' "$status_file" >/dev/null
                    verify_v18_failed_core_claim_bindings "$verified_core"
                fi
                test "$(/usr/bin/sha256sum \
                    "$verified_core/build-closure/replay-candidate/bench_leopard2" | \
                    /usr/bin/cut -d' ' -f1)" = \
                    "$(/usr/bin/jq -er '.candidate_binary_sha256 | strings' \
                        "$verified_core/manifest.json")"
                test "$(/usr/bin/sha256sum \
                    "$verified_core/build-closure/replay-baseline/leopard_main_benchmark" | \
                    /usr/bin/cut -d' ' -f1)" = \
                    "$(/usr/bin/jq -er '.baseline_binary_sha256 | strings' \
                        "$verified_core/manifest.json")"
                failure_path="$verified_core/campaign/failure.json"
                retained_failure_sha="$(/usr/bin/jq -er \
                    '.failure_sha256 | if . == null then "null" else strings end' \
                    "$verified_core/manifest.json")"
                if [[ -f "$failure_path" ]]; then
                    test "$retained_failure_sha" != null
                    test "$(/usr/bin/sha256sum "$failure_path" | \
                        /usr/bin/cut -d' ' -f1)" = "$retained_failure_sha"
                    failed_verify_root="$(/usr/bin/mktemp -d \
                        /tmp/leopard-v17-failed-replay.XXXXXX)"
                    reconstruct_owner_only_campaign_tree \
                        "$verified_core/campaign" "$failed_verify_root/campaign"
                    /usr/bin/mkdir -p \
                        "$failed_verify_root/controller/experiments/leopard2/main_compare" \
                        "$failed_verify_root/controller/experiments/leopard2/decoder_dispatch" \
                        "$failed_verify_root/controller/tools"
                    /usr/bin/cp --reflink=never "$verified_core/run_abba.py" \
                        "$failed_verify_root/controller/experiments/leopard2/main_compare/run_abba.py"
                    /usr/bin/cp --reflink=never "$verified_core/git_capture.py" \
                        "$failed_verify_root/controller/experiments/leopard2/main_compare/git_capture.py"
                    /usr/bin/cp --reflink=never \
                        "$verified_core/balanced_evidence_common.py" \
                        "$failed_verify_root/controller/experiments/leopard2/decoder_dispatch/balanced_evidence_common.py"
                    /usr/bin/cp --reflink=never \
                        "$verified_core/leopard2_build_provenance.py" \
                        "$failed_verify_root/controller/tools/leopard2_build_provenance.py"
                    /usr/bin/find "$failed_verify_root/controller" -type d \
                        -exec /usr/bin/chmod 0700 {} +
                    /usr/bin/find "$failed_verify_root/controller" -type f \
                        -exec /usr/bin/chmod 0600 {} +
                    set +e
                    /usr/bin/python3 -I -S -B \
                        "$failed_verify_root/controller/experiments/leopard2/main_compare/run_abba.py" \
                        verify-failure \
                        --failure "$failed_verify_root/campaign/failure.json" \
                        > "$failed_verify_root/verification.log" 2>&1
                    replay_failure_status=$?
                    set -e
                    test "$replay_failure_status" -eq \
                        "$(/usr/bin/jq -er '.failure_verify_status | numbers' \
                            "$verified_core/manifest.json")"
                    if [[ "$(/usr/bin/jq -er \
                            '.failure_verified | booleans | tostring' \
                            "$verified_core/manifest.json")" == true ]]; then
                        test "$replay_failure_status" -eq 0
                    else
                        test "$replay_failure_status" -ne 0
                    fi
                else
                    test "$retained_failure_sha" = null
                    test "$(/usr/bin/jq -er \
                        '.failure_verified | booleans | tostring' \
                        "$verified_core/manifest.json")" = false
                    test "$(/usr/bin/jq -er '.failure_verify_status | numbers' \
                        "$verified_core/manifest.json")" -eq 1
                fi
            else
                if [[ "$evidence_generation" == failed-v1 ]]; then
                    /usr/bin/jq -e '
                    keys == (["baseline_commit","campaign_exit_status",
                        "core_sha256sums_sha256","core_sha256sums_verified",
                        "promotion_passed","schema","source_commit","source_tree",
                        "stage","status"] | sort) and
                    .status == "failed" and .promotion_passed == false and
                    .campaign_exit_status != 0 and
                    (.core_sha256sums_verified | type == "boolean") and
                    (.stage | type == "string" and length > 0)
                    ' "$status_file" >/dev/null
                else
                    /usr/bin/jq -e '
                        keys == (["acquisition_generation","attempt","attempt_budget",
                            "attempt_lineage_sha256",
                            "baseline_commit","campaign_exit_status",
                            "core_sha256sums_sha256","core_sha256sums_verified",
                            "promotion_passed","schema","source_commit","source_tree",
                            "stage","status"] | sort) and
                        .acquisition_generation == "passive-v2" and
                        .attempt_budget == 3 and .status == "failed" and
                        .promotion_passed == false and
                        (.campaign_exit_status | type) == "number" and
                        .campaign_exit_status ==
                            (.campaign_exit_status | floor) and
                        .campaign_exit_status >= 1 and
                        .campaign_exit_status <= 255 and
                        .baseline_commit ==
                            "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198" and
                        (.core_sha256sums_verified | type) == "boolean" and
                        (.stage | type) == "string" and (.stage | length) > 0
                    ' "$status_file" >/dev/null
                fi
                if [[ "$(/usr/bin/jq -er \
                        '.core_sha256sums_verified | booleans | tostring' \
                        "$status_file")" == true ]]; then
                    test -f "$verified_core/manifest.json"
                    failed_core_schema="$(/usr/bin/jq -er '.schema | strings' \
                        "$verified_core/manifest.json")"
                    if [[ "$evidence_generation" == failed-v1 ]]; then
                        test "$failed_core_schema" = \
                            leopard2-v17-gfni-main-core-manifest/v3 || \
                        test "$failed_core_schema" = \
                            leopard2-v17-gfni-main-passive-core-manifest/v4 || \
                        test "$failed_core_schema" = \
                            leopard2-v17-gfni-main-failed-core-manifest/v1
                    else
                        test "$failed_core_schema" = \
                            leopard2-v18-gfni-main-passive-core-manifest/v1 || \
                        test "$failed_core_schema" = \
                            leopard2-v18-gfni-main-passive-failed-core-manifest/v1
                        if [[ "$failed_core_schema" == \
                                leopard2-v18-gfni-main-passive-core-manifest/v1 ]]; then
                            /usr/bin/jq -e '
                                .stage as $stage |
                                ["seal_core", "independent_postseal_audit",
                                 "publish_not_promoted_envelope"] |
                                index($stage) != null
                            ' "$status_file" >/dev/null
                            verify_v18_complete_core_claim_bindings \
                                "$verified_core"
                        else
                            /usr/bin/jq -e \
                                '.stage == "verify_and_seal_failed_campaign"' \
                                "$status_file" >/dev/null
                            verify_v18_failed_core_claim_bindings \
                                "$verified_core"
                        fi
                        /usr/bin/jq -e --slurpfile status "$status_file" '
                            .acquisition_generation == "passive-v2" and
                            (.attempt | type) == "number" and
                            .attempt == (.attempt | floor) and
                            .attempt >= 1 and .attempt <= 3 and
                            .attempt == $status[0].attempt and
                            .attempt_budget == 3 and
                            .attempt_budget == $status[0].attempt_budget and
                            .attempt_lineage_sha256 ==
                                $status[0].attempt_lineage_sha256 and
                            .baseline_commit ==
                                "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
                        ' "$verified_core/manifest.json" >/dev/null
                    fi
                    test "$(/usr/bin/jq -er '.source_commit | strings' \
                        "$verified_core/manifest.json")" = \
                        "$(/usr/bin/jq -er '.source_commit | strings' \
                            "$status_file")"
                    test "$(/usr/bin/jq -er '.source_tree | strings' \
                        "$verified_core/manifest.json")" = \
                        "$(/usr/bin/jq -er '.source_tree | strings' \
                            "$status_file")"
                fi
            fi
        fi
    fi
    /usr/bin/printf 'authoritative v17/v18 envelope verified: %s\n' \
        "$verified_envelope"
}

if [[ $# -eq 1 && $1 == --self-test-passive-contract ]]; then
    test "$passive_mode" = false
    passive_contract_self_test
    exit 0
fi

if [[ $# -eq 4 && $1 == --verify-v18-complete-core-semantics &&
      $2 == /* &&
      $3 == /tmp/leopard-v18-complete-core-replay.*/NOT_PROMOTED.json &&
      $4 == /tmp/leopard-v18-complete-core-replay.*/postseal-audit.json ]]; then
    test "$passive_mode" = false
    verify_v18_complete_core_claim_bindings_preflight "$2/core"
    verify_envelope "$2" "$3" "$4"
    exit 0
fi

if [[ $# -eq 2 && $1 == --verify && $2 == /* ]]; then
    if [[ "$passive_mode" == true ]]; then
        /usr/bin/printf \
            '%s cannot be combined with --verify; generation comes from the sealed envelope\n' \
            --passive-shared-host >&2
        exit 2
    fi
    verify_envelope "$2"
    exit 0
fi

if [[ $# -ne 1 || $1 != /* ]]; then
    /usr/bin/printf 'usage: %s --passive-shared-host --attempt N --attempt-budget 3 /absolute/repository/.research/envelope\n' \
        "$0" >&2
    /usr/bin/printf '       %s --verify /absolute/repository/.research/envelope\n' \
        "$0" >&2
    /usr/bin/printf '       %s --self-test-passive-contract\n' "$0" >&2
    exit 2
fi
if [[ "$passive_mode" != true ]]; then
    /usr/bin/printf \
        'fresh active-v17 acquisition is disabled because the live runner now emits passive-only v18 evidence; retained v17 envelopes remain verifiable\n' >&2
    exit 2
fi
if ! validate_attempt_contract "$attempt" "$attempt_budget"; then
    /usr/bin/printf \
        'passive v18 acquisition requires --attempt N (1..3) and --attempt-budget 3\n' >&2
    exit 2
fi
envelope=$1
case "$envelope" in
    "$repo"/.research/*) ;;
    *)
        /usr/bin/printf 'artifact envelope must be below %s/.research\n' \
            "$repo" >&2
        exit 2
        ;;
esac
test "$(/usr/bin/readlink -m "$envelope")" = "$envelope"

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
lane=
commit=unknown
tree=unknown
attempt_lineage_json=
attempt_lineage_sha256=

next_stage()
{
    stage=$1
    /usr/bin/printf 'AUTHORITATIVE_STAGE %s\n' "$stage" | \
        /usr/bin/tee -a "$lane/wrapper-stage.log"
}

failure_record()
{
    local status=$?
    local lineage_ready=false
    local failed_verification=
    local failed_verification_status=1
    trap - ERR
    set +e
    if [[ -n "$lane" && -d "$envelope" ]]; then
        /usr/bin/chmod u+w "$envelope" 2>/dev/null
        if [[ ! -d "$lane" ]]; then
            /usr/bin/mkdir -m 0700 "$lane" 2>/dev/null
        fi
        /usr/bin/chmod u+w "$lane" 2>/dev/null
        if [[ ! -f "$lane/run-authoritative.sh" ]]; then
            /usr/bin/cp --reflink=never "$0" \
                "$lane/run-authoritative.sh" 2>/dev/null
            /usr/bin/chmod 0555 "$lane/run-authoritative.sh" 2>/dev/null
        fi
        if write_v18_attempt_lineage; then
            lineage_ready=true
        fi
    fi
    if [[ -n "$lane" && -d "$lane" && -w "$lane" &&
          -f "$lane/attempt-lineage.json" && "$lineage_ready" == true ]]; then
        /usr/bin/jq -n \
            --arg stage "$stage" \
            --argjson exit_status "$status" \
            --argjson attempt "$attempt" \
            --argjson attempt_budget "$attempt_budget" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
            '{schema:"leopard2-v18-gfni-main-passive-wrapper-failure/v1",status:"failed",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,promotion_passed:false,stage:$stage,exit_status:$exit_status,source_commit:$commit,source_tree:$tree}' \
            > "$lane/WRAPPER_FAILURE.json"
        (
            cd "$lane" || exit
            /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
                /usr/bin/sort -z | \
                /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
        )
        /usr/bin/find "$lane" -type f -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        /usr/bin/find "$lane" -type d -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)"
        /usr/bin/jq -n \
            --argjson exit_status "$status" \
            --arg stage "$stage" \
            --argjson attempt "$attempt" \
            --argjson attempt_budget "$attempt_budget" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
            --arg core_sha256sums_sha256 "$core_sha" \
            '{schema:"leopard2-v18-gfni-main-failed-envelope/v1",status:"failed",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,promotion_passed:false,campaign_exit_status:$exit_status,stage:$stage,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256}' \
            > "$envelope/FAILED.json"
        /usr/bin/printf '' > "$envelope/SHA256SUMS"
        write_tree_metadata_manifest "$envelope"
        (
            cd "$envelope" || exit
            /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
                /usr/bin/sort -z | \
                /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
        )
        /usr/bin/find "$envelope" -type f -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        /usr/bin/find "$envelope" -type d -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        failed_verification="$(/usr/bin/mktemp -d \
            /tmp/leopard-v18-wrapper-failed-verify.XXXXXX)"
        /usr/bin/bash "$lane/run-authoritative.sh" --verify "$envelope" \
            > "$failed_verification/verification.txt" \
            2> "$failed_verification/verification-stderr.log"
        failed_verification_status=$?
        if [[ "$failed_verification_status" -eq 0 ]]; then
            /usr/bin/printf \
                'sealed_failed_envelope=%s\nexternal_verification=%s\n' \
                "$envelope" "$failed_verification"
        else
            /usr/bin/printf \
                'sealed FAILED envelope did not verify; diagnostic=%s\n' \
                "$failed_verification" >&2
        fi
    fi
    /usr/bin/printf 'authoritative wrapper failed in stage %s with status %s\n' \
        "$stage" "$status" >&2
    exit "$status"
}

postseal_failure_record()
{
    local status=$?
    local failed_core_sha=
    local failed_core_verified=false
    local displaced_status=
    local failed_verification=
    local failed_verification_status=1
    trap - ERR
    set +e
    exec 8>&-
    if [[ -d "$envelope" && -f "$lane/SHA256SUMS" ]]; then
        /usr/bin/chmod u+w "$envelope" 2>/dev/null
        if [[ -e "$envelope/SHA256SUMS" ]]; then
            /usr/bin/chmod u+w "$envelope/SHA256SUMS" 2>/dev/null
        fi
        if (
            cd "$lane" || exit
            /usr/bin/sha256sum -c SHA256SUMS
        ) >/dev/null 2>&1; then
            failed_core_verified=true
        fi
        failed_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)"
        attempt_lineage_sha256="$(/usr/bin/sha256sum \
            "$lane/attempt-lineage.json" | /usr/bin/cut -d' ' -f1)"
        for displaced_status in \
            COMPLETED.json NOT_PROMOTED.json FAILED.json TREE-METADATA.json; do
            if [[ -e "$envelope/$displaced_status" ]]; then
                /usr/bin/mv "$envelope/$displaced_status" \
                    "$envelope/UNCOMMITTED-$displaced_status"
            fi
        done
        /usr/bin/jq -n \
            --argjson exit_status "$status" \
            --argjson core_sha256sums_verified "$failed_core_verified" \
            --argjson attempt "$attempt" \
            --argjson attempt_budget "$attempt_budget" \
            --arg stage "$stage" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg baseline_commit "$main_commit" \
            --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
            --arg core_sha256sums_sha256 "$failed_core_sha" \
            '{schema:"leopard2-v18-gfni-main-failed-envelope/v1",status:"failed",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,promotion_passed:false,campaign_exit_status:$exit_status,stage:$stage,source_commit:$commit,source_tree:$tree,baseline_commit:$baseline_commit,core_sha256sums_verified:$core_sha256sums_verified,core_sha256sums_sha256:$core_sha256sums_sha256}' \
            > "$envelope/FAILED.json"
        /usr/bin/printf '' > "$envelope/SHA256SUMS"
        write_tree_metadata_manifest "$envelope"
        (
            cd "$envelope" || exit
            /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
                /usr/bin/sort -z | \
                /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
        )
        /usr/bin/find "$envelope" -type f -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        /usr/bin/find "$envelope" -type d -perm /222 \
            -exec /usr/bin/chmod a-w {} +
        failed_verification="$(/usr/bin/mktemp -d \
            /tmp/leopard-v18-postseal-failed-verify.XXXXXX)"
        /usr/bin/bash "$lane/run-authoritative.sh" --verify "$envelope" \
            > "$failed_verification/verification.txt" \
            2> "$failed_verification/verification-stderr.log"
        failed_verification_status=$?
        if [[ "$failed_verification_status" -eq 0 ]]; then
            /usr/bin/printf \
                'sealed_failed_envelope=%s\nexternal_verification=%s\n' \
                "$envelope" "$failed_verification"
        else
            /usr/bin/printf \
                'postseal FAILED envelope did not verify; diagnostic=%s\n' \
                "$failed_verification" >&2
        fi
    fi
    /usr/bin/printf \
        'authoritative wrapper failed after core sealing in stage %s with status %s\n' \
        "$stage" "$status" >&2
    exit "$status"
}

preregistered_commit="$(/usr/bin/git -C "$repo" rev-parse HEAD)"
preregistered_tree="$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')"
[[ "$preregistered_commit" =~ ^[0-9a-f]{40}$ ]]
[[ "$preregistered_tree" =~ ^[0-9a-f]{40}$ ]]
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$repo/$relative_wrapper" -ef "$0"
commit=$preregistered_commit
tree=$preregistered_tree
expected_attempt_envelope="$repo/.research/leopard-79h/${preregistered_commit:0:7}-v18-passive-main-a${attempt}"
if [[ "$envelope" != "$expected_attempt_envelope" ]]; then
    /usr/bin/printf \
        'passive v18 attempt envelope must be exactly %s\n' \
        "$expected_attempt_envelope" >&2
    exit 2
fi
prior_attempt_records='[]'
if (( attempt > 1 )); then
    for ((prior_attempt = 1; prior_attempt < attempt; ++prior_attempt)); do
        prior_envelope="$repo/.research/leopard-79h/${preregistered_commit:0:7}-v18-passive-main-a${prior_attempt}"
        test -d "$prior_envelope"
        prior_sha_before="$(/usr/bin/sha256sum \
            "$prior_envelope/SHA256SUMS" | /usr/bin/cut -d' ' -f1)"
        verify_envelope "$prior_envelope"
        prior_sha_after="$(/usr/bin/sha256sum \
            "$prior_envelope/SHA256SUMS" | /usr/bin/cut -d' ' -f1)"
        test "$prior_sha_after" = "$prior_sha_before"
        if [[ ! -f "$prior_envelope/FAILED.json" ]]; then
            /usr/bin/printf \
                'refusing attempt %s because prior attempt %s produced valid observational evidence\n' \
                "$attempt" "$prior_attempt" >&2
            exit 2
        fi
        if v18_envelope_has_observational_output "$prior_envelope"; then
            /usr/bin/printf \
                'refusing attempt %s because prior attempt %s reached observational output before its later failure\n' \
                "$attempt" "$prior_attempt" >&2
            exit 2
        fi
        /usr/bin/jq -e --argjson prior_attempt "$prior_attempt" '
            .schema == "leopard2-v18-gfni-main-failed-envelope/v1" and
            .acquisition_generation == "passive-v2" and
            .attempt == $prior_attempt and .attempt_budget == 3
        ' "$prior_envelope/FAILED.json" >/dev/null
        test "$(/usr/bin/jq -er '.source_commit | strings' \
            "$prior_envelope/FAILED.json")" = "$preregistered_commit"
        test "$(/usr/bin/jq -er '.source_tree | strings' \
            "$prior_envelope/FAILED.json")" = "$preregistered_tree"
        /usr/bin/cmp "$prior_envelope/core/run-authoritative.sh" \
            "$repo/$relative_wrapper"
        prior_attempt_records="$(/usr/bin/jq -cS -n \
            --argjson records "$prior_attempt_records" \
            --argjson prior_attempt "$prior_attempt" \
            --arg prior_envelope "$prior_envelope" \
            --arg prior_sha "$prior_sha_before" \
            '$records + [{attempt:$prior_attempt,envelope:$prior_envelope,terminal:"FAILED.json",terminal_schema:"leopard2-v18-gfni-main-failed-envelope/v1",envelope_sha256sums_sha256:$prior_sha}]')"
    done
fi
for ((future_attempt = attempt + 1; future_attempt <= attempt_budget; \
        ++future_attempt)); do
    future_envelope="$repo/.research/leopard-79h/${preregistered_commit:0:7}-v18-passive-main-a${future_attempt}"
    require_path_absent "$future_envelope"
done
attempt_lineage_json="$(/usr/bin/jq -cS -n \
    --argjson attempt "$attempt" --argjson attempt_budget "$attempt_budget" \
    --arg commit "$preregistered_commit" --arg tree "$preregistered_tree" \
    --argjson prior_attempts "$prior_attempt_records" \
    '{schema:"leopard2-v18-gfni-main-attempt-lineage/v1",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,source_commit:$commit,source_tree:$tree,prior_attempts:$prior_attempts}')"
if ! require_path_absent "$envelope"; then
    /usr/bin/printf 'refusing to reuse artifact envelope: %s\n' \
        "$envelope" >&2
    exit 2
fi
/usr/bin/mkdir -p "$(/usr/bin/dirname "$envelope")"
test "$(/usr/bin/readlink -f "$(/usr/bin/dirname "$envelope")")" = \
    "$(/usr/bin/dirname "$envelope")"
/usr/bin/mkdir -m 0700 "$envelope"
lane="$envelope/core"
trap failure_record ERR
test "$(/usr/bin/readlink -f "$envelope")" = "$envelope"
/usr/bin/mkdir -m 0700 "$lane"
write_v18_attempt_lineage
/usr/bin/cp --reflink=never "$0" "$lane/run-authoritative.sh"
/usr/bin/chmod 0555 "$lane/run-authoritative.sh"
test ! -L "$lane/run-authoritative.sh"
test "$(/usr/bin/stat -c %h "$lane/run-authoritative.sh")" = 1

next_stage canonical_lock
exec 9>> "$lock"
/usr/bin/flock -n 9
test -f "$lock"
test ! -L "$lock"
test "$(/usr/bin/stat -c %h "$lock")" = 1

next_stage pre_tool_closure
/usr/bin/sha256sum "${tools_to_hash[@]}" > "$lane/pre-tool-SHA256SUMS"
/usr/bin/python3 -I -S -B -c \
    'import json,sys; assert sys.flags.isolated == 1 and sys.flags.no_site == 1 and sys.flags.ignore_environment == 1 and sys.dont_write_bytecode; print(json.dumps({"executable":sys.executable,"flags":{"dont_write_bytecode":sys.dont_write_bytecode,"ignore_environment":sys.flags.ignore_environment,"isolated":sys.flags.isolated,"no_site":sys.flags.no_site,"no_user_site":sys.flags.no_user_site}},sort_keys=True,separators=(",",":")))' \
    > "$lane/python-isolation.json"

next_stage build_order_normalizer_protocol_self_test
install_build_order_normalizer "$lane/build-order-normalizer.py"
normalizer_reconstruction_root="$(/usr/bin/mktemp -d \
    /tmp/leopard-v17-normalizer-reconstruction.XXXXXX)"
install_build_order_normalizer \
    "$normalizer_reconstruction_root/build-order-normalizer.py"
/usr/bin/cmp "$lane/build-order-normalizer.py" \
    "$normalizer_reconstruction_root/build-order-normalizer.py"
/usr/bin/python3 -I -S -B "$lane/build-order-normalizer.py" self-test \
    > "$lane/build-order-normalizer-self-test.log" 2>&1
normalized_evidence_closure_self_test \
    > "$lane/normalized-evidence-closure-self-test.log" 2>&1
candidate_commit_binding_self_test \
    > "$lane/candidate-commit-binding-self-test.log" 2>&1
/usr/bin/jq -n \
    '{schema:"leopard2-v17-build-order-normalizer-reconstruction/v3",embedded_wrapper_reconstruction_byte_identical:true,self_test_passed:true,normalized_evidence_closure_self_test_passed:true,candidate_commit_binding_self_test_passed:true,timing_performed:false,profiles:["candidate-tests","candidate-timing"],profile_sha256s:{"candidate-tests":"bead85b797c0e9966d4f09ba5b18b48739a350f550f831d1afa5bd2ae68dfa69","candidate-timing":"d2a07561d2e5919051c2b7d7464f6aa844752e42e92099190feeea004f2716b3"},closed_world_profile_templates:true,normalized_evidence_exact_file_count:11,template_path_tokens:{build:"${BUILD}",commit:"${COMMIT}",source:"${SOURCE}"}}' \
    > "$lane/build-order-normalizer-reconstruction.json"

next_stage tree_metadata_protocol_self_test
tree_self_test="$lane/tree-metadata-protocol-self-test"
/usr/bin/mkdir -m 0700 "$tree_self_test"
/usr/bin/mkdir -m 0700 "$tree_self_test/child"
/usr/bin/printf 'canonical-tree-protocol-v1\n' \
    > "$tree_self_test/child/payload.txt"
/usr/bin/printf '' > "$tree_self_test/SHA256SUMS"
write_tree_metadata_manifest "$tree_self_test"
(
    cd "$tree_self_test"
    /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
        /usr/bin/sort -z | \
        /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
)
/usr/bin/find "$tree_self_test" -type f -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$tree_self_test" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
verify_sealed_tree "$tree_self_test"
verify_tree_metadata_manifest "$tree_self_test"
/usr/bin/jq -n \
    '{schema:"leopard2-authoritative-tree-metadata-self-test/v1",passed:true,timing_performed:false}' \
    > "$lane/tree-metadata-protocol-self-test.json"

next_stage clean_source_preflight
commit="$(/usr/bin/git -C "$repo" rev-parse HEAD)"
tree="$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')"
test "$commit" = "$preregistered_commit"
test "${#commit}" = 40
test "${#tree}" = 40
test "$(/usr/bin/git -C "$repo" rev-parse --show-toplevel)" = "$repo"
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/readlink -f "$repo/$relative_wrapper")" = \
    "$(/usr/bin/readlink -f "$0")"
test "$(/usr/bin/grep -Fxc \
    'static std::atomic<uint32_t> g_auto_gf16_gfni_encode_mode(1U);' \
    "$repo/leopard2.cpp")" = 1
test "$(/usr/bin/grep -Fc 'g_auto_gf16_gfni_encode_mode(' \
    "$repo/leopard2.cpp")" = 1
candidate_submodule_commit="$(/usr/bin/git -C "$repo" ls-tree \
    "$commit" sse2neon | /usr/bin/cut -d' ' -f3 | \
    /usr/bin/cut -f1)"
baseline_submodule_commit="$(/usr/bin/git -C "$repo" ls-tree \
    "$main_commit" sse2neon | /usr/bin/cut -d' ' -f3 | \
    /usr/bin/cut -f1)"
test "${#candidate_submodule_commit}" = 40
test "$candidate_submodule_commit" = "$baseline_submodule_commit"
test "$(/usr/bin/git -C "$repo/sse2neon" rev-parse HEAD)" = \
    "$candidate_submodule_commit"
require_empty_output /usr/bin/git -C "$repo/sse2neon" status \
    --porcelain=v1 --untracked-files=normal
/usr/bin/printf '%s\n' "$commit" > "$lane/source-commit.txt"
/usr/bin/printf '%s\n' "$tree" > "$lane/source-tree.txt"
/usr/bin/printf '%s\n' "$main_commit" > "$lane/leopard1-commit.txt"
/usr/bin/printf '%s\n' \
    "$(/usr/bin/cat /sys/devices/system/cpu/cpu52/topology/thread_siblings_list)" \
    > "$lane/cpu52-thread-siblings.txt"
test "$(/usr/bin/cat "$lane/cpu52-thread-siblings.txt")" = 52,116
/usr/bin/python3 -I -S -B -c '
import json, os
allowed = sorted(os.sched_getaffinity(0))
if not (52 in allowed and 116 in allowed and
        len(set(allowed) - {52, 116}) > 0):
    raise SystemExit("wrapper launch affinity lacks the reserved pair or housekeeping")
print(json.dumps({"allowed_cpus": allowed}, sort_keys=True, separators=(",", ":")))
' > "$lane/wrapper-launch-affinity.json"

next_stage fresh_source_and_build_roots
work_root="$(/usr/bin/mktemp -d /tmp/leopard-v17-gfni-main.XXXXXX)"
candidate_source="$work_root/candidate-source"
baseline_source="$work_root/leopard1-source"
candidate_build="$work_root/candidate-build"
baseline_build="$work_root/baseline-build"
candidate_test_build="$work_root/candidate-test-build"
first_candidate_test_build="$work_root/first-candidate-test-build"
first_candidate_build="$work_root/first-candidate-build"
first_baseline_build="$work_root/first-baseline-build"
/usr/bin/jq -n \
    --arg work_root "$work_root" \
    --arg candidate_source "$candidate_source" \
    --arg baseline_source "$baseline_source" \
    --arg candidate_build "$candidate_build" \
    --arg baseline_build "$baseline_build" \
    --arg candidate_test_build "$candidate_test_build" \
    --arg first_candidate_test_build "$first_candidate_test_build" \
    --arg first_candidate_build "$first_candidate_build" \
    --arg first_baseline_build "$first_baseline_build" \
    '{work_root:$work_root,candidate_source:$candidate_source,baseline_source:$baseline_source,candidate_build:$candidate_build,baseline_build:$baseline_build,candidate_test_build:$candidate_test_build,first_candidate_test_build:$first_candidate_test_build,first_candidate_build:$first_candidate_build,first_baseline_build:$first_baseline_build,replay_policy:"two clean timing builds and two clean candidate correctness builds at identical absolute source/build paths; each first build retained before the second configure"}' \
    > "$lane/build-paths.json"

next_stage detached_source_clones
/usr/bin/git clone --no-hardlinks --no-checkout "$repo" "$candidate_source" \
    > "$lane/candidate-clone.log" 2>&1
/usr/bin/git -C "$candidate_source" checkout --detach "$commit" \
    > "$lane/candidate-checkout.log" 2>&1
/usr/bin/git -C "$candidate_source" config submodule.sse2neon.url \
    "$repo/sse2neon"
/usr/bin/git -c protocol.file.allow=always -C "$candidate_source" \
    submodule update --init --recursive \
    > "$lane/candidate-submodule.log" 2>&1
/usr/bin/git clone --no-hardlinks --no-checkout "$repo" "$baseline_source" \
    > "$lane/baseline-clone.log" 2>&1
/usr/bin/git -C "$baseline_source" checkout --detach "$main_commit" \
    > "$lane/baseline-checkout.log" 2>&1
/usr/bin/git -C "$baseline_source" config submodule.sse2neon.url \
    "$repo/sse2neon"
/usr/bin/git -c protocol.file.allow=always -C "$baseline_source" \
    submodule update --init --recursive \
    > "$lane/baseline-submodule.log" 2>&1
test "$(/usr/bin/git -C "$candidate_source" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$candidate_source" rev-parse 'HEAD^{tree}')" = "$tree"
test "$(/usr/bin/git -C "$baseline_source" rev-parse HEAD)" = "$main_commit"
test "$(/usr/bin/git -C "$candidate_source/sse2neon" rev-parse HEAD)" = \
    "$candidate_submodule_commit"
test "$(/usr/bin/git -C "$baseline_source/sse2neon" rev-parse HEAD)" = \
    "$baseline_submodule_commit"
test -z "$(/usr/bin/git -C "$candidate_source" symbolic-ref -q HEAD || true)"
test -z "$(/usr/bin/git -C "$baseline_source" symbolic-ref -q HEAD || true)"
require_empty_output /usr/bin/git -C "$candidate_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
require_empty_output /usr/bin/git -C "$baseline_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/grep -Fxc \
    'static std::atomic<uint32_t> g_auto_gf16_gfni_encode_mode(1U);' \
    "$candidate_source/leopard2.cpp")" = 1

next_stage canonical_source_archives
write_canonical_git_archive "$repo" "$commit" candidate-source \
    "$lane/candidate-source.tar"
write_canonical_git_archive "$candidate_source" "$commit" candidate-source \
    "$work_root/replayed-candidate-source.tar"
/usr/bin/cmp "$lane/candidate-source.tar" \
    "$work_root/replayed-candidate-source.tar"
write_canonical_git_archive "$repo" "$main_commit" leopard1-source \
    "$lane/leopard1-source.tar"
write_canonical_git_archive "$baseline_source" "$main_commit" \
    leopard1-source "$work_root/replayed-leopard1-source.tar"
/usr/bin/cmp "$lane/leopard1-source.tar" \
    "$work_root/replayed-leopard1-source.tar"
write_canonical_git_archive "$repo/sse2neon" \
    "$candidate_submodule_commit" sse2neon-source \
    "$lane/sse2neon-source.tar"
write_canonical_git_archive "$candidate_source/sse2neon" \
    "$candidate_submodule_commit" sse2neon-source \
    "$work_root/replayed-candidate-sse2neon-source.tar"
write_canonical_git_archive "$baseline_source/sse2neon" \
    "$baseline_submodule_commit" sse2neon-source \
    "$work_root/replayed-baseline-sse2neon-source.tar"
/usr/bin/cmp "$lane/sse2neon-source.tar" \
    "$work_root/replayed-candidate-sse2neon-source.tar"
/usr/bin/cmp "$lane/sse2neon-source.tar" \
    "$work_root/replayed-baseline-sse2neon-source.tar"

next_stage freeze_committed_controllers
/usr/bin/cmp "$candidate_source/$relative_wrapper" \
    "$lane/run-authoritative.sh"
/usr/bin/cp --reflink=never "$candidate_source/$relative_auditor" \
    "$lane/audit_v17_gfni_main_compare.py"
/usr/bin/cp --reflink=never "$candidate_source/$relative_runner" \
    "$lane/run_abba.py"
/usr/bin/cp --reflink=never "$candidate_source/$relative_git_capture" \
    "$lane/git_capture.py"
/usr/bin/cp --reflink=never "$candidate_source/$relative_helper" \
    "$lane/balanced_evidence_common.py"
/usr/bin/cp --reflink=never "$candidate_source/$relative_build_provenance" \
    "$lane/leopard2_build_provenance.py"
/usr/bin/cp --reflink=never "$candidate_source/$relative_supervisor" \
    "$lane/leopard2_affinity_supervisor.py"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/cp --reflink=never "$candidate_source/$relative_census" \
        "$lane/passive_environment_census.py"
fi
/usr/bin/chmod 0555 \
    "$lane/audit_v17_gfni_main_compare.py" \
    "$lane/run_abba.py" \
    "$lane/leopard2_affinity_supervisor.py"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/chmod 0555 "$lane/passive_environment_census.py"
fi
/usr/bin/chmod 0444 \
    "$lane/git_capture.py" \
    "$lane/balanced_evidence_common.py" \
    "$lane/leopard2_build_provenance.py" \
    "$lane/candidate-source.tar" \
    "$lane/leopard1-source.tar" \
    "$lane/sse2neon-source.tar"
for frozen_controller in \
    "$lane/run-authoritative.sh" \
    "$lane/audit_v17_gfni_main_compare.py" \
    "$lane/run_abba.py" \
    "$lane/git_capture.py" \
    "$lane/balanced_evidence_common.py" \
    "$lane/leopard2_build_provenance.py" \
    "$lane/build-order-normalizer.py" \
    "$lane/leopard2_affinity_supervisor.py"; do
    test ! -L "$frozen_controller"
    test "$(/usr/bin/stat -c %h "$frozen_controller")" = 1
done
if [[ "$passive_mode" == true ]]; then
    test ! -L "$lane/passive_environment_census.py"
    test "$(/usr/bin/stat -c %h \
        "$lane/passive_environment_census.py")" = 1
fi

candidate_configure=(
    -G "Unix Makefiles"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    -DLEO2_BUILD_TESTS=OFF
    -DLEO2_BUILD_BENCHMARKS=ON
    -DLEO2_BUILD_FUZZERS=OFF
    -DLEO2_ENABLE_CUDA=OFF
    -DLEO2_BACKEND_VARIANT=auto
    -DLEOPARD_ENABLE_GF8=ON
    -DLEOPARD_ENABLE_GF16=ON
    -DENABLE_OPENMP=ON
    -DLEO2_FLAG_MAVX512F=FALSE
    -DLEO2_FLAG_MAVX512BW=FALSE
    -DLEO2_FLAG_MAVX512VL=FALSE
    -DLEO2_FLAG_MPREFER_VECTOR_WIDTH_256=FALSE
    -DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git
)
candidate_test_configure=(
    -G "Unix Makefiles"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    -DLEO2_BUILD_TESTS=ON
    -DLEO2_BUILD_BENCHMARKS=ON
    -DLEO2_BUILD_FUZZERS=OFF
    -DLEO2_ENABLE_CUDA=OFF
    -DLEO2_BACKEND_VARIANT=auto
    -DLEOPARD_ENABLE_GF8=ON
    -DLEOPARD_ENABLE_GF16=ON
    -DENABLE_OPENMP=ON
    -DLEO2_FLAG_MAVX512F=FALSE
    -DLEO2_FLAG_MAVX512BW=FALSE
    -DLEO2_FLAG_MAVX512VL=FALSE
    -DLEO2_FLAG_MPREFER_VECTOR_WIDTH_256=FALSE
    -DLEO2_BENCHMARK_GIT_EXECUTABLE=/usr/bin/git
)
baseline_configure=(
    -G "Unix Makefiles"
    -DCMAKE_BUILD_TYPE=Release
    -DLEOPARD_MAIN_SOURCE_DIR="$baseline_source"
    -DLEO_MAIN_PURE_AVX2=OFF
)

assert_common_cache()
{
    local cache=$1/CMakeCache.txt
    /usr/bin/grep -Fx 'CMAKE_BUILD_TYPE:STRING=Release' "$cache"
    /usr/bin/grep -Fx 'CMAKE_CXX_COMPILER:FILEPATH=/usr/bin/c++' "$cache"
    /usr/bin/grep -Fx 'CMAKE_AR:FILEPATH=/usr/bin/ar' "$cache"
    /usr/bin/grep -Fx 'CMAKE_RANLIB:FILEPATH=/usr/bin/ranlib' "$cache"
    /usr/bin/grep -Fx 'CMAKE_MAKE_PROGRAM:FILEPATH=/usr/bin/gmake' "$cache"
    /usr/bin/grep -Fx 'CMAKE_GENERATOR:INTERNAL=Unix Makefiles' "$cache"
}

assert_candidate_cache()
{
    local cache=$1/CMakeCache.txt
    assert_common_cache "$1"
    for expected in \
        'ENABLE_OPENMP:BOOL=ON' \
        'LEO2_BACKEND_VARIANT:STRING=auto' \
        'LEO2_BUILD_BENCHMARKS:BOOL=ON' \
        'LEO2_BUILD_FUZZERS:BOOL=OFF' \
        'LEO2_ENABLE_CUDA:BOOL=OFF' \
        'LEOPARD_ENABLE_GF8:BOOL=ON' \
        'LEOPARD_ENABLE_GF16:BOOL=ON' \
        'LEO2_EXPERIMENT_HIGH_T16_Q2_B64_FUSED:BOOL=ON' \
        'LEO2_EXPERIMENT_SMALL_DUAL_LOCATOR_TERMS:BOOL=ON' \
        'LEO2_EXPERIMENT_SMALL_DUAL_REGULAR_FALLBACK:BOOL=ON' \
        'LEO2_FLAG_MGFNI:INTERNAL=1' \
        'LEO2_FLAG_MAVX512F:UNINITIALIZED=FALSE' \
        'LEO2_FLAG_MAVX512BW:UNINITIALIZED=FALSE' \
        'LEO2_FLAG_MAVX512VL:UNINITIALIZED=FALSE' \
        'LEO2_FLAG_MPREFER_VECTOR_WIDTH_256:UNINITIALIZED=FALSE' \
        'LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA:INTERNAL=leopard2-benchmark-build-configuration/v13'; do
        /usr/bin/grep -Fx "$expected" "$cache"
    done
}

assert_baseline_cache()
{
    local cache=$1/CMakeCache.txt
    assert_common_cache "$1"
    /usr/bin/grep -Fx "LEOPARD_MAIN_SOURCE_DIR:PATH=$baseline_source" "$cache"
    /usr/bin/grep -Fx 'LEO_MAIN_PURE_AVX2:BOOL=OFF' "$cache"
    /usr/bin/grep -Fx 'LEO_MAIN_HAS_MARCH_NATIVE:INTERNAL=1' "$cache"
}

candidate_test_selected_files=(
    bench_leopard2
    bench_leopard2_prevalidated_batch
    leopard2_auto_gf16_gfni_production_test
    leopard2_backend_failures_test
    leopard2_legacy_golden_test
    libleopard.a
    libleopard_test_hooks.a
    generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h
    generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt
)

copy_candidate_test_closure()
{
    local build=$1
    local destination=$2
    local source_file=
    local relative=
    /usr/bin/mkdir -m 0700 "$destination"
    for relative in \
        CMakeCache.txt compile_commands.json Makefile cmake_install.cmake \
        CTestTestfile.cmake CMakeFiles/Makefile.cmake CMakeFiles/Makefile2 \
        CMakeFiles/TargetDirectories.txt; do
        test -f "$build/$relative"
        /usr/bin/mkdir -p "$destination/$(/usr/bin/dirname "$relative")"
        /usr/bin/cp --reflink=never "$build/$relative" \
            "$destination/$relative"
    done
    while IFS= read -r -d '' source_file; do
        relative=${source_file#"$build"/}
        /usr/bin/mkdir -p "$destination/$(/usr/bin/dirname "$relative")"
        /usr/bin/cp --reflink=never "$source_file" "$destination/$relative"
    done < <(/usr/bin/find "$build/CMakeFiles" -type f \
        \( -name '*.o' -o -name link.txt -o -name flags.make \
           -o -name build.make -o -name DependInfo.cmake \
        \) -print0 | /usr/bin/sort -z)
    for relative in "${candidate_test_selected_files[@]}"; do
        test -f "$build/$relative"
        /usr/bin/mkdir -p \
            "$destination/$(/usr/bin/dirname "$relative")"
        /usr/bin/cp --reflink=never "$build/$relative" \
            "$destination/$relative"
    done
}

configure_and_build_candidate_tests()
{
    local generation=$1
    /usr/bin/cmake -S "$candidate_source" -B "$candidate_test_build" \
        "${candidate_test_configure[@]}" \
        > "$lane/$generation-candidate-test-configure.log" 2>&1
    assert_candidate_cache "$candidate_test_build" \
        >> "$lane/$generation-candidate-test-configure.log"
    /usr/bin/grep -Fx 'LEO2_BUILD_TESTS:BOOL=ON' \
        "$candidate_test_build/CMakeCache.txt" \
        >> "$lane/$generation-candidate-test-configure.log"
    /usr/bin/cmake --build "$candidate_test_build" --parallel 2 --target \
        bench_leopard2 \
        bench_leopard2_prevalidated_batch \
        leopard2_auto_gf16_gfni_production_test \
        leopard2_backend_failures_test \
        leopard2_legacy_golden_test \
        > "$lane/$generation-candidate-test-build.log" 2>&1
}

next_stage first_candidate_correctness_build
/usr/bin/mkdir -m 0700 "$lane/build-closure"
configure_and_build_candidate_tests first
copy_candidate_test_closure "$candidate_test_build" \
    "$lane/build-closure/live-candidate-tests"

next_stage retire_first_candidate_correctness_build
/usr/bin/mv "$candidate_test_build" "$first_candidate_test_build"
test ! -e "$candidate_test_build"

next_stage replay_candidate_correctness_build
configure_and_build_candidate_tests replay
copy_candidate_test_closure "$candidate_test_build" \
    "$lane/build-closure/replay-candidate-tests"

next_stage freeze_candidate_correctness_closure_inputs
/usr/bin/find \
    "$lane/build-closure/live-candidate-tests" \
    "$lane/build-closure/replay-candidate-tests" \
    -type f -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find \
    "$lane/build-closure/live-candidate-tests" \
    "$lane/build-closure/replay-candidate-tests" \
    -type d -perm /222 -exec /usr/bin/chmod a-w {} +

next_stage candidate_correctness_build_qualified_closure
compare_candidate_build_closures \
    "$lane/build-order-normalizer.py" candidate-tests \
    "$lane/build-closure/live-candidate-tests" \
    "$lane/build-closure/replay-candidate-tests" \
    "$candidate_test_build" "$candidate_source" "$commit" \
    "$lane/build-closure/candidate-test-order-normalization"
/usr/bin/cp --reflink=never \
    "$lane/build-closure/candidate-test-order-normalization/raw-tree.diff" \
    "$lane/candidate-test-build-byte-diff.txt"
/usr/bin/cp --reflink=never \
    "$lane/build-closure/candidate-test-order-normalization/raw-tree.diff.status" \
    "$lane/candidate-test-build-byte-diff.status"
(
    cd "$lane/build-closure/replay-candidate-tests"
    /usr/bin/sha256sum "${candidate_test_selected_files[@]}" \
        > "$lane/build-closure/candidate-test-selected-SHA256SUMS"
)
(
    cd "$lane/build-closure/live-candidate-tests"
    /usr/bin/sha256sum -c \
        "$lane/build-closure/candidate-test-selected-SHA256SUMS"
) > "$lane/candidate-test-live-selected-sha256-check.txt"
(
    cd "$lane/build-closure/replay-candidate-tests"
    /usr/bin/sha256sum -c \
        "$lane/build-closure/candidate-test-selected-SHA256SUMS"
) > "$lane/candidate-test-replay-selected-sha256-check.txt"
for reproduced_test_executable in \
    bench_leopard2_prevalidated_batch leopard2_backend_failures_test; do
    /usr/bin/cmp \
        "$candidate_test_build/$reproduced_test_executable" \
        "$lane/build-closure/live-candidate-tests/$reproduced_test_executable"
    test "$(/usr/bin/stat -c %i \
        "$candidate_test_build/$reproduced_test_executable")" != \
        "$(/usr/bin/stat -c %i \
        "$lane/build-closure/live-candidate-tests/$reproduced_test_executable")"
done

candidate_test_regex='^leopard2_(auto_gf16_gfni_production|portable_isa|portable_isa_registration|backend_failure_scalar_ff8_allocation|backend_failure_scalar_ff16_allocation|backend_failure_scalar_kat|backend_failure_ssse3_ff8_allocation|backend_failure_ssse3_ff16_allocation|backend_failure_ssse3_kat|backend_failure_avx2_ff8_allocation|backend_failure_avx2_ff16_allocation|backend_failure_avx2_kat|backend_failure_gfni_ff8_allocation|backend_failure_gfni_ff16_allocation|backend_failure_gfni_kat|backend_auto_gfni_encode_disabled_inert|backend_auto_gfni_encode_ineligible_inert|backend_auto_gfni_encode_kat_fallback|backend_auto_gfni_encode_ff16_allocation_fallback|backend_auto_gfni_encode_ff8_allocation_fallback|legacy_golden|benchmark_json_regression)$'
next_stage candidate_correctness_census
/usr/bin/ctest --test-dir "$candidate_test_build" -N \
    -R "$candidate_test_regex" --show-only=json-v1 \
    > "$lane/candidate-test-census.json"
/usr/bin/jq -e '
    [.tests[].name] == [
        "leopard2_auto_gf16_gfni_production",
        "leopard2_portable_isa",
        "leopard2_portable_isa_registration",
        "leopard2_backend_failure_scalar_ff8_allocation",
        "leopard2_backend_failure_scalar_ff16_allocation",
        "leopard2_backend_failure_scalar_kat",
        "leopard2_backend_failure_ssse3_ff8_allocation",
        "leopard2_backend_failure_ssse3_ff16_allocation",
        "leopard2_backend_failure_ssse3_kat",
        "leopard2_backend_failure_avx2_ff8_allocation",
        "leopard2_backend_failure_avx2_ff16_allocation",
        "leopard2_backend_failure_avx2_kat",
        "leopard2_backend_failure_gfni_ff8_allocation",
        "leopard2_backend_failure_gfni_ff16_allocation",
        "leopard2_backend_failure_gfni_kat",
        "leopard2_backend_auto_gfni_encode_disabled_inert",
        "leopard2_backend_auto_gfni_encode_ineligible_inert",
        "leopard2_backend_auto_gfni_encode_kat_fallback",
        "leopard2_backend_auto_gfni_encode_ff16_allocation_fallback",
        "leopard2_backend_auto_gfni_encode_ff8_allocation_fallback",
        "leopard2_legacy_golden",
        "leopard2_benchmark_json_regression"
    ]
' "$lane/candidate-test-census.json" >/dev/null

next_stage candidate_correctness_tests
/usr/bin/ctest --test-dir "$candidate_test_build" -j1 --output-on-failure \
    --no-tests=error -R "$candidate_test_regex" \
    > "$lane/candidate-focused-tests.log" 2>&1

next_stage candidate_posttest_build_closure
posttest_candidate_test_closure="$work_root/posttest-candidate-test-closure"
copy_candidate_test_closure "$candidate_test_build" \
    "$posttest_candidate_test_closure"
/usr/bin/diff -qr "$lane/build-closure/replay-candidate-tests" \
    "$posttest_candidate_test_closure" \
    > "$lane/candidate-posttest-build-byte-diff.txt"
(
    cd "$candidate_test_build"
    /usr/bin/sha256sum -c \
        "$lane/build-closure/candidate-test-selected-SHA256SUMS"
) > "$lane/candidate-posttest-selected-sha256-check.txt"

next_stage isolated_controller_self_tests
/usr/bin/python3 -I -S -B \
    "$candidate_source/$relative_runner" --help \
    > "$lane/runner-help.txt" 2> "$lane/runner-help-stderr.log"
/usr/bin/python3 -I -S -B \
    "$candidate_source/$relative_auditor" --self-test \
    > "$lane/auditor-self-test.log" 2>&1
if [[ "$passive_mode" == true ]]; then
    /usr/bin/python3 -I -S -B \
        "$candidate_source/$relative_census" self-test \
        > "$lane/passive-census-self-test.log" 2>&1
    test ! -e "$lane/supervisor-self-test.log"
else
    /usr/bin/python3 -I -S -B \
        "$candidate_source/$relative_supervisor" self-test \
        > "$lane/supervisor-self-test.log" 2>&1
fi

copy_timing_closure()
{
    local build=$1
    local destination=$2
    local role=$3
    local source_file=
    local relative=
    test "$role" = candidate || test "$role" = baseline
    /usr/bin/mkdir -m 0700 "$destination"
    for relative in \
        CMakeCache.txt compile_commands.json Makefile cmake_install.cmake \
        CMakeFiles/Makefile.cmake CMakeFiles/Makefile2 \
        CMakeFiles/TargetDirectories.txt; do
        test -f "$build/$relative"
        /usr/bin/mkdir -p "$destination/$(/usr/bin/dirname "$relative")"
        /usr/bin/cp --reflink=never "$build/$relative" \
            "$destination/$relative"
    done
    while IFS= read -r -d '' source_file; do
        relative=${source_file#"$build"/}
        /usr/bin/mkdir -p "$destination/$(/usr/bin/dirname "$relative")"
        /usr/bin/cp --reflink=never "$source_file" "$destination/$relative"
    done < <(/usr/bin/find "$build/CMakeFiles" -type f \
        \( -name '*.o' -o -name link.txt -o -name flags.make \
           -o -name build.make -o -name DependInfo.cmake \) -print0 | \
        /usr/bin/sort -z)
    if [[ "$role" == candidate ]]; then
        for relative in \
            bench_leopard2 libleopard.a \
            generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h \
            generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt; do
            test -f "$build/$relative"
            /usr/bin/mkdir -p \
                "$destination/$(/usr/bin/dirname "$relative")"
            /usr/bin/cp --reflink=never "$build/$relative" \
                "$destination/$relative"
        done
    else
        for relative in leopard_main_benchmark libleopard_main_exact.a; do
            test -f "$build/$relative"
            /usr/bin/cp --reflink=never "$build/$relative" \
                "$destination/$relative"
        done
    fi
}

configure_and_build_timing_pair()
{
    local generation=$1
    /usr/bin/cmake -S "$candidate_source" -B "$candidate_build" \
        "${candidate_configure[@]}" \
        > "$lane/$generation-candidate-configure.log" 2>&1
    assert_candidate_cache "$candidate_build" \
        >> "$lane/$generation-candidate-configure.log"
    /usr/bin/grep -Fx 'LEO2_BUILD_TESTS:BOOL=OFF' \
        "$candidate_build/CMakeCache.txt" \
        >> "$lane/$generation-candidate-configure.log"
    /usr/bin/cmake --build "$candidate_build" --parallel 2 --target \
        bench_leopard2 \
        > "$lane/$generation-candidate-build.log" 2>&1

    /usr/bin/cmake \
        -S "$candidate_source/experiments/leopard2/main_compare" \
        -B "$baseline_build" "${baseline_configure[@]}" \
        > "$lane/$generation-baseline-configure.log" 2>&1
    assert_baseline_cache "$baseline_build" \
        >> "$lane/$generation-baseline-configure.log"
    /usr/bin/cmake --build "$baseline_build" --parallel 2 --target \
        leopard_main_benchmark \
        > "$lane/$generation-baseline-build.log" 2>&1
}

next_stage first_canonical_timing_builds
configure_and_build_timing_pair first
copy_timing_closure "$candidate_build" \
    "$lane/build-closure/live-candidate" candidate
copy_timing_closure "$baseline_build" \
    "$lane/build-closure/live-baseline" baseline

next_stage retire_first_builds_without_path_change
/usr/bin/mv "$candidate_build" "$first_candidate_build"
/usr/bin/mv "$baseline_build" "$first_baseline_build"
test ! -e "$candidate_build"
test ! -e "$baseline_build"

next_stage replay_canonical_timing_builds
configure_and_build_timing_pair replay
copy_timing_closure "$candidate_build" \
    "$lane/build-closure/replay-candidate" candidate
copy_timing_closure "$baseline_build" \
    "$lane/build-closure/replay-baseline" baseline

next_stage freeze_candidate_timing_closure_inputs
/usr/bin/find \
    "$lane/build-closure/live-candidate" \
    "$lane/build-closure/replay-candidate" \
    -type f -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find \
    "$lane/build-closure/live-candidate" \
    "$lane/build-closure/replay-candidate" \
    -type d -perm /222 -exec /usr/bin/chmod a-w {} +

next_stage timing_live_replay_qualified_closure
compare_candidate_build_closures \
    "$lane/build-order-normalizer.py" candidate-timing \
    "$lane/build-closure/live-candidate" \
    "$lane/build-closure/replay-candidate" \
    "$candidate_build" "$candidate_source" "$commit" \
    "$lane/build-closure/candidate-timing-order-normalization"
/usr/bin/cp --reflink=never \
    "$lane/build-closure/candidate-timing-order-normalization/raw-tree.diff" \
    "$lane/candidate-build-byte-diff.txt"
/usr/bin/cp --reflink=never \
    "$lane/build-closure/candidate-timing-order-normalization/raw-tree.diff.status" \
    "$lane/candidate-build-byte-diff.status"
/usr/bin/diff -qr "$lane/build-closure/live-baseline" \
    "$lane/build-closure/replay-baseline" \
    > "$lane/baseline-build-byte-diff.txt"
/usr/bin/cmp "$candidate_build/bench_leopard2" \
    "$lane/build-closure/live-candidate/bench_leopard2"
/usr/bin/cmp "$candidate_build/libleopard.a" \
    "$lane/build-closure/live-candidate/libleopard.a"
/usr/bin/cmp "$baseline_build/leopard_main_benchmark" \
    "$lane/build-closure/live-baseline/leopard_main_benchmark"
/usr/bin/cmp "$baseline_build/libleopard_main_exact.a" \
    "$lane/build-closure/live-baseline/libleopard_main_exact.a"
test "$(/usr/bin/stat -c %i "$candidate_build/bench_leopard2")" != \
    "$(/usr/bin/stat -c %i \
        "$lane/build-closure/replay-candidate/bench_leopard2")"
test "$(/usr/bin/stat -c %i "$baseline_build/leopard_main_benchmark")" != \
    "$(/usr/bin/stat -c %i \
        "$lane/build-closure/replay-baseline/leopard_main_benchmark")"

next_stage baseline_test_census
/usr/bin/ctest --test-dir "$baseline_build" -N \
    -R '^leopard_main_(benchmark_smoke|compare_runner_selftest)$' \
    --show-only=json-v1 > "$lane/baseline-test-census.json"
/usr/bin/jq -e '
    [.tests[].name] == [
        "leopard_main_benchmark_smoke",
        "leopard_main_compare_runner_selftest"
    ]
' "$lane/baseline-test-census.json" >/dev/null

next_stage baseline_focused_tests
/usr/bin/ctest --test-dir "$baseline_build" -j1 --output-on-failure \
    --no-tests=error \
    -R '^leopard_main_(benchmark_smoke|compare_runner_selftest)$' \
    > "$lane/baseline-focused-tests.log" 2>&1

next_stage retain_build_analysis
/usr/bin/objdump -drwC -Mintel \
    "$candidate_build/CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNI.cpp.o" \
    > "$lane/build-closure/candidate-gfni-disassembly.txt"
/usr/bin/readelf -Ws \
    "$candidate_build/CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNI.cpp.o" \
    > "$lane/build-closure/candidate-gfni-symbols.txt"
/usr/bin/nm -S --size-sort \
    "$candidate_build/CMakeFiles/leopard2_backend_gfni.dir/Leopard2BackendGFNI.cpp.o" \
    > "$lane/build-closure/candidate-gfni-nm-sizes.txt"
/usr/bin/objdump -drwC -Mintel \
    "$baseline_build/CMakeFiles/leopard_main_exact.dir/${baseline_source#/}/LeopardFF16.cpp.o" \
    > "$lane/build-closure/leopard1-ff16-disassembly.txt"
/usr/bin/c++ --version > "$lane/build-closure/compiler-version.txt"
/usr/bin/cmake --version > "$lane/build-closure/cmake-version.txt"
/usr/bin/gmake --version > "$lane/build-closure/gmake-version.txt"
/usr/bin/git --version > "$lane/build-closure/git-version.txt"
/usr/bin/python3 --version > "$lane/build-closure/python-version.txt" 2>&1
/usr/bin/lscpu --json > "$lane/build-closure/lscpu.json"
/usr/bin/uname -a > "$lane/build-closure/uname.txt"

next_stage committed_controller_closure
/usr/bin/git -C "$repo" show "$commit:$relative_wrapper" \
    > "$lane/build-closure/committed-wrapper.sh"
/usr/bin/git -C "$repo" show "$commit:$relative_runner" \
    > "$lane/build-closure/committed-runner.py"
/usr/bin/git -C "$repo" show "$commit:$relative_git_capture" \
    > "$lane/build-closure/committed-git-capture.py"
/usr/bin/git -C "$repo" show "$commit:$relative_helper" \
    > "$lane/build-closure/committed-evidence-helper.py"
/usr/bin/git -C "$repo" show "$commit:$relative_build_provenance" \
    > "$lane/build-closure/committed-build-provenance.py"
/usr/bin/git -C "$repo" show "$commit:$relative_auditor" \
    > "$lane/build-closure/committed-auditor.py"
/usr/bin/git -C "$repo" show "$commit:$relative_supervisor" \
    > "$lane/build-closure/committed-supervisor.py"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/git -C "$repo" show "$commit:$relative_census" \
        > "$lane/build-closure/committed-passive-census.py"
fi
/usr/bin/cmp "$lane/run-authoritative.sh" \
    "$lane/build-closure/committed-wrapper.sh"
/usr/bin/cmp "$lane/run_abba.py" \
    "$lane/build-closure/committed-runner.py"
/usr/bin/cmp "$lane/git_capture.py" \
    "$lane/build-closure/committed-git-capture.py"
/usr/bin/cmp "$lane/balanced_evidence_common.py" \
    "$lane/build-closure/committed-evidence-helper.py"
/usr/bin/cmp "$lane/leopard2_build_provenance.py" \
    "$lane/build-closure/committed-build-provenance.py"
/usr/bin/cmp "$lane/audit_v17_gfni_main_compare.py" \
    "$lane/build-closure/committed-auditor.py"
/usr/bin/cmp "$lane/leopard2_affinity_supervisor.py" \
    "$lane/build-closure/committed-supervisor.py"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/cmp "$lane/passive_environment_census.py" \
        "$lane/build-closure/committed-passive-census.py"
fi
/usr/bin/cmp "$candidate_source/$relative_wrapper" \
    "$lane/run-authoritative.sh"
/usr/bin/cmp "$candidate_source/$relative_auditor" \
    "$lane/audit_v17_gfni_main_compare.py"

candidate_binary_hash="$(/usr/bin/sha256sum \
    "$candidate_build/bench_leopard2" | /usr/bin/cut -d' ' -f1)"
candidate_archive_hash="$(/usr/bin/sha256sum \
    "$candidate_build/libleopard.a" | /usr/bin/cut -d' ' -f1)"
baseline_binary_hash="$(/usr/bin/sha256sum \
    "$baseline_build/leopard_main_benchmark" | /usr/bin/cut -d' ' -f1)"
baseline_archive_hash="$(/usr/bin/sha256sum \
    "$baseline_build/libleopard_main_exact.a" | /usr/bin/cut -d' ' -f1)"
backend_failures_test_hash="$(/usr/bin/sha256sum \
    "$candidate_test_build/leopard2_backend_failures_test" | \
    /usr/bin/cut -d' ' -f1)"
prevalidated_batch_test_hash="$(/usr/bin/sha256sum \
    "$candidate_test_build/bench_leopard2_prevalidated_batch" | \
    /usr/bin/cut -d' ' -f1)"
candidate_test_selected_sha256sums_hash="$(/usr/bin/sha256sum \
    "$lane/build-closure/candidate-test-selected-SHA256SUMS" | \
    /usr/bin/cut -d' ' -f1)"
runner_hash="$(/usr/bin/sha256sum "$lane/run_abba.py" | \
    /usr/bin/cut -d' ' -f1)"
auditor_hash="$(/usr/bin/sha256sum \
    "$lane/audit_v17_gfni_main_compare.py" | /usr/bin/cut -d' ' -f1)"
supervisor_hash="$(/usr/bin/sha256sum \
    "$lane/leopard2_affinity_supervisor.py" | /usr/bin/cut -d' ' -f1)"
passive_census_hash=null
if [[ "$passive_mode" == true ]]; then
    passive_census_hash="$(/usr/bin/sha256sum \
        "$lane/passive_environment_census.py" | /usr/bin/cut -d' ' -f1)"
fi
wrapper_hash="$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | \
    /usr/bin/cut -d' ' -f1)"
build_order_normalizer_hash="$(/usr/bin/sha256sum \
    "$lane/build-order-normalizer.py" | /usr/bin/cut -d' ' -f1)"
candidate_timing_normalization_report_hash="$(/usr/bin/sha256sum \
    "$lane/build-closure/candidate-timing-order-normalization/normalized/report.json" | \
    /usr/bin/cut -d' ' -f1)"
candidate_timing_raw_tree_byte_identical="$(/usr/bin/jq -er \
    '.contract.raw_tree_file_bytes_identical | booleans | tostring' \
    "$lane/build-closure/candidate-timing-order-normalization/normalized/report.json")"
candidate_test_normalization_report_hash="$(/usr/bin/sha256sum \
    "$lane/build-closure/candidate-test-order-normalization/normalized/report.json" | \
    /usr/bin/cut -d' ' -f1)"
candidate_test_raw_tree_byte_identical="$(/usr/bin/jq -er \
    '.contract.raw_tree_file_bytes_identical | booleans | tostring' \
    "$lane/build-closure/candidate-test-order-normalization/normalized/report.json")"
candidate_source_archive_hash="$(/usr/bin/sha256sum \
    "$lane/candidate-source.tar" | /usr/bin/cut -d' ' -f1)"
baseline_source_archive_hash="$(/usr/bin/sha256sum \
    "$lane/leopard1-source.tar" | /usr/bin/cut -d' ' -f1)"
sse2neon_source_archive_hash="$(/usr/bin/sha256sum \
    "$lane/sse2neon-source.tar" | /usr/bin/cut -d' ' -f1)"

/usr/bin/jq -n \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg main_commit "$main_commit" \
    --arg candidate_source "$candidate_source" \
    --arg baseline_source "$baseline_source" \
    --arg candidate_build "$candidate_build" \
    --arg baseline_build "$baseline_build" \
    --arg candidate_test_build "$candidate_test_build" \
    --arg candidate_binary_sha256 "$candidate_binary_hash" \
    --arg candidate_archive_sha256 "$candidate_archive_hash" \
    --arg baseline_binary_sha256 "$baseline_binary_hash" \
    --arg baseline_archive_sha256 "$baseline_archive_hash" \
    --arg backend_failures_test_sha256 "$backend_failures_test_hash" \
    --arg prevalidated_batch_test_sha256 "$prevalidated_batch_test_hash" \
    --arg candidate_test_selected_sha256sums_sha256 \
        "$candidate_test_selected_sha256sums_hash" \
    --arg runner_sha256 "$runner_hash" \
    --arg auditor_sha256 "$auditor_hash" \
    --arg supervisor_sha256 "$supervisor_hash" \
    --arg wrapper_sha256 "$wrapper_hash" \
    --arg build_order_normalizer_sha256 "$build_order_normalizer_hash" \
    --arg candidate_timing_normalization_report_sha256 \
        "$candidate_timing_normalization_report_hash" \
    --argjson candidate_timing_raw_tree_byte_identical \
        "$candidate_timing_raw_tree_byte_identical" \
    --arg candidate_test_normalization_report_sha256 \
        "$candidate_test_normalization_report_hash" \
    --argjson candidate_test_raw_tree_byte_identical \
        "$candidate_test_raw_tree_byte_identical" \
    --arg candidate_source_archive_sha256 "$candidate_source_archive_hash" \
    --arg baseline_source_archive_sha256 "$baseline_source_archive_hash" \
    --arg sse2neon_commit "$candidate_submodule_commit" \
    --arg sse2neon_source_archive_sha256 "$sse2neon_source_archive_hash" \
    '{
        schema:"leopard2-v17-gfni-main-build-closure/v3",
        candidate:{
            commit:$commit,tree:$tree,source:$candidate_source,
            build:$candidate_build,
            profile:"standard AUTO Release; tests off; benchmarks on; GF8/GF16 on",
            binary_sha256:$candidate_binary_sha256,
            archive_sha256:$candidate_archive_sha256,
            source_archive_sha256:$candidate_source_archive_sha256,
            reproduction:{
                qualified_semantic_equal:true,
                normalizer_profile:"candidate-timing",
                raw_tree_file_bytes_identical:$candidate_timing_raw_tree_byte_identical,
                all_nonexception_file_bytes_identical:true,
                order_normalized_exception_paths:["compile_commands.json","CMakeFiles/Makefile2"],
                compile_commands_exact_entry_count:30,
                compile_commands_compiler_counts:{"/usr/bin/c++":30},
                makefile2_exact_validated_block_count:4,
                makefile2_exact_normalized_block_count:4,
                closed_world_template_sha256s:{
                    compile_commands:"6ad3a97ffe31f4eaf5c8f2ff5459dd310635e6e430181c460d9803f11cc5fb28",
                    compile_output_list:"a39c7e87cc0a506148faaa4023b11f96c2e26b22c103584a62b62ae91f2387b7",
                    makefile2:"5ab77166d931546b1cb998cfdf4eb5c7dd1850994b2eef65d8c5e627a3899a62"
                },
                node_census_count:131,
                node_census_sha256:"d3194ea1135fbdc0067cd3af72693cbf27d1e28a4734804b52bf5be1a3f37aca",
                profile_sha256:"d2a07561d2e5919051c2b7d7464f6aa844752e42e92099190feeea004f2716b3",
                normalized_evidence_exact_file_count:11,
                normalized_evidence_relative_canonical_sha256sums:true,
                normalization_report_sha256:$candidate_timing_normalization_report_sha256
            }
        },
        candidate_tests:{
            build:$candidate_test_build,
            profile:"standard AUTO Release; tests and benchmarks on; GF8/GF16 on",
            selected_files:["bench_leopard2","bench_leopard2_prevalidated_batch","leopard2_auto_gf16_gfni_production_test","leopard2_backend_failures_test","leopard2_legacy_golden_test","libleopard.a","libleopard_test_hooks.a","generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h","generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt"],
            selected_sha256sums_sha256:$candidate_test_selected_sha256sums_sha256,
            backend_failures_test_sha256:$backend_failures_test_sha256,
            prevalidated_batch_test_sha256:$prevalidated_batch_test_sha256,
            reproduction:{
                qualified_semantic_equal:true,
                normalizer_profile:"candidate-tests",
                raw_tree_file_bytes_identical:$candidate_test_raw_tree_byte_identical,
                all_nonexception_file_bytes_identical:true,
                order_normalized_exception_paths:["compile_commands.json","CMakeFiles/Makefile2"],
                compile_commands_exact_entry_count:170,
                compile_commands_compiler_counts:{"/usr/bin/c++":169,"/usr/bin/cc":1},
                single_c_record:{
                    source:"tests/leopard2/test_codec_options_abi.c",
                    output:"CMakeFiles/leopard2_codec_options_abi_test.dir/tests/leopard2/test_codec_options_abi.c.o"
                },
                makefile2_exact_validated_block_count:6,
                makefile2_exact_normalized_block_count:6,
                closed_world_template_sha256s:{
                    compile_commands:"dda169bcd458db8b17c09c1e8873c9342731a4c5639874849ad87d3c1003badb",
                    compile_output_list:"8aeea7120bbbe5e3f5e71059b29007b4fb3e01b153c55f0aa874c45de9e48574",
                    makefile2:"c56eb276850e7acbc5242e91c26b8fecf993040feeeb88b56596a189924a8479"
                },
                node_census_count:556,
                node_census_sha256:"2fa17fb9bb1458ea2bb7c361de077323e6f0a55359ba137eb445d231dd13d24b",
                profile_sha256:"bead85b797c0e9966d4f09ba5b18b48739a350f550f831d1afa5bd2ae68dfa69",
                normalized_evidence_exact_file_count:11,
                normalized_evidence_relative_canonical_sha256sums:true,
                normalization_report_sha256:$candidate_test_normalization_report_sha256
            },
            two_clean_builds_qualified_semantic_equal:true,
            two_clean_builds_raw_tree_file_bytes_identical:$candidate_test_raw_tree_byte_identical,
            posttest_raw_byte_identical:true,
            postcampaign_revalidation_required:true,
            complete_object_link_cache_closure:true
        },
        baseline:{
            commit:$main_commit,source:$baseline_source,build:$baseline_build,
            profile:"canonical Leopard1 native Release (-march=native; LEO_MAIN_PURE_AVX2=OFF)",
            binary_sha256:$baseline_binary_sha256,
            archive_sha256:$baseline_archive_sha256,
            source_archive_sha256:$baseline_source_archive_sha256,
            two_clean_builds_raw_byte_identical:true
        },
        sse2neon:{
            commit:$sse2neon_commit,archive_prefix:"sse2neon-source/",
            source_archive_sha256:$sse2neon_source_archive_sha256,
            reproduced_from_candidate_and_baseline_clones:true
        },
        generator:"Unix Makefiles",compiler:"/usr/bin/c++",
        build_reproduction:{
            candidate_timing_qualified_semantic_equal:true,
            candidate_timing_raw_tree_file_bytes_identical:$candidate_timing_raw_tree_byte_identical,
            candidate_tests_qualified_semantic_equal:true,
            candidate_tests_raw_tree_file_bytes_identical:$candidate_test_raw_tree_byte_identical,
            baseline_raw_byte_identical:true,
            identical_absolute_source_and_build_paths:true,
            objects_archives_binaries_link_recipes_cache_and_nonordering_generated_inputs_raw_byte_identical:true
        },
        controllers:{
            runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,
            supervisor_sha256:$supervisor_sha256,wrapper_sha256:$wrapper_sha256,
            build_order_normalizer_sha256:$build_order_normalizer_sha256,
            build_order_normalizer_origin:"embedded deterministic helper from committed wrapper"
        },
        audit_boundary:"independent campaign auditor covers retained campaign semantics; build-order qualification is separately recomputed by the frozen wrapper normalizer during durable verification",
        canonical_lock:"/tmp/leopard-gf8-authoritative.lock",
        python_controller:["/usr/bin/python3","-I","-S","-B"]
    }' \
    > "$lane/build-closure.json"

next_stage freeze_build_and_source_inputs
/usr/bin/find "$lane/build-closure" -type f -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$lane/build-closure" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/chmod a-w "$lane/build-closure.json"
/usr/bin/find "$candidate_build" "$baseline_build" \
    "$candidate_test_build" -type f -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$candidate_build" "$baseline_build" \
    "$candidate_test_build" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
# Keep Git metadata readable and available to the source-capture protocol, but
# make every worktree directory and tracked leaf non-writable so the supervised
# isolated interpreter cannot leave __pycache__ bytes beside committed source.
/usr/bin/find "$candidate_source" -path "$candidate_source/.git" -prune -o \
    \( -type f -o -type d \) -perm /222 -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$baseline_source" -path "$baseline_source/.git" -prune -o \
    \( -type f -o -type d \) -perm /222 -exec /usr/bin/chmod a-w {} +
require_empty_output /usr/bin/git -C "$candidate_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
require_empty_output /usr/bin/git -C "$baseline_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none

next_stage coordinator_reservation
reservation="$lane/cpu-reservation.json"
/usr/bin/python3 -I -S -B -c '
import json, os, sys
path = sys.argv[1]
payload = {
    "benchmark_cpu": 52,
    "nonce": os.urandom(32).hex(),
    "owner": f"v18 passive exact-main authoritative wrapper attempt {sys.argv[2]}/{sys.argv[3]}",
    "reserved_sibling": 116,
    "schema": "leopard2-cpu-reservation/v1",
    "status": "held",
}
data = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 0o600)
try:
    view = memoryview(data)
    while view:
        count = os.write(fd, view)
        if count <= 0:
            raise RuntimeError("short reservation write")
        view = view[count:]
    os.fsync(fd)
finally:
    os.close(fd)
parent = os.open(os.path.dirname(path), os.O_RDONLY | os.O_DIRECTORY)
try:
    os.fsync(parent)
finally:
    os.close(parent)
' "$reservation" "$attempt" "$attempt_budget"
test "$(/usr/bin/stat -c %a "$reservation")" = 600

original_allowed_csv=
housekeeping_csv=
if [[ "$passive_mode" == true ]]; then
    next_stage passive_controller_affinity
    original_allowed_csv="$(/usr/bin/jq -er \
        '.allowed_cpus | map(tostring) | join(",")' \
        "$lane/wrapper-launch-affinity.json")"
    housekeeping_csv="$(/usr/bin/jq -er \
        "$passive_housekeeping_jq" \
        "$lane/wrapper-launch-affinity.json")"
    current_allowed="$(/usr/bin/python3 -I -S -B -c '
import json, os, sys
print(json.dumps({"allowed_cpus": sorted(os.sched_getaffinity(int(sys.argv[1])))}, sort_keys=True, separators=(",", ":")))
' "$$")"
    test "$current_allowed" = \
        "$(/usr/bin/jq -cS . "$lane/wrapper-launch-affinity.json")"
    wrapper_children="$(<"/proc/$$/task/$$/children")"
    test -z "${wrapper_children//[[:space:]]/}"
    /usr/bin/taskset -pc "$housekeeping_csv" "$$" \
        > /dev/null
    /usr/bin/python3 -I -S -B -c '
import json, os, sys
source, output, original_csv, housekeeping_csv, wrapper_pid = sys.argv[1:]
def require(condition, message):
    if not condition:
        raise RuntimeError(message)
with open(source, "r", encoding="utf-8") as handle:
    before = json.load(handle)["allowed_cpus"]
after = sorted(os.sched_getaffinity(int(wrapper_pid)))
def parse(text):
    return [int(item) for item in text.split(",")]
require(before == parse(original_csv), "original affinity serialization differs")
require(after == parse(housekeeping_csv), "housekeeping affinity differs")
require(52 not in after and 116 not in after and bool(after),
        "housekeeping affinity includes the reserved pair or is empty")
payload = {
    "schema": "leopard2-v18-passive-controller-affinity/v1",
    "acquisition_generation": "passive-v2",
    "wrapper_pid": int(wrapper_pid),
    "before_allowed_cpus": before,
    "after_allowed_cpus": after,
    "runner_launch_allowed_cpus": before,
    "benchmark_cpu": 52,
    "reserved_sibling": 116,
    "affinity_mutation_scope": "wrapper-process-and-owned-descendants-only",
    "active_affinity_supervisor_executed": False,
}
data = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode() + b"\n"
fd = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, 0o600)
try:
    view = memoryview(data)
    while view:
        count = os.write(fd, view)
        require(count > 0, "short controller-affinity write")
        view = view[count:]
    os.fsync(fd)
finally:
    os.close(fd)
directory = os.open(os.path.dirname(output), os.O_RDONLY | os.O_DIRECTORY)
try:
    os.fsync(directory)
finally:
    os.close(directory)
' "$lane/wrapper-launch-affinity.json" "$lane/controller-affinity.json" \
        "$original_allowed_csv" "$housekeeping_csv" "$$"
    next_stage passive_environment_census_pre
    /usr/bin/python3 -I -S -B "$lane/passive_environment_census.py" capture \
        --phase pre --reserved-cpus 52,116 \
        --output "$lane/environment-census-pre.json"
    /usr/bin/python3 -I -S -B "$lane/passive_environment_census.py" verify \
        --input "$lane/environment-census-pre.json" --phase pre \
        --generation passive-v2
fi

if [[ "$passive_mode" == true ]]; then
    next_stage passive_v18_windowed_campaign
else
    next_stage supervised_campaign
fi
campaign_dir="$lane/campaign"
affinity_report="$lane/affinity-report.json"
campaign_command=()
if [[ "$passive_mode" == true ]]; then
    campaign_command=(
        /usr/bin/taskset -c "$original_allowed_csv"
        /usr/bin/python3 -I -S -B
        "$candidate_source/$relative_runner"
    )
else
    campaign_command=(
        /usr/bin/python3 -I -S -B
        "$candidate_source/$relative_supervisor"
        run
        --report "$affinity_report"
        --reserved-cpus 52,116
        --
        # Supervisor v14 requires the child Python source as direct argv[1].
        /usr/bin/python3
        "$candidate_source/$relative_runner"
    )
fi
campaign_command+=(
    run
    --baseline "$baseline_build/leopard_main_benchmark"
    --candidate "$candidate_build/bench_leopard2"
    --baseline-archive "$baseline_build/libleopard_main_exact.a"
    --candidate-archive "$candidate_build/libleopard.a"
    --baseline-build-dir "$baseline_build"
    --candidate-build-dir "$candidate_build"
    --baseline-source-root "$baseline_source"
    --candidate-source-root "$candidate_source"
    --candidate-commit "$commit"
    --baseline-native
    --candidate-mode auto
    --reservation-file "$reservation"
    --output "$campaign_dir"
    --cpu 52
    --reserved-sibling 116
    --taskset /usr/bin/taskset
    --ldd /usr/bin/ldd
    --preset v17-gfni-encode
    --reuse 8
    --iterations 9
    --warmup 2
    --timeout 600
)
/usr/bin/jq -n --args '$ARGS.positional' -- "${campaign_command[@]}" \
    > "$lane/campaign-command.json"
trap - ERR
set +e
"${campaign_command[@]}" > "$lane/campaign-summary.log" \
    2> "$lane/campaign-stderr.log"
campaign_status=$?
set -e
trap failure_record ERR

if [[ "$passive_mode" == true ]]; then
    next_stage passive_environment_census_post
    /usr/bin/python3 -I -S -B "$lane/passive_environment_census.py" capture \
        --phase post --reserved-cpus 52,116 \
        --output "$lane/environment-census-post.json"
fi
if [[ "$passive_mode" == true && "$campaign_status" -eq 0 ]]; then
    next_stage passive_environment_policy_v2
    /usr/bin/python3 -I -S -B "$lane/passive_environment_census.py" compare \
        --pre "$lane/environment-census-pre.json" \
        --post "$lane/environment-census-post.json" \
        --raw "$campaign_dir/raw.json" \
        --controller "$lane/controller-affinity.json" \
        --output "$lane/passive-environment-policy.json"
    next_stage passive_windowed_contamination_gate
    /usr/bin/jq -e '
        .schema == "leopard2-passive-shared-host-policy/v2" and
        .windowed_contamination.gated == true and
        .windowed_contamination.retained_window_count == 12 and
        .windowed_contamination.all_benchmark_cpu_excess_zero == true and
        .windowed_contamination.all_reserved_sibling_nonidle_zero == true and
        .windowed_contamination.windowed.benchmark_cpu_nonidle_excess_jiffies == 0 and
        .windowed_contamination.windowed.reserved_sibling_nonidle_jiffies == 0 and
        .outer_disclosure.gated == false
    ' "$lane/passive-environment-policy.json" >/dev/null
fi

if [[ "$campaign_status" -ne 0 ]]; then
    next_stage verify_and_seal_failed_campaign
    failure_verified=false
    failure_verify_status=1
    if [[ -f "$campaign_dir/failure.json" ]]; then
        trap - ERR
        set +e
        /usr/bin/python3 -I -S -B \
            "$candidate_source/$relative_runner" verify-failure \
            --failure "$campaign_dir/failure.json" \
            > "$lane/failure-verification.log" \
            2> "$lane/failure-verification-stderr.log"
        failure_verify_status=$?
        set -e
        trap failure_record ERR
        if [[ "$failure_verify_status" -eq 0 ]]; then
            failure_verified=true
        fi
    fi
    failure_sha=null
    if [[ -f "$campaign_dir/failure.json" ]]; then
        failure_sha="$(/usr/bin/sha256sum "$campaign_dir/failure.json" | \
            /usr/bin/cut -d' ' -f1)"
    fi
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --argjson failure_verify_status "$failure_verify_status" \
        --argjson failure_verified "$failure_verified" \
        --argjson attempt "$attempt" \
        --argjson attempt_budget "$attempt_budget" \
        --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg main_commit "$main_commit" \
        --arg candidate_binary_sha256 "$candidate_binary_hash" \
        --arg baseline_binary_sha256 "$baseline_binary_hash" \
        --arg failure_sha256 "$failure_sha" \
        '{schema:"leopard2-v18-gfni-main-passive-failed-core-manifest/v1",status:"failed",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,promotion_passed:false,campaign_exit_status:$campaign_exit_status,failure_verify_status:$failure_verify_status,failure_verified:$failure_verified,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,candidate_binary_sha256:$candidate_binary_sha256,baseline_binary_sha256:$baseline_binary_sha256,failure_sha256:($failure_sha256 | if . == "null" then null else . end),canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$lane/manifest.json"
    (
        cd "$lane"
        /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
            /usr/bin/sort -z | \
            /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    )
    trap postseal_failure_record ERR
    /usr/bin/find "$lane" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$lane" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    failed_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --argjson failure_verified "$failure_verified" \
        --argjson attempt "$attempt" \
        --argjson attempt_budget "$attempt_budget" \
        --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg core_sha256sums_sha256 "$failed_core_sha" \
        '{schema:"leopard2-v18-gfni-main-failed-envelope/v1",status:"failed",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,promotion_passed:false,campaign_exit_status:$campaign_exit_status,failure_verified:$failure_verified,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256}' \
        > "$envelope/FAILED.json"
    /usr/bin/printf '' > "$envelope/SHA256SUMS"
    write_tree_metadata_manifest "$envelope"
    (
        cd "$envelope"
        /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
            /usr/bin/sort -z | \
            /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    )
    /usr/bin/find "$envelope" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$envelope" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    failed_verification="$(/usr/bin/mktemp -d \
        /tmp/leopard-v17-gfni-main-failed-verify.XXXXXX)"
    (
        cd "$envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) > "$failed_verification/sha256-check.txt"
    "$lane/run-authoritative.sh" --verify "$envelope" \
        > "$failed_verification/verification.txt" \
        2> "$failed_verification/verification-stderr.log"
    /usr/bin/printf 'sealed_failed_envelope=%s\nexternal_verification=%s\n' \
        "$envelope" "$failed_verification"
    trap - ERR
    exit "$campaign_status"
fi

affinity_binding=
if [[ "$passive_mode" == true ]]; then
    next_stage passive_null_supervision_contract
    test "$(/usr/bin/jq -er '.supervision == null' \
        "$campaign_dir/raw.json")" = true
    test "$(/usr/bin/jq -er '.supervision == null' \
        "$campaign_dir/manifest.json")" = true
    for forbidden_affinity_artifact in \
        "$affinity_report" \
        "$lane/affinity-report.json.accepted.json" \
        "$lane/affinity-report.json.ambiguous" \
        "$lane/affinity-binding.json" \
        "$lane/affinity-report-verification.log" \
        "$lane/affinity-binding-create.log" \
        "$lane/affinity-binding-verification.log" \
        "$lane/supervisor-self-test.log" \
        "$lane/passive-auditor-self-test.log" \
        "$lane/passive-controller-taskset.log"; do
        test ! -e "$forbidden_affinity_artifact"
    done
else
    next_stage accepted_supervisor_report
    test -f "$affinity_report"
    test "$(/usr/bin/jq -er '.schema | strings' "$affinity_report")" = \
        leopard2-affinity-supervisor/v14
    test "$(/usr/bin/jq -er '.accepted | booleans | tostring' \
        "$affinity_report")" = true
    /usr/bin/python3 -I -S -B \
        "$candidate_source/$relative_supervisor" verify-report \
        --report "$affinity_report" \
        > "$lane/affinity-report-verification.log" 2>&1

    next_stage bind_affinity_to_manifest
    test -f "$campaign_dir/manifest.json"
    affinity_binding="$lane/affinity-binding.json"
    /usr/bin/python3 -I -S -B \
        "$candidate_source/$relative_supervisor" bind \
        --report "$affinity_report" \
        --manifest "$campaign_dir/manifest.json" \
        --output "$affinity_binding" \
        > "$lane/affinity-binding-create.log" 2>&1
    test "$(/usr/bin/jq -er '.schema | strings' "$affinity_binding")" = \
        leopard2-affinity-main-binding/v8
    /usr/bin/python3 -I -S -B \
        "$candidate_source/$relative_supervisor" verify-binding \
        --binding "$affinity_binding" \
        --manifest "$campaign_dir/manifest.json" \
        --manifest-sha256 \
        "$(/usr/bin/sha256sum "$campaign_dir/manifest.json" | \
            /usr/bin/cut -d' ' -f1)" \
        > "$lane/affinity-binding-verification.log" 2>&1
fi

next_stage producer_bundle_verification
producer_verify_args=(
    /usr/bin/python3 -I -S -B
    "$candidate_source/$relative_runner" verify
    --manifest "$campaign_dir/manifest.json"
)
if [[ "$passive_mode" == false ]]; then
    producer_verify_args+=(--affinity-binding "$affinity_binding")
fi
"${producer_verify_args[@]}" > "$lane/producer-verification.log" \
    2> "$lane/producer-verification-stderr.log"

next_stage independent_preseal_audit
audit_args=(
    /usr/bin/python3 -I -S -B
    "$lane/audit_v17_gfni_main_compare.py"
    --manifest "$campaign_dir/manifest.json"
)
if [[ "$passive_mode" == true ]]; then
    audit_args+=(--supervision windowed)
fi
"${audit_args[@]}" \
    --output "$lane/audit.json" \
    > "$lane/audit-summary.json" 2> "$lane/audit-stderr.log"
test "$(/usr/bin/jq -er '.audit_passed | booleans | tostring' \
    "$lane/audit.json")" = true
test "$(/usr/bin/jq -er '.timing_performed | booleans | tostring' \
    "$lane/audit.json")" = false
test "$(/usr/bin/jq -er '.benchmark_executed | booleans | tostring' \
    "$lane/audit.json")" = false
test "$(/usr/bin/jq -er \
    '.reporting_policy.combined_or_stacked_ratio_emitted | booleans | tostring' \
    "$lane/audit.json")" = false
if [[ "$passive_mode" == true ]]; then
    test "$(/usr/bin/jq -er '.schema | strings' "$lane/audit.json")" = \
        leopard2-main-compare-v18-passive-independent-audit/v1
    /usr/bin/jq -e '
        .audit_mode == "passive-windowed-shared-host-observation" and
        .acquisition_generation == "passive-v2" and
        .evidence_class == "passive-windowed-shared-host-observation/v1" and
        .affinity_supervisor_binding_verified == false and
        .promotion_eligible == false and .promotion_passed == false and
        .causal_performance_claim_eligible == false and
        .cpu_pair_exclusive == false and
        .isolation_claim.benchmark_cpu_exclusive_ownership_claimed == false and
        .isolation_claim.foreign_process_affinity_mutation_claimed == false and
        .isolation_claim.foreign_process_signalling_claimed == false and
        .isolation_claim.every_retained_window_benchmark_cpu_excess_zero == true and
        .isolation_claim.every_retained_window_reserved_sibling_nonidle_zero == true and
        .isolation_claim.whole_campaign_and_out_of_window_activity_gated == false and
        .isolation_claim.windowed_screen ==
            "per retained benchmark invocation" and
        .isolation_claim.out_of_window_activity_gated == false and
        .isolation_claim.reserved_sibling_zero_nonidle_jiffies_in_every_retained_window == true and
        .isolation_claim.promotion_eligible == false and
        .isolation_claim.promotion_passed == false and
        .windowed_observation.retained_window_count == 12 and
        .windowed_observation.windowed.benchmark_cpu_nonidle_excess_jiffies == 0 and
        .windowed_observation.windowed.reserved_sibling_nonidle_jiffies == 0 and
        .windowed_observation.out_of_window.gated == false and
        .gate_results.reservation_pair_lease_null_supervision_windowed_closure == true and
        .gate_results.benchmark_cpu_nonidle_excess_zero_in_every_retained_window == true and
        .gate_results.reserved_sibling_nonidle_zero_in_every_retained_window == true and
        .gate_results.out_of_window_activity_disclosed_not_gated == true and
        .gate_results.benchmark_cpu_exclusivity_not_claimed == true
    ' "$lane/audit.json" >/dev/null
fi

next_stage post_campaign_closure
postcampaign_candidate_test_closure=\
"$work_root/postcampaign-candidate-test-closure"
copy_candidate_test_closure "$candidate_test_build" \
    "$postcampaign_candidate_test_closure"
/usr/bin/diff -qr "$lane/build-closure/replay-candidate-tests" \
    "$postcampaign_candidate_test_closure" \
    > "$lane/candidate-postcampaign-build-byte-diff.txt"
test "$(/usr/bin/sha256sum \
    "$lane/build-closure/candidate-test-selected-SHA256SUMS" | \
    /usr/bin/cut -d' ' -f1)" = \
    "$candidate_test_selected_sha256sums_hash"
(
    cd "$candidate_test_build"
    /usr/bin/sha256sum -c \
        "$lane/build-closure/candidate-test-selected-SHA256SUMS"
) > "$lane/candidate-postcampaign-selected-sha256-check.txt"
/usr/bin/jq -n \
    --arg selected_sha256sums_sha256 \
        "$candidate_test_selected_sha256sums_hash" \
    --arg normalization_report_sha256 \
        "$candidate_test_normalization_report_hash" \
    --argjson two_clean_builds_raw_tree_file_bytes_identical \
        "$candidate_test_raw_tree_byte_identical" \
    '{schema:"leopard2-v17-gfni-main-candidate-test-temporal-closure/v2",selected_sha256sums_sha256:$selected_sha256sums_sha256,normalization_report_sha256:$normalization_report_sha256,two_clean_builds_qualified_semantic_equal:true,two_clean_builds_raw_tree_file_bytes_identical:$two_clean_builds_raw_tree_file_bytes_identical,posttest_byte_identical:true,postcampaign_byte_identical:true,canonical_test_build_frozen_during_campaign:true}' \
    > "$lane/candidate-test-temporal-closure.json"
/usr/bin/cmp "$candidate_build/bench_leopard2" \
    "$lane/build-closure/replay-candidate/bench_leopard2"
/usr/bin/cmp "$candidate_build/libleopard.a" \
    "$lane/build-closure/replay-candidate/libleopard.a"
/usr/bin/cmp "$baseline_build/leopard_main_benchmark" \
    "$lane/build-closure/replay-baseline/leopard_main_benchmark"
/usr/bin/cmp "$baseline_build/libleopard_main_exact.a" \
    "$lane/build-closure/replay-baseline/libleopard_main_exact.a"
test "$(/usr/bin/sha256sum "$candidate_build/bench_leopard2" | \
    /usr/bin/cut -d' ' -f1)" = "$candidate_binary_hash"
test "$(/usr/bin/sha256sum "$candidate_build/libleopard.a" | \
    /usr/bin/cut -d' ' -f1)" = "$candidate_archive_hash"
test "$(/usr/bin/sha256sum "$baseline_build/leopard_main_benchmark" | \
    /usr/bin/cut -d' ' -f1)" = "$baseline_binary_hash"
test "$(/usr/bin/sha256sum "$baseline_build/libleopard_main_exact.a" | \
    /usr/bin/cut -d' ' -f1)" = "$baseline_archive_hash"
require_empty_output /usr/bin/find "$candidate_build" "$baseline_build" \
    "$candidate_test_build" \
    -type f -perm /222 -print -quit
require_empty_output /usr/bin/find "$candidate_build" "$baseline_build" \
    "$candidate_test_build" \
    -type d -perm /222 -print -quit
/usr/bin/cmp "$candidate_source/$relative_wrapper" \
    "$lane/run-authoritative.sh"
/usr/bin/cmp "$candidate_source/$relative_runner" "$lane/run_abba.py"
/usr/bin/cmp "$candidate_source/$relative_git_capture" "$lane/git_capture.py"
/usr/bin/cmp "$candidate_source/$relative_helper" \
    "$lane/balanced_evidence_common.py"
/usr/bin/cmp "$candidate_source/$relative_build_provenance" \
    "$lane/leopard2_build_provenance.py"
/usr/bin/cmp "$candidate_source/$relative_auditor" \
    "$lane/audit_v17_gfni_main_compare.py"
/usr/bin/cmp "$candidate_source/$relative_supervisor" \
    "$lane/leopard2_affinity_supervisor.py"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/cmp "$candidate_source/$relative_census" \
        "$lane/passive_environment_census.py"
fi
test "$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | \
    /usr/bin/cut -d' ' -f1)" = "$wrapper_hash"
test "$(/usr/bin/sha256sum "$lane/run_abba.py" | \
    /usr/bin/cut -d' ' -f1)" = "$runner_hash"
test "$(/usr/bin/sha256sum "$lane/audit_v17_gfni_main_compare.py" | \
    /usr/bin/cut -d' ' -f1)" = "$auditor_hash"
test "$(/usr/bin/sha256sum "$lane/leopard2_affinity_supervisor.py" | \
    /usr/bin/cut -d' ' -f1)" = "$supervisor_hash"
if [[ "$passive_mode" == true ]]; then
    test "$(/usr/bin/sha256sum "$lane/passive_environment_census.py" | \
        /usr/bin/cut -d' ' -f1)" = "$passive_census_hash"
fi
test "$(/usr/bin/sha256sum "$lane/build-order-normalizer.py" | \
    /usr/bin/cut -d' ' -f1)" = "$build_order_normalizer_hash"
require_empty_output /usr/bin/git -C "$candidate_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
require_empty_output /usr/bin/git -C "$baseline_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/git -C "$candidate_source" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$candidate_source" rev-parse 'HEAD^{tree}')" = "$tree"
test "$(/usr/bin/git -C "$baseline_source" rev-parse HEAD)" = "$main_commit"
write_canonical_git_archive "$candidate_source" "$commit" candidate-source \
    "$work_root/post-candidate-source.tar"
write_canonical_git_archive "$baseline_source" "$main_commit" \
    leopard1-source "$work_root/post-leopard1-source.tar"
write_canonical_git_archive "$candidate_source/sse2neon" \
    "$candidate_submodule_commit" sse2neon-source \
    "$work_root/post-candidate-sse2neon-source.tar"
write_canonical_git_archive "$baseline_source/sse2neon" \
    "$baseline_submodule_commit" sse2neon-source \
    "$work_root/post-baseline-sse2neon-source.tar"
/usr/bin/cmp "$lane/candidate-source.tar" \
    "$work_root/post-candidate-source.tar"
/usr/bin/cmp "$lane/leopard1-source.tar" \
    "$work_root/post-leopard1-source.tar"
/usr/bin/cmp "$lane/sse2neon-source.tar" \
    "$work_root/post-candidate-sse2neon-source.tar"
/usr/bin/cmp "$lane/sse2neon-source.tar" \
    "$work_root/post-baseline-sse2neon-source.tar"
test "$(/usr/bin/git -C "$repo" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$repo" rev-parse 'HEAD^{tree}')" = "$tree"
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal
require_empty_output /usr/bin/git -C "$repo" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/git -C "$repo/sse2neon" rev-parse HEAD)" = \
    "$candidate_submodule_commit"
require_empty_output /usr/bin/git -C "$repo/sse2neon" status \
    --porcelain=v1 --untracked-files=normal
test "$(/usr/bin/sha256sum "$lane/sse2neon-source.tar" | \
    /usr/bin/cut -d' ' -f1)" = "$sse2neon_source_archive_hash"
/usr/bin/sha256sum "${tools_to_hash[@]}" > "$lane/post-tool-SHA256SUMS"
/usr/bin/cmp "$lane/pre-tool-SHA256SUMS" "$lane/post-tool-SHA256SUMS"

next_stage result_policy
test "$(/usr/bin/jq -er '.schema | strings' \
    "$campaign_dir/manifest.json")" = leopard2-main-compare-manifest/v18
test "$(/usr/bin/jq -er '.valid | booleans | tostring' \
    "$campaign_dir/manifest.json")" = true
performance_gate_passed="$(/usr/bin/jq -er '
    [
        .analysis["gf16-high-full"].encode
            .promotion_lower_bound_at_least_1_05,
        .analysis["gf16-high-full"].one_shot_encode
            .promotion_lower_bound_at_least_1_05
    ] |
    if all(.[]; type == "boolean") then all else error("invalid gates") end |
    tostring
' "$campaign_dir/manifest.json")"
if [[ "$passive_mode" == true ]]; then
    /usr/bin/jq --slurpfile audit "$lane/audit.json" \
        --slurpfile lineage "$lane/attempt-lineage.json" \
        --slurpfile policy "$lane/passive-environment-policy.json" \
        --argjson attempt "$attempt" --argjson attempt_budget "$attempt_budget" \
        --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
        --arg current_envelope "$envelope" '{
        schema:"leopard2-v18-gfni-main-passive-result-summary/v1",
        acquisition_generation:"passive-v2",
        attempt:$attempt,
        attempt_budget:$attempt_budget,
        attempt_lineage_sha256:$attempt_lineage_sha256,
        attempt_statement:("attempt \($attempt) of at most \($attempt_budget)"),
        sealed_attempt_envelopes:
            ([$lineage[0].prior_attempts[].envelope] + [$current_envelope]),
        active_affinity_supervisor_executed:false,
        evidence_class:"passive-windowed-shared-host-observation/v1",
        promotion_eligible:false,
        promotion_passed:false,
        causal_performance_claim_eligible:false,
        cpu_pair_exclusive:false,
        ratio_orientation:
            "exact_leopard1_native_time_over_leopard2_candidate_time",
        ratios_are_separate_correlated_and_must_not_be_multiplied:true,
        clustered_t_interval_reported_as_nominal_under_shared_host_load:true,
        same_binary_gfni_over_avx2_is_a_separate_campaign:true,
        isolation_claim:$audit[0].isolation_claim,
        windowed_observation:$audit[0].windowed_observation,
        outer_disclosure:$policy[0].outer_disclosure,
        shared_host_exposure:$policy[0].shared_host_exposure,
        ordinary_one_item_batch:.analysis["gf16-high-full"].encode,
        one_shot_encode:.analysis["gf16-high-full"].one_shot_encode
    }' "$campaign_dir/manifest.json" > "$lane/result-summary.json"
else
    /usr/bin/jq '{
        schema:"leopard2-v17-gfni-main-result-summary/v1",
        ratio_orientation:
            "exact_leopard1_native_time_over_leopard2_candidate_time",
        ratios_are_separate_correlated_and_must_not_be_multiplied:true,
        same_binary_gfni_over_avx2_is_a_separate_campaign:true,
        ordinary_one_item_batch:.analysis["gf16-high-full"].encode,
        one_shot_encode:.analysis["gf16-high-full"].one_shot_encode
    }' "$campaign_dir/manifest.json" > "$lane/result-summary.json"
fi

campaign_manifest_hash="$(/usr/bin/sha256sum \
    "$campaign_dir/manifest.json" | /usr/bin/cut -d' ' -f1)"
campaign_raw_hash="$(/usr/bin/sha256sum "$campaign_dir/raw.json" | \
    /usr/bin/cut -d' ' -f1)"
affinity_report_hash=null
affinity_binding_hash=null
if [[ "$passive_mode" == false ]]; then
    affinity_report_hash="$(/usr/bin/sha256sum "$affinity_report" | \
        /usr/bin/cut -d' ' -f1)"
    affinity_binding_hash="$(/usr/bin/sha256sum "$affinity_binding" | \
        /usr/bin/cut -d' ' -f1)"
fi
audit_hash="$(/usr/bin/sha256sum "$lane/audit.json" | \
    /usr/bin/cut -d' ' -f1)"

next_stage final_core_manifest
if [[ "$passive_mode" == true ]]; then
    controller_affinity_hash="$(/usr/bin/sha256sum \
        "$lane/controller-affinity.json" | /usr/bin/cut -d' ' -f1)"
    census_pre_hash="$(/usr/bin/sha256sum \
        "$lane/environment-census-pre.json" | /usr/bin/cut -d' ' -f1)"
    census_post_hash="$(/usr/bin/sha256sum \
        "$lane/environment-census-post.json" | /usr/bin/cut -d' ' -f1)"
    passive_policy_hash="$(/usr/bin/sha256sum \
        "$lane/passive-environment-policy.json" | /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n --slurpfile audit "$lane/audit.json" \
        --argjson campaign_exit_status "$campaign_status" \
        --argjson performance_gate_passed "$performance_gate_passed" \
        --argjson attempt "$attempt" \
        --argjson attempt_budget "$attempt_budget" \
        --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg main_commit "$main_commit" \
        --arg candidate_binary_sha256 "$candidate_binary_hash" \
        --arg candidate_archive_sha256 "$candidate_archive_hash" \
        --arg baseline_binary_sha256 "$baseline_binary_hash" \
        --arg baseline_archive_sha256 "$baseline_archive_hash" \
        --arg runner_sha256 "$runner_hash" \
        --arg auditor_sha256 "$auditor_hash" \
        --arg supervisor_sha256 "$supervisor_hash" \
        --arg passive_census_sha256 "$passive_census_hash" \
        --arg wrapper_sha256 "$wrapper_hash" \
        --arg build_order_normalizer_sha256 "$build_order_normalizer_hash" \
        --arg candidate_timing_normalization_report_sha256 \
            "$candidate_timing_normalization_report_hash" \
        --argjson candidate_timing_raw_tree_byte_identical \
            "$candidate_timing_raw_tree_byte_identical" \
        --arg candidate_test_normalization_report_sha256 \
            "$candidate_test_normalization_report_hash" \
        --argjson candidate_test_raw_tree_byte_identical \
            "$candidate_test_raw_tree_byte_identical" \
        --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
        --arg campaign_raw_sha256 "$campaign_raw_hash" \
        --arg audit_sha256 "$audit_hash" \
        --arg controller_affinity_sha256 "$controller_affinity_hash" \
        --arg environment_census_pre_sha256 "$census_pre_hash" \
        --arg environment_census_post_sha256 "$census_post_hash" \
        --arg passive_environment_policy_sha256 "$passive_policy_hash" \
        --arg sse2neon_commit "$candidate_submodule_commit" \
        --arg sse2neon_source_archive_sha256 "$sse2neon_source_archive_hash" \
        '{schema:"leopard2-v18-gfni-main-passive-core-manifest/v1",status:"complete",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,campaign_exit_status:$campaign_exit_status,evidence_valid:true,evidence_class:"passive-windowed-shared-host-observation/v1",performance_gate_passed:$performance_gate_passed,promotion_eligible:false,promotion_passed:false,causal_performance_claim_eligible:false,cpu_pair_exclusive:false,windowed_benchmark_cpu_nonidle_excess_jiffies:$audit[0].windowed_observation.windowed.benchmark_cpu_nonidle_excess_jiffies,windowed_reserved_sibling_nonidle_jiffies:$audit[0].windowed_observation.windowed.reserved_sibling_nonidle_jiffies,out_of_window_benchmark_cpu_nonidle_jiffies:$audit[0].windowed_observation.out_of_window.benchmark_cpu_nonidle_jiffies,out_of_window_reserved_sibling_nonidle_jiffies:$audit[0].windowed_observation.out_of_window.reserved_sibling_nonidle_jiffies,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,sse2neon_commit:$sse2neon_commit,sse2neon_source_archive_sha256:$sse2neon_source_archive_sha256,candidate_binary_sha256_pre_post:$candidate_binary_sha256,candidate_archive_sha256_pre_post:$candidate_archive_sha256,baseline_binary_sha256_pre_post:$baseline_binary_sha256,baseline_archive_sha256_pre_post:$baseline_archive_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,supervisor_sha256:$supervisor_sha256,supervisor_role:"retained-active-v1-verifier-only",active_affinity_supervisor_executed:false,passive_census_sha256:$passive_census_sha256,wrapper_sha256:$wrapper_sha256,build_order_normalizer_sha256:$build_order_normalizer_sha256,candidate_timing_normalization_report_sha256:$candidate_timing_normalization_report_sha256,candidate_test_normalization_report_sha256:$candidate_test_normalization_report_sha256,campaign_manifest_sha256:$campaign_manifest_sha256,campaign_raw_sha256:$campaign_raw_sha256,audit_sha256:$audit_sha256,controller_affinity_sha256:$controller_affinity_sha256,environment_census_pre_sha256:$environment_census_pre_sha256,environment_census_post_sha256:$environment_census_post_sha256,passive_environment_policy_sha256:$passive_environment_policy_sha256,isolation_claim:$audit[0].isolation_claim,build_reproduction:{candidate_timing_qualified_semantic_equal:true,candidate_timing_raw_tree_file_bytes_identical:$candidate_timing_raw_tree_byte_identical,candidate_tests_qualified_semantic_equal:true,candidate_tests_raw_tree_file_bytes_identical:$candidate_test_raw_tree_byte_identical,baseline_raw_byte_identical:true},producer_verification_passed:true,producer_verification_mode:"manifest-without-affinity-binding",independent_preseal_audit_passed:true,independent_auditor_supervision_mode:"windowed",independent_auditor_scope:"campaign semantics only; build-order qualification is separately recomputed by the frozen wrapper normalizer",canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116,ratio_policy:{ordinary_and_one_shot_are_separate:true,ratios_are_separate_correlated_and_must_not_be_multiplied:true,combined_or_stacked_ratio_emitted:false,same_binary_ratio_is_another_campaign:true,clustered_t_interval_reported_as_nominal_under_shared_host_load:true},postseal_policy:"passive-v2 shared-host observations always publish NOT_PROMOTED.json, independent of the observed speed gate; at most three preregistered attempts and every outcome is retained"}' \
        > "$lane/manifest.json"
else
/usr/bin/jq -n \
    --argjson campaign_exit_status "$campaign_status" \
    --argjson performance_gate_passed "$performance_gate_passed" \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg main_commit "$main_commit" \
    --arg candidate_binary_sha256 "$candidate_binary_hash" \
    --arg candidate_archive_sha256 "$candidate_archive_hash" \
    --arg baseline_binary_sha256 "$baseline_binary_hash" \
    --arg baseline_archive_sha256 "$baseline_archive_hash" \
    --arg runner_sha256 "$runner_hash" \
    --arg auditor_sha256 "$auditor_hash" \
    --arg supervisor_sha256 "$supervisor_hash" \
    --arg wrapper_sha256 "$wrapper_hash" \
    --arg build_order_normalizer_sha256 "$build_order_normalizer_hash" \
    --arg candidate_timing_normalization_report_sha256 \
        "$candidate_timing_normalization_report_hash" \
    --argjson candidate_timing_raw_tree_byte_identical \
        "$candidate_timing_raw_tree_byte_identical" \
    --arg candidate_test_normalization_report_sha256 \
        "$candidate_test_normalization_report_hash" \
    --argjson candidate_test_raw_tree_byte_identical \
        "$candidate_test_raw_tree_byte_identical" \
    --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
    --arg campaign_raw_sha256 "$campaign_raw_hash" \
    --arg affinity_report_sha256 "$affinity_report_hash" \
    --arg affinity_binding_sha256 "$affinity_binding_hash" \
    --arg audit_sha256 "$audit_hash" \
    --arg sse2neon_commit "$candidate_submodule_commit" \
    --arg sse2neon_source_archive_sha256 "$sse2neon_source_archive_hash" \
    '{schema:"leopard2-v17-gfni-main-core-manifest/v3",status:"complete",campaign_exit_status:$campaign_exit_status,evidence_valid:true,performance_gate_passed:$performance_gate_passed,promotion_passed:false,promotion_requires_completion_envelope:true,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,sse2neon_commit:$sse2neon_commit,sse2neon_source_archive_sha256:$sse2neon_source_archive_sha256,candidate_binary_sha256_pre_post:$candidate_binary_sha256,candidate_archive_sha256_pre_post:$candidate_archive_sha256,baseline_binary_sha256_pre_post:$baseline_binary_sha256,baseline_archive_sha256_pre_post:$baseline_archive_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,supervisor_sha256:$supervisor_sha256,wrapper_sha256:$wrapper_sha256,build_order_normalizer_sha256:$build_order_normalizer_sha256,candidate_timing_normalization_report_sha256:$candidate_timing_normalization_report_sha256,candidate_test_normalization_report_sha256:$candidate_test_normalization_report_sha256,campaign_manifest_sha256:$campaign_manifest_sha256,campaign_raw_sha256:$campaign_raw_sha256,affinity_report_sha256:$affinity_report_sha256,affinity_binding_sha256:$affinity_binding_sha256,audit_sha256:$audit_sha256,build_reproduction:{candidate_timing_qualified_semantic_equal:true,candidate_timing_raw_tree_file_bytes_identical:$candidate_timing_raw_tree_byte_identical,candidate_tests_qualified_semantic_equal:true,candidate_tests_raw_tree_file_bytes_identical:$candidate_test_raw_tree_byte_identical,baseline_raw_byte_identical:true},producer_verification_passed:true,independent_preseal_audit_passed:true,independent_auditor_scope:"campaign semantics only; build-order qualification is separately recomputed by the frozen wrapper normalizer",canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116,ratio_policy:{ordinary_and_one_shot_are_separate:true,ratios_are_separate_correlated_and_must_not_be_multiplied:true,combined_or_stacked_ratio_emitted:false,same_binary_ratio_is_another_campaign:true},postseal_policy:"promotion requires the enclosing COMPLETED.json written only after byte-identical independent pre/post audits, qualified-semantic build closure verification, core SHA verification, clean-source recheck, and a zero campaign exit"}' \
    > "$lane/manifest.json"
fi

next_stage seal_core
(
    cd "$lane"
    /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
        /usr/bin/sort -z | \
        /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
)
trap postseal_failure_record ERR
/usr/bin/find "$lane" -type f -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$lane" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
verify_sealed_tree "$lane"
(
    cd "$lane"
    /usr/bin/sha256sum -c SHA256SUMS
) > "$envelope/core-sha256-check.txt"

stage=independent_postseal_audit
/usr/bin/printf 'AUTHORITATIVE_STAGE independent_postseal_audit\n'
postseal_audit_args=(
    /usr/bin/python3 -I -S -B
    "$lane/audit_v17_gfni_main_compare.py"
    --manifest "$campaign_dir/manifest.json"
)
if [[ "$passive_mode" == true ]]; then
    postseal_audit_args+=(--supervision windowed)
fi
"${postseal_audit_args[@]}" \
    --output "$envelope/postseal-audit.json" \
    > "$envelope/postseal-audit-summary.json" \
    2> "$envelope/postseal-audit-stderr.log"
/usr/bin/cmp "$lane/audit.json" "$envelope/postseal-audit.json"
(
    cd "$lane"
    /usr/bin/sha256sum -c SHA256SUMS
) > "$envelope/postaudit-core-sha256-check.txt"
verify_sealed_tree "$lane"

core_manifest_sha="$(/usr/bin/sha256sum "$lane/manifest.json" | \
    /usr/bin/cut -d' ' -f1)"
core_sha256sums_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
    /usr/bin/cut -d' ' -f1)"
postseal_audit_sha="$(/usr/bin/sha256sum \
    "$envelope/postseal-audit.json" | /usr/bin/cut -d' ' -f1)"

if [[ "$passive_mode" == true || "$performance_gate_passed" != true ]]; then
    stage=publish_not_promoted_envelope
    /usr/bin/printf 'AUTHORITATIVE_STAGE publish_not_promoted_envelope\n'
    if [[ "$passive_mode" == true ]]; then
        /usr/bin/jq -n --slurpfile audit "$lane/audit.json" \
            --argjson campaign_exit_status "$campaign_status" \
            --argjson performance_gate_passed "$performance_gate_passed" \
            --argjson attempt "$attempt" \
            --argjson attempt_budget "$attempt_budget" \
            --arg attempt_lineage_sha256 "$attempt_lineage_sha256" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg main_commit "$main_commit" \
            --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
            --arg audit_sha256 "$audit_hash" \
            --arg postseal_audit_sha256 "$postseal_audit_sha" \
            --arg core_manifest_sha256 "$core_manifest_sha" \
            --arg core_sha256sums_sha256 "$core_sha256sums_sha" \
            '{schema:"leopard2-v18-gfni-main-passive-not-promoted-envelope/v1",status:"complete",acquisition_generation:"passive-v2",attempt:$attempt,attempt_budget:$attempt_budget,attempt_lineage_sha256:$attempt_lineage_sha256,evidence_class:"passive-windowed-shared-host-observation/v1",promotion_eligible:false,promotion_passed:false,causal_performance_claim_eligible:false,cpu_pair_exclusive:false,windowed_benchmark_cpu_nonidle_excess_jiffies:$audit[0].windowed_observation.windowed.benchmark_cpu_nonidle_excess_jiffies,windowed_reserved_sibling_nonidle_jiffies:$audit[0].windowed_observation.windowed.reserved_sibling_nonidle_jiffies,out_of_window_benchmark_cpu_nonidle_jiffies:$audit[0].windowed_observation.out_of_window.benchmark_cpu_nonidle_jiffies,out_of_window_reserved_sibling_nonidle_jiffies:$audit[0].windowed_observation.out_of_window.reserved_sibling_nonidle_jiffies,performance_gate_passed:$performance_gate_passed,evidence_valid:true,campaign_exit_status:$campaign_exit_status,preseal_audit_passed:true,postseal_audit_passed:true,postseal_audit_byte_identical:true,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,campaign_manifest_sha256:$campaign_manifest_sha256,audit_sha256:$audit_sha256,postseal_audit_sha256:$postseal_audit_sha256,core_manifest_sha256:$core_manifest_sha256,core_sha256sums_sha256:$core_sha256sums_sha256,isolation_claim:$audit[0].isolation_claim}' \
            > "$envelope/NOT_PROMOTED.json"
    else
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg main_commit "$main_commit" \
        --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
        --arg audit_sha256 "$audit_hash" \
        --arg postseal_audit_sha256 "$postseal_audit_sha" \
        --arg core_manifest_sha256 "$core_manifest_sha" \
        --arg core_sha256sums_sha256 "$core_sha256sums_sha" \
        '{schema:"leopard2-v17-gfni-main-not-promoted-envelope/v1",status:"complete",promotion_passed:false,performance_gate_passed:false,evidence_valid:true,campaign_exit_status:$campaign_exit_status,preseal_audit_passed:true,postseal_audit_passed:true,postseal_audit_byte_identical:true,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,campaign_manifest_sha256:$campaign_manifest_sha256,audit_sha256:$audit_sha256,postseal_audit_sha256:$postseal_audit_sha256,core_manifest_sha256:$core_manifest_sha256,core_sha256sums_sha256:$core_sha256sums_sha256}' \
        > "$envelope/NOT_PROMOTED.json"
    fi
    /usr/bin/printf '' > "$envelope/SHA256SUMS"
    write_tree_metadata_manifest "$envelope"
    (
        cd "$envelope"
        /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
            /usr/bin/sort -z | \
            /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    )
    /usr/bin/find "$envelope" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$envelope" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    verification_root="$(/usr/bin/mktemp -d \
        /tmp/leopard-v17-gfni-main-not-promoted-verify.XXXXXX)"
    "$lane/run-authoritative.sh" --verify "$envelope" \
        > "$verification_root/verification.txt" \
        2> "$verification_root/verification-stderr.log"
    /usr/bin/printf \
        'sealed_not_promoted_envelope=%s\nexternal_verification=%s\n' \
        "$envelope" "$verification_root"
    trap - ERR
    exit 0
fi

stage=publish_completion_envelope
/usr/bin/printf 'AUTHORITATIVE_STAGE publish_completion_envelope\n'
test "$campaign_status" -eq 0
completion_json="$(/usr/bin/jq -c -n \
    --argjson campaign_exit_status "$campaign_status" \
    --arg commit "$commit" \
    --arg tree "$tree" \
    --arg main_commit "$main_commit" \
    --arg candidate_binary_sha256 "$candidate_binary_hash" \
    --arg baseline_binary_sha256 "$baseline_binary_hash" \
    --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
    --arg audit_sha256 "$audit_hash" \
    --arg postseal_audit_sha256 "$postseal_audit_sha" \
    --arg core_manifest_sha256 "$core_manifest_sha" \
    --arg core_sha256sums_sha256 "$core_sha256sums_sha" \
    '{schema:"leopard2-v17-gfni-main-completion-envelope/v1",status:"complete",promotion_passed:true,performance_gate_passed:true,evidence_valid:true,campaign_exit_status:$campaign_exit_status,preseal_audit_passed:true,postseal_audit_passed:true,postseal_audit_byte_identical:true,core_sha256sums_verified:true,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,candidate_binary_sha256:$candidate_binary_sha256,baseline_binary_sha256:$baseline_binary_sha256,campaign_manifest_sha256:$campaign_manifest_sha256,audit_sha256:$audit_sha256,postseal_audit_sha256:$postseal_audit_sha256,core_manifest_sha256:$core_manifest_sha256,core_sha256sums_sha256:$core_sha256sums_sha256,ratios_are_separate_correlated_and_must_not_be_multiplied:true,verification_command:["core/run-authoritative.sh","--verify","<absolute-envelope-path>"]}')"
completion_hash="$(/usr/bin/printf '%s\n' "$completion_json" | \
    /usr/bin/sha256sum | /usr/bin/cut -d' ' -f1)"
exec 8> "$envelope/COMPLETED.json"
/usr/bin/printf '' > "$envelope/SHA256SUMS"
write_tree_metadata_manifest "$envelope"
(
    cd "$envelope"
    /usr/bin/find . -type f ! -path './SHA256SUMS' \
        ! -path './COMPLETED.json' -print0 | \
        /usr/bin/sort -z | \
        /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    /usr/bin/printf '%s  ./COMPLETED.json\n' "$completion_hash" \
        >> SHA256SUMS
)
/usr/bin/find "$envelope" -type f -perm /222 \
    ! -path "$envelope/COMPLETED.json" -exec /usr/bin/chmod a-w {} +
/usr/bin/find "$envelope" -type d -perm /222 \
    -exec /usr/bin/chmod a-w {} +
/usr/bin/python3 -I -S -B -c '
import os, sys
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
    os.close(directory)
' "$completion_json" "$envelope"
exec 8>&-
/usr/bin/chmod a-w "$envelope/COMPLETED.json"

/usr/bin/printf 'AUTHORITATIVE_STAGE verify_completion_envelope\n'
verification_root="$(/usr/bin/mktemp -d \
    /tmp/leopard-v17-gfni-main-envelope-verify.XXXXXX)"
"$lane/run-authoritative.sh" --verify "$envelope" \
    > "$verification_root/verification.txt" \
    2> "$verification_root/verification-stderr.log"

/usr/bin/printf 'AUTHORITATIVE_STAGE complete\n'
/usr/bin/jq '{
    ratio_orientation,
    ratios_are_separate_correlated_and_must_not_be_multiplied,
    ordinary_one_item_batch,
    one_shot_encode
}' "$lane/result-summary.json"
/usr/bin/printf 'campaign_manifest_sha256=%s\n' "$campaign_manifest_hash"
/usr/bin/printf 'sealed_envelope=%s\nexternal_verification=%s\n' \
    "$envelope" "$verification_root"
trap - ERR
exit 0
