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
relative_supervisor=tools/leopard2_affinity_supervisor.py
lock=/tmp/leopard-gf8-authoritative.lock

tools_to_hash=(
    /usr/bin/bash /usr/bin/env /usr/bin/c++ /usr/bin/cc /usr/bin/ld
    /usr/bin/ar /usr/bin/ranlib /usr/bin/cmake /usr/bin/gmake /usr/bin/git
    /usr/bin/objdump /usr/bin/readelf /usr/bin/nm /usr/bin/python3
    /usr/bin/taskset /usr/bin/ldd /usr/bin/flock /usr/bin/jq
    /usr/bin/sha256sum /usr/bin/find /usr/bin/grep /usr/bin/cmp /usr/bin/cp
    /usr/bin/stat /usr/bin/chmod /usr/bin/sort /usr/bin/xargs /usr/bin/cut
    /usr/bin/tee /usr/bin/ctest /usr/bin/cat /usr/bin/dirname /usr/bin/lscpu
    /usr/bin/mkdir /usr/bin/mktemp /usr/bin/printf /usr/bin/readlink
    /usr/bin/uname /usr/bin/mv /usr/bin/diff
)

require_empty_output()
{
    local observed_output
    observed_output="$("$@")"
    test -z "$observed_output"
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
    local verified_core="$verified_envelope/core"
    local status_file=
    local status_schema=
    local status_value=
    local promotion_value=
    local campaign_exit_value=
    local campaign_manifest_sha=
    local replay_campaign=
    local replay_controller_root=
    local verifier_tmp=

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
        cd "$verified_core"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null
    (
        cd "$verified_envelope"
        /usr/bin/sha256sum -c SHA256SUMS
    ) >/dev/null
    verify_sealed_tree "$verified_envelope"
    verify_tree_metadata_manifest "$verified_envelope"

    if [[ -f "$verified_envelope/COMPLETED.json" ]]; then
        test ! -e "$verified_envelope/NOT_PROMOTED.json"
        test ! -e "$verified_envelope/FAILED.json"
        status_file="$verified_envelope/COMPLETED.json"
        status_schema=leopard2-v17-gfni-main-completion-envelope/v1
        status_value=complete
        promotion_value=true
    elif [[ -f "$verified_envelope/NOT_PROMOTED.json" ]]; then
        test ! -e "$verified_envelope/COMPLETED.json"
        test ! -e "$verified_envelope/FAILED.json"
        status_file="$verified_envelope/NOT_PROMOTED.json"
        status_schema=leopard2-v17-gfni-main-not-promoted-envelope/v1
        status_value=complete
        promotion_value=false
    else
        test -f "$verified_envelope/FAILED.json"
        test ! -e "$verified_envelope/COMPLETED.json"
        test ! -e "$verified_envelope/NOT_PROMOTED.json"
        status_file="$verified_envelope/FAILED.json"
        status_schema=leopard2-v17-gfni-main-failed-envelope/v1
        status_value=failed
        promotion_value=false
    fi
    test "$(/usr/bin/jq -er '.schema | strings' "$status_file")" = \
        "$status_schema"
    test "$(/usr/bin/jq -er '.status | strings' "$status_file")" = \
        "$status_value"
    test "$(/usr/bin/jq -er \
        '.promotion_passed | booleans | tostring' "$status_file")" = \
        "$promotion_value"
    campaign_exit_value="$(/usr/bin/jq -er \
        '.campaign_exit_status | numbers' "$status_file")"
    test "$(/usr/bin/sha256sum "$verified_core/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)" = \
        "$(/usr/bin/jq -er '.core_sha256sums_sha256 | strings' \
            "$status_file")"

    if [[ "$status_value" == complete ]]; then
        test "$campaign_exit_value" -eq 0
        test -f "$verified_core/manifest.json"
        test -f "$verified_core/campaign/manifest.json"
        test -f "$verified_core/audit.json"
        test -f "$verified_core/affinity-report.json"
        test -f "$verified_core/affinity-binding.json"
        test -f "$verified_core/run_abba.py"
        test -f "$verified_core/git_capture.py"
        test -f "$verified_core/balanced_evidence_common.py"
        test -f "$verified_core/leopard2_build_provenance.py"
        test -f "$verified_core/leopard2_affinity_supervisor.py"
        test -f "$verified_core/sse2neon-source.tar"
        test -f "$verified_core/build-closure.json"
        test -f "$verified_core/candidate-test-temporal-closure.json"
        test -f \
            "$verified_core/build-closure/candidate-test-selected-SHA256SUMS"
        test -f \
            "$verified_core/build-closure/live-candidate-tests/leopard2_backend_failures_test"
        test -f \
            "$verified_core/build-closure/replay-candidate-tests/leopard2_backend_failures_test"
        test -f \
            "$verified_core/build-closure/live-candidate-tests/bench_leopard2_prevalidated_batch"
        test -f \
            "$verified_core/build-closure/replay-candidate-tests/bench_leopard2_prevalidated_batch"
        test -f "$verified_envelope/postseal-audit.json"
        test -f "$verified_core/audit_v17_gfni_main_compare.py"
        /usr/bin/cmp "$verified_core/audit.json" \
            "$verified_envelope/postseal-audit.json"
        /usr/bin/cmp "$verified_core/run-authoritative.sh" \
            "$verified_core/build-closure/committed-wrapper.sh"
        /usr/bin/cmp "$verified_core/audit_v17_gfni_main_compare.py" \
            "$verified_core/build-closure/committed-auditor.py"
        /usr/bin/cmp "$verified_core/leopard2_build_provenance.py" \
            "$verified_core/build-closure/committed-build-provenance.py"
        /usr/bin/cmp \
            "$verified_core/build-closure/live-candidate-tests/leopard2_backend_failures_test" \
            "$verified_core/build-closure/replay-candidate-tests/leopard2_backend_failures_test"
        /usr/bin/cmp \
            "$verified_core/build-closure/live-candidate-tests/bench_leopard2_prevalidated_batch" \
            "$verified_core/build-closure/replay-candidate-tests/bench_leopard2_prevalidated_batch"
        /usr/bin/diff -qr \
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
            '.candidate_tests.two_clean_builds | booleans | tostring' \
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
            .candidate_tests.posttest_byte_identical == true and
            .candidate_tests.postcampaign_revalidation_required == true
        ' "$verified_core/build-closure.json" >/dev/null
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/candidate-test-temporal-closure.json")" = \
            leopard2-v17-gfni-main-candidate-test-temporal-closure/v1
        test "$(/usr/bin/jq -er \
            '.selected_sha256sums_sha256 | strings' \
            "$verified_core/candidate-test-temporal-closure.json")" = \
            "$(/usr/bin/jq -er \
            '.candidate_tests.selected_sha256sums_sha256 | strings' \
            "$verified_core/build-closure.json")"
        /usr/bin/jq -e '
            .two_clean_builds_byte_identical == true and
            .posttest_byte_identical == true and
            .postcampaign_byte_identical == true and
            .canonical_test_build_frozen_during_campaign == true
        ' "$verified_core/candidate-test-temporal-closure.json" >/dev/null
        test "$(/usr/bin/jq -er '.evidence_valid | booleans | tostring' \
            "$status_file")" = true
        test "$(/usr/bin/jq -er \
            '.performance_gate_passed | booleans | tostring' \
            "$status_file")" = "$promotion_value"
        test "$(/usr/bin/jq -er '.schema | strings' \
            "$verified_core/manifest.json")" = \
            leopard2-v17-gfni-main-core-manifest/v1
        test "$(/usr/bin/jq -er '.campaign_exit_status | numbers' \
            "$verified_core/manifest.json")" -eq 0
        test "$(/usr/bin/jq -er '.evidence_valid | booleans | tostring' \
            "$verified_core/manifest.json")" = true
        test "$(/usr/bin/jq -er \
            '.promotion_requires_completion_envelope | booleans | tostring' \
            "$verified_core/manifest.json")" = true
        test "$(/usr/bin/jq -er '.promotion_passed | booleans | tostring' \
            "$verified_core/manifest.json")" = false
        test "$(/usr/bin/jq -er \
            '.performance_gate_passed | booleans | tostring' \
            "$verified_core/manifest.json")" = "$promotion_value"
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
        test "$(/usr/bin/sha256sum \
            "$verified_envelope/postseal-audit.json" | \
            /usr/bin/cut -d' ' -f1)" = \
            "$(/usr/bin/jq -er '.postseal_audit_sha256 | strings' \
                "$status_file")"
        verifier_tmp="$(/usr/bin/mktemp -d \
            /tmp/leopard-v17-gfni-main-verify.XXXXXX)"
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
        /usr/bin/python3 -I -S -B \
            "$verified_core/audit_v17_gfni_main_compare.py" \
            --manifest "$verified_core/campaign/manifest.json" \
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
    fi
    /usr/bin/printf 'authoritative v17 envelope verified: %s\n' \
        "$verified_envelope"
}

if [[ $# -eq 2 && $1 == --verify && $2 == /* ]]; then
    verify_envelope "$2"
    exit 0
fi

if [[ $# -ne 1 || $1 != /* ]]; then
    /usr/bin/printf 'usage: %s /absolute/repository/.research/envelope\n' \
        "$0" >&2
    /usr/bin/printf '       %s --verify /absolute/repository/.research/envelope\n' \
        "$0" >&2
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

next_stage()
{
    stage=$1
    /usr/bin/printf 'AUTHORITATIVE_STAGE %s\n' "$stage" | \
        /usr/bin/tee -a "$lane/wrapper-stage.log"
}

failure_record()
{
    local status=$?
    trap - ERR
    set +e
    if [[ -n "$lane" && -d "$lane" && -w "$lane" ]]; then
        /usr/bin/jq -n \
            --arg stage "$stage" \
            --argjson exit_status "$status" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            '{schema:"leopard2-v17-gfni-main-wrapper-failure/v1",status:"failed",promotion_passed:false,stage:$stage,exit_status:$exit_status,source_commit:$commit,source_tree:$tree}' \
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
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg core_sha256sums_sha256 "$core_sha" \
            '{schema:"leopard2-v17-gfni-main-failed-envelope/v1",status:"failed",promotion_passed:false,campaign_exit_status:$exit_status,stage:$stage,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256}' \
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
    trap - ERR
    set +e
    exec 8>&-
    if [[ -d "$envelope" && -w "$envelope" && \
          -f "$lane/SHA256SUMS" ]]; then
        if (
            cd "$lane" || exit
            /usr/bin/sha256sum -c SHA256SUMS
        ) >/dev/null 2>&1; then
            failed_core_verified=true
        fi
        failed_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
            /usr/bin/cut -d' ' -f1)"
        for displaced_status in \
            COMPLETED.json NOT_PROMOTED.json TREE-METADATA.json; do
            if [[ -e "$envelope/$displaced_status" ]]; then
                /usr/bin/mv "$envelope/$displaced_status" \
                    "$envelope/UNCOMMITTED-$displaced_status"
            fi
        done
        /usr/bin/jq -n \
            --argjson exit_status "$status" \
            --argjson core_sha256sums_verified "$failed_core_verified" \
            --arg stage "$stage" \
            --arg commit "$commit" \
            --arg tree "$tree" \
            --arg baseline_commit "$main_commit" \
            --arg core_sha256sums_sha256 "$failed_core_sha" \
            '{schema:"leopard2-v17-gfni-main-failed-envelope/v1",status:"failed",promotion_passed:false,campaign_exit_status:$exit_status,stage:$stage,source_commit:$commit,source_tree:$tree,baseline_commit:$baseline_commit,core_sha256sums_verified:$core_sha256sums_verified,core_sha256sums_sha256:$core_sha256sums_sha256}' \
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
    fi
    /usr/bin/printf \
        'authoritative wrapper failed after core sealing in stage %s with status %s\n' \
        "$stage" "$status" >&2
    exit "$status"
}

if [[ -e "$envelope" ]]; then
    /usr/bin/printf 'refusing to reuse artifact envelope: %s\n' \
        "$envelope" >&2
    exit 2
fi
/usr/bin/mkdir -p "$(/usr/bin/dirname "$envelope")"
test "$(/usr/bin/readlink -f "$(/usr/bin/dirname "$envelope")")" = \
    "$(/usr/bin/dirname "$envelope")"
/usr/bin/mkdir -m 0700 "$envelope"
test "$(/usr/bin/readlink -f "$envelope")" = "$envelope"
lane="$envelope/core"
/usr/bin/mkdir -m 0700 "$lane"
/usr/bin/cp --reflink=never "$0" "$lane/run-authoritative.sh"
/usr/bin/chmod 0555 "$lane/run-authoritative.sh"
test ! -L "$lane/run-authoritative.sh"
test "$(/usr/bin/stat -c %h "$lane/run-authoritative.sh")" = 1
trap failure_record ERR

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
assert 52 in allowed and 116 in allowed and len(set(allowed) - {52, 116}) > 0
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
/usr/bin/git -C "$repo" archive --format=tar \
    --prefix=candidate-source/ --output="$lane/candidate-source.tar" "$commit"
/usr/bin/git -C "$candidate_source" archive --format=tar \
    --prefix=candidate-source/ \
    --output="$work_root/replayed-candidate-source.tar" "$commit"
/usr/bin/cmp "$lane/candidate-source.tar" \
    "$work_root/replayed-candidate-source.tar"
/usr/bin/git -C "$repo" archive --format=tar \
    --prefix=leopard1-source/ --output="$lane/leopard1-source.tar" \
    "$main_commit"
/usr/bin/git -C "$baseline_source" archive --format=tar \
    --prefix=leopard1-source/ \
    --output="$work_root/replayed-leopard1-source.tar" "$main_commit"
/usr/bin/cmp "$lane/leopard1-source.tar" \
    "$work_root/replayed-leopard1-source.tar"
/usr/bin/git -C "$repo/sse2neon" archive --format=tar \
    --prefix=sse2neon-source/ --output="$lane/sse2neon-source.tar" \
    "$candidate_submodule_commit"
/usr/bin/git -C "$candidate_source/sse2neon" archive --format=tar \
    --prefix=sse2neon-source/ \
    --output="$work_root/replayed-candidate-sse2neon-source.tar" \
    "$candidate_submodule_commit"
/usr/bin/git -C "$baseline_source/sse2neon" archive --format=tar \
    --prefix=sse2neon-source/ \
    --output="$work_root/replayed-baseline-sse2neon-source.tar" \
    "$baseline_submodule_commit"
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
/usr/bin/chmod 0555 \
    "$lane/audit_v17_gfni_main_compare.py" \
    "$lane/run_abba.py" \
    "$lane/leopard2_affinity_supervisor.py"
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
    "$lane/leopard2_affinity_supervisor.py"; do
    test ! -L "$frozen_controller"
    test "$(/usr/bin/stat -c %h "$frozen_controller")" = 1
done

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

next_stage candidate_correctness_build_byte_closure
/usr/bin/diff -qr "$lane/build-closure/live-candidate-tests" \
    "$lane/build-closure/replay-candidate-tests" \
    > "$lane/candidate-test-build-byte-diff.txt"
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
/usr/bin/python3 -I -S -B \
    "$candidate_source/$relative_supervisor" self-test \
    > "$lane/supervisor-self-test.log" 2>&1

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

next_stage live_replay_byte_closure
/usr/bin/diff -qr "$lane/build-closure/live-candidate" \
    "$lane/build-closure/replay-candidate" \
    > "$lane/candidate-build-byte-diff.txt"
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
wrapper_hash="$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | \
    /usr/bin/cut -d' ' -f1)"
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
    --arg candidate_source_archive_sha256 "$candidate_source_archive_hash" \
    --arg baseline_source_archive_sha256 "$baseline_source_archive_hash" \
    --arg sse2neon_commit "$candidate_submodule_commit" \
    --arg sse2neon_source_archive_sha256 "$sse2neon_source_archive_hash" \
    '{schema:"leopard2-v17-gfni-main-build-closure/v1",candidate:{commit:$commit,tree:$tree,source:$candidate_source,build:$candidate_build,profile:"standard AUTO Release; tests off; benchmarks on; GF8/GF16 on",binary_sha256:$candidate_binary_sha256,archive_sha256:$candidate_archive_sha256,source_archive_sha256:$candidate_source_archive_sha256},candidate_tests:{build:$candidate_test_build,profile:"standard AUTO Release; tests and benchmarks on; GF8/GF16 on",selected_files:["bench_leopard2","bench_leopard2_prevalidated_batch","leopard2_auto_gf16_gfni_production_test","leopard2_backend_failures_test","leopard2_legacy_golden_test","libleopard.a","libleopard_test_hooks.a","generated/leopard2-benchmark-attestation/leopard2_benchmark_source_attestation.h","generated/leopard2-benchmark-attestation/leopard2_benchmark_build_configuration.txt"],selected_sha256sums_sha256:$candidate_test_selected_sha256sums_sha256,backend_failures_test_sha256:$backend_failures_test_sha256,prevalidated_batch_test_sha256:$prevalidated_batch_test_sha256,two_clean_builds:true,posttest_byte_identical:true,postcampaign_revalidation_required:true,complete_object_link_cache_closure:true},baseline:{commit:$main_commit,source:$baseline_source,build:$baseline_build,profile:"canonical Leopard1 native Release (-march=native; LEO_MAIN_PURE_AVX2=OFF)",binary_sha256:$baseline_binary_sha256,archive_sha256:$baseline_archive_sha256,source_archive_sha256:$baseline_source_archive_sha256},sse2neon:{commit:$sse2neon_commit,archive_prefix:"sse2neon-source/",source_archive_sha256:$sse2neon_source_archive_sha256,reproduced_from_candidate_and_baseline_clones:true},generator:"Unix Makefiles",compiler:"/usr/bin/c++",byte_identical:{candidate_two_clean_builds:true,candidate_test_two_clean_builds:true,baseline_two_clean_builds:true,identical_absolute_source_and_build_paths:true,selected_complete_object_link_cache_closure:true},controllers:{runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,supervisor_sha256:$supervisor_sha256,wrapper_sha256:$wrapper_sha256},canonical_lock:"/tmp/leopard-gf8-authoritative.lock",python_controller:["/usr/bin/python3","-I","-S","-B"]}' \
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
    "owner": "v17 exact-main authoritative wrapper",
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
' "$reservation"
test "$(/usr/bin/stat -c %a "$reservation")" = 600

next_stage supervised_campaign
campaign_dir="$lane/campaign"
affinity_report="$lane/affinity-report.json"
campaign_command=(
    /usr/bin/python3 -I -S -B
    "$candidate_source/$relative_supervisor"
    run
    --report "$affinity_report"
    --reserved-cpus 52,116
    --
    # Supervisor v14 requires the child Python source as direct argv[1].  It
    # seals both files and launches its bootstrap with isolated/no-site flags;
    # the clean environment's PYTHONDONTWRITEBYTECODE=1 supplies -B semantics.
    /usr/bin/python3
    "$candidate_source/$relative_runner"
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
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg main_commit "$main_commit" \
        --arg candidate_binary_sha256 "$candidate_binary_hash" \
        --arg baseline_binary_sha256 "$baseline_binary_hash" \
        --arg failure_sha256 "$failure_sha" \
        '{schema:"leopard2-v17-gfni-main-failed-core-manifest/v1",status:"failed",promotion_passed:false,campaign_exit_status:$campaign_exit_status,failure_verify_status:$failure_verify_status,failure_verified:$failure_verified,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,candidate_binary_sha256:$candidate_binary_sha256,baseline_binary_sha256:$baseline_binary_sha256,failure_sha256:($failure_sha256 | if . == "null" then null else . end),canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116}' \
        > "$lane/manifest.json"
    trap - ERR
    (
        cd "$lane"
        /usr/bin/find . -type f ! -name SHA256SUMS -print0 | \
            /usr/bin/sort -z | \
            /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    )
    /usr/bin/find "$lane" -type f -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    /usr/bin/find "$lane" -type d -perm /222 \
        -exec /usr/bin/chmod a-w {} +
    failed_core_sha="$(/usr/bin/sha256sum "$lane/SHA256SUMS" | \
        /usr/bin/cut -d' ' -f1)"
    /usr/bin/jq -n \
        --argjson campaign_exit_status "$campaign_status" \
        --argjson failure_verified "$failure_verified" \
        --arg commit "$commit" \
        --arg tree "$tree" \
        --arg core_sha256sums_sha256 "$failed_core_sha" \
        '{schema:"leopard2-v17-gfni-main-failed-envelope/v1",status:"failed",promotion_passed:false,campaign_exit_status:$campaign_exit_status,failure_verified:$failure_verified,source_commit:$commit,source_tree:$tree,core_sha256sums_sha256:$core_sha256sums_sha256}' \
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
    exit "$campaign_status"
fi

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

next_stage producer_bundle_verification
/usr/bin/python3 -I -S -B \
    "$candidate_source/$relative_runner" verify \
    --manifest "$campaign_dir/manifest.json" \
    --affinity-binding "$affinity_binding" \
    > "$lane/producer-verification.log" \
    2> "$lane/producer-verification-stderr.log"

next_stage independent_preseal_audit
/usr/bin/python3 -I -S -B \
    "$lane/audit_v17_gfni_main_compare.py" \
    --manifest "$campaign_dir/manifest.json" \
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
    '{schema:"leopard2-v17-gfni-main-candidate-test-temporal-closure/v1",selected_sha256sums_sha256:$selected_sha256sums_sha256,two_clean_builds_byte_identical:true,posttest_byte_identical:true,postcampaign_byte_identical:true,canonical_test_build_frozen_during_campaign:true}' \
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
test "$(/usr/bin/sha256sum "$lane/run-authoritative.sh" | \
    /usr/bin/cut -d' ' -f1)" = "$wrapper_hash"
test "$(/usr/bin/sha256sum "$lane/run_abba.py" | \
    /usr/bin/cut -d' ' -f1)" = "$runner_hash"
test "$(/usr/bin/sha256sum "$lane/audit_v17_gfni_main_compare.py" | \
    /usr/bin/cut -d' ' -f1)" = "$auditor_hash"
test "$(/usr/bin/sha256sum "$lane/leopard2_affinity_supervisor.py" | \
    /usr/bin/cut -d' ' -f1)" = "$supervisor_hash"
require_empty_output /usr/bin/git -C "$candidate_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
require_empty_output /usr/bin/git -C "$baseline_source" status \
    --porcelain=v1 --untracked-files=normal --ignore-submodules=none
test "$(/usr/bin/git -C "$candidate_source" rev-parse HEAD)" = "$commit"
test "$(/usr/bin/git -C "$candidate_source" rev-parse 'HEAD^{tree}')" = "$tree"
test "$(/usr/bin/git -C "$baseline_source" rev-parse HEAD)" = "$main_commit"
/usr/bin/git -C "$candidate_source" archive --format=tar \
    --prefix=candidate-source/ \
    --output="$work_root/post-candidate-source.tar" "$commit"
/usr/bin/git -C "$baseline_source" archive --format=tar \
    --prefix=leopard1-source/ \
    --output="$work_root/post-leopard1-source.tar" "$main_commit"
/usr/bin/git -C "$candidate_source/sse2neon" archive --format=tar \
    --prefix=sse2neon-source/ \
    --output="$work_root/post-candidate-sse2neon-source.tar" \
    "$candidate_submodule_commit"
/usr/bin/git -C "$baseline_source/sse2neon" archive --format=tar \
    --prefix=sse2neon-source/ \
    --output="$work_root/post-baseline-sse2neon-source.tar" \
    "$baseline_submodule_commit"
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
    "$campaign_dir/manifest.json")" = leopard2-main-compare-manifest/v17
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
/usr/bin/jq '{
    schema:"leopard2-v17-gfni-main-result-summary/v1",
    ratio_orientation:
        "exact_leopard1_native_time_over_leopard2_candidate_time",
    ratios_are_separate_correlated_and_must_not_be_multiplied:true,
    same_binary_gfni_over_avx2_is_a_separate_campaign:true,
    ordinary_one_item_batch:.analysis["gf16-high-full"].encode,
    one_shot_encode:.analysis["gf16-high-full"].one_shot_encode
}' "$campaign_dir/manifest.json" > "$lane/result-summary.json"

campaign_manifest_hash="$(/usr/bin/sha256sum \
    "$campaign_dir/manifest.json" | /usr/bin/cut -d' ' -f1)"
campaign_raw_hash="$(/usr/bin/sha256sum "$campaign_dir/raw.json" | \
    /usr/bin/cut -d' ' -f1)"
affinity_report_hash="$(/usr/bin/sha256sum "$affinity_report" | \
    /usr/bin/cut -d' ' -f1)"
affinity_binding_hash="$(/usr/bin/sha256sum "$affinity_binding" | \
    /usr/bin/cut -d' ' -f1)"
audit_hash="$(/usr/bin/sha256sum "$lane/audit.json" | \
    /usr/bin/cut -d' ' -f1)"

next_stage final_core_manifest
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
    --arg campaign_manifest_sha256 "$campaign_manifest_hash" \
    --arg campaign_raw_sha256 "$campaign_raw_hash" \
    --arg affinity_report_sha256 "$affinity_report_hash" \
    --arg affinity_binding_sha256 "$affinity_binding_hash" \
    --arg audit_sha256 "$audit_hash" \
    --arg sse2neon_commit "$candidate_submodule_commit" \
    --arg sse2neon_source_archive_sha256 "$sse2neon_source_archive_hash" \
    '{schema:"leopard2-v17-gfni-main-core-manifest/v1",status:"complete",campaign_exit_status:$campaign_exit_status,evidence_valid:true,performance_gate_passed:$performance_gate_passed,promotion_passed:false,promotion_requires_completion_envelope:true,source_commit:$commit,source_tree:$tree,baseline_commit:$main_commit,sse2neon_commit:$sse2neon_commit,sse2neon_source_archive_sha256:$sse2neon_source_archive_sha256,candidate_binary_sha256_pre_post:$candidate_binary_sha256,candidate_archive_sha256_pre_post:$candidate_archive_sha256,baseline_binary_sha256_pre_post:$baseline_binary_sha256,baseline_archive_sha256_pre_post:$baseline_archive_sha256,runner_sha256:$runner_sha256,auditor_sha256:$auditor_sha256,supervisor_sha256:$supervisor_sha256,wrapper_sha256:$wrapper_sha256,campaign_manifest_sha256:$campaign_manifest_sha256,campaign_raw_sha256:$campaign_raw_sha256,affinity_report_sha256:$affinity_report_sha256,affinity_binding_sha256:$affinity_binding_sha256,audit_sha256:$audit_sha256,build_replay_byte_identical:true,producer_verification_passed:true,independent_preseal_audit_passed:true,canonical_lock:"/tmp/leopard-gf8-authoritative.lock",cpu:52,sibling:116,ratio_policy:{ordinary_and_one_shot_are_separate:true,ratios_are_separate_correlated_and_must_not_be_multiplied:true,combined_or_stacked_ratio_emitted:false,same_binary_ratio_is_another_campaign:true},postseal_policy:"promotion requires the enclosing COMPLETED.json written only after byte-identical independent pre/post audits, core SHA verification, clean-source recheck, and a zero campaign exit"}' \
    > "$lane/manifest.json"

next_stage seal_core
trap - ERR
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
/usr/bin/python3 -I -S -B \
    "$lane/audit_v17_gfni_main_compare.py" \
    --manifest "$campaign_dir/manifest.json" \
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

if [[ "$performance_gate_passed" != true ]]; then
    stage=publish_not_promoted_envelope
    /usr/bin/printf 'AUTHORITATIVE_STAGE publish_not_promoted_envelope\n'
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
    /usr/bin/printf '' > "$envelope/SHA256SUMS"
    write_tree_metadata_manifest "$envelope"
    (
        cd "$envelope"
        /usr/bin/find . -type f ! -path './SHA256SUMS' -print0 | \
            /usr/bin/sort -z | \
            /usr/bin/xargs -0 /usr/bin/sha256sum > SHA256SUMS
    )
    trap - ERR
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
trap - ERR
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
exit 0
