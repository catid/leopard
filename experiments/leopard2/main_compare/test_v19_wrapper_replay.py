#!/usr/bin/python3

"""Synthetic sealed-envelope replay for dormant conditioned-passive v19."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[3]
MAIN_COMPARE = ROOT / "experiments/leopard2/main_compare"
WRAPPER = MAIN_COMPARE / "run_authoritative_v17_gfni_main_compare.sh"
REPLAY_PATH = MAIN_COMPARE / "conditioned_v19_wrapper_replay.py"


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load synthetic v19 dependency: {path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


fixtures = load_module(
    "v19_wrapper_replay_fixtures", MAIN_COMPARE / "test_v19_end_to_end.py")
replay = load_module("conditioned_v19_wrapper_replay_tested", REPLAY_PATH)
runner = fixtures.runner
census = fixtures.census
auditor = fixtures.auditor


def command_output(*arguments: str) -> bytes:
    completed = subprocess.run(
        list(arguments), stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        timeout=60)
    if completed.returncode != 0 or completed.stderr:
        raise RuntimeError(
            f"fixture command failed: {arguments!r}: "
            f"{completed.stderr.decode(errors='replace')}")
    return completed.stdout


def git_identity(expression: str) -> str:
    return command_output(
        "/usr/bin/git", "-C", str(ROOT), "rev-parse", "--verify",
        expression).decode("ascii").strip()


CONTROLLER_SOURCE_COMMIT = git_identity("HEAD^{commit}")
CONTROLLER_SOURCE_TREE = git_identity("HEAD^{tree}")
CANDIDATE_SOURCE_COMMIT = "cf7a7056e0bd7f54b8da436a39cae857beab10c1"
CANDIDATE_SOURCE_TREE = "499697b4ee71ce7b2a12831f7d72d0e423cd0806"


def write_bytes(path: Path, data: bytes, mode: int = 0o600) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor = os.open(
        path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC, mode)
    try:
        os.fchmod(descriptor, mode)
        view = memoryview(data)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise RuntimeError(f"short fixture write: {path}")
            view = view[written:]
    finally:
        os.close(descriptor)


def write_json(path: Path, value: Any, mode: int = 0o600) -> None:
    write_bytes(path, replay.canonical_json_bytes(value), mode)


def replace_json(path: Path, value: Any) -> None:
    path.write_bytes(replay.canonical_json_bytes(value))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_source(source: Path, destination: Path, mode: int) -> None:
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    shutil.copyfile(source, destination, follow_symlinks=False)
    destination.chmod(mode)


def common_sources(core: Path) -> None:
    sources = (
        (WRAPPER, core / "run-authoritative.sh", 0o700),
        (REPLAY_PATH, core / "conditioned_v19_wrapper_replay.py", 0o700),
        (MAIN_COMPARE / "run_abba.py", core / "run_abba.py", 0o600),
        (MAIN_COMPARE / "git_capture.py", core / "git_capture.py", 0o600),
        (ROOT / "tools/leopard2_build_provenance.py",
         core / "leopard2_build_provenance.py", 0o600),
        (ROOT / "experiments/leopard2/decoder_dispatch/"
         "balanced_evidence_common.py",
         core / "balanced_evidence_common.py", 0o600),
        (MAIN_COMPARE / "audit_v17_gfni_main_compare.py",
         core / "audit_v17_gfni_main_compare.py", 0o600),
        (MAIN_COMPARE / "passive_environment_census.py",
         core / "passive_environment_census.py", 0o600),
    )
    for source, destination, mode in sources:
        copy_source(source, destination, mode)
    for lane_path, source_path in replay.CONTROLLER_BINDINGS:
        copy_source(ROOT / source_path, core / lane_path, 0o600)


def wrapper_records(core: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    preregistration_data = command_output(
        str(WRAPPER), "--print-conditioned-v19-preregistration")
    closure_data = command_output(
        str(WRAPPER), "--print-conditioned-v19-controller-closure")
    write_bytes(core / "conditioned-v19-preregistration.json",
                preregistration_data)
    write_bytes(core / "conditioned-v19-controller-closure.json",
                closure_data)
    return json.loads(preregistration_data), json.loads(closure_data)


def emitted_record(
    kind: str, envelope: Path, request: dict[str, Any] | None = None,
) -> dict[str, Any]:
    arguments = [
        str(WRAPPER), "--emit-conditioned-v19-record", kind, str(envelope),
    ]
    request_path = envelope.parent / f"{envelope.name}-record-request.json"
    if request is not None:
        write_json(request_path, request)
        arguments.append(str(request_path))
    try:
        output = command_output(*arguments)
    finally:
        if request_path.exists():
            request_path.unlink()
    return json.loads(output)


def attempt_two_not_acquired(prior_failure_sha256: str) -> dict[str, Any]:
    return census.v19_attempt_record(
        attempt=2, prior_attempt_failure_sha256=prior_failure_sha256,
        prior_attempt_selection_status="not-acquired")


def attempt_lineage(
    envelope: Path, authority: dict[str, Any],
    prior_attempts: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    previous = [] if prior_attempts is None else copy.deepcopy(prior_attempts)
    return emitted_record("lineage", envelope, {
        "controller_source_commit": CONTROLLER_SOURCE_COMMIT,
        "controller_source_tree": CONTROLLER_SOURCE_TREE,
        "prior_attempts": previous,
        "external_attempt": copy.deepcopy(authority),
    })


def envelope_path(root: Path, attempt: int) -> Path:
    return root / \
        f"{CANDIDATE_SOURCE_COMMIT[:7]}-v19-conditioned-main-a{attempt}"


def pair_verdict(raw: dict[str, Any]) -> dict[str, Any]:
    qualification = raw["pair_qualification"]
    contract = runner.pair_v19.contract
    return runner.pair_v19.verifier.require_accepted_pair_qualification_bundle(
        contract.canonical_json_bytes(qualification["acquisition"]),
        contract.canonical_json_bytes(qualification["bridge"]),
        expected_policy=qualification["policy"],
        expected_policy_sha256=qualification["policy_sha256"],
        expected_frozen_pair=
            qualification["attempt"]["frozen_pair_from_prior_attempt"],
        expected_acquisition_window_count=
            runner.pair_v19.QUALIFICATION_WINDOW_COUNT,
        expected_acquisition_nominal_window_ns=
            runner.pair_v19.QUALIFICATION_NOMINAL_WINDOW_NS,
        expected_bridge_geometry=runner.pair_v19.bridge_geometry_record())


def core_checksum_paths(core: Path) -> list[str]:
    return sorted(
        (path.relative_to(core).as_posix() for path in core.rglob("*")
         if path.is_file() and path.name != "SHA256SUMS"),
        key=lambda value: value.encode("utf-8"))


def checksum_lines(root: Path, paths: list[str]) -> bytes:
    return "".join(
        f"{sha256(root / relative)}  ./{relative}\n" for relative in paths
    ).encode("ascii")


def tree_metadata(root: Path) -> dict[str, Any]:
    root_stat = root.lstat()
    records = []
    for path in (root, *sorted(root.rglob("*"),
                               key=lambda item: item.as_posix())):
        relative = "." if path == root else path.relative_to(root).as_posix()
        if relative == "TREE-METADATA.json":
            continue
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            kind = "directory"
        elif stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1:
            kind = "file"
        else:
            raise RuntimeError(f"unsafe synthetic tree node: {path}")
        records.append({
            "gid": metadata.st_gid,
            "mode": format(stat.S_IMODE(metadata.st_mode) & ~0o222, "04o"),
            "nlink": metadata.st_nlink,
            "path": relative,
            "type": kind,
            "uid": metadata.st_uid,
        })
    records.sort(key=lambda record: record["path"])
    return {
        "entries": records,
        "excluded_paths": ["TREE-METADATA.json"],
        "final_mode_policy": "observed mode with all write bits removed",
        "root": ".",
        "schema": "leopard2-authoritative-tree-metadata/v1",
        "self_policy": {
            "gid": root_stat.st_gid,
            "mode": "0400",
            "nlink": 1,
            "sha256_binding":
                "exactly one ./TREE-METADATA.json checksum entry",
            "type": "file",
            "uid": root_stat.st_uid,
        },
        "uid_gid_policy": {
            "gid": root_stat.st_gid,
            "rule": "every retained node has the invoking effective uid and gid",
            "uid": root_stat.st_uid,
        },
    }


def unseal(envelope: Path) -> None:
    for path in (envelope, *envelope.rglob("*")):
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            path.chmod(0o700)
        elif stat.S_ISREG(metadata.st_mode):
            executable = bool(stat.S_IMODE(metadata.st_mode) & 0o111)
            path.chmod(0o700 if executable else 0o600)


def seal(envelope: Path, *, failure_exit_status: int | None = None) -> None:
    core = envelope / "core"
    for path in (envelope / "TREE-METADATA.json",
                 envelope / "SHA256SUMS", core / "SHA256SUMS"):
        if path.exists():
            path.unlink()
    success_path = envelope / "NOT_PROMOTED.json"
    failure_path = envelope / "FAILED.json"
    if success_path.exists() and failure_path.exists():
        raise RuntimeError("synthetic envelope has two terminals")
    terminal_path = success_path if success_path.exists() else failure_path
    terminal = (
        json.loads(terminal_path.read_bytes()) if terminal_path.exists() else
        None)
    write_bytes(core / "SHA256SUMS",
                checksum_lines(core, core_checksum_paths(core)))
    if terminal is None:
        manifest = json.loads((core / "manifest.json").read_bytes())
        if manifest["status"] == "complete":
            terminal_path = success_path
            terminal = emitted_record("success-terminal", envelope)
        else:
            if failure_exit_status is None:
                raise RuntimeError("synthetic failure lacks its shell exit status")
            terminal_path = failure_path
            terminal = emitted_record(
                "failure-terminal", envelope,
                {"campaign_exit_status": failure_exit_status})
        write_json(terminal_path, terminal)
    else:
        terminal["core_manifest_sha256"] = sha256(core / "manifest.json")
        terminal["core_sha256sums_sha256"] = sha256(core / "SHA256SUMS")
        replace_json(terminal_path, terminal)
    outer_paths = [terminal_path.name]
    if (envelope / "postseal-audit.json").exists():
        outer_paths.append("postseal-audit.json")
    outer_paths.append("core/SHA256SUMS")
    write_bytes(envelope / "SHA256SUMS", checksum_lines(envelope, outer_paths))
    write_json(envelope / "TREE-METADATA.json", tree_metadata(envelope))
    with (envelope / "SHA256SUMS").open("ab") as stream:
        stream.write(
            f"{sha256(envelope / 'TREE-METADATA.json')}  "
            "./TREE-METADATA.json\n".encode("ascii"))
    for path in sorted(envelope.rglob("*"), reverse=True):
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode):
            path.chmod(stat.S_IMODE(metadata.st_mode) & ~0o222)
        elif stat.S_ISREG(metadata.st_mode):
            path.chmod(stat.S_IMODE(metadata.st_mode) & ~0o222)
    envelope.chmod(stat.S_IMODE(envelope.lstat().st_mode) & ~0o222)


def reseal_with_retained_tree_metadata(envelope: Path) -> None:
    terminal = (
        "NOT_PROMOTED.json" if (envelope / "NOT_PROMOTED.json").exists()
        else "FAILED.json")
    outer_paths = [terminal]
    if (envelope / "postseal-audit.json").exists():
        outer_paths.append("postseal-audit.json")
    outer_paths.extend(("core/SHA256SUMS", "TREE-METADATA.json"))
    sums = envelope / "SHA256SUMS"
    if sums.exists():
        sums.unlink()
    write_bytes(sums, checksum_lines(envelope, outer_paths))
    for path in sorted(envelope.rglob("*"), reverse=True):
        metadata = path.lstat()
        if stat.S_ISDIR(metadata.st_mode) or stat.S_ISREG(metadata.st_mode):
            path.chmod(stat.S_IMODE(metadata.st_mode) & ~0o222)
    envelope.chmod(stat.S_IMODE(envelope.lstat().st_mode) & ~0o222)


def success_core_manifest(envelope: Path) -> dict[str, Any]:
    return emitted_record("complete-core", envelope)


def build_success(envelope: Path, raw: dict[str, Any] | None = None) -> Path:
    envelope.mkdir(mode=0o700)
    core = envelope / "core"
    core.mkdir(mode=0o700)
    common_sources(core)
    _preregistration, _closure = wrapper_records(core)
    value = fixtures.canonical_raw() if raw is None else copy.deepcopy(raw)
    authority = copy.deepcopy(value["pair_qualification"]["attempt"])
    campaign_manifest = \
        fixtures.replay_fixtures.materialize_windowed_auditor_bundle(
            core / "campaign", value)
    retained_raw = json.loads((core / "campaign/raw.json").read_bytes())
    pre, post, controller = fixtures.census_witnesses(retained_raw)
    policy = fixtures.compare_census(retained_raw, authority, (pre, post, controller))
    audit = auditor.replay(
        campaign_manifest, supervision_mode="conditioned",
        expected_v19_attempt=authority)
    write_json(core / "environment-census-pre.json", pre)
    write_json(core / "environment-census-post.json", post)
    write_json(core / "controller-affinity.json", controller)
    write_json(core / "passive-environment-policy.json", policy)
    write_json(core / "audit.json", audit)
    write_json(envelope / "postseal-audit.json", audit)
    qualification = retained_raw["pair_qualification"]
    write_json(core / "pair-qualification-policy.json",
               qualification["policy"])
    write_json(core / "pair-qualification-acquisition.json",
               qualification["acquisition"])
    write_json(core / "pair-qualification-bridge.json",
               qualification["bridge"])
    write_json(core / "pair-qualification-verdict.json",
               pair_verdict(retained_raw))
    lineage = attempt_lineage(envelope, authority)
    write_json(core / "attempt-lineage.json", lineage)
    manifest = success_core_manifest(envelope)
    write_json(core / "manifest.json", manifest)
    seal(envelope)
    return envelope


def build_failure(
    envelope: Path, *, authority: dict[str, Any] | None = None,
    prior_attempts: list[dict[str, Any]] | None = None,
    campaign_exit_status: int = 1,
    selected_with_streams: bool = False,
) -> Path:
    if authority is None:
        prior_envelope = envelope_path(envelope.parent, 1)
        build_failure(
            prior_envelope,
            authority=census.v19_attempt_record(attempt=1),
            prior_attempts=[])
        prior_terminal = prior_envelope / "FAILED.json"
        prior_attempts = [{
            "attempt": 1,
            "envelope": str(prior_envelope),
            "terminal": "FAILED.json",
            "terminal_schema": replay.FAILURE_TERMINAL_SCHEMA,
            "envelope_sha256sums_sha256":
                sha256(prior_envelope / "SHA256SUMS"),
            "terminal_sha256": sha256(prior_terminal),
        }]
        authority = attempt_two_not_acquired(sha256(prior_terminal))
    envelope.mkdir(mode=0o700)
    core = envelope / "core"
    core.mkdir(mode=0o700)
    campaign = core / "campaign"
    campaign.mkdir(mode=0o700)
    common_sources(core)
    _preregistration, _closure = wrapper_records(core)
    fixture_module = fixtures.replay_fixtures.fixtures
    if selected_with_streams:
        raw = fixture_module.synthetic_raw(raw_schema=runner.RAW_SCHEMA_V19)
        failure = fixture_module.synthetic_failure(runner.RAW_SCHEMA_V19)
        invocation = copy.deepcopy(raw["invocations"][0])
        stdout_data = runner.canonical_bytes(invocation["result"])
        stderr_data = b""
        write_bytes(campaign / "first.stdout", stdout_data)
        write_bytes(campaign / "first.stderr", stderr_data)
        invocation["stdout"] = {
            "path": "first.stdout", "size": len(stdout_data),
            "sha256": runner.sha256_bytes(stdout_data),
        }
        invocation["stderr"] = {
            "path": "first.stderr", "size": len(stderr_data),
            "sha256": runner.sha256_bytes(stderr_data),
        }
        failure["invocations"] = [invocation]
        isolation = raw["isolation"]
        failure["isolation"] = runner.isolation_record_v2(
            failure["campaign"]["benchmark_cpu"],
            failure["campaign"]["reserved_sibling"],
            failure["pair_lease"], isolation["before"]["monotonic_ns"],
            isolation["after"]["monotonic_ns"],
            isolation["before"]["benchmark_cpu"],
            isolation["after"]["benchmark_cpu"],
            isolation["before"]["reserved_sibling"],
            isolation["after"]["reserved_sibling"],
            [invocation["cpu_window"]])
        failure["retained_files"] = runner.retained_file_records(campaign)
    else:
        failure = fixture_module.synthetic_unselected_v19_failure(
            acquired=False)
    failure["pair_qualification"]["attempt"] = copy.deepcopy(authority)
    failure = fixture_module.resign(failure)
    write_json(campaign / "failure.json", failure)
    lineage = attempt_lineage(envelope, authority, prior_attempts)
    write_json(core / "attempt-lineage.json", lineage)
    manifest = emitted_record(
        "failed-core", envelope,
        {"campaign_exit_status": campaign_exit_status})
    write_json(core / "manifest.json", manifest)
    seal(envelope, failure_exit_status=campaign_exit_status)
    return envelope


def run_wrapper(envelope: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(WRAPPER), "--verify-conditioned-v19-campaign-core",
         str(envelope)],
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, text=True, timeout=120)


def mutate_reseal(envelope: Path, mutation: Callable[[Path], None]) -> None:
    unseal(envelope)
    mutation(envelope)
    seal(envelope)


def mutate_campaign_raw(
    envelope: Path, mutation: Callable[[dict[str, Any]], None],
) -> None:
    campaign = envelope / "core/campaign"
    raw = json.loads((campaign / "raw.json").read_bytes())
    mutation(raw)
    shutil.rmtree(campaign)
    fixtures.replay_fixtures.materialize_windowed_auditor_bundle(
        campaign, raw)
    manifest_path = envelope / "core/manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["campaign_manifest_sha256"] = sha256(campaign / "manifest.json")
    manifest["campaign_raw_sha256"] = sha256(campaign / "raw.json")
    replace_json(manifest_path, manifest)


class V19WrapperReplayTests(unittest.TestCase):
    def test_synthetic_success_and_failure_replay(self) -> None:
        for kind, builder in (("success", build_success),
                              ("failure", build_failure)):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory(
                    prefix="leopard-v19-wrapper-replay.") as temporary:
                if kind == "success":
                    envelope = builder(envelope_path(Path(temporary), 1))
                else:
                    envelope = builder(envelope_path(Path(temporary), 2))
                completed = run_wrapper(envelope)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                self.assertEqual(
                    completed.stdout,
                    f"conditioned v19 campaign core verified: {envelope}\n")
                self.assertEqual(completed.stderr, "")

    def test_wrapper_level_splices_fail_closed(self) -> None:
        def authority(envelope: Path) -> None:
            path = envelope / "core/attempt-lineage.json"
            value = json.loads(path.read_bytes())
            value["external_attempt"] = attempt_two_not_acquired("0" * 64)
            replace_json(path, value)
            digest = sha256(path)
            for record_path in (
                    envelope / "core/manifest.json",
                    envelope / "NOT_PROMOTED.json"):
                record = json.loads(record_path.read_bytes())
                record["attempt_lineage_sha256"] = digest
                replace_json(record_path, record)

        def controller(envelope: Path) -> None:
            path = envelope / "core/controller-affinity.json"
            value = json.loads(path.read_bytes())
            value["active_affinity_supervisor_executed"] = True
            replace_json(path, value)
            manifest_path = envelope / "core/manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            manifest["controller_affinity_sha256"] = sha256(path)
            replace_json(manifest_path, manifest)

        def pair_splice(envelope: Path) -> None:
            for path in (envelope / "core/manifest.json",
                         envelope / "NOT_PROMOTED.json"):
                value = json.loads(path.read_bytes())
                value["selected_pair"]["reserved_sibling"] = 66
                replace_json(path, value)

        def bridge(envelope: Path) -> None:
            path = envelope / "core/pair-qualification-bridge.json"
            value = json.loads(path.read_bytes())
            value["bridge_tail_sha256"] = "0" * 64
            replace_json(path, value)
            manifest_path = envelope / "core/manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            manifest["pair_qualification_bridge_sha256"] = sha256(path)
            replace_json(manifest_path, manifest)

        def affinity(envelope: Path) -> None:
            path = envelope / "core/controller-affinity.json"
            value = json.loads(path.read_bytes())
            value["before_allowed_cpus"].append(128)
            value["runner_launch_allowed_cpus"].append(128)
            replace_json(path, value)
            manifest_path = envelope / "core/manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            manifest["controller_affinity_sha256"] = sha256(path)
            replace_json(manifest_path, manifest)

        def resource(envelope: Path) -> None:
            path = envelope / "core/conditioned-v19-preregistration.json"
            value = json.loads(path.read_bytes())
            value["resource_envelope"]["release_max_jobs"] = 2
            replace_json(path, value)

        def handoff(envelope: Path) -> None:
            def splice(raw: dict[str, Any]) -> None:
                raw["pair_qualification"]["first_window_handoff"] \
                    ["handoff_elapsed_ns"] = 0

            mutate_campaign_raw(envelope, splice)

        def lineage(envelope: Path) -> None:
            path = envelope / "core/attempt-lineage.json"
            value = json.loads(path.read_bytes())
            value["v18_failure_lineage_sha256"] = "0" * 64
            replace_json(path, value)
            digest = sha256(path)
            for record_path in (
                    envelope / "core/manifest.json",
                    envelope / "NOT_PROMOTED.json"):
                record = json.loads(record_path.read_bytes())
                record["attempt_lineage_sha256"] = digest
                replace_json(record_path, record)

        def inner_v18_lineage(envelope: Path) -> None:
            def splice(raw: dict[str, Any]) -> None:
                raw["pair_qualification"]["v18_failure_lineage"] \
                    ["attempts"][0]["envelope_sha256sums_sha256"] = \
                    "0" * 64

            mutate_campaign_raw(envelope, splice)

        def claim(envelope: Path) -> None:
            for path in (envelope / "core/manifest.json",
                         envelope / "NOT_PROMOTED.json"):
                value = json.loads(path.read_bytes())
                value["shared_host_claim_ceiling"]["promotion_eligible"] = \
                    True
                replace_json(path, value)

        def terminal(envelope: Path) -> None:
            path = envelope / "NOT_PROMOTED.json"
            value = json.loads(path.read_bytes())
            value["unexpected"] = False
            replace_json(path, value)

        def hash_binding(envelope: Path) -> None:
            path = envelope / "core/manifest.json"
            value = json.loads(path.read_bytes())
            value["audit_sha256"] = "0" * 64
            replace_json(path, value)

        def exit_status_type(envelope: Path) -> None:
            for path in (envelope / "core/manifest.json",
                         envelope / "NOT_PROMOTED.json"):
                value = json.loads(path.read_bytes())
                value["campaign_exit_status"] = False
                replace_json(path, value)

        def extra_campaign_directory(envelope: Path) -> None:
            (envelope / "core/campaign/unbound-empty").mkdir(mode=0o700)

        def extra_campaign_file(envelope: Path) -> None:
            write_bytes(envelope / "core/campaign/unbound.bin", b"unbound")

        def candidate_source_cross_splice(envelope: Path) -> None:
            lineage_path = envelope / "core/attempt-lineage.json"
            lineage = json.loads(lineage_path.read_bytes())
            lineage["source_commit"] = CONTROLLER_SOURCE_COMMIT
            lineage["source_tree"] = CONTROLLER_SOURCE_TREE
            replace_json(lineage_path, lineage)
            lineage_hash = sha256(lineage_path)
            for record_path in (
                    envelope / "core/manifest.json",
                    envelope / "NOT_PROMOTED.json"):
                record = json.loads(record_path.read_bytes())
                record["source_commit"] = CONTROLLER_SOURCE_COMMIT
                record["source_tree"] = CONTROLLER_SOURCE_TREE
                record["attempt_lineage_sha256"] = lineage_hash
                replace_json(record_path, record)

        cases = (
            ("authority", authority),
            ("controller", controller),
            ("pair", pair_splice),
            ("bridge", bridge),
            ("affinity", affinity),
            ("handoff", handoff),
            ("lineage", lineage),
            ("inner-v18-lineage", inner_v18_lineage),
            ("claim", claim),
            ("resource", resource),
            ("terminal", terminal),
            ("hash", hash_binding),
            ("exit-status-type", exit_status_type),
            ("extra-campaign-directory", extra_campaign_directory),
            ("extra-campaign-file", extra_campaign_file),
            ("candidate-source-cross-splice", candidate_source_cross_splice),
        )
        for label, mutation in cases:
            with self.subTest(splice=label), tempfile.TemporaryDirectory(
                    prefix="leopard-v19-wrapper-replay.") as temporary:
                envelope = build_success(envelope_path(Path(temporary), 1))
                mutate_reseal(envelope, mutation)
                completed = run_wrapper(envelope)
                self.assertNotEqual(completed.returncode, 0)
                self.assertEqual(completed.stdout, "")

    def test_failure_join_and_mode_splices_fail_closed(self) -> None:
        def attempt_join(envelope: Path) -> None:
            terminal_path = envelope / "FAILED.json"
            terminal = json.loads(terminal_path.read_bytes())
            terminal["attempt"] = 1
            replace_json(terminal_path, terminal)

        def exit_status_type(envelope: Path) -> None:
            for path in (envelope / "core/manifest.json",
                         envelope / "FAILED.json"):
                value = json.loads(path.read_bytes())
                value["campaign_exit_status"] = True
                replace_json(path, value)

        def verify_status_type(envelope: Path) -> None:
            path = envelope / "core/manifest.json"
            value = json.loads(path.read_bytes())
            value["failure_verify_status"] = False
            replace_json(path, value)

        def prior_checksum(envelope: Path) -> None:
            path = envelope / "core/attempt-lineage.json"
            value = json.loads(path.read_bytes())
            value["prior_attempts"][0]["envelope_sha256sums_sha256"] = \
                "2" * 64
            replace_json(path, value)
            digest = sha256(path)
            for record_path in (
                    envelope / "core/manifest.json",
                    envelope / "FAILED.json"):
                record = json.loads(record_path.read_bytes())
                record["attempt_lineage_sha256"] = digest
                replace_json(record_path, record)

        def prior_attempt_type(envelope: Path) -> None:
            path = envelope / "core/attempt-lineage.json"
            value = json.loads(path.read_bytes())
            value["prior_attempts"][0]["attempt"] = True
            replace_json(path, value)
            digest = sha256(path)
            for record_path in (
                    envelope / "core/manifest.json",
                    envelope / "FAILED.json"):
                record = json.loads(record_path.read_bytes())
                record["attempt_lineage_sha256"] = digest
                replace_json(record_path, record)

        def mismatched_exit_status(envelope: Path) -> None:
            path = envelope / "FAILED.json"
            value = json.loads(path.read_bytes())
            value["campaign_exit_status"] = 8
            replace_json(path, value)

        def bounded_exit_status(envelope: Path, value: int) -> None:
            for path in (envelope / "core/manifest.json",
                         envelope / "FAILED.json"):
                record = json.loads(path.read_bytes())
                record["campaign_exit_status"] = value
                replace_json(path, record)

        for label, mutation in (
                ("attempt-join", attempt_join),
                ("exit-status-type", exit_status_type),
                ("verify-status-type", verify_status_type),
                ("prior-checksum", prior_checksum),
                ("prior-attempt-type", prior_attempt_type),
                ("mismatched-exit-status", mismatched_exit_status),
                ("zero-exit-status",
                 lambda envelope: bounded_exit_status(envelope, 0)),
                ("oversized-exit-status",
                 lambda envelope: bounded_exit_status(envelope, 256))):
            with self.subTest(splice=label), tempfile.TemporaryDirectory(
                    prefix="leopard-v19-wrapper-replay.") as temporary:
                envelope = build_failure(envelope_path(Path(temporary), 2))
                mutate_reseal(envelope, mutation)
                completed = run_wrapper(envelope)
                self.assertNotEqual(completed.returncode, 0)
                self.assertEqual(completed.stdout, "")

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            envelope = build_success(envelope_path(Path(temporary), 1))
            target = envelope / "core/pair_qualification_contract.py"
            target.chmod(0o440)
            completed = run_wrapper(envelope)
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, "")

    def test_selected_failure_stream_inventory_and_exit_status(self) -> None:
        def selected_failure(root: Path) -> Path:
            return build_failure(
                envelope_path(root, 1),
                authority=census.v19_attempt_record(attempt=1),
                prior_attempts=[], campaign_exit_status=7,
                selected_with_streams=True)

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            envelope = selected_failure(Path(temporary))
            completed = run_wrapper(envelope)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            terminal = json.loads((envelope / "FAILED.json").read_bytes())
            manifest = json.loads(
                (envelope / "core/manifest.json").read_bytes())
            self.assertEqual(terminal["campaign_exit_status"], 7)
            self.assertEqual(manifest["campaign_exit_status"], 7)
            self.assertEqual(terminal["selected_pair"], {
                "benchmark_cpu": 1, "reserved_sibling": 65,
            })

        def missing_listed(envelope: Path) -> None:
            (envelope / "core/campaign/first.stdout").unlink()

        def changed_listed(envelope: Path) -> None:
            path = envelope / "core/campaign/first.stdout"
            path.write_bytes(path.read_bytes() + b"changed")

        def unlisted_file(envelope: Path) -> None:
            write_bytes(envelope / "core/campaign/unlisted.bin", b"unlisted")

        def unlisted_directory(envelope: Path) -> None:
            (envelope / "core/campaign/unlisted-empty").mkdir(mode=0o700)

        for label, mutation in (
                ("missing-listed", missing_listed),
                ("changed-listed", changed_listed),
                ("unlisted-file", unlisted_file),
                ("unlisted-directory", unlisted_directory)):
            with self.subTest(splice=label), tempfile.TemporaryDirectory(
                    prefix="leopard-v19-wrapper-replay.") as temporary:
                envelope = selected_failure(Path(temporary))
                mutate_reseal(envelope, mutation)
                completed = run_wrapper(envelope)
                self.assertNotEqual(completed.returncode, 0)
                self.assertEqual(completed.stdout, "")

    def test_prior_root_and_tree_metadata_splices_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            root = Path(temporary)
            envelope = build_failure(envelope_path(root, 2))
            prior = envelope_path(root, 1)
            displaced = root / "displaced-prior-envelope"
            prior.rename(displaced)
            prior.symlink_to(displaced, target_is_directory=True)
            completed = run_wrapper(envelope)
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, "")

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            root = Path(temporary)
            envelope = build_failure(envelope_path(root, 2))
            prior = envelope_path(root, 1)
            unseal(prior)
            metadata_path = prior / "TREE-METADATA.json"
            metadata = json.loads(metadata_path.read_bytes())
            metadata["entries"][0]["mode"] = "0700"
            replace_json(metadata_path, metadata)
            reseal_with_retained_tree_metadata(prior)
            completed = run_wrapper(envelope)
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, "")

    def test_conditioned_failure_paths_dominate_legacy_emitters(self) -> None:
        source = WRAPPER.read_text(encoding="utf-8")
        for function_name in ("failure_record", "postseal_failure_record"):
            start = source.index(f"{function_name}()")
            legacy = source.index(
                'leopard2-v18-gfni-main-failed-envelope/v1', start)
            refusal = source.index(
                "refusing legacy v18 failure emission", start, legacy)
            self.assertLess(refusal, legacy)
        campaign = source.index('if [[ "$campaign_status" -ne 0 ]]')
        next_legacy = source.index(
            'leopard2-v18-gfni-main-passive-failed-core-manifest/v1', campaign)
        refusal = source.index(
            "refusing legacy v18 failure emission", campaign, next_legacy)
        self.assertLess(refusal, next_legacy)

    def test_controller_source_authority_is_exact(self) -> None:
        def verify(core: Path, commit: str, tree: str) -> \
                subprocess.CompletedProcess[str]:
            return subprocess.run(
                [str(WRAPPER),
                 "--verify-conditioned-v19-core-source-authority",
                 str(core), commit, tree],
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False, text=True, timeout=60)

        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            envelope = build_success(envelope_path(Path(temporary), 1))
            core = envelope / "core"
            current = verify(
                core, CONTROLLER_SOURCE_COMMIT, CONTROLLER_SOURCE_TREE)
            self.assertEqual(current.returncode, 0, current.stderr)
            self.assertEqual(current.stdout, "")

            wrong_tree = verify(
                core, CONTROLLER_SOURCE_COMMIT, CANDIDATE_SOURCE_TREE)
            self.assertNotEqual(wrong_tree.returncode, 0)
            wrong_commit = verify(
                core, CANDIDATE_SOURCE_COMMIT, CANDIDATE_SOURCE_TREE)
            self.assertNotEqual(wrong_commit.returncode, 0)

            retained = core / "run_abba.py"
            retained.chmod(0o600)
            retained.write_bytes(retained.read_bytes() + b"\n")
            retained.chmod(0o400)
            changed = verify(
                core, CONTROLLER_SOURCE_COMMIT, CONTROLLER_SOURCE_TREE)
            self.assertNotEqual(changed.returncode, 0)

    def test_retained_source_is_data_not_executable(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard-v19-wrapper-replay.") as temporary:
            root = Path(temporary)
            envelope = build_success(envelope_path(root, 1))
            marker = root / "retained-source-executed"

            def replace_controller(target: Path) -> None:
                source = target / "core/pair_qualified_v19_contract.py"
                source.write_text(
                    "from pathlib import Path\n"
                    f"Path({str(marker)!r}).write_text('executed')\n",
                    encoding="utf-8")
                closure_path = target / \
                    "core/conditioned-v19-controller-closure.json"
                closure = json.loads(closure_path.read_bytes())
                for entry in closure["files"]:
                    if entry["lane_path"] == source.name:
                        entry["sha256"] = sha256(source)
                replace_json(closure_path, closure)
                closure_hash = sha256(closure_path)
                lineage_path = target / "core/attempt-lineage.json"
                lineage = json.loads(lineage_path.read_bytes())
                lineage["controller_closure_sha256"] = closure_hash
                replace_json(lineage_path, lineage)
                lineage_hash = sha256(lineage_path)
                for record_path in (
                        target / "core/manifest.json",
                        target / "NOT_PROMOTED.json"):
                    record = json.loads(record_path.read_bytes())
                    record["controller_closure_sha256"] = closure_hash
                    record["attempt_lineage_sha256"] = lineage_hash
                    replace_json(record_path, record)

            mutate_reseal(envelope, replace_controller)
            completed = run_wrapper(envelope)
            self.assertNotEqual(completed.returncode, 0)
            self.assertEqual(completed.stdout, "")
            self.assertFalse(marker.exists())

    def test_preflight_identity_join_is_exact(self) -> None:
        preregistration = json.loads(command_output(
            str(WRAPPER), "--print-conditioned-v19-preregistration"))
        preflight = preregistration["build_preflight"]
        identity = {
            "candidate_source": {
                "head": preflight["source_commit"],
                "tree": preflight["source_tree"],
            },
            "candidate_executable": {
                "sha256": preflight["candidate_binary_sha256"]},
            "candidate_archive": {
                "sha256": preflight["candidate_archive_sha256"]},
            "baseline_source": {
                "head": preflight["baseline_commit"],
                "tree": preflight["baseline_tree"],
            },
            "baseline_executable": {
                "sha256": preflight["baseline_binary_sha256"]},
            "baseline_archive": {
                "sha256": preflight["baseline_archive_sha256"]},
        }
        raw = {
            "identities_initial": copy.deepcopy(identity),
            "identities_final": copy.deepcopy(identity),
            "invocations": [{
                "identity_before": copy.deepcopy(identity),
                "identity_after": copy.deepcopy(identity),
            }],
        }
        campaign_manifest = {"identities": copy.deepcopy(identity)}
        replay.validate_preflight_identity(
            raw, campaign_manifest, preregistration)
        raw["identities_final"]["candidate_executable"]["sha256"] = \
            "0" * 64
        with self.assertRaises(replay.ReplayError):
            replay.validate_preflight_identity(
                raw, campaign_manifest, preregistration)


if __name__ == "__main__":
    unittest.main()
