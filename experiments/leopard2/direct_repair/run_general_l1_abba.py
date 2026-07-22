#!/usr/bin/env python3
"""Authoritative same-source screen for the general GF8/AVX2 L=1 repair.

The control and candidate are compiled from one clean commit.  By default,
their only permitted compile-command and object difference is the
source-scoped LEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1 definition on
leopard2.cpp.  --pair-fusion additionally permits its one definition on the
backend qualification and AVX2 implementation sources, and no other delta.
--pair-attribution instead compiles direct repair into both executables and
permits only the pair-fusion backend sources to differ.  This isolates the
two-source kernel and exposes any byte range where a runtime threshold is
needed.
The runner freezes both executables before timing, forces and attests AVX2,
uses serial ABBA rounds on one physical core, and reports decode-plan setup
separately from byte-heavy execution.  An exact Leopard-main executable is
invoked only for valid legacy-high cells that first clear the same-source
promotion gate.

This is deliberately a dedicated experiment runner.  It does not change the
historical Leopard-main or expanded-final-map campaigns.
"""

from __future__ import annotations

import argparse
import copy
import contextlib
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Any, Iterable, Mapping, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPOSITORY_ROOT / "tools"))

from leopard2_build_provenance import (  # noqa: E402
    candidate_build_provenance,
)


SCHEMA = "leopard2-general-direct-l1-abba/v2"
MATRIX_SCHEMA = "leopard2-general-direct-l1-matrix/v1"
COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
EXPERIMENT = "LEO2_EXPERIMENTAL_GF8_AVX2_GENERAL_DIRECT_L1"
PAIR_FUSION = "LEO2_EXPERIMENTAL_GF8_AVX2_DIRECT_L1_PAIR_FUSION"
TRANSFORM_VS_SIMPLE = "transform-vs-simple"
TRANSFORM_VS_PAIR = "transform-vs-pair"
SIMPLE_VS_PAIR = "simple-vs-pair"
ABBA = ("control", "candidate", "candidate", "control")
BASE_BYTES = (64, 2048, 4096, 65536, 1048576)
BASE_REUSE = (1, 8, 64)
BYTE_NEIGHBORS = (
    63, 65, 2047, 2049, 4095, 4097,
    65535, 65537, 1048575, 1048577,
)
DEFAULT_LOCK = "/tmp/leopard-gf8-authoritative.lock"

# Required cells from the experiment bead.
HIGH_CORE = (
    (240, 16), (224, 32), (192, 64), (200, 30),
)
LOW_CORE = (
    (17, 31), (31, 200), (32, 224), (127, 128),
)

# Stay on the GF8 side of each dyadic boundary.  Some requested core cells
# sit exactly at N=256, so their +K or +R neighbor is intentionally invalid.
HIGH_NEIGHBORS = (
    (239, 16), (240, 15),
    (223, 32), (224, 31),
    (191, 64), (192, 63),
    (199, 30), (201, 30), (200, 29), (200, 31),
    # Explicit legacy-high is legal above rate 1/2 even though AUTO chooses
    # LOW there.  These cover the two shortened high-profile correctness
    # shapes selected by the focused API test; exact Leopard main is invalid.
    (17, 33), (17, 128),
)
LOW_NEIGHBORS = (
    (18, 31), (17, 30), (17, 32),
    (30, 200), (32, 200), (31, 199), (31, 201),
    (31, 224), (32, 223),
    (126, 128),
)
LOW_BALANCED_NEIGHBORS = (
    (17, 17), (31, 31), (32, 32), (127, 127), (128, 128),
)


class CampaignError(RuntimeError):
    """The campaign cannot produce valid evidence."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CampaignError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path | str, value: Any) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(
        destination.name + ".tmp.%d" % os.getpid())
    temporary.write_bytes(canonical_bytes(value) + b"\n")
    os.replace(str(temporary), str(destination))


def read_json(path: Path | str) -> Any:
    with open(path, "r", encoding="utf-8") as stream:
        return json.load(stream)


def ceil_pow2(value: int) -> int:
    require(value > 0, "ceil_pow2 requires a positive value")
    return 1 << (value - 1).bit_length()


def parent_count(profile: str, k: int, r: int) -> int:
    if profile == "high":
        return ceil_pow2(k + ceil_pow2(r))
    require(profile == "low", "unknown profile")
    return ceil_pow2(ceil_pow2(k) + r)


def make_cell(
    profile: str,
    k: int,
    r: int,
    shard_bytes: int,
    reuse: int,
    tier: str,
) -> dict[str, Any]:
    require(k > 16 and r > 1 and shard_bytes > 0 and reuse > 0,
            "matrix contains an ineligible direct-repair cell")
    parent = parent_count(profile, k, r)
    require(parent <= 256, "matrix cell escapes the GF8 parent")
    if profile == "low" and tier != "balanced-neighbor":
        require(r > k, "non-balanced LOW candidate must have R > K")
    key = "%s:k=%d:r=%d:b=%d:u=%d" % (
        profile, k, r, shard_bytes, reuse)
    seed = int.from_bytes(hashlib.sha256(key.encode()).digest()[:8], "big")
    return {
        "id": "%s-%s-k%d-r%d-b%d-u%d" % (
            profile, tier, k, r, shard_bytes, reuse),
        "profile": profile,
        "tier": tier,
        "K": k,
        "R": r,
        "parent_count": parent,
        "bytes": shard_bytes,
        "loss": 1,
        "batch": 1,
        "reuse": reuse,
        "seed": seed,
        "candidate_expected_path": (
            "transform" if profile == "low" and r == k else "direct"),
    }


def build_matrix() -> dict[str, Any]:
    cells: dict[tuple[str, int, int, int, int], dict[str, Any]] = {}

    def add(
        profile: str,
        shapes: Iterable[tuple[int, int]],
        byte_values: Iterable[int],
        reuse_values: Iterable[int],
        tier: str,
    ) -> None:
        for k, r in shapes:
            for shard_bytes in byte_values:
                for reuse in reuse_values:
                    key = (profile, k, r, shard_bytes, reuse)
                    require(key not in cells,
                            "matrix tiers overlap at %r" % (key,))
                    cells[key] = make_cell(
                        profile, k, r, shard_bytes, reuse, tier)

    add("high", HIGH_CORE, BASE_BYTES, BASE_REUSE, "core")
    add("low", LOW_CORE, BASE_BYTES, BASE_REUSE, "core")
    # Balanced LOW is a negative selector control and therefore receives the
    # complete size/reuse matrix rather than a token spot check.
    add("low", LOW_BALANCED_NEIGHBORS, BASE_BYTES, BASE_REUSE,
        "balanced-neighbor")
    # Shape-neighbor screens are compact: two representative working-set
    # sizes and the middle reuse point.  Full expansion follows only if one
    # of these exposes a boundary regression.
    add("high", HIGH_NEIGHBORS, (4096, 65536), (8,), "shape-neighbor")
    add("low", LOW_NEIGHBORS, (4096, 65536), (8,), "shape-neighbor")
    # Exercise vector/tile tails around every requested byte boundary while
    # holding reuse at the representative amortization point.
    add("high", HIGH_CORE, BYTE_NEIGHBORS, (8,), "byte-neighbor")
    add("low", LOW_CORE, BYTE_NEIGHBORS, (8,), "byte-neighbor")

    high = sorted(
        (cell for cell in cells.values() if cell["profile"] == "high"),
        key=lambda cell: cell["id"])
    low = sorted(
        (cell for cell in cells.values() if cell["profile"] == "low"),
        key=lambda cell: cell["id"])
    document: dict[str, Any] = {
        "schema": MATRIX_SCHEMA,
        "definitions": {
            "field": "gf8",
            "backend": "avx2",
            "loss_count": 1,
            "base_bytes": list(BASE_BYTES),
            "base_reuse": list(BASE_REUSE),
            "byte_neighbors": list(BYTE_NEIGHBORS),
        },
        "high": high,
        "low": low,
        "cell_count": len(high) + len(low),
    }
    document["matrix_sha256"] = sha256_bytes(canonical_bytes(document))
    return document


def validate_matrix(document: Mapping[str, Any]) -> list[dict[str, Any]]:
    require(document.get("schema") == MATRIX_SCHEMA,
            "matrix schema differs")
    retained = dict(document)
    digest = retained.pop("matrix_sha256", None)
    require(digest == sha256_bytes(canonical_bytes(retained)),
            "matrix digest differs")
    high = document.get("high")
    low = document.get("low")
    require(isinstance(high, list) and isinstance(low, list),
            "matrix profiles are malformed")
    cells = high + low
    require(all(isinstance(cell, dict) and cell.get("profile") == "high"
                for cell in high),
            "HIGH matrix contains a cross-labeled cell")
    require(all(isinstance(cell, dict) and cell.get("profile") == "low"
                for cell in low),
            "LOW matrix contains a cross-labeled cell")
    require(document.get("cell_count") == len(cells) and cells,
            "matrix cell count differs")
    require(len({cell.get("id") for cell in cells}) == len(cells),
            "matrix has duplicate cell ids")
    for cell in cells:
        require(isinstance(cell, dict), "matrix cell is malformed")
        expected = make_cell(
            cell["profile"], cell["K"], cell["R"], cell["bytes"],
            cell["reuse"], cell["tier"])
        require(cell == expected, "matrix cell differs from generator")
    return cells


def checked(
    command: Sequence[str],
    *,
    cwd: Path | str | None = None,
    timeout: int = 1800,
    maximum_bytes: int = 64 << 20,
) -> bytes:
    completed = subprocess.run(
        list(command), cwd=cwd, stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env={"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C"},
        timeout=timeout, check=False)
    require(completed.returncode == 0,
            "command failed (%d): %s\n%s" % (
                completed.returncode, " ".join(command),
                completed.stderr.decode(errors="replace")))
    require(len(completed.stdout) <= maximum_bytes and
            len(completed.stderr) <= maximum_bytes,
            "command output exceeds retention bound")
    return completed.stdout


def git_text(source: Path, *arguments: str) -> str:
    return checked(
        ("/usr/bin/git", "-C", str(source), *arguments),
        timeout=120).decode("utf-8", errors="strict").strip()


def source_identity(source: Path, expected_commit: str) -> dict[str, Any]:
    source = source.resolve(strict=True)
    require(COMMIT_PATTERN.fullmatch(expected_commit) is not None,
            "--commit must be a full lowercase SHA")
    require(Path(git_text(source, "rev-parse", "--show-toplevel")) == source,
            "source is not the Git top level")
    require(not git_text(
        source, "for-each-ref", "--format=%(refname)", "refs/replace"),
        "replace refs make source identity ambiguous")
    head = git_text(source, "rev-parse", "HEAD")
    require(head == expected_commit, "source HEAD differs from --commit")
    for flag in ("-v", "-f"):
        records = [record for record in git_text(
            source, "ls-files", flag, "-z").split("\0") if record]
        require(records and all(record.startswith("H ") for record in records),
                "source index uses assume-unchanged, skip-worktree, "
                "fsmonitor-valid, or another non-default flag")
    branch = git_text(source, "branch", "--show-current")
    require(branch not in {"main", "master", ""},
            "authoritative experiment must use a topic branch")
    status = git_text(
        source, "status", "--porcelain=v1", "--untracked-files=normal")
    require(not status, "source tree is not clean:\n" + status)
    return {
        "path": str(source), "head": head,
        "tree": git_text(source, "rev-parse", "HEAD^{tree}"),
        "branch": branch, "status": "clean",
    }


@contextlib.contextmanager
def campaign_lock(path: str):
    descriptor = os.open(path, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        print("waiting for exclusive campaign lock %s" % path, flush=True)
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        print("acquired campaign lock", flush=True)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def prepare_owned_directory(root: Path, identity: Mapping[str, Any]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    marker = root / ".leopard2-general-l1-owner.json"
    if marker.exists():
        require(read_json(marker) == identity,
                "build root belongs to another campaign")
    else:
        require(not any(root.iterdir()),
                "refusing to claim a nonempty unmarked build root")
        atomic_json(marker, identity)


def configure_and_build(
    source: Path,
    build: Path,
    enabled: bool,
    pair_fusion: bool,
    jobs: int,
) -> None:
    require(not pair_fusion or enabled,
            "pair fusion requires the general direct-L1 candidate")
    cmake = "/usr/bin/cmake"
    configure = (
        cmake, "-S", str(source), "-B", str(build),
        "-G", "Unix Makefiles",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_CXX_FLAGS=",
        "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG",
        "-DENABLE_OPENMP=ON",
        "-DLEO2_BACKEND_VARIANT=avx2",
        "-DLEO2_BUILD_ALLK_DIAGNOSTIC=OFF",
        "-DLEO2_BUILD_BENCHMARKS=ON",
        "-DLEO2_BUILD_FUZZERS=OFF",
        "-DLEO2_BUILD_TESTS=OFF",
        "-DLEO2_ENABLE_CUDA=OFF",
        "-DLEOPARD_ENABLE_GF8=ON",
        "-DLEOPARD_ENABLE_GF16=ON",
        # Cache these checks false so no AVX-512 object is compiled or linked.
        "-DLEO2_FLAG_MAVX512F=FALSE",
        "-DLEO2_FLAG_MAVX512BW=FALSE",
        "-DLEO2_FLAG_MAVX512VL=FALSE",
        "-D%s=%s" % (PAIR_FUSION, "ON" if pair_fusion else "OFF"),
        "-D%s=%s" % (EXPERIMENT, "ON" if enabled else "OFF"),
    )
    checked(configure, timeout=300)
    checked((
        cmake, "--build", str(build), "--target", "bench_leopard2",
        "--parallel", str(jobs),
    ), timeout=1800)


def cache_values(build: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in (build / "CMakeCache.txt").read_text(
            encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name_type, value = line.split("=", 1)
        if ":" in name_type:
            result[name_type.split(":", 1)[0]] = value
    return result


def closure_by_source(provenance: Mapping[str, Any]) -> dict[str, Any]:
    source = Path(str(provenance["source_root"]))
    result = {}
    for record in provenance["source_object_compile_closure"]:
        path = Path(record["source"]["path"])
        relative = path.relative_to(source).as_posix()
        require(relative not in result, "duplicate source in compile closure")
        result[relative] = record
    return result


def normalized_arguments(record: Mapping[str, Any], build: Path) -> list[str]:
    result = []
    build_text = str(build.resolve(strict=True))
    for argument in record["compile_entry"]["arguments"]:
        result.append(argument.replace(build_text, "<BUILD>"))
    return result


def compare_build_closures(
    source: Path,
    control_build: Path,
    candidate_build: Path,
    control: Mapping[str, Any],
    candidate: Mapping[str, Any],
    control_general: bool = False,
    control_pair: bool = False,
    candidate_general: bool = True,
    candidate_pair: bool = False,
) -> dict[str, Any]:
    require(not control_pair or control_general,
            "control pair fusion requires general direct L1")
    require(not candidate_pair or candidate_general,
            "candidate pair fusion requires general direct L1")
    require((control_general, control_pair) !=
            (candidate_general, candidate_pair),
            "control and candidate configurations are identical")
    control_cache = cache_values(control_build)
    candidate_cache = cache_values(candidate_build)
    require(control_cache.get(EXPERIMENT) ==
                ("ON" if control_general else "OFF") and
            control_cache.get(PAIR_FUSION) ==
                ("ON" if control_pair else "OFF") and
            candidate_cache.get(EXPERIMENT) ==
                ("ON" if candidate_general else "OFF") and
            candidate_cache.get(PAIR_FUSION) ==
                ("ON" if candidate_pair else "OFF"),
            "control/candidate experiment cache identity differs")
    left = closure_by_source(control)
    right = closure_by_source(candidate)
    require(set(left) == set(right), "control/candidate source closure differs")
    comparison = []
    changed_objects = []
    macro = "-D%s=1" % EXPERIMENT
    pair_macro = "-D%s=1" % PAIR_FUSION
    def configured_macros(general: bool, pair: bool) -> dict[str, str]:
        result = {"leopard2.cpp": macro} if general else {}
        if pair:
            result.update({
                "Leopard2Backend.cpp": pair_macro,
                "Leopard2BackendAVX2.cpp": pair_macro,
            })
        return result

    control_macros = configured_macros(control_general, control_pair)
    candidate_macros = configured_macros(candidate_general, candidate_pair)
    expected_changed = sorted(
        relative for relative in set(control_macros) | set(candidate_macros)
        if control_macros.get(relative) != candidate_macros.get(relative))
    for relative in sorted(left):
        control_record = left[relative]
        candidate_record = right[relative]
        require(control_record["source"]["sha256"] ==
                candidate_record["source"]["sha256"],
                "source bytes differ for " + relative)
        require(control_record["flag_profile"] ==
                candidate_record["flag_profile"],
                "ISA flag profile differs for " + relative)
        control_args = normalized_arguments(control_record, control_build)
        candidate_args = normalized_arguments(candidate_record, candidate_build)
        control_expected = control_macros.get(relative)
        candidate_expected = candidate_macros.get(relative)
        for arguments, expected, role in (
                (control_args, control_expected, "control"),
                (candidate_args, candidate_expected, "candidate")):
            for definition in (macro, pair_macro):
                require(arguments.count(definition) ==
                        (1 if definition == expected else 0),
                        role + " definition scope differs for " + relative)
        control_base = [item for item in control_args
                        if item not in (macro, pair_macro)]
        candidate_base = [item for item in candidate_args
                          if item not in (macro, pair_macro)]
        require(candidate_base == control_base,
                relative + " compile recipe differs beyond approved macros")
        same_object = (
            control_record["object"]["sha256"] ==
            candidate_record["object"]["sha256"] and
            control_record["object"]["size"] ==
            candidate_record["object"]["size"])
        if same_object:
            require(relative not in expected_changed,
                    "candidate macro did not change " + relative + " object")
        else:
            changed_objects.append(relative)
            require(relative in expected_changed,
                    "unrelated object bytes differ for " + relative)
        comparison.append({
            "source": relative,
            "source_sha256": control_record["source"]["sha256"],
            "object_equal": same_object,
            "control_object_sha256": control_record["object"]["sha256"],
            "candidate_object_sha256": candidate_record["object"]["sha256"],
            "flag_profile": control_record["flag_profile"],
        })
    require(changed_objects == expected_changed,
            "object delta differs from the permitted source scope")
    require(not any(name.startswith("Leopard2BackendAVX512") for name in left),
            "pure AVX2 closure contains an AVX-512 source")
    return {
        "schema": "leopard2-general-l1-source-scoped-closure/v2",
        "source_root": str(source),
        "control_macros": sorted(control_macros.items()),
        "candidate_macros": sorted(candidate_macros.items()),
        "allowed_changed_sources": expected_changed,
        "changed_objects": changed_objects,
        "records": comparison,
        "control_archive_sha256": control["archive"]["sha256"],
        "candidate_archive_sha256": candidate["archive"]["sha256"],
        "control_executable_sha256": control["executable"]["sha256"],
        "candidate_executable_sha256": candidate["executable"]["sha256"],
    }


def analyze_disassembly(disassembly: str) -> dict[str, int]:
    instruction = re.compile(r"^\s*[0-9a-f]+:\s+((?:[0-9a-f]{2}\s+)+)")
    evex = 0
    ymm = 0
    for line in disassembly.splitlines():
        match = instruction.match(line)
        if match and match.group(1).split()[0] == "62":
            evex += 1
        if re.search(r"\bymm[0-9]+\b", line):
            ymm += 1
    require(evex == 0, "pure AVX2 executable contains EVEX instructions")
    require(ymm > 0, "pure AVX2 executable contains no visible YMM work")
    return {
        "evex_prefixed_instruction_count": evex,
        "ymm_operand_instruction_count": ymm,
    }


def inspect_pure_avx2(executable: Path) -> dict[str, Any]:
    completed = subprocess.run(
        ("/usr/bin/objdump", "-d", "-M", "intel", str(executable)),
        stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, text=True, timeout=180, check=False)
    require(completed.returncode == 0,
            "objdump failed: " + completed.stderr.strip())
    result: dict[str, Any] = analyze_disassembly(completed.stdout)
    result.update({
        "objdump_version": checked(
            ("/usr/bin/objdump", "--version"), timeout=30).decode(
                errors="replace").splitlines()[0],
    })
    return result


def inspect_object_footprints(
    provenance: Mapping[str, Any], sources: Sequence[str],
) -> dict[str, Any]:
    records = closure_by_source(provenance)
    result: dict[str, Any] = {}
    for relative in sources:
        require(relative in records,
                "footprint source is absent from compile closure: " + relative)
        path = Path(records[relative]["object"]["path"])
        require(path.is_file(), "footprint object is missing: " + str(path))
        completed = subprocess.run(
            ("/usr/bin/size", "-A", "-d", str(path)),
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, timeout=60, check=False)
        require(completed.returncode == 0,
                "size failed for %s: %s" % (
                    relative, completed.stderr.strip()))
        sections: dict[str, int] = {}
        for line in completed.stdout.splitlines():
            columns = line.split()
            if len(columns) >= 2 and columns[0].startswith("."):
                try:
                    sections[columns[0]] = int(columns[1], 10)
                except ValueError:
                    continue
        require(any(name.startswith(".text") for name in sections),
                "size emitted no text section for " + relative)
        result[relative] = {
            "path": str(path.resolve(strict=True)),
            "sha256": sha256_file(path),
            "file_bytes": path.stat().st_size,
            "text_bytes": sum(value for name, value in sections.items()
                              if name.startswith(".text")),
            "rodata_bytes": sum(value for name, value in sections.items()
                                if name.startswith(".rodata")),
            "sections": sections,
        }
    return {
        "size_version": checked(
            ("/usr/bin/size", "--version"), timeout=30).decode(
                errors="replace").splitlines()[0],
        "objects": result,
    }


def inspect_pair_object(
    provenance: Mapping[str, Any], expected: bool,
) -> dict[str, Any]:
    records = closure_by_source(provenance)
    relative = "Leopard2BackendAVX2.cpp"
    require(relative in records,
            "AVX2 backend is absent from the compile closure")
    path = Path(records[relative]["object"]["path"])
    names = checked((
        "/usr/bin/nm", "-C", "--defined-only", str(path)),
        timeout=60).decode(errors="replace")
    pair_names = sorted(line.strip() for line in names.splitlines()
                        if "AVX2FF8LinearCombination2" in line)
    require(bool(pair_names) == expected,
            "pair-fusion symbols do not match the configured option")
    disassembly = checked((
        "/usr/bin/objdump", "-d", "-C", "-M", "intel", str(path)),
        timeout=180).decode(errors="replace")
    pair_lines = [line for line in disassembly.splitlines()
                  if "AVX2FF8LinearCombination2" in line]
    functions: dict[str, list[str]] = {}
    active_name = None
    for line in disassembly.splitlines():
        label = re.match(r"^\s*[0-9a-f]+ <(.*)>:$", line)
        if label:
            active_name = label.group(1)
            functions.setdefault(active_name, [])
        elif active_name is not None:
            functions[active_name].append(line)
    loop_specializations = {}
    loop_pattern = re.compile(
        r"AVX2FF8LinearCombination2Loop<"
        r"(true|false), (true|false), (true|false)>")
    for name, lines in functions.items():
        match = loop_pattern.search(name)
        if not match:
            continue
        add, identity0, identity1 = match.groups()
        key = "add=%s,identity0=%s,identity1=%s" % (
            add, identity0, identity1)
        body = "\n".join(lines)
        loop_specializations[key] = {
            "symbol": name,
            "vbroadcast_count": len(re.findall(
                r"\bvbroadcast(?:i128|f128|ss|sd)\b", body)),
            "vpshufb_count": len(re.findall(r"\bvpshufb\b", body)),
            "body_sha256": sha256_bytes(body.encode()),
        }
    identity_loops = [value for key, value in loop_specializations.items()
                      if "identity0=true,identity1=true" in key]
    if identity_loops:
        require(all(value["vbroadcast_count"] == 0 and
                    value["vpshufb_count"] == 0
                    for value in identity_loops),
                "both-identity pair specialization retained table work")
    if expected:
        require(any("AVX2FF8LinearCombination2(" in name
                    for name in pair_names),
                "pair-fusion callback symbol is missing")
        require(re.search(r"\bymm[0-9]+\b", disassembly) is not None and
                "vpshufb" in disassembly,
                "pair-fusion object has no visible AVX2 table work")
    return {
        "object_path": str(path.resolve(strict=True)),
        "object_sha256": sha256_file(path),
        "pair_symbol_count": len(pair_names),
        "pair_symbols": pair_names,
        "pair_reference_line_count": len(pair_lines),
        "loop_specializations": loop_specializations,
        "both_identity_table_work_elided": (
            bool(identity_loops) and all(
                value["vbroadcast_count"] == 0 and
                value["vpshufb_count"] == 0
                for value in identity_loops)),
        "complete_disassembly_sha256": sha256_bytes(disassembly.encode()),
        "complete_vbroadcast_count": len(re.findall(
            r"\bvbroadcast(?:i128|f128|ss|sd)\b", disassembly)),
        "complete_vpshufb_count": len(re.findall(
            r"\bvpshufb\b", disassembly)),
    }


def freeze_executable(source: Path, destination: Path) -> dict[str, Any]:
    source = source.resolve(strict=True)
    source_digest = sha256_file(source)
    source_size = source.stat().st_size
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        temporary = destination.with_name(
            destination.name + ".tmp.%d" % os.getpid())
        shutil.copyfile(source, temporary)
        os.chmod(temporary, 0o555)
        os.replace(str(temporary), str(destination))
    retained_status = destination.lstat()
    require(stat.S_ISREG(retained_status.st_mode) and
            not destination.is_symlink() and retained_status.st_mode & 0o111,
            "frozen executable is missing or not executable")
    require(retained_status.st_size == source_size and
            sha256_file(destination) == source_digest,
            "frozen executable differs from its build output")
    return {
        "build_path": str(source),
        "frozen_path": str(destination.resolve(strict=True)),
        "sha256": source_digest,
        "size": source_size,
        "mode": oct(retained_status.st_mode & 0o777),
    }


def verify_frozen_artifact(artifact: Mapping[str, Any], role: str) -> None:
    path = Path(str(artifact.get("frozen_path", "")))
    retained_status = path.lstat()
    require(stat.S_ISREG(retained_status.st_mode) and
            not path.is_symlink() and retained_status.st_mode & 0o111,
            role + " frozen executable is missing or not executable")
    require(retained_status.st_size == artifact.get("size") and
            sha256_file(path) == artifact.get("sha256"),
            role + " frozen executable digest differs")


def parse_cpu_list(text: str) -> set[int]:
    result: set[int] = set()
    for part in text.strip().split(","):
        if not part:
            continue
        if "-" in part:
            first, last = (int(value) for value in part.split("-", 1))
            result.update(range(first, last + 1))
        else:
            result.add(int(part))
    return result


def cpu_topology(requested: int | None) -> dict[str, Any]:
    allowed = set(os.sched_getaffinity(0))
    require(bool(allowed), "process affinity is empty")
    candidates = [requested] if requested is not None else sorted(allowed)
    for cpu in candidates:
        require(cpu is not None, "invalid CPU")
        if cpu not in allowed:
            continue
        path = Path(
            "/sys/devices/system/cpu/cpu%d/topology/thread_siblings_list" % cpu)
        siblings = parse_cpu_list(path.read_text(encoding="ascii"))
        if siblings.issubset(allowed):
            reserved = sorted(siblings - {cpu})
            return {
                "benchmark_cpu": cpu,
                "thread_siblings": sorted(siblings),
                "reserved_sibling": reserved[0] if reserved else None,
                "allowed_cpus": sorted(allowed),
            }
    raise CampaignError(
        "no requested physical core has its complete sibling set available")


def cpu_stat(cpu: int) -> dict[str, int]:
    prefix = "cpu%d " % cpu
    with open("/proc/stat", encoding="ascii") as stream:
        for line in stream:
            if line.startswith(prefix):
                values = [int(value) for value in line.split()[1:]]
                while len(values) < 8:
                    values.append(0)
                names = (
                    "user", "nice", "system", "idle", "iowait",
                    "irq", "softirq", "steal",
                )
                return dict(zip(names, values[:8]))
    raise CampaignError("CPU%d is absent from /proc/stat" % cpu)


def nonidle_delta(before: Mapping[str, int], after: Mapping[str, int]) -> int:
    return sum(after[name] - before[name] for name in (
        "user", "nice", "system", "irq", "softirq", "steal"))


def benchmark_environment() -> dict[str, str]:
    return {
        "PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
        "TZ": "UTC", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
        "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE",
    }


def iterations_for(shard_bytes: int) -> int:
    if shard_bytes >= 1048575:
        return 5
    if shard_bytes >= 65535:
        return 7
    if shard_bytes >= 4095:
        return 9
    return 11


def benchmark_command(
    executable: Path,
    cell: Mapping[str, Any],
    iterations: int,
    *,
    leopard2: bool,
) -> list[str]:
    result = [
        str(executable),
        "--k", str(cell["K"]), "--r", str(cell["R"]),
        "--bytes", str(cell["bytes"]), "--loss", str(cell["loss"]),
        "--batch", "1", "--reuse", str(cell["reuse"]),
        "--iterations", str(iterations), "--warmup", "2",
        "--threads", "1", "--seed", str(cell["seed"]),
        "--json", "-",
    ]
    if leopard2:
        result.extend((
            "--profile", cell["profile"], "--field", "gf8",
            "--backend", "avx2", "--skip-legacy", "--retain-samples",
            "--report-decode-path",
        ))
    return result


def run_invocation(
    role: str,
    executable: Path,
    executable_digest: str,
    cell: Mapping[str, Any],
    topology: Mapping[str, Any],
    invocation_dir: Path,
    iterations: int,
    maximum_attempts: int,
    *,
    leopard2: bool,
) -> dict[str, Any]:
    command = benchmark_command(
        executable, cell, iterations, leopard2=leopard2)
    taskset_command = [
        "/usr/bin/taskset", "-c", str(topology["benchmark_cpu"]), *command]
    reserved = topology["reserved_sibling"]
    rejected = []
    timeout = 1200 if cell["bytes"] >= 1048575 else 300
    for attempt in range(maximum_attempts):
        require(sha256_file(executable) == executable_digest,
                role + " frozen executable changed before invocation")
        before = cpu_stat(reserved) if reserved is not None else None
        started_ns = time.monotonic_ns()
        try:
            completed = subprocess.run(
                taskset_command, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env=benchmark_environment(), timeout=timeout, check=False)
        except subprocess.TimeoutExpired as error:
            raise CampaignError(
                "%s benchmark timed out after %ds: %s" % (
                    role, timeout, (error.stderr or b"").decode(
                        errors="replace"))) from error
        ended_ns = time.monotonic_ns()
        after = cpu_stat(reserved) if reserved is not None else None
        require(sha256_file(executable) == executable_digest,
                role + " frozen executable changed after invocation")
        attempt_dir = invocation_dir / ("attempt-%d" % attempt)
        attempt_dir.mkdir(parents=True, exist_ok=True)
        (attempt_dir / (role + ".stdout.json")).write_bytes(completed.stdout)
        (attempt_dir / (role + ".stderr")).write_bytes(completed.stderr)
        require(completed.returncode == 0,
                "%s benchmark failed (%d); see %s" % (
                    role, completed.returncode, attempt_dir))
        sibling_work = (
            nonidle_delta(before, after)
            if before is not None and after is not None else 0)
        if sibling_work:
            rejected.append({
                "attempt": attempt,
                "reason": "reserved sibling performed non-idle work",
                "sibling_nonidle_jiffies": sibling_work,
                "started_ns": started_ns,
                "ended_ns": ended_ns,
            })
            time.sleep(0.05)
            continue
        try:
            document = json.loads(completed.stdout)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise CampaignError(
                "%s benchmark emitted invalid JSON; see %s" % (
                    role, attempt_dir)) from error
        return {
            "role": role,
            "attempt": attempt,
            "rejected_attempts": rejected,
            "started_ns": started_ns,
            "ended_ns": ended_ns,
            "command": taskset_command,
            "stdout_sha256": sha256_bytes(completed.stdout),
            "stderr_sha256": sha256_bytes(completed.stderr),
            "document": document,
        }
    raise CampaignError(
        "%s exhausted attempts because the reserved sibling was busy" % role)


def digest_tuple(document: Mapping[str, Any]) -> tuple[str, str, str]:
    digests = document["workload_digests"]
    require(digests.get("algorithm") == "fnv1a64",
            "unexpected workload digest algorithm")
    return (
        digests["original_data"], digests["transmitted_parity"],
        digests["recovered_originals"],
    )


def missing_tuple(
    document: Mapping[str, Any], cell: Mapping[str, Any], role: str,
) -> tuple[int, ...]:
    values = document["parameters"].get("missing_original_indices")
    require(isinstance(values, list) and
            len(values) == cell["loss"] and
            all(type(value) is int and 0 <= value < cell["K"]
                for value in values) and
            values == sorted(set(values)),
            role + " missing-original identity is malformed")
    return tuple(values)


def validate_leopard2_document(
    role: str,
    document: Mapping[str, Any],
    cell: Mapping[str, Any],
    direct_control: bool = False,
) -> None:
    require(document.get("schema") == "leopard2-benchmark-v3",
            role + " did not emit path-report schema v3")
    parameters = document["parameters"]
    expected_profile = (
        "legacy_high_v1" if cell["profile"] == "high" else "low_v1")
    expected_iterations = int(cell.get(
        "benchmark_iterations", iterations_for(cell["bytes"])))
    require(parameters["K"] == cell["K"] and
            parameters["R"] == cell["R"] and
            parameters["shard_bytes"] == cell["bytes"] and
            parameters["loss_count"] == cell["loss"] and
            parameters["batch"] == 1 and
            parameters["reuse"] == cell["reuse"] and
            parameters["iterations"] == expected_iterations and
            parameters["warmup"] == 2 and
            parameters["thread_count"] == 1 and
            parameters["seed"] == cell["seed"],
            role + " parameters differ from the cell")
    require(parameters["requested_profile"] == expected_profile and
            parameters["requested_field"] == "gf8" and
            parameters["requested_backend"] == "avx2",
            role + " request attestation differs")
    require(parameters["force_generic_decode"] is False and
            parameters["force_specialized_decode"] is False and
            parameters["force_tiled_decode"] is False and
            parameters["force_materialized_decode"] is False and
            parameters["skip_legacy"] is True and
            parameters["retain_samples"] is True and
            parameters["report_decode_path"] is True,
            role + " benchmark-mode attestation differs")
    missing_tuple(document, cell, role)
    resolved = document["resolved"]
    require(resolved["profile"] == expected_profile and
            resolved["field"] == "gf8" and
            resolved["backend"] == "avx2" and
            resolved["thread_count"] == 1 and
            resolved["parent_count"] == cell["parent_count"],
            role + " resolved identity differs")
    require(document["correctness"]["leopard2_round_trip"] is True,
            role + " round trip did not pass")
    path = resolved["selected_decode_path"]
    direct_counts = (
        resolved["decode_direct_term_count"],
        resolved["decode_direct_unit_term_count"],
        resolved["decode_direct_pair_count"],
        resolved["decode_direct_pair_with_unit_count"],
    )
    if role == "control" and not direct_control:
        require(path != "direct",
                "control unexpectedly selected direct repair")
        require(direct_counts == (0, 0, 0, 0),
                "control exposed direct-plan terms")
    elif cell["candidate_expected_path"] == "direct":
        require(path == "direct",
                "candidate did not select direct repair")
        require(resolved["selected_decode_rule"] == "direct",
                "candidate direct path has a non-direct rule")
        require(direct_counts[0] == cell["K"],
                "one-loss direct plan did not retain exactly K terms")
        require(0 <= direct_counts[1] <= direct_counts[0] and
                direct_counts[2] == direct_counts[0] // 2 and
                0 <= direct_counts[3] <= direct_counts[2] and
                direct_counts[3] <= direct_counts[1],
                "direct-plan pair/unit accounting is inconsistent")
    else:
        require(path != "direct",
                "negative selector control unexpectedly selected direct repair")
        require(direct_counts == (0, 0, 0, 0),
                "negative selector control exposed direct-plan terms")


def direct_plan_mix(document: Mapping[str, Any]) -> dict[str, int]:
    resolved = document["resolved"]
    return {
        "term_count": int(resolved["decode_direct_term_count"]),
        "unit_term_count": int(resolved["decode_direct_unit_term_count"]),
        "pair_count": int(resolved["decode_direct_pair_count"]),
        "pair_with_unit_count": int(
            resolved["decode_direct_pair_with_unit_count"]),
    }


def timings(document: Mapping[str, Any], reuse: int) -> dict[str, float]:
    metrics = document["metrics"]
    codec_setup = float(metrics["codec_setup"]["median_us"])
    setup = float(metrics["decode_plan_setup"]["median_us"])
    execution = float(
        metrics["decode_execution"]["median_us_per_batch_call"])
    amortized = float(
        metrics["decode_amortized_at_reuse"][
            "derived_median_us_per_batch_call"])
    require(codec_setup > 0 and setup > 0 and
            execution > 0 and amortized > 0,
            "benchmark emitted a nonpositive timing")
    require(abs(amortized - (execution + setup / reuse)) <=
            max(1e-5, amortized * 2e-6),
            "amortized timing does not match setup/execution")
    encode = float(metrics["encode_execution"]["median_us_per_batch_call"])
    return {
        "codec_setup_us": codec_setup,
        "plan_setup_us": setup,
        "decode_execution_us": execution,
        "decode_amortized_us": amortized,
        "decode_cold_codec_plan_us": codec_setup + setup + execution,
        "encode_execution_us": encode,
    }


def confidence(log_contrasts: Sequence[float], orientation: str) -> dict[str, Any]:
    require(len(log_contrasts) >= 2,
            "at least two ABBA rounds are required")
    critical = {
        2: 12.706204736, 3: 4.302652730, 4: 3.182446305,
        5: 2.776445105, 6: 2.570581836, 7: 2.446911851,
        8: 2.364624252, 9: 2.306004135, 10: 2.262157163,
    }.get(len(log_contrasts), 1.96)
    mean = statistics.mean(log_contrasts)
    half = critical * statistics.stdev(log_contrasts) / math.sqrt(
        len(log_contrasts))
    return {
        "ratio": math.exp(mean),
        "ci95_low": math.exp(mean - half),
        "ci95_high": math.exp(mean + half),
        "round_log_contrasts": list(log_contrasts),
        "orientation": orientation,
    }


def same_source_metrics(
    calls: Sequence[Mapping[str, Any]], rounds: int, reuse: int,
) -> dict[str, Any]:
    require(len(calls) == rounds * len(ABBA), "ABBA call count differs")
    result = {}
    for metric in (
        "codec_setup_us", "plan_setup_us", "decode_execution_us",
        "decode_amortized_us", "decode_cold_codec_plan_us",
        "encode_execution_us",
    ):
        contrasts = []
        for round_index in range(rounds):
            group = calls[round_index * 4:(round_index + 1) * 4]
            control = [
                math.log(call["timings"][metric]) for call in group
                if call["role"] == "control"]
            candidate = [
                math.log(call["timings"][metric]) for call in group
                if call["role"] == "candidate"]
            require(len(control) == 2 and len(candidate) == 2,
                    "ABBA round role count differs")
            contrasts.append(statistics.mean(control) -
                             statistics.mean(candidate))
        result[metric] = confidence(
            contrasts, "control_time_over_candidate_time")
    result["reuse_count"] = reuse
    return result


def validate_main_document(
    document: Mapping[str, Any],
    cell: Mapping[str, Any],
    main_commit: str,
) -> None:
    require(document.get("schema") == "leopard-main-benchmark-v1",
            "exact main benchmark schema differs")
    require(document["build"]["main_source_commit"] == main_commit,
            "exact main benchmark commit differs")
    parameters = document["parameters"]
    require(parameters["K"] == cell["K"] and
            parameters["R"] == cell["R"] and
            parameters["shard_bytes"] == cell["bytes"] and
            parameters["loss_count"] == 1 and
            parameters["batch"] == 1 and
            parameters["reuse"] == cell["reuse"] and
            parameters["iterations"] == iterations_for(cell["bytes"]) and
            parameters["warmup"] == 2 and
            parameters["thread_count"] == 1 and
            parameters["seed"] == cell["seed"],
            "exact main parameters differ")
    missing_tuple(document, cell, "exact main")
    require(document["resolved"]["profile"] == "legacy_high_v1" and
            document["resolved"]["field"] == "gf8" and
            document["resolved"]["thread_count"] == 1,
            "exact main resolved identity differs")
    require(document["correctness"]["round_trip"] is True,
            "exact main round trip did not pass")


def run_exact_main_followup(
    cell: Mapping[str, Any],
    candidate: Path,
    candidate_digest: str,
    main: Path,
    main_digest: str,
    main_commit: str,
    topology: Mapping[str, Any],
    cell_dir: Path,
    rounds: int,
    maximum_attempts: int,
) -> dict[str, Any]:
    require(cell["profile"] == "high" and cell["R"] <= cell["K"] and
            cell["bytes"] % 64 == 0,
            "exact main follow-up requested for an invalid cell")
    calls = []
    for round_index in range(rounds):
        for slot, role in enumerate(("main", "candidate", "candidate", "main")):
            executable = main if role == "main" else candidate
            digest = main_digest if role == "main" else candidate_digest
            invocation = run_invocation(
                role, executable, digest, cell, topology,
                cell_dir / "exact-main" / ("round-%d-slot-%d" % (
                    round_index, slot)), iterations_for(cell["bytes"]),
                maximum_attempts, leopard2=role == "candidate")
            document = invocation["document"]
            if role == "main":
                validate_main_document(document, cell, main_commit)
                value = float(document["metrics"]["decode_including_setup"][
                    "median_us_per_batch_call"])
            else:
                validate_leopard2_document(role, document, cell)
                value = timings(document, cell["reuse"])[
                    "decode_amortized_us"]
            require(value > 0, "exact-main comparison timing is nonpositive")
            invocation["comparison_us"] = value
            calls.append(invocation)
    expected_digest = digest_tuple(calls[0]["document"])
    require(all(digest_tuple(call["document"]) == expected_digest
                for call in calls),
            "exact main and candidate workload digests differ")
    expected_missing = missing_tuple(calls[0]["document"], cell, "exact main")
    require(all(missing_tuple(call["document"], cell, call["role"]) ==
                expected_missing for call in calls),
            "exact main and candidate missing-original identities differ")
    contrasts = []
    for round_index in range(rounds):
        group = calls[round_index * 4:(round_index + 1) * 4]
        main_values = [math.log(call["comparison_us"]) for call in group
                       if call["role"] == "main"]
        candidate_values = [math.log(call["comparison_us"]) for call in group
                            if call["role"] == "candidate"]
        contrasts.append(statistics.mean(main_values) -
                         statistics.mean(candidate_values))
    return {
        "gate": "same-source decode-amortized CI95 low >= promotion threshold",
        "semantics": (
            "main includes setup in every decode; candidate is execution plus "
            "decode-plan setup amortized at the cell reuse count"),
        "metrics": confidence(
            contrasts, "exact_main_time_over_candidate_time"),
        "calls": calls,
    }


def l2_fallback_cells() -> list[dict[str, Any]]:
    result = []
    for profile, shapes in (("high", HIGH_CORE), ("low", LOW_CORE)):
        for k, r in shapes:
            cell = make_cell(profile, k, r, 65, 1, "core")
            cell.update({
                "id": "%s-l2-fallback-k%d-r%d-b65-u1" % (
                    profile, k, r),
                "tier": "l2-fallback",
                "loss": 2,
                "benchmark_iterations": 1,
                "candidate_expected_path": "transform",
            })
            result.append(cell)
    return result


def run_l2_fallback_probes(
    frozen_candidate: Mapping[str, Any],
    topology: Mapping[str, Any],
    run_dir: Path,
    identity_sha256: str,
    maximum_attempts: int,
) -> dict[str, Any]:
    path = run_dir / "preflight" / "l2-fallback.json"
    if path.exists():
        retained = read_json(path)
        require(retained.get("identity_sha256") == identity_sha256 and
                retained.get("status") == "complete",
                "resumed L=2 fallback probe identity differs")
        verify_frozen_artifact(frozen_candidate, "candidate")
        return retained
    calls = []
    for cell in l2_fallback_cells():
        invocation = run_invocation(
            "candidate", Path(frozen_candidate["frozen_path"]),
            frozen_candidate["sha256"], cell, topology,
            run_dir / "invocations" / cell["id"], 1,
            maximum_attempts, leopard2=True)
        validate_leopard2_document("candidate", invocation["document"], cell)
        calls.append({
            "cell": cell,
            "selected_decode_path": invocation["document"]["resolved"][
                "selected_decode_path"],
            "workload_digests": invocation["document"]["workload_digests"],
            "invocation": invocation,
        })
    result = {
        "schema": SCHEMA,
        "status": "complete",
        "identity_sha256": identity_sha256,
        "purpose": "L=2 must fall back from experimental L=1 direct repair",
        "calls": calls,
    }
    atomic_json(path, result)
    return result


def completed_cell(
    path: Path, identity_sha256: str, cell: Mapping[str, Any],
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    result = read_json(path)
    require(result.get("status") == "complete" and
            result.get("identity_sha256") == identity_sha256 and
            result.get("cell") == cell,
            "resumed cell identity differs: " + str(path))
    return result


def run_cell(
    cell: Mapping[str, Any],
    frozen: Mapping[str, Mapping[str, Any]],
    topology: Mapping[str, Any],
    run_dir: Path,
    identity_sha256: str,
    rounds: int,
    maximum_attempts: int,
    promotion_threshold: float,
    main_commit: str | None,
    comparison_mode: str,
) -> dict[str, Any]:
    verify_frozen_artifact(frozen["control"], "control")
    verify_frozen_artifact(frozen["candidate"], "candidate")
    if frozen.get("main") is not None:
        verify_frozen_artifact(frozen["main"], "main")
    result_path = run_dir / "cells" / (cell["id"] + ".json")
    resumed = completed_cell(result_path, identity_sha256, cell)
    if resumed is not None:
        return resumed
    calls = []
    cell_dir = run_dir / "invocations" / cell["id"]
    for round_index in range(rounds):
        for slot, role in enumerate(ABBA):
            artifact = frozen[role]
            invocation = run_invocation(
                role, Path(artifact["frozen_path"]), artifact["sha256"],
                cell, topology,
                cell_dir / "same-source" / ("round-%d-slot-%d" % (
                    round_index, slot)), iterations_for(cell["bytes"]),
                maximum_attempts, leopard2=True)
            validate_leopard2_document(
                role, invocation["document"], cell,
                direct_control=comparison_mode == SIMPLE_VS_PAIR)
            invocation["timings"] = timings(
                invocation["document"], cell["reuse"])
            calls.append(invocation)
    expected_digest = digest_tuple(calls[0]["document"])
    require(all(digest_tuple(call["document"]) == expected_digest
                for call in calls),
            "control/candidate workload digests differ")
    expected_missing = missing_tuple(calls[0]["document"], cell, "control")
    require(all(missing_tuple(call["document"], cell, call["role"]) ==
                expected_missing for call in calls),
            "control/candidate missing-original identities differ")
    expected_plan_mix = {
        role: direct_plan_mix(next(
            call for call in calls if call["role"] == role)["document"])
        for role in ("control", "candidate")
    }
    require(all(direct_plan_mix(call["document"]) ==
                expected_plan_mix[call["role"]] for call in calls),
            "a repeated direct-plan coefficient mix differs")
    if comparison_mode == SIMPLE_VS_PAIR:
        require(expected_plan_mix["control"] ==
                expected_plan_mix["candidate"],
                "pair attribution changed the direct-plan coefficients")
    metrics = same_source_metrics(calls, rounds, cell["reuse"])
    winner = (
        cell["candidate_expected_path"] == "direct" and
        metrics["decode_amortized_us"]["ci95_low"] >= promotion_threshold)
    exact_main = None
    if frozen.get("main") is not None:
        # The exact-main executable is never used as a discovery oracle.
        # Invalid, LOW, and same-source non-winning cells stop here.
        if (winner and cell["profile"] == "high" and
                cell["R"] <= cell["K"] and cell["bytes"] % 64 == 0):
            require(main_commit is not None,
                    "exact main commit was not supplied")
            exact_main = run_exact_main_followup(
                cell,
                Path(frozen["candidate"]["frozen_path"]),
                frozen["candidate"]["sha256"],
                Path(frozen["main"]["frozen_path"]),
                frozen["main"]["sha256"],
                main_commit, topology, cell_dir, rounds, maximum_attempts)
    result = {
        "schema": SCHEMA,
        "status": "complete",
        "identity_sha256": identity_sha256,
        "cell": cell,
        "workload_digests": {
            "algorithm": "fnv1a64",
            "original_data": expected_digest[0],
            "transmitted_parity": expected_digest[1],
            "recovered_originals": expected_digest[2],
        },
        "same_source": {
            "comparison_mode": comparison_mode,
            "control_configuration": (
                "general-on,pair-off" if comparison_mode == SIMPLE_VS_PAIR
                else "general-off,pair-off"),
            "candidate_configuration": (
                "general-on,pair-on" if comparison_mode in
                    (TRANSFORM_VS_PAIR, SIMPLE_VS_PAIR)
                else "general-on,pair-off"),
            "direct_plan_mix": expected_plan_mix,
            "metrics": metrics,
            "memory": {
                "control": next(call for call in calls
                                if call["role"] == "control")["document"][
                                    "memory"],
                "candidate": next(call for call in calls
                                  if call["role"] == "candidate")[
                                      "document"]["memory"],
            },
            "calls": calls,
        },
        "same_source_candidate_win": winner,
        "exact_main": exact_main,
        "completed_at": time.strftime(
            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_json(result_path, result)
    return result


def summarize(
    results: Sequence[Mapping[str, Any]], promotion_threshold: float,
    comparison_mode: str,
) -> dict[str, Any]:
    profiles: dict[str, list[dict[str, Any]]] = {"high": [], "low": []}
    credible_regressions = []
    for result in results:
        cell = result["cell"]
        metrics = result["same_source"]["metrics"]
        compact = {
            "id": cell["id"], "tier": cell["tier"],
            "K": cell["K"], "R": cell["R"],
            "bytes": cell["bytes"], "reuse": cell["reuse"],
            "candidate_expected_path": cell["candidate_expected_path"],
            "execution": metrics["decode_execution_us"],
            "codec_setup": metrics["codec_setup_us"],
            "plan_setup": metrics["plan_setup_us"],
            "amortized": metrics["decode_amortized_us"],
            "cold_codec_plan": metrics["decode_cold_codec_plan_us"],
            "memory": result["same_source"]["memory"],
            "direct_plan_mix": result["same_source"]["direct_plan_mix"],
            "same_source_candidate_win": result[
                "same_source_candidate_win"],
            "exact_main": (
                result["exact_main"]["metrics"]
                if result["exact_main"] is not None else None),
        }
        profiles[cell["profile"]].append(compact)
        if (metrics["decode_execution_us"]["ci95_high"] < 0.98 or
                metrics["decode_amortized_us"]["ci95_high"] < 0.98):
            credible_regressions.append(cell["id"])
    for values in profiles.values():
        values.sort(key=lambda value: value["id"])
    result = {
        "schema": SCHEMA,
        "promotion_threshold": promotion_threshold,
        "regression_threshold": 0.98,
        "ratio_orientation": "control_time_over_candidate_time",
        "comparison_mode": comparison_mode,
        "high": profiles["high"], "low": profiles["low"],
        "same_source_winner_count": sum(
            bool(result["same_source_candidate_win"]) for result in results),
        "credible_regressions": sorted(credible_regressions),
        "exact_main_followup_count": sum(
            result["exact_main"] is not None for result in results),
    }
    coefficient_mix: dict[str, dict[str, int]] = {}
    for cell_result in results:
        cell = cell_result["cell"]
        if cell["candidate_expected_path"] != "direct":
            continue
        key = "%s:k=%d:r=%d" % (
            cell["profile"], cell["K"], cell["R"])
        mix = cell_result["same_source"]["direct_plan_mix"]["candidate"]
        if key in coefficient_mix:
            require(coefficient_mix[key] == mix,
                    "direct coefficient mix changed with bytes or reuse")
        else:
            coefficient_mix[key] = mix
    result["candidate_direct_coefficient_mix_by_shape"] = coefficient_mix

    if comparison_mode == SIMPLE_VS_PAIR:
        direct_results = [cell_result for cell_result in results
                          if cell_result["cell"][
                              "candidate_expected_path"] == "direct"]

        def observed_threshold(
            selected: Sequence[Mapping[str, Any]],
        ) -> dict[str, Any]:
            byte_values = sorted({entry["cell"]["bytes"]
                                  for entry in selected})
            accepted = None
            for byte_count in byte_values:
                tail = [entry for entry in selected
                        if entry["cell"]["bytes"] >= byte_count]
                if tail and all(entry["same_source"]["metrics"][
                        "decode_execution_us"]["ci95_low"] >=
                        promotion_threshold for entry in tail):
                    accepted = byte_count
                    break
            return {
                "observed_execution_threshold_bytes": accepted,
                "observed_byte_values": byte_values,
                "execution_winner_ids": sorted(
                    entry["cell"]["id"] for entry in selected
                    if entry["same_source"]["metrics"][
                        "decode_execution_us"]["ci95_low"] >=
                        promotion_threshold),
                "execution_regression_ids": sorted(
                    entry["cell"]["id"] for entry in selected
                    if entry["same_source"]["metrics"][
                        "decode_execution_us"]["ci95_high"] < 0.98),
            }

        by_shape_reuse: dict[str, Any] = {}
        keys = sorted({(
            entry["cell"]["profile"], entry["cell"]["K"],
            entry["cell"]["R"], entry["cell"]["reuse"])
            for entry in direct_results})
        for profile, k, r, reuse in keys:
            selected = [entry for entry in direct_results
                        if (entry["cell"]["profile"], entry["cell"]["K"],
                            entry["cell"]["R"], entry["cell"]["reuse"]) ==
                            (profile, k, r, reuse)]
            by_shape_reuse["%s:k=%d:r=%d:u=%d" % (
                profile, k, r, reuse)] = observed_threshold(selected)
        result["pair_fusion_threshold_screen"] = {
            "scope": (
                "observed cells only; a production boundary still requires "
                "neighbor confirmation"),
            "gate": "decode execution CI95 low >= promotion threshold",
            "all_direct_cells": observed_threshold(direct_results),
            "by_shape_and_reuse": by_shape_reuse,
        }
    return result


def machine_identity(topology: Mapping[str, Any]) -> dict[str, Any]:
    model = "unknown"
    with open("/proc/cpuinfo", encoding="ascii", errors="replace") as stream:
        for line in stream:
            if line.startswith("model name") and ":" in line:
                model = line.split(":", 1)[1].strip()
                break
    return {
        "platform": platform.platform(),
        "uname": list(platform.uname()),
        "python": sys.version,
        "cpu_model": model,
        "topology": topology,
    }


def expect_campaign_error(callable_value, label: str) -> None:
    try:
        callable_value()
    except CampaignError:
        return
    raise CampaignError("adversarial self-test was accepted: " + label)


def synthetic_document(
    cell: Mapping[str, Any], path: str,
) -> dict[str, Any]:
    profile = "legacy_high_v1" if cell["profile"] == "high" else "low_v1"
    return {
        "schema": "leopard2-benchmark-v3",
        "parameters": {
            "K": cell["K"], "R": cell["R"],
            "shard_bytes": cell["bytes"],
            "loss_count": cell["loss"],
            "missing_original_indices": list(range(cell["loss"])),
            "batch": 1, "reuse": cell["reuse"],
            "iterations": int(cell.get(
                "benchmark_iterations", iterations_for(cell["bytes"]))),
            "warmup": 2,
            "thread_count": 1, "seed": cell["seed"],
            "requested_profile": profile,
            "requested_field": "gf8", "requested_backend": "avx2",
            "force_generic_decode": False,
            "force_specialized_decode": False,
            "force_tiled_decode": False,
            "force_materialized_decode": False,
            "skip_legacy": True, "retain_samples": True,
            "report_decode_path": True,
        },
        "resolved": {
            "profile": profile, "field": "gf8", "backend": "avx2",
            "thread_count": 1, "parent_count": cell["parent_count"],
            "selected_decode_path": path,
            "selected_decode_rule": "direct" if path == "direct" else "auto",
            "decode_direct_term_count": (
                cell["K"] if path == "direct" else 0),
            "decode_direct_unit_term_count": (
                1 if path == "direct" else 0),
            "decode_direct_pair_count": (
                cell["K"] // 2 if path == "direct" else 0),
            "decode_direct_pair_with_unit_count": (
                1 if path == "direct" and cell["K"] >= 2 else 0),
        },
        "correctness": {"leopard2_round_trip": True},
        "workload_digests": {
            "algorithm": "fnv1a64", "original_data": "01",
            "transmitted_parity": "02", "recovered_originals": "03",
        },
    }


def synthetic_build_closure(
    source: Path, control_build: Path, candidate_build: Path,
    control_general: bool = False,
    control_pair: bool = False,
    candidate_general: bool = True,
    candidate_pair: bool = False,
) -> tuple[dict[str, Any], dict[str, Any]]:
    sources = (
        ("leopard2.cpp", "portable-core"),
        ("Leopard2Backend.cpp", "portable-core"),
        ("Leopard2BackendAVX2.cpp", "avx2-no-avx512"),
        ("bench/leopard2/benchmark.cpp", "portable-core"),
    )

    def provenance(
        build: Path, label: str, general: bool, pair: bool,
    ) -> dict[str, Any]:
        records = []
        for index, (relative, flag_profile) in enumerate(sources):
            source_path = source / relative
            object_path = build / ("object-%d.o" % index)
            arguments = [
                "/usr/bin/c++", "-O3", "-c", str(source_path),
                "-o", str(object_path),
            ]
            if general and relative == "leopard2.cpp":
                arguments.insert(2, "-D%s=1" % EXPERIMENT)
            if pair and relative in {
                    "Leopard2Backend.cpp", "Leopard2BackendAVX2.cpp"}:
                arguments.insert(2, "-D%s=1" % PAIR_FUSION)
            relevant = (
                "general=%d" % general if relative == "leopard2.cpp" else
                "pair=%d" % pair if relative in {
                    "Leopard2Backend.cpp", "Leopard2BackendAVX2.cpp"} else
                "same")
            object_hash = "%s:%s" % (relative, relevant)
            records.append({
                "role": (
                    "benchmark" if relative.startswith("bench/")
                    else "archive"),
                "source": {
                    "path": str(source_path),
                    "sha256": "source-%d" % index,
                },
                "object": {
                    "path": str(object_path), "sha256": object_hash,
                    "size": 100 + index,
                },
                "compile_entry": {"arguments": arguments},
                "flag_profile": flag_profile,
            })
        return {
            "source_root": str(source),
            "source_object_compile_closure": records,
            "archive": {"sha256": "archive-%s" % label},
            "executable": {"sha256": "executable-%s" % label},
        }

    (control_build / "CMakeCache.txt").write_text(
        "%s:BOOL=%s\n%s:BOOL=%s\n" % (
            EXPERIMENT, "ON" if control_general else "OFF",
            PAIR_FUSION, "ON" if control_pair else "OFF"),
        encoding="utf-8")
    (candidate_build / "CMakeCache.txt").write_text(
        "%s:BOOL=ON\n%s:BOOL=%s\n" % (
            EXPERIMENT, PAIR_FUSION, "ON" if candidate_pair else "OFF"),
        encoding="utf-8")
    if not candidate_general:
        (candidate_build / "CMakeCache.txt").write_text(
            "%s:BOOL=OFF\n%s:BOOL=%s\n" % (
                EXPERIMENT, PAIR_FUSION,
                "ON" if candidate_pair else "OFF"), encoding="utf-8")
    return (
        provenance(control_build, "control", control_general, control_pair),
        provenance(candidate_build, "candidate",
                   candidate_general, candidate_pair),
    )


def adversarial_self_tests() -> dict[str, str]:
    with tempfile.TemporaryDirectory(
            prefix="leopard2-general-l1-self-test-") as temporary:
        root = Path(temporary)
        source = root / "source"
        control_build = root / "control"
        candidate_build = root / "candidate"
        source.mkdir()
        control_build.mkdir()
        candidate_build.mkdir()
        control, candidate = synthetic_build_closure(
            source, control_build, candidate_build)
        compare_build_closures(
            source, control_build, candidate_build, control, candidate)

        leakage = copy.deepcopy(candidate)
        leakage["source_object_compile_closure"][1]["compile_entry"][
            "arguments"].insert(2, "-D%s=1" % EXPERIMENT)
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, control, leakage),
            "option source-scope leakage")

        stale = copy.deepcopy(candidate)
        stale["source_object_compile_closure"][0]["object"]["sha256"] = \
            control["source_object_compile_closure"][0]["object"]["sha256"]
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, control, stale),
            "stale candidate leopard2.cpp object")

        mixed = copy.deepcopy(candidate)
        mixed["source_object_compile_closure"][1]["object"]["sha256"] = \
            "mixed-off-on-object"
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, control, mixed),
            "mixed unrelated OFF/ON object")

        pair_control, pair_candidate = synthetic_build_closure(
            source, control_build, candidate_build, candidate_pair=True)
        compare_build_closures(source, control_build, candidate_build,
            pair_control, pair_candidate, candidate_pair=True)

        missing_pair_macro = copy.deepcopy(pair_candidate)
        missing_pair_macro["source_object_compile_closure"][2][
            "compile_entry"]["arguments"] = [argument for argument in
                missing_pair_macro["source_object_compile_closure"][2][
                    "compile_entry"]["arguments"]
                if argument != "-D%s=1" % PAIR_FUSION]
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, pair_control,
            missing_pair_macro, candidate_pair=True),
            "missing AVX2 pair-fusion macro")

        leaked_pair_macro = copy.deepcopy(pair_candidate)
        leaked_pair_macro["source_object_compile_closure"][3][
            "compile_entry"]["arguments"].insert(
                2, "-D%s=1" % PAIR_FUSION)
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, pair_control,
            leaked_pair_macro, candidate_pair=True),
            "pair-fusion macro leakage")

        stale_pair_object = copy.deepcopy(pair_candidate)
        stale_pair_object["source_object_compile_closure"][1]["object"][
            "sha256"] = pair_control["source_object_compile_closure"][1][
                "object"]["sha256"]
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build, pair_control,
            stale_pair_object, candidate_pair=True),
            "stale pair-fusion qualification object")

        direct_control, fused_candidate = synthetic_build_closure(
            source, control_build, candidate_build,
            control_general=True, candidate_pair=True)
        compare_build_closures(
            source, control_build, candidate_build,
            direct_control, fused_candidate,
            control_general=True, candidate_pair=True)
        changed_direct_object = copy.deepcopy(fused_candidate)
        changed_direct_object["source_object_compile_closure"][0]["object"][
            "sha256"] = "unexpected-direct-object-drift"
        expect_campaign_error(lambda: compare_build_closures(
            source, control_build, candidate_build,
            direct_control, changed_direct_object,
            control_general=True, candidate_pair=True),
            "pair attribution changed direct executor object")

        analyze_disassembly(
            "  0: c5 fd ef c0          vpxor ymm0,ymm0,ymm0\n")
        expect_campaign_error(lambda: analyze_disassembly(
            "  0: 62 f1 7d 28 ef c0    vpxord ymm0,ymm0,ymm0\n"),
            "EVEX injection")
        expect_campaign_error(lambda: analyze_disassembly(
            "  0: 66 0f ef c0          pxor xmm0,xmm0\n"),
            "missing AVX2/YMM instructions")

        cell = make_cell("high", 240, 16, 64, 1, "core")
        valid = synthetic_document(cell, "direct")
        validate_leopard2_document("candidate", valid, cell)
        for section, name, bad_value in (
            ("resolved", "profile", "low_v1"),
            ("resolved", "field", "gf16"),
            ("resolved", "backend", "avx512"),
            ("parameters", "requested_profile", "low_v1"),
            ("parameters", "requested_field", "gf16"),
            ("parameters", "requested_backend", "scalar"),
        ):
            malformed = copy.deepcopy(valid)
            malformed[section][name] = bad_value
            expect_campaign_error(lambda malformed=malformed:
                validate_leopard2_document("candidate", malformed, cell),
                "wrong %s.%s" % (section, name))

        l2_cell = l2_fallback_cells()[0]
        l2_document = synthetic_document(l2_cell, "tiled")
        validate_leopard2_document("candidate", l2_document, l2_cell)
        l2_direct = copy.deepcopy(l2_document)
        l2_direct["resolved"]["selected_decode_path"] = "direct"
        l2_direct["resolved"]["selected_decode_rule"] = "direct"
        expect_campaign_error(lambda: validate_leopard2_document(
            "candidate", l2_direct, l2_cell), "L=2 direct-path leak")

        cross_labeled = copy.deepcopy(build_matrix())
        cross_labeled["low"].append(cross_labeled["high"].pop())
        retained = dict(cross_labeled)
        retained.pop("matrix_sha256")
        cross_labeled["matrix_sha256"] = sha256_bytes(
            canonical_bytes(retained))
        expect_campaign_error(lambda: validate_matrix(cross_labeled),
                              "HIGH/LOW cross-labeling")

        executable = root / "frozen-benchmark"
        executable.write_bytes(b"first")
        executable.chmod(0o555)
        artifact = {
            "frozen_path": str(executable), "sha256": sha256_file(executable),
            "size": executable.stat().st_size,
        }
        verify_frozen_artifact(artifact, "self-test")
        executable.chmod(0o755)
        executable.write_bytes(b"second")
        executable.chmod(0o555)
        expect_campaign_error(
            lambda: verify_frozen_artifact(artifact, "resumed"),
            "resumed artifact hash mismatch")

    return {
        "option_source_scope_leakage": "rejected",
        "pair_fusion_source_scope": "accepted_exact_rejected_drift",
        "pair_attribution_source_scope": "accepted_exact_rejected_drift",
        "stale_or_mixed_objects": "rejected",
        "evex_injection": "rejected",
        "wrong_profile_field_backend": "rejected",
        "l2_direct_path": "rejected",
        "high_low_cross_label": "rejected",
        "resumed_artifact_hash_mismatch": "rejected",
    }


def run_campaign(args: argparse.Namespace) -> None:
    require(2 <= args.rounds <= 20, "--rounds must be in 2..20")
    require(1 <= args.maximum_attempts <= 10,
            "--maximum-attempts must be in 1..10")
    require(1.0 < args.promotion_threshold < 2.0,
            "--promotion-threshold must be in (1,2)")
    comparison_mode = (
        SIMPLE_VS_PAIR if args.pair_attribution else
        TRANSFORM_VS_PAIR if args.pair_fusion else
        TRANSFORM_VS_SIMPLE)
    main_arguments = (
        args.main_benchmark, args.main_commit, args.main_sha256)
    require(all(main_arguments) or not any(main_arguments),
            "exact main arguments must be supplied together")
    require(comparison_mode != SIMPLE_VS_PAIR or not any(main_arguments),
            "pair attribution does not invoke the exact-main follow-up")
    jobs = args.jobs if args.jobs is not None else min(
        128, len(os.sched_getaffinity(0)))
    require(1 <= jobs <= 128, "--jobs must be in 1..128")
    source = Path(args.source).resolve(strict=True)
    build_root = Path(args.build_root).resolve(strict=False)
    run_dir = Path(args.run_dir).resolve(strict=False)
    require(not build_root.is_relative_to(source) and
            not run_dir.is_relative_to(source),
            "authoritative build and run roots must be outside the source tree")
    require(not build_root.is_relative_to(run_dir) and
            not run_dir.is_relative_to(build_root),
            "authoritative build and result roots must be disjoint")
    matrix = build_matrix()
    cells = validate_matrix(matrix)
    selected_tiers = set(args.tiers.split(","))
    known_tiers = {
        "core", "balanced-neighbor", "shape-neighbor", "byte-neighbor"}
    require(selected_tiers and selected_tiers.issubset(known_tiers),
            "--tiers contains an unknown or empty tier")
    cells = [cell for cell in cells if cell["tier"] in selected_tiers]
    if args.cell_regex:
        try:
            cell_pattern = re.compile(args.cell_regex)
        except re.error as error:
            raise CampaignError("--cell-regex is invalid: %s" % error)
        cells = [cell for cell in cells if cell_pattern.search(cell["id"])]
    require(cells, "selected matrix is empty")
    source_record = source_identity(source, args.commit)
    topology = cpu_topology(args.cpu)

    owner = {
        "schema": SCHEMA, "source_head": args.commit,
        "runner": str(Path(__file__).resolve(strict=True)),
        "comparison_mode": comparison_mode,
    }
    with campaign_lock(args.lock):
        prepare_owned_directory(build_root, {**owner, "role": "build"})
        prepare_owned_directory(run_dir, {**owner, "role": "results"})
        control_build = build_root / "control"
        candidate_build = build_root / "candidate"
        control_general = comparison_mode == SIMPLE_VS_PAIR
        candidate_pair = comparison_mode in (
            TRANSFORM_VS_PAIR, SIMPLE_VS_PAIR)
        configure_and_build(
            source, control_build, control_general, False, jobs)
        configure_and_build(
            source, candidate_build, True, candidate_pair, jobs)
        # The build and timing phases share one exclusive lease so another
        # worker cannot replace either binary between provenance and timing.
        control_provenance = candidate_build_provenance(
            control_build, source, control_build / "bench_leopard2",
            "bench_leopard2")
        candidate_provenance = candidate_build_provenance(
            candidate_build, source, candidate_build / "bench_leopard2",
            "bench_leopard2")
        closure = compare_build_closures(
            source, control_build, candidate_build,
            control_provenance, candidate_provenance,
            control_general=control_general,
            candidate_pair=candidate_pair)
        require(closure["control_executable_sha256"] !=
                closure["candidate_executable_sha256"],
                "control and candidate executable bytes are unexpectedly equal")
        isa = {
            "control": inspect_pure_avx2(control_build / "bench_leopard2"),
            "candidate": inspect_pure_avx2(
                candidate_build / "bench_leopard2"),
        }
        footprint_sources = sorted(
            set(dict(closure["control_macros"])) |
            set(dict(closure["candidate_macros"])))
        object_footprints = {
            "control": inspect_object_footprints(
                control_provenance, footprint_sources),
            "candidate": inspect_object_footprints(
                candidate_provenance, footprint_sources),
        }
        pair_object_inspection = {
            "control": inspect_pair_object(
                control_provenance, expected=False),
            "candidate": inspect_pair_object(
                candidate_provenance, expected=candidate_pair),
        }
        footprint_deltas = {}
        for relative in footprint_sources:
            control_object = object_footprints["control"]["objects"][relative]
            candidate_object = object_footprints["candidate"]["objects"][relative]
            footprint_deltas[relative] = {
                "file_bytes": candidate_object["file_bytes"] -
                    control_object["file_bytes"],
                "text_bytes": candidate_object["text_bytes"] -
                    control_object["text_bytes"],
                "rodata_bytes": candidate_object["rodata_bytes"] -
                    control_object["rodata_bytes"],
            }
        frozen = {
            "control": freeze_executable(
                control_build / "bench_leopard2",
                run_dir / "artifacts" / "control-bench_leopard2"),
            "candidate": freeze_executable(
                candidate_build / "bench_leopard2",
                run_dir / "artifacts" / "candidate-bench_leopard2"),
        }
        main_commit = None
        if args.main_benchmark:
            require(args.main_commit is not None and
                    COMMIT_PATTERN.fullmatch(args.main_commit) is not None,
                    "--main-benchmark requires a full --main-commit")
            require(args.main_sha256 is not None and
                    SHA256_PATTERN.fullmatch(args.main_sha256) is not None,
                    "--main-benchmark requires an exact --main-sha256")
            main_commit = args.main_commit
            frozen["main"] = freeze_executable(
                Path(args.main_benchmark),
                run_dir / "artifacts" / "exact-main-benchmark")
            require(frozen["main"]["sha256"] == args.main_sha256,
                    "exact main executable differs from --main-sha256")

        identity = {
            "schema": SCHEMA,
            "source": source_record,
            "runner": {
                "path": str(Path(__file__).resolve(strict=True)),
                "sha256": sha256_file(__file__),
            },
            "matrix_sha256": matrix["matrix_sha256"],
            "selected_tiers": sorted(selected_tiers),
            "cell_regex": args.cell_regex,
            "selected_cell_ids": [cell["id"] for cell in cells],
            "rounds": args.rounds,
            "maximum_attempts": args.maximum_attempts,
            "promotion_threshold": args.promotion_threshold,
            "comparison_mode": comparison_mode,
            "build_closure": closure,
            "pure_avx2": isa,
            "object_footprints": object_footprints,
            "object_footprint_candidate_minus_control": footprint_deltas,
            "pair_object_inspection": pair_object_inspection,
            "frozen_executables": frozen,
            "main_commit": main_commit,
            "l2_fallback_probe_ids": [
                cell["id"] for cell in l2_fallback_cells()],
            "machine": machine_identity(topology),
            "lock": os.path.realpath(args.lock),
        }
        identity_sha256 = sha256_bytes(canonical_bytes(identity))
        identity_path = run_dir / "identity.json"
        if identity_path.exists():
            require(read_json(identity_path) == identity,
                    "run directory identity differs")
        else:
            atomic_json(identity_path, identity)
            atomic_json(run_dir / "manifest.json", matrix)
            atomic_json(run_dir / "provenance" / "control.json",
                        control_provenance)
            atomic_json(run_dir / "provenance" / "candidate.json",
                        candidate_provenance)
        (run_dir / "cells").mkdir(exist_ok=True)
        run_l2_fallback_probes(
            frozen["candidate"], topology, run_dir, identity_sha256,
            args.maximum_attempts)
        results = []
        for index, cell in enumerate(cells, 1):
            print("[%d/%d] %s" % (index, len(cells), cell["id"]),
                  flush=True)
            results.append(run_cell(
                cell, frozen, topology, run_dir, identity_sha256,
                args.rounds, args.maximum_attempts,
                args.promotion_threshold, main_commit, comparison_mode))
        # Re-attest immutable source and frozen binary bytes after all timing.
        require(source_identity(source, args.commit) == source_record,
                "source identity changed during campaign")
        for role, artifact in frozen.items():
            require(sha256_file(artifact["frozen_path"]) == artifact["sha256"],
                    role + " frozen executable changed during campaign")
        require(inspect_object_footprints(
                    control_provenance, footprint_sources) ==
                object_footprints["control"] and
                inspect_object_footprints(
                    candidate_provenance, footprint_sources) ==
                object_footprints["candidate"],
                "source object footprint or bytes changed during campaign")
        require(inspect_pair_object(
                    control_provenance, expected=False) ==
                pair_object_inspection["control"] and
                inspect_pair_object(
                    candidate_provenance, expected=candidate_pair) ==
                pair_object_inspection["candidate"],
                "pair object or assembly changed during campaign")
        summary = summarize(
            results, args.promotion_threshold, comparison_mode)
        summary["identity_sha256"] = identity_sha256
        atomic_json(run_dir / "summary.json", summary)
        print(json.dumps({
            "status": "complete", "run_dir": str(run_dir),
            "cells": len(results),
            "same_source_winners": summary["same_source_winner_count"],
            "credible_regressions": len(summary["credible_regressions"]),
            "exact_main_followups": summary["exact_main_followup_count"],
        }, sort_keys=True))


def self_test() -> None:
    matrix = build_matrix()
    cells = validate_matrix(matrix)
    high = [cell for cell in cells if cell["profile"] == "high"]
    low = [cell for cell in cells if cell["profile"] == "low"]
    require({(cell["K"], cell["R"]) for cell in high
             if cell["tier"] == "core"} == set(HIGH_CORE),
            "HIGH core shapes differ")
    require({(cell["K"], cell["R"]) for cell in low
             if cell["tier"] == "core"} == set(LOW_CORE),
            "LOW core shapes differ")
    require(all(cell["candidate_expected_path"] == "transform"
                for cell in low if cell["tier"] == "balanced-neighbor"),
            "balanced LOW is not a negative selector control")
    require(all(cell["parent_count"] <= 256 for cell in cells),
            "matrix escaped GF8")
    expected_core = (len(HIGH_CORE) + len(LOW_CORE)) * len(BASE_BYTES) * \
        len(BASE_REUSE)
    require(sum(cell["tier"] == "core" for cell in cells) == expected_core,
            "core matrix cardinality differs")
    adversarial = adversarial_self_tests()
    print(json.dumps({
        "self_test": "passed", "matrix_sha256": matrix["matrix_sha256"],
        "cells": len(cells), "high_cells": len(high),
        "low_cells": len(low), "core_cells": expected_core,
        "l2_fallback_probes": len(l2_fallback_cells()),
        "adversarial": adversarial,
    }, sort_keys=True))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    matrix = subparsers.add_parser("make-matrix")
    matrix.add_argument("--output", required=True)
    subparsers.add_parser("self-test")
    run = subparsers.add_parser("run")
    run.add_argument("--source", required=True)
    run.add_argument("--commit", required=True)
    run.add_argument("--build-root", required=True)
    run.add_argument("--run-dir", required=True)
    run.add_argument("--main-benchmark")
    run.add_argument("--main-commit")
    run.add_argument("--main-sha256")
    run.add_argument("--cpu", type=int)
    run.add_argument("--jobs", type=int)
    run.add_argument("--rounds", type=int, default=5)
    run.add_argument("--maximum-attempts", type=int, default=3)
    run.add_argument("--promotion-threshold", type=float, default=1.05)
    comparison = run.add_mutually_exclusive_group()
    comparison.add_argument(
        "--pair-fusion", action="store_true",
        help="enable the nested two-source/one-output AVX2 candidate")
    comparison.add_argument(
        "--pair-attribution", action="store_true",
        help="compare general direct repair with pair fusion OFF versus ON")
    run.add_argument(
        "--tiers",
        default="core,balanced-neighbor,shape-neighbor,byte-neighbor",
        help="comma-separated subset of core, balanced-neighbor, "
             "shape-neighbor, byte-neighbor")
    run.add_argument(
        "--cell-regex",
        help="optional regular expression selecting cell IDs after tiers")
    run.add_argument("--lock", default=DEFAULT_LOCK)
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    if args.command == "make-matrix":
        atomic_json(args.output, build_matrix())
    elif args.command == "self-test":
        self_test()
    else:
        run_campaign(args)


if __name__ == "__main__":
    try:
        main()
    except (CampaignError, OSError, subprocess.SubprocessError) as error:
        print("error: %s" % error, file=sys.stderr)
        raise SystemExit(1)
