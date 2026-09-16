#!/usr/bin/env python3
"""Clock-free native parity checks; run inside a serialized 256-MiB scope."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

BASELINE = "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198"
CELLS = (
    ("copy", 1, 1, 4096), ("small", 16, 8, 64),
    ("gf8-high", 240, 16, 65536), ("gf8-balanced", 128, 128, 65536),
    ("gf16-inflation", 200, 50, 65536),
    ("gf16-gfni-region", 1000, 200, 65536),
    ("gf16-explicit-avx2", 1000, 200, 65536),
    ("gf16-large", 4096, 512, 4096),
)
ARTIFACTS = ("main", "current", "main-guard", "current-guard", "guard-control",
             "main.a", "current.a")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
        os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
    return digest.hexdigest()


def validate(record, index, implementation, commit):
    name, k, r, size = CELLS[index]
    baseline = implementation == "main"
    exact = {
        "schema": "leopard-native-release-encode-check/v1",
        "codec_commit": commit,
        "implementation": "leopard1-native" if baseline else "leopard2",
        "cell": index, "id": name, "k": k, "r": r, "bytes": size,
        "requested_backend": "native" if baseline else "avx2" if index == 6 else "auto",
        "threads": 1, "input_bytes": k * size, "parity_bytes": r * size,
        "separate_output_bytes": 0 if baseline else r * size,
        "output_layout": "parity_is_first_r_work_rows" if baseline else
                         "separate_parity_and_scratch",
        "input_unchanged": True, "repeated_encode_equal": True,
        "public_encode_calls": 2,
    }
    require(type(record) is dict and set(record) == set(exact) |
            {"workspace_bytes", "context_backend", "input_hash", "parity_hash"},
            "unexpected record fields")
    for key, expected in exact.items():
        require(type(record[key]) is type(expected) and record[key] == expected,
                "record mismatch: " + key)
    workspace = record["workspace_bytes"]
    require(type(workspace) is int and 0 <= workspace < 256 * 1024 * 1024,
            "invalid workspace")
    backend = record["context_backend"]
    require(type(backend) is int and
            (backend == -1 if baseline else backend == 3 if index == 6 else 1 <= backend <= 6),
            "backend mismatch")
    if baseline:
        count = r if k == 1 else 1 if r == 1 else 2 * (1 << (r - 1).bit_length())
        require(workspace == count * size, "native work geometry mismatch")
    for key in ("input_hash", "parity_hash"):
        require(type(record[key]) is str and re.fullmatch("[0-9a-f]{16}", record[key]),
                "invalid digest")


def equal_files(left, right, size):
    require(left.stat().st_size == right.stat().st_size == size, "parity size mismatch")
    with left.open("rb") as a, right.open("rb") as b:
        while True:
            block = a.read(65536)
            require(block == b.read(65536), "full parity bytes differ")
            if not block:
                break
        for stream in (a, b):
            os.posix_fadvise(stream.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)


def run(bundle, output, candidate_commit):
    require(re.fullmatch("[0-9a-f]{40}", candidate_commit), "invalid candidate commit")
    manifest = bundle / "SHA256SUMS"
    pins = {}
    for line in manifest.read_text().splitlines():
        digest, name = line.split()
        require(re.fullmatch("[0-9a-f]{64}", digest) and name not in pins,
                "invalid or duplicate artifact hash")
        pins[name] = digest
    require(set(pins) == set(ARTIFACTS), "artifact inventory mismatch")
    manifest_digest = sha256(manifest)
    identities = {}

    def verify():
        require(sha256(manifest) == manifest_digest, "manifest changed")
        for name in ARTIFACTS:
            path = bundle / name
            info = path.lstat()
            require(path.is_file() and not path.is_symlink() and not info.st_mode & 0o222,
                    "artifact is not a frozen regular file: " + name)
            require(sha256(path) == pins[name], "artifact bytes changed: " + name)
            identity = (info.st_dev, info.st_ino, info.st_size, info.st_mode,
                        info.st_mtime_ns, info.st_ctime_ns)
            require(name not in identities or identities[name] == identity,
                    "artifact metadata changed: " + name)
            identities[name] = identity

    verify()
    output.mkdir(mode=0o700)  # No resume or overwrite.
    env = {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C",
           "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE", "OMP_THREAD_LIMIT": "1"}
    records = []
    try:
        def invoke(name, arguments, label, expected_code):
            verify()
            with (output / (label + ".stdout")).open("xb") as stdout, \
                    (output / (label + ".stderr")).open("xb") as stderr:
                child = subprocess.run([str(bundle / name), *arguments], env=env,
                                       stdout=stdout, stderr=stderr, timeout=60)
            verify()
            require(child.returncode == expected_code, "unexpected exit: " + label)
            return output / (label + ".stdout"), output / (label + ".stderr")

        for clock in ("clock_gettime", "gettimeofday", "clock"):
            stdout, stderr = invoke("guard-control", [clock], "guard-" + clock, 86)
            require(stdout.stat().st_size == 0 and stderr.read_text() ==
                    "unexpected timing clock in clock-free probe\n", "guard did not fire")
        for name in ARTIFACTS[:4]:
            for index, arguments in enumerate(([], ["--measure", "0"], ["--check", "8"],
                                              ["--check", "00"], ["--check", "-1"])):
                stdout, stderr = invoke(name, arguments, f"bad-{name}-{index}", 1)
                require(stdout.stat().st_size == 0 and stderr.stat().st_size > 0,
                        "invalid CLI accepted")
        for index, (_, _, r, size) in enumerate(CELLS):
            cell_records = {}
            for name in ARTIFACTS[:4]:
                parity = output / f"cell-{index}-{name}.parity"
                stdout, stderr = invoke(name, ["--check", str(index), str(parity)],
                                        f"cell-{index}-{name}", 0)
                require(stderr.stat().st_size == 0, "unexpected diagnostic")
                record = json.loads(stdout.read_text())
                implementation = name.split("-")[0]
                validate(record, index, implementation,
                         BASELINE if implementation == "main" else candidate_commit)
                cell_records[name] = record
                require(parity.stat().st_size == r * size, "incorrect parity file size")
            for name in ARTIFACTS[1:4]:
                equal_files(output / f"cell-{index}-main.parity",
                            output / f"cell-{index}-{name}.parity", r * size)
                require(all(cell_records[name][key] == cell_records["main"][key]
                            for key in ("input_hash", "parity_hash")), "digest mismatch")
            for name in ("main", "current"):
                require(cell_records[name] == cell_records[name + "-guard"],
                        "guard changed check record")
            records.append(cell_records)
        verify()
    finally:
        (output / "records.json").write_text(json.dumps(records, indent=2) + "\n")
    (output / "result.json").write_text(json.dumps({
        "complete": True, "cells": len(records), "positive_processes": 32,
        "clock_guard_controls": 3, "cli_rejections": 20, "artifact_sha256": pins,
        "candidate_commit": candidate_commit, "baseline_commit": BASELINE,
        "timings": False, "performance_conclusion": None,
    }, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidate-commit", required=True)
    args = parser.parse_args()
    run(args.bundle.resolve(), args.output.resolve(), args.candidate_commit)
