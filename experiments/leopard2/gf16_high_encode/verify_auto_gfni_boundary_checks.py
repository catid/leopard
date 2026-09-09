#!/usr/bin/env python3
"""Retained-only implementation checkpoint audit; no native execution or clocks."""
import hashlib
import json
from pathlib import Path
import sys


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def same(actual, expected):
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), "typed record differs")


def archive_members(path):
    """Independently hash GNU ar payloads, without invoking ar or the builder."""
    result, long_names = {}, b""
    with path.open("rb") as stream:
        require(stream.read(8) == b"!<arch>\n", "regular ar required")
        while True:
            header = stream.read(60)
            if not header:
                break
            require(len(header) == 60 and header[58:] == b"`\n", "ar header")
            size = int(header[48:58].strip())
            require(0 <= size < 64 * 1024 * 1024, "ar member bound")
            payload = stream.read(size)
            require(len(payload) == size, "ar payload length")
            if size % 2:
                require(stream.read(1) == b"\n", "ar padding")
            name = header[:16].decode("ascii").strip()
            if name in ("/", "/SYM64/"):
                continue
            if name == "//":
                long_names = payload
                continue
            if name.startswith("/"):
                offset = int(name[1:])
                require(0 <= offset < len(long_names), "ar long-name offset")
                name = long_names[offset:].split(b"/\n", 1)[0].decode("ascii")
            else:
                name = name.rstrip("/")
            require(name and Path(name).name == name and name not in result,
                    "unambiguous member names")
            result[name] = hashlib.sha256(payload).hexdigest()
    require(len(result) == 24, "complete object inventory")
    return result


def resources(path, peak, maximum):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1 and "\tExit status: 0" in lines, "scope exit")
    same(lines[lines.index("memory.peak"):], ["memory.peak", str(peak),
        "memory.max", str(maximum), "memory.events", "low 0", "high 0", "max 0",
        "oom 0", "oom_kill 0", "oom_group_kill 0", "memory.swap.current", "0",
        "memory.swap.max", "0"])


def audit(root, reference):
    fixed = {
        "codec-source-v2/leopard2.cpp": "0e4cdc4485e96e6d1ca711ad0b5192925df01e7a4951bca3d40b35e580c061b7",
        "codec-source-v2/Leopard2Direct.h": "0c6ff3efdfbb8754cd5abebaf520f05241f16e4f5487f27b093aae4dbd2ca065",
        "release/candidate.a": "80bffc9e873585d9a18fcf6a294413b8c2a76d84fa6e7f9a437fb96ca8458533",
        "sanitize/candidate.a": "9a87a6baa349195e02c810fab8131fe5d4edfaf86de5a0efd1b7cafedb207c38",
        "release/check": "e7c7a4817553dee4c585eef14d52e7f4e3554b3e2663d679647096a4188f787f",
        "sanitize/check": "b255217014ce05957440c944aff9daa58778a3b30649b79ccebda08ce898c99e",
        "release/production-check": "1772b435ad9f4377bfcd483bc5855890da60bab29e986c9ff07f49fab7a015c8",
        "sanitize/production-check": "5bbf266d3e20769b5985deda83f701a6bb980ea2d4bf49a557d0512bbefbc5e5",
        "drivers/gfni_boundary_screen.cpp": "3e26af733e27413bb4f4f732c3be65748080fc16e8649ded88a7cd112b5b833e",
        "drivers/test_gfni_boundary.cpp": "7fd60aed56c51f7ccb88b8c1bb11e531fc763f91685ec7775b73b979d6931062",
        "drivers/test_auto_gfni_boundary.cpp": "77a1c57162e2c6b9ba17a917be695d20fda1e2756e1dd566fd0c608b465e0f30",
        "drivers/test_auto_gf16_gfni_production.cpp": "cd83bb52a97a4ab747994c3af5490d10b5f741cb80e8a724f9e3e94a06c9e063",
        "build-archive-v2.log": "9831aa5d64c5242e81597bd44eaefd372a06f1d3c291c5bb05714bbf71f5e000",
        "build-checks.log": "e130f6c0f4a748264fa076d7cdfee0dacb820b48a36bdb0718cd0f02fc17949b",
        "checks.log": "aa14291609e73a19b7d47d51d60a1ffe34b6fbd490b98294617268c757f4b157",
    }
    for name, digest in fixed.items():
        require(sha(root / name) == digest, "fixed evidence: " + name)
    for name, digest in {
        "frozen/current.a": "259c2c9aa3f51b1f941eac84270ad37c39c2a5b3cc06b1676a414f8088ad88e6",
        "sanitizer/libleopard.a": "1cac86a52866496279b4b8a08158ebc71438e14b062332f3516ed6beece88b58",
    }.items():
        require(sha(reference / name) == digest, "reference archive: " + name)
    source_names = [Path(name).name for name in
                    (reference / "current-source-inputs.txt").read_text().splitlines()]
    same(sorted(path.name for path in (root / "codec-source-v2").iterdir()), sorted(source_names))
    for name in source_names:
        path = root / "codec-source-v2" / name
        require(not path.is_symlink() and not path.stat().st_mode & 0o222, "immutable source")
        if name not in ("leopard2.cpp", "Leopard2Direct.h"):
            require(sha(path) == sha(reference / "codec-current" / name), "unrelated source changed")
    for profile, original_name in (("release", "frozen/current.a"),
                                   ("sanitize", "sanitizer/libleopard.a")):
        before = archive_members(reference / original_name)
        after = archive_members(root / profile / "candidate.a")
        same(sorted(before), sorted(after))
        require([name for name in before if before[name] != after[name]] == ["leopard2.cpp.o"],
                "only core object may change")
        same(after["leopard2.cpp.o"], sha(root / profile / "leopard2.cpp.o"))
        same(json.loads((root / profile / "members.json").read_text()),
             {name: {"before": before[name], "after": after[name]} for name in before})
    count = 0
    names = []
    for profile in ("release", "sanitize"):
        records = {"production": None, "routes": "--routes"}
        for cell in range(2):
            for kind in ("api", "concurrent"):
                records[f"{kind}-{cell}"] = "--" + kind
            for fault in ("host", "unavailable", "oom", "kat"):
                records[f"fault-{cell}-{fault}"] = "--fault"
        for cell in range(8):
            for mode in range(2):
                records[f"guards-{cell}-{mode}"] = "--guards"
        for suffix, case in records.items():
            name = profile + "-" + suffix
            names.append(name)
            require((root / "checks" / (name + ".stderr")).stat().st_size == 0, "check stderr")
            lines = (root / "checks" / (name + ".stdout")).read_text().splitlines()
            if case is None:
                same(lines, ["Production AUTO GF16 GFNI route passed"])
            else:
                expected = {"schema": "leopard-auto-gfni-boundary-check/v1",
                            "case": case, "passed": True, "timed": False}
                same(json.loads(lines[-1]), expected)
                require(len(lines) == (2 if case == "--guards" else 1), "record count")
                if case == "--guards":
                    cell = int(suffix.split("-")[1])
                    same(json.loads(lines[0]), {
                        "schema": "leopard-gfni-boundary-guards/v1", "cell": cell,
                        "k": 1000 if cell < 6 else 17,
                        "r": (199 if cell % 2 else 200) if cell < 6 else 7,
                        "bytes": ((65536 if cell % 2 else 32768) + (2 if cell >= 4 else 0))
                            if cell < 6 else (65 if cell == 6 else 66),
                        "misalignment": 0 if cell < 2 else 1 if cell < 4 or cell >= 6 else 2,
                        "field": 1 if cell == 6 else 2, "subset_masks": 6,
                        # Ragged FF16 buffers retain the pre-existing 64,000-byte
                        # additional tail staging, also pinned by the earlier
                        # explicit-backend guarded checkpoint.
                        "scratch_bytes": 16808512 if cell < 4 else 16872512 if cell < 6 else 2816,
                        "timed": False})
            count += 1
    same(sorted(path.stem for path in (root / "checks").glob("*.stdout")), sorted(names))
    same(sorted(path.stem for path in (root / "checks").glob("*.stderr")), sorted(names))
    resources(root / "build-archive-v2.log", 181243904, 536870912)
    resources(root / "build-checks.log", 94679040, 536870912)
    resources(root / "checks.log", 175464448, 268435456)
    rejected = 0
    expected = {"passed": True, "timed": False, "cell": 0}
    for key, value in (("passed", False), ("passed", 1), ("timed", True), ("cell", False)):
        changed = dict(expected)
        changed[key] = value
        try:
            same(changed, expected)
        except ValueError:
            rejected += 1
    require(rejected == 4, "typed mutation escaped")
    return {"native_checks": count, "archive_objects_per_build": 24,
            "unchanged_objects_per_build": 23, "changed_object": "leopard2.cpp.o",
            "source_inputs": len(source_names), "rejected_record_mutations": rejected,
            "candidate_default_enabled": False, "timed_invocations": 0,
            "performance_qualified": False, "production_promotion": False}


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_auto_gfni_boundary_checks.py ROOT REFERENCE")
    print(json.dumps(audit(Path(sys.argv[1]), Path(sys.argv[2])), sort_keys=True))
