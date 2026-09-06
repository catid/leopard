#!/usr/bin/env python3
"""Offline replay of the one failed .38.5.4.10 diagnostic; executes no codec.

Independent of collector imports. This proves retained failure/parity evidence,
not performance, arbitrary bundle security, or v19 build/runtime qualification.
"""
import hashlib
import json
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(65536), b""):
            result.update(block)
    return result.hexdigest()


def read_json(path):
    require(path.stat().st_size < 1048576, "oversized JSON")
    return json.loads(path.read_text())


def replay(root):
    frozen = root / "frozen"
    fixed = {
        "attempt1/attempt.json": "f2089e56a1e787bc0dd41ce7fd3c6ec6b54bc89d0b1ef759d403f060fbf4af51",
        "attempt1.log": "78510b95f0b065afc906813e65d3f2e89dcc9b1cd00a14ede5dcc0b69977c3f2",
        "frozen/pins.json": "51e5417317d37ae99a7f16174da9239f17870541371d688299995aa8376d0b69",
    }
    for name, expected in fixed.items():
        require(digest(root / name) == expected, "fixed evidence SHA: " + name)
    pins = read_json(frozen / "pins.json")
    require(len(pins["files"]) == 13, "frozen inventory size")
    for name, expected in pins["files"].items():
        require(Path(name).name == name, "non-flat frozen name")
        require(digest(frozen / name) == expected, "artifact SHA: " + name)
    attempt = read_json(root / "attempt1/attempt.json")
    require(attempt["pins"] == pins, "journal pins")
    require(attempt["plan_sha256"] == pins["files"]["current_route_screen_plan.json"],
            "journal plan")
    require(attempt["complete"] is False and attempt["invocations"] == [] and
            "analysis" not in attempt, "failure cannot support timings")
    require(attempt["failure"] == "ValueError: passive sibling activity; attempt stopped",
            "terminal reason")
    require(attempt["passive"] == {"before": 191252, "after": 191260,
            "elapsed_ns": 10000058885}, "passive activity/window")
    expected = read_json(frozen / "expected.json")
    require(len(attempt["preflight"]) == 12 and
            set(expected) == {"main", "current"} and
            all(len(records) == 6 for records in expected.values()), "check count")
    parity_bytes = 0
    for cell in range(6):
        for slot, variant in enumerate(("main", "current")):
            record = dict(expected[variant][cell], samples_ns=[])
            require(attempt["preflight"][2 * cell + slot] == record, "journal check")
            prefix = root / "attempt1" / ("check-%d-%s" % (cell, variant))
            require(read_json(prefix.with_suffix(".stdout")) == record, "raw check")
            require(prefix.with_suffix(".stderr").read_bytes() == b"", "raw stderr")
            local = root / "preflight" / ("%d-%s.json" % (cell, variant))
            require(read_json(local) == record, "local release check")
        current = root / "preflight" / ("%d-current.json" % cell)
        sanitizer = root / "preflight" / ("%d-sanitize.json" % cell)
        require(current.read_bytes() == sanitizer.read_bytes(), "sanitizer check")
        count = expected["current"][cell]["r"] * expected["current"][cell]["bytes"]
        paths = [root / "preflight" / ("%d-%s.parity" % (cell, v))
                 for v in ("main", "current")]
        require(all(path.stat().st_size == count for path in paths), "parity length")
        with paths[0].open("rb") as main, paths[1].open("rb") as candidate:
            while True:
                block = main.read(65536)
                require(block == candidate.read(65536), "full parity bytes")
                if not block:
                    break
        parity_bytes += count
    require(parity_bytes == 61014016, "total parity bytes")
    return {"checks": 12, "sanitizer_checks": 6, "parity_pairs": 6,
            "compared_parity_bytes": parity_bytes, "timed_invocations": 0,
            "sibling_nonidle_jiffies": 8, "performance_inference": False}


if __name__ == "__main__":
    require(len(sys.argv) == 2, "usage: replay_current_route_failure.py ARCHIVE")
    print(json.dumps(replay(Path(sys.argv[1])), sort_keys=True))
