#!/usr/bin/env python3
"""Clock-free schedule model and retained-only callback/parity replay (.12)."""
from collections import Counter
import json
from pathlib import Path
import sys

from verify_gfni_source_stage import equal, require, resource_peak

CELLS = ((1000,200,65536,6,32768,2), (1000,200,65536,3,32768,2),
         (1000,200,65536,5,65536,1), (1000,200,32768,3,32768,1),
         (1000,199,65536,3,32768,2), (4096,512,4096,3,4096,1))
OPS = ("multiply", "multiply_add", "xor", "xor_2to1", "xor4", "copy", "ifft2", "fft2",
       "fft2_out", "ifft2_xor", "ifft4", "fft4", "ifft4_out", "fft4_out", "ifft4_range", "fft4_range")


def model(cell):
    return model_shape(CELLS[cell])


def model_shape(shape):
    """Structural high-profile traversal, not a model of CPU time or traffic.

    The supported fixture shapes have complete four-source groups, nonzero inverse
    skews, and payloads above the field's 128-byte fused threshold. The first
    forward group has zero-skew mask5; all remaining groups have mask0.
    A radix-four group represents four logical radix-two edges per lane.
    """
    k, r, public_bytes, kind, tile, passes = shape
    side = 1 << (r - 1).bit_length()
    require(k % 4 == 0 and tile > 128 and kind in (3,5,6), "model scope")
    counts = Counter()

    def add(op, calls, distance=1, mask=0):
        if calls:
            counts[op, distance, mask, False, tile] += calls

    for _ in range(passes):
        for block, start in enumerate(range(0, k, side)):
            available = min(side, k - start)
            staged = public_bytes <= 16384
            distance = 1
            if staged:
                add("ifft4_out", available // 4)
                distance = 4
            while 4 * distance <= side:
                groups = (available + 4 * distance - 1) // (4 * distance)
                if 4 * distance == side and block:
                    add("ifft2", groups * 2 * distance)
                    add("ifft2_xor", groups * 2 * distance)
                elif distance == 1:
                    add("ifft2", groups * 4)
                else:
                    add("ifft4_range", groups, distance)
                distance *= 4
            if distance < side:
                add("ifft2_xor" if block else "ifft2", side // 2)
        distance = side // 4
        while distance:
            groups = (r + 4 * distance - 1) // (4 * distance)
            if distance == 1:
                add("xor", 3)
                add("fft2", 4 * groups - 3)
            else:
                add("fft4_range", 1, distance, 5)
                add("fft4_range", groups - 1, distance)
            distance //= 4
        if (side.bit_length() - 1) % 2:
            add("xor", 1)
            add("fft2", (r + 1) // 2 - 1)
    return counts


def validate_counts(record, cell):
    return validate_shape_counts(record, CELLS[cell])


def validate_shape_counts(record, shape):
    require(set(record) == {"schema", "timed", "calls", "passes", "buckets"} and
            record["schema"] == "gf16-callback-counts/v1" and record["timed"] is False,
            "callback schema")
    k,r,public_bytes,kind,tile,passes = shape
    side = 1 << (r - 1).bit_length()
    equal(record["passes"], [{"kind":kind,"k":k,"r":r,"requested":r,"side":side,
        "sparse_blocks":0,"bytes":tile,"source_policy":public_bytes}] * passes)
    require(type(record["calls"]) is int and 0 < record["calls"] <= 1000000 and
            type(record["buckets"]) is list and 0 < len(record["buckets"]) <= 256, "count bounds")
    actual = Counter()
    for bucket in record["buckets"]:
        require(set(bucket) == {"op","distance","zero_mask","prefer_fused","bytes","calls"}, "bucket schema")
        require(bucket["op"] in OPS and type(bucket["prefer_fused"]) is bool and
                all(type(bucket[key]) is int for key in ("distance","zero_mask","bytes","calls")) and
                bucket["distance"] > 0 and 0 <= bucket["zero_mask"] < 8 and
                bucket["bytes"] > 0 and bucket["calls"] > 0, "bucket values")
        key = tuple(bucket[name] for name in ("op","distance","zero_mask","prefer_fused","bytes"))
        require(key not in actual, "duplicate bucket")
        actual[key] = bucket["calls"]
    require(actual == model_shape(shape) and sum(actual.values()) == record["calls"], "schedule model mismatch")
    return {op: {"calls":sum(n for key,n in actual.items() if key[0] == op),
                 "lane_groups":sum(key[1]*n for key,n in actual.items() if key[0] == op)}
            for op in OPS if any(key[0] == op for key in actual)}


def replay(root, baseline):
    peaks = []
    compared = 0
    summaries = []
    for flavor in ("release", "sanitize"):
        for cell in range(6):
            prefix = root / "checks" / ("%s-%d" % (flavor, cell))
            peaks.append(resource_peak(prefix.with_suffix(".log")))
            require(prefix.with_suffix(".json").read_bytes() ==
                    (baseline / ("%d-current.json" % cell)).read_bytes(), "public workload")
            record = json.loads(prefix.with_suffix(".callbacks.json").read_text())
            summary = validate_counts(record, cell)
            if flavor == "release":
                summaries.append({"cell":cell,"calls":record["calls"],"operations":summary})
                reference = baseline / ("%d-main.parity" % cell)
                parity = prefix.with_suffix(".parity")
                require(parity.stat().st_size == reference.stat().st_size == CELLS[cell][1]*CELLS[cell][2], "parity size")
                with parity.open("rb") as actual, reference.open("rb") as original:
                    while True:
                        block = actual.read(65536)
                        require(block == original.read(65536), "exact Leopard1 parity")
                        if not block:
                            break
                compared += parity.stat().st_size
            else:
                require(prefix.with_suffix(".callbacks.json").read_bytes() ==
                        (root / "checks" / ("release-%d.callbacks.json" % cell)).read_bytes(), "sanitizer trace")
        unit = root / ("unit-%s.log" % flavor)
        peaks.append(resource_peak(unit))
        require(unit.read_text().splitlines().count("callback observer: 16 exact delegations, metadata/null preservation, masks/hints, aggregation, reentry and both bounds passed") == 1, "unit checks")
    require(compared == 61014016, "parity total")
    return {"cells":summaries,"full_parity_bytes":compared,"validated_scopes":len(peaks),
            "peak_bytes":max(peaks),"timed":False,"performance_inference":False}


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_gf16_callback_probe.py ROOT BASELINE")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])), sort_keys=True))
