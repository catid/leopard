#!/usr/bin/env python3
"""Retained-only .13 structural contrast: no codec execution or speed inference."""
import copy
import json
from pathlib import Path
import sys
from verify_gf16_callback_probe import CELLS, validate_counts
from verify_gfni_source_stage import equal, require, resource_peak


def validate(record, cell, enabled):
    # Undo exactly the selected first-stage substitution, then require the
    # independently derived original traversal in every remaining field/bucket.
    original = copy.deepcopy(record)
    if cell == 0 and enabled:
        k, _, _, _, tile, passes = CELLS[cell]
        groups = passes * k // 4
        expected = dict(op="ifft4_range",distance=1,zero_mask=0,
                        prefer_fused=False,bytes=tile,calls=groups)
        candidates = [b for b in original["buckets"] if b.get("op") == "ifft4_range" and b.get("distance") == 1]
        require(len(candidates) == 1, "unique first-stage range bucket")
        equal(candidates[0], expected)
        original["buckets"].remove(candidates[0])
        pairs = [b for b in original["buckets"] if b.get("op") == "ifft2" and b.get("distance") == 1 and b.get("zero_mask") == 0]
        require(len(pairs) == 1 and type(pairs[0]["calls"]) is int and
                type(original["calls"]) is int, "pair/count types")
        pairs[0]["calls"] += 4 * groups
        original["calls"] += 3 * groups
    validate_counts(original, cell)


def replay(root, baseline):
    parity_bytes = 0
    peaks = []
    for flavor in ("release", "sanitize"):
        for enabled in (False, True):
            for cell, (k,r,size,kind,tile,passes) in enumerate(CELLS):
                path = root / "checks" / ("%s-%d-%d" % (flavor,enabled,cell))
                peaks.append(resource_peak(path.with_suffix(".log")))
                require(path.with_suffix(".json").read_bytes() ==
                        (baseline / ("%d-current.json" % cell)).read_bytes(), "public workload unchanged")
                lines = path.with_suffix(".trace.jsonl").read_text().splitlines()
                require(len(lines) == 2, "exact hook and callback records")
                hook, callbacks = map(json.loads, lines)
                matched = passes if cell == 0 else 0
                equal(hook, {"schema":"gfni-inverse-distance1/v1","timed":False,
                    "enabled":enabled,"calls":passes,"matches":matched,"changed":matched if enabled else 0,
                    "records":[dict(kind=kind,k=k,r=r,requested=r,side=512 if cell==5 else 256,
                        sparse_blocks=0,bytes=tile,source_policy=size,sparse_present=True)] * passes})
                validate(callbacks, cell, enabled)
                if flavor == "sanitize":
                    require(path.with_suffix(".trace.jsonl").read_bytes() ==
                            (root / "checks" / ("release-%d-%d.trace.jsonl" % (enabled,cell))).read_bytes(),
                            "sanitizer trace identity")
                else:
                    reference = baseline / ("%d-main.parity" % cell)
                    parity = path.with_suffix(".parity")
                    require(parity.stat().st_size == reference.stat().st_size == r*size, "parity length")
                    with parity.open("rb") as actual, reference.open("rb") as expected:
                        while True:
                            block = actual.read(65536)
                            require(block == expected.read(65536), "exact Leopard1 parity")
                            if not block: break
                    parity_bytes += r*size
    require(parity_bytes == 122028032, "full parity total")
    directed = {}
    for flavor in ("release", "sanitize"):
        for item in ("kernel", *(str(i) for i in range(14))):
            path = root / "directed" / ("%s-%s.log" % (flavor,item))
            peaks.append(resource_peak(path))
            lines = [line for line in path.read_text().splitlines() if line.startswith("distance-one ")]
            if item == "kernel":
                equal(lines, ["distance-one hook: 17 negative predicates, both modes, descriptor identity and atomic 16-pass bound passed",
                    "distance-one GFNI kernel: 480 exact-end/unaligned/zero-skew/hint scalar-pair comparisons passed"])
            else:
                index = int(item)
                if index in (1,2):
                    equal(lines,["distance-one public %d: native odd rejected in both modes" % index])
                else:
                    changes = 2 if index in (0,3,4) else 0
                    passes = 3 if index in (4,13) else 1 if index in (10,11) else 2
                    require(len(lines) == 1 and lines[0].startswith("distance-one public %d: K" % index) and
                        lines[0].endswith(" calls%d changed%d parity exact" % (passes,changes)), "directed result/counts")
            if flavor == "release": directed[item] = lines
            else: equal(lines,directed[item])
    return dict(public_records=24,full_parity_bytes=parity_bytes,peak_bytes=max(peaks),
                target_callbacks_off=4132,target_callbacks_on=2632,
                first_inverse_pairs_replaced=2000,first_inverse_range_calls=500,
                directed_public_cases_per_build=14,kernel_cases_per_build=480,validated_scopes=len(peaks),
                timed=False,performance_inference=False)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_gfni_inverse_distance1.py ROOT BASELINE")
    print(json.dumps(replay(*(Path(v) for v in sys.argv[1:])), sort_keys=True))
