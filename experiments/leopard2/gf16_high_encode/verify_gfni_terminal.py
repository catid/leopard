#!/usr/bin/env python3
"""Retained-only .14 replay: exact terminal substitution, no speed inference."""
import copy
import json
from pathlib import Path
import sys
from verify_gf16_callback_probe import CELLS, validate_counts
from verify_gfni_source_stage import equal, require, resource_peak


def terminal_shape(cell):
    k, r, _, _, tile, passes = CELLS[cell]
    side = 1 << (r - 1).bit_length()
    calls = passes * ((k + side - 1) // side - 1)
    groups = calls * (side // 4)
    return calls, groups, groups * 4 * tile * 2


def validate(record, kernel, cell, enabled):
    selected = cell == 0 and enabled
    calls, groups, _ = terminal_shape(cell)
    equal(kernel, dict(schema="gfni-terminal-kernel-counts/v1", timed=False,
                       calls=calls if selected else 0, lane_groups=groups if selected else 0))
    # Restore only the removed terminal pair buckets, then check every bucket
    # against the independently derived original traversal. First stage stays
    # split; the external kernel's calls are not part of the Ops total.
    original = copy.deepcopy(record)
    if selected:
        pairs = [b for b in original["buckets"] if b.get("op") == "ifft2"]
        require(len(pairs) == 1 and type(pairs[0]["calls"]) is int and
                type(original["calls"]) is int, "pair/count types")
        pairs[0]["calls"] += 2 * groups
        original["buckets"].append(dict(op="ifft2_xor", distance=1,
            zero_mask=0, prefer_fused=False, bytes=CELLS[cell][4], calls=2*groups))
        original["calls"] += 4 * groups
    validate_counts(original, cell)


PUBLIC = (
    (1000,200,65536,0,0,2,2), (1000,200,65535,0,0,0,0),
    (1000,200,65537,0,0,0,0), (1000,200,65536,6,0,2,2),
    (1000,200,65538,6,0,3,2), (1000,200,65536,0,1,2,0),
    (1000,200,65536,6,1,2,0), (1000,199,65536,0,0,2,0),
    (999,200,65536,6,0,2,0), (1000,201,65536,6,0,2,0),
    (1000,200,32768,6,0,1,0), (1000,200,65536,5,0,1,0),
    (1000,200,65534,0,0,2,0), (1000,200,65538,0,0,3,0))


def resource_result(path, status=0, cap=268435456):
    lines = path.read_text().splitlines()
    require(lines.count("memory.peak") == 1, "resource block count")
    index = lines.index("memory.peak")
    peak = int(lines[index+1])
    require(0 < peak <= cap and lines[index+2:] == [
        "memory.max", str(cap), "memory.events", "low 0", "high 0", "max 0",
        "oom 0", "oom_kill 0", "oom_group_kill 0", "memory.swap.current", "0",
        "memory.swap.max", "0"], "resource envelope")
    equal([line for line in lines if line.startswith("\tExit status:")], ["\tExit status: %d" % status])
    return peak


def ancillary(root):
    peaks = []
    for name in ("isa-release.log", "isa-sanitize-modified.log"):
        path = root / name
        peaks.append(resource_result(path))
        equal([line for line in path.read_text().splitlines() if line.startswith("portable ISA check:")],
            ["portable ISA check: PASS (SSE2 baseline; named SSSE3/AVX2/GFNI/AVX-512VL/probe members isolated)"])
    # Retain, do not turn into a pass, the separately tracked full-sanitizer
    # qualification failure. Baseline and candidate have the same unchanged
    # offending member; provenance.sh additionally compares its exact bytes.
    for name in ("isa-sanitize.log", "isa-sanitize-baseline.log"):
        path = root / name
        peaks.append(resource_result(path, 1))
        equal([line for line in path.read_text().splitlines() if line.startswith("portable ISA check:")],
            ["portable ISA check: AVX-512VL member widened beyond YMM: Leopard2BackendAVX512.cpp.o"])
    guard = root / "guards.log"
    peaks.append(resource_result(guard))
    require(guard.read_text().splitlines().count(
        "14 malformed or timing requests rejected before workload execution") == 1, "guard result")
    provenance = root / "provenance.log"
    peaks.append(resource_result(provenance))
    require(provenance.read_text().splitlines().count(
        "All frozen build inputs/artifacts unchanged; field overlay applies to exact base; production unchanged") == 1,
        "post-check provenance result")
    build_peak = max(resource_result(root / name, cap=536870912)
                     for name in ("build.log", "build-unit.log"))
    return dict(release_archive_isa="pass",modified_sanitizer_members_isa="pass",
                whole_sanitizer_archive_isa="fails_unchanged_AVX512_width",
                sanitizer_isa_followup="leopard-79h.38.5.4.15",full_build_metadata_qualification=False,
                ancillary_peak_bytes=max(peaks),build_peak_bytes=build_peak)


def replay(root, baseline):
    parity_bytes, peaks = 0, []
    for flavor in ("release", "sanitize"):
        for enabled in (False, True):
            for cell, (k,r,size,kind,tile,passes) in enumerate(CELLS):
                path = root / "checks" / ("%s-%d-%d" % (flavor,enabled,cell))
                peaks.append(resource_peak(path.with_suffix(".log")))
                require(path.with_suffix(".json").read_bytes() ==
                        (baseline / ("%d-current.json" % cell)).read_bytes(), "public workload unchanged")
                lines = path.with_suffix(".trace.jsonl").read_text().splitlines()
                require(len(lines) == 3, "exact hook, callback and terminal records")
                hook, callbacks, kernel = map(json.loads, lines)
                matched = passes if cell == 0 else 0
                equal(hook, {"schema":"gfni-terminal/v1","timed":False,
                    "enabled":enabled,"calls":passes,"matches":matched,"changed":matched if enabled else 0,
                    "records":[dict(kind=kind,k=k,r=r,requested=r,side=512 if cell==5 else 256,
                        sparse_blocks=0,bytes=tile,source_policy=size,sparse_present=True)] * passes})
                validate(callbacks, kernel, cell, enabled)
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
    for flavor in ("release", "sanitize"):
        for item in ("kernel", *(str(i) for i in range(14))):
            path = root / "directed" / ("%s-%s.log" % (flavor,item))
            peaks.append(resource_peak(path))
            lines = [line for line in path.read_text().splitlines() if line.startswith("terminal ")]
            if item == "kernel":
                expected = ["terminal hook: 17 negative predicates, both modes, descriptor identity and atomic 16-pass bound passed",
                    "terminal GFNI kernel: 512 scalar-layer/XOR, zero-skew, exact-end, unaligned, readonly-input and cancellation cases passed"]
            elif int(item) in (1,2):
                expected = ["terminal public %s: native odd rejected in both modes" % item]
            else:
                expected = ["terminal public %d: K%d/R%d/B%d backend%d partial%d calls%d changed%d parity exact" %
                            (int(item), *PUBLIC[int(item)])]
            equal(lines, expected)
    return dict(public_records=24,full_parity_bytes=parity_bytes,peak_bytes=max(peaks),
                target_callbacks_off=4132,target_callbacks_on=2596,
                terminal_pairs_removed=768,terminal_accumulating_pairs_removed=768,
                terminal_kernel_calls=6,terminal_four_way_groups=384,
                first_stage_pairs_retained=2000,logical_scratch_bytes_avoided=terminal_shape(0)[2],
                directed_public_cases_per_build=14,kernel_cases_per_build=512,validated_scopes=len(peaks),
                timed=False,performance_inference=False,qualification=ancillary(root))


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_gfni_terminal.py ROOT BASELINE")
    print(json.dumps(replay(*(Path(v) for v in sys.argv[1:])), sort_keys=True))
