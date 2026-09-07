#!/usr/bin/env python3
"""Retained-only four-mode correctness replay for leopard-79h.38.5.4.16.

Derives the selected schedule from the independent original traversal model.
Does not execute a codec, import a collector, or estimate performance.
"""
import json
from pathlib import Path
import sys

from verify_gf16_callback_probe import CELLS, model
from verify_gfni_source_stage import equal, require
from verify_gfni_terminal import PUBLIC, resource_result


def scope(cell, mode):
    require(type(cell) is int and 0 <= cell < len(CELLS), "cell scope")
    require(type(mode) is int and 0 <= mode <= 3, "mode scope")


def expected_hook(cell, mode):
    scope(cell, mode)
    k, r, size, kind, tile, passes = CELLS[cell]
    matched = passes if cell == 0 else 0
    return dict(schema="gfni-combined/v1", timed=False, mode=mode,
        calls=passes, matches=matched, first=matched if mode & 1 else 0,
        terminal=matched if mode & 2 else 0,
        records=[dict(kind=kind, k=k, r=r, requested=r, side=1 << (r-1).bit_length(),
            sparse_blocks=0, bytes=tile, source_policy=size, sparse_present=True)] * passes)


def validate_hook(hook, cell, mode):
    equal(hook, expected_hook(cell, mode))


def validate(record, kernel, cell, mode):
    scope(cell, mode)
    k, r, size, kind, tile, passes = CELLS[cell]
    side = 1 << (r-1).bit_length()
    counts = model(cell)
    terminal_calls = terminal_groups = 0
    pair = ("ifft2", 1, 0, False, tile)
    if cell == 0:
        if mode & 1:
            # These fixed shapes have complete first-stage four-source groups.
            groups = passes * k // 4
            counts[pair] -= 4 * groups
            counts["ifft4_range", 1, 0, False, tile] += groups
        if mode & 2:
            terminal_calls = passes * ((k + side - 1) // side - 1)
            terminal_groups = terminal_calls * (side // 4)
            counts[pair] -= 2 * terminal_groups
            counts["ifft2_xor", 1, 0, False, tile] -= 2 * terminal_groups
    require(all(value >= 0 for value in counts.values()), "negative modeled count")
    equal(kernel, dict(schema="gfni-terminal-kernel-counts/v1", timed=False,
        calls=terminal_calls, lane_groups=terminal_groups))
    expected = dict(schema="gf16-callback-counts/v1", timed=False,
        calls=sum(counts.values()), passes=[dict(kind=kind, k=k, r=r, requested=r,
            side=side, sparse_blocks=0, bytes=tile, source_policy=size)] * passes,
        buckets=[dict(op=op, distance=distance, zero_mask=mask, prefer_fused=hint,
            bytes=length, calls=count)
            for (op, distance, mask, hint, length), count in counts.items() if count])
    # Ignore only bucket order. Compare serialized types/keys exactly; do not
    # normalize booleans into integers, merge duplicates or discard zero rows.
    require(type(record) is dict and type(record.get("buckets")) is list,
            "callback record/bucket types")
    canonical = lambda value: json.dumps(value, sort_keys=True)
    actual = dict(record, buckets=sorted(record["buckets"], key=canonical))
    expected["buckets"].sort(key=canonical)
    equal(actual, expected)


def replay(root, baseline):
    peaks, parity_bytes = [], 0
    for flavor in ("release", "sanitize"):
        for mode in range(4):
            for cell, (_, r, size, _, _, _) in enumerate(CELLS):
                path = root / "checks" / ("%s-%d-%d" % (flavor, mode, cell))
                peaks.append(resource_result(path.with_suffix(".log")))
                require(path.with_suffix(".json").read_bytes() ==
                    (baseline / ("%d-current.json" % cell)).read_bytes(), "public workload identity")
                lines = path.with_suffix(".trace.jsonl").read_text().splitlines()
                require(len(lines) == 3, "exact hook/callback/kernel records")
                hook, callbacks, kernel = map(json.loads, lines)
                validate_hook(hook, cell, mode)
                validate(callbacks, kernel, cell, mode)
                if flavor == "sanitize":
                    require(path.with_suffix(".trace.jsonl").read_bytes() ==
                        (root / "checks" / ("release-%d-%d.trace.jsonl" % (mode, cell))).read_bytes(),
                        "sanitizer trace identity")
                else:
                    reference = baseline / ("%d-main.parity" % cell)
                    parity = path.with_suffix(".parity")
                    require(parity.stat().st_size == reference.stat().st_size == r * size,
                            "full parity length")
                    with parity.open("rb") as actual, reference.open("rb") as expected:
                        while True:
                            block = actual.read(65536)
                            require(block == expected.read(65536), "exact Leopard1 parity")
                            if not block:
                                break
                    parity_bytes += r * size
    require(parity_bytes == 244056064, "full parity total")
    for flavor in ("release", "sanitize"):
        for item in ("first", "kernel", *(str(i) for i in range(14))):
            path = root / "directed" / ("%s-%s.log" % (flavor, item))
            peaks.append(resource_result(path))
            lines = [line for line in path.read_text().splitlines()
                     if line.startswith(("combined ", "distance-one "))]
            if item == "first":
                expected = ["distance-one hook: 17 negative predicates, both modes, descriptor identity and atomic 16-pass bound passed",
                    "distance-one GFNI kernel: 480 exact-end/unaligned/zero-skew/hint scalar-pair comparisons passed"]
            elif item == "kernel":
                expected = ["combined hook: 17 negative predicates, four masks, atomic 16-pass bound/reset and 9 CLI guards passed",
                    "combined GFNI kernel: 512 scalar-layer/XOR, zero-skew, exact-end, unaligned, readonly-input and cancellation cases passed"]
            elif int(item) in (1, 2):
                expected = ["combined public %s: native odd rejected in all four modes" % item]
            else:
                k, r, size, backend, partial, passes, changes = PUBLIC[int(item)]
                expected = ["combined public %d: K%d/R%d/B%d backend%d partial%d calls%d first%d terminal%d modes4 parity exact" %
                    (int(item), k, r, size, backend, partial, passes, changes, changes)]
            equal(lines, expected)
    ancillary_peaks = []
    for name, message in (
        ("isa-release.log", "portable ISA check: PASS (SSE2 baseline; named SSSE3/AVX2/GFNI/AVX-512VL/probe members isolated)"),
        ("guards.log", "24 malformed or timing requests rejected before workload execution"),
        ("provenance.log", "All frozen build inputs/artifacts unchanged; combined overlay applies to exact base; production unchanged")):
        path = root / name
        ancillary_peaks.append(resource_result(path))
        require(path.read_text().splitlines().count(message) == 1, "ancillary result")
    return dict(bead="leopard-79h.38.5.4.16", public_records=48,
        full_parity_files=24, full_parity_bytes=parity_bytes,
        modes={"0":"neither", "1":"first_only", "2":"terminal_only", "3":"both"},
        target_ops_callbacks=[4132, 2632, 2596, 1096],
        target_external_terminal_calls=[0, 0, 6, 6],
        kernel_cases_per_build=dict(first=480, terminal=512),
        directed_public_cases_per_build=14, public_modes_per_case=4,
        validated_native_scopes=len(peaks), native_peak_bytes=max(peaks),
        ancillary_peak_bytes=max(ancillary_peaks),
        build_peak_bytes=resource_result(root / "build.log", cap=536870912),
        release_archive_isa="pass", sanitizer_correctness="pass",
        full_build_metadata_qualification=False, timed=False,
        performance_inference=False, production_promotion=False)


if __name__ == "__main__":
    require(len(sys.argv) == 3, "usage: verify_gfni_combined.py ROOT BASELINE")
    print(json.dumps(replay(*(Path(value) for value in sys.argv[1:])), sort_keys=True))
