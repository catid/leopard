#!/usr/bin/env python3
"""Merge repeated GFNI-affine CSV runs into deterministic comparisons."""

import csv
import pathlib
import statistics
import sys


def load(path):
    rows = {}
    with path.open(newline="") as source:
        for row in csv.DictReader(source):
            key = (row["kind"], int(row["bytes"]), int(row["multiplier"]), row["kernel"])
            rows[key] = float(row["median_ns"])
    return rows


def main(argv):
    if len(argv) < 2:
        raise SystemExit("usage: analyze.py RESULTS.csv [RESULTS2.csv ...]")
    runs = [load(pathlib.Path(name)) for name in argv[1:]]
    keys = set.intersection(*(set(run) for run in runs))
    merged = {key: statistics.median(run[key] for run in runs) for key in keys}

    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow([
        "kind", "bytes", "avx512_nibble_median_ns", "gfni_zmm_median_ns",
        "gfni_zmm_improvement_percent", "gfni_ymm_improvement_percent",
        "best_kernel", "meets_5_percent_cell_threshold",
    ])
    for kind in ("mul", "chain"):
        multiplier = 0x5D if kind == "mul" else -1
        sizes = sorted({key[1] for key in keys if key[0] == kind and key[2] == multiplier})
        for size in sizes:
            values = {
                key[3]: value for key, value in merged.items()
                if key[:3] == (kind, size, multiplier)
            }
            baseline = values["avx512_nibble"]
            zmm_gain = (baseline / values["gfni_zmm"] - 1.0) * 100.0
            ymm_gain = (baseline / values["gfni_ymm"] - 1.0) * 100.0
            best = min(values, key=values.get)
            writer.writerow([
                kind, size, f"{baseline:.3f}", f"{values['gfni_zmm']:.3f}",
                f"{zmm_gain:.3f}", f"{ymm_gain:.3f}", best,
                "yes" if zmm_gain >= 5.0 else "no",
            ])


if __name__ == "__main__":
    main(sys.argv)
