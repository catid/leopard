#!/usr/bin/env python3
"""Merge two retained AVX-512 experiment runs deterministically."""
import csv
import statistics
import sys


def load(path):
    rows = {}
    with open(path, newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = (row["operation"], int(row["bytes"]), row["variant"])
            if key in rows:
                raise SystemExit(f"duplicate row {key} in {path}")
            if int(row["samples"]) != 9:
                raise SystemExit(f"unexpected sample count in {path}: {row}")
            rows[key] = row
    return rows


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: analyze.py results_run1.csv results_run2.csv")
    runs = [load(path) for path in sys.argv[1:]]
    if runs[0].keys() != runs[1].keys():
        raise SystemExit("run key sets differ")
    fields = ["operation", "bytes", "avx2_ns", "avx512bw_ns",
              "avx512vbmi_ns", "bw_vs_avx2_pct", "vbmi_vs_avx2_pct",
              "bw_run_spread_pct", "vbmi_run_spread_pct", "best"]
    out = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
    out.writeheader()
    operations = ["copy", "xor", "mul", "radix4", "radix8", "codec"]
    byte_counts = sorted({key[1] for key in runs[0]})
    for operation in operations:
        for byte_count in byte_counts:
            values = {}
            spreads = {}
            for variant in ("avx2", "avx512bw", "avx512vbmi"):
                samples = [float(run[(operation, byte_count, variant)]["median_ns"])
                           for run in runs]
                values[variant] = statistics.median(samples)
                spreads[variant] = abs(samples[0] - samples[1]) / values[variant] * 100.0
            avx2 = values["avx2"]
            out.writerow({
                "operation": operation,
                "bytes": byte_count,
                "avx2_ns": f"{avx2:.3f}",
                "avx512bw_ns": f"{values['avx512bw']:.3f}",
                "avx512vbmi_ns": f"{values['avx512vbmi']:.3f}",
                "bw_vs_avx2_pct": f"{(avx2 / values['avx512bw'] - 1) * 100:.2f}",
                "vbmi_vs_avx2_pct": f"{(avx2 / values['avx512vbmi'] - 1) * 100:.2f}",
                "bw_run_spread_pct": f"{spreads['avx512bw']:.2f}",
                "vbmi_run_spread_pct": f"{spreads['avx512vbmi']:.2f}",
                "best": min(values, key=values.get),
            })


if __name__ == "__main__":
    main()
