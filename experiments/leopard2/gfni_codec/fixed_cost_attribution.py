#!/usr/bin/env python3
"""Separate fixed per-call cost from byte-proportional cost, L2 vs exact main."""
import json
import subprocess
import sys

ENV = {"LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
       "PATH": "/usr/bin:/bin", "TZ": "UTC"}
MAIN = "/tmp/leopard2-final-gf8-main-baseline-a2515d1-20260722/leopard_main_benchmark"
L2 = "/tmp/l2-head-build/bench_leopard2"
CPU = 14


def enc(binary, k, r, nbytes, extra=()):
    cmd = ["taskset", "-c", str(CPU), binary, "--k", str(k), "--r", str(r),
           "--bytes", str(nbytes), "--loss", "1", "--batch", "1",
           "--reuse", "64", "--iterations", "9", "--warmup", "3",
           "--threads", "1", "--seed", "1", "--json", "-"]
    cmd.extend(extra)
    out = subprocess.run(cmd, capture_output=True, env=ENV, text=True)
    if out.returncode != 0:
        return None
    d = json.loads(out.stdout)
    return d["metrics"]["encode_execution"]["median_us_per_batch_call"]


def fit(points):
    # least squares on (bytes, us)
    n = len(points)
    sx = sum(p[0] for p in points)
    sy = sum(p[1] for p in points)
    sxx = sum(p[0] * p[0] for p in points)
    sxy = sum(p[0] * p[1] for p in points)
    denom = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return intercept, slope


SIZES = [64, 128, 256, 512, 1024, 2048, 4096]
print("%-12s %8s %10s %10s %10s" % ("shape", "bytes", "main_us", "l2_us", "ratio"))
for (k, r) in [(2, 2), (4, 4), (8, 8), (16, 4), (16, 16), (32, 32)]:
    pts_m, pts_l = [], []
    for nbytes in SIZES:
        m = enc(MAIN, k, r, nbytes)
        l = enc(L2, k, r, nbytes, ("--skip-legacy",))
        if m is None or l is None:
            continue
        pts_m.append((nbytes, m))
        pts_l.append((nbytes, l))
        print("%-12s %8d %10.4f %10.4f %10.3f"
              % ("K=%d R=%d" % (k, r), nbytes, m, l, m / l))
    im, sm = fit(pts_m)
    il, sl = fit(pts_l)
    print("  -> fixed: main %.4f us, leopard2 %.4f us  (delta %.4f us)"
          % (im, il, il - im))
    print("  -> per-KiB: main %.4f us, leopard2 %.4f us"
          % (sm * 1024, sl * 1024))
    sys.stdout.flush()
