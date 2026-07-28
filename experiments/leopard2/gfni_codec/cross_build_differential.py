#!/usr/bin/env python3
"""Randomized cross-build differential: stock Leopard2 vs GFNI Leopard2.

Compares original/parity/recovered FNV digests and the round-trip flag for
random (K, R, bytes, loss, field, backend) shapes.  Any digest difference is a
wire-compatibility break in the GFNI kernels.
"""
import argparse
import json
import random
import subprocess
import sys

ENV = {"LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
       "PATH": "/usr/bin:/bin", "TZ": "UTC"}


def run(binary, cpu, shape):
    k, r, nbytes, loss, field, backend, threads, batch = shape
    cmd = ["taskset", "-c", str(cpu), binary,
           "--k", str(k), "--r", str(r), "--bytes", str(nbytes),
           "--loss", str(loss), "--batch", str(batch), "--reuse", "1",
           "--iterations", "1", "--warmup", "0", "--threads", str(threads),
           "--seed", str(k * 7919 + r * 131 + loss), "--skip-legacy",
           "--backend", backend, "--json", "-"]
    if field:
        cmd += ["--field", field]
    out = subprocess.run(cmd, capture_output=True, env=ENV, text=True)
    if out.returncode != 0:
        return {"error": out.stderr.strip()[-200:]}
    d = json.loads(out.stdout)
    return {
        "digests": d["workload_digests"],
        "round_trip": d["correctness"]["leopard2_round_trip"],
        "resolved": d["resolved"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stock", required=True)
    ap.add_argument("--gfni", required=True)
    ap.add_argument("--count", type=int, default=300)
    ap.add_argument("--cpu", type=int, default=14)
    ap.add_argument("--seed", type=int, default=20260724)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    sizes = [1, 2, 3, 31, 32, 33, 63, 64, 65, 100, 127, 128, 255, 256, 500,
             1000, 1024, 4096, 4097, 8192, 65536]
    backends = ["auto", "scalar", "ssse3", "avx2", "avx512"]
    mismatches = 0
    errors = 0
    for index in range(args.count):
        k = rng.choice([1, 2, 3, 4, 5, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65,
                        100, 127, 128, 129, 200, 240, 255, 256, 300, 512,
                        1000, 2000])
        r = rng.choice([1, 2, 3, 4, 5, 8, 16, 17, 32, 64, 65, 100, 128, 129,
                        200, 255, 256, 300])
        if k + r > 65536:
            continue
        loss = rng.randint(0, min(r, k))
        nbytes = rng.choice(sizes)
        field = rng.choice([None, None, "gf8", "gf16"])
        backend = rng.choice(backends)
        threads = rng.choice([1, 1, 1, 2])
        batch = rng.choice([1, 1, 1, 2])
        shape = (k, r, nbytes, loss, field, backend, threads, batch)
        a = run(args.stock, args.cpu, shape)
        b = run(args.gfni, args.cpu, shape)
        if "error" in a and "error" in b:
            # Both reject the shape identically; that is agreement.
            if a["error"] != b["error"]:
                print("ERROR TEXT DIFFERS", shape, a["error"], b["error"])
                mismatches += 1
            errors += 1
            continue
        if ("error" in a) != ("error" in b):
            print("ACCEPTANCE DIFFERS", shape, a, b)
            mismatches += 1
            continue
        if a != b:
            print("MISMATCH", shape)
            print("  stock", a)
            print("  gfni ", b)
            mismatches += 1
        if (index + 1) % 50 == 0:
            sys.stderr.write("%d/%d checked, %d rejected shapes, %d mismatches\n"
                             % (index + 1, args.count, errors, mismatches))
            sys.stderr.flush()
    print("checked=%d rejected_shapes=%d mismatches=%d"
          % (args.count, errors, mismatches))
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
