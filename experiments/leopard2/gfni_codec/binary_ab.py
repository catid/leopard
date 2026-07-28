#!/usr/bin/env python3
"""Clustered-ABBA A/B between two Leopard2 benchmark binaries (diagnostic)."""
import argparse
import json
import math
import subprocess
import sys

ENV = {"LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
       "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin",
       "TZ": "UTC"}


def run(binary, cpu, cell, cfg):
    k, r, nbytes, loss, field = cell
    cmd = ["taskset", "-c", str(cpu), binary, "--k", str(k), "--r", str(r),
           "--bytes", str(nbytes), "--loss", str(loss), "--batch", "1",
           "--reuse", str(cfg["reuse"]), "--iterations", str(cfg["iters"]),
           "--warmup", str(cfg["warmup"]), "--threads", "1", "--seed", "1",
           "--skip-legacy", "--json", "-"]
    if field:
        cmd += ["--field", field]
    out = subprocess.run(cmd, capture_output=True, env=ENV, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-300:])
    return json.loads(out.stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--cells", required=True)
    ap.add_argument("--cpu", type=int, default=14)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--reuse", type=int, default=8)
    ap.add_argument("--iters", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--output", default="")
    args = ap.parse_args()
    cfg = {"reuse": args.reuse, "iters": args.iters, "warmup": args.warmup}
    cells = json.load(open(args.cells))
    results = []
    print("%-34s %8s %8s %11s %11s %s" % (
        "cell", "encode", "decode", "base_us", "cand_us", "digest"))
    for cell in cells:
        enc, dec = [], []
        ok = True
        base_us = cand_us = 0.0
        for _ in range(args.rounds):
            a1 = run(args.baseline, args.cpu, cell, cfg)
            b1 = run(args.candidate, args.cpu, cell, cfg)
            b2 = run(args.candidate, args.cpu, cell, cfg)
            a2 = run(args.baseline, args.cpu, cell, cfg)
            key = "encode_execution"
            ea = 0.5 * (math.log(a1["metrics"][key]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"][key]["median_us_per_batch_call"]))
            eb = 0.5 * (math.log(b1["metrics"][key]["median_us_per_batch_call"]) +
                        math.log(b2["metrics"][key]["median_us_per_batch_call"]))
            enc.append(ea - eb)
            key = "decode_execution"
            da = 0.5 * (math.log(a1["metrics"][key]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"][key]["median_us_per_batch_call"]))
            db = 0.5 * (math.log(b1["metrics"][key]["median_us_per_batch_call"]) +
                        math.log(b2["metrics"][key]["median_us_per_batch_call"]))
            dec.append(da - db)
            base_us = a1["metrics"]["encode_execution"]["median_us_per_batch_call"]
            cand_us = b1["metrics"]["encode_execution"]["median_us_per_batch_call"]
            for k2 in ("original_data", "transmitted_parity", "recovered_originals"):
                if a1["workload_digests"][k2] != b1["workload_digests"][k2]:
                    ok = False
        k, r, nbytes, loss, field = cell
        entry = {"cell": {"k": k, "r": r, "bytes": nbytes, "loss": loss,
                          "field": field},
                 "encode_ratio": math.exp(sum(enc) / len(enc)),
                 "decode_ratio": math.exp(sum(dec) / len(dec)),
                 "encode_rounds": [math.exp(x) for x in enc],
                 "digests_match": ok}
        results.append(entry)
        print("%-34s %8.4f %8.4f %11.2f %11.2f %s" % (
            "K=%d R=%d %s l=%d %s" % (k, r, nbytes, loss, field or ""),
            entry["encode_ratio"], entry["decode_ratio"], base_us, cand_us,
            "OK" if ok else "MISMATCH"))
        sys.stdout.flush()
    if args.output:
        json.dump(results, open(args.output, "w"), indent=1)


if __name__ == "__main__":
    main()
