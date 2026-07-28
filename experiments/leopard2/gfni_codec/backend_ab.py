#!/usr/bin/env python3
"""Within-binary backend A/B for the Leopard2 benchmark (non-authoritative)."""
import argparse
import json
import math
import subprocess
import sys

ENV = {"LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
       "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin",
       "TZ": "UTC"}


def run(binary, cpu, cell, backend, cfg):
    k, r, nbytes, loss, field = cell
    cmd = ["taskset", "-c", str(cpu), binary, "--k", str(k), "--r", str(r),
           "--bytes", str(nbytes), "--loss", str(loss), "--batch", "1",
           "--reuse", str(cfg["reuse"]), "--iterations", str(cfg["iters"]),
           "--warmup", str(cfg["warmup"]), "--threads", "1", "--seed", "1",
           "--skip-legacy", "--backend", backend, "--json", "-"]
    if field:
        cmd += ["--field", field]
    out = subprocess.run(cmd, capture_output=True, env=ENV, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-300:])
    return json.loads(out.stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", required=True)
    ap.add_argument("--baseline", default="avx2")
    ap.add_argument("--candidate", default="avx512")
    ap.add_argument("--cells", required=True)
    ap.add_argument("--cpu", type=int, default=14)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--reuse", type=int, default=8)
    ap.add_argument("--iters", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=2)
    args = ap.parse_args()
    cfg = {"reuse": args.reuse, "iters": args.iters, "warmup": args.warmup}
    cells = json.load(open(args.cells))

    print("%-36s %8s %8s %11s %11s %s" % (
        "cell", "encode", "decode", "base_enc_us", "cand_enc_us", "digest"))
    for cell in cells:
        enc, dec = [], []
        digest_ok = True
        base_us = cand_us = 0.0
        for _ in range(args.rounds):
            a1 = run(args.binary, args.cpu, cell, args.baseline, cfg)
            b1 = run(args.binary, args.cpu, cell, args.candidate, cfg)
            b2 = run(args.binary, args.cpu, cell, args.candidate, cfg)
            a2 = run(args.binary, args.cpu, cell, args.baseline, cfg)
            ea = 0.5 * (math.log(a1["metrics"]["encode_execution"]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"]["encode_execution"]["median_us_per_batch_call"]))
            eb = 0.5 * (math.log(b1["metrics"]["encode_execution"]["median_us_per_batch_call"]) +
                        math.log(b2["metrics"]["encode_execution"]["median_us_per_batch_call"]))
            enc.append(ea - eb)
            da = 0.5 * (math.log(a1["metrics"]["decode_execution"]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"]["decode_execution"]["median_us_per_batch_call"]))
            db = 0.5 * (math.log(b1["metrics"]["decode_execution"]["median_us_per_batch_call"]) +
                        math.log(b2["metrics"]["decode_execution"]["median_us_per_batch_call"]))
            dec.append(da - db)
            base_us = a1["metrics"]["encode_execution"]["median_us_per_batch_call"]
            cand_us = b1["metrics"]["encode_execution"]["median_us_per_batch_call"]
            for key in ("original_data", "transmitted_parity", "recovered_originals"):
                if a1["workload_digests"][key] != b1["workload_digests"][key]:
                    digest_ok = False
            if not (a1["correctness"]["leopard2_round_trip"] and
                    b1["correctness"]["leopard2_round_trip"]):
                digest_ok = False
        k, r, nbytes, loss, field = cell
        print("%-36s %8.3f %8.3f %11.2f %11.2f %s" % (
            "K=%d R=%d %s l=%d %s" % (k, r, nbytes, loss, field or ""),
            math.exp(sum(enc) / len(enc)), math.exp(sum(dec) / len(dec)),
            base_us, cand_us, "OK" if digest_ok else "MISMATCH"))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
