#!/usr/bin/env python3
"""Lightweight NON-AUTHORITATIVE screening harness: leopard2 vs exact-main.

Saturates the allowed physical cores (one logical CPU per core, SMT sibling
left idle) and runs clustered ABBA rounds per cell.  This is a diagnostic to
locate regions worth an isolated campaign, not promotion evidence.
"""
import argparse
import concurrent.futures
import json
import math
import os
import subprocess
import sys

MAIN = "/tmp/leopard2-final-gf8-main-baseline-a2515d1-20260722/leopard_main_benchmark"
L2 = "/tmp/l2-sess-prod/bench_leopard2"

ENV = {
    "LANG": "C", "LC_ALL": "C", "OMP_DYNAMIC": "FALSE", "OMP_NUM_THREADS": "1",
    "OMP_PLACES": "cores", "OMP_PROC_BIND": "TRUE", "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}


def run(binary, cpu, k, r, nbytes, loss, reuse, iters, warmup, seed, extra=()):
    cmd = ["taskset", "-c", str(cpu), binary,
           "--k", str(k), "--r", str(r), "--bytes", str(nbytes),
           "--loss", str(loss), "--batch", "1", "--reuse", str(reuse),
           "--iterations", str(iters), "--warmup", str(warmup),
           "--threads", "1", "--seed", str(seed), "--json", "-"]
    cmd.extend(extra)
    out = subprocess.run(cmd, capture_output=True, env=ENV, text=True)
    if out.returncode != 0:
        raise RuntimeError("%s failed: %s" % (binary, out.stderr[-400:]))
    return json.loads(out.stdout)


def cell_worker(job):
    cpu, cell, cfg = job
    k, r, nbytes, loss, field = cell
    rounds = cfg["rounds"]
    enc_ratios, dec_ratios = [], []
    digest = None
    l2_meta = {}
    try:
        extra_l2 = ["--skip-legacy"]
        if field:
            extra_l2 += ["--field", field]
        for _ in range(rounds):
            a1 = run(MAIN, cpu, k, r, nbytes, loss, cfg["reuse"], cfg["iters"],
                     cfg["warmup"], cfg["seed"])
            b1 = run(L2, cpu, k, r, nbytes, loss, cfg["reuse"], cfg["iters"],
                     cfg["warmup"], cfg["seed"], extra_l2)
            b2 = run(L2, cpu, k, r, nbytes, loss, cfg["reuse"], cfg["iters"],
                     cfg["warmup"], cfg["seed"], extra_l2)
            a2 = run(MAIN, cpu, k, r, nbytes, loss, cfg["reuse"], cfg["iters"],
                     cfg["warmup"], cfg["seed"])
            for x in (a1, a2, b1, b2):
                pass
            ma = 0.5 * (math.log(a1["metrics"]["encode_execution"]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"]["encode_execution"]["median_us_per_batch_call"]))
            mb = 0.5 * (math.log(b1["metrics"]["encode_execution"]["median_us_per_batch_call"]) +
                        math.log(b2["metrics"]["encode_execution"]["median_us_per_batch_call"]))
            enc_ratios.append(ma - mb)
            da = 0.5 * (math.log(a1["metrics"]["decode_including_setup"]["median_us_per_batch_call"]) +
                        math.log(a2["metrics"]["decode_including_setup"]["median_us_per_batch_call"]))
            db = 0.5 * (math.log(b1["metrics"]["decode_amortized_at_reuse"]["derived_median_us_per_batch_call"]) +
                        math.log(b2["metrics"]["decode_amortized_at_reuse"]["derived_median_us_per_batch_call"]))
            dec_ratios.append(da - db)
            dg_a = a1["workload_digests"]
            dg_b = b1["workload_digests"]
            same = all(dg_a[key] == dg_b[key] for key in
                       ("original_data", "transmitted_parity", "recovered_originals"))
            digest = same if digest is None else (digest and same)
            l2_meta = {
                "field": b1["resolved"]["field"],
                "profile": b1["resolved"]["profile"],
                "backend": b1["resolved"]["backend"],
                "parent": b1["resolved"]["parent_count"],
                "enc_us": b1["metrics"]["encode_execution"]["median_us_per_batch_call"],
                "main_enc_us": a1["metrics"]["encode_execution"]["median_us_per_batch_call"],
                "dec_us": b1["metrics"]["decode_execution"]["median_us_per_batch_call"],
                "plan_us": b1["metrics"]["decode_plan_setup"]["median_us"],
                "main_dec_us": a1["metrics"]["decode_including_setup"]["median_us_per_batch_call"],
            }
        return {
            "cell": {"k": k, "r": r, "bytes": nbytes, "loss": loss, "field": field},
            "encode_ratio": math.exp(sum(enc_ratios) / len(enc_ratios)),
            "decode_ratio": math.exp(sum(dec_ratios) / len(dec_ratios)),
            "encode_rounds": [math.exp(x) for x in enc_ratios],
            "decode_rounds": [math.exp(x) for x in dec_ratios],
            "digests_match": digest,
            "meta": l2_meta,
        }
    except Exception as exc:  # noqa: BLE001
        return {"cell": {"k": k, "r": r, "bytes": nbytes, "loss": loss, "field": field},
                "error": str(exc)}


def main():
    global MAIN, L2
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="JSON file with cell list")
    ap.add_argument("--output", required=True)
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--reuse", type=int, default=8)
    ap.add_argument("--iters", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--main", default=MAIN)
    ap.add_argument("--l2", default=L2)
    ap.add_argument("--cpus", default="", help="comma list; default all cores < 15")
    args = ap.parse_args()

    MAIN, L2 = args.main, args.l2

    with open(args.cells) as handle:
        cells = json.load(handle)

    allowed = sorted(os.sched_getaffinity(0))
    if args.cpus:
        cpus = [int(x) for x in args.cpus.split(",")]
    else:
        cpus = [c for c in allowed if c < 15]
    cfg = {"rounds": args.rounds, "reuse": args.reuse, "iters": args.iters,
           "warmup": args.warmup, "seed": args.seed}
    jobs = [(cpus[i % len(cpus)], tuple(cell), cfg) for i, cell in enumerate(cells)]
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        for res in pool.map(cell_worker, jobs):
            results.append(res)
            sys.stderr.write(".")
            sys.stderr.flush()
    sys.stderr.write("\n")
    with open(args.output, "w") as handle:
        json.dump(results, handle, indent=1)

    ok = [r for r in results if "error" not in r]
    ok.sort(key=lambda r: r["encode_ratio"])
    print("%-34s %10s %10s %6s" % ("cell", "encode", "decode", "digest"))
    for r in ok:
        c = r["cell"]
        print("%-34s %10.4f %10.4f %6s" % (
            "K=%d R=%d %s B loss=%d %s" % (c["k"], c["r"], c["bytes"], c["loss"],
                                           r["meta"].get("field", "?")),
            r["encode_ratio"], r["decode_ratio"], r["digests_match"]))
    for r in results:
        if "error" in r:
            print("ERROR", r["cell"], r["error"][:200])


if __name__ == "__main__":
    main()
