#!/usr/bin/env python3
"""Run bounded, clock-free L1 attribution checks (leopard-79h.38.5.4.18).

Run under the canonical lock and a 256MiB/no-swap scope. This is an attribution
comparator, not a replacement for native L1, a benchmark, or L2 qualification.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

from avx2_isa_attribution import sha, summarize


def require(value, message):
    if not value:
        raise ValueError(message)


def manifest(directory):
    for line in (directory / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        require(sha(directory / name) == digest, "manifest mismatch: " + name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    repo = args.repo.resolve(strict=True)
    reference = repo / ".research/leopard-79h/gf16-current-route-failed.STc10h"
    artifacts = root / "check-artifacts"
    out = args.out.absolute()
    out.mkdir()
    manifest(artifacts)
    report = {"bead": "leopard-79h.38.5.4.18", "timings": False,
              "source_commit": "6e5725ebdf9da4370b0bcc4f70fa8eb66f4e6198",
              "source_pins": {}, "objects": {}, "records": [], "comparisons": [],
              "negative_commands": []}
    source = Path("/tmp/leopard-gf16-current-routes.BMj72w/main-source")
    (out / "source").mkdir()
    for stem in ("leopard", "LeopardCommon", "LeopardFF8", "LeopardFF16"):
        for suffix in (".cpp", ".h"):
            name = stem + suffix
            expected = subprocess.check_output(["git", "-C", str(repo), "show",
                                                report["source_commit"] + ":" + name])
            require((source / name).read_bytes() == expected ==
                    (reference / "codec-main" / name).read_bytes(), "source drift")
            shutil.copyfile(source / name, out / "source" / name)
            report["source_pins"][name] = sha(source / name)
    recipes = json.loads((root / "pure-build/compile_commands.json").read_text())
    report["compile_recipes"] = [r for r in recipes if
        r["output"].startswith("CMakeFiles/leopard_main_exact.dir/")]
    require(len(report["compile_recipes"]) == 4, "unexpected build recipe count")
    for recipe in report["compile_recipes"]:
        require("-march=x86-64 -mtune=generic -mavx2 -mno-avx512f" in recipe["command"]
                and "-march=native" not in recipe["command"], "unexpected ISA flags")
    members = subprocess.check_output(["ar", "t", str(artifacts / "pure.a")], text=True).splitlines()
    require(sorted(members) == sorted(n + ".cpp.o" for n in
            ("leopard", "LeopardCommon", "LeopardFF8", "LeopardFF16")), "archive members")
    for member in members:
        obj = out / member
        with obj.open("xb") as stream:
            subprocess.run(["ar", "p", str(artifacts / "pure.a"), member], stdout=stream, check=True)
        raw = subprocess.check_output(["objdump", "-d", "-w", "-C", str(obj)], text=True)
        (out / (member + ".disassembly")).write_text(raw)
        result = summarize(raw)
        require(not any(result["counts"][k] for k in
                        ("evex", "high_ymm", "zmm", "ternary_logic", "gfni")),
                "ISA ceiling violated")
        report["objects"][member] = {"sha256": sha(obj), **result}
    retained_pins = dict((str(Path(name)), digest) for digest, name in
        (line.split("  ", 1) for line in (reference / "SHA256SUMS").read_text().splitlines()))
    for cell in range(6):
        name = f"preflight/{cell}-main.parity"
        require(sha(reference / name) == retained_pins[name], "old parity pin changed")
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    for profile, cells in (("native", (0, 6, 7)), ("pure", range(8))):
        for cell in cells:
            parity = out / f"{profile}-{cell}.parity"
            command = ["prlimit", "--cpu=30:30", str(artifacts / f"check-{profile}"),
                       "--check", str(cell), str(parity)]
            child = subprocess.run(command, env=env, capture_output=True, text=True)
            (out / f"{profile}-{cell}.stdout").write_text(child.stdout)
            (out / f"{profile}-{cell}.stderr").write_text(child.stderr)
            require(child.returncode == 0 and child.stderr == "", "native check failed")
            record = json.loads(child.stdout)
            require(record["timings"] is False and record["cell"] == cell and
                    record["profile"] == profile + "_l1" and record["outer_guards"] is True,
                    "native record mismatch")
            require(parity.stat().st_size == record["parity_bytes"], "short parity")
            report["records"].append({"argv": command, "result": record, "parity_sha256": sha(parity)})
            expected = reference / f"preflight/{cell}-main.parity" if cell < 6 else out / f"native-{cell}.parity"
            if parity != expected:
                subprocess.run(["cmp", str(parity), str(expected)], check=True)
                report["comparisons"].append({"actual": str(parity), "expected": str(expected),
                                              "bytes": parity.stat().st_size,
                                              "sha256": sha(expected)})
        for bad in ([], ["--check"], ["--measure", "0", "unused"],
                    ["--check", "8", "unused"], ["--check", "00", "unused"],
                    ["--check", "-1", "unused"], ["--check", "x", "unused"]):
            command = [str(artifacts / f"check-{profile}"), *bad]
            child = subprocess.run(command, env=env, capture_output=True, text=True)
            require(child.returncode == 1 and not child.stdout and child.stderr,
                    "invalid CLI accepted")
            report["negative_commands"].append({"argv": command, "stderr": child.stderr})
    manifest(artifacts)
    for name, digest in report["source_pins"].items():
        require(sha(source / name) == digest, "post-run source drift")
    report["full_comparison_bytes"] = sum(c["bytes"] for c in report["comparisons"])
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"records": len(report["records"]),
                      "full_comparison_bytes": report["full_comparison_bytes"],
                      "negative_commands": len(report["negative_commands"]),
                      "objects": {k: v["counts"] for k, v in report["objects"].items()}}, indent=2))


if __name__ == "__main__":
    main()
