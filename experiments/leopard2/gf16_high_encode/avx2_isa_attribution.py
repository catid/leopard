#!/usr/bin/env python3
"""Clock-free object-code attribution for leopard-79h.38.5.4.18.

Static instruction counts are not dynamic counts, cost estimates, or timings.
Never execute the input codecs. Preserve the original native-L1 comparator.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


FUNCTION = re.compile(r"^[0-9a-f]+ <(.+)>:$")
HIGH_YMM = re.compile(r"%ymm(?:1[6-9]|2[0-9]|3[01])\b")


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def summarize(disassembly):
    functions = {}
    counts = Counter()
    current = None
    for line in disassembly.splitlines():
        match = FUNCTION.fullmatch(line)
        if match:
            current = match.group(1)
            if current in functions:
                raise ValueError("duplicate function label")
            functions[current] = {"counts": Counter(), "examples": {}}
            continue
        columns = line.split("\t")
        if len(columns) < 3 or not re.fullmatch(r"\s*[0-9a-f]+:", columns[0]):
            continue
        if current is None:
            raise ValueError("instruction before function label")
        raw = bytes.fromhex(columns[1])
        assembly = columns[2].strip()
        if not raw or not assembly:
            raise ValueError("missing instruction bytes or assembly")
        mnemonic = assembly.split()[0]
        flags = {
            "instructions": True,
            "evex": raw[0] == 0x62,
            "high_ymm": bool(HIGH_YMM.search(assembly)),
            "zmm": bool(re.search(r"%zmm[0-9]+\b", assembly)),
            "ternary_logic": mnemonic in ("vpternlogd", "vpternlogq"),
            "byte_shuffle": mnemonic == "vpshufb",
            "gfni": mnemonic.startswith("vgf2p8"),
        }
        record = functions[current]
        for name, present in flags.items():
            counts[name] += int(present)
            record["counts"][name] += int(present)
            if present and name != "instructions":
                record["examples"].setdefault(name, line.strip())
    if not counts["instructions"]:
        raise ValueError("no instructions parsed")
    return {"counts": dict(counts), "functions": functions}


def select_recipe(recipes, source, object_target):
    selected = [r for r in recipes if Path(r["file"]).name == source
                and r["output"].startswith(f"CMakeFiles/{object_target}.dir/")]
    if len(selected) != 1:
        raise ValueError("ambiguous compile recipe")
    return selected[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve(strict=True)
    args.out.mkdir()  # Refuse to overwrite an earlier audit.
    out = args.out.resolve(strict=True)
    original = repo / ".research/leopard-79h/gf16-current-route-failed.STc10h"
    production = repo / ".research/leopard-79h/auto-gfni-boundary-production.B5NV9R"
    lanes = {
        "native_l1": (
            original / "frozen/main.a",
            "3f83c55599cbfbeea5d15adf3627405bcde1cbed2ddf4a672fbc75b186f5c4c1",
            "LeopardFF16.cpp.o",
            original / "build-metadata/main-build/compile_commands.json",
            "LeopardFF16.cpp", "leopard_main_exact",
        ),
        "measured_l2_avx2": (
            original / "frozen/current.a",
            "259c2c9aa3f51b1f941eac84270ad37c39c2a5b3cc06b1676a414f8088ad88e6",
            "Leopard2BackendAVX2.cpp.o",
            original / "build-metadata/current-build/compile_commands.json",
            "Leopard2BackendAVX2.cpp", "leopard2_backend_avx2",
        ),
        "production_l2_avx2": (
            production / "release/libleopard.a",
            "d8eb0985e7347e491ea069e2e39d64ec256a35808422d02c1df36641fb315c8a",
            "Leopard2BackendAVX2.cpp.o",
            production / "release/compile_commands.json",
            "Leopard2BackendAVX2.cpp", "leopard2_backend_avx2",
        ),
    }
    report = {"bead": "leopard-79h.38.5.4.18", "timings": False, "lanes": {}}
    report["tools"] = {
        tool: subprocess.check_output([tool, "--version"], text=True).splitlines()[0]
        for tool in ("ar", "objdump")
    }
    for lane, (archive, expected, member, recipes, source, target) in lanes.items():
        if sha(archive) != expected:
            raise ValueError(f"original archive changed: {lane}")
        directory = out / lane
        directory.mkdir()
        frozen = directory / "codec.a"
        shutil.copyfile(archive, frozen)
        frozen.chmod(0o444)
        if sha(frozen) != expected:
            raise ValueError(f"copy changed: {lane}")
        members = subprocess.check_output(["ar", "t", str(frozen)], text=True).splitlines()
        if members.count(member) != 1:
            raise ValueError(f"ambiguous member: {lane}")
        obj = directory / member
        with obj.open("xb") as stream:
            subprocess.run(["ar", "p", str(frozen), member], stdout=stream, check=True)
        obj.chmod(0o444)
        obj_hash = sha(obj)
        command = ["objdump", "-d", "-w", "-C", str(obj)]
        raw = subprocess.check_output(command, text=True)
        (directory / "disassembly.txt").write_text(raw)
        shutil.copyfile(recipes, directory / "compile_commands.json")
        selected = select_recipe(json.loads(recipes.read_text()), source, target)
        if sha(obj) != obj_hash or sha(archive) != expected or sha(frozen) != expected:
            raise ValueError(f"input drift: {lane}")
        report["lanes"][lane] = {
            "archive": str(archive), "archive_sha256": expected,
            "member": member, "member_sha256": obj_hash,
            "recipe": selected, "objdump_argv": command,
            **summarize(raw),
        }
    current = report["lanes"]["measured_l2_avx2"]
    latest = report["lanes"]["production_l2_avx2"]
    if current["member_sha256"] != latest["member_sha256"]:
        raise ValueError("production AVX2 changed since measured comparison")
    for lane in (current, latest):
        if any(lane["counts"][k] for k in ("evex", "high_ymm", "zmm", "ternary_logic", "gfni")):
            raise ValueError("restricted AVX2 object has unexpected instructions")
    (out / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    shutil.copyfile(__file__, out / Path(__file__).name)
    manifest = "".join(f"{sha(p)}  {p.relative_to(out)}\n" for p in sorted(out.rglob("*")) if p.is_file())
    (out / "SHA256SUMS").write_text(manifest)
    print(json.dumps({k: {"sha256": v["member_sha256"], **v["counts"]} for k, v in report["lanes"].items()}, indent=2))


if __name__ == "__main__":
    main()
