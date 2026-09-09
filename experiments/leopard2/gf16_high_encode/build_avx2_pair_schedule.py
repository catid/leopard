#!/usr/bin/env python3
"""Untimed, serial codegen experiment for leopard-79h.38.5.4.18.2.

Run inside a 512 MiB/no-swap scope under the canonical build/test lock.
Input is a private copy of production sources with the retained pair overlay.
This does not execute a codec and cannot justify a performance claim.
"""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import shlex
import shutil
import subprocess


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", type=Path)
    args = parser.parse_args()
    root = args.workspace.resolve(strict=True)
    repo = Path(__file__).resolve().parents[3]
    production = repo / ".research/leopard-79h/auto-gfni-boundary-production.B5NV9R"
    source = root / "source"
    original = production / "source"
    member = "Leopard2BackendAVX2.cpp.o"
    cpp = "Leopard2BackendAVX2.cpp"
    header = Path(__file__).with_name("avx2_pair_schedule.h")
    if set(p.name for p in source.iterdir()) != set(p.name for p in original.iterdir()):
        raise ValueError("source inventory mismatch")
    for path in original.iterdir():
        tracked = subprocess.check_output(
            ["git", "show", "3a2f064:" + path.name], cwd=repo)
        checked = path if path.name == cpp else source / path.name
        if hashlib.sha256(tracked).hexdigest() != sha(checked):
            raise ValueError("production source identity: " + path.name)
    out = root / "codegen"
    out.mkdir()  # Never overwrite a prior preparation.
    shutil.copyfile(header, source / header.name)
    shutil.copyfile(__file__, out / Path(__file__).name)
    patch = "".join(difflib.unified_diff(
        (original / cpp).read_text().splitlines(keepends=True),
        (source / cpp).read_text().splitlines(keepends=True),
        fromfile="a/" + cpp, tofile="b/" + cpp))
    (out / "source.patch").write_text(patch)
    recipes_path = production / "release/compile_commands.json"
    recipes = json.loads(recipes_path.read_text())
    choices = [r for r in recipes if r["output"] ==
               "CMakeFiles/leopard2_backend_avx2.dir/" + member]
    if len(choices) != 1:
        raise ValueError("ambiguous production recipe")
    recipe = choices[0]
    command = shlex.split(recipe["command"])
    command[command.index("-o") + 1] = "REPLACE_OUTPUT"
    command[command.index("-c") + 1] = str(source / cpp)
    command[command.index("-I" + str(repo))] = "-I" + str(source)
    inputs = {p.name: sha(p) for p in source.iterdir()}
    state = {"bead": "leopard-79h.38.5.4.18.2", "timed": False,
             "source_base": "3a2f064", "source_pins": inputs,
             "recipe": recipe, "recipe_sha256": sha(recipes_path),
             "compiler": subprocess.check_output([command[0], "--version"], text=True),
             "modes": {}}
    try:
        for mode in range(3):
            directory = out / str(mode)
            directory.mkdir()
            obj = directory / member
            argv = [str(obj) if arg == "REPLACE_OUTPUT" else arg for arg in command]
            argv.append("-DLEO_AVX2_PAIR_SCHEDULE=" + str(mode))
            state["modes"][str(mode)] = {"command": argv}
            with (directory / "compile.stdout").open("wb") as stdout, \
                    (directory / "compile.stderr").open("wb") as stderr:
                subprocess.run(argv, cwd=directory, stdout=stdout, stderr=stderr, check=True)
            state["modes"][str(mode)]["object_sha256"] = sha(obj)
            with (directory / "disassembly.txt").open("wb") as output:
                subprocess.run(["objdump", "-drwC", str(obj)], stdout=output, check=True)
            if mode == 0 and sha(obj) != \
                    "bf615cd0bb211ea3b8fd8a2c7e6c8a9651f5317857cd8ddcd910c03e74fe1e5d":
                raise ValueError("default-off object differs from production")
            print(json.dumps({"mode": mode, "object_sha256": sha(obj)}), flush=True)
        if inputs != {p.name: sha(p) for p in source.iterdir()}:
            raise ValueError("source changed during build")
    finally:
        (out / "build.json").write_text(json.dumps(state, indent=2) + "\n")


if __name__ == "__main__":
    main()
