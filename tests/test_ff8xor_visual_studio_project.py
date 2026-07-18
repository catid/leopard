#!/usr/bin/env python3

"""Keep mandatory FF8 XOR fallback objects in the Visual Studio project."""

from pathlib import Path
import sys
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT / "proj" / "Leopard.vcxproj"
FILTERS = ROOT / "proj" / "Leopard.vcxproj.filters"
MSBUILD = "{http://schemas.microsoft.com/developer/msbuild/2003}"


def includes(path: Path, item: str):
    root = ET.parse(str(path)).getroot()
    return {
        node.attrib["Include"]
        for node in root.iter(MSBUILD + item)
        if "Include" in node.attrib
    }


def main() -> int:
    required = {
        "ClInclude": {
            r"..\LeopardFF8XorAVX512Four.h",
            r"..\LeopardFF8XorDerivative.h",
        },
        "ClCompile": {r"..\LeopardFF8XorAVX512Four.cpp"},
    }
    failures = []
    for item, expected in required.items():
        for path in (PROJECT, FILTERS):
            missing = expected - includes(path, item)
            if missing:
                failures.append(
                    "%s is missing %s %s"
                    % (path.relative_to(ROOT), item, sorted(missing)))

    if failures:
        for failure in failures:
            print("FAIL:", failure, file=sys.stderr)
        print(
            "LeopardFF8Xor.cpp directly references the four-buffer namespace; "
            "the separate translation unit supplies required portable stubs "
            "when AVX-512 generation is disabled.",
            file=sys.stderr,
        )
        return 1

    print("Visual Studio FF8 XOR fallback manifest: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
