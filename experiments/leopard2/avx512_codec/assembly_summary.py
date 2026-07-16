#!/usr/bin/env python3
"""Report code size and conservative stack-vector-reference counts."""
import re
import subprocess
import sys


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: assembly_summary.py BINARY")
    binary = sys.argv[1]
    nm = subprocess.check_output(["nm", "-S", "--defined-only", "-C", binary], text=True)
    sizes = {}
    for line in nm.splitlines():
        parts = line.split(maxsplit=3)
        match = re.search(r"::(Avx(?:2|512)\w+)\(", parts[3]) if len(parts) == 4 else None
        if match and re.search(r"(Copy|Xor|Mul|Radix|Codec)", match.group(1)):
            sizes[match.group(1)] = int(parts[1], 16)
    assembly = subprocess.check_output(["objdump", "-d", "--no-show-raw-insn", "-C", binary], text=True)
    current = None
    stats = {name: [0, 0] for name in sizes}
    for line in assembly.splitlines():
        if re.match(r"^[0-9a-f]+ <.*>:$", line):
            header = re.search(r"<\(anonymous namespace\)::(Avx(?:2|512)\w+)\(", line)
            current = header.group(1) if header and header.group(1) in stats else None
            continue
        if current and re.match(r"\s*[0-9a-f]+:", line):
            stats[current][0] += 1
            if re.search(r"\bvmov\w*\s+.*\(%r(?:sp|bp)\)|\bvmov\w*\s+.*\(%r(?:sp|bp)\).*,", line):
                stats[current][1] += 1
    print("function,code_bytes,instructions,vector_stack_refs")
    for name in sorted(sizes):
        print(f"{name},{sizes[name]},{stats[name][0]},{stats[name][1]}")


if __name__ == "__main__":
    main()
