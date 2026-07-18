#!/usr/bin/env python3
"""Compile and inspect representative generated AVX-512 radix-4 kernels.

The 32-wire maps exactly fill the architectural ZMM register file.  Ordinary
intrinsic wrappers let GCC/Clang reassociate the circuit DAG and were observed
to spill.  This inspection intentionally uses destructive ``+v`` inline-asm
XOR wrappers, the lowering required by a future production integration.
"""

from __future__ import print_function

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "generated" / "FF8XorFourBufferCircuits.json"


def select_representatives(metadata):
    result = []
    for family in ("fft", "ifft"):
        records = metadata["families"][family]["records"]
        # The common case is selected by the exhaustive API-shape census, but
        # excludes the unusually cheap sentinel tuple.  Retain the maximum-op
        # stress case too so a clean representative cannot conceal pressure at
        # the worst edge of the corpus.
        common = max(
            (record for record in records if 255 not in record["tuple"]),
            key=lambda record: (
                record["parameter_census_calls"], -record["tuple_index"]))
        for lowering in ("xor2", "xor3"):
            field = lowering
            stress = max(
                records,
                key=lambda record: (
                    record[lowering + "_selected"][field][
                        "instruction_count"],
                    -record["tuple_index"]))
            for selection, record in (
                    ("frequency-representative", common),
                    ("max-instruction-stress", stress)):
                result.append({
                    "symbol": "%s4_%s_%s" % (
                        family, lowering,
                        "representative" if selection.startswith("frequency")
                        else "stress"),
                    "selection": selection,
                    "family": family,
                    "lowering": lowering,
                    "tuple_index": record["tuple_index"],
                    "tuple": record["tuple"],
                    "parameter_census_calls": record[
                        "parameter_census_calls"],
                    "gate_count": record[lowering + "_selected"][
                        "gate_count"],
                    "source_instruction_count": record[
                        lowering + "_selected"][field]["instruction_count"],
                })
    return result


def emit_kernel(record):
    prefix = "FFT4" if record["family"] == "fft" else "IFFT4"
    circuit = "%sCircuitXor%s<%d>" % (
        prefix, "2" if record["lowering"] == "xor2" else "3",
        record["tuple_index"])
    names = ["%s%d" % (letter, bit)
             for letter in ("a", "b", "c", "d") for bit in range(8)]
    pointer_for_letter = {"a": "buffer0", "b": "buffer1",
                          "c": "buffer2", "d": "buffer3"}
    lines = [
        "extern \"C\" __attribute__((noinline)) void %s(" % record["symbol"],
        "    unsigned char* __restrict buffer0,",
        "    unsigned char* __restrict buffer1,",
        "    unsigned char* __restrict buffer2,",
        "    unsigned char* __restrict buffer3,",
        "    unsigned long long plane_bytes)",
        "{",
        "    for (unsigned long long offset = 0; offset + 64 <= plane_bytes; offset += 64)",
        "    {",
    ]
    for name in names:
        pointer = pointer_for_letter[name[0]]
        plane = int(name[1])
        lines.append(
            "        __m512i %s = _mm512_loadu_si512((const void*)(%s + %d * plane_bytes + offset));" %
            (name, pointer, plane))
    lines.append("        %s::Apply(" % circuit)
    for index, name in enumerate(names):
        suffix = "," if index + 1 != len(names) else ");"
        lines.append("            %s%s" % (name, suffix))
    for name in names:
        pointer = pointer_for_letter[name[0]]
        plane = int(name[1])
        lines.append(
            "        _mm512_storeu_si512((void*)(%s + %d * plane_bytes + offset), %s);" %
            (pointer, plane, name))
    lines.extend(["    }", "}", ""])
    return "\n".join(lines)


def emit_source(representatives):
    lines = [
        "#include <immintrin.h>",
        "#define LEO_FORCE_INLINE inline __attribute__((always_inline))",
        "namespace leopard { namespace ff8xor {",
        "static LEO_FORCE_INLINE __m512i XorValue(__m512i a, __m512i b)",
        "{ __asm__(\"vpxord %1, %0, %0\" : \"+v\"(a) : \"v\"(b)); return a; }",
        "static LEO_FORCE_INLINE __m512i Xor3Value(__m512i a, __m512i b, __m512i c)",
        "{ __asm__(\"vpternlogd $0x96, %2, %1, %0\" : \"+v\"(a) : \"v\"(b), \"v\"(c)); return a; }",
        "}}",
        "#include \"generated/LeopardFF8XorFourBufferCircuits.inl\"",
        "using namespace leopard::ff8xor::generated4;",
        "",
    ]
    for record in representatives:
        lines.append(emit_kernel(record))
    return "\n".join(lines)


def compiler_version(compiler):
    result = subprocess.run(
        [compiler, "--version"], check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return result.stdout.splitlines()[0]


def compiler_supports_flag(compiler, flag):
    result = subprocess.run(
        [compiler, "-Werror", flag, "-x", "c++", "-fsyntax-only", "-"],
        input="int main() { return 0; }\n", text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode == 0


def symbol_size(nm, object_path, symbol):
    result = subprocess.run(
        [nm, "-S", "--defined-only", str(object_path)], check=True,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pattern = re.compile(
        r"^[0-9a-fA-F]+\s+([0-9a-fA-F]+)\s+\S\s+" +
        re.escape(symbol) + r"$", re.MULTILINE)
    match = pattern.search(result.stdout)
    if not match:
        raise RuntimeError("cannot find symbol size for %s" % symbol)
    return int(match.group(1), 16)


def extract_symbol(disassembly, symbol):
    marker = re.search(
        r"^[0-9a-fA-F]+ <" + re.escape(symbol) + r">:\n",
        disassembly, re.MULTILINE)
    if not marker:
        raise RuntimeError("cannot find disassembly for %s" % symbol)
    following = disassembly[marker.end():]
    next_symbol = re.search(r"^[0-9a-fA-F]+ <[^>]+>:\n", following,
                            re.MULTILINE)
    return following[:next_symbol.start()] if next_symbol else following


def inspect_body(body):
    instruction_lines = []
    mnemonics = []
    for line in body.splitlines():
        match = re.match(r"^\s*[0-9a-fA-F]+:\s+(\S+)(?:\s+(.*))?$", line)
        if not match:
            continue
        instruction_lines.append(line.strip())
        mnemonics.append(match.group(1).lower())
    stack_vector_references = [
        line for line in instruction_lines
        if re.search(r"%(?:zmm|ymm|xmm)\d+", line) and
        # Both compilers omit the frame pointer at -O3 and freely use RBP as
        # one of the four payload-plane address bases.  Only RSP is therefore
        # unambiguously stack storage in this deliberately fixed command.
        re.search(r"\(%rsp", line)]
    calls = [line for line, mnemonic in zip(instruction_lines, mnemonics)
             if mnemonic.startswith("call")]
    forbidden = [
        line for line, mnemonic in zip(instruction_lines, mnemonics)
        if ("pshufb" in mnemonic or "gf2p8" in mnemonic or
            "pclmul" in mnemonic)]
    return {
        "instruction_count": len(instruction_lines),
        "vector_xor_count": sum(
            mnemonic.startswith("vpxor") for mnemonic in mnemonics),
        "vpternlog_count": sum(
            mnemonic.startswith("vpternlog") for mnemonic in mnemonics),
        "vector_move_count": sum(
            mnemonic.startswith("vmov") for mnemonic in mnemonics),
        "call_count": len(calls),
        "calls": calls,
        "vector_stack_reference_count": len(stack_vector_references),
        "vector_stack_references": stack_vector_references,
        "forbidden_field_instruction_count": len(forbidden),
        "forbidden_field_instructions": forbidden,
        "instruction_shape": sorted(set(mnemonics)),
    }


def inspect_compiler(
        compiler, source, representatives, temporary_directory, nm,
        objdump):
    compiler_name = Path(compiler).name
    source_path = temporary_directory / (compiler_name + ".cpp")
    object_path = temporary_directory / (compiler_name + ".o")
    source_path.write_text(source, encoding="utf-8")
    command = [
        compiler, "-O3", "-std=c++11", "-fno-exceptions", "-fno-rtti",
        "-fno-stack-protector", "-march=x86-64", "-mavx2", "-mavx512f",
        "-mavx512vl",
    ]
    # This matches CMake's production FF8XOR_AVX512_COMPILE_FLAGS contract:
    # GCC gets the supported reassociation guard while Clang does not accept
    # that GCC-specific switch.
    if compiler_supports_flag(compiler, "-fno-tree-reassoc"):
        command.append("-fno-tree-reassoc")
    command.extend([
        "-I", str(ROOT), "-c", str(source_path), "-o", str(object_path),
    ])
    subprocess.run(command, check=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, text=True)
    disassembly = subprocess.run(
        [objdump, "-d", "--no-show-raw-insn", str(object_path)], check=True,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    kernels = []
    for record in representatives:
        body = extract_symbol(disassembly, record["symbol"])
        observation = dict(record)
        observation["code_bytes"] = symbol_size(
            nm, object_path, record["symbol"])
        observation.update(inspect_body(body))
        kernels.append(observation)
    return {
        "compiler": os.path.realpath(compiler),
        "compiler_version": compiler_version(compiler),
        "command": command,
        "kernels": kernels,
    }


def validate_strict(results):
    failures = []
    for compiler in results:
        for kernel in compiler["kernels"]:
            label = "%s:%s" % (Path(compiler["compiler"]).name,
                                kernel["symbol"])
            if kernel["call_count"]:
                failures.append("%s contains calls" % label)
            if kernel["vector_stack_reference_count"]:
                failures.append("%s contains vector stack spills" % label)
            if kernel["forbidden_field_instruction_count"]:
                failures.append("%s contains forbidden field instructions" % label)
            if kernel["lowering"] == "xor3" and not kernel["vpternlog_count"]:
                failures.append("%s did not emit VPTERNLOG" % label)
            if kernel["lowering"] == "xor2" and not kernel["vector_xor_count"]:
                failures.append("%s did not emit vector XOR" % label)
    return failures


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", action="append", default=[])
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))
    representatives = select_representatives(metadata)
    source = emit_source(representatives)
    compilers = arguments.compiler
    if not compilers:
        gcc = shutil.which("g++")
        clang = (shutil.which("clang++") or shutil.which("clang++-18") or
                 shutil.which("clang++-17") or shutil.which("clang++-16"))
        compilers = [compiler for compiler in (gcc, clang) if compiler]
    if not compilers:
        print("no GCC or Clang C++ compiler found", file=sys.stderr)
        return 2
    nm = arguments.nm or shutil.which("llvm-nm") or shutil.which("nm")
    objdump = (arguments.objdump or shutil.which("llvm-objdump") or
               shutil.which("objdump"))
    if not nm or not objdump:
        print("nm and objdump (GNU or LLVM) are required", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory(prefix="ff8xor-four-avx512-") as directory:
        temporary_directory = Path(directory)
        results = [inspect_compiler(
            compiler, source, representatives, temporary_directory,
            nm, objdump)
            for compiler in compilers]

    failures = validate_strict(results)
    if arguments.json:
        print(json.dumps({"results": results, "failures": failures},
                         sort_keys=True, indent=2))
    else:
        for compiler in results:
            print(compiler["compiler_version"])
            for kernel in compiler["kernels"]:
                print("  %s tuple=%s gates=%d source_ops=%d code=%dB "
                      "insns=%d xor=%d ternlog=%d moves=%d spills=%d calls=%d" % (
                          kernel["symbol"], kernel["tuple"],
                          kernel["gate_count"],
                          kernel["source_instruction_count"],
                          kernel["code_bytes"], kernel["instruction_count"],
                          kernel["vector_xor_count"], kernel["vpternlog_count"],
                          kernel["vector_move_count"],
                          kernel["vector_stack_reference_count"],
                          kernel["call_count"]))
        for failure in failures:
            print("FAIL: %s" % failure, file=sys.stderr)
    return 1 if arguments.strict and failures else 0


if __name__ == "__main__":
    sys.exit(main())
