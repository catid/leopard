#!/usr/bin/env python3
"""Inspect every generated four-buffer specialization in a built archive.

Unlike ``inspect_ff8xor_four_buffer_avx512.py``, which compiles a small source
harness while developing the generated circuit corpus, this tool reads the
actual ``LeopardFF8XorAVX512Four.cpp.o`` member linked into ``libleopard``.  It
checks all 64 tuples in both directions and both XOR2/XOR3 lowerings.

The optimized kernels are allowed to spill scalar plane-address bookkeeping.
They are not allowed to spill a payload vector, call from the offset loop, or
perform table-based field arithmetic.  RBP is only treated as stack storage
when the disassembly establishes it from RSP; Clang otherwise commonly uses
RBP as one of the caller-provided payload bases.
"""

from __future__ import print_function

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "generated" / "FF8XorFourBufferCircuits.json"
OBJECT_BASENAME = "LeopardFF8XorAVX512Four.cpp.o"
SPECIALIZATION_RE = re.compile(
    r"ApplyWhole<(\d+)u?,\s*(true|false),\s*(true|false)>\(")
DISASSEMBLY_HEADER_RE = re.compile(r"^[0-9a-fA-F]+ <(.+)>:$")
INSTRUCTION_RE = re.compile(
    r"^\s*[0-9a-fA-F]+:\s+([a-zA-Z][a-zA-Z0-9_.]*)(?:\s+(.*))?$")
VECTOR_RE = re.compile(r"%(?:zmm|ymm|xmm)\d+", re.IGNORECASE)
RSP_MEMORY_RE = re.compile(r"\([^)]*%rsp(?:,|\))|\[\s*rsp(?:\s*[+\]])",
                           re.IGNORECASE)
RBP_MEMORY_RE = re.compile(r"\([^)]*%rbp(?:,|\))|\[\s*rbp(?:\s*[+\]])",
                           re.IGNORECASE)
RIP_MEMORY_RE = re.compile(r"\([^)]*%rip(?:,|\))|\[\s*rip(?:\s*[+\]])",
                           re.IGNORECASE)


def run_text(command):
    return subprocess.run(
        command, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout


def select_archive_member(member_names):
    """Return the unique production object member from ``ar t`` output."""
    exact = [name for name in member_names
             if Path(name).name == OBJECT_BASENAME]
    if len(exact) == 1:
        return exact[0]
    # A few generators decorate object members while retaining the source
    # stem.  Accept that form only when it is still unambiguous.
    stem = "LeopardFF8XorAVX512Four.cpp"
    fallback = [name for name in member_names
                if stem in Path(name).name and name.endswith((".o", ".obj"))]
    if len(fallback) == 1:
        return fallback[0]
    if not exact and not fallback:
        raise RuntimeError("archive does not contain %s" % OBJECT_BASENAME)
    raise RuntimeError("archive contains ambiguous four-buffer objects: %s" %
                       ", ".join(sorted(exact or fallback)))


def extract_archive_member(ar, archive, destination):
    members = [line.strip() for line in run_text([ar, "t", str(archive)]).
               splitlines() if line.strip()]
    member = select_archive_member(members)
    result = subprocess.run(
        [ar, "p", str(archive), member], check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    destination.write_bytes(result.stdout)
    return member


def parse_nm_sizes(output):
    """Parse ``nm -S -C --defined-only`` into demangled-name sizes."""
    result = {}
    for line in output.splitlines():
        parts = line.strip().split(None, 3)
        if len(parts) != 4:
            continue
        try:
            size = int(parts[1], 16)
        except ValueError:
            continue
        result[parts[3]] = size
    return result


def parse_disassembly(output):
    """Return demangled symbol bodies without assuming templates lack ``>``."""
    result = {}
    current = None
    lines = []
    for line in output.splitlines():
        match = DISASSEMBLY_HEADER_RE.match(line.strip())
        if match:
            if current is not None:
                result[current] = "\n".join(lines)
            current = match.group(1)
            lines = []
        elif current is not None:
            lines.append(line)
    if current is not None:
        result[current] = "\n".join(lines)
    return result


def instruction_records(body):
    records = []
    for line in body.splitlines():
        match = INSTRUCTION_RE.match(line)
        if match:
            records.append({
                "line": line.strip(),
                "mnemonic": match.group(1).lower(),
                "operands": (match.group(2) or "").strip(),
            })
    return records


def establishes_rbp_frame(records):
    for record in records[:24]:
        if not record["mnemonic"].startswith("mov"):
            continue
        operands = re.sub(r"\s+", "", record["operands"].lower())
        # GNU/LLVM objdump's default AT&T form, followed by Intel form for
        # parser tests and downstream override tools.
        if operands.startswith("%rsp,%rbp") or operands.startswith("rbp,rsp"):
            return True
    return False


def split_operands(operands):
    """Split assembler operands while retaining commas inside addresses."""
    result = []
    start = 0
    depth = 0
    for index, character in enumerate(operands):
        if character in "([":
            depth += 1
        elif character in ")]":
            depth = max(0, depth - 1)
        elif character == "," and depth == 0:
            result.append(operands[start:index].strip())
            start = index + 1
    result.append(operands[start:].strip())
    return result


def is_memory_operand(operand):
    return "(" in operand or "[" in operand


def ternary_immediate_is_xor3(operands):
    """Recognize the actual AT&T/Intel immediate, not an address displacement."""
    immediates = []
    for operand in split_operands(operands):
        text = operand.strip().lower()
        if text.startswith("$"):
            text = text[1:]
        if not re.fullmatch(r"(?:0x[0-9a-f]+|[0-9]+)", text):
            continue
        try:
            immediates.append(int(text, 0))
        except ValueError:
            return False
    return immediates == [0x96]


def memory_vector_direction(record):
    """Classify a vector move as payload load/store, or return ``None``."""
    mnemonic = record["mnemonic"]
    operands = record["operands"]
    if not mnemonic.startswith("vmov") or not VECTOR_RE.search(operands):
        return None
    pieces = split_operands(operands)
    if len(pieces) < 2:
        return None
    first_memory = "(" in pieces[0] or "[" in pieces[0]
    last_memory = "(" in pieces[-1] or "[" in pieces[-1]
    if first_memory and not last_memory:
        return "load"
    if last_memory and not first_memory:
        return "store"
    return None


def inspect_body(body):
    records = instruction_records(body)
    rbp_frame = establishes_rbp_frame(records)
    stack_vectors = []
    payload_rbp_vectors = []
    calls = []
    forbidden = []
    non_xor3_ternary = []
    non_move_vector_memory = []
    rip_relative = []
    loads = 0
    stores = 0
    mnemonics = []
    for record in records:
        line = record["line"]
        mnemonic = record["mnemonic"]
        operands = record["operands"]
        mnemonics.append(mnemonic)
        has_vector = VECTOR_RE.search(operands) is not None
        has_memory = any(
            is_memory_operand(operand) for operand in split_operands(operands))
        if mnemonic.startswith("call"):
            calls.append(line)
        if has_vector and RSP_MEMORY_RE.search(operands):
            stack_vectors.append(line)
        if has_vector and RBP_MEMORY_RE.search(operands):
            if rbp_frame:
                stack_vectors.append(line)
            else:
                payload_rbp_vectors.append(line)
        if RIP_MEMORY_RE.search(operands):
            rip_relative.append(line)
        if ("pshufb" in mnemonic or "gf2p8" in mnemonic or
                "pclmul" in mnemonic or mnemonic.startswith("vpmull") or
                mnemonic.startswith(("vpand", "pand")) or
                mnemonic in ("andps", "andpd", "andnps", "andnpd",
                             "vandps", "vandpd")):
            forbidden.append(line)
        if (mnemonic.startswith("vpternlog") and
                not ternary_immediate_is_xor3(operands)):
            non_xor3_ternary.append(line)
        if has_vector and has_memory and not mnemonic.startswith("vmov"):
            non_move_vector_memory.append(line)
        direction = memory_vector_direction(record)
        if direction == "load":
            loads += 1
        elif direction == "store":
            stores += 1
    return {
        "instruction_count": len(records),
        "vector_xor_count": sum(
            mnemonic.startswith("vpxor") for mnemonic in mnemonics),
        "vpternlog_count": sum(
            mnemonic.startswith("vpternlog") for mnemonic in mnemonics),
        "vector_move_count": sum(
            mnemonic.startswith("vmov") for mnemonic in mnemonics),
        "payload_vector_load_count": loads,
        "payload_vector_store_count": stores,
        "call_count": len(calls),
        "calls": calls,
        "rbp_is_frame_pointer": rbp_frame,
        "vector_stack_reference_count": len(stack_vectors),
        "vector_stack_references": stack_vectors,
        "payload_rbp_vector_reference_count": len(payload_rbp_vectors),
        "forbidden_field_instruction_count": len(forbidden),
        "forbidden_field_instructions": forbidden,
        "non_xor3_ternary_count": len(non_xor3_ternary),
        "non_xor3_ternary_instructions": non_xor3_ternary,
        "non_move_vector_memory_count": len(non_move_vector_memory),
        "non_move_vector_memory_instructions": non_move_vector_memory,
        "rip_relative_memory_reference_count": len(rip_relative),
        "rip_relative_memory_references": rip_relative,
        "instruction_shape": sorted(set(mnemonics)),
    }


def specialization_key(symbol):
    match = SPECIALIZATION_RE.search(symbol)
    if not match:
        return None
    return (int(match.group(1)), match.group(2) == "true",
            match.group(3) == "true")


def expected_operation_counts(metadata, key):
    tuple_index, inverse, use_xor3 = key
    family = "ifft" if inverse else "fft"
    record = metadata["families"][family]["records"][tuple_index]
    selected = record["xor3_selected" if use_xor3 else "xor2_selected"]
    schedule = selected["xor3" if use_xor3 else "xor2"]
    return schedule["xor2_count"], schedule["xor3_count"]


def inspect_object(object_path, nm, objdump, metadata):
    nm_output = run_text([
        nm, "-S", "-C", "--defined-only", str(object_path)])
    sizes = parse_nm_sizes(nm_output)
    disassembly = run_text([
        objdump, "-d", "-C", "--no-show-raw-insn", str(object_path)])
    bodies = parse_disassembly(disassembly)
    kernels = []
    seen = set()
    failures = []
    for symbol in sorted(bodies):
        key = specialization_key(symbol)
        if key is None:
            continue
        if key in seen:
            failures.append("duplicate ApplyWhole specialization %s" % (key,))
            continue
        seen.add(key)
        observation = {
            "symbol": symbol,
            "tuple_index": key[0],
            "family": "ifft" if key[1] else "fft",
            "lowering": "xor3" if key[2] else "xor2",
            "code_bytes": sizes.get(symbol, 0),
        }
        observation.update(inspect_body(bodies[symbol]))
        expected_xor2, expected_xor3 = expected_operation_counts(metadata, key)
        observation["expected_vector_xor_count"] = expected_xor2
        observation["expected_vpternlog_count"] = expected_xor3
        kernels.append(observation)

    expected = set((index, inverse, xor3)
                   for index in range(64)
                   for inverse in (False, True)
                   for xor3 in (False, True))
    for key in sorted(expected - seen):
        failures.append("missing ApplyWhole specialization %s" % (key,))
    for key in sorted(seen - expected):
        failures.append("unexpected ApplyWhole specialization %s" % (key,))
    return kernels, failures


def validate_kernels(kernels, initial_failures):
    failures = list(initial_failures)
    for kernel in kernels:
        label = "%s/%s/tuple%d" % (
            kernel["family"], kernel["lowering"], kernel["tuple_index"])
        if not kernel["code_bytes"]:
            failures.append("%s has no nm symbol size" % label)
        if kernel["call_count"]:
            failures.append("%s contains calls" % label)
        if kernel["vector_stack_reference_count"]:
            failures.append("%s contains payload-vector stack references" % label)
        if kernel["forbidden_field_instruction_count"]:
            failures.append("%s contains forbidden field instructions" % label)
        if kernel["non_xor3_ternary_count"]:
            failures.append("%s contains non-XOR3 VPTERNLOG immediates" % label)
        if kernel["non_move_vector_memory_count"]:
            failures.append("%s contains non-move vector memory operands" % label)
        if kernel["rip_relative_memory_reference_count"]:
            failures.append("%s contains RIP-relative memory references" % label)
        if kernel["payload_vector_load_count"] != 32:
            failures.append("%s has %d payload loads, expected 32" %
                            (label, kernel["payload_vector_load_count"]))
        # A compiler may prove that sentinel-map outputs are unchanged and
        # omit their write-back.  More than 32 stores is always redundant;
        # fewer stores remain covered by exhaustive runtime map equivalence.
        if kernel["payload_vector_store_count"] > 32:
            failures.append("%s has %d payload stores, expected at most 32" %
                            (label, kernel["payload_vector_store_count"]))
        if kernel["vector_xor_count"] != \
                kernel["expected_vector_xor_count"]:
            failures.append("%s emitted %d vector XORs, expected %d" % (
                label, kernel["vector_xor_count"],
                kernel["expected_vector_xor_count"]))
        if kernel["vpternlog_count"] != kernel["expected_vpternlog_count"]:
            failures.append("%s emitted %d VPTERNLOGs, expected %d" % (
                label, kernel["vpternlog_count"],
                kernel["expected_vpternlog_count"]))
    return failures


def summarize(kernels):
    result = {}
    for family in ("fft", "ifft"):
        for lowering in ("xor2", "xor3"):
            selected = [kernel for kernel in kernels
                        if kernel["family"] == family and
                        kernel["lowering"] == lowering]
            key = family + "_" + lowering
            code = [kernel["code_bytes"] for kernel in selected]
            instructions = [kernel["instruction_count"] for kernel in selected]
            moves = [kernel["vector_move_count"] for kernel in selected]
            loads = [kernel["payload_vector_load_count"] for kernel in selected]
            stores = [kernel["payload_vector_store_count"] for kernel in selected]
            result[key] = {
                "specialization_count": len(selected),
                "code_bytes_min": min(code) if code else 0,
                "code_bytes_max": max(code) if code else 0,
                "code_bytes_average": (sum(code) / float(len(code)))
                    if code else 0.0,
                "instruction_count_min": min(instructions) if instructions else 0,
                "instruction_count_max": max(instructions) if instructions else 0,
                "instruction_count_average": (
                    sum(instructions) / float(len(instructions)))
                    if instructions else 0.0,
                "vector_move_count_min": min(moves) if moves else 0,
                "vector_move_count_max": max(moves) if moves else 0,
                "payload_vector_load_count_min": min(loads) if loads else 0,
                "payload_vector_load_count_max": max(loads) if loads else 0,
                "payload_vector_store_count_min": min(stores) if stores else 0,
                "payload_vector_store_count_max": max(stores) if stores else 0,
                "excess_payload_reload_count": sum(
                    max(0, count - 32) for count in loads),
                "clean_memory_traffic_specialization_count": sum(
                    load == 32 and store <= 32
                    for load, store in zip(loads, stores)),
            }
    return result


def select_observations(kernels, metadata):
    by_key = {(kernel["family"], kernel["lowering"],
               kernel["tuple_index"]): kernel for kernel in kernels}
    result = []
    for family in ("fft", "ifft"):
        records = metadata["families"][family]["records"]
        common = max(
            (record for record in records if 255 not in record["tuple"]),
            key=lambda record: (
                record["parameter_census_calls"], -record["tuple_index"]))
        for lowering in ("xor2", "xor3"):
            schedule_name = lowering
            stress = max(records, key=lambda record: (
                record[lowering + "_selected"][schedule_name][
                    "instruction_count"], -record["tuple_index"]))
            for selection, record in (("frequency-representative", common),
                                      ("max-instruction-stress", stress)):
                kernel = by_key.get((family, lowering, record["tuple_index"]))
                if kernel:
                    result.append({
                        "selection": selection,
                        "family": family,
                        "lowering": lowering,
                        "tuple_index": record["tuple_index"],
                        "tuple": record["tuple"],
                        "parameter_census_calls": record[
                            "parameter_census_calls"],
                        "code_bytes": kernel["code_bytes"],
                        "instruction_count": kernel["instruction_count"],
                        "vector_move_count": kernel["vector_move_count"],
                        "payload_vector_load_count": kernel[
                            "payload_vector_load_count"],
                        "payload_vector_store_count": kernel[
                            "payload_vector_store_count"],
                        "excess_payload_reload_count": max(
                            0, kernel["payload_vector_load_count"] - 32),
                        "vector_stack_reference_count": kernel[
                            "vector_stack_reference_count"],
                        "call_count": kernel["call_count"],
                    })
    return result


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True,
                        help="CMake-built libleopard static archive")
    parser.add_argument("--metadata", default=str(DEFAULT_METADATA))
    parser.add_argument("--ar")
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    archive = Path(arguments.archive).resolve()
    ar = arguments.ar or shutil.which("llvm-ar") or shutil.which("ar")
    nm = arguments.nm or shutil.which("llvm-nm") or shutil.which("nm")
    objdump = (arguments.objdump or shutil.which("llvm-objdump") or
               shutil.which("objdump"))
    if not archive.is_file():
        print("archive does not exist: %s" % archive, file=sys.stderr)
        return 2
    if not ar or not nm or not objdump:
        print("ar, nm, and objdump are required", file=sys.stderr)
        return 2
    metadata = json.loads(Path(arguments.metadata).read_text(encoding="utf-8"))

    with tempfile.TemporaryDirectory(
            prefix="ff8xor-four-production-") as directory:
        object_path = Path(directory) / OBJECT_BASENAME
        member = extract_archive_member(ar, archive, object_path)
        kernels, parse_failures = inspect_object(
            object_path, nm, objdump, metadata)
        object_bytes = object_path.read_bytes()

    failures = validate_kernels(kernels, parse_failures)
    excess_payload_reloads = sum(
        max(0, kernel["payload_vector_load_count"] - 32)
        for kernel in kernels)
    report = {
        "archive": str(archive),
        "archive_member": member,
        "object_bytes": len(object_bytes),
        "object_sha256": hashlib.sha256(object_bytes).hexdigest(),
        "specialization_count": len(kernels),
        "status": "clean" if not failures else "failed",
        "excess_payload_reload_count": excess_payload_reloads,
        "summary": summarize(kernels),
        "representative_and_stress": select_observations(kernels, metadata),
        "failures": failures,
    }
    if arguments.json:
        print(json.dumps(report, sort_keys=True, indent=2))
    else:
        print("archive member: %s (%d bytes, sha256=%s)" % (
            member, report["object_bytes"], report["object_sha256"]))
        print("ApplyWhole specializations: %d" % len(kernels))
        for key in sorted(report["summary"]):
            item = report["summary"][key]
            print("  %s: count=%d code=%d..%dB avg=%.1fB "
                  "insns=%d..%d avg=%.1f moves=%d..%d "
                  "loads=%d..%d stores=%d..%d excess_reloads=%d" % (
                      key, item["specialization_count"],
                      item["code_bytes_min"], item["code_bytes_max"],
                      item["code_bytes_average"],
                      item["instruction_count_min"],
                      item["instruction_count_max"],
                      item["instruction_count_average"],
                      item["vector_move_count_min"],
                      item["vector_move_count_max"],
                      item["payload_vector_load_count_min"],
                      item["payload_vector_load_count_max"],
                      item["payload_vector_store_count_min"],
                      item["payload_vector_store_count_max"],
                      item["excess_payload_reload_count"]))
        for item in report["representative_and_stress"]:
            print("  %s %s %s tuple=%d %s: code=%dB insns=%d moves=%d "
                  "loads=%d stores=%d excess_reloads=%d spills=%d calls=%d" % (
                      item["family"], item["lowering"], item["selection"],
                      item["tuple_index"], item["tuple"], item["code_bytes"],
                      item["instruction_count"], item["vector_move_count"],
                      item["payload_vector_load_count"],
                      item["payload_vector_store_count"],
                      item["excess_payload_reload_count"],
                      item["vector_stack_reference_count"],
                      item["call_count"]))
        for failure in failures:
            print("FAIL: %s" % failure, file=sys.stderr)
    return 1 if arguments.strict and failures else 0


if __name__ == "__main__":
    sys.exit(main())
