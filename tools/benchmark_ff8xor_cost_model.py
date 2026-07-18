#!/usr/bin/env python3
"""Offline AVX2 calibration for the FF8 ISA-aware circuit cost model.

The tool emits and compiles one named-register multiplier kernel per logarithmic
coefficient, measures every kernel in balanced deterministic rounds, and joins
the samples to the repository's assembly-census parser.  It is an offline
developer tool: neither the normal build nor runtime dispatch invokes it.
"""

from __future__ import print_function

import argparse
import json
import math
import os
import platform
import re
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import ff8_xor_cost_model as cost_model  # noqa: E402
import generate_ff8_xor_circuits as circuits  # noqa: E402
import inspect_ff8xor_assembly as inspector  # noqa: E402


NM_CALIBRATION_SYMBOL = re.compile(
    r"^\s*([0-9a-fA-F]+)\s+([0-9a-fA-F]+)\s+[A-Za-z]\s+.*"
    r"avx2::Multiply<(\d+)u?>")
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1


def run(command, cwd=None):
    try:
        return subprocess.check_output(
            command, cwd=str(cwd) if cwd else None,
            stderr=subprocess.STDOUT, universal_newlines=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError("command failed (%s):\n%s" %
                           (" ".join(command), error.output))


def run_stdout(command, cwd=None):
    """Capture a successful program's data stream without diagnostics."""
    process = subprocess.Popen(
        command, cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError("command failed (%s):\n%s%s" %
                           (" ".join(command), stdout, stderr))
    if stderr and not all(line.startswith("sink=")
                          for line in stderr.splitlines()):
        raise RuntimeError("unexpected calibration diagnostic:\n%s" % stderr)
    return stdout


def compiler_identity(compiler):
    first = run([compiler, "--version"]).splitlines()[0].strip()
    return first


def cpu_model():
    try:
        for line in Path("/proc/cpuinfo").read_text(
                encoding="utf-8", errors="replace").splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def generated_source():
    lines = [
        "#include <algorithm>",
        "#include <chrono>",
        "#include <cmath>",
        "#include <cstdint>",
        "#include <cstdlib>",
        "#include <cstdio>",
        "#include <vector>",
        "#include <immintrin.h>",
        "#if defined(__linux__)",
        "# include <sched.h>",
        "#endif",
        "static inline __m256i XorValue(__m256i a, __m256i b)",
        "{ return _mm256_xor_si256(a, b); }",
        "#define LEO_FORCE_INLINE inline __attribute__((always_inline))",
        '#include "generated/LeopardFF8XorCircuits.inl"',
        "#undef LEO_FORCE_INLINE",
        "namespace leopard { namespace ff8xor { namespace avx2 {",
        "template <unsigned Coefficient>",
        "__attribute__((noinline)) void Multiply(",
        "    uint64_t iterations, unsigned char* buffer)",
        "{",
        "    __m256i x0 = _mm256_load_si256((const __m256i*)(buffer + 0 * 32));",
        "    __m256i x1 = _mm256_load_si256((const __m256i*)(buffer + 1 * 32));",
        "    __m256i x2 = _mm256_load_si256((const __m256i*)(buffer + 2 * 32));",
        "    __m256i x3 = _mm256_load_si256((const __m256i*)(buffer + 3 * 32));",
        "    __m256i x4 = _mm256_load_si256((const __m256i*)(buffer + 4 * 32));",
        "    __m256i x5 = _mm256_load_si256((const __m256i*)(buffer + 5 * 32));",
        "    __m256i x6 = _mm256_load_si256((const __m256i*)(buffer + 6 * 32));",
        "    __m256i x7 = _mm256_load_si256((const __m256i*)(buffer + 7 * 32));",
        "    for (uint64_t iteration = 0; iteration < iterations; ++iteration)",
        "    {",
        "        generated::MultiplyCircuit<Coefficient>::Apply(",
        "            x0, x1, x2, x3, x4, x5, x6, x7);",
        "        __asm__ volatile(\"\" : \"+x\"(x0), \"+x\"(x1),",
        "            \"+x\"(x2), \"+x\"(x3), \"+x\"(x4), \"+x\"(x5),",
        "            \"+x\"(x6), \"+x\"(x7));",
        "    }",
        "    _mm256_store_si256((__m256i*)(buffer + 0 * 32), x0);",
        "    _mm256_store_si256((__m256i*)(buffer + 1 * 32), x1);",
        "    _mm256_store_si256((__m256i*)(buffer + 2 * 32), x2);",
        "    _mm256_store_si256((__m256i*)(buffer + 3 * 32), x3);",
        "    _mm256_store_si256((__m256i*)(buffer + 4 * 32), x4);",
        "    _mm256_store_si256((__m256i*)(buffer + 5 * 32), x5);",
        "    _mm256_store_si256((__m256i*)(buffer + 6 * 32), x6);",
        "    _mm256_store_si256((__m256i*)(buffer + 7 * 32), x7);",
        "}",
        "#define LEO_INSTANTIATE(Coefficient) \\",
        "    template void Multiply<Coefficient>(uint64_t, unsigned char*);",
        "LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_INSTANTIATE)",
        "#undef LEO_INSTANTIATE",
        "}}} // namespace leopard::ff8xor::avx2",
        "typedef void (*Kernel)(uint64_t, unsigned char*);",
        "static const Kernel Kernels[256] = {",
        "#define LEO_ENTRY(Coefficient) \\",
        "    &leopard::ff8xor::avx2::Multiply<Coefficient>,",
        "LEOPARD_FF8XOR_FOR_EACH_LOG(LEO_ENTRY)",
        "#undef LEO_ENTRY",
        "};",
        "static volatile uint64_t Sink = 0;",
        "static double Median(std::vector<double> values)",
        "{",
        "    std::sort(values.begin(), values.end());",
        "    const size_t middle = values.size() / 2;",
        "    return values.size() & 1 ? values[middle] :",
        "        (values[middle - 1] + values[middle]) * 0.5;",
        "}",
        "int main(int argc, char** argv)",
        "{",
        "    if (argc != 4) return 2;",
        "    const uint64_t iterations = std::strtoull(argv[1], 0, 10);",
        "    const unsigned rounds = (unsigned)std::strtoul(argv[2], 0, 10);",
        "    const unsigned cpu = (unsigned)std::strtoul(argv[3], 0, 10);",
        "    if (!iterations || rounds < 3) return 2;",
        "#if defined(__linux__)",
        "    if (cpu >= CPU_SETSIZE) return 3;",
        "    cpu_set_t set; CPU_ZERO(&set); CPU_SET(cpu, &set);",
        "    if (sched_setaffinity(0, sizeof(set), &set) != 0) return 3;",
        "#else",
        "    (void)cpu;",
        "#endif",
        "    alignas(64) unsigned char buffer[8 * 32];",
        "    for (unsigned i = 0; i < sizeof(buffer); ++i)",
        "        buffer[i] = (unsigned char)(i * 29u + 17u);",
        "    std::vector<std::vector<double> > samples(256);",
        "    std::vector<unsigned> order(256);",
        "    for (unsigned i = 0; i < 256; ++i) { order[i] = i; Kernels[i](17, buffer); }",
        "    uint64_t random = UINT64_C(0x243f6a8885a308d3);",
        "    for (unsigned round = 0; round < rounds; ++round)",
        "    {",
        "        for (unsigned i = 256; i > 1; --i)",
        "        {",
        "            random ^= random >> 12; random ^= random << 25; random ^= random >> 27;",
        "            const unsigned other = (unsigned)((random * UINT64_C(2685821657736338717)) % i);",
        "            std::swap(order[i - 1], order[other]);",
        "        }",
        "        if (round & 1) std::reverse(order.begin(), order.end());",
        "        for (unsigned position = 0; position < 256; ++position)",
        "        {",
        "            const unsigned coefficient = order[position];",
        "            const std::chrono::steady_clock::time_point begin =",
        "                std::chrono::steady_clock::now();",
        "            Kernels[coefficient](iterations, buffer);",
        "            const std::chrono::steady_clock::time_point end =",
        "                std::chrono::steady_clock::now();",
        "            const double nanoseconds =",
        "                std::chrono::duration<double, std::nano>(end - begin).count();",
        "            samples[coefficient].push_back(nanoseconds / iterations);",
        "            Sink ^= buffer[(coefficient * 37u + round) & 255u];",
        "        }",
        "    }",
        "    std::printf(\"checksum,%s\\n\",",
        "        leopard::ff8xor::generated::kCircuitChecksum);",
        "    std::printf(\"coefficient,median_ns,best_ns,mad_ns\\n\");",
        "    for (unsigned coefficient = 0; coefficient < 256; ++coefficient)",
        "    {",
        "        const double median = Median(samples[coefficient]);",
        "        std::vector<double> deviations;",
        "        for (size_t i = 0; i < samples[coefficient].size(); ++i)",
        "            deviations.push_back(std::abs(samples[coefficient][i] - median));",
        "        std::printf(\"%u,%.9f,%.9f,%.9f\\n\", coefficient, median,",
        "            *std::min_element(samples[coefficient].begin(),",
        "                samples[coefficient].end()), Median(deviations));",
        "    }",
        "    std::fprintf(stderr, \"sink=%llu\\n\",",
        "        (unsigned long long)Sink);",
        "    return 0;",
        "}",
    ]
    return "\n".join(lines) + "\n"


def parse_timings(text):
    # run() deliberately preserves stderr with failed-command diagnostics; the
    # benchmark's volatile sink is written there and may therefore share this
    # captured stream on a successful run.
    lines = [line for line in text.splitlines()
             if not line.startswith("sink=")]
    if len(lines) != 258 or not lines[0].startswith("checksum,") or \
            lines[1] != "coefficient,median_ns,best_ns,mad_ns":
        raise RuntimeError("unexpected calibration benchmark output")
    result = {}
    for line in lines[2:]:
        fields = line.split(",")
        if len(fields) != 4:
            raise RuntimeError("malformed timing row")
        coefficient = int(fields[0])
        if coefficient in result or not (0 <= coefficient < 256):
            raise RuntimeError("duplicate or invalid timing coefficient")
        timing = {
            "median_ns": float(fields[1]),
            "best_ns": float(fields[2]),
            "mad_ns": float(fields[3]),
        }
        if (not math.isfinite(timing["median_ns"]) or
                not math.isfinite(timing["best_ns"]) or
                not math.isfinite(timing["mad_ns"]) or
                timing["median_ns"] <= 0 or timing["best_ns"] <= 0 or
                timing["mad_ns"] < 0):
            raise RuntimeError("invalid calibration timing value")
        result[coefficient] = timing
    if set(result) != set(range(256)):
        raise RuntimeError("calibration output omitted coefficients")
    return lines[0].split(",", 1)[1], result


def symbol_addresses(nm_output):
    result = {}
    for line in nm_output.splitlines():
        match = NM_CALIBRATION_SYMBOL.match(line)
        if not match:
            continue
        coefficient = int(match.group(3))
        candidate = (int(match.group(1), 16), int(match.group(2), 16))
        if coefficient not in result or candidate[1] > result[coefficient][1]:
            result[coefficient] = candidate
    return result


def compile_and_measure(arguments, temporary, circuit_data):
    source = temporary / "ff8xor_cost_calibration.cpp"
    artifact = temporary / "ff8xor_cost_calibration"
    source.write_text(generated_source(), encoding="utf-8")
    flags = (
        "-std=c++11", "-O3", "-DNDEBUG", "-march=x86-64", "-mavx2",
        "-fno-tree-reassoc", "-fno-omit-frame-pointer")
    run([arguments.compiler] + list(flags) + [
        "-I", str(ROOT), str(source), "-o", str(artifact)])
    timing_text = run_stdout([
        str(artifact), str(arguments.iterations), str(arguments.rounds),
        str(arguments.cpu)])
    measured_checksum, timings = parse_timings(timing_text)
    expected_checksum = circuits.circuit_checksum(circuit_data)
    if measured_checksum != expected_checksum:
        raise RuntimeError("compiled circuit checksum does not match generator")

    nm_tool = inspector.find_tool(arguments.nm, ("nm", "llvm-nm"))
    objdump_tool = inspector.find_tool(
        arguments.objdump, ("objdump", "llvm-objdump"))
    nm_output = run([nm_tool, "-S", "-C", "--defined-only", str(artifact)])
    sizes = inspector.parse_sizes(nm_output)
    addresses = symbol_addresses(nm_output)
    disassembly = run([
        objdump_tool, "-drwC", "--no-show-raw-insn", str(artifact)])
    functions = inspector.parse_disassembly(disassembly, sizes)
    by_coefficient = {}
    for function in functions:
        if function["family"] != "multiply_avx2":
            continue
        coefficient = function["coefficient"]
        if (coefficient not in by_coefficient or
                function["code_bytes"] >
                by_coefficient[coefficient]["code_bytes"]):
            by_coefficient[coefficient] = function
    if set(by_coefficient) != set(range(256)) or \
            set(addresses) != set(range(256)):
        raise RuntimeError("assembly census did not find all 256 kernels")

    records = []
    for coefficient, gates in enumerate(circuit_data["multiply_circuits"]):
        function = by_coefficient[coefficient]
        address, symbol_size = addresses[coefficient]
        code_bytes = function["code_bytes"] or symbol_size
        line_offset = address % cost_model.ICACHE_LINE_BYTES
        icache_lines = (line_offset + code_bytes +
                        cost_model.ICACHE_LINE_BYTES - 1) // \
            cost_model.ICACHE_LINE_BYTES
        source_metrics = cost_model.source_metrics(
            gates, circuits.WIRE_COUNT_MULTIPLY)
        function_with_placement = dict(function)
        function_with_placement["icache_lines"] = icache_lines
        compiled_metrics = cost_model.assembly_metrics(
            source_metrics, function_with_placement)
        records.append({
            "coefficient": coefficient,
            "source_gate_count": len(gates),
            "source_depth": circuits.circuit_depth(
                gates, circuits.WIRE_COUNT_MULTIPLY),
            "synthesis_variant": circuit_data["multiply_variants"][coefficient],
            "source_metrics": source_metrics,
            "compiled_metrics": compiled_metrics,
            "symbol_address_mod_64": line_offset,
            "timing": timings[coefficient],
        })

    evidence = {
        "schema": cost_model.CALIBRATION_SCHEMA,
        "circuit_checksum_sha256": expected_checksum,
        "target_isa": "avx2",
        "compiler": compiler_identity(arguments.compiler),
        "compiler_flags": list(flags),
        "cpu_model": cpu_model(),
        "operating_system": "%s %s" % (platform.system(), platform.release()),
        "pinned_cpu": arguments.cpu,
        "measurement_order": (
            "deterministically shuffled coefficient order; alternate reversal"),
        "clock": "std::chrono::steady_clock",
        "rounds": arguments.rounds,
        "circuit_iterations_per_sample": arguments.iterations,
        "kernel_shape": (
            "eight named YMM wires loaded once, repeated literal circuit with "
            "compiler barrier, eight stores"),
        "assembly_tools": {
            "nm": os.path.basename(nm_tool),
            "objdump": os.path.basename(objdump_tool),
        },
        "record_count": len(records),
        "records": records,
    }
    evidence["checksum_sha256"] = cost_model.checksum(evidence)
    return evidence


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--rounds", type=int, default=15)
    parser.add_argument("--iterations", type=int, default=100000)
    parser.add_argument("--nm")
    parser.add_argument("--objdump")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def validate_options(iterations, rounds, cpu):
    if (any(not isinstance(value, int) or isinstance(value, bool)
            for value in (iterations, rounds, cpu)) or
            not 1 <= iterations <= UINT64_MAX or
            not 3 <= rounds <= UINT32_MAX or
            not 0 <= cpu <= UINT32_MAX):
        raise ValueError("invalid calibration count or CPU")


def main():
    arguments = parse_arguments()
    try:
        validate_options(
            arguments.iterations, arguments.rounds, arguments.cpu)
    except ValueError as error:
        raise SystemExit(str(error))
    circuit_data = circuits.build_circuits()
    with tempfile.TemporaryDirectory(prefix="leopard-ff8xor-cost-") as path:
        evidence = compile_and_measure(arguments, Path(path), circuit_data)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(evidence, sort_keys=True, indent=2) + "\n",
        encoding="utf-8")
    nonidentity = [record for record in evidence["records"]
                   if record["coefficient"] not in (0, 255)]
    gates = [record["source_gate_count"] for record in nonidentity]
    timings = [record["timing"]["median_ns"] for record in nonidentity]
    print("Generated %s" % arguments.output)
    print("Calibration checksum: %s" % evidence["checksum_sha256"])
    print("Raw gate-count Spearman: %.6f" %
          cost_model.spearman(gates, timings))
    return 0


if __name__ == "__main__":
    sys.exit(main())
