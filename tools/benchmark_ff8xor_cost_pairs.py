#!/usr/bin/env python3
"""Same-process AVX2 timing of equivalent portable/machine FF8 schedules."""

from __future__ import print_function

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import evaluate_ff8xor_cost_candidate as evaluation  # noqa: E402
import ff8_xor_cost_model as cost_model  # noqa: E402
from benchmark_ff8xor_cost_model import cpu_model  # noqa: E402


UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1


def run(command):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError("command failed (%s):\n%s%s" %
                           (" ".join(command), stdout, stderr))
    return stdout, stderr


def wire_name(wire):
    return ("x%d" % wire) if wire < 8 else ("y%d" % (wire - 8))


def function_name(profile, family, coefficient):
    return "FF8Cost_%s_%s_%d" % (profile, family, coefficient)


def emit_kernel(lines, profile, family, coefficient, gates):
    width = 8 if family == "multiply" else 16
    name = function_name(profile, family, coefficient)
    lines.extend([
        'extern "C" __attribute__((noinline,noclone)) void %s(' % name,
        "    uint64_t iterations, unsigned char* x_buffer,",
        "    unsigned char* y_buffer)",
        "{",
    ])
    for wire in range(8):
        lines.append(
            "    __m256i x%d = _mm256_load_si256(" % wire +
            "(const __m256i*)(x_buffer + %d * 32));" % wire)
    if width == 16:
        for wire in range(8):
            lines.append(
                "    __m256i y%d = _mm256_load_si256(" % wire +
                "(const __m256i*)(y_buffer + %d * 32));" % wire)
    else:
        lines.append("    (void)y_buffer;")
    lines.extend([
        "    for (uint64_t iteration = 0; iteration < iterations; ++iteration)",
        "    {",
    ])
    for destination, source in gates:
        destination_name = wire_name(destination)
        source_name = wire_name(source)
        lines.append("        %s = _mm256_xor_si256(%s, %s);" % (
            destination_name, destination_name, source_name))
    # Split sixteen read/write operands across two barriers: GCC counts a '+'
    # operand twice against its thirty-operand limit.  Both barriers remain in
    # each loop iteration, preventing cross-iteration algebraic folding.
    lines.append(
        '        __asm__ volatile("" : "+x"(x0), "+x"(x1), "+x"(x2), '
        '"+x"(x3), "+x"(x4), "+x"(x5), "+x"(x6), "+x"(x7));')
    if width == 16:
        lines.append(
            '        __asm__ volatile("" : "+x"(y0), "+x"(y1), "+x"(y2), '
            '"+x"(y3), "+x"(y4), "+x"(y5), "+x"(y6), "+x"(y7));')
    lines.append("    }")
    for wire in range(8):
        lines.append(
            "    _mm256_store_si256((__m256i*)(x_buffer + %d * 32), x%d);" %
            (wire, wire))
    if width == 16:
        for wire in range(8):
            lines.append(
                "    _mm256_store_si256((__m256i*)(y_buffer + %d * 32), y%d);" %
                (wire, wire))
    lines.extend(["}", ""])


def generated_source(portable, machine):
    pairs = []
    lines = [
        "#include <algorithm>",
        "#include <chrono>",
        "#include <cmath>",
        "#include <cstdint>",
        "#include <cstdio>",
        "#include <cstdlib>",
        "#include <cstring>",
        "#include <vector>",
        "#include <immintrin.h>",
        "#if defined(__linux__)",
        "# include <sched.h>",
        "#endif",
        "",
    ]
    for family in evaluation.FAMILIES:
        portable_circuits = portable[family]
        machine_circuits = machine[family]
        for coefficient in range(256):
            if portable_circuits[coefficient] == machine_circuits[coefficient]:
                continue
            emit_kernel(lines, "portable", family, coefficient,
                        portable_circuits[coefficient])
            emit_kernel(lines, "machine", family, coefficient,
                        machine_circuits[coefficient])
            pairs.append((family, coefficient))

    lines.extend([
        "typedef void (*Kernel)(uint64_t, unsigned char*, unsigned char*);",
        "struct Pair { const char* family; unsigned coefficient; Kernel portable; Kernel machine; };",
        "static const Pair Pairs[] = {",
    ])
    for family, coefficient in pairs:
        lines.append('    { "%s", %d, &%s, &%s },' % (
            family, coefficient,
            function_name("portable", family, coefficient),
            function_name("machine", family, coefficient)))
    lines.extend([
        "};",
        "static volatile uint64_t Sink = 0;",
        "static double Median(std::vector<double> values)",
        "{",
        "    std::sort(values.begin(), values.end());",
        "    const size_t middle = values.size() / 2;",
        "    return values.size() & 1 ? values[middle] :",
        "        (values[middle - 1] + values[middle]) * 0.5;",
        "}",
        "static double Time(Kernel kernel, uint64_t iterations,",
        "                   unsigned char* x, unsigned char* y)",
        "{",
        "    const std::chrono::steady_clock::time_point begin =",
        "        std::chrono::steady_clock::now();",
        "    kernel(iterations, x, y);",
        "    const std::chrono::steady_clock::time_point end =",
        "        std::chrono::steady_clock::now();",
        "    Sink ^= x[17] ^ y[29];",
        "    return std::chrono::duration<double, std::nano>(end - begin).count() / iterations;",
        "}",
        "int main(int argc, char** argv)",
        "{",
        "    if (argc != 4) return 2;",
        "    const uint64_t base_iterations = std::strtoull(argv[1], 0, 10);",
        "    const unsigned rounds = (unsigned)std::strtoul(argv[2], 0, 10);",
        "    const unsigned cpu = (unsigned)std::strtoul(argv[3], 0, 10);",
        "#if defined(__linux__)",
        "    if (cpu >= CPU_SETSIZE) return 3;",
        "    cpu_set_t set; CPU_ZERO(&set); CPU_SET(cpu, &set);",
        "    if (sched_setaffinity(0, sizeof(set), &set) != 0) return 3;",
        "#else",
        "    (void)cpu;",
        "#endif",
        "    alignas(64) unsigned char px[256], py[256], mx[256], my[256];",
        '    std::printf("family,coefficient,portable_median_ns,machine_median_ns,portable_best_ns,machine_best_ns,portable_mad_ns,machine_mad_ns\\n");',
        "    for (size_t pair_index = 0; pair_index < sizeof(Pairs) / sizeof(Pairs[0]); ++pair_index)",
        "    {",
        "        const Pair& pair = Pairs[pair_index];",
        "        for (unsigned i = 0; i < 256; ++i)",
        "            px[i] = py[i] = mx[i] = my[i] = (unsigned char)(i * 29u + pair.coefficient * 7u + pair_index);",
        "        pair.portable(1, px, py); pair.machine(1, mx, my);",
        "        if (std::memcmp(px, mx, 256) || std::memcmp(py, my, 256)) return 4;",
        "        const uint64_t iterations = pair.family[0] == 'm' ? base_iterations : base_iterations / 2;",
        "        std::vector<double> portable_samples, machine_samples;",
        "        pair.portable(31, px, py); pair.machine(31, mx, my);",
        "        for (unsigned round = 0; round < rounds; ++round)",
        "        {",
        "            portable_samples.push_back(Time(pair.portable, iterations, px, py));",
        "            machine_samples.push_back(Time(pair.machine, iterations, mx, my));",
        "            machine_samples.push_back(Time(pair.machine, iterations, mx, my));",
        "            portable_samples.push_back(Time(pair.portable, iterations, px, py));",
        "        }",
        "        const double portable_median = Median(portable_samples);",
        "        const double machine_median = Median(machine_samples);",
        "        std::vector<double> portable_deviation, machine_deviation;",
        "        for (size_t i = 0; i < portable_samples.size(); ++i)",
        "        { portable_deviation.push_back(std::abs(portable_samples[i] - portable_median));",
        "          machine_deviation.push_back(std::abs(machine_samples[i] - machine_median)); }",
        '        std::printf("%s,%u,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\\n",',
        "            pair.family, pair.coefficient, portable_median, machine_median,",
        "            *std::min_element(portable_samples.begin(), portable_samples.end()),",
        "            *std::min_element(machine_samples.begin(), machine_samples.end()),",
        "            Median(portable_deviation), Median(machine_deviation));",
        "    }",
        '    std::fprintf(stderr, "sink=%llu\\n", (unsigned long long)Sink);',
        "    return 0;",
        "}",
    ])
    return "\n".join(lines) + "\n", pairs


def parse_timings(text, pairs):
    lines = text.splitlines()
    header = ("family,coefficient,portable_median_ns,machine_median_ns,"
              "portable_best_ns,machine_best_ns,portable_mad_ns,machine_mad_ns")
    if not lines or lines[0] != header or len(lines) != len(pairs) + 1:
        raise RuntimeError("unexpected equivalent-pair benchmark output")
    result = {}
    expected = set(pairs)
    for line in lines[1:]:
        fields = line.split(",")
        if len(fields) != 8:
            raise RuntimeError("malformed equivalent-pair timing row")
        key = (fields[0], int(fields[1]))
        if key in result or key not in expected:
            raise RuntimeError("duplicate or unknown equivalent-pair timing")
        timing = {
            "portable_median_ns": float(fields[2]),
            "machine_median_ns": float(fields[3]),
            "portable_best_ns": float(fields[4]),
            "machine_best_ns": float(fields[5]),
            "portable_mad_ns": float(fields[6]),
            "machine_mad_ns": float(fields[7]),
        }
        values = list(timing.values())
        if (any(not math.isfinite(value) for value in values) or
                any(timing[name] <= 0 for name in (
                    "portable_median_ns", "machine_median_ns",
                    "portable_best_ns", "machine_best_ns")) or
                timing["portable_mad_ns"] < 0 or
                timing["machine_mad_ns"] < 0):
            raise RuntimeError("invalid equivalent-pair timing value")
        result[key] = timing
    if set(result) != expected:
        raise RuntimeError("equivalent-pair benchmark omitted a circuit")
    return result


def compiler_identity(compiler):
    return run([compiler, "--version"])[0].splitlines()[0]


def validate_options(iterations, rounds, cpu):
    if any(not isinstance(value, int) or isinstance(value, bool)
           for value in (iterations, rounds, cpu)):
        raise ValueError("benchmark counts and CPU must be integers")
    if iterations <= 1 or iterations > UINT64_MAX:
        raise ValueError(
            "--iterations must fit uint64 and be greater than one "
            "(butterflies use half)")
    if rounds < 3 or rounds > UINT32_MAX:
        raise ValueError("--rounds must fit unsigned and be at least three")
    if cpu < 0 or cpu > UINT32_MAX:
        raise ValueError("--cpu must fit unsigned and be non-negative")


def validate_pairs(pairs, workload):
    if not pairs:
        raise ValueError("portable and machine corpora contain no changed pairs")
    total_weight = sum(workload[family].get(coefficient, 0)
                       for family, coefficient in pairs)
    if total_weight <= 0:
        raise ValueError("changed circuit pairs have zero workload weight")
    for family in evaluation.FAMILIES:
        selected = [coefficient for pair_family, coefficient in pairs
                    if pair_family == family]
        if selected and sum(workload[family].get(coefficient, 0)
                            for coefficient in selected) <= 0:
            raise ValueError("changed %s pairs have zero workload weight" %
                             family)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--portable-circuits", type=Path, required=True)
    parser.add_argument("--machine-circuits", type=Path, required=True)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--rounds", type=int, default=15)
    parser.add_argument("--iterations", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    validate_options(arguments.iterations, arguments.rounds, arguments.cpu)

    profile_artifact = cost_model.load_profile_artifact(
        ROOT / "generated" / "FF8XorCostProfiles.json")
    machine_profile = cost_model.find_profile(
        profile_artifact, "amd-zen5-gcc13-avx2-avx512-v1")
    portable_text = arguments.portable_circuits.read_text(encoding="utf-8")
    machine_text = arguments.machine_circuits.read_text(encoding="utf-8")
    portable = {family: evaluation.parse_circuits(portable_text, family)
                for family in evaluation.FAMILIES}
    machine = {family: evaluation.parse_circuits(machine_text, family)
               for family in evaluation.FAMILIES}
    source, pairs = generated_source(portable, machine)
    workload = evaluation.weights(profile_artifact)
    validate_pairs(pairs, workload)
    flags = ("-std=c++11", "-O3", "-DNDEBUG", "-march=x86-64", "-mavx2",
             "-fno-tree-reassoc", "-fno-omit-frame-pointer")
    with tempfile.TemporaryDirectory(prefix="leopard-ff8xor-pairs-") as path:
        source_path = Path(path) / "pairs.cpp"
        artifact = Path(path) / "pairs"
        source_path.write_text(source, encoding="utf-8")
        run([arguments.compiler] + list(flags) + [str(source_path), "-o", str(artifact)])
        stdout, stderr = run([
            str(artifact), str(arguments.iterations), str(arguments.rounds),
            str(arguments.cpu)])
        if stderr and not all(line.startswith("sink=") for line in stderr.splitlines()):
            raise RuntimeError("unexpected pair benchmark diagnostic: %s" % stderr)
        timings = parse_timings(stdout, pairs)

    records = []
    for family, coefficient in pairs:
        width = 8 if family == "multiply" else 16
        portable_metrics = cost_model.source_metrics(
            portable[family][coefficient], width, machine_profile)
        machine_metrics = cost_model.source_metrics(
            machine[family][coefficient], width, machine_profile)
        portable_score = cost_model.score_metrics(
            portable_metrics, machine_profile)
        machine_score = cost_model.score_metrics(machine_metrics, machine_profile)
        timing = timings[(family, coefficient)]
        records.append({
            "family": family,
            "coefficient": coefficient,
            "workload_weight": workload[family].get(coefficient, 0),
            "portable_metrics": portable_metrics,
            "machine_metrics": machine_metrics,
            "predicted_saving": portable_score - machine_score,
            "portable_score": portable_score,
            "machine_score": machine_score,
            "timing": timing,
            "measured_speedup": (
                timing["portable_median_ns"] /
                timing["machine_median_ns"]),
        })

    predicted = [record["predicted_saving"] for record in records]
    raw_gate_savings = [
        record["portable_metrics"]["literal_xor2_count"] -
        record["machine_metrics"]["literal_xor2_count"]
        for record in records]
    speedups = [record["measured_speedup"] for record in records]
    weights = [record["workload_weight"] for record in records]
    total_weight = sum(weights)
    weighted_log_speedup = sum(
        weight * math.log(speedup)
        for weight, speedup in zip(weights, speedups)) / total_weight
    by_family = {}
    for family in evaluation.FAMILIES:
        selected = [record for record in records if record["family"] == family]
        family_weight = sum(record["workload_weight"] for record in selected)
        by_family[family] = {
            "pair_count": len(selected),
            "median_win_count": sum(
                record["measured_speedup"] > 1 for record in selected),
            "median_geomean_speedup": evaluation.geometric_mean(
                [record["measured_speedup"] for record in selected]),
            "workload_weighted_geomean_speedup": math.exp(sum(
                record["workload_weight"] * math.log(record["measured_speedup"])
                for record in selected) / family_weight),
            "raw_gate_saving_vs_measured_speedup_spearman":
                cost_model.spearman([
                    record["portable_metrics"]["literal_xor2_count"] -
                    record["machine_metrics"]["literal_xor2_count"]
                    for record in selected], [
                        record["measured_speedup"] for record in selected]),
            "predicted_saving_vs_measured_speedup_spearman":
                cost_model.spearman([
                    record["predicted_saving"] for record in selected], [
                        record["measured_speedup"] for record in selected]),
        }
    result = {
        "schema": "leopard.ff8xor.equivalent-pair-calibration.v1",
        "portable_circuit_file_sha256":
            evaluation.file_checksum(arguments.portable_circuits),
        "machine_circuit_file_sha256":
            evaluation.file_checksum(arguments.machine_circuits),
        "machine_profile_checksum_sha256":
            machine_profile["checksum_sha256"],
        "workload_schedule_checksum_sha256":
            profile_artifact["workload"]["schedule_checksum_sha256"],
        "compiler": compiler_identity(arguments.compiler),
        "compiler_flags": list(flags),
        "cpu_model": cpu_model(),
        "pinned_cpu": arguments.cpu,
        "rounds": arguments.rounds,
        "base_iterations": arguments.iterations,
        "measurement_order": "same-process ABBA for each equivalent pair",
        "semantic_validation": "portable and machine outputs compared before timing",
        "pair_count": len(records),
        "median_win_count": sum(speedup > 1 for speedup in speedups),
        "median_geomean_speedup": evaluation.geometric_mean(speedups),
        "workload_weighted_geomean_speedup": math.exp(weighted_log_speedup),
        "predicted_saving_vs_measured_speedup_spearman":
            cost_model.spearman(predicted, speedups),
        "raw_gate_saving_vs_measured_speedup_spearman":
            cost_model.spearman(raw_gate_savings, speedups),
        "by_family": by_family,
        "records": records,
    }
    result["checksum_sha256"] = cost_model.checksum(result)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print("Generated %s" % arguments.output)
    print("Pair checksum: %s" % result["checksum_sha256"])
    print("Pairs: %d wins: %d geomean: %.6f weighted: %.6f correlation: %.6f" % (
        result["pair_count"], result["median_win_count"],
        result["median_geomean_speedup"],
        result["workload_weighted_geomean_speedup"],
        result["predicted_saving_vs_measured_speedup_spearman"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
