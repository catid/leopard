#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
compiler="${CXX:-c++}"
cpu="${LEO2_AVX512_CPU:-$(awk -F: '/Cpus_allowed_list/ {gsub(/[[:space:]]/,"",$2); print $2}' /proc/self/status | cut -d, -f1 | cut -d- -f1)}"
common=(-std=c++17 -march=x86-64-v2 -Wall -Wextra -Wpedantic
        -Wconversion -Wsign-conversion -Wshadow -Werror)
release="/tmp/leo2-avx512-codec-$UID"
sanitized="/tmp/leo2-avx512-codec-asan-$UID"
trap 'rm -f "$release" "$sanitized"' EXIT

"$compiler" -O3 -g "${common[@]}" avx512_codec.cpp -o "$release"
"$compiler" -O1 -g "${common[@]}" -fsanitize=address,undefined \
    -fno-omit-frame-pointer avx512_codec.cpp -o "$sanitized"
"$compiler" --version | head -1 > compiler_version.txt
lscpu > lscpu.txt
{
    for item in scaling_driver scaling_governor cpuinfo_min_freq cpuinfo_max_freq; do
        path="/sys/devices/system/cpu/cpu${cpu}/cpufreq/${item}"
        printf '%s=' "$item"
        if [[ -r "$path" ]]; then cat "$path"; else echo unavailable; fi
    done
    printf 'perf_event_paranoid='
    if [[ -r /proc/sys/kernel/perf_event_paranoid ]]; then
        cat /proc/sys/kernel/perf_event_paranoid
    else
        echo unavailable
    fi
} > frequency_environment.txt
taskset -c "$cpu" "$release" --validate > validate.stdout.txt 2> validate.stderr.txt
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
taskset -c "$cpu" "$sanitized" --validate \
    > sanitizer.stdout.txt 2> sanitizer.stderr.txt
taskset -c "$cpu" "$release" --benchmark --cpu "$cpu" \
    > results_run1.csv 2> benchmark_run1.stderr.txt
taskset -c "$cpu" "$release" --benchmark --cpu "$cpu" \
    > results_run2.csv 2> benchmark_run2.stderr.txt
python3 analyze.py results_run1.csv results_run2.csv > comparison.csv
python3 assembly_summary.py "$release" > assembly_summary.csv
sha256sum avx512_codec.cpp analyze.py assembly_summary.py run.sh \
    results_run1.csv results_run2.csv comparison.csv assembly_summary.csv \
    validate.stdout.txt validate.stderr.txt sanitizer.stdout.txt \
    sanitizer.stderr.txt benchmark_run1.stderr.txt benchmark_run2.stderr.txt \
    compiler_version.txt lscpu.txt \
    frequency_environment.txt ../../../docs/leopard2_avx512_codec.md \
    > SHA256SUMS
echo "AVX-512 experiment completed on CPU $cpu"
