#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if [[ -z "${GFNI_CPU:-}" ]]; then
    allowed="$(awk -F: '/Cpus_allowed_list/ { gsub(/[[:space:]]/, "", $2); print $2 }' /proc/self/status)"
    first_group="${allowed%%,*}"
    GFNI_CPU="${first_group%%-*}"
fi

compiler="${CXX:-c++}"
flags=(-O3 -g -std=c++17 -march=native -Wall -Wextra -Werror)

"$compiler" "${flags[@]}" gfni_affine.cpp -o gfni_affine
taskset -c "$GFNI_CPU" ./gfni_affine --validate \
    > validate.stdout.txt 2> validate.stderr.txt
taskset -c "$GFNI_CPU" ./gfni_affine --benchmark \
    > results.csv 2> benchmark.stderr.txt
taskset -c "$GFNI_CPU" ./gfni_affine --benchmark \
    > results_run2.csv 2> benchmark_run2.stderr.txt
python3 analyze.py results.csv results_run2.csv > comparison.csv
sha256sum gfni_affine.cpp analyze.py run.sh results.csv results_run2.csv \
    comparison.csv validate.stderr.txt benchmark.stderr.txt \
    benchmark_run2.stderr.txt > SHA256SUMS

echo "GFNI affine experiment completed on CPU $GFNI_CPU"
