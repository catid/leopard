#!/usr/bin/env bash
# Bounded-job resource evidence, not a benchmark timer.
set -uo pipefail
"$@"
check_status=$?
printf 'TOWER_CHILD_EXIT=%s\n' "$check_status"
resource_group=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
for name in memory.peak memory.max memory.events memory.swap.current memory.swap.max; do
    printf '%s\n' "$name"
    cat "/sys/fs/cgroup$resource_group/$name"
done
exit "$check_status"
