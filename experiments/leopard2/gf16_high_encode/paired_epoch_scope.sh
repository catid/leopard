#!/usr/bin/env bash
# Bind the resource footer to this exact scope and controller command.
set -euo pipefail
epoch_group=$(awk -F: '$1 == "0" {print $3}' /proc/self/cgroup)
printf 'EPOCH_SCOPE=%s\n' "$epoch_group"
printf 'EPOCH_WRAPPER=%s\n' "${BASH_SOURCE[0]}"
printf 'EPOCH_COMMAND_BEGIN\n'
printf '%s\n' "$@"
printf 'EPOCH_COMMAND_END\n'
epoch_directory=$(dirname -- "${BASH_SOURCE[0]}")
exec /bin/bash "$epoch_directory/tower_public_scope.sh" "$@"
