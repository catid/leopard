#!/usr/bin/env bash
set -euo pipefail
check_service() {
    local condition_unit="$1"
    shift
    local condition_text
    condition_text=$(systemctl "$@" show "$condition_unit" -p Id -p ActiveState -p MainPID -p UnitFileState)
    printf '%s\n' "$condition_text"
    grep -Fxq 'MainPID=0' <<< "$condition_text"
    grep -Fxq 'ActiveState=inactive' <<< "$condition_text"
    grep -Fxq 'UnitFileState=disabled' <<< "$condition_text"
}
check_service slipgate-catstream.service --user
check_service slipgate-obs.service --user
check_service slipgate-headless-x.service
for condition_container in 3fad3dfb247b3d4f3991e0a6801fba0805aae0901bdc89018841d4707a2d888c 002ba0a05bf860b12f2d799b49da1976aed663cb1494b2f8236751e18dc50182; do
    condition_state=$(docker inspect --format '{{.State.Status}}|{{.State.Running}}|{{.HostConfig.RestartPolicy.Name}}' "$condition_container")
    printf '%s %s\n' "$condition_container" "$condition_state"
    test "$condition_state" = 'exited|false|no'
done
