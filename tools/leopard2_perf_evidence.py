#!/usr/bin/env python3
"""Canonical Linux ``perf stat`` evidence parsing for Leopard2 tools.

Both the lab runner's resume path and the benchmark-matrix collector import
this module.  Keeping the derivation here prevents either consumer from
trusting measurements copied into result JSON instead of the retained raw
``perf-stat.txt`` bytes.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


PERF_EVENT_CANONICAL_ALIASES = {
    "branch-instructions": "branches",
    "branches": "branches",
    "cpu-cycles": "cycles",
    "cycles": "cycles",
}


def canonical_perf_event(event: str) -> str:
    """Return an explicit generic-event alias, never a positional guess."""
    return PERF_EVENT_CANONICAL_ALIASES.get(event, event)


def perf_events_match(requested: str, reported: str) -> bool:
    return canonical_perf_event(requested) == canonical_perf_event(reported)


def parse_counter_value(raw_value: str) -> float | None:
    value = raw_value.strip().replace(",", "")
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed) or parsed < 0.0:
        return None
    return parsed


def parse_perf_stat(
    path: Path | str, requested_events: Sequence[str],
) -> tuple[list[dict], str | None]:
    """Derive canonical measurements from perf's delimiter-format output."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        return [], "cannot read perf-stat output: {}".format(error)

    return parse_perf_stat_text(text, requested_events)


def parse_perf_stat_text(
    text: str, requested_events: Sequence[str],
) -> tuple[list[dict], str | None]:
    """Derive canonical measurements from already captured raw text."""
    lines = text.splitlines()

    rows = []
    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [field.strip() for field in line.split(";")]
        if len(fields) >= 3:
            rows.append(fields)

    measurements = []
    unused = list(range(len(rows)))
    for requested in requested_events:
        selected = None
        for row_index in unused:
            if perf_events_match(requested, rows[row_index][2]):
                selected = row_index
                break
        if selected is None:
            measurements.append({"event": requested, "status": "missing"})
            continue

        unused.remove(selected)
        fields = rows[selected]
        parsed = parse_counter_value(fields[0])
        measurement = {
            "event": requested,
            "reported_event": fields[2],
            "raw_value": fields[0],
            "unit": fields[1],
            "status": "counted" if parsed is not None else "not-counted",
        }
        if parsed is not None:
            measurement["value"] = parsed
        if len(fields) > 3 and fields[3]:
            measurement["runtime"] = fields[3]
        if len(fields) > 4 and fields[4]:
            parsed_percentage = parse_counter_value(fields[4].rstrip("%"))
            if parsed_percentage is not None:
                measurement["running_percentage"] = parsed_percentage
        measurements.append(measurement)

    counted = sum(
        measurement["status"] == "counted" for measurement in measurements)
    if counted == len(requested_events):
        return measurements, None
    return measurements, "{} of {} requested events were counted".format(
        counted, len(requested_events))


def derive_perf_stat(
    path: Path | str, requested_events: Sequence[str],
) -> tuple[list[dict], str, str | None]:
    """Return measurements, availability status, and canonical detail."""
    measurements, detail = parse_perf_stat(path, requested_events)
    status = perf_measurement_status(measurements, requested_events)
    return measurements, status, detail


def perf_measurement_status(
    measurements: Sequence[Mapping[str, object]],
    requested_events: Sequence[str],
) -> str:
    counted = sum(
        measurement.get("status") == "counted"
        for measurement in measurements)
    if counted == len(requested_events):
        return "available"
    if counted:
        return "partial"
    return "unavailable"


def read_perf_stat_evidence(
    path: Path | str, requested_events: Sequence[str],
) -> tuple[dict, list[dict], str, str | None]:
    """Hash and parse one immutable byte snapshot of retained perf output."""
    data = Path(path).read_bytes()
    measurements, detail = parse_perf_stat_text(
        data.decode("utf-8", errors="replace"), requested_events)
    identity = {
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
    }
    return (
        identity, measurements,
        perf_measurement_status(measurements, requested_events), detail)


def perf_probe_command(
    perf_executable: str, events: Iterable[str], probe_python: str,
) -> list[str]:
    """Build the exact manifest-bound preflight argv."""
    command = [
        perf_executable, "stat", "--no-big-num", "-x", ";",
    ]
    for event in events:
        command.extend(("-e", event))
    command.extend(("--", probe_python, "-c", "pass"))
    return command


def probe_command_matches_request(
    request: Mapping[str, object], command: object,
) -> bool:
    """Check the complete manifest-bound probe argv and its fixed workload."""
    recorded = request.get("probe_command")
    executable = request.get("executable")
    events = request.get("events")
    if (not isinstance(recorded, list) or command != recorded or
            not all(isinstance(argument, str) and argument
                    for argument in recorded) or
            len(recorded) < 9 or recorded[-4] != "--" or
            recorded[-2:] != ["-c", "pass"] or
            not isinstance(executable, Mapping) or
            not isinstance(executable.get("path"), str) or
            not isinstance(events, list) or
            not all(isinstance(event, str) for event in events)):
        return False
    return recorded == perf_probe_command(
        executable["path"], events, recorded[-3])
