#!/usr/bin/env python3
"""Self and mutation tests for the exact-main ABBA evidence runner."""

from __future__ import annotations

import argparse
import base64
from contextlib import ExitStack
import copy
import errno
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock
from dataclasses import asdict
from pathlib import Path
from typing import Mapping


MODULE_PATH = Path(__file__).with_name("run_abba.py")
SPEC = importlib.util.spec_from_file_location("main_compare_run_abba", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load the main-comparison evidence runner")
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def load_plan_runner():
    path = MODULE_PATH.parents[1] / "decoder_dispatch" / \
        "plan_balanced_promotion.py"
    module_name = "main_compare_test_balanced_promotion"
    existing = sys.modules.get(module_name)
    if existing is not None:
        return existing
    module_spec = importlib.util.spec_from_file_location(module_name, path)
    if module_spec is None or module_spec.loader is None:
        raise RuntimeError("cannot load the balanced-promotion planner")
    module = importlib.util.module_from_spec(module_spec)
    decoder_directory = str(path.parent)
    inserted = decoder_directory not in sys.path
    if inserted:
        sys.path.insert(0, decoder_directory)
    try:
        sys.modules[module_name] = module
        module_spec.loader.exec_module(module)
    finally:
        if inserted:
            sys.path.remove(decoder_directory)
    return module


BASELINE_TREE = "b7c8830d96a978f6ec14fe747095f066e351ae72"
BASELINE_TREE_BASE64 = (
    "MTAwNjQ0IC5naXRpZ25vcmUAd4ZzgeN0hGZMvhC3SNIlEZKrNxwxMDA2NDQgLmdpdG1v"
    "ZHVsZXMAoR0P8aCA8DMN4ErrDW1SZUrZ7gkxMDA2NDQgLnRyYXZpcy55bWwAe2QAdEdU"
    "NkTzN2pl35FU2Znb8OMxMDA2NDQgQmVuY2htYXJrcy5tZAAy22E4N5ZlfvNeS/jjtv5b"
    "RRRoMDEwMDY0NCBDTWFrZUxpc3RzLnR4dABe9Z2BR/FnaxWoF0pFeZG0oty+DjEwMDY0"
    "NCBMZW9wYXJkQ29tbW9uLmNwcACWMpFU4u8AMsbBmIJ84vzZ+EZfxzEwMDY0NCBMZW9w"
    "YXJkQ29tbW9uLmgARbrD1dirgtfN8BKjfnZzqBSp0IYxMDA2NDQgTGVvcGFyZEZGMTYu"
    "Y3BwADJBxm+Ilm/icdojn1SlGjwYKEEWMTAwNjQ0IExlb3BhcmRGRjE2LmgAb4umwORH"
    "9+eWxuwkMXY6R0VeWXsxMDA2NDQgTGVvcGFyZEZGOC5jcHAAL5Qfg/tbWzOL9zJBr89n"
    "WYPGSIcxMDA2NDQgTGVvcGFyZEZGOC5oAP4Z7VXYbeW4HWdOTIx865F/JKNCMTAwNjQ0"
    "IExpY2Vuc2UubWQAd9VDahPkRvTJTSvq02FqlNLYnoIxMDA2NDQgUkVBRE1FLm1kAKBF"
    "mPF6SZtcTi19PM8bHB90tEZ2MTAwNjQ0IGFwcHZleW9yLnltbAC7Sa/AtwInE9m6HIdg"
    "ztpNYm7mRzQwMDAwIGRvY3MApWTCpdc70ujVI5sIb1+E9J9u38YxMDA2NDQgbGVvcGFy"
    "ZC5jcHAADAhH8XYqtMLMBvNUyHhUc5uy1UUxMDA2NDQgbGVvcGFyZC5oACqPp5kbKTIE"
    "bw5nnrWLqOFlbLkZNDAwMDAgcHJvagDi+Iv1XmaaeUq1e8tkoeJBkViXgDE2MDAwMCBz"
    "c2UybmVvbgDK1RipOzJvD2RLeXLUiNBOqisEdTQwMDAwIHRlc3RzAEStx5PlLFmpccsI"
    "T+QhDlDeAC7e")
BASELINE_COMMIT_BASE64 = (
    "dHJlZSBiN2M4ODMwZDk2YTk3OGY2ZWMxNGZlNzQ3MDk1ZjA2NmUzNTFhZTcyCnBhcmVudCAy"
    "MmRkYzc4MDQ5OThkMzFjOGYxYTI2MTdlZTcyMGUwNjNiMWZhNmNkCnBhcmVudCAzNjQyN2Rk"
    "MjViZjY2NWY0MTA1MjVhMDNhMWYwYTBlYTkyNzUxNTBiCmF1dGhvciBDaHJpcyBUYXlsb3Ig"
    "PG1yY2F0aWRAZ21haWwuY29tPiAxNzExMTkyODE2IC0wNzAwCmNvbW1pdHRlciBHaXRIdWIg"
    "PG5vcmVwbHlAZ2l0aHViLmNvbT4gMTcxMTE5MjgxNiAtMDcwMApncGdzaWcgLS0tLS1CRUdJ"
    "TiBQR1AgU0lHTkFUVVJFLS0tLS0KIAogd3NGY0JBQUJDQUFRQlFKbC9ycndDUkMxYVE3dXU1"
    "VWhsQUFBWUVVUUFKQkh1ZjBtUWRxQW9mclZTR3ZHT2ZRZwogSllQQVNZTmh5azhoUEs0N3E3"
    "REhPWVI0UUN3eHZySC9VaVpYcGMyZndSK3FlNUxSRkhqUjFaNE43K3AvK0J6cgogMnd6cXBZ"
    "NmR1SUt5NzVvdjFNR0RGa0FocmpPdlBXbDFCR2RlcTNCYnRQL0NzQ0JVSVEyOCtNZG5OVGRq"
    "dGprVwogMkJMWUFYN3VnU2pKWGQ1OEE0eTJ2MHJvUDZKeVRWbHR3b0NRcUtlYUlTVkNDa0Y0"
    "amcvcmJTZHRUck9WeE1oZAogenlPdmtTMzlzcGVOUk9UaUlsTTB4N3VESkR6ZFZja2draEV2"
    "N2dDN3Z6VWdzRFRHYUh6WXkwTU5ISDlmaHhBaAogcEdNQXEzV3p4ZGgxdjJLQTF5cExvZ1VV"
    "ekwwSzNFRFpvRjREL1RZMG96SDdESm1TUytzV1BVRWpUSHNXQzg5eAogdFR2bk8vSzhLNjZj"
    "SkZDQ1N1LzBvSm4yVFdhQnZUdldnNXJ5WFN2RlNibmVjNGRzQlAzczBZVS9odGJRbDN4Tgog"
    "TTkvMmw3M21ra2NJNHIxOVRHQU0zV0JzYTVxKzBrUHk5QnE2SHloeTM5ampOV2pkTnVzZFFG"
    "a3p0dXNSVFJsNwogYUx5SlNRTXB4cWVzeHRacWVEMjRmQmZYUmFRMlhhYkxFUytsTThBd0ly"
    "dytMUFFScWcrQ20rb3Rxb21BeG9QaQogemVmcUU3OGxXVkJRenVESkpLZnNHaERGY1BFaWJJ"
    "UXYyN0hOT3BsaTBtWkFJbGFnckJiOEF4cjRkbytVS0d6cAogMkhYYlprTjZFTWVGYVRKNzY4"
    "MUhLZHlQRFdiWUM3STVPRFBJNktkUW1kaWFJNmZldUJ4QTgybjU2R3laQ2w1OQogclFHZnFM"
    "dGU3VFVvRnVMTE53c1YKID00MTBSCiAtLS0tLUVORCBQR1AgU0lHTkFUVVJFLS0tLS0KIAoK"
    "TWVyZ2UgcHVsbCByZXF1ZXN0ICMyMyBmcm9tIGdibGV0cjQyL21hc3RlcgoKSGFuZGxlIHRo"
    "ZSBjYXNlIHdoZW4gbm8gb3JpZ2luYWwgZGF0YSB3YXMgbG9zdA==")
BASELINE_TREE_OBJECTS_BASE64 = {
    "44adc793e52c59a971cb084fe4210e50de002ede": (
        "MTAwNjQ0IGJlbmNobWFyay5jcHAA8/MtWl7o7Bdvt5BBzp3e8MnZso8xMDA2NDQg"
        "ZXhwZXJpbWVudHMuY3BwANpEfrpL5nWbdjkbETceAS2Zfe2MNDAwMDAgcHJvagDH"
        "GQZEGb/M9gs9CA8VILkuUuImZw=="),
    "a564c2a5d73bd2e8d5239b086f5f84f49f6edfc6": (
        "MTAwNjQ0IEhpZ2hSYXRlRGVjb2Rlci5wZGYAbOUFRKRKJQRMcXPJeirIKJ3xtTkx"
        "MDA2NDQgTG93UmF0ZURlY29kZXIucGRmAJO6Zey0wawNQZu5uKeEJzPIKeJtMTAw"
        "NjQ0IE5vdmVsUG9seW5vbWlhbEJhc2lzRkZUMjAxNi5wZGYARCnf9d5bbKTOwgWl"
        "skYv3s+YW3MxMDA2NDQgcGxhbmstZmFzdDEzLnBkZgCll+moGI1Ke5iDLwrlqsQv"
        "Rrr4azEwMDY0NCB2ZWN0b3JfZndodF80LnR4dABMrTpvJY4Ip+BsCuyWAPerCc5d"
        "Ew=="),
    "b7c8830d96a978f6ec14fe747095f066e351ae72": BASELINE_TREE_BASE64,
    "c719064419bfccf60b3d080f1520b92e52e22667": (
        "MTAwNjQ0IEJlbmNobWFyay52Y3hwcm9qAEFYP/u2S3vdPw/ORCj4lXwVuC3eMTAw"
        "NjQ0IEJlbmNobWFyay52Y3hwcm9qLmZpbHRlcnMAUKBd1CCxdNe8IEDFjKmkk4JN"
        "3LgxMDA2NDQgRXhwZXJpbWVudHMuZmlsdGVycwBQoF3UILF017wgQMWMqaSTgk3c"
        "uDEwMDY0NCBFeHBlcmltZW50cy52Y3hwcm9qABh9gEr0hJWl9kijPSYnpI97x8Po"),
    "e2f88bf55e669a794ab57bcb64a1e24191589780": (
        "MTAwNjQ0IExlb3BhcmQuc2xuANqp9YiW4MOj8tOJnz2jikBjEnvXMTAwNjQ0IExl"
        "b3BhcmQudmN4cHJvagCM+38iiyYCY9Fr/p3w6R/Wnzz6ajEwMDY0NCBMZW9wYXJk"
        "LnZjeHByb2ouZmlsdGVycwCQ3j+mbqr7xNmN4YBs3wtHtOkyrw=="),
}
SSE2NEON_COMMIT = "cad518a93b326f0f644b7972d488d04eaa2b0475"
SSE2NEON_COMMIT_BASE64 = (
    "dHJlZSA2ODEyZmJlYTE1ZTY5OTU0ZmU1OWFlZjJiMWNkNTc3MjYyNjZlMjcyCnBh"
    "cmVudCAxMGYxMzU3MzcxYjI1MzFjM2JkYTViN2VhYTlmNmI4YTZkMTc1MjNkCmF1"
    "dGhvciBKaW0gSHVhbmcgPGpzZXJ2QGJpaWxhYnMuaW8+IDE2NTI2NzU1MjUgKzA4"
    "MDAKY29tbWl0dGVyIEppbSBIdWFuZyA8anNlcnZAYmlpbGFicy5pbz4gMTY1MjY3"
    "NTUyNSArMDgwMAoKUmVmZXIgdG8gdGhlIGFydGljbGVzIGF0IEFybSB3ZWIgc2l0"
    "ZQo=")
SSE2NEON_TREE_OBJECTS_BASE64 = {
    "4b448767151979de0ef817e2026ba5a016634cdd": (
        "MTAwNjQ0IENPREVPV05FUlMAWMrhK58Gr6ty4fodJyBlLFmN7Zg0MDAwMCB3b3Jr"
        "Zmxvd3MAfWoEvyNBMGScZvrvPsnjETqFX4I="),
    "6812fbea15e69954fe59aef2b1cd57726266e272": (
        "NDAwMDAgLmNpAH7KmO0VXKQe+D8/RRccPsfv15Q9MTAwNjQ0IC5jbGFuZy1mb3Jt"
        "YXQATnzbSQ1fdF7CsKCtW4btqOiJ+48xMDA2NDQgLmdpdGF0dHJpYnV0ZXMAkHZ3"
        "Sd+BUchIeD+9V6ZK6GQXUS80MDAwMCAuZ2l0aHViAEtEh2cVGXneDvgX4gJrpaAW"
        "Y0zdMTAwNjQ0IC5naXRpZ25vcmUAZNPdmM/1yeABFZEbrRrIcBniF1ExMDA2NDQg"
        "Q09OVFJJQlVUSU5HLm1kAHz/2CGl6uU4IROqipur9gXCqf7IMTAwNjQ0IExJQ0VO"
        "U0UAnPEGJyrDtWsMTIAhjo/BCmZMpfQxMDA2NDQgTWFrZWZpbGUAmkNEh+9h8uNA"
        "TRGbDeDBwGvrpmwxMDA2NDQgUkVBRE1FLm1kALhGuXRZw1/Pr8G1Jioz2ovZhBGE"
        "MTAwNjQ0IHNzZTJuZW9uLmgAK98k8nat31DZJNyGEZwCbSkMVZU0MDAwMCB0ZXN0"
        "cwDGTMDmzMETNLVGz+gUXjRFw9eJKQ=="),
    "7d6a04bf234130649c66faef3ec9e3113a855f82": (
        "MTAwNjQ0IGdpdGh1Yl9hY3Rpb25zLnltbAA+TRnw73FLSp1KGKo/gGxfgx7qvg=="),
    "7eca98ed155ca41ef83f3f45171c3ec7efd7943d": (
        "MTAwNzU1IGNoZWNrLWZvcm1hdC5zaACm/cHJIWVjjxXsT5QPk0kvo1d/fTEwMDY0"
        "NCBjb21tb24uc2gAwDlgYrqgjZEhaRwVXoBtPzJX4BgxMDA3NTUgY3Jvc3MtY2hl"
        "Y2suc2gAXvuf786rNckEX+0mOwUlNWSOMjMxMDA3NTUgY3Jvc3MtdG9vbC5zaADx"
        "SeprMfpfcQOTKacmEew5kmIHRw=="),
    "c64cc0e6ccc11334b546cfe8145e3445c3d78929": (
        "MTAwNjQ0IFJFQURNRS5tZACvsfn4CLTf3rC1J6V8ZknMO7NNYzEwMDY0NCBiaW5k"
        "aW5nLmNwcABOnq3bo9AY1OuKRoNQmpCXxFYTGjEwMDY0NCBiaW5kaW5nLmgAZKn6"
        "UG3AS9RE1nzoRDPURTqYjrsxMDA2NDQgY29tbW9uLmNwcADQ26UD3ExXR6VaQp4p"
        "sYOcoDmfXDEwMDY0NCBjb21tb24uaACI+gFnWIAR+k2eDBcafOOynHq/xDEwMDY0"
        "NCBpbXBsLmNwcAAbGU+kfkn7tBPzNd2+Lcun0dN0NTEwMDY0NCBpbXBsLmgADneU"
        "9+bwCeGKb1LI5e8NyFMNSXsxMDA2NDQgbWFpbi5jcHAAlT1MNrURJCYPKuV8PcWD"
        "SPQioUk="),
}
CANDIDATE_TREE = "4b825dc642cb6eb9a060e54bf8d69288fbee4904"
CANDIDATE_COMMIT_RAW = (
    f"tree {CANDIDATE_TREE}\nauthor Fixture <fixture@example.com> 1 +0000\n"
    "committer Fixture <fixture@example.com> 1 +0000\n\nfixture\n"
).encode("utf-8")
CANDIDATE_COMMIT = hashlib.sha1(
    f"commit {len(CANDIDATE_COMMIT_RAW)}\0".encode("ascii") +
    CANDIDATE_COMMIT_RAW).hexdigest()


CELL = runner.Cell("fixture", 8, 4, 64, 2, 7)
SPECIFICATION = {
    "runner": "/fixture/run_abba.py",
    "taskset": "/usr/bin/taskset",
    "ldd": "/usr/bin/ldd",
    "baseline_executable": "/fixture/baseline",
    "candidate_executable": "/fixture/candidate",
    "baseline_archive": "/fixture/baseline.a",
    "candidate_archive": "/fixture/candidate-build/libleopard.a",
    "baseline_build_dir": "/fixture/baseline-build",
    "candidate_build_dir": "/fixture/candidate-build",
    "baseline_source_root": "/fixture/main",
    "candidate_source_root": "/fixture/candidate-source",
    "candidate_commit": CANDIDATE_COMMIT,
}
CAMPAIGN = {
    "rounds": runner.ROUNDS,
    "order": list(runner.ORDER),
    "cells": [asdict(CELL)],
    "candidate_mode": "auto",
    "batch": 1,
    "reuse": 8,
    "iterations": 3,
    "warmup": 2,
    "threads": 1,
    "child_environment": copy.deepcopy(runner.CHILD_ENVIRONMENT),
    "benchmark_cpu": 0,
    "reserved_sibling": 1,
    "allowed_cpu_set_at_launch": [0, 1, 2],
    "timeout_seconds": 10.0,
    "statistics": runner.statistics_policy(),
}
DIGESTS = {
    "algorithm": "fnv1a64",
    "original_data": "0123456789abcdef",
    "transmitted_parity": "1111111111111111",
    # This digest covers only missing originals, not the full original_data
    # sequence.  A partial-loss workload must not require the two to match.
    "recovered_originals": "2222222222222222",
}
RESERVATION_PAYLOAD = {
    "benchmark_cpu": 0,
    "nonce": "fixture-nonce",
    "owner": "unit test",
    "reserved_sibling": 1,
    "schema": runner.RESERVATION_SCHEMA,
    "status": "held",
}
RESERVATION = {
    "path": "/fixture/reservation.json",
    "sha256": runner.sha256_bytes(runner.canonical_bytes(RESERVATION_PAYLOAD)),
    "payload": RESERVATION_PAYLOAD,
    "lock": "exclusive_nonblocking",
}


def initialize_git_fixture(root: Path, content: str = "one\n") -> str:
    subprocess.run(["/usr/bin/git", "init", "-q", str(root)], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                    "user.name", "Leopard2 Test"], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "config",
                    "user.email", "test@example.invalid"], check=True)
    (root / "tracked.txt").write_text(content, encoding="utf-8")
    subprocess.run(["/usr/bin/git", "-C", str(root), "add", "."], check=True)
    subprocess.run(["/usr/bin/git", "-C", str(root), "commit", "-qm",
                    "fixture"], check=True)
    return subprocess.check_output(
        ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
        text=True).strip()


def host_cpu(cpu: int) -> dict:
    return {
        "cpu": cpu,
        "online": "1",
        "cpuinfo": {"processor": str(cpu), "model name": "fixture cpu"},
        "topology": {
            "core_id": "0", "physical_package_id": "0", "die_id": "0",
            "cluster_id": None, "thread_siblings_list": "0-1",
            "core_siblings_list": "0-2",
        },
        "frequency_policy": {
            "scaling_driver": "fixture", "scaling_governor": "performance",
            "energy_performance_preference": "performance",
            "scaling_min_freq": "1", "scaling_max_freq": "2",
            "cpuinfo_min_freq": "1", "cpuinfo_max_freq": "2",
        },
        "cache_hierarchy": [
            {
                "index": 0, "level": "1", "type": "Data", "size": "32K",
                "coherency_line_size": "64", "number_of_sets": "64",
                "ways_of_associativity": "8", "physical_line_partition": "1",
                "shared_cpu_list": "0-1", "shared_cpu_map": "00000003",
                "allocation_policy": None, "write_policy": "WriteBack",
            },
            {
                "index": 1, "level": "3", "type": "Unified", "size": "8M",
                "coherency_line_size": "64", "number_of_sets": "16384",
                "ways_of_associativity": "8", "physical_line_partition": "1",
                "shared_cpu_list": "0-2", "shared_cpu_map": "00000007",
                "allocation_policy": None, "write_policy": "WriteBack",
            },
        ],
        "cache_index_inventory": ["index0", "index1"],
        "cache_directory_inventory_text": runner.exact_text_content(
            "index0\nindex1\n", "fixture cache-directory inventory"),
        "numa_nodes": [0],
        "numa_node_inventory": ["node0"],
        "core_class": {"core_type": None, "cpu_capacity": None},
    }


HOST = {
    "system": {
        "system": "Linux", "node": "fixture", "release": "fixture",
        "version": "fixture", "machine": "x86_64", "python": "3",
        "page_size": 4096,
    },
    "allowed_cpu_set_at_launch": [0, 1, 2],
    "online_cpu_set": [0, 1, 2],
    "online_cpu_list_text": runner.exact_text_content(
        "0-2\n", "fixture online CPU list"),
    "online_node_list_text": runner.exact_text_content(
        "0\n", "fixture online node list"),
    "benchmark_cpu": host_cpu(0),
    "reserved_sibling": host_cpu(1),
    "turbo_and_pstate": {
        "/sys/devices/system/cpu/intel_pstate/no_turbo": "0",
        "/sys/devices/system/cpu/amd_pstate/status": None,
        "/sys/devices/system/cpu/cpufreq/boost": None,
    },
}


def cpu_stat(cpu: int, *, user: int, idle: int) -> dict:
    fields = {
        "user": user, "nice": 0, "system": 0, "idle": idle,
        "iowait": 0, "irq": 0, "softirq": 0, "steal": 0,
    }
    return {"cpu": cpu, "fields": fields, "total_jiffies": sum(fields.values())}


PAIR_LEASE_PAYLOAD = runner.pair_lease_payload(0, 1)
PAIR_LEASE = {
    "device": 1,
    "directory_device": 1,
    "directory_inode": 2,
    "inode": 3,
    "lock": "exclusive_nonblocking_pair_wide",
    "path": str(runner.pair_lease_directory() / runner.pair_lease_name(0, 1)),
    "payload": PAIR_LEASE_PAYLOAD,
    "sha256": runner.sha256_bytes(runner.canonical_bytes(PAIR_LEASE_PAYLOAD)),
}
ISOLATION = runner.isolation_record(
    0, 1, PAIR_LEASE, 1_000, 2_000,
    cpu_stat(0, user=100, idle=100), cpu_stat(0, user=110, idle=110),
    cpu_stat(1, user=100, idle=100), cpu_stat(1, user=100, idle=120),
)
SUPERVISION = runner.supervision_record(
    "ab" * 32, 900, 2_100, CAMPAIGN, RESERVATION, ISOLATION)


def summary(samples: list[float], setup: bool = False) -> dict:
    middle = sorted(samples)[len(samples) // 2]
    deviations = sorted(abs(value - middle) for value in samples)
    suffix = "" if setup else "_per_batch_call"
    result = {
        f"median_us{suffix}": middle,
        f"mad_us{suffix}": deviations[len(deviations) // 2],
        f"minimum_us{suffix}": min(samples),
        f"maximum_us{suffix}": max(samples),
        "samples_us" if setup else "samples_us_per_batch_call": samples,
    }
    return result


def common_parameters() -> dict:
    return {
        "K": CELL.k,
        "R": CELL.r,
        "shard_bytes": CELL.shard_bytes,
        "loss_count": CELL.losses,
        "missing_original_indices": runner.expected_losses(CELL),
        "batch": 1,
        "reuse": CAMPAIGN["reuse"],
        "iterations": CAMPAIGN["iterations"],
        "warmup": CAMPAIGN["warmup"],
        "thread_count": 1,
        "seed": CELL.seed,
    }


def baseline_result(scale: float = 1.0) -> dict:
    encode = summary([10.0 * scale, 11.0 * scale, 12.0 * scale])
    decode = summary([20.0 * scale, 21.0 * scale, 22.0 * scale])
    encode.update({"input_GB_per_s": 1.0, "parity_output_GB_per_s": 0.5})
    decode.update({"offered_received_GB_per_s": 1.0,
                   "repaired_output_GB_per_s": 0.5})
    return {
        "schema": "leopard-main-benchmark-v1",
        "build": {"main_source_commit": runner.MAIN_COMMIT, "cplusplus": 201103},
        "parameters": common_parameters(),
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8",
            "parent_count": 16, "padded_side": 4, "thread_count": 1,
        },
        "correctness": {"round_trip": True},
        "workload_digests": copy.deepcopy(DIGESTS),
        "metrics": {
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": encode,
            "decode_including_setup": decode,
            "rate_semantics": "fixture",
        },
    }


def candidate_result(
    scale: float = 0.8, raw_schema: str = runner.RAW_SCHEMA,
    campaign: Mapping[str, object] = CAMPAIGN,
) -> dict:
    parameters = common_parameters()
    parameters.update({
        "requested_profile": "legacy_high_v1",
        "requested_field": "auto",
        "requested_backend": "auto",
        "skip_legacy": True,
        "retain_samples": True,
    })
    flags = runner.candidate_mode_flags(
        runner.candidate_mode_for_campaign(campaign))
    parameters.update({
        "force_generic_decode": flags["force_generic_decode"],
        "force_specialized_decode": flags["force_specialized_decode"],
    })
    if raw_schema in runner.WORKSPACE_SELECTOR_SCHEMAS:
        parameters.update({name: flags[name] for name in (
            "force_tiled_decode", "force_materialized_decode")})
    codec = summary([3.0, 3.1, 3.2], setup=True)
    plan = summary([4.0, 4.1, 4.2], setup=True)
    encode = summary([10.0 * scale, 11.0 * scale, 12.0 * scale])
    decode = summary([12.0 * scale, 13.0 * scale, 14.0 * scale])
    encode.update({"input_GB_per_s": 1.0, "parity_output_GB_per_s": 0.5})
    decode.update({"offered_received_GB_per_s": 1.0,
                   "repaired_output_GB_per_s": 0.5})
    amortized = decode["median_us_per_batch_call"] + \
        plan["median_us"] / CAMPAIGN["reuse"]
    return {
        "schema": "leopard2-benchmark-v2",
        "build": {"compiler": "fixture"},
        "parameters": parameters,
        "resolved": {
            "profile": "legacy_high_v1", "field": "gf8", "backend": "avx2",
            "thread_count": 1, "parent_count": 16, "padded_side": 4,
        },
        "correctness": {"leopard2_round_trip": True, "legacy_comparison": None},
        "workload_digests": copy.deepcopy(DIGESTS),
        "metrics": {
            "codec_setup": codec,
            "encode_execution": encode,
            "decode_plan_setup": plan,
            "decode_execution": decode,
            "decode_amortized_at_reuse": {
                "reuse_count": CAMPAIGN["reuse"],
                "derived_median_us_per_batch_call": amortized,
                "offered_received_GB_per_s": 1.0,
                "repaired_output_GB_per_s": 0.5,
            },
            "rate_semantics": "fixture",
        },
        "legacy": {
            "available": False,
            "unavailable_reason": "disabled by --skip-legacy",
            "codec_setup": None,
            "decode_timing_includes_setup": True,
            "encode_execution": None,
            "decode_including_setup": None,
        },
    }


def archive_recipe_fixture_text(cmake: dict | Mapping[str, str]) -> str:
    return (
        f"/usr/bin/ar qc {cmake['archive']} "
        f"CMakeFiles/{cmake['target_directory']}/LeopardCommon.cpp.o "
        "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o\n"
        f"/usr/bin/ranlib {cmake['archive']}\n"
    )


def complete_artifact(path: str, kind: str, character: str,
                      *, content: str | None = None) -> dict:
    encoded = None if content is None else content.encode("utf-8")
    return {
        "path": path,
        "kind": kind,
        "size": 1 if encoded is None else len(encoded),
        "mode": 0o755 if kind in {"executable", "compiler", "archiver",
                                   "ranlib", "build_tool"} else 0o644,
        "mtime_ns": 1,
        "sha256": (character * 64 if encoded is None else
                   runner.sha256_bytes(encoded)),
    }


def complete_ninja_graph_fixture(build_dir: str) -> dict:
    contents = {
        "build-Release.ninja":
            "include CMakeFiles/impl-Release.ninja\n",
        "CMakeFiles/impl-Release.ninja":
            "include CMakeFiles/common.ninja\n"
            "build Release/libleopard.a: phony\n"
            "build Release/bench_leopard2: phony Release/libleopard.a\n",
        "CMakeFiles/common.ninja":
            "rule fixture_noop\n"
            "  command = :\n",
    }
    files = []
    for index, relative in enumerate(sorted(contents)):
        content = runner.exact_text_content(
            contents[relative], f"fixture Ninja graph {relative}")
        files.append({
            "relative_path": relative,
            "artifact": complete_artifact(
                f"{build_dir}/{relative}", "ninja_graph_input",
                format(index + 1, "x"), content=contents[relative]),
            "content": content,
        })
    return {
        "schema": runner.NINJA_GRAPH_CLOSURE_SCHEMA,
        "entrypoint": "build-Release.ninja",
        "files": files,
    }


def replace_ninja_graph_fixture_text(
    build: dict, relative_path: str, text: str,
) -> None:
    graph = build["multi_config_ninja_graph"]
    record = next(
        item for item in graph["files"]
        if item["relative_path"] == relative_path)
    content = runner.exact_text_content(
        text, f"mutated Ninja graph {relative_path}")
    record["content"] = content
    record["artifact"]["size"] = content["size"]
    record["artifact"]["sha256"] = content["sha256"]
    if relative_path == "CMakeFiles/impl-Release.ninja":
        for name in ("archive_link_recipe", "executable_link_recipe"):
            build[name]["size"] = content["size"]
            build[name]["sha256"] = content["sha256"]


def compile_commands_fixture(root: Path, role: str) -> tuple[Path, dict, Path]:
    """Materialize one canonical compile database without invoking a compiler."""
    baseline_root = (root / "baseline-source").resolve()
    candidate_root = (root / "candidate-source").resolve()
    build = (root / f"{role}-build").resolve()
    compiler = (root / "toolchain" / "c++").resolve()
    for directory in (baseline_root, candidate_root, build, compiler.parent):
        directory.mkdir(parents=True, exist_ok=True)
    compiler.write_bytes(b"fixture compiler\n")
    compiler.chmod(0o755)
    specification = {
        "baseline_source_root": str(baseline_root),
        "candidate_source_root": str(candidate_root),
        "candidate_commit": CANDIDATE_COMMIT,
        f"{role}_build_dir": str(build),
    }
    build_configuration = None
    if role == "candidate":
        configuration_entries = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_CONFIGURATION_TYPES": "",
            "CMAKE_CXX_COMPILER": str(compiler),
            "CMAKE_CXX_FLAGS": " -Wall -Wextra -fopenmp",
            "CMAKE_CXX_FLAGS_DEBUG": "-g -O0",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
            "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
            "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
            "ENABLE_OPENMP": "ON",
            "LEOPARD_ENABLE_GF8": "ON",
            "LEOPARD_ENABLE_GF16": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
            "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
            "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
        }
        material = runner.build_configuration_material(configuration_entries)
        digest = runner.sha256_bytes(material)
        sidecar = build / runner.BUILD_CONFIGURATION_RELATIVE_PATH
        sidecar.parent.mkdir(parents=True, exist_ok=True)
        sidecar.write_bytes(
            (f"schema={runner.BUILD_CONFIGURATION_FILE_SCHEMA}\n"
             f"sha256={digest}\n").encode("ascii") + material)
        helper = candidate_root / \
            runner.BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH
        helper.parent.mkdir(parents=True, exist_ok=True)
        helper.write_text("# fixture helper\n", encoding="utf-8")
        cache_values = {
            **configuration_entries,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                runner.BUILD_CONFIGURATION_FILE_SCHEMA,
            "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256": digest,
        }
        (build / "CMakeCache.txt").write_text("".join(
            f"{name}:{next(iter(runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES[name]))}"
            f"={value}\n"
            for name, value in cache_values.items()), encoding="utf-8")
        build_configuration = runner.capture_candidate_build_configuration(
            specification, cache_values)
    else:
        cache_values = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_GENERATOR": "Unix Makefiles",
            "CMAKE_CONFIGURATION_TYPES": "",
            "CMAKE_CXX_COMPILER": str(compiler),
        }
        (build / "CMakeCache.txt").write_text("".join(
            f"{name}:{next(iter(runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES[name]))}"
            f"={value}\n"
            for name, value in cache_values.items()), encoding="utf-8")
    if role == "baseline":
        sources = [baseline_root / name for name in runner.BASELINE_LIBRARY_SOURCES]
        adapter = candidate_root / \
            "experiments/leopard2/main_compare/legacy_main_benchmark.cpp"
        sources.append(adapter)
    else:
        sources = [
            candidate_root / name for name in runner.CANDIDATE_CONFIGURED_SOURCES
        ]
    entries = []
    for source in sources:
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text("// fixture source\n", encoding="utf-8")
        output = runner.expected_compile_output(role, source, specification)
        output_path = build / output
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(b"fixture object\n")
        arguments = runner.expected_compile_argv(
            role, source, specification, str(compiler),
            build_configuration=build_configuration)
        entry = {
            "directory": str(build), "file": str(source), "output": output,
        }
        # Exercise both representations while retaining one exact token stream.
        if len(entries) % 2:
            entry["arguments"] = arguments
        else:
            entry["command"] = runner.shlex.join(arguments)
        entries.append(entry)
    path = build / "compile_commands.json"
    path.write_text(json.dumps(entries), encoding="utf-8")
    return path, specification, compiler


def complete_build_fixture(
    role: str, raw_schema: str = runner.RAW_SCHEMA, *,
    multi_config: bool = False,
) -> dict:
    if multi_config and raw_schema not in runner.BUILD_CLOSURE_V7_SCHEMAS:
        raise ValueError("multi-config fixture requires the current schema")
    baseline = role == "baseline"
    build_dir = SPECIFICATION[f"{role}_build_dir"]
    source_root = SPECIFICATION[f"{role}_source_root"]
    selected_configuration = "Release" if multi_config else None
    configuration_types = "Debug;Release" if multi_config else ""
    generator = "Ninja Multi-Config" if multi_config else "Unix Makefiles"
    if baseline:
        archive_name = "libleopard_main_exact.a"
        executable_name = "leopard_main_benchmark"
        target_directory = "leopard_main_exact.dir"
        benchmark_source = (SPECIFICATION["candidate_source_root"] +
            "/experiments/leopard2/main_compare/legacy_main_benchmark.cpp")
        cache = {
            "CMAKE_BUILD_TYPE": "" if multi_config else "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "LEO_MAIN_HAS_MARCH_NATIVE": "1",
        }
        isa_policy = "whole-build -march=native"
        library_names = runner.BASELINE_LIBRARY_SOURCES
        entry_count = 5
    else:
        archive_name = runner.CANONICAL_CMAKE_IDENTITY["archive"]
        executable_name = "bench_leopard2"
        target_directory = runner.CANONICAL_CMAKE_IDENTITY["target_directory"]
        benchmark_source = source_root + "/bench/leopard2/benchmark.cpp"
        cache = {
            "CMAKE_BUILD_TYPE": "" if multi_config else "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        if raw_schema in runner.BUILD_CLOSURE_V7_SCHEMAS:
            (configuration_record_schema, configuration_file_schema,
             configuration_variables) = \
                runner.build_configuration_contract_for_raw_schema(raw_schema)
            cache.update({
                "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
                "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
                **({
                    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
                } if raw_schema == runner.RAW_SCHEMA else {}),
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
            })
            configuration_entries = {
                "CMAKE_BUILD_TYPE": "" if multi_config else "Release",
                "CMAKE_GENERATOR": generator,
                "CMAKE_CONFIGURATION_TYPES": configuration_types,
                "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
                "CMAKE_CXX_FLAGS": " -Wall -Wextra -fopenmp",
                "CMAKE_CXX_FLAGS_DEBUG": "-g -O0",
                "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG -O3",
                "CMAKE_CXX_FLAGS_RELWITHDEBINFO": "-O2 -g -DNDEBUG",
                "CMAKE_CXX_FLAGS_MINSIZEREL": "-Os -DNDEBUG",
                "ENABLE_OPENMP": "ON",
                "LEOPARD_ENABLE_GF8": "ON",
                "LEOPARD_ENABLE_GF16": "ON",
                "LEO2_BACKEND_VARIANT": "auto",
                "LEO2_BENCHMARK_GIT_EXECUTABLE": "/usr/bin/git",
                "LEO2_BUILD_BENCHMARKS": "ON",
                "LEO2_BUILD_TESTS": "OFF",
                "LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN": "OFF",
                "LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE": "OFF",
                **({
                    "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT": "OFF",
                } if raw_schema == runner.RAW_SCHEMA else {}),
                "LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE": "0",
            }
            configuration_material = runner.build_configuration_material(
                configuration_entries, configuration_variables)
            configuration_digest = runner.sha256_bytes(
                configuration_material)
            configuration_text = (
                f"schema={configuration_file_schema}\n"
                f"sha256={configuration_digest}\n").encode("ascii") + \
                configuration_material
            configuration_text_value = configuration_text.decode("utf-8")
            effective_configuration = {
                "schema": configuration_record_schema,
                "artifact": complete_artifact(
                    build_dir + "/" +
                    runner.BUILD_CONFIGURATION_RELATIVE_PATH,
                    "generated_build_configuration", "4",
                    content=configuration_text_value),
                "content": runner.exact_text_content(
                    configuration_text_value,
                    "fixture effective build configuration"),
                "configuration_schema":
                    configuration_file_schema,
                "configuration_sha256": configuration_digest,
                "entries": configuration_entries,
                "embedded_build_type": "Release",
                "helper_source": complete_artifact(
                    source_root + "/" +
                    runner.BENCHMARK_ATTESTATION_HELPER_RELATIVE_PATH,
                    "source_file", "5"),
            }
            cache.update({
                "CMAKE_CONFIGURATION_TYPES": configuration_types,
                "CMAKE_GENERATOR": generator,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SCHEMA":
                    configuration_file_schema,
                "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256":
                    configuration_digest,
            })
        else:
            effective_configuration = None
        isa_policy = (
            "portable core with ISA flags only on SSSE3, AVX2, and "
            "AVX-512VL translation units")
        library_names = runner.CANDIDATE_LIBRARY_SOURCES
        entry_count = runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT
    if multi_config and not baseline:
        object_library_sources = (
            "Leopard2BackendSSSE3.cpp",
            "Leopard2BackendAVX2.cpp",
            "Leopard2BackendAVX2Xor.cpp",
            "Leopard2BackendAVX512.cpp",
            "Leopard2BackendGFNI.cpp",
        )
        library_names = (
            *object_library_sources,
            *(name for name in library_names
              if name not in object_library_sources),
        )
    if raw_schema in runner.BUILD_CLOSURE_V7_SCHEMAS:
        cache.update({
            "CMAKE_CONFIGURATION_TYPES": configuration_types,
            "CMAKE_GENERATOR": generator,
        })
        if multi_config:
            cache["CMAKE_MAKE_PROGRAM"] = "/usr/bin/ninja"
    entry_count *= 2 if multi_config else 1
    benchmark_object_relative = runner.expected_compile_output(
        role, Path(benchmark_source), SPECIFICATION, selected_configuration)
    benchmark_object = build_dir + "/" + benchmark_object_relative
    library_pairs = []
    archive_objects = []
    for index, name in enumerate(library_names):
        source_path = source_root + "/" + name
        object_relative = runner.expected_compile_output(
            role, Path(source_path), SPECIFICATION, selected_configuration)
        archive_objects.append(object_relative)
        library_pairs.append({
            "source": complete_artifact(
                source_path, "source_file",
                format((index % 9) + 1, "x")),
            "object": complete_artifact(
                build_dir + "/" + object_relative, "object_file",
                format(((index + 3) % 9) + 1, "x")),
        })
    pairs = sorted([
        *library_pairs,
        {
            "source": complete_artifact(
                benchmark_source, "source_file", "5" if baseline else "6"),
            "object": complete_artifact(
                benchmark_object, "object_file", "7" if baseline else "8"),
        },
    ], key=lambda record: record["source"]["path"])
    output_prefix = (
        f"{selected_configuration}/" if selected_configuration else "")
    archive_recipe_name = output_prefix + archive_name
    executable_recipe_name = output_prefix + executable_name
    archive_recipe_text = (
        f"/usr/bin/ar qc {archive_recipe_name} "
        + " ".join(archive_objects) + "\n"
        + f"/usr/bin/ranlib {archive_recipe_name}\n")
    archive_recipe_content = runner.exact_text_content(
        archive_recipe_text, f"fixture {role} archive recipe")
    executable_recipe_text = (
        f"/usr/bin/compiler -O3 -fopenmp {archive_recipe_name} "
        f"{Path(benchmark_object).relative_to(Path(build_dir)).as_posix()} "
        f"-o {executable_recipe_name} "
        "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so "
        "/usr/lib/x86_64-linux-gnu/libpthread.a\n")
    executable_recipe_content = runner.exact_text_content(
        executable_recipe_text, f"fixture {role} executable recipe")
    archive_path = build_dir + "/" + archive_recipe_name
    executable_path = build_dir + "/" + executable_recipe_name
    if selected_configuration:
        executable_recipe_path = (
            build_dir + f"/CMakeFiles/impl-{selected_configuration}.ninja")
        archive_recipe_path = executable_recipe_path
    else:
        executable_recipe_path = (
            build_dir + f"/CMakeFiles/{executable_name}.dir/link.txt")
        archive_recipe_path = (
            build_dir + f"/CMakeFiles/{target_directory}/link.txt")
    ninja_graph = (
        complete_ninja_graph_fixture(build_dir) if multi_config else None)
    if ninja_graph is not None:
        impl_artifact = copy.deepcopy(next(
            record["artifact"] for record in ninja_graph["files"]
            if record["relative_path"] ==
                "CMakeFiles/impl-Release.ninja"))
        impl_artifact["kind"] = "build_metadata"
        executable_recipe_artifact = copy.deepcopy(impl_artifact)
        archive_recipe_artifact = copy.deepcopy(impl_artifact)
    else:
        executable_recipe_artifact = complete_artifact(
            executable_recipe_path, "build_metadata", "b",
            content=executable_recipe_text)
        archive_recipe_artifact = complete_artifact(
            archive_recipe_path, "build_metadata", "c",
            content=archive_recipe_text)
    compiler_text = "fixture compiler 1.0\n"
    external_link_inputs = runner.capture_external_link_inputs(
        runner.shlex.split(executable_recipe_text),
        f"fixture {role} executable recipe")
    validated_executable = complete_artifact(
        executable_path, "executable", "0" if baseline else "1")
    validated_executable["mtime_ns"] = 1 + max(
        record["artifact"]["mtime_ns"] for record in external_link_inputs)
    compile_identity = {
        "schema": runner.compile_commands_schema_for_raw_schema(raw_schema),
        "implementation": role,
        "profile": runner.compile_profile_for_implementation(role, raw_schema),
        "entry_count": entry_count,
        "required_sources": sorted(pair["source"]["path"] for pair in pairs),
        "validated_optimization": "-O3",
        "validated_openmp": True,
        "required_source_object_pairs": pairs,
        "required_entries": sorted([{
            "directory": build_dir,
            "file": pair["source"]["path"],
            "output": runner.expected_compile_output(
                role, Path(pair["source"]["path"]), SPECIFICATION,
                selected_configuration),
            "arguments": runner.expected_compile_argv(
                role, Path(pair["source"]["path"]), SPECIFICATION,
                "/usr/bin/compiler", raw_schema,
                CANDIDATE_TREE if role == "candidate" else None,
                effective_configuration if not baseline else None,
                selected_configuration),
        } for pair in pairs], key=lambda entry: entry["file"]),
        "isa_policy": isa_policy,
    }
    if raw_schema in (
            runner.RAW_SCHEMA_V6, runner.RAW_SCHEMA_V7,
            runner.RAW_SCHEMA_V8, runner.RAW_SCHEMA):
        if baseline:
            compile_identity["generated_attestation_header"] = None
        else:
            header_bytes = runner.benchmark_attestation_header_bytes(
                CANDIDATE_COMMIT, CANDIDATE_TREE, False)
            header_text = header_bytes.decode("ascii")
            header_path = (
                build_dir +
                "/generated/leopard2-benchmark-attestation/"
                "leopard2_benchmark_source_attestation.h")
            compile_identity["generated_attestation_header"] = {
                "artifact": complete_artifact(
                    header_path, "generated_compile_input", "4",
                    content=header_text),
                "content": runner.exact_text_content(
                    header_text, "fixture generated attestation"),
                "source_commit": CANDIDATE_COMMIT,
                "source_tree": CANDIDATE_TREE,
                "source_tracked_dirty": False,
            }
    if raw_schema in runner.BUILD_CLOSURE_V7_SCHEMAS:
        compile_identity["effective_build_configuration"] = (
            None if baseline else effective_configuration)
    result = {
        "build_dir": build_dir,
        "cmake_cache": complete_artifact(
            build_dir + "/CMakeCache.txt", "build_metadata", "9"),
        "compile_commands": complete_artifact(
            build_dir + "/compile_commands.json", "build_metadata", "a"),
        "executable_link_recipe": executable_recipe_artifact,
        "archive_link_recipe": archive_recipe_artifact,
        "compiler": complete_artifact("/usr/bin/compiler", "compiler", "d"),
        "compiler_invocation": {
            "invocation": "/usr/bin/compiler",
            "resolved_path": "/usr/bin/compiler",
        },
        "compiler_version_stdout": {
            "sha256": runner.sha256_bytes(compiler_text.encode("utf-8")),
            "text": compiler_text,
        },
        "archiver": complete_artifact("/usr/bin/ar", "archiver", "e"),
        "ranlib": complete_artifact("/usr/bin/ranlib", "ranlib", "f"),
        "validated_archive_members": [Path(name).name + ".o"
                                       for name in library_names],
        "validated_executable": validated_executable,
        "validated_archive": complete_artifact(
            archive_path, "archive", "2" if baseline else "3"),
        "validated_cache": cache,
        "validated_compile_commands": compile_identity,
        "archive_link_recipe_content": archive_recipe_content,
        "executable_link_recipe_content": executable_recipe_content,
        "archive_link_tool_invocations": {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                           "resolved_path": "/usr/bin/ranlib"},
        },
        "validated_external_link_inputs": external_link_inputs,
    }
    if raw_schema in runner.BUILD_CLOSURE_V7_SCHEMAS:
        ninja_text = "1.11.1\n"
        result["multi_config_build_tool"] = (
            complete_artifact(
                runner.CANONICAL_NINJA_PATH, "build_tool", "7")
            if multi_config else None)
        result["multi_config_build_tool_version_stdout"] = (
            {
                "sha256": runner.sha256_bytes(ninja_text.encode("utf-8")),
                "text": ninja_text,
            } if multi_config else None)
        result["multi_config_ninja_graph"] = ninja_graph
    return result


def complete_runtime_fixture(
    executable: str, character: str, raw_schema: str,
) -> dict:
    raw = (
        "linux-vdso.so.1 (0x00000000)\n"
        "libc.so.6 => /lib/libc.so.6 (0x00000000)\n"
        "libm.so.6 => /lib/libm.so.6 (0x00000000)\n"
        "/lib64/ld-linux-x86-64.so.2 (0x00000000)\n")
    result = {
        "executable": executable,
        "dependencies": [
            {
                "soname": "ld-linux-x86-64.so.2",
                "loader_path": "/lib64/ld-linux-x86-64.so.2",
                "file": complete_artifact(
                    "/lib64/ld-linux-x86-64.so.2",
                    "dynamic_loader", character),
            },
            {
                "soname": "libc.so.6",
                "loader_path": "/lib/libc.so.6",
                "file": complete_artifact(
                    "/lib/libc.so.6",
                    "shared_library", character),
            },
            {
                "soname": "libm.so.6",
                "loader_path": "/lib/libm.so.6",
                "file": complete_artifact(
                    "/lib/libm.so.6", "shared_library",
                    hex((int(character, 16) + 1) & 15)[2:]),
            },
            {"soname": "linux-vdso.so.1", "virtual": True},
        ],
    }
    if raw_schema == runner.RAW_SCHEMA_V5:
        result["raw_ldd_output"] = runner.exact_text_content(
            raw, "fixture raw ldd output")
    elif raw_schema in (
            runner.RAW_SCHEMA_V6, runner.RAW_SCHEMA_V7,
            runner.RAW_SCHEMA_V8, runner.RAW_SCHEMA):
        result["canonical_ldd_output"] = runner.canonical_ldd_output(
            raw, "fixture raw ldd output")
    else:
        raise ValueError("complete runtime fixture requires schema v5, v6, or v7")
    return result


def retained_ldd_text(closure: dict) -> str:
    key = ("canonical_ldd_output" if "canonical_ldd_output" in closure else
           "raw_ldd_output")
    return closure[key]["text"]


def replace_retained_ldd_text(closure: dict, text: str, label: str) -> None:
    if "canonical_ldd_output" in closure:
        identity = runner.exact_text_content(text, label)
        closure["canonical_ldd_output"] = {
            "schema": runner.CANONICAL_LDD_SCHEMA,
            "normalization": runner.CANONICAL_LDD_NORMALIZATION,
            **identity,
        }
    else:
        closure["raw_ldd_output"] = runner.exact_text_content(text, label)


def git_object_fixture(kind: str, raw: bytes) -> dict:
    object_id = hashlib.sha1(
        f"{kind} {len(raw)}\0".encode("ascii") + raw).hexdigest()
    return {
        "encoding": "base64", "size": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "object_id": object_id,
        "base64": base64.b64encode(raw).decode("ascii"),
    }


def rich_git_source_fixture(
    path: str,
    head: str,
    commit_raw: bytes,
    tree_objects_base64: Mapping[str, str],
    *,
    detached: bool,
    nested: Mapping[str, dict] | None = None,
) -> dict:
    nested = {} if nested is None else dict(nested)
    commit_object = git_object_fixture("commit", commit_raw)
    if commit_object["object_id"] != head:
        raise AssertionError("fixture commit bytes do not match HEAD")
    tree = commit_raw.split(b"\n", 1)[0].decode("ascii").split(" ", 1)[1]
    tree_objects = []
    decoded_trees = {}
    for object_id, encoded in sorted(tree_objects_base64.items()):
        raw = base64.b64decode(encoded, validate=True)
        identity = git_object_fixture("tree", raw)
        if identity["object_id"] != object_id:
            raise AssertionError("fixture tree bytes do not match object ID")
        tree_objects.append(identity)
        decoded_trees[object_id] = raw
    leaves = runner.git_capture._flatten_tree_objects(tree, decoded_trees)
    records = []
    submodules = []
    for entry in leaves:
        if entry["git_type"] == "blob":
            records.append({**entry, "kind": "regular"})
            continue
        identity = nested.get(entry["path"])
        if identity is None:
            raise AssertionError("fixture submodule identity is absent")
        identity_sha256 = runner.git_capture._digest(identity)
        records.append({
            **entry, "kind": "submodule",
            "identity_sha256": identity_sha256,
        })
        submodules.append({
            "path": entry["path"],
            "object_id": entry["object_id"],
            "identity_sha256": identity_sha256,
            "identity": identity,
        })
    tree_listing = b"".join(
        (
            f"{record['git_mode']} {record['git_type']} "
            f"{record['object_id']}\t{record['path']}\0"
        ).encode("utf-8")
        for record in records)
    index_stage = b"".join(
        (
            f"{record['git_mode']} {record['object_id']} 0\t"
            f"{record['path']}\0"
        ).encode("utf-8")
        for record in records)
    default_flags = b"".join(
        f"H {record['path']}\0".encode("utf-8")
        for record in records)
    empty_bytes = {
        "size": 0, "sha256": hashlib.sha256(b"").hexdigest()}
    git_sha = "a" * 64
    gitdir = path + "/.git"
    return {
        "schema": runner.git_capture.SCHEMA,
        "path": path,
        "head": head,
        "tree": tree,
        "detached": detached,
        "head_ref": None if detached else "refs/heads/fixture",
        "superproject_worktree": None,
        "tracked_tree_listing_sha256":
            hashlib.sha256(tree_listing).hexdigest(),
        "tracked_status": "clean",
        "commit_object": commit_object,
        "tree_objects": tree_objects,
        "git_executable": {
            "source": {
                "path": "/usr/bin/git", "size": 1, "mode": 0o755,
                "sha256": git_sha,
            },
            "sealed": {
                "protocol": runner.git_capture.GIT_EXECUTABLE_PROTOCOL,
                "size": 1, "mode": 0o755, "sha256": git_sha,
                "seals": runner.git_capture.REQUIRED_SEALS,
                "source_sha256": git_sha,
            },
        },
        "git_metadata": {
            "layout": "ordinary", "gitdir": gitdir, "commondir": gitdir,
            "guarded_components": [gitdir],
            "guard_policy": runner.git_capture.METADATA_GUARD_POLICY,
            "guarded_file_count": 0,
            "guarded_files_sha256": hashlib.sha256(
                runner.canonical_bytes([])).hexdigest(),
        },
        "worktree_guard_policy":
            runner.git_capture.WORKTREE_GUARD_POLICY,
        "config": copy.deepcopy(empty_bytes),
        "index": {
            "entry_count": len(records),
            "stage": runner.git_capture._byte_identity(index_stage),
            "flags_v": runner.git_capture._byte_identity(default_flags),
            "flags_f": runner.git_capture._byte_identity(default_flags),
        },
        "tracked_files": records,
        "tracked_files_sha256": runner.git_capture._digest(records),
        "submodules": submodules,
    }


def complete_source_fixture(
        role: str, raw_schema: str = runner.RAW_SCHEMA) -> dict:
    baseline = role == "baseline"
    raw = (base64.b64decode(BASELINE_COMMIT_BASE64)
           if baseline else CANDIDATE_COMMIT_RAW)
    head = runner.MAIN_COMMIT if baseline else CANDIDATE_COMMIT
    tree = BASELINE_TREE if baseline else CANDIDATE_TREE
    commit_object = runner.git_commit_object_identity(raw, head)
    runner.validate_git_commit_object_identity(
        commit_object, head, tree, f"fixture {role}")
    legacy = {
        "path": SPECIFICATION[f"{role}_source_root"],
        "head": head, "tree": tree, "detached": baseline,
        "tracked_tree_listing_sha256": "7" * 64 if baseline else "9" * 64,
        "tracked_status": "clean", "commit_object": commit_object,
    }
    if raw_schema not in runner.SEALED_EXECUTABLE_SCHEMAS:
        return legacy
    source_root = SPECIFICATION[f"{role}_source_root"]
    if not baseline:
        return rich_git_source_fixture(
            source_root, head, raw, {CANDIDATE_TREE: ""},
            detached=False)
    nested_path = source_root + "/sse2neon"
    nested = rich_git_source_fixture(
        nested_path, SSE2NEON_COMMIT,
        base64.b64decode(SSE2NEON_COMMIT_BASE64, validate=True),
        SSE2NEON_TREE_OBJECTS_BASE64, detached=True)
    return rich_git_source_fixture(
        source_root, head, raw, BASELINE_TREE_OBJECTS_BASE64,
        detached=True, nested={"sse2neon": nested})


def cmake_fixture_identity(
    raw_schema: str, *, multi_config: bool = False,
) -> tuple[dict, dict]:
    if raw_schema in runner.COMPLETE_EVIDENCE_SCHEMAS:
        baseline_build = complete_build_fixture(
            "baseline", raw_schema, multi_config=multi_config)
        candidate_build = complete_build_fixture(
            "candidate", raw_schema, multi_config=multi_config)
        baseline_executable = baseline_build["validated_executable"]
        candidate_executable = candidate_build["validated_executable"]
        baseline_archive = baseline_build["validated_archive"]
        candidate_archive = candidate_build["validated_archive"]
        compiler = baseline_build["compiler"]
        candidate_build["compiler"] = copy.deepcopy(compiler)
        candidate_build["compiler_version_stdout"] = copy.deepcopy(
            baseline_build["compiler_version_stdout"])
        identity = {
            "runner": complete_artifact(
                SPECIFICATION["runner"], "file", "1"),
            "taskset": complete_artifact(
                SPECIFICATION["taskset"], "executable", "2"),
            "ldd": complete_artifact(
                SPECIFICATION["ldd"], "executable", "3"),
            "baseline_executable": copy.deepcopy(baseline_executable),
            "candidate_executable": copy.deepcopy(candidate_executable),
            "baseline_archive": copy.deepcopy(baseline_archive),
            "candidate_archive": copy.deepcopy(candidate_archive),
            "baseline_build": baseline_build,
            "candidate_build": candidate_build,
            "baseline_runtime_closure": complete_runtime_fixture(
                baseline_executable["path"], "4", raw_schema),
            "candidate_runtime_closure": complete_runtime_fixture(
                candidate_executable["path"], "5", raw_schema),
            "baseline_source": complete_source_fixture(
                "baseline", raw_schema),
            "candidate_source": complete_source_fixture(
                "candidate", raw_schema),
        }
        if raw_schema in runner.BUILD_CLOSURE_V7_SCHEMAS:
            identity["evidence_helper"] = complete_artifact(
                SPECIFICATION["candidate_source_root"] + "/" +
                runner.EVIDENCE_HELPER_RELATIVE_PATH, "file", "6")
        specification = copy.deepcopy(SPECIFICATION)
        specification.update({
            "baseline_executable": baseline_executable["path"],
            "candidate_executable": candidate_executable["path"],
            "baseline_archive": baseline_archive["path"],
            "candidate_archive": candidate_archive["path"],
        })
        return identity, specification

    cmake = runner.cmake_identity_for_raw_schema(raw_schema)
    archive_path = f"/fixture/candidate-build/{cmake['archive']}"
    archive = {"path": archive_path, "sha256": "a" * 64}
    recipe_text = archive_recipe_fixture_text(cmake)
    recipe_content = runner.exact_text_content(
        recipe_text, "fixture archive link recipe")
    recipe = {
        "path": "/fixture/candidate-build/CMakeFiles/"
                f"{cmake['target_directory']}/link.txt",
        "size": recipe_content["size"],
        "sha256": recipe_content["sha256"],
    }
    build = {
        "build_dir": "/fixture/candidate-build",
        "validated_archive": copy.deepcopy(archive),
        "archive_link_recipe": recipe,
        "archiver": {"path": "/usr/bin/ar"},
        "validated_archive_members": [
            "LeopardCommon.cpp.o", "Leopard2BackendAVX2.cpp.o"],
        "validated_compile_commands": {
            "required_source_object_pairs": [
                {
                    "source": {"path": "/fixture/source/LeopardCommon.cpp"},
                    "object": {"path": "/fixture/candidate-build/CMakeFiles/"
                               f"{cmake['target_directory']}/LeopardCommon.cpp.o"},
                },
                {
                    "source": {"path": "/fixture/source/Leopard2BackendAVX2.cpp"},
                    "object": {"path": "/fixture/candidate-build/CMakeFiles/"
                               "leopard2_backend_avx2.dir/"
                               "Leopard2BackendAVX2.cpp.o"},
                },
            ],
        },
    }
    if raw_schema in runner.HARDENED_BUILD_SCHEMAS:
        build["archive_link_recipe_content"] = recipe_content
        build["ranlib"] = {"path": "/usr/bin/ranlib"}
        build["archive_link_tool_invocations"] = {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                       "resolved_path": "/usr/bin/ranlib"},
        }
    identity = {
        "fixture": {"sha256": "a" * 64},
        "candidate_archive": copy.deepcopy(archive),
        "candidate_build": build,
    }
    specification = copy.deepcopy(SPECIFICATION)
    specification["candidate_archive"] = archive_path
    return identity, specification


def sealed_executable_fixtures(identity: Mapping[str, object]) -> dict:
    result = {}
    for role in ("baseline", "candidate"):
        source = copy.deepcopy(identity[f"{role}_executable"])
        closure = copy.deepcopy(identity[f"{role}_runtime_closure"])
        closure["executable"] = runner.SEALED_EXECUTABLE_COMMAND[role]
        result[role] = {
            "source": source,
            "snapshot": {
                "protocol": runner.SEALED_EXECUTABLE_PROTOCOL,
                "size": source["size"],
                "mode": source["mode"],
                "sha256": source["sha256"],
                "seals": runner.LINUX_REQUIRED_EXECUTABLE_SEALS,
                "elf": True,
            },
            "runtime_closure": closure,
        }
    return result


def synthetic_raw(
    candidate_scale: float = 0.8, raw_schema: str = runner.RAW_SCHEMA,
    candidate_mode: str = "auto",
    *, multi_config: bool = False,
) -> dict:
    identity, specification = cmake_fixture_identity(
        raw_schema, multi_config=multi_config)
    executable_snapshots = (
        sealed_executable_fixtures(identity)
        if raw_schema in runner.SEALED_EXECUTABLE_SCHEMAS else None)
    campaign = copy.deepcopy(CAMPAIGN)
    if raw_schema in runner.CANDIDATE_MODE_SCHEMAS:
        campaign["candidate_mode"] = candidate_mode
    else:
        campaign.pop("candidate_mode", None)
    invocations = []
    for round_index in range(runner.ROUNDS):
        for slot, implementation in enumerate(runner.ORDER):
            result = (baseline_result() if implementation == "baseline"
                      else candidate_result(
                          candidate_scale, raw_schema, campaign))
            normalized = runner.validate_result(
                implementation, result, CELL, campaign, raw_schema)
            invocations.append({
                "cell_id": CELL.identifier,
                "round": round_index,
                "slot": slot,
                "implementation": implementation,
                "command": [
                    specification["taskset"], "-c", "0",
                    *runner.benchmark_arguments(
                        implementation,
                        (Path(runner.SEALED_EXECUTABLE_COMMAND[implementation])
                         if executable_snapshots is not None else
                         Path(specification[f"{implementation}_executable"])),
                        CELL, campaign),
                ],
                "environment": copy.deepcopy(runner.CHILD_ENVIRONMENT),
                "pinned_cpu": 0,
                "started_utc": "2026-07-16T00:00:00Z",
                "duration_ns": 1,
                "returncode": 0,
                "stdout": {"path": "unused", "size": 0, "sha256": "0" * 64},
                "stderr": {"path": "unused", "size": 0, "sha256": "0" * 64},
                "result": result,
                "normalized": normalized,
                "identity_before": identity,
                "identity_after": identity,
                "reservation_before": RESERVATION,
                "reservation_after": RESERVATION,
            })
            if executable_snapshots is not None:
                invocations[-1]["execution_protocol"] = \
                    runner.SEALED_EXECUTABLE_PROTOCOL
                invocations[-1]["executable_snapshot"] = copy.deepcopy(
                    executable_snapshots[implementation])
    analysis = runner.analyze(invocations, campaign)
    payload = {
        "schema": raw_schema,
        "created_utc": "2026-07-16T00:00:00Z",
        "validity_is_independent_of_speed": True,
        "campaign": campaign,
        "host_initial": copy.deepcopy(HOST),
        "isolation": copy.deepcopy(ISOLATION),
        "reservation": RESERVATION,
        "input_specification": specification,
        "identities_initial": identity,
        "invocations": invocations,
        "identities_final": identity,
        "host_final": copy.deepcopy(HOST),
        "analysis": analysis,
    }
    if raw_schema in runner.SUPERVISION_SCHEMAS:
        payload["supervision"] = runner.supervision_record(
            "ab" * 32, 900, 2_100, campaign, RESERVATION, ISOLATION)
    if executable_snapshots is not None:
        payload["executable_snapshots"] = executable_snapshots
    return runner.signed(payload)


def resign(value: dict) -> dict:
    payload = copy.deepcopy(value)
    payload.pop("digest", None)
    return runner.signed(payload)


def synchronize_identity(value: dict) -> None:
    value["identities_final"] = copy.deepcopy(value["identities_initial"])
    for invocation in value["invocations"]:
        invocation["identity_before"] = copy.deepcopy(value["identities_initial"])
        invocation["identity_after"] = copy.deepcopy(value["identities_initial"])


def write_complete_evidence_bundle(
    root: Path, value: dict, manifest_schema: str = runner.MANIFEST_SCHEMA,
) -> Path:
    for index, invocation in enumerate(value["invocations"]):
        stdout = json.dumps(invocation["result"]).encode("utf-8")
        stdout_path = root / f"{index}.stdout"
        stderr_path = root / f"{index}.stderr"
        stdout_path.write_bytes(stdout)
        stderr_path.write_bytes(b"")
        if manifest_schema in (
                runner.MANIFEST_SCHEMA_V7, runner.MANIFEST_SCHEMA_V8,
                runner.MANIFEST_SCHEMA):
            stdout_path.chmod(0o600)
            stderr_path.chmod(0o600)
        invocation["stdout"] = {
            "path": stdout_path.name,
            "size": len(stdout),
            "sha256": runner.sha256_bytes(stdout),
        }
        invocation["stderr"] = {
            "path": stderr_path.name,
            "size": 0,
            "sha256": runner.sha256_bytes(b""),
        }
    value = resign(value)
    raw_path = root / "raw.json"
    runner.write_json_exclusive(raw_path, value)
    manifest_payload = {
        "schema": manifest_schema,
        "created_utc": "2026-07-16T00:00:00Z",
        "valid": True,
        "validity_is_independent_of_speed": True,
        "raw": {
            "path": raw_path.name,
            "size": raw_path.stat().st_size,
            "sha256": runner.sha256_file(raw_path),
            "payload_digest": value["digest"],
        },
        "campaign": value["campaign"],
        "host": value["host_initial"],
        "isolation": value["isolation"],
        "reservation": value["reservation"],
        "identities": value["identities_initial"],
        "analysis": value["analysis"],
    }
    if manifest_schema in (
        runner.MANIFEST_SCHEMA_V5, runner.MANIFEST_SCHEMA_V6,
        runner.MANIFEST_SCHEMA_V7, runner.MANIFEST_SCHEMA_V8,
        runner.MANIFEST_SCHEMA,
    ):
        manifest_payload["supervision"] = value["supervision"]
    if manifest_schema in (
            runner.MANIFEST_SCHEMA_V8, runner.MANIFEST_SCHEMA):
        manifest_payload["executable_snapshots"] = \
            value["executable_snapshots"]
    manifest = runner.signed(manifest_payload)
    manifest_path = root / "manifest.json"
    runner.write_json_exclusive(manifest_path, manifest)
    return manifest_path


def recursively_replace_strings(value: object, replacements: tuple[tuple[str, str], ...]
                                ) -> object:
    if isinstance(value, str):
        result = value
        for before, after in replacements:
            result = result.replace(before, after)
        return result
    if isinstance(value, list):
        return [recursively_replace_strings(item, replacements) for item in value]
    if isinstance(value, dict):
        return {
            key: recursively_replace_strings(item, replacements)
            for key, item in value.items()
        }
    return value


def attach_recipe_content(value: object, content: dict) -> None:
    if isinstance(value, list):
        for item in value:
            attach_recipe_content(item, content)
        return
    if not isinstance(value, dict):
        return
    build = value.get("candidate_build")
    if isinstance(build, dict) and isinstance(build.get("archive_link_recipe"), dict):
        build["archive_link_recipe"]["size"] = content["size"]
        build["archive_link_recipe"]["sha256"] = content["sha256"]
        build["archive_link_recipe_content"] = copy.deepcopy(content)
        build["ranlib"] = {"path": "/usr/bin/ranlib"}
        build["archive_link_tool_invocations"] = {
            "archiver": {"invocation": "/usr/bin/ar",
                         "resolved_path": "/usr/bin/ar"},
            "ranlib": {"invocation": "/usr/bin/ranlib",
                       "resolved_path": "/usr/bin/ranlib"},
        }
    for item in value.values():
        attach_recipe_content(item, content)


def replace_current_recipe_text(value: dict, text: str) -> None:
    content = runner.exact_text_content(text, "mutated archive link recipe")
    build = value["identities_initial"]["candidate_build"]
    build["archive_link_recipe_content"] = content
    build["archive_link_recipe"]["size"] = content["size"]
    build["archive_link_recipe"]["sha256"] = content["sha256"]
    value["identities_final"] = copy.deepcopy(value["identities_initial"])
    for invocation in value["invocations"]:
        invocation["identity_before"] = copy.deepcopy(value["identities_initial"])
        invocation["identity_after"] = copy.deepcopy(value["identities_initial"])


def replace_current_executable_recipe_text(value: dict, text: str) -> None:
    content = runner.exact_text_content(
        text, "mutated executable link recipe")
    build = value["identities_initial"]["candidate_build"]
    build["executable_link_recipe_content"] = content
    build["executable_link_recipe"]["size"] = content["size"]
    build["executable_link_recipe"]["sha256"] = content["sha256"]
    synchronize_identity(value)


def synthetic_failure(raw_schema: str) -> dict:
    raw = synthetic_raw(raw_schema=raw_schema)
    failure_schema = {
        runner.RAW_SCHEMA_V2: runner.FAILURE_SCHEMA_V2,
        runner.RAW_SCHEMA_V3: runner.FAILURE_SCHEMA_V3,
        runner.RAW_SCHEMA_V4: runner.FAILURE_SCHEMA_V4,
        runner.RAW_SCHEMA_V5: runner.FAILURE_SCHEMA_V5,
        runner.RAW_SCHEMA_V6: runner.FAILURE_SCHEMA_V6,
        runner.RAW_SCHEMA_V7: runner.FAILURE_SCHEMA_V7,
        runner.RAW_SCHEMA_V8: runner.FAILURE_SCHEMA_V8,
        runner.RAW_SCHEMA: runner.FAILURE_SCHEMA,
    }[raw_schema]
    payload = {
        "schema": failure_schema,
        "created_utc": "2026-07-16T00:00:00Z",
        "status": "failed",
        "valid": False,
        "error_type": "EvidenceError",
        "error": "fixture failure",
        "campaign": copy.deepcopy(raw["campaign"]),
        "host_initial": copy.deepcopy(raw["host_initial"]),
        "reservation": copy.deepcopy(raw["reservation"]),
        "pair_lease": copy.deepcopy(PAIR_LEASE),
        "isolation": copy.deepcopy(raw["isolation"]),
        "input_specification": copy.deepcopy(raw["input_specification"]),
        "identities_initial": copy.deepcopy(raw["identities_initial"]),
        "invocations": [],
        "retained_files": [],
        "traceback": "fixture traceback",
    }
    if raw_schema == runner.RAW_SCHEMA:
        payload["evidence_contract"] = runner.FAILURE_EVIDENCE_CONTRACT
    if raw_schema in runner.SUPERVISION_SCHEMAS:
        payload["supervision"] = copy.deepcopy(raw["supervision"])
    if raw_schema in runner.SEALED_EXECUTABLE_SCHEMAS:
        payload["executable_snapshots"] = copy.deepcopy(
            raw["executable_snapshots"])
    return runner.signed(payload)


class MainCompareRunnerTests(unittest.TestCase):
    def assert_rejected(self, value: dict) -> None:
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(
                resign(value), None, check_files=False, check_current_inputs=False)

    def test_direct_script_help_loads_shared_git_capture(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-direct-script-") as directory:
            completed = subprocess.run(
                [sys.executable, str(MODULE_PATH), "--help"],
                cwd=directory, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={**os.environ, "PYTHONWARNINGS":
                     "error::ResourceWarning"},
                timeout=30, check=False)
        self.assertEqual(
            completed.returncode, 0,
            completed.stderr.decode("utf-8", errors="replace"))

    def test_exact_main_pure_avx2_compile_profile(self) -> None:
        specification = dict(SPECIFICATION)
        specification["baseline_pure_avx2"] = True
        source = Path(specification["baseline_source_root"]) / \
            "LeopardFF8.cpp"
        arguments = runner.expected_compile_argv(
            "baseline", source, specification, "/usr/bin/c++")
        first = arguments.index("-march=x86-64")
        self.assertEqual(arguments[first:first + 4], [
            "-march=x86-64", "-mtune=generic", "-mavx2",
            "-mno-avx512f",
        ])
        self.assertNotIn("-march=native", arguments)
        self.assertEqual(
            runner.baseline_compile_profile(specification),
            runner.BASELINE_PURE_AVX2_COMPILE_PROFILE)
        self.assertEqual(
            runner.baseline_isa_policy(specification),
            "whole-build -march=x86-64 -mtune=generic -mavx2 "
            "-mno-avx512f")
        specification["baseline_pure_avx2"] = 1
        with self.assertRaises(runner.EvidenceError):
            runner.expected_compile_argv(
                "baseline", source, specification, "/usr/bin/c++")

    @unittest.skipUnless(Path("/proc/self/fd").is_dir(),
                         "descriptor-count regression needs procfs")
    def test_git_symlink_type_failure_closes_untransferred_descriptor(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-symlink-fd-") as directory, \
                ExitStack() as stack:
            root = Path(directory)
            (root / "victim").write_text("regular", encoding="ascii")
            tree = runner.git_capture._OpenDirectoryTree(
                root, stack, "symlink descriptor fixture")
            before = len(os.listdir("/proc/self/fd"))
            with mock.patch.object(
                    runner.git_capture.os, "readlink",
                    return_value="forged-target"), \
                 self.assertRaisesRegex(
                        runner.git_capture.GitCaptureError,
                        "not a symlink"):
                runner.git_capture._open_symlink(tree, "victim", stack)
            self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    @unittest.skipUnless(Path("/proc/self/fd").is_dir(),
                         "descriptor-count regression needs procfs")
    def test_git_capture_symlink_return_trace_cannot_strand_descriptor(
            self) -> None:
        class InjectedTraceFailure(BaseException):
            pass

        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-symlink-trace-") as directory:
            root = Path(directory)
            initialize_git_fixture(root)
            (root / "tracked-link").symlink_to("tracked.txt")
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "add", "tracked-link"],
                check=True)
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "commit", "-qm",
                 "tracked symlink"], check=True)
            commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
                text=True).strip()
            before = len(os.listdir("/proc/self/fd"))
            triggered = False

            def fail_at_symlink_return(frame, event, argument):
                nonlocal triggered
                del argument
                if event == "return" and frame.f_code is \
                        runner.git_capture._open_symlink.__code__ and \
                        frame.f_locals.get("relative") == "tracked-link":
                    triggered = True
                    sys.settrace(None)
                    raise InjectedTraceFailure(
                        "injected symlink-return interruption")
                return fail_at_symlink_return

            sys.settrace(fail_at_symlink_return)
            try:
                with self.assertRaises(InjectedTraceFailure):
                    runner.git_identity(
                        root, commit, False, include_commit_object=True,
                        rich=True)
            finally:
                sys.settrace(None)
            self.assertTrue(triggered)
            self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_git_capture_accepts_detached_main_and_attached_candidate(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-clean-") as directory:
            root = Path(directory)
            main = root / "main"
            candidate = root / "candidate"
            main.mkdir()
            candidate.mkdir()
            main_commit = initialize_git_fixture(main, "main\n")
            candidate_commit = initialize_git_fixture(
                candidate, "candidate\n")
            subprocess.run(
                ["/usr/bin/git", "-C", str(main), "checkout", "-q",
                 "--detach", main_commit], check=True)
            main_identity = runner.git_identity(
                main, main_commit, True, include_commit_object=True,
                rich=True)
            candidate_identity = runner.git_identity(
                candidate, candidate_commit, False,
                include_commit_object=True, rich=True)
            self.assertTrue(main_identity["detached"])
            self.assertFalse(candidate_identity["detached"])
            self.assertEqual(
                main_identity["git_executable"],
                candidate_identity["git_executable"])
            runner.git_capture.validate_git_capture(
                main_identity, str(main.resolve()), main_commit,
                require_detached=True)
            runner.git_capture.validate_git_capture(
                candidate_identity, str(candidate.resolve()),
                candidate_commit, require_detached=False)

    def test_git_capture_tree_closure_rejects_coherent_transcript_forgery(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-tree-proof-") as directory:
            root = Path(directory)
            (root / "nested").mkdir()
            commit = initialize_git_fixture(root)
            (root / "nested" / "second.txt").write_text(
                "second\n", encoding="utf-8")
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "add", "."], check=True)
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "commit", "-qm",
                 "nested"], check=True)
            commit = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"],
                text=True).strip()
            identity = runner.git_identity(
                root, commit, False, include_commit_object=True, rich=True)
            self.assertGreater(len(identity["tree_objects"]), 1)
            self.assertTrue(all(
                set(record) == {
                    "path", "git_mode", "git_type", "object_id", "kind"}
                for record in identity["tracked_files"]))

            def reseal_transcripts(value: dict) -> None:
                records = value["tracked_files"]
                value["tracked_files_sha256"] = \
                    runner.git_capture._digest(records)
                tree_listing = b"".join(
                    (
                        f"{record['git_mode']} {record['git_type']} "
                        f"{record['object_id']}\t{record['path']}\0"
                    ).encode("utf-8")
                    for record in records)
                index_stage = b"".join(
                    (
                        f"{record['git_mode']} {record['object_id']} 0\t"
                        f"{record['path']}\0"
                    ).encode("utf-8")
                    for record in records)
                default_flags = b"".join(
                    f"H {record['path']}\0".encode("utf-8")
                    for record in records)
                value["tracked_tree_listing_sha256"] = \
                    hashlib.sha256(tree_listing).hexdigest()
                value["index"]["entry_count"] = len(records)
                value["index"]["stage"] = \
                    runner.git_capture._byte_identity(index_stage)
                value["index"]["flags_v"] = \
                    runner.git_capture._byte_identity(default_flags)
                value["index"]["flags_f"] = \
                    runner.git_capture._byte_identity(default_flags)

            forged_record = copy.deepcopy(identity)
            forged_record["tracked_files"][0]["object_id"] = "0" * 40
            reseal_transcripts(forged_record)
            with self.assertRaisesRegex(
                    runner.git_capture.GitCaptureError,
                    "differs from recursive tree objects"):
                runner.git_capture.validate_git_capture(
                    forged_record, str(root.resolve()), commit,
                    require_detached=False)

            omitted_record = copy.deepcopy(identity)
            omitted_record["tracked_files"].pop()
            reseal_transcripts(omitted_record)
            with self.assertRaises(runner.git_capture.GitCaptureError):
                runner.git_capture.validate_git_capture(
                    omitted_record, str(root.resolve()), commit,
                    require_detached=False)

            forged_tree = copy.deepcopy(identity)
            subtree = next(
                record for record in forged_tree["tree_objects"]
                if record["object_id"] != forged_tree["tree"])
            content = bytearray(base64.b64decode(
                subtree["base64"], validate=True))
            content[-1] ^= 1
            replacement = git_object_fixture("tree", bytes(content))
            forged_tree["tree_objects"][
                forged_tree["tree_objects"].index(subtree)] = replacement
            forged_tree["tree_objects"].sort(
                key=lambda record: record["object_id"])
            with self.assertRaises(runner.git_capture.GitCaptureError):
                runner.git_capture.validate_git_capture(
                    forged_tree, str(root.resolve()), commit,
                    require_detached=False)

            extra_field = copy.deepcopy(identity)
            extra_field["tracked_files"][0]["sha256"] = "f" * 64
            reseal_transcripts(extra_field)
            with self.assertRaisesRegex(
                    runner.git_capture.GitCaptureError,
                    "record shape differs"):
                runner.git_capture.validate_git_capture(
                    extra_field, str(root.resolve()), commit,
                    require_detached=False)

    def test_git_tree_object_rejects_blob_tree_duplicate_basename(self) -> None:
        raw = (
            b"100644 foo\0" + b"\x11" * 20 +
            b"40000 foo\0" + b"\x22" * 20
        )
        with self.assertRaisesRegex(
                runner.git_capture.GitCaptureError,
                "duplicate entry names"):
            runner.git_capture._parse_tree_object(raw, "fixture tree")

    def test_git_capture_rejects_deterministic_mixed_command_state(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-mixed-") as directory:
            root = Path(directory)
            first = initialize_git_fixture(root)
            (root / "tracked.txt").write_text("two\n", encoding="utf-8")
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "commit", "-qam",
                 "second"], check=True)
            alternate_tree = subprocess.check_output(
                ["/usr/bin/git", "-C", str(root), "rev-parse",
                 "HEAD^{tree}"], text=True).strip()
            subprocess.run(
                ["/usr/bin/git", "-C", str(root), "reset", "--hard", "-q",
                 first], check=True)
            invoke = runner.git_capture._invoke_git

            def mixed(*args, **kwargs):
                arguments = args[3]
                if tuple(arguments) == (
                        "rev-parse", "--verify", "HEAD^{tree}"):
                    return (alternate_tree + "\n").encode("ascii")
                return invoke(*args, **kwargs)

            with mock.patch.object(
                    runner.git_capture, "_invoke_git", side_effect=mixed), \
                    self.assertRaises(runner.EvidenceError):
                runner.git_identity(
                    root, first, False, include_commit_object=True, rich=True)

    def test_git_capture_rejects_actual_git_directory_aba(self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leo2-main-git-aba-") as directory:
            parent = Path(directory)
            root = parent / "source"
            root.mkdir()
            commit = initialize_git_fixture(root)
            alternate = parent / "alternate.git"
            saved = parent / "saved.git"
            shutil.copytree(root / ".git", alternate)
            invoke = runner.git_capture._invoke_git
            triggered = False

            def aba(*args, **kwargs):
                nonlocal triggered
                output = invoke(*args, **kwargs)
                if not triggered:
                    triggered = True
                    (root / ".git").rename(saved)
                    alternate.rename(root / ".git")
                    (root / ".git").rename(alternate)
                    saved.rename(root / ".git")
                return output

            with mock.patch.object(
                    runner.git_capture, "_invoke_git", side_effect=aba), \
                    self.assertRaises(runner.EvidenceError):
                runner.git_identity(
                    root, commit, False,
                    include_commit_object=True, rich=True)
            self.assertTrue(triggered)

    def test_valid_fixture_and_analysis(self) -> None:
        value = synthetic_raw()
        analysis = runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        for metric in ("encode", "decode_first_use", "decode_reuse_amortized"):
            result = analysis[CELL.identifier][metric]
            self.assertEqual(result["independent_round_count"], 3)
            self.assertEqual(result["constituent_pair_count"], 6)
            self.assertEqual(result["degrees_of_freedom"], 2)
            self.assertTrue(math.isfinite(result["ci95_lower"]))
            self.assertTrue(result["performance_result_does_not_affect_evidence_validity"])

    def test_v8_exact_shapes_and_timestamps_fail_closed(self) -> None:
        mutations = []

        def add_raw_field(value: dict) -> None:
            value["unbound"] = True

        def remove_raw_time(value: dict) -> None:
            value.pop("created_utc")

        def noncanonical_raw_time(value: dict) -> None:
            value["created_utc"] = "2026-07-16T00:00:00.1Z"

        def invalid_raw_time(value: dict) -> None:
            value["created_utc"] = "2026-02-30T00:00:00.000000Z"

        def add_invocation_field(value: dict) -> None:
            value["invocations"][0]["unbound"] = True

        def add_campaign_field(value: dict) -> None:
            value["campaign"]["unversioned_claim"] = {"speedup": 999}

        def add_reservation_field_everywhere(value: dict) -> None:
            reservation = copy.deepcopy(value["reservation"])
            reservation["unversioned_claim"] = {"lease": "forged"}
            value["reservation"] = copy.deepcopy(reservation)
            for invocation in value["invocations"]:
                invocation["reservation_before"] = copy.deepcopy(reservation)
                invocation["reservation_after"] = copy.deepcopy(reservation)

        def add_invocation_reservation_field(value: dict) -> None:
            reservation = copy.deepcopy(
                value["invocations"][0]["reservation_before"])
            reservation["unversioned_claim"] = {"lease": "forged"}
            value["invocations"][0]["reservation_before"] = reservation

        def remove_invocation_time(value: dict) -> None:
            value["invocations"][0].pop("started_utc")

        def noncanonical_invocation_time(value: dict) -> None:
            value["invocations"][0]["started_utc"] = \
                "2026-07-16T00:00:00+00:00"

        def add_stream_field(value: dict) -> None:
            value["invocations"][0]["stdout"]["unbound"] = True

        def remove_stream_digest(value: dict) -> None:
            value["invocations"][0]["stdout"].pop("sha256")

        mutations.extend((
            add_raw_field, remove_raw_time, noncanonical_raw_time,
            invalid_raw_time,
            add_invocation_field, add_campaign_field,
            add_reservation_field_everywhere,
            add_invocation_reservation_field, remove_invocation_time,
            noncanonical_invocation_time, add_stream_field,
            remove_stream_digest,
        ))
        for mutate in mutations:
            with self.subTest(mutation=mutate.__name__):
                value = synthetic_raw()
                mutate(value)
                self.assert_rejected(value)

        historical = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7)
        historical["unbound_historical_field"] = True
        historical["campaign"]["unbound_historical_field"] = True
        historical["invocations"][0]["unbound_historical_field"] = True
        supervision = historical["supervision"]
        historical["supervision"] = runner.supervision_record(
            supervision["execution_nonce"],
            supervision["runner_started_monotonic_ns"],
            supervision["runner_finished_monotonic_ns"],
            historical["campaign"], historical["reservation"],
            historical["isolation"])
        historical["supervision"]["runner_pid"] = supervision["runner_pid"]
        runner.validate_raw(
            resign(historical), None, check_files=False,
            check_current_inputs=False)

    def test_v8_manifest_and_raw_identity_shapes_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original_path = write_complete_evidence_bundle(
                root, synthetic_raw())
            original = json.loads(original_path.read_text(encoding="utf-8"))
            mutations = (
                ("manifest-extra",
                 lambda value: value.__setitem__("unbound", True)),
                ("manifest-missing-time",
                 lambda value: value.pop("created_utc")),
                ("manifest-noncanonical-time",
                 lambda value: value.__setitem__(
                     "created_utc", "2026-07-16T00:00:00.1Z")),
                ("manifest-invalid-time",
                 lambda value: value.__setitem__(
                     "created_utc", "2026-02-30T00:00:00.000000Z")),
                ("raw-info-extra",
                 lambda value: value["raw"].__setitem__("unbound", True)),
                ("raw-info-missing-digest",
                 lambda value: value["raw"].pop("payload_digest")),
            )
            for label, mutate in mutations:
                with self.subTest(label=label):
                    changed = copy.deepcopy(original)
                    mutate(changed)
                    changed = resign(changed)
                    original_path.unlink()
                    runner.write_json_exclusive(original_path, changed)
                    with self.assertRaises(runner.EvidenceError):
                        runner.verify_campaign(argparse.Namespace(
                            manifest=original_path,
                            no_current_input_check=True))
            original_path.unlink()
            runner.write_json_exclusive(original_path, original)

    def test_failure_input_specification_is_exact_and_never_keyerrors(
        self,
    ) -> None:
        keys = sorted(runner.INPUT_SPECIFICATION_KEYS)
        for key in keys:
            with self.subTest(missing=key):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["input_specification"].pop(key)
                with self.assertRaises(runner.EvidenceError) as captured:
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)
                self.assertNotIsInstance(captured.exception, KeyError)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["input_specification"]["unbound"] = "/tmp/unbound"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["created_utc"] = "2026-07-16T00:00:00.1Z"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["campaign"]["unversioned_claim"] = {"speedup": 999}
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["reservation"]["unversioned_claim"] = {"lease": "forged"}
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)

    def test_v8_failed_invocation_and_stream_shapes_fail_closed(self) -> None:
        raw = synthetic_raw()
        original_invocation = copy.deepcopy(raw["invocations"][0])

        def failure_with_invocation() -> dict:
            failure = synthetic_failure(runner.RAW_SCHEMA)
            invocation = copy.deepcopy(original_invocation)
            invocation["stdout"]["path"] = "unused.stdout"
            invocation["stderr"]["path"] = "unused.stderr"
            failure["invocations"] = [invocation]
            failure["retained_files"] = [
                copy.deepcopy(invocation["stdout"]),
                copy.deepcopy(invocation["stderr"]),
            ]
            return failure

        failure = failure_with_invocation()
        runner.validate_failure(
            resign(failure), Path("/unused"), check_files=False)

        failure = failure_with_invocation()
        failure["invocations"][0]["unbound"] = True
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)

        failure = failure_with_invocation()
        failure["invocations"][0]["stderr"].pop("sha256")
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)

    def test_count_and_cpu_fields_require_bounded_non_boolean_integers(self) -> None:
        for key, replacement in (
            ("rounds", True), ("batch", True), ("threads", True),
            ("iterations", True), ("iterations", runner.MAX_CAMPAIGN_COUNT + 1),
            ("reuse", True), ("reuse", runner.MAX_CAMPAIGN_COUNT + 1),
            ("warmup", True), ("warmup", runner.MAX_CAMPAIGN_COUNT + 1),
            ("benchmark_cpu", False),
            ("benchmark_cpu", runner.MAX_CPU_ID + 1),
            ("reserved_sibling", True),
            ("reserved_sibling", runner.MAX_CPU_ID + 1),
        ):
            with self.subTest(kind="raw-campaign", key=key, value=replacement):
                value = synthetic_raw()
                value["campaign"][key] = replacement
                self.assert_rejected(value)

        for key, replacement in (
            ("rounds", True), ("batch", True), ("threads", True),
            ("iterations", True), ("iterations", runner.MAX_CAMPAIGN_COUNT + 1),
            ("reuse", True), ("reuse", runner.MAX_CAMPAIGN_COUNT + 1),
            ("warmup", True), ("warmup", runner.MAX_CAMPAIGN_COUNT + 1),
            ("benchmark_cpu", False),
            ("benchmark_cpu", runner.MAX_CPU_ID + 1),
            ("reserved_sibling", True),
            ("reserved_sibling", runner.MAX_CPU_ID + 1),
        ):
            with self.subTest(kind="failed-campaign", key=key, value=replacement):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["campaign"][key] = replacement
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)

        for invalid in (False, runner.MAX_CPU_ID + 1):
            with self.subTest(kind="topology", value=invalid), \
                 self.assertRaises(runner.EvidenceError):
                runner.validate_topology(invalid, 1)
            with self.subTest(kind="reservation-request", value=invalid), \
                 self.assertRaises(runner.EvidenceError):
                runner.parse_reservation(
                    runner.canonical_bytes(RESERVATION_PAYLOAD), invalid, 1)

        for key in ("benchmark_cpu", "reserved_sibling"):
            for invalid in (False, runner.MAX_CPU_ID + 1):
                with self.subTest(kind="reservation-payload", key=key,
                                  value=invalid):
                    payload = copy.deepcopy(RESERVATION_PAYLOAD)
                    payload[key] = invalid
                    with self.assertRaises(runner.EvidenceError):
                        runner.parse_reservation(
                            runner.canonical_bytes(payload), 0, 1)

        for field in ("allowed_cpu_set_at_launch", "online_cpu_set"):
            for invalid in (False, runner.MAX_CPU_ID + 1):
                with self.subTest(kind="host", field=field, value=invalid):
                    host = copy.deepcopy(HOST)
                    host[field][0] = invalid
                    with self.assertRaises(runner.EvidenceError):
                        runner.validate_host_record(
                            host, 0, 1, CAMPAIGN["allowed_cpu_set_at_launch"],
                            runner.RAW_SCHEMA)

        for invalid in (False, runner.MAX_CPU_ID + 1):
            with self.subTest(kind="campaign-allowed", value=invalid):
                value = synthetic_raw()
                value["campaign"]["allowed_cpu_set_at_launch"][0] = invalid
                self.assert_rejected(value)
            with self.subTest(kind="failed-campaign-allowed", value=invalid):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["campaign"]["allowed_cpu_set_at_launch"][0] = invalid
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)
            with self.subTest(kind="failed-host-online", value=invalid):
                failure = synthetic_failure(runner.RAW_SCHEMA)
                failure["host_initial"]["online_cpu_set"][0] = invalid
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_failure(
                        resign(failure), Path("/unused"), check_files=False)

        for field, replacement in (
            ("k", True), ("r", True), ("shard_bytes", True),
            ("shard_bytes", runner.MAX_SHARD_BYTES + 64),
            ("losses", True), ("seed", True), ("seed", runner.MASK64 + 1),
        ):
            with self.subTest(kind="cell", field=field, value=replacement):
                values = asdict(CELL)
                values[field] = replacement
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_cell(runner.Cell(**values))

    def test_uniformly_incomplete_complete_identities_fail_offline_verify_and_promotion(
        self,
    ) -> None:
        plan = load_plan_runner()

        with tempfile.TemporaryDirectory() as directory:
            valid_manifest = write_complete_evidence_bundle(
                Path(directory),
                synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7),
                runner.MANIFEST_SCHEMA_V7)
            runner.verify_campaign(argparse.Namespace(
                manifest=valid_manifest, no_current_input_check=True))
            document, scope, _ = plan.verify_exact_manifest(valid_manifest)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA_V7)
            self.assertEqual(
                plan.validate_evidence_scope(scope)["schema"],
                plan.EVIDENCE_SCOPE_SCHEMA_V4)

        def missing_sources(value: dict) -> None:
            for role in ("baseline", "candidate"):
                value["identities_initial"][f"{role}_source"].pop("tree")
            synchronize_identity(value)

        def incomplete_runtime_files(value: dict) -> None:
            for role in ("baseline", "candidate"):
                dependency = value["identities_initial"][
                    f"{role}_runtime_closure"]["dependencies"][0]
                dependency["file"] = {"path": dependency["file"]["path"]}
            synchronize_identity(value)

        def path_only_compile_pairs(value: dict) -> None:
            for role in ("baseline", "candidate"):
                pairs = value["identities_initial"][f"{role}_build"][
                    "validated_compile_commands"]["required_source_object_pairs"]
                for pair in pairs:
                    pair["source"] = {"path": pair["source"]["path"]}
                    pair["object"] = {"path": pair["object"]["path"]}
            synchronize_identity(value)

        def empty_toolchain_records(value: dict) -> None:
            for role in ("baseline", "candidate"):
                build = value["identities_initial"][f"{role}_build"]
                for name in ("compiler", "cmake_cache", "compile_commands",
                             "executable_link_recipe", "archive_link_recipe"):
                    build[name] = {}
            synchronize_identity(value)

        def reduced_outputs(value: dict) -> None:
            identities = value["identities_initial"]
            for role in ("baseline", "candidate"):
                for output in ("archive", "executable"):
                    key = f"{role}_{output}"
                    old = identities[key]
                    reduced = {
                        "path": old["path"], "kind": old["kind"],
                        "sha256": old["sha256"],
                    }
                    identities[key] = copy.deepcopy(reduced)
                    identities[f"{role}_build"][f"validated_{output}"] = \
                        copy.deepcopy(reduced)
            synchronize_identity(value)

        def topology_only_host(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                record = value["host_initial"][name]
                value["host_initial"][name] = {
                    "cpu": record["cpu"], "topology": record["topology"]}
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def truncated_tu_closure(value: dict) -> None:
            for role, names in (("baseline", runner.BASELINE_LIBRARY_SOURCES),
                                ("candidate", runner.CANDIDATE_LIBRARY_SOURCES)):
                build = value["identities_initial"][f"{role}_build"]
                semantics = build["validated_compile_commands"]
                suffix = "/" + names[-1]
                pair = next(item for item in semantics[
                    "required_source_object_pairs"]
                    if item["source"]["path"].endswith(suffix))
                semantics["required_source_object_pairs"].remove(pair)
                semantics["required_sources"].remove(pair["source"]["path"])
                semantics["entry_count"] -= 1
                member = Path(pair["object"]["path"]).name
                build["validated_archive_members"].remove(member)
                relative = Path(pair["object"]["path"]).relative_to(
                    Path(build["build_dir"])).as_posix()
                text = build["archive_link_recipe_content"]["text"].replace(
                    " " + relative, "", 1)
                content = runner.exact_text_content(
                    text, f"truncated {role} archive recipe")
                build["archive_link_recipe_content"] = content
                build["archive_link_recipe"]["size"] = content["size"]
                build["archive_link_recipe"]["sha256"] = content["sha256"]
            synchronize_identity(value)

        def multiline_executable_recipe(value: dict) -> None:
            for role in ("baseline", "candidate"):
                build = value["identities_initial"][f"{role}_build"]
                text = build["executable_link_recipe_content"]["text"] + \
                    "\n-O3\n"
                content = runner.exact_text_content(
                    text, f"multiline {role} executable recipe")
                build["executable_link_recipe_content"] = content
                build["executable_link_recipe"]["size"] = content["size"]
                build["executable_link_recipe"]["sha256"] = content["sha256"]
            synchronize_identity(value)

        def missing_dynamic_loader(value: dict) -> None:
            for role in ("baseline", "candidate"):
                closure = value["identities_initial"][f"{role}_runtime_closure"]
                closure["dependencies"] = [item for item in closure["dependencies"]
                    if item["soname"] != "ld-linux-x86-64.so.2"]
                text = "\n".join(line for line in
                    retained_ldd_text(closure).splitlines()
                    if not line.startswith("/lib64/ld-linux")) + "\n"
                replace_retained_ldd_text(
                    closure, text, f"truncated {role} ldd output")
            synchronize_identity(value)

        def swapped_runtime_file_records(value: dict) -> None:
            for role in ("baseline", "candidate"):
                dependencies = value["identities_initial"][
                    f"{role}_runtime_closure"]["dependencies"]
                libc = next(item for item in dependencies
                            if item["soname"] == "libc.so.6")
                libm = next(item for item in dependencies
                            if item["soname"] == "libm.so.6")
                libc["file"], libm["file"] = libm["file"], libc["file"]
            synchronize_identity(value)

        def noncanonical_runtime_loader_paths(value: dict) -> None:
            for role in ("baseline", "candidate"):
                closure = value["identities_initial"][f"{role}_runtime_closure"]
                libc = next(item for item in closure["dependencies"]
                            if item["soname"] == "libc.so.6")
                libc["loader_path"] = "/lib/./libc.so.6"
                libc["file"]["path"] = libc["loader_path"]
                text = retained_ldd_text(closure).replace(
                    "/lib/libc.so.6", libc["loader_path"])
                replace_retained_ldd_text(
                    closure, text, f"noncanonical {role} ldd output")
            synchronize_identity(value)

        def truncated_cache_inventory(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                value["host_initial"][name]["cache_hierarchy"].pop()
                value["host_initial"][name]["cache_index_inventory"].pop()
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def empty_numa_summary(value: dict) -> None:
            for name in ("benchmark_cpu", "reserved_sibling"):
                value["host_initial"][name]["numa_nodes"] = []
                value["host_initial"][name]["numa_node_inventory"] = []
            value["host_final"] = copy.deepcopy(value["host_initial"])

        def source_tree_commit_mismatch(value: dict) -> None:
            for role in ("baseline", "candidate"):
                value["identities_initial"][f"{role}_source"]["tree"] = "f" * 40
            synchronize_identity(value)

        def boolean_campaign_counts(value: dict) -> None:
            value["campaign"]["iterations"] = True

        mutations = {
            "sources": missing_sources,
            "runtime-files": incomplete_runtime_files,
            "path-only-compile-pairs": path_only_compile_pairs,
            "empty-toolchain-records": empty_toolchain_records,
            "reduced-outputs": reduced_outputs,
            "topology-only-host": topology_only_host,
            "truncated-tu-closure": truncated_tu_closure,
            "multiline-executable-recipe": multiline_executable_recipe,
            "missing-dynamic-loader": missing_dynamic_loader,
            "swapped-runtime-file-records": swapped_runtime_file_records,
            "noncanonical-runtime-loader-paths":
                noncanonical_runtime_loader_paths,
            "truncated-cache-inventory": truncated_cache_inventory,
            "empty-numa-summary": empty_numa_summary,
            "source-tree-commit-mismatch": source_tree_commit_mismatch,
            "boolean-campaign-counts": boolean_campaign_counts,
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                value = synthetic_raw()
                mutate(value)
                manifest_path = write_complete_evidence_bundle(root, value)
                with self.assertRaises(runner.EvidenceError):
                    runner.verify_campaign(argparse.Namespace(
                        manifest=manifest_path, no_current_input_check=True))
                with self.assertRaises(plan.PlanError):
                    plan.verify_exact_manifest(manifest_path)

    def test_promotion_accepts_v9_and_rejects_schema_downgrades(self) -> None:
        plan = load_plan_runner()
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), synthetic_raw())
            document, scope, _ = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA)
            self.assertEqual(
                plan.validate_evidence_scope(scope)["schema"],
                plan.EVIDENCE_SCOPE_SCHEMA)
        for manifest_schema, raw_schema in (
            (runner.MANIFEST_SCHEMA_V8, runner.RAW_SCHEMA),
            (runner.MANIFEST_SCHEMA, runner.RAW_SCHEMA_V8),
            (runner.MANIFEST_SCHEMA_V7, runner.RAW_SCHEMA),
            (runner.MANIFEST_SCHEMA, runner.RAW_SCHEMA_V7),
            (runner.MANIFEST_SCHEMA_V6, runner.RAW_SCHEMA_V5),
        ):
            with self.subTest(
                    manifest_schema=manifest_schema,
                    raw_schema=raw_schema), self.assertRaises(plan.PlanError):
                plan._validate_exact_schema_pair(
                    {"schema": manifest_schema}, {"schema": raw_schema})

    def test_coherent_ordinary_runtime_rewrite_is_internal_consistency_only(
        self,
    ) -> None:
        plan = load_plan_runner()
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7)
        for role in ("baseline", "candidate"):
            closure = value["identities_initial"][f"{role}_runtime_closure"]
            closure["dependencies"] = [
                dependency for dependency in closure["dependencies"]
                if dependency["soname"] != "libm.so.6"]
            text = "\n".join(
                line for line in retained_ldd_text(closure).splitlines()
                if not line.startswith("libm.so.6 =>")) + "\n"
            replace_retained_ldd_text(
                closure, text, f"coherently rewritten {role} ldd output")
        synchronize_identity(value)
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), value, runner.MANIFEST_SCHEMA_V7)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)), 0)
            document, scope, _ = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA_V7)
            plan.validate_evidence_scope(scope)

    def test_coherent_cache_inventory_rewrite_is_internal_consistency_only(
        self,
    ) -> None:
        plan = load_plan_runner()
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7)
        for name in ("benchmark_cpu", "reserved_sibling"):
            record = value["host_initial"][name]
            record["cache_hierarchy"].pop()
            record["cache_index_inventory"].pop()
            text = "".join(
                f"{entry}\n" for entry in record["cache_index_inventory"])
            record["cache_directory_inventory_text"] = \
                runner.exact_text_content(
                    text, f"coherently rewritten {name} cache inventory")
        value["host_final"] = copy.deepcopy(value["host_initial"])
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), value, runner.MANIFEST_SCHEMA_V7)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)), 0)
            document, scope, _ = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA_V7)
            plan.validate_evidence_scope(scope)

    def test_promotion_consumes_the_exact_verified_manifest_and_raw_snapshots(
        self,
    ) -> None:
        plan = load_plan_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = write_complete_evidence_bundle(
                root, synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7),
                runner.MANIFEST_SCHEMA_V7)
            raw_path = root / "raw.json"
            exact_runner = plan.load_exact_main_runner()

            class InterleavingRunner:
                accepted = None

                def verified_campaign_bundle(
                    self, path: Path, no_current_input_check: bool = False,
                ):
                    self.accepted = exact_runner.verified_campaign_bundle(
                        path, no_current_input_check)
                    path.write_text("{}\n", encoding="utf-8")
                    raw_path.write_text("{}\n", encoding="utf-8")
                    return self.accepted

            interleaving = InterleavingRunner()
            with mock.patch.object(
                plan, "load_exact_main_runner", return_value=interleaving,
            ):
                document, scope, snapshot = plan.verify_exact_manifest(manifest_path)
            self.assertIsNotNone(interleaving.accepted)
            self.assertEqual(document, interleaving.accepted[0])
            self.assertEqual(snapshot, interleaving.accepted[3])
            plan.validate_evidence_scope(scope)
            self.assertEqual(manifest_path.read_text(encoding="utf-8"), "{}\n")
            self.assertEqual(raw_path.read_text(encoding="utf-8"), "{}\n")

    def test_survivor_selection_references_the_exact_verified_snapshot(
        self,
    ) -> None:
        plan = load_plan_runner()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_path = write_complete_evidence_bundle(
                root, synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7),
                runner.MANIFEST_SCHEMA_V7)
            exact_runner = plan.load_exact_main_runner()

            class InterleavingRunner:
                accepted = None

                def verified_campaign_bundle(
                    self, path: Path, no_current_input_check: bool = False,
                ):
                    self.accepted = exact_runner.verified_campaign_bundle(
                        path, no_current_input_check)
                    path.write_text("{}\n", encoding="utf-8")
                    return self.accepted

            interleaving = InterleavingRunner()
            captured = {}

            def capture_selection(plan_root, manifests, references, scopes):
                captured["manifest"] = manifests[0]
                captured["reference"] = references[0]
                captured["scope"] = scopes[0]
                return {"reference": references[0]}

            with mock.patch.object(
                    plan, "load_exact_main_runner", return_value=interleaving), \
                 mock.patch.object(
                    plan, "derive_survivors", side_effect=capture_selection):
                plan.select_survivors(
                    root, [manifest_path], root / "survivors.json")

            self.assertIsNotNone(interleaving.accepted)
            self.assertEqual(captured["manifest"], interleaving.accepted[0])
            self.assertEqual(
                {key: captured["reference"][key]
                 for key in ("size", "sha256")},
                interleaving.accepted[3])
            self.assertEqual(
                captured["reference"]["payload_digest"],
                interleaving.accepted[0]["digest"])
            self.assertNotEqual(
                captured["reference"]["sha256"],
                runner.sha256_file(manifest_path))
            plan.validate_evidence_scope(captured["scope"])

    def test_legacy_v1_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V1)
        value.pop("isolation")
        value = resign(value)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v2_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v3_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V3)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v4_raw_fixture_remains_replayable(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V4)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_legacy_v4_manifest_fixture_remains_replayable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = write_complete_evidence_bundle(
                Path(directory), synthetic_raw(raw_schema=runner.RAW_SCHEMA_V4),
                runner.MANIFEST_SCHEMA_V4)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest, no_current_input_check=True)), 0)

    def test_legacy_v5_raw_and_manifest_fixtures_remain_replayable(self) -> None:
        plan = load_plan_runner()
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V5)
        entries = value["identities_initial"]["candidate_build"][
            "validated_compile_commands"]["required_entries"]
        benchmark_entry = next(
            entry for entry in entries
            if entry["file"].endswith("/bench/leopard2/benchmark.cpp"))
        self.assertIn(
            f'-DLEO2_BENCHMARK_SOURCE_COMMIT="{CANDIDATE_COMMIT}"',
            benchmark_entry["arguments"])
        self.assertIn(
            f'-DLEO2_BENCHMARK_SOURCE_TREE="{CANDIDATE_TREE}"',
            benchmark_entry["arguments"])
        self.assertIn(
            "-DLEO2_BENCHMARK_SOURCE_TRACKED_DIRTY=0",
            benchmark_entry["arguments"])
        self.assertFalse(any(
            "LEO2_BENCHMARK_SOURCE_ATTESTATION_HEADER" in token or
            "generated/leopard2-benchmark-attestation" in token
            for token in benchmark_entry["arguments"]))
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        with tempfile.TemporaryDirectory() as directory:
            manifest = write_complete_evidence_bundle(
                Path(directory), value, runner.MANIFEST_SCHEMA_V5)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest, no_current_input_check=True)), 0)
            document, scope, _ = plan.verify_exact_manifest(manifest)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA_V5)
            plan.validate_evidence_scope(scope)

    def test_legacy_v5_current_input_comparison_ignores_only_addresses(self) -> None:
        retained = synthetic_raw(
            raw_schema=runner.RAW_SCHEMA_V5)["identities_initial"]
        current = copy.deepcopy(retained)
        for role in ("baseline", "candidate"):
            closure = current[f"{role}_runtime_closure"]
            text = closure["raw_ldd_output"]["text"].replace(
                "0x00000000", "0xabcdef12")
            closure["raw_ldd_output"] = runner.exact_text_content(
                text, f"changed-address {role} raw ldd output")
        self.assertNotEqual(current, retained)
        self.assertTrue(runner.input_snapshots_equal(
            current, retained, runner.RAW_SCHEMA_V5))

        changed_path = copy.deepcopy(current)
        closure = changed_path["candidate_runtime_closure"]
        text = closure["raw_ldd_output"]["text"].replace(
            "/lib/libc.so.6", "/lib/libc-forged.so.6")
        closure["raw_ldd_output"] = runner.exact_text_content(
            text, "changed-path raw ldd output")
        self.assertFalse(runner.input_snapshots_equal(
            changed_path, retained, runner.RAW_SCHEMA_V5))

    @unittest.skipUnless(
        sys.platform.startswith("linux") and Path("/usr/bin/ldd").is_file() and
        Path("/usr/bin/env").is_file(), "requires Linux /usr/bin/ldd and env")
    def test_unmocked_same_binary_runtime_snapshots_are_aslr_stable(self) -> None:
        executable = Path("/usr/bin/env")
        first = runner.runtime_closure(
            Path("/usr/bin/ldd"), executable, runner.RAW_SCHEMA)
        second = runner.runtime_closure(
            Path("/usr/bin/ldd"), executable, runner.RAW_SCHEMA)
        self.assertEqual(first, second)
        self.assertIn("canonical_ldd_output", first)
        self.assertNotIn("raw_ldd_output", first)
        transcript = first["canonical_ldd_output"]
        self.assertEqual(transcript["schema"], runner.CANONICAL_LDD_SCHEMA)
        self.assertIn(runner.CANONICAL_LDD_ADDRESS, transcript["text"])
        runner.validate_complete_runtime_closure(
            first, "unmocked /usr/bin/env", str(executable.resolve()),
            runner.RAW_SCHEMA)

    def test_canonical_ldd_normalization_changes_only_terminal_addresses(
        self,
    ) -> None:
        first = (
            "linux-vdso.so.1 (0x1)\n"
            "libc.so.6 => /lib/libc.so.6 (0x2)\n"
            "/lib64/ld-linux-x86-64.so.2 (0x3)\n")
        second = first.replace("0x1", "0x111111").replace(
            "0x2", "0x222222").replace("0x3", "0x333333")
        self.assertEqual(
            runner.canonical_ldd_output(first, "first"),
            runner.canonical_ldd_output(second, "second"))
        changed_path = second.replace("/lib/libc.so.6", "/lib/libc-alt.so.6")
        self.assertNotEqual(
            runner.canonical_ldd_output(first, "first"),
            runner.canonical_ldd_output(changed_path, "changed path"))
        for label, malformed in {
            "missing address": first.replace(" (0x2)", "", 1),
            "address-to-path": first.replace("0x2", "0x2/forged", 1),
            "suffix after address": first.replace("(0x2)", "(0x2) /forged", 1),
        }.items():
            with self.subTest(label=label), self.assertRaises(runner.EvidenceError):
                runner.canonical_ldd_output(malformed, label)

    def test_canonical_ldd_transcript_and_summary_attacks_fail_closed(self) -> None:
        def reject(mutate) -> None:
            value = synthetic_raw()
            closure = value["identities_initial"]["candidate_runtime_closure"]
            mutate(closure)
            synchronize_identity(value)
            self.assert_rejected(value)

        def edit_text(closure: dict, transform) -> None:
            replace_retained_ldd_text(
                closure, transform(retained_ldd_text(closure)),
                "adversarial canonical ldd output")

        mutations = {
            "library-line-removal": lambda closure: edit_text(
                closure, lambda text: "\n".join(
                    line for line in text.splitlines()
                    if not line.startswith("libm.so.6 =>")) + "\n"),
            "loader-line-removal": lambda closure: edit_text(
                closure, lambda text: "\n".join(
                    line for line in text.splitlines()
                    if not line.startswith("/lib64/ld-linux")) + "\n"),
            "soname-relabel": lambda closure: edit_text(
                closure, lambda text: text.replace("libc.so.6 =>", "libx.so.6 =>")),
            "path-relabel": lambda closure: edit_text(
                closure, lambda text: text.replace(
                    "/lib/libc.so.6", "/lib/libc-forged.so.6")),
            "address-marker-removal": lambda closure: edit_text(
                closure, lambda text: text.replace(
                    f" ({runner.CANONICAL_LDD_ADDRESS})", "", 1)),
            "summary-removal": lambda closure: closure["dependencies"].pop(1),
            "normalization-relabel": lambda closure: closure[
                "canonical_ldd_output"].update({"normalization": "anything/v2"}),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                reject(mutate)

    def test_complete_schema_supervision_semantics_are_bound(self) -> None:
        value = synthetic_raw()
        for name, replacement in (
            ("execution_nonce", "bad"),
            ("launch_cpus", [0, 1]),
            ("reserved_cpus", [0, 2]),
            ("campaign_sha256", "f" * 64),
            ("reservation_sha256", "e" * 64),
            ("runner_started_monotonic_ns", 1_001),
            ("runner_finished_monotonic_ns", 1_999),
        ):
            edited = copy.deepcopy(value)
            edited["supervision"][name] = replacement
            self.assert_rejected(edited)
        edited = copy.deepcopy(value)
        edited["campaign"]["reuse"] += 1
        self.assert_rejected(edited)
        unsupervised = copy.deepcopy(value)
        unsupervised["supervision"] = None
        runner.validate_raw(
            resign(unsupervised), None, check_files=False,
            check_current_inputs=False)

    def test_candidate_modes_bind_flags_and_exact_argv(self) -> None:
        expected_arguments = {
            "auto": (),
            "generic": ("--force-generic",),
            "materialized": ("--force-specialized", "--force-materialized"),
            "tiled": ("--force-specialized", "--force-tiled"),
        }
        for mode in runner.CANDIDATE_MODES:
            with self.subTest(mode=mode):
                value = synthetic_raw(candidate_mode=mode)
                runner.validate_raw(
                    value, None, check_files=False,
                    check_current_inputs=False)
                flags = runner.candidate_mode_flags(mode)
                invocation = next(item for item in value["invocations"]
                                  if item["implementation"] == "candidate")
                parameters = invocation["result"]["parameters"]
                for name, expected in flags.items():
                    self.assertIs(parameters[name], expected)
                command = invocation["command"]
                for argument in expected_arguments[mode]:
                    self.assertIn(argument, command)

        value = synthetic_raw(candidate_mode="generic")
        value["campaign"]["candidate_mode"] = "auto"
        self.assert_rejected(value)

        value = synthetic_raw(candidate_mode="generic")
        invocation = next(item for item in value["invocations"]
                          if item["implementation"] == "candidate")
        invocation["command"].remove("--force-generic")
        self.assert_rejected(value)

        value = synthetic_raw(candidate_mode="tiled")
        invocation = next(item for item in value["invocations"]
                          if item["implementation"] == "candidate")
        invocation["result"]["parameters"]["force_tiled_decode"] = False
        self.assert_rejected(value)

        value = synthetic_raw()
        value["campaign"].pop("candidate_mode")
        self.assert_rejected(value)

        value = synthetic_raw()
        value["campaign"]["candidate_mode"] = ["auto"]
        self.assert_rejected(value)

        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V3)
        value["campaign"]["candidate_mode"] = "auto"
        self.assert_rejected(value)

    def test_workspace_selector_schema_boundary_is_fail_closed(self) -> None:
        for name in ("force_tiled_decode", "force_materialized_decode"):
            value = synthetic_raw()
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"].pop(name)
            self.assert_rejected(value)

            value = synthetic_raw()
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"][name] = True
            self.assert_rejected(value)

        for raw_schema in (runner.RAW_SCHEMA_V1, runner.RAW_SCHEMA_V2):
            value = synthetic_raw(raw_schema=raw_schema)
            for invocation in value["invocations"]:
                if invocation["implementation"] == "candidate":
                    invocation["result"]["parameters"].update({
                        "force_tiled_decode": False,
                        "force_materialized_decode": False,
                    })
            self.assert_rejected(value)

    def test_cmake_identity_and_cross_schema_relabels_are_rejected(self) -> None:
        value = synthetic_raw()
        value["input_specification"]["candidate_archive"] = \
            "/fixture/candidate-build/liblibleopard.a"
        self.assert_rejected(value)

        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe"]["path"] = \
            "/fixture/candidate-build/CMakeFiles/libleopard.dir/link.txt"
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["schema"] = runner.RAW_SCHEMA_V2
        self.assert_rejected(value)

        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        value["schema"] = runner.RAW_SCHEMA
        self.assert_rejected(value)

        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V5)
        value["schema"] = runner.RAW_SCHEMA
        self.assert_rejected(value)

        value = synthetic_raw()
        value["schema"] = runner.RAW_SCHEMA_V5
        self.assert_rejected(value)

    def test_private_hardened_historical_build_schema_is_not_evidence(self) -> None:
        private = runner.HARDENED_HISTORICAL_BUILD_SCHEMA
        self.assertIn(private, runner.BUILD_SCHEMA_TO_CMAKE_IDENTITY)
        self.assertIn(private, runner.HARDENED_BUILD_SCHEMAS)
        self.assertNotIn(private, runner.RAW_TO_CMAKE_IDENTITY)
        self.assertNotIn(private, runner.MANIFEST_TO_RAW_SCHEMA.values())
        self.assertNotIn(private, runner.FAILURE_TO_RAW_SCHEMA.values())
        self.assertEqual(
            runner.HISTORICAL_CMAKE_IDENTITY,
            runner.cmake_identity_for_build_schema(private))
        with self.assertRaises(runner.EvidenceError):
            runner.cmake_identity_for_raw_schema(private)

        value = synthetic_raw()
        value["schema"] = private
        self.assert_rejected(value)

    def test_coherent_historical_recipe_relabel_is_rejected(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V2)
        value = recursively_replace_strings(value, (
            ("liblibleopard.a", "libleopard.a"),
            ("CMakeFiles/libleopard.dir", "CMakeFiles/leopard.dir"),
        ))
        self.assertIsInstance(value, dict)
        value["schema"] = runner.RAW_SCHEMA
        historical = runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA_V2)
        old_content = runner.exact_text_content(
            archive_recipe_fixture_text(historical),
            "historical fixture archive link recipe")
        attach_recipe_content(value, old_content)
        # Every path and every surrounding digest can be coherently relabeled,
        # but the retained old recipe bytes still describe the old CMake target.
        self.assert_rejected(value)

    def test_recipe_content_and_identity_mutations_are_rejected(self) -> None:
        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe_content"]["text"] += "\n"
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

    def test_generated_attestation_header_is_a_canonical_compile_input(
        self,
    ) -> None:
        def retained(value: dict) -> dict:
            return value["identities_initial"]["candidate_build"][
                "validated_compile_commands"]["generated_attestation_header"]

        def benchmark_entry(value: dict) -> dict:
            entries = value["identities_initial"]["candidate_build"][
                "validated_compile_commands"]["required_entries"]
            return next(
                entry for entry in entries
                if entry["file"].endswith("/bench/leopard2/benchmark.cpp"))

        value = synthetic_raw()
        record = retained(value)
        injected = record["content"]["text"].replace(
            "#endif\n", '#include "stale_or_hostile.cpp"\n#endif\n')
        record["content"] = runner.exact_text_content(
            injected, "mutated generated attestation")
        record["artifact"]["size"] = record["content"]["size"]
        record["artifact"]["sha256"] = record["content"]["sha256"]
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        retained(value)["artifact"]["path"] = (
            SPECIFICATION["candidate_build_dir"] +
            "/foreign/leopard2_benchmark_source_attestation.h")
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        entry = benchmark_entry(value)
        expected_path = (
            SPECIFICATION["candidate_build_dir"] +
            "/generated/leopard2-benchmark-attestation/"
            "leopard2_benchmark_source_attestation.h")
        source_local_path = (
            SPECIFICATION["candidate_source_root"] +
            "/bench/leopard2/leopard2_benchmark_source_attestation.h")
        original = list(entry["arguments"])
        entry["arguments"] = [
            token.replace(expected_path, source_local_path)
            for token in entry["arguments"]]
        self.assertNotEqual(entry["arguments"], original)
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        record = retained(value)
        benchmark_pair = next(
            pair for pair in value["identities_initial"]["candidate_build"][
                "validated_compile_commands"]["required_source_object_pairs"]
            if pair["source"]["path"].endswith(
                "/bench/leopard2/benchmark.cpp"))
        record["artifact"]["mtime_ns"] = \
            benchmark_pair["object"]["mtime_ns"] + 1
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        record = retained(value)
        alternate_tree = "a" * 40
        header = runner.benchmark_attestation_header_bytes(
            CANDIDATE_COMMIT, alternate_tree, False).decode("ascii")
        record["source_tree"] = alternate_tree
        record["content"] = runner.exact_text_content(
            header, "coherently relabeled generated attestation")
        record["artifact"]["size"] = record["content"]["size"]
        record["artifact"]["sha256"] = record["content"]["sha256"]
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "validated_compile_commands"].pop("generated_attestation_header")
        synchronize_identity(value)
        self.assert_rejected(value)

    def test_effective_build_configuration_is_a_canonical_compile_input(
        self,
    ) -> None:
        def candidate_build(value: dict) -> dict:
            return value["identities_initial"]["candidate_build"]

        def retained(value: dict) -> dict:
            return candidate_build(value)["validated_compile_commands"][
                "effective_build_configuration"]

        def benchmark_entry(value: dict) -> dict:
            return next(
                entry for entry in candidate_build(value)[
                    "validated_compile_commands"]["required_entries"]
                if entry["file"].endswith("/bench/leopard2/benchmark.cpp"))

        value = synthetic_raw()
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

        for label, mutate in (
            ("configuration path", lambda value: retained(value)["artifact"].update({
                "path": SPECIFICATION["candidate_build_dir"] +
                        "/foreign/build_configuration.txt"})),
            ("helper path", lambda value: retained(value)["helper_source"].update({
                "path": SPECIFICATION["candidate_source_root"] +
                        "/cmake/foreign.cmake"})),
            ("embedded type", lambda value: retained(value).update({
                "embedded_build_type": "Debug"})),
            ("missing record", lambda value: candidate_build(value)[
                "validated_compile_commands"].pop(
                    "effective_build_configuration")),
            ("evidence helper path", lambda value: value[
                "identities_initial"]["evidence_helper"].update({
                    "path": SPECIFICATION["candidate_source_root"] +
                            "/experiments/leopard2/main_compare/run_abba.py"})),
            ("missing evidence helper", lambda value: value[
                "identities_initial"].pop("evidence_helper")),
        ):
            with self.subTest(label=label):
                edited = synthetic_raw()
                mutate(edited)
                synchronize_identity(edited)
                self.assert_rejected(edited)

        value = synthetic_raw()
        entry = benchmark_entry(value)
        entry["arguments"] = [
            ('-DLEO2_BENCHMARK_BUILD_TYPE="Debug"'
             if token == '-DLEO2_BENCHMARK_BUILD_TYPE="Release"' else token)
            for token in entry["arguments"]]
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        record = retained(value)
        benchmark_pair = next(
            pair for pair in candidate_build(value)[
                "validated_compile_commands"]["required_source_object_pairs"]
            if pair["source"]["path"].endswith("/bench/leopard2/benchmark.cpp"))
        record["artifact"]["mtime_ns"] = \
            benchmark_pair["object"]["mtime_ns"] + 1
        synchronize_identity(value)
        self.assert_rejected(value)

        # Even a coherent digest/cache/compiler-argv rewrite cannot turn a
        # noncanonical selector into acceptable Release/AUTO evidence.
        for variable, invalid_value in (
            ("LEO2_EXPERIMENT_DIRECT_SOURCE_PLAN", "ON"),
            ("LEO2_EXPERIMENT_HIGH_DIRECT_ENCODE", "ON"),
            ("LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT", "ON"),
            ("LEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE", "1"),
        ):
            with self.subTest(selector=variable):
                value = synthetic_raw()
                record = retained(value)
                record["entries"][variable] = invalid_value
                material = runner.build_configuration_material(
                    record["entries"])
                digest = runner.sha256_bytes(material)
                text = (
                    f"schema={runner.BUILD_CONFIGURATION_FILE_SCHEMA}\n"
                    f"sha256={digest}\n").encode("ascii") + material
                record["content"] = runner.exact_text_content(
                    text.decode("utf-8"),
                    "coherently rewritten configuration")
                record["artifact"].update({
                    "size": record["content"]["size"],
                    "sha256": record["content"]["sha256"],
                })
                record["configuration_sha256"] = digest
                build = candidate_build(value)
                build["validated_cache"][variable] = invalid_value
                build["validated_cache"][
                    "LEO2_BENCHMARK_EFFECTIVE_CONFIGURATION_SHA256"] = digest
                entry = benchmark_entry(value)
                entry["arguments"] = [
                    (f'-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256="{digest}"'
                     if token.startswith(
                         "-DLEO2_BENCHMARK_BUILD_CONFIGURATION_SHA256=")
                     else token)
                    for token in entry["arguments"]]
                synchronize_identity(value)
                self.assert_rejected(value)

    def test_effective_build_type_single_and_multi_config_semantics(self) -> None:
        entries = complete_build_fixture("candidate")[
            "validated_compile_commands"]["effective_build_configuration"][
                "entries"]
        self.assertEqual(
            runner.validate_embedded_build_type(
                entries, "Release", authoritative=True),
            "Release")

        single_mismatch = copy.deepcopy(entries)
        single_mismatch["CMAKE_BUILD_TYPE"] = "Debug"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_embedded_build_type(
                single_mismatch, "Release", authoritative=True)

        multi = copy.deepcopy(entries)
        multi.update({
            "CMAKE_BUILD_TYPE": "",
            "CMAKE_GENERATOR": "Ninja Multi-Config",
            "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
        })
        self.assertEqual(
            runner.validate_embedded_build_type(
                multi, "Release", authoritative=True),
            "Release")
        self.assertEqual(
            runner.validate_canonical_build_configuration_entries(
                multi, {
                    "CMAKE_BUILD_TYPE": "",
                    "CMAKE_GENERATOR": "Ninja Multi-Config",
                    "CMAKE_CONFIGURATION_TYPES": "Debug;Release",
                    "CMAKE_CXX_COMPILER": multi["CMAKE_CXX_COMPILER"],
                }),
            multi)
        for encoded_types, embedded in (
            ("Debug;RelWithDebInfo", "Release"),
            ("Debug;Release;Release", "Release"),
            ("Debug;Release", "Debug"),
        ):
            with self.subTest(types=encoded_types, embedded=embedded):
                edited = copy.deepcopy(multi)
                edited["CMAKE_CONFIGURATION_TYPES"] = encoded_types
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_embedded_build_type(
                        edited, embedded, authoritative=True)

    def test_cmake_cache_selected_keys_have_unique_canonical_types(self) -> None:
        values = {
            name: (
                "Release" if name == "CMAKE_BUILD_TYPE" else
                "Unix Makefiles" if name == "CMAKE_GENERATOR" else
                "/tmp/value" if next(iter(types)) in {"FILEPATH", "PATH"} else
                "ON" if next(iter(types)) == "BOOL" else "fixture")
            for name, types in runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES.items()
        }

        def render(overrides: Mapping[str, str] | None = None) -> str:
            selected = dict(overrides or {})
            return "".join(
                f"{name}:{selected.get(name, next(iter(types)))}="
                f"{values[name]}\n"
                for name, types in
                runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES.items())

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "CMakeCache.txt"
            path.write_text(render(), encoding="utf-8")
            parsed = runner.parse_cmake_cache(path)
            self.assertEqual(
                set(parsed), set(runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES))

            for name, types in runner.CMAKE_CACHE_REQUIRED_ENTRY_TYPES.items():
                canonical_type = next(iter(types))
                wrong_type = next(
                    item for item in sorted(runner.CMAKE_CACHE_ENTRY_TYPES)
                    if item != canonical_type)
                with self.subTest(key=name, wrong_type=wrong_type):
                    path.write_text(
                        render({name: wrong_type}), encoding="utf-8")
                    with self.assertRaises(runner.EvidenceError):
                        runner.parse_cmake_cache(path)

            duplicate = render() + "CMAKE_GENERATOR:INTERNAL=Ninja\n"
            path.write_text(duplicate, encoding="utf-8")
            with self.assertRaises(runner.EvidenceError):
                runner.parse_cmake_cache(path)

            path.write_text(
                render().replace(
                    "CMAKE_GENERATOR:INTERNAL=",
                    "CMAKE_GENERATOR=", 1),
                encoding="utf-8")
            with self.assertRaises(runner.EvidenceError):
                runner.parse_cmake_cache(path)

            path.write_text(
                render().replace(
                    "CMAKE_GENERATOR:INTERNAL=",
                    "CMAKE_GENERATOR:FOREIGN=", 1),
                encoding="utf-8")
            with self.assertRaises(runner.EvidenceError):
                runner.parse_cmake_cache(path)

    def test_exact_lf_and_strict_utf8_text_boundaries(self) -> None:
        entries = complete_build_fixture("candidate")[
            "validated_compile_commands"]["effective_build_configuration"][
                "entries"]
        entries = copy.deepcopy(entries)
        entries["CMAKE_CXX_COMPILER"] = "/usr/bin/c++\u2028fixture"
        material = runner.build_configuration_material(entries)
        digest = runner.sha256_bytes(material)
        retained = (
            f"schema={runner.BUILD_CONFIGURATION_FILE_SCHEMA}\n"
            f"sha256={digest}\n").encode("ascii") + material
        self.assertEqual(
            runner.parse_build_configuration_bytes(retained)["entries"],
            entries)
        plan = load_plan_runner()
        self.assertEqual(
            plan._parse_build_configuration_bytes(retained)["entries"],
            entries)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "CMakeCache.txt"
            path.write_text(
                "CMAKE_CXX_FLAGS:STRING=-O3\u2028-fno-omit-frame-pointer\n",
                encoding="utf-8")
            self.assertEqual(
                runner.parse_cmake_cache(path)["CMAKE_CXX_FLAGS"],
                "-O3\u2028-fno-omit-frame-pointer")

        with self.assertRaises(runner.EvidenceError):
            runner.exact_text_content("\ud800", "lone surrogate")
        with self.assertRaises(runner.EvidenceError):
            runner.canonical_bytes({"text": "\ud800"})
        with self.assertRaises(runner.EvidenceError):
            runner.build_configuration_material({
                **entries, "CMAKE_CXX_FLAGS": "\ud800"})
        with self.assertRaises(plan.PlanError):
            plan.canonical_bytes({"text": "\ud800"})
        with self.assertRaises(plan.PlanError):
            plan._validate_scope_text({
                "encoding": "utf-8", "text": "\ud800",
                "size": 1, "sha256": "0" * 64,
            }, "lone surrogate")

    def test_single_and_multi_config_complete_closures_are_accepted(self) -> None:
        runner.validate_raw(
            synthetic_raw(), None,
            check_files=False, check_current_inputs=False)
        multi = synthetic_raw(multi_config=True)
        runner.validate_raw(
            multi, None, check_files=False, check_current_inputs=False)
        for role in ("baseline", "candidate"):
            build = multi["identities_initial"][f"{role}_build"]
            self.assertEqual(
                build["validated_compile_commands"]["entry_count"],
                (runner.BASELINE_EXPECTED_COMPILE_COMMAND_COUNT
                 if role == "baseline" else
                 runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT) * 2)
            self.assertIn(
                "/Release/",
                build["validated_compile_commands"]["required_entries"][0][
                    "output"])
            self.assertTrue(
                build["validated_executable"]["path"].endswith(
                    "/Release/" + (
                        "leopard_main_benchmark" if role == "baseline"
                        else "bench_leopard2")))
            self.assertEqual(
                build["multi_config_build_tool"]["path"],
                runner.CANONICAL_NINJA_PATH)
            self.assertEqual(
                build["multi_config_build_tool_version_stdout"]["text"],
                "1.11.1\n")
            self.assertEqual(
                [record["relative_path"]
                 for record in build["multi_config_ninja_graph"]["files"]],
                [
                    "CMakeFiles/common.ninja",
                    "CMakeFiles/impl-Release.ninja",
                    "build-Release.ninja",
                ])

        plan = load_plan_runner()
        identities = multi["identities_initial"]
        replacements = (
            (SPECIFICATION["baseline_source_root"], "$BASELINE_SOURCE"),
            (SPECIFICATION["candidate_source_root"], "$CANDIDATE_SOURCE"),
            (SPECIFICATION["baseline_build_dir"], "$BASELINE_BUILD"),
            (SPECIFICATION["candidate_build_dir"], "$CANDIDATE_BUILD"),
        )
        for role in ("baseline", "candidate"):
            source = plan._normalize_bound_paths(
                identities[f"{role}_source"], replacements)
            build = plan._normalize_bound_paths(
                identities[f"{role}_build"], replacements)
            plan._validate_scope_build(build, role, source)

        def mutate_top_level_include(value: dict) -> None:
            replace_ninja_graph_fixture_text(
                value["identities_initial"]["candidate_build"],
                "build-Release.ninja",
                "include CMakeFiles/substituted.ninja\n")

        def mutate_transitive_include(value: dict) -> None:
            replace_ninja_graph_fixture_text(
                value["identities_initial"]["candidate_build"],
                "CMakeFiles/impl-Release.ninja",
                "include CMakeFiles/common.ninja\n"
                "include CMakeFiles/late.ninja\n"
                "build Release/libleopard.a: phony\n"
                "build Release/bench_leopard2: phony "
                "Release/libleopard.a\n")

        for label, mutate in (
            ("configuration inventory", lambda value:
                value["identities_initial"]["baseline_build"][
                    "validated_cache"].update({
                        "CMAKE_CONFIGURATION_TYPES":
                            "Debug;Release;RelWithDebInfo"})),
            ("missing make program", lambda value:
                value["identities_initial"]["candidate_build"][
                    "validated_cache"].pop("CMAKE_MAKE_PROGRAM")),
            ("substituted Ninja", lambda value:
                value["identities_initial"]["candidate_build"][
                    "multi_config_build_tool"].update({
                        "path": "/tmp/ninja"})),
            ("Ninja version", lambda value:
                value["identities_initial"]["candidate_build"][
                    "multi_config_build_tool_version_stdout"].update({
                        "text": "9.9.9\n",
                        "sha256": runner.sha256_bytes(b"9.9.9\n"),
                    })),
            ("top-level Ninja include substitution", mutate_top_level_include),
            ("transitive Ninja include mutation", mutate_transitive_include),
            ("Ninja graph artifact substitution", lambda value:
                value["identities_initial"]["candidate_build"][
                    "multi_config_ninja_graph"]["files"][0]["artifact"].update({
                        "path":
                            f"{SPECIFICATION['candidate_build_dir']}/"
                            "CMakeFiles/substituted.ninja"})),
            ("entry count", lambda value:
                value["identities_initial"]["candidate_build"][
                    "validated_compile_commands"].update({
                        "entry_count":
                            runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT})),
            ("selected compiler output", lambda value:
                value["identities_initial"]["candidate_build"][
                    "validated_compile_commands"]["required_entries"][0].update({
                        "output": value["identities_initial"]["candidate_build"][
                            "validated_compile_commands"]["required_entries"][0][
                                "output"].replace(
                                    "/Release/", "/Debug/", 1)})),
            ("selected executable", lambda value:
                value["identities_initial"]["candidate_build"][
                    "validated_executable"].update({
                        "path":
                            f"{SPECIFICATION['candidate_build_dir']}/"
                            "Debug/bench_leopard2"})),
        ):
            with self.subTest(mutation=label):
                edited = synthetic_raw(multi_config=True)
                mutate(edited)
                synchronize_identity(edited)
                self.assert_rejected(edited)

    def test_ninja_graph_capture_rejects_substitution_and_extraction_race(
            self) -> None:
        with tempfile.TemporaryDirectory(
                prefix="leopard2-ninja-graph-race-") as directory:
            build = Path(directory).resolve()
            (build / "CMakeFiles").mkdir()
            (build / "build-Release.ninja").write_text(
                "include CMakeFiles/impl-Release.ninja\n",
                encoding="utf-8")
            (build / "CMakeFiles/impl-Release.ninja").write_text(
                "include CMakeFiles/common.ninja\n",
                encoding="utf-8")
            common = build / "CMakeFiles/common.ninja"
            common.write_text(
                "rule fixture_noop\n  command = :\n", encoding="utf-8")
            runner.capture_ninja_graph_closure(build, "Release")

            replacement = build / "replacement.ninja"
            replacement.write_text(
                "rule replacement\n  command = false\n", encoding="utf-8")
            common.unlink()
            common.symlink_to(replacement)
            with self.assertRaises(runner.EvidenceError):
                runner.capture_ninja_graph_closure(build, "Release")
            common.unlink()
            common.write_text(
                "rule fixture_noop\n  command = :\n", encoding="utf-8")

            tool_identity = complete_artifact(
                runner.CANONICAL_NINJA_PATH, "build_tool", "a")
            tool_version = {
                "sha256": runner.sha256_bytes(b"1.11.1\n"),
                "text": "1.11.1\n",
            }

            def mutate_during_extraction(*_args, **_kwargs):
                common.write_text(
                    "rule raced\n  command = false\n", encoding="utf-8")
                return ("archive recipe\n", "executable recipe\n")

            cache = {
                "CMAKE_GENERATOR": "Ninja Multi-Config",
                "CMAKE_MAKE_PROGRAM": runner.CANONICAL_NINJA_PATH,
            }
            with mock.patch.object(
                    runner, "canonical_ninja_identity",
                    return_value=(tool_identity, tool_version)), \
                    mock.patch.object(
                        runner, "ninja_multi_link_recipes",
                        side_effect=mutate_during_extraction), \
                    self.assertRaises(runner.EvidenceError):
                runner.stable_ninja_multi_link_recipes(
                    cache, build, "Release", "Release/bench_leopard2",
                    "Release/libleopard.a", "Release/bench_leopard2")

    def test_normalized_multi_config_ninja_link_identity_tracks_graph(
            self) -> None:
        build = complete_build_fixture("candidate", multi_config=True)
        impl_relative = "CMakeFiles/impl-Release.ninja"
        impl = next(
            record for record in build["multi_config_ninja_graph"]["files"]
            if record["relative_path"] == impl_relative)
        replace_ninja_graph_fixture_text(
            build, impl_relative,
            impl["content"]["text"] +
            "# absolute source " +
            SPECIFICATION["candidate_source_root"] +
            "/leopard2.cpp\n")
        runner.validate_complete_build_identity(
            build, "candidate", SPECIFICATION, runner.RAW_SCHEMA,
            CANDIDATE_TREE)

        plan = load_plan_runner()
        replacements = (
            (SPECIFICATION["baseline_source_root"], "$BASELINE_SOURCE"),
            (SPECIFICATION["candidate_source_root"], "$CANDIDATE_SOURCE"),
            (SPECIFICATION["baseline_build_dir"], "$BASELINE_BUILD"),
            (SPECIFICATION["candidate_build_dir"], "$CANDIDATE_BUILD"),
        )
        normalized_build = plan._normalize_bound_paths(build, replacements)
        normalized_source = plan._normalize_bound_paths(
            complete_source_fixture("candidate"), replacements)
        normalized_impl = next(
            record for record in
                normalized_build["multi_config_ninja_graph"]["files"]
            if record["relative_path"] == impl_relative)
        graph_artifact = normalized_impl["artifact"]
        self.assertIn(
            "$CANDIDATE_SOURCE/leopard2.cpp",
            normalized_impl["content"]["text"])
        self.assertNotEqual(
            graph_artifact["sha256"], impl["artifact"]["sha256"])
        for name in ("archive_link_recipe", "executable_link_recipe"):
            with self.subTest(identity=name):
                identity = normalized_build[name]
                self.assertEqual(
                    {key: identity[key]
                     for key in ("path", "size", "mode", "sha256")},
                    {key: graph_artifact[key]
                     for key in ("path", "size", "mode", "sha256")})
        self.assertNotEqual(
            graph_artifact["sha256"],
            normalized_build["archive_link_recipe_content"]["sha256"])
        plan._validate_scope_build(
            normalized_build, "candidate", normalized_source)

    def test_real_single_and_ninja_multi_config_build_provenance(self) -> None:
        if shutil.which("cmake") is None or shutil.which("ninja") is None or \
                shutil.which("make") is None:
            self.skipTest("CMake, Ninja, and Make are required")
        source_checkout = MODULE_PATH.resolve().parents[3]
        tracked = subprocess.run(
            ["git", "-C", str(source_checkout), "ls-files", "-z"],
            check=True, stdout=subprocess.PIPE).stdout.split(b"\0")
        with tempfile.TemporaryDirectory(
                prefix="leopard2-real-cmake-provenance-") as directory:
            root = Path(directory)
            source = root / "source"
            source.mkdir()
            for encoded in tracked:
                if not encoded:
                    continue
                relative = os.fsdecode(encoded)
                original = source_checkout / relative
                destination = source / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                if original.is_dir():
                    # A gitlink such as sse2neon is not a regular tracked file;
                    # the x86 provenance fixtures do not consume its checkout.
                    continue
                if original.is_symlink():
                    destination.symlink_to(os.readlink(original))
                else:
                    shutil.copy2(original, destination)
            subprocess.run(["git", "init", "-q", str(source)], check=True)
            subprocess.run(
                ["git", "-C", str(source), "config", "user.email",
                 "leopard2-fixture@example.invalid"], check=True)
            subprocess.run(
                ["git", "-C", str(source), "config", "user.name",
                 "Leopard2 Fixture"], check=True)
            subprocess.run(
                ["git", "-C", str(source), "add", "-A"], check=True)
            subprocess.run(
                ["git", "-C", str(source), "commit", "-q", "-m",
                 "fixture"], check=True)
            commit = subprocess.run(
                ["git", "-C", str(source), "rev-parse", "HEAD"],
                check=True, stdout=subprocess.PIPE, text=True).stdout.strip()
            tree = subprocess.run(
                ["git", "-C", str(source), "rev-parse", "HEAD^{tree}"],
                check=True, stdout=subprocess.PIPE, text=True).stdout.strip()

            common = [
                "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
                "-DENABLE_OPENMP=ON",
                "-DLEOPARD_ENABLE_GF8=ON",
                "-DLEOPARD_ENABLE_GF16=ON",
                "-DLEO2_BACKEND_VARIANT=auto",
                "-DLEO2_BUILD_BENCHMARKS=ON",
                "-DLEO2_BUILD_FUZZERS=OFF",
                "-DLEO2_BUILD_TESTS=OFF",
                "-DLEO2_ENABLE_CUDA=OFF",
                "-DLEO2_EXPERIMENT_DIRECT_SOURCE_PLAN=OFF",
                "-DLEO2_EXPERIMENT_HIGH_DIRECT_ENCODE=OFF",
                "-DLEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT=OFF",
                "-DLEO2_EXPERIMENT_GF8_SMALL_DIRECT_MODE=0",
            ]
            for label, generator, configure_extra, output_prefix in (
                ("single", "Unix Makefiles",
                 ["-DCMAKE_BUILD_TYPE=Release"], Path()),
                ("multi", "Ninja Multi-Config",
                 ["-DCMAKE_CONFIGURATION_TYPES=Debug;Release"],
                 Path("Release")),
            ):
                with self.subTest(layout=label):
                    build = root / f"build-{label}"
                    subprocess.run(
                        ["cmake", "-S", str(source), "-B", str(build),
                         "-G", generator, *configure_extra, *common],
                        check=True, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, timeout=120)
                    build_command = [
                        "cmake", "--build", str(build), "--parallel", "16",
                        "--target", "bench_leopard2",
                    ]
                    if label == "multi":
                        build_command.extend(["--config", "Release"])
                    subprocess.run(
                        build_command, check=True, stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT, timeout=120)
                    executable = build / output_prefix / "bench_leopard2"
                    archive = build / output_prefix / "libleopard.a"
                    specification = {
                        "baseline_source_root": str(source),
                        "candidate_source_root": str(source),
                        "candidate_build_dir": str(build),
                        "candidate_executable": str(executable),
                        "candidate_archive": str(archive),
                        "candidate_commit": commit,
                    }
                    provenance = runner.build_provenance(
                        "candidate", specification, runner.RAW_SCHEMA)
                    runner.validate_complete_build_identity(
                        provenance, "candidate", specification,
                        runner.RAW_SCHEMA, tree)
                    self.assertEqual(
                        provenance["validated_compile_commands"][
                            "entry_count"],
                        runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT *
                        (2 if label == "multi" else 1))
                    if label == "multi":
                        self.assertEqual(
                            provenance["multi_config_build_tool"]["path"],
                            runner.CANONICAL_NINJA_PATH)
                        self.assertTrue(
                            provenance[
                                "multi_config_build_tool_version_stdout"][
                                    "text"])
                        self.assertGreaterEqual(
                            len(provenance[
                                "multi_config_ninja_graph"]["files"]), 2)
                    else:
                        self.assertIsNone(
                            provenance["multi_config_build_tool"])
                        self.assertIsNone(
                            provenance[
                                "multi_config_build_tool_version_stdout"])
                        self.assertIsNone(
                            provenance["multi_config_ninja_graph"])

    def test_recipe_command_semantic_mutations_are_rejected(self) -> None:
        canonical = archive_recipe_fixture_text(
            runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA))
        mutations = {
            "noncanonical output path": canonical.replace(
                "ar qc libleopard.a", "ar qc nested/libleopard.a", 1),
            "different archiver": canonical.replace(
                "/usr/bin/ar", "/tmp/ar", 1),
            "response file": canonical.replace(
                "LeopardCommon.cpp.o ", "LeopardCommon.cpp.o @objects.rsp ", 1),
            "object order": canonical.replace(
                "CMakeFiles/leopard.dir/LeopardCommon.cpp.o "
                "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o",
                "CMakeFiles/leopard2_backend_avx2.dir/Leopard2BackendAVX2.cpp.o "
                "CMakeFiles/leopard.dir/LeopardCommon.cpp.o"),
            "different ranlib": canonical.replace(
                "/usr/bin/ranlib", "/tmp/ranlib", 1),
        }
        for label, recipe in mutations.items():
            with self.subTest(label=label):
                value = synthetic_raw()
                replace_current_recipe_text(value, recipe)
                self.assert_rejected(value)

        value = synthetic_raw()
        value["identities_initial"]["candidate_build"][
            "archive_link_recipe"]["sha256"] = "b" * 64
        value["identities_final"] = copy.deepcopy(value["identities_initial"])
        for invocation in value["invocations"]:
            invocation["identity_before"] = copy.deepcopy(
                value["identities_initial"])
            invocation["identity_after"] = copy.deepcopy(
                value["identities_initial"])
        self.assert_rejected(value)

    def test_executable_recipe_requires_the_declared_archive_operand(self) -> None:
        canonical = complete_build_fixture("candidate")[
            "executable_link_recipe_content"]["text"]
        mutations = {
            "second flag command": canonical + "\n-O3\n",
            "leading flag command": "-O3\r\n" + canonical,
            "absolute same-basename archive": canonical.replace(
                " libleopard.a ", " /tmp/libleopard.a ", 1),
            "sibling same-basename archive": canonical.replace(
                " libleopard.a ", " ../foreign/libleopard.a ", 1),
            "duplicate same-basename archive": canonical.replace(
                " libleopard.a ",
                " libleopard.a /tmp/libleopard.a ", 1),
            "foreign differently-named archive": canonical.replace(
                " libleopard.a ",
                " libleopard.a /tmp/evil.a ", 1),
            "opaque library switch": canonical.replace(
                " libleopard.a ", " libleopard.a -levil ", 1),
            "linker response file": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,@evil.rsp ", 1),
            "driver response file": canonical.replace(
                " libleopard.a ", " libleopard.a @evil.rsp ", 1),
            "fused linker script": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,--script=/tmp/evil.ld ", 1),
            "comma linker script": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,-T,/tmp/evil.ld ", 1),
            "fused library search": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,-L,/tmp,-levil ", 1),
            "compiler specs": canonical.replace(
                " libleopard.a ",
                " libleopard.a -specs=/tmp/evil.specs ", 1),
            "alternate tool root": canonical.replace(
                " libleopard.a ", " libleopard.a -B/tmp/toolchain ", 1),
            "compiler plugin": canonical.replace(
                " libleopard.a ",
                " libleopard.a -fplugin=/tmp/evil.so ", 1),
            "linker plugin": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,--plugin,/tmp/evil.so ", 1),
            "alternate linker": canonical.replace(
                " libleopard.a ",
                " libleopard.a -fuse-ld=/tmp/evil-ld ", 1),
            "alternate runtime loader": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Wl,--dynamic-linker,/tmp/evil-ld.so ", 1),
            "xlinker control": canonical.replace(
                " libleopard.a ",
                " libleopard.a -Xlinker /tmp/evil.ld ", 1),
            "undeclared shared input": canonical.replace(
                " libleopard.a ",
                " libleopard.a /tmp/libevil.so ", 1),
            "duplicate OpenMP runtime": canonical.replace(
                " /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so ",
                " /usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so "
                "/usr/lib/gcc/x86_64-linux-gnu/13/libgomp.so ", 1),
            "duplicate pthread archive": canonical.replace(
                " /usr/lib/x86_64-linux-gnu/libpthread.a",
                " /usr/lib/x86_64-linux-gnu/libpthread.a "
                "/usr/lib/x86_64-linux-gnu/libpthread.a", 1),
        }
        for label, recipe in mutations.items():
            with self.subTest(label=label):
                value = synthetic_raw()
                replace_current_executable_recipe_text(value, recipe)
                self.assert_rejected(value)

        value = synthetic_raw()
        runner.validate_raw(
            resign(value), None, check_files=False, check_current_inputs=False)

    def test_executable_recipe_parser_preserves_command_boundaries(self) -> None:
        canonical = complete_build_fixture("candidate")[
            "executable_link_recipe_content"]["text"]
        self.assertEqual(
            runner.parse_single_executable_recipe(
                canonical, "canonical executable recipe"),
            runner.shlex.split(canonical, posix=True))
        for label, recipe in (
            ("second LF command", canonical + "\n-O3\n"),
            ("leading CRLF command", "-O3\r\n" + canonical),
            ("second Unicode line command", canonical + "\u2028-O3"),
        ):
            with self.subTest(label=label), \
                 self.assertRaises(runner.EvidenceError):
                runner.parse_single_executable_recipe(
                    recipe, "multiline executable recipe")

    def test_external_link_inputs_bind_paths_roles_and_bytes(self) -> None:
        canonical = complete_build_fixture("candidate")[
            "executable_link_recipe_content"]["text"]

        for label, replacement in (
            ("temporary pthread root", "/tmp/libpthread.a"),
            ("alternate system pthread root",
             "/usr/lib/aarch64-linux-gnu/libpthread.a"),
        ):
            with self.subTest(label=label):
                value = synthetic_raw()
                build = value["identities_initial"]["candidate_build"]
                pthread = build["validated_external_link_inputs"][1]
                pthread["operand"] = replacement
                pthread["artifact"]["path"] = replacement
                replace_current_executable_recipe_text(
                    value, canonical.replace(
                        "/usr/lib/x86_64-linux-gnu/libpthread.a",
                        replacement, 1))
                self.assert_rejected(value)

        value = synthetic_raw()
        value["invocations"][0]["identity_before"] = copy.deepcopy(
            value["identities_initial"])
        value["invocations"][0]["identity_before"]["candidate_build"][
            "validated_external_link_inputs"][1]["artifact"]["sha256"] = "e" * 64
        # The recipe still names the same lexical archive, but one retained
        # identity copy no longer agrees with the signed campaign snapshots.
        self.assert_rejected(value)

        for mutation in ("missing", "duplicate", "swapped"):
            with self.subTest(mutation=mutation):
                value = synthetic_raw()
                records = value["identities_initial"]["candidate_build"][
                    "validated_external_link_inputs"]
                if mutation == "missing":
                    records.pop()
                elif mutation == "duplicate":
                    records.append(copy.deepcopy(records[-1]))
                else:
                    records.reverse()
                synchronize_identity(value)
                self.assert_rejected(value)

        for label, mutate in (
            ("foreign resolved pthread path", lambda record:
                record["artifact"].update({"path": "/tmp/libpthread.a"})),
            ("foreign resolved OpenMP root", lambda record:
                record["artifact"].update({
                    "path": "/usr/lib/aarch64-linux-gnu/libgomp.so.1.0.0"})),
            ("wrong external role", lambda record:
                record.update({"role": "openmp_runtime_shared"})),
            ("truncated pthread archive identity", lambda record:
                record["artifact"].update({"size": 1})),
            ("malformed pthread archive digest", lambda record:
                record["artifact"].update({"sha256": "0" * 63})),
        ):
            with self.subTest(label=label):
                value = synthetic_raw()
                records = value["identities_initial"]["candidate_build"][
                    "validated_external_link_inputs"]
                mutate(records[0] if "OpenMP" in label else records[1])
                synchronize_identity(value)
                self.assert_rejected(value)

        value = synthetic_raw()
        candidate = value["identities_initial"]["candidate_build"]
        candidate["validated_external_link_inputs"][1]["artifact"][
            "sha256"] = "e" * 64
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        for role in ("baseline_build", "candidate_build"):
            value["identities_initial"][role][
                "validated_external_link_inputs"][0]["artifact"]["path"] = \
                "/usr/lib/x86_64-linux-gnu/libgomp.so.999.999"
        synchronize_identity(value)
        self.assert_rejected(value)

        value = synthetic_raw()
        candidate = value["identities_initial"]["candidate_build"]
        candidate["validated_external_link_inputs"][1]["artifact"][
            "mtime_ns"] = candidate["validated_executable"]["mtime_ns"] + 1
        synchronize_identity(value)
        self.assert_rejected(value)

        for label, mutate in (
            ("external file mode", lambda record:
                record["artifact"].update({
                    "mode": record["artifact"]["mode"] ^ 0o100})),
            ("lexical link target", lambda record:
                record["lexical_symlink_chain"][0].update({
                    "target": "../../../x86_64-linux-gnu/libgomp.so.999"})),
            ("lexical link mode", lambda record:
                record["lexical_symlink_chain"][0].update({"mode": 0o700})),
            ("truncated lexical link chain", lambda record:
                record["lexical_symlink_chain"].pop()),
        ):
            with self.subTest(external_chain=label):
                value = synthetic_raw()
                record = value["identities_initial"]["candidate_build"][
                    "validated_external_link_inputs"][0]
                mutate(record)
                synchronize_identity(value)
                self.assert_rejected(value)

        value = synthetic_raw()
        pthread = value["identities_initial"]["candidate_build"][
            "validated_external_link_inputs"][1]
        pthread["lexical_symlink_chain"] = [{
            "path": pthread["operand"], "target": "libpthread-evil.a",
            "mode": 0o777,
        }]
        synchronize_identity(value)
        self.assert_rejected(value)

        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "external.bin"
            operand = Path(directory) / "external-link"
            old_bytes = b"0123456789abcdef"
            new_bytes = b"fedcba9876543210"
            target.write_bytes(old_bytes)
            operand.symlink_to(target.name)
            metadata, digest, snapshot, resolved, chain = \
                runner.link_common.current_external_file_snapshot(
                    operand, "immutable external fixture")
            self.assertEqual(snapshot, old_bytes)
            self.assertEqual(resolved, target)
            self.assertEqual(chain[0]["path"], str(operand))
            target.write_bytes(new_bytes)
            os.utime(target, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))
            _, new_digest, new_snapshot, _, _ = \
                runner.link_common.current_external_file_snapshot(
                    operand, "rewritten external fixture")
            self.assertEqual(snapshot, old_bytes)
            self.assertEqual(digest, hashlib.sha256(old_bytes).hexdigest())
            self.assertEqual(new_snapshot, new_bytes)
            self.assertNotEqual(new_digest, digest)

            resolver = runner.link_common._resolve_external_lexical_path
            calls = 0

            def retarget_after_read(path: Path, label: str):
                nonlocal calls
                calls += 1
                resolved_path, public, private = resolver(path, label)
                if calls == 2:
                    public = copy.deepcopy(public)
                    public[0]["target"] = "retargeted-after-read"
                return resolved_path, public, private

            with mock.patch.object(
                    runner.link_common, "_resolve_external_lexical_path",
                    side_effect=retarget_after_read), \
                 self.assertRaises(runner.link_common.EvidenceError):
                runner.link_common.current_external_file_snapshot(
                    operand, "retarget race fixture")

    def test_effective_release_and_executable_flags_fail_closed(self) -> None:
        for label, flags in (
            ("final optimization downgrade", "-O3 -O0"),
            ("bare optimization downgrade", "-O3 -O"),
            ("unknown optimization spelling", "-O4 -O3"),
            ("sanitizer after optimization", "-O3 -fsanitize=address"),
            ("profile after optimization", "-O3 -fprofile-generate"),
            ("LTO after optimization", "-O3 -flto"),
            ("instrumentation after optimization",
             "-O3 -finstrument-functions"),
            ("vector disable after optimization",
             "-O3 -fno-tree-vectorize"),
            ("coverage after optimization", "-O3 --coverage"),
            ("long optimize alias", "-O3 --optimize=0"),
            ("long sanitizer alias", "-O3 --sanitize=address"),
            ("long instrumentation alias", "-O3 --instrument-functions"),
            ("long LTO alias", "-O3 --lto"),
            ("long profile alias", "-O3 --profile"),
            ("inline disable", "-O3 -fno-inline"),
            ("GCC pass disable", "-O3 -fdisable-tree-vect"),
            ("Release response file", "-O3 @evil.rsp"),
        ):
            with self.subTest(label=label):
                value = synthetic_raw()
                value["identities_initial"]["candidate_build"][
                    "validated_cache"]["CMAKE_CXX_FLAGS_RELEASE"] = flags
                synchronize_identity(value)
                self.assert_rejected(value)

        canonical = complete_build_fixture("candidate")[
            "executable_link_recipe_content"]["text"]
        for label, recipe in (
            ("recipe -O0", canonical.replace(" -O3 ", " -O0 ", 1)),
            ("recipe trailing -O0", canonical.replace(
                " -O3 ", " -O3 -O0 ", 1)),
            ("recipe sanitizer", canonical.replace(
                " -O3 ", " -O3 -fsanitize=address ", 1)),
            ("recipe long optimize alias", canonical.replace(
                " -O3 ", " -O3 --optimize=0 ", 1)),
            ("recipe inline disable", canonical.replace(
                " -O3 ", " -O3 -fno-inline ", 1)),
        ):
            with self.subTest(label=label):
                value = synthetic_raw()
                replace_current_executable_recipe_text(value, recipe)
                self.assert_rejected(value)

        # The producer-side compile database is independently hostile input.
        # Exercise its complete argv grammar directly; these cases run under
        # both ordinary Python and ``python -O`` with the rest of this suite.
        for role in ("baseline", "candidate"):
            with self.subTest(compile_profile=role), \
                 tempfile.TemporaryDirectory() as directory:
                path, specification, compiler = compile_commands_fixture(
                    Path(directory), role)
                proof = runner.validate_compile_commands(
                    path, role, specification, compiler,
                    compiler_invocation=str(compiler))
                self.assertEqual(proof["schema"], runner.COMPILE_COMMANDS_SCHEMA)
                self.assertEqual(
                    proof["implementation"], role)
                self.assertEqual(
                    proof["profile"],
                    runner.compile_profile_for_implementation(role))
                self.assertEqual(
                    len(proof["required_entries"]),
                    len(proof["required_source_object_pairs"]))
                path.write_bytes(
                    b'[{"file":"first.cpp","file":"second.cpp"}]')
                with self.assertRaisesRegex(
                        runner.EvidenceError, "duplicate key"):
                    runner.validate_compile_commands(
                        path, role, specification, compiler,
                        compiler_invocation=str(compiler))
                path.write_bytes(
                    b'[{"untrusted":' + b"9" * 5000 + b"}]")
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_compile_commands(
                        path, role, specification, compiler,
                        compiler_invocation=str(compiler))

        compile_mutations = (
            ("response file", "candidate", "/leopard.cpp", "response"),
            ("second positional source", "candidate", "/leopard.cpp", "extra_source"),
            ("stdin positional source", "candidate", "/leopard.cpp", "stdin_source"),
            ("different positional source", "candidate", "/leopard.cpp", "wrong_source"),
            ("missing positional source", "candidate", "/leopard.cpp", "missing_source"),
            ("duplicate positional source", "candidate", "/leopard.cpp", "duplicate_source"),
            ("wrong output operand", "candidate", "/leopard.cpp", "wrong_output"),
            ("coherent wrong output metadata", "candidate", "/leopard.cpp", "coherent_output"),
            ("duplicate output option", "candidate", "/leopard.cpp", "duplicate_output"),
            ("missing compile option", "candidate", "/leopard.cpp", "missing_compile"),
            ("duplicate compile option", "candidate", "/leopard.cpp", "duplicate_compile"),
            ("extra definition", "candidate", "/leopard.cpp", "extra_definition"),
            ("reordered definitions", "candidate", "/leopard.cpp", "reorder_definitions"),
            ("extra include", "candidate", "/leopard.cpp", "extra_include"),
            ("different include", "candidate", "/leopard.cpp", "different_include"),
            ("last language option wins", "candidate", "/leopard.cpp", "language_last"),
            ("last optimization option wins", "candidate", "/leopard.cpp", "optimization_last"),
            ("last AVX2 option wins", "candidate", "/Leopard2BackendAVX2.cpp", "avx2_last"),
            ("missing AVX2 isolation", "candidate", "/Leopard2BackendAVX2.cpp", "missing_avx2"),
            ("last native ISA option wins", "baseline", "/leopard.cpp", "baseline_isa_last"),
            ("compiler wrapper", "candidate", "/leopard.cpp", "compiler_wrapper"),
            ("ambiguous argv representations", "candidate", "/leopard.cpp", "both_forms"),
            ("unrelated target response file", "candidate",
             "/tests/benchmark.cpp", "response"),
        )
        for label, role, suffix, mutation in compile_mutations:
            with self.subTest(compiler_argv=label), \
                 tempfile.TemporaryDirectory() as directory:
                path, specification, compiler = compile_commands_fixture(
                    Path(directory), role)
                entries = json.loads(path.read_text(encoding="utf-8"))
                entry = next(item for item in entries
                             if item["file"].endswith(suffix))
                tokens = runner.compile_command_tokens(entry)
                entry.pop("command", None)
                entry["arguments"] = tokens
                output_index = tokens.index("-o")
                compile_index = tokens.index("-c")
                source_index = tokens.index(entry["file"])
                other = next(item for item in entries
                             if item["file"] != entry["file"])
                if mutation == "response":
                    tokens.insert(source_index, "@evil.rsp")
                elif mutation == "extra_source":
                    tokens.insert(source_index, other["file"])
                elif mutation == "stdin_source":
                    tokens.insert(source_index, "-")
                elif mutation == "wrong_source":
                    tokens[source_index] = other["file"]
                elif mutation == "missing_source":
                    tokens.pop(source_index)
                elif mutation == "duplicate_source":
                    tokens.append(entry["file"])
                elif mutation == "wrong_output":
                    tokens[output_index + 1] = other["output"]
                elif mutation == "coherent_output":
                    entry["output"] = other["output"]
                    tokens[output_index + 1] = other["output"]
                elif mutation == "duplicate_output":
                    tokens[compile_index:compile_index] = [
                        "-o", entry["output"]]
                elif mutation == "missing_compile":
                    tokens.pop(compile_index)
                elif mutation == "duplicate_compile":
                    tokens.insert(compile_index, "-c")
                elif mutation == "extra_definition":
                    tokens.insert(1, "-DEVIL=1")
                elif mutation == "reorder_definitions":
                    definitions = [index for index, token in enumerate(tokens)
                                   if token.startswith("-D")]
                    tokens[definitions[0]], tokens[definitions[1]] = \
                        tokens[definitions[1]], tokens[definitions[0]]
                elif mutation == "extra_include":
                    include = next(index for index, token in enumerate(tokens)
                                   if token.startswith("-I"))
                    tokens.insert(include + 1, "-I/tmp/evil")
                elif mutation == "different_include":
                    include = next(index for index, token in enumerate(tokens)
                                   if token.startswith("-I"))
                    tokens[include] = "-I/tmp/evil"
                elif mutation == "language_last":
                    language = tokens.index("-std=gnu++11")
                    tokens.insert(language + 1, "-std=gnu++20")
                elif mutation == "optimization_last":
                    tokens.insert(output_index, "-O0")
                elif mutation == "avx2_last":
                    tokens.insert(output_index, "-mno-avx2")
                elif mutation == "missing_avx2":
                    tokens.remove("-mavx2")
                elif mutation == "baseline_isa_last":
                    tokens.insert(output_index, "-march=x86-64")
                elif mutation == "compiler_wrapper":
                    tokens.insert(1, "--compiler-wrapper=/tmp/evil")
                elif mutation == "both_forms":
                    entry["command"] = runner.shlex.join(tokens)
                else:  # pragma: no cover - table and implementation are local
                    self.fail(f"unknown compiler mutation: {mutation}")
                path.write_text(json.dumps(entries), encoding="utf-8")
                with self.assertRaises(runner.EvidenceError):
                    runner.validate_compile_commands(
                        path, role, specification, compiler,
                        compiler_invocation=str(compiler))

    def test_unrelated_compile_output_may_be_absent_but_argv_stays_closed(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path, specification, compiler = compile_commands_fixture(
                Path(directory), "candidate")
            canonical = json.loads(path.read_text(encoding="utf-8"))
            extra_index = next(
                index for index, entry in enumerate(canonical)
                if entry["file"].endswith("/tests/benchmark.cpp"))
            extra = canonical[extra_index]
            output = Path(specification["candidate_build_dir"]) / extra["output"]
            output.unlink()

            proof = runner.validate_compile_commands(
                path, "candidate", specification, compiler,
                compiler_invocation=str(compiler))
            self.assertEqual(
                proof["entry_count"],
                runner.CANDIDATE_EXPECTED_COMPILE_COMMAND_COUNT)
            self.assertFalse(output.exists())

            outside = Path(directory) / "outside"
            outside.mkdir()
            external_source = outside / "external.cpp"
            external_source.write_text(
                "// external substitute\n", encoding="utf-8")
            escape = Path(specification["candidate_build_dir"]) / "escape"
            escape.symlink_to(outside, target_is_directory=True)
            for label, mutation in (
                ("response file", "response"),
                ("compiler", "compiler"),
                ("second source", "source"),
                ("missing compile option", "compile"),
                ("duplicate output option", "output"),
                ("output path escape", "escape"),
                ("extra definition", "definition"),
                ("coherent external source", "coherent_source"),
                ("coherent output", "coherent_output"),
            ):
                with self.subTest(unrelated_entry=label):
                    entries = copy.deepcopy(canonical)
                    entry = entries[extra_index]
                    tokens = runner.compile_command_tokens(entry)
                    entry.pop("command", None)
                    entry["arguments"] = tokens
                    if mutation == "response":
                        tokens.append("@evil.rsp")
                    elif mutation == "compiler":
                        tokens[0] = "/tmp/unapproved-c++"
                    elif mutation == "source":
                        tokens.append(canonical[0]["file"])
                    elif mutation == "compile":
                        tokens.remove("-c")
                    elif mutation == "output":
                        tokens[1:1] = ["-o", extra["output"]]
                    elif mutation == "escape":
                        entry["output"] = "escape/evil.o"
                        tokens[tokens.index("-o") + 1] = entry["output"]
                    elif mutation == "definition":
                        tokens.insert(1, "-DEVIL=1")
                    elif mutation == "coherent_source":
                        source_index = tokens.index(entry["file"])
                        entry["file"] = str(external_source)
                        tokens[source_index] = entry["file"]
                    elif mutation == "coherent_output":
                        entry["output"] = (
                            "CMakeFiles/bench_leopard.dir/tests/substitute.cpp.o")
                        tokens[tokens.index("-o") + 1] = entry["output"]
                    else:  # pragma: no cover - local table is exhaustive
                        self.fail(f"unknown unrelated-entry mutation: {mutation}")
                    path.write_text(json.dumps(entries), encoding="utf-8")
                    with self.assertRaises(runner.EvidenceError):
                        runner.validate_compile_commands(
                            path, "candidate", specification, compiler,
                            compiler_invocation=str(compiler))

            path.write_text(json.dumps(canonical), encoding="utf-8")
            required = next(
                entry for entry in canonical
                if entry["file"].endswith("/leopard.cpp"))
            required_output = Path(
                specification["candidate_build_dir"]) / required["output"]
            required_output.unlink()
            with self.assertRaises((OSError, runner.EvidenceError)):
                runner.validate_compile_commands(
                    path, "candidate", specification, compiler,
                    compiler_invocation=str(compiler))

    def test_coherent_failed_historical_recipe_relabel_is_rejected(self) -> None:
        historical = synthetic_failure(runner.RAW_SCHEMA_V2)
        runner.validate_failure(historical, Path("/unused"), check_files=False)
        current = synthetic_failure(runner.RAW_SCHEMA)
        runner.validate_failure(current, Path("/unused"), check_files=False)

        relabeled = recursively_replace_strings(historical, (
            ("liblibleopard.a", "libleopard.a"),
            ("CMakeFiles/libleopard.dir", "CMakeFiles/leopard.dir"),
        ))
        self.assertIsInstance(relabeled, dict)
        relabeled["schema"] = runner.FAILURE_SCHEMA
        old_identity = runner.cmake_identity_for_raw_schema(runner.RAW_SCHEMA_V2)
        attach_recipe_content(relabeled, runner.exact_text_content(
            archive_recipe_fixture_text(old_identity),
            "historical failed fixture archive link recipe"))
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(relabeled), Path("/unused"), check_files=False)

    def test_partial_failure_v8_v9_relabels_are_rejected(self) -> None:
        def before_topology(raw_schema: str) -> dict:
            value = synthetic_failure(raw_schema)
            value.pop("digest")
            value["campaign"].pop("allowed_cpu_set_at_launch")
            for field in (
                    "host_initial", "reservation", "supervision", "pair_lease",
                    "isolation", "identities_initial", "executable_snapshots"):
                value[field] = None
            value["invocations"] = []
            value["retained_files"] = []
            return runner.signed(value)

        historical = before_topology(runner.RAW_SCHEMA_V8)
        current = before_topology(runner.RAW_SCHEMA)
        runner.validate_failure(
            historical, Path("/unused"), check_files=False)
        runner.validate_failure(
            current, Path("/unused"), check_files=False)

        upgraded = copy.deepcopy(historical)
        upgraded["schema"] = runner.FAILURE_SCHEMA
        with self.assertRaisesRegex(
                runner.EvidenceError, "unexpected or missing fields"):
            runner.validate_failure(
                resign(upgraded), Path("/unused"), check_files=False)

        downgraded = copy.deepcopy(current)
        downgraded["schema"] = runner.FAILURE_SCHEMA_V8
        with self.assertRaisesRegex(
                runner.EvidenceError, "unexpected or missing fields"):
            runner.validate_failure(
                resign(downgraded), Path("/unused"), check_files=False)

        missing_contract = copy.deepcopy(current)
        missing_contract.pop("evidence_contract")
        with self.assertRaisesRegex(
                runner.EvidenceError, "unexpected or missing fields"):
            runner.validate_failure(
                resign(missing_contract), Path("/unused"), check_files=False)

        wrong_contract = copy.deepcopy(current)
        wrong_contract["evidence_contract"] += "-forged"
        with self.assertRaisesRegex(
                runner.EvidenceError, "evidence contract differs"):
            runner.validate_failure(
                resign(wrong_contract), Path("/unused"), check_files=False)

    def test_non_string_schema_values_fail_closed(self) -> None:
        value = synthetic_raw()
        value["schema"] = {"unexpected": "object"}
        self.assert_rejected(value)
        failure = synthetic_failure(runner.RAW_SCHEMA)
        failure["schema"] = [runner.FAILURE_SCHEMA]
        with self.assertRaises(runner.EvidenceError):
            runner.validate_failure(
                resign(failure), Path("/unused"), check_files=False)

    def test_slower_candidate_is_valid_evidence(self) -> None:
        value = synthetic_raw(candidate_scale=2.0)
        analysis = runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        result = analysis[CELL.identifier]["encode"]
        self.assertLess(result["geometric_speedup"], 1.0)
        self.assertFalse(result["faster_ci_excludes_one"])

    def test_parameter_and_missing_mutations_rejected(self) -> None:
        for name, replacement in (
            ("K", 9),
            ("missing_original_indices", [0, 1]),
            ("requested_profile", "low_v1"),
        ):
            with self.subTest(name=name):
                value = synthetic_raw()
                value["invocations"][1]["result"]["parameters"][name] = replacement
                self.assert_rejected(value)

    def test_type_coerced_benchmark_evidence_is_rejected(self) -> None:
        mutations = {
            "integer-parameter-as-float": lambda value: value["invocations"][0]
                ["result"]["parameters"].__setitem__("K", float(CELL.k)),
            "boolean-option-as-integer": lambda value: value["invocations"][1]
                ["result"]["parameters"].__setitem__("skip_legacy", 1),
            "resolved-parent-as-float": lambda value: value["invocations"][0]
                ["result"]["resolved"].__setitem__(
                    "parent_count",
                    float(value["invocations"][0]["result"]["resolved"]["parent_count"])),
            "resolved-thread-as-boolean": lambda value: value["invocations"][0]
                ["result"]["resolved"].__setitem__("thread_count", True),
            "returncode-as-boolean": lambda value: value["invocations"][0]
                .__setitem__("returncode", False),
            "pinned-cpu-as-boolean": lambda value: value["invocations"][0]
                .__setitem__("pinned_cpu", False),
            "normalized-sample-as-integer": lambda value: value["invocations"][0]
                ["normalized"]["encode"].__setitem__(
                    0, int(value["invocations"][0]["normalized"]["encode"][0])),
            "analysis-count-as-float": lambda value: value["analysis"][CELL.identifier]
                ["encode"].__setitem__("independent_round_count", float(runner.ROUNDS)),
            "statistics-count-as-float": lambda value: value["campaign"]["statistics"]
                .__setitem__("degrees_of_freedom", float(runner.ROUNDS - 1)),
        }
        for name, mutation in mutations.items():
            with self.subTest(name=name):
                value = synthetic_raw()
                mutation(value)
                self.assert_rejected(value)

    def test_type_coerced_retained_text_sizes_are_rejected(self) -> None:
        for role, content_name in (
                ("candidate", "archive_link_recipe_content"),
                ("candidate", "executable_link_recipe_content"),
                ("baseline", "archive_link_recipe_content"),
                ("baseline", "executable_link_recipe_content")):
            with self.subTest(role=role, content=content_name):
                value = synthetic_raw()
                content = value["identities_initial"][f"{role}_build"][content_name]
                content["size"] = float(content["size"])
                synchronize_identity(value)
                self.assert_rejected(value)

    def test_speed_independent_validity_requires_exact_true(self) -> None:
        for replacement in (False, 0, 1, None):
            with self.subTest(replacement=replacement):
                value = synthetic_raw()
                value["validity_is_independent_of_speed"] = replacement
                self.assert_rejected(value)

    def test_successful_round_trip_requires_cross_codec_recovered_digest(self) -> None:
        value = synthetic_raw()
        value["invocations"][0]["result"]["workload_digests"][
            "recovered_originals"] = "f" * 16
        self.assert_rejected(value)

    def test_strict_json_parser_rejects_ambiguous_or_nonfinite_values(self) -> None:
        for data in (
            b'{"duplicate":1,"duplicate":2}',
            b'{"overflow":1e9999}',
            b'{"extension":NaN}',
            b'{"extension":Infinity}',
            b'{"oversized_integer":' + b"9" * 5000 + b"}",
        ):
            with self.subTest(data=data):
                with self.assertRaises(runner.EvidenceError):
                    runner.strict_json_loads(data, "fixture")
        self.assertEqual(
            runner.strict_json_loads(b'{"exact":1,"finite":1.5}', "fixture"),
            {"exact": 1, "finite": 1.5},
        )

    def test_1000_digit_json_number_has_bounded_evidence_error(self) -> None:
        parsed = runner.strict_json_loads(
            b'{"value":' + b"9" * 1000 + b"}", "1000-digit fixture")
        self.assertIs(type(parsed["value"]), int)
        with self.assertRaisesRegex(
                runner.EvidenceError, "outside the finite float range"):
            runner.finite_number(
                parsed["value"], "1000-digit fixture")

    def test_every_fnv_digest_mutation_rejected(self) -> None:
        for name in ("original_data", "transmitted_parity", "recovered_originals"):
            with self.subTest(name=name):
                value = synthetic_raw()
                value["invocations"][1]["result"]["workload_digests"][name] = "f" * 16
                self.assert_rejected(value)

    def test_sample_mutations_rejected(self) -> None:
        mutations = (
            lambda samples: samples.pop(),
            lambda samples: samples.__setitem__(0, 0.0),
            lambda samples: samples.__setitem__(0, float("nan")),
        )
        for mutation in mutations:
            value = synthetic_raw()
            samples = value["invocations"][0]["result"]["metrics"] \
                ["encode_execution"]["samples_us_per_batch_call"]
            mutation(samples)
            self.assert_rejected(value)

    def test_order_identity_and_analysis_mutations_rejected(self) -> None:
        value = synthetic_raw()
        value["invocations"][0]["implementation"] = "candidate"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["invocations"][0]["identity_after"] = {"fixture": {"sha256": "b" * 64}}
        self.assert_rejected(value)
        value = synthetic_raw()
        value["analysis"][CELL.identifier]["encode"]["geometric_speedup"] += 1.0
        self.assert_rejected(value)

    def test_environment_host_and_build_spec_mutations_rejected(self) -> None:
        value = synthetic_raw()
        value["invocations"][0]["environment"]["LD_PRELOAD"] = "/tmp/injected.so"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["campaign"]["child_environment"]["OMP_NUM_THREADS"] = "2"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["host_final"]["benchmark_cpu"]["frequency_policy"] \
            ["scaling_governor"] = "powersave"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["input_specification"].pop("candidate_build_dir")
        self.assert_rejected(value)
        value = synthetic_raw()
        value["invocations"][5]["result"]["resolved"]["backend"] = "scalar"
        value["invocations"][5]["normalized"]["backend"] = "scalar"
        self.assert_rejected(value)
        value = synthetic_raw()
        value["reservation"]["sha256"] = "f" * 64
        for invocation in value["invocations"]:
            invocation["reservation_before"]["sha256"] = "f" * 64
            invocation["reservation_after"]["sha256"] = "f" * 64
        self.assert_rejected(value)

    def test_isolation_mutations_and_active_sibling_rejected(self) -> None:
        value = synthetic_raw()
        value["isolation"]["delta"]["reserved_sibling"]["idle_jiffies"] += 1
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["pair_lease"]["path"] = "/tmp/wrong-pair.lock"
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["accepted"] = False
        self.assert_rejected(value)

        value = synthetic_raw()
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        active_after = copy.deepcopy(after["reserved_sibling"])
        active_after["fields"]["user"] += 1
        active_after["total_jiffies"] += 1
        value["isolation"] = runner.isolation_record(
            0, 1, value["isolation"]["pair_lease"],
            before["monotonic_ns"], after["monotonic_ns"],
            before["benchmark_cpu"], after["benchmark_cpu"],
            before["reserved_sibling"], active_after)
        self.assertFalse(value["isolation"]["accepted"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["before"]["benchmark_cpu"]["cpu"] = 98
        value["isolation"]["after"]["benchmark_cpu"]["cpu"] = 98
        value["isolation"]["delta"]["benchmark_cpu"]["cpu"] = 98
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["before"]["benchmark_cpu"] = []
        self.assert_rejected(value)

        value = synthetic_raw()
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        before["benchmark_cpu"], before["reserved_sibling"] = (
            before["reserved_sibling"], before["benchmark_cpu"])
        after["benchmark_cpu"], after["reserved_sibling"] = (
            after["reserved_sibling"], after["benchmark_cpu"])
        value["isolation"]["delta"]["benchmark_cpu"], \
            value["isolation"]["delta"]["reserved_sibling"] = (
                value["isolation"]["delta"]["reserved_sibling"],
                value["isolation"]["delta"]["benchmark_cpu"])
        self.assert_rejected(value)

        value = synthetic_raw()
        value["invocations"][0]["duration_ns"] = 0
        self.assert_rejected(value)

        value = synthetic_raw()
        value["isolation"]["after"]["monotonic_ns"] = 1_001
        self.assert_rejected(value)

    def test_isolation_replay_does_not_require_original_uid(self) -> None:
        value = synthetic_raw()
        payload = runner.pair_lease_payload(0, 1, uid=12345)
        lease = {
            "device": 10,
            "directory_device": 11,
            "directory_inode": 12,
            "inode": 13,
            "lock": "exclusive_nonblocking_pair_wide",
            "path": str(runner.pair_lease_directory(12345) /
                runner.pair_lease_name(0, 1, uid=12345)),
            "payload": payload,
            "sha256": runner.sha256_bytes(runner.canonical_bytes(payload)),
        }
        before = value["isolation"]["before"]
        after = value["isolation"]["after"]
        value["isolation"] = runner.isolation_record(
            0, 1, lease, before["monotonic_ns"], after["monotonic_ns"],
            before["benchmark_cpu"], after["benchmark_cpu"],
            before["reserved_sibling"], after["reserved_sibling"])
        value = resign(value)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)

    def test_signature_mutation_rejected(self) -> None:
        value = synthetic_raw()
        value["created_utc"] = "changed"
        with self.assertRaises(runner.EvidenceError):
            runner.validate_raw(
                value, None, check_files=False, check_current_inputs=False)

    def test_retained_stdout_mutation_rejected(self) -> None:
        value = synthetic_raw()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, invocation in enumerate(value["invocations"]):
                stdout = json.dumps(invocation["result"]).encode("utf-8")
                stderr = b""
                stdout_path = root / f"{index}.stdout"
                stderr_path = root / f"{index}.stderr"
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(stderr)
                invocation["stdout"] = {
                    "path": stdout_path.name,
                    "size": len(stdout),
                    "sha256": runner.sha256_bytes(stdout),
                }
                invocation["stderr"] = {
                    "path": stderr_path.name,
                    "size": 0,
                    "sha256": runner.sha256_bytes(stderr),
                }
            value = resign(value)
            runner.validate_raw(
                value, root, check_files=True, check_current_inputs=False)
            (root / "0.stdout").write_bytes(b"{}")
            with self.assertRaises(runner.EvidenceError):
                runner.validate_raw(
                    value, root, check_files=True, check_current_inputs=False)

    def test_v7_manifest_binds_complete_identity_and_replays_portably(self) -> None:
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V7)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, invocation in enumerate(value["invocations"]):
                stdout = json.dumps(invocation["result"]).encode("utf-8")
                stdout_path = root / f"{index}.stdout"
                stderr_path = root / f"{index}.stderr"
                stdout_path.write_bytes(stdout)
                stderr_path.write_bytes(b"")
                stdout_path.chmod(0o600)
                stderr_path.chmod(0o600)
                invocation["stdout"] = {
                    "path": stdout_path.name,
                    "size": len(stdout),
                    "sha256": runner.sha256_bytes(stdout),
                }
                invocation["stderr"] = {
                    "path": stderr_path.name,
                    "size": 0,
                    "sha256": runner.sha256_bytes(b""),
                }
            value = resign(value)
            raw_path = root / "raw.json"
            runner.write_json_exclusive(raw_path, value)
            manifest = runner.signed({
                "schema": runner.MANIFEST_SCHEMA_V7,
                "created_utc": "2026-07-16T00:00:00Z",
                "valid": True,
                "validity_is_independent_of_speed": True,
                "raw": {
                    "path": raw_path.name,
                    "size": raw_path.stat().st_size,
                    "sha256": runner.sha256_file(raw_path),
                    "payload_digest": value["digest"],
                },
                "campaign": value["campaign"],
                "host": value["host_initial"],
                "isolation": value["isolation"],
                "reservation": value["reservation"],
                "supervision": value["supervision"],
                "identities": value["identities_initial"],
                "analysis": value["analysis"],
            })
            manifest_path = root / "manifest.json"
            runner.write_json_exclusive(manifest_path, manifest)
            options = argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)
            self.assertEqual(runner.verify_campaign(options), 0)

            options.affinity_binding = root / "affinity-binding.json"
            completed = subprocess.CompletedProcess(
                [], 0, stdout=b"verified\n", stderr=b"")
            with mock.patch.object(
                    runner, "run_process_bounded",
                    return_value=completed) as bounded:
                self.assertEqual(runner.verify_campaign(options), 0)
            verifier_argv = bounded.call_args.args[0]
            supervisor = MODULE_PATH.resolve().parents[3] / \
                "tools/leopard2_affinity_supervisor.py"
            self.assertEqual(
                verifier_argv[:4],
                [sys.executable, "-I", "-S", str(supervisor)])
            self.assertEqual(
                verifier_argv[4:7],
                ["verify-binding", "--binding", str(options.affinity_binding)])

            options.affinity_binding = root / "missing-affinity-binding.json"
            with self.assertRaises(runner.EvidenceError):
                runner.verify_campaign(options)

            options.affinity_binding = None

            edited = copy.deepcopy(manifest)
            edited["isolation"]["accepted"] = False
            edited = resign(edited)
            manifest_path.unlink()
            runner.write_json_exclusive(manifest_path, edited)
            with self.assertRaises(runner.EvidenceError):
                runner.verify_campaign(options)

            edited = copy.deepcopy(manifest)
            edited["validity_is_independent_of_speed"] = 1
            edited = resign(edited)
            manifest_path.unlink()
            runner.write_json_exclusive(manifest_path, edited)
            with self.assertRaises(runner.EvidenceError):
                runner.verify_campaign(options)

            cross_schema = copy.deepcopy(manifest)
            cross_schema["schema"] = runner.MANIFEST_SCHEMA_V2
            cross_schema = resign(cross_schema)
            manifest_path.unlink()
            runner.write_json_exclusive(manifest_path, cross_schema)
            with self.assertRaises(runner.EvidenceError):
                runner.verify_campaign(options)

    def test_legacy_v6_raw_and_manifest_fixtures_remain_replayable(self) -> None:
        plan = load_plan_runner()
        value = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V6)
        runner.validate_raw(
            value, None, check_files=False, check_current_inputs=False)
        runner.validate_failure(
            synthetic_failure(runner.RAW_SCHEMA_V6), Path("."),
            check_files=False)
        with tempfile.TemporaryDirectory() as directory:
            manifest_path = write_complete_evidence_bundle(
                Path(directory), value, runner.MANIFEST_SCHEMA_V6)
            self.assertEqual(runner.verify_campaign(argparse.Namespace(
                manifest=manifest_path, no_current_input_check=True)), 0)
            document, scope, _ = plan.verify_exact_manifest(manifest_path)
            self.assertEqual(document["schema"], runner.MANIFEST_SCHEMA_V6)
            self.assertEqual(
                plan.validate_evidence_scope(scope)["schema"],
                plan.EVIDENCE_SCOPE_SCHEMA_V3)

    def test_literal_v5_v6_candidate_cache_shape_remains_replayable(
            self) -> None:
        historical_cache = {
            "CMAKE_BUILD_TYPE": "Release",
            "CMAKE_CXX_COMPILER": "/usr/bin/compiler",
            "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG",
            "ENABLE_OPENMP": "ON",
            "LEO2_BACKEND_VARIANT": "auto",
            "LEO2_BUILD_BENCHMARKS": "ON",
            "LEO2_BUILD_FUZZERS": "OFF",
            "LEO2_BUILD_TESTS": "OFF",
            "LEO2_ENABLE_CUDA": "OFF",
        }
        for raw_schema in (runner.RAW_SCHEMA_V5, runner.RAW_SCHEMA_V6):
            with self.subTest(raw_schema=raw_schema):
                value = synthetic_raw(raw_schema=raw_schema)
                self.assertEqual(
                    value["identities_initial"]["candidate_build"][
                        "validated_cache"],
                    historical_cache)
                runner.validate_raw(
                    value, None, check_files=False,
                    check_current_inputs=False)

    def test_failed_bundle_binds_retained_invocation_files(self) -> None:
        value = synthetic_raw()
        invocation = copy.deepcopy(value["invocations"][0])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stdout = root / "first.stdout"
            stderr = root / "first.stderr"
            stdout.write_text(json.dumps(invocation["result"]), encoding="utf-8")
            stderr.write_bytes(b"diagnostic")
            stdout.chmod(0o600)
            stderr.chmod(0o600)
            invocation["stdout"] = {
                "path": stdout.name, "size": stdout.stat().st_size,
                "sha256": runner.sha256_file(stdout),
            }
            invocation["stderr"] = {
                "path": stderr.name, "size": stderr.stat().st_size,
                "sha256": runner.sha256_file(stderr),
            }
            failure = runner.signed({
                "schema": runner.FAILURE_SCHEMA,
                "evidence_contract": runner.FAILURE_EVIDENCE_CONTRACT,
                "created_utc": "2026-07-16T00:00:00Z",
                "status": "failed", "valid": False,
                "error_type": "EvidenceError", "error": "fixture failure",
                "campaign": copy.deepcopy(value["campaign"]),
                "host_initial": copy.deepcopy(HOST),
                "reservation": copy.deepcopy(RESERVATION),
                "pair_lease": copy.deepcopy(PAIR_LEASE),
                "isolation": copy.deepcopy(ISOLATION),
                "supervision": copy.deepcopy(value["supervision"]),
                "input_specification": copy.deepcopy(value["input_specification"]),
                "identities_initial": copy.deepcopy(invocation["identity_before"]),
                "executable_snapshots": copy.deepcopy(
                    value["executable_snapshots"]),
                "invocations": [invocation],
                "retained_files": runner.retained_file_records(root),
                "traceback": "fixture traceback",
            })
            failure_path = root / "failure.json"
            runner.write_json_exclusive(failure_path, failure)
            runner.validate_failure(failure, root, check_files=True)
            self.assertEqual(runner.verify_failed_campaign(
                argparse.Namespace(failure=failure_path)), 0)

            duplicate_path = root / "duplicate-failure.json"
            duplicate_path.write_bytes(b'{"schema":"first","schema":"second"}')
            duplicate_path.chmod(0o600)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "duplicate key"):
                runner.verify_failed_campaign(
                    argparse.Namespace(failure=duplicate_path))
            huge_integer_path = root / "huge-integer-failure.json"
            huge_integer_path.write_bytes(
                b'{"untrusted":' + b"9" * 5000 + b"}")
            huge_integer_path.chmod(0o600)
            with self.assertRaises(runner.EvidenceError):
                runner.verify_failed_campaign(
                    argparse.Namespace(failure=huge_integer_path))

            stdout.write_bytes(b"{}")
            semantic_mismatch = copy.deepcopy(failure)
            replacement = {
                "path": stdout.name, "size": stdout.stat().st_size,
                "sha256": runner.sha256_file(stdout),
            }
            semantic_mismatch["invocations"][0]["stdout"] = replacement
            for index, retained in enumerate(semantic_mismatch["retained_files"]):
                if retained["path"] == stdout.name:
                    semantic_mismatch["retained_files"][index] = replacement
            semantic_mismatch = resign(semantic_mismatch)
            with self.assertRaises(runner.EvidenceError):
                runner.validate_failure(semantic_mismatch, root, check_files=True)

            stdout.write_bytes(b"edited")
            with self.assertRaises(runner.EvidenceError):
                runner.validate_failure(failure, root, check_files=True)

            legacy = copy.deepcopy(failure)
            legacy["schema"] = runner.FAILURE_SCHEMA_V2
            legacy.pop("evidence_contract")
            legacy["campaign"].pop("candidate_mode", None)
            legacy.pop("supervision", None)
            legacy.pop("executable_snapshots", None)
            old_identity, old_specification = cmake_fixture_identity(
                runner.RAW_SCHEMA_V2)
            legacy["input_specification"] = old_specification
            legacy["identities_initial"] = old_identity
            legacy["invocations"] = []
            legacy["retained_files"] = runner.retained_file_records(root)
            legacy = resign(legacy)
            runner.validate_failure(legacy, root, check_files=True)

    def test_sealed_executable_survives_same_inode_and_path_replacement(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "benchmark"
            shutil.copyfile("/usr/bin/echo", source)
            source.chmod(0o755)
            expected = runner.artifact_identity(source, "executable")
            descriptor, record = runner.capture_sealed_executable(
                source, expected, "fixture benchmark")
            duplicate = -1
            try:
                original = source.read_bytes()
                with source.open("r+b", buffering=0) as mutable:
                    first = mutable.read(1)
                    mutable.seek(0)
                    mutable.write(bytes((first[0] ^ 1,)))
                    mutable.seek(0)
                    mutable.write(first)
                replacement = root / "replacement"
                shutil.copyfile("/usr/bin/false", replacement)
                replacement.chmod(0o755)
                os.replace(replacement, source)

                self.assertEqual(
                    runner.sealed_executable_identity(
                        descriptor, "fixture benchmark"),
                    record["snapshot"])
                duplicate = runner.duplicate_snapshot_for_execution(
                    descriptor, "fixture benchmark")
                completed = runner.run_process_bounded(
                    [f"/proc/self/fd/{duplicate}", "sealed"],
                    inherited_descriptors=(duplicate,))
                self.assertEqual(completed.returncode, 0)
                self.assertEqual(completed.stdout, b"sealed\n")
                self.assertNotEqual(source.read_bytes(), original)
            finally:
                if duplicate >= 0:
                    os.close(duplicate)
                os.close(descriptor)

    def test_sealed_capture_copy_seal_and_owner_cleanup_faults(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "benchmark"
            shutil.copyfile("/usr/bin/true", source)
            source.chmod(0o755)
            expected = runner.artifact_identity(source, "executable")
            native_create = runner.linux_executable_memfd

            for operation in ("copy", "seal"):
                created = []

                def tracked_create(name: str) -> int:
                    descriptor = native_create(name)
                    created.append(descriptor)
                    return descriptor

                if operation == "copy":
                    fault = mock.patch.object(
                        runner.os, "write",
                        side_effect=OSError(errno.EIO, "injected copy fault"))
                else:
                    native_fcntl = runner.fcntl.fcntl

                    def seal_fault(descriptor, command, argument=0):
                        if command == runner.LINUX_F_ADD_SEALS:
                            raise OSError(errno.EIO, "injected seal fault")
                        return native_fcntl(descriptor, command, argument)

                    fault = mock.patch.object(
                        runner.fcntl, "fcntl", side_effect=seal_fault)
                with self.subTest(operation=operation), \
                     mock.patch.object(
                         runner, "linux_executable_memfd",
                         side_effect=tracked_create), fault, \
                     self.assertRaises(OSError):
                    runner.capture_sealed_executable(
                        source, expected, "faulted benchmark")
                self.assertEqual(len(created), 1)
                with self.assertRaises(OSError):
                    os.fstat(created[0])

        first = os.open("/dev/null", os.O_RDONLY)
        second = os.open("/dev/null", os.O_RDONLY)
        owner = runner.ExecutableSnapshotOwner()
        owner.descriptors = {"baseline": first, "candidate": second}
        native_close = runner.os.close
        closed = []

        def close_fault(descriptor: int) -> None:
            closed.append(descriptor)
            native_close(descriptor)
            if descriptor == first:
                raise OSError(errno.EIO, "injected close fault")

        with mock.patch.object(runner.os, "close", side_effect=close_fault), \
             self.assertRaisesRegex(
                 runner.EvidenceError, "descriptor cleanup failed"):
            owner.close()
        self.assertCountEqual(closed, [first, second])
        self.assertEqual(owner.descriptors, {})
        for descriptor in (first, second):
            with self.assertRaises(OSError):
                os.fstat(descriptor)

    def test_snapshot_owner_close_retries_close_before_effect(self) -> None:
        descriptor = os.open("/dev/null", os.O_RDONLY)
        owner = runner.ExecutableSnapshotOwner()
        owner.descriptors = {"baseline": descriptor}
        try:
            with mock.patch.object(
                    runner.os, "close",
                    side_effect=KeyboardInterrupt("before close")), \
                 self.assertRaisesRegex(
                     runner.EvidenceError, "descriptor cleanup failed"):
                owner.close()
            self.assertEqual(owner.descriptors, {"baseline": descriptor})
            os.fstat(descriptor)

            owner.close()
            self.assertEqual(owner.descriptors, {})
            with self.assertRaises(OSError):
                os.fstat(descriptor)
        finally:
            if owner.descriptors:
                os.close(owner.descriptors.pop("baseline"))

    def test_snapshot_owner_retry_refuses_recycled_descriptor(self) -> None:
        descriptor = os.open("/dev/null", os.O_RDONLY)
        replacement = -1
        owner = runner.ExecutableSnapshotOwner()
        owner.descriptors = {"baseline": descriptor}
        native_close = runner.os.close
        native_fstat = runner.os.fstat
        fstat_calls = 0

        def close_then_raise(value: int) -> None:
            native_close(value)
            raise OSError(errno.EIO, "close completed before diagnostic")

        def interrupt_post_close_probe(value: int):
            nonlocal fstat_calls
            fstat_calls += 1
            if fstat_calls == 2:
                raise KeyboardInterrupt("post-close probe interrupted")
            return native_fstat(value)

        try:
            with mock.patch.object(
                    runner.os, "close", side_effect=close_then_raise), \
                 mock.patch.object(
                     runner.os, "fstat",
                     side_effect=interrupt_post_close_probe), \
                 self.assertRaisesRegex(
                     runner.EvidenceError, "descriptor state probe failed"):
                owner.close()
            self.assertEqual(owner.descriptors, {"baseline": descriptor})
            self.assertIn("baseline", owner.descriptor_identities)

            replacement = os.open("/dev/zero", os.O_RDONLY)
            self.assertEqual(replacement, descriptor)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "identity changed before cleanup"):
                owner.close()
            self.assertEqual(owner.descriptors, {})
            self.assertEqual(owner.descriptor_identities, {})
            os.fstat(replacement)
        finally:
            if replacement >= 0:
                os.close(replacement)
            elif owner.descriptors:
                try:
                    os.close(owner.descriptors.pop("baseline"))
                except OSError:
                    pass

    def test_snapshot_capture_interrupt_after_owner_transfer_is_not_leaked(
            self) -> None:
        class InterruptAfterInsert(dict):
            def __setitem__(self, key, value) -> None:
                super().__setitem__(key, value)
                raise KeyboardInterrupt("after ownership transfer")

        source = Path("/usr/bin/true")
        expected = runner.artifact_identity(source, "executable")
        owner = runner.ExecutableSnapshotOwner()
        owner.descriptors = InterruptAfterInsert()
        descriptor = -1
        try:
            with self.assertRaisesRegex(
                    KeyboardInterrupt, "after ownership transfer"):
                owner.capture("baseline", source, expected)
            descriptor = owner.descriptors["baseline"]
            runner.sealed_executable_identity(
                descriptor, "transferred benchmark")
            owner.close()
            self.assertEqual(owner.descriptors, {})
            with self.assertRaises(OSError):
                os.fstat(descriptor)
        finally:
            if owner.descriptors:
                os.close(owner.descriptors.pop("baseline"))

    def test_snapshot_capture_interrupt_after_source_open_is_not_leaked(
            self) -> None:
        source = Path("/usr/bin/true")
        expected = runner.artifact_identity(source, "executable")
        native_open = runner.os.open
        native_fstat = runner.os.fstat
        opened = []

        def tracked_open(*arguments, **keywords):
            descriptor = native_open(*arguments, **keywords)
            opened.append(descriptor)
            return descriptor

        def interrupt_first_fstat(descriptor: int):
            if opened and descriptor == opened[0]:
                raise KeyboardInterrupt("after source open")
            return native_fstat(descriptor)

        try:
            with mock.patch.object(
                    runner.os, "open", side_effect=tracked_open), \
                 mock.patch.object(
                     runner.os, "fstat", side_effect=interrupt_first_fstat), \
                 self.assertRaisesRegex(
                     KeyboardInterrupt, "after source open"):
                runner.capture_sealed_executable(
                    source, expected, "interrupted benchmark")
            self.assertEqual(len(opened), 1)
            with self.assertRaises(OSError):
                os.fstat(opened[0])
        finally:
            for descriptor in opened:
                try:
                    os.close(descriptor)
                except OSError:
                    pass

    def test_run_child_closes_duplicate_when_argument_setup_interrupts(
            self) -> None:
        retained_descriptor = os.open("/dev/null", os.O_RDONLY)
        execution_descriptor = os.open("/dev/null", os.O_RDONLY)
        owner = runner.ExecutableSnapshotOwner()
        owner.descriptors = {"baseline": retained_descriptor}
        retained_snapshot = {"snapshot": {}}
        try:
            with mock.patch.object(
                    runner, "validate_sealed_executable_record"), \
                 mock.patch.object(owner, "inspect", return_value={}), \
                 mock.patch.object(
                     runner, "duplicate_snapshot_for_execution",
                     return_value=execution_descriptor), \
                 mock.patch.object(
                     runner, "benchmark_arguments",
                     side_effect=KeyboardInterrupt("argument setup")), \
                 self.assertRaisesRegex(KeyboardInterrupt, "argument setup"):
                runner.run_child(
                    "baseline", CELL, 0, 0, CAMPAIGN,
                    {"taskset": "/usr/bin/taskset"}, {}, RESERVATION,
                    Path("/unused"), mock.Mock(), 0, 1.0,
                    {"baseline": retained_snapshot}, owner)
            with self.assertRaises(OSError):
                os.fstat(execution_descriptor)
            os.fstat(retained_descriptor)
        finally:
            try:
                os.close(execution_descriptor)
            except OSError:
                pass
            os.close(retained_descriptor)

    def test_campaign_keyboard_interrupt_closes_both_snapshots_and_reports_cleanup(
        self,
    ) -> None:
        class FakeContext:
            def __init__(self, *args, **kwargs) -> None:
                pass

            def __enter__(self):
                return copy.deepcopy(PAIR_LEASE)

            def __exit__(self, exc_type, exc, traceback_value) -> bool:
                return False

            def validate_current(self) -> None:
                pass

        class FakeEvidence:
            def __init__(self, path: Path) -> None:
                self.path = path

            def validate_current(self) -> None:
                pass

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = argparse.Namespace(
                taskset=Path("/usr/bin/taskset"),
                ldd=Path("/usr/bin/ldd"),
                baseline=Path("/usr/bin/true"),
                candidate=Path("/usr/bin/true"),
                baseline_archive=Path("/usr/bin/true"),
                candidate_archive=Path("/usr/bin/true"),
                baseline_build_dir=root,
                candidate_build_dir=root,
                baseline_source_root=root,
                candidate_source_root=root,
                candidate_commit="1" * 40,
                candidate_mode="auto",
                cpu=0,
                reserved_sibling=1,
                reservation_file=root / "reservation.json",
                reuse=1,
                iterations=3,
                warmup=1,
                timeout=1.0,
                cell=None,
                preset="smoke",
            )
            captured = []

            def capture(_specification, _initial, owner):
                for role in ("baseline", "candidate"):
                    descriptor = os.open("/dev/null", os.O_RDONLY)
                    captured.append(descriptor)
                    owner.descriptors[role] = descriptor
                return {"baseline": {}, "candidate": {}}

            def close_with_diagnostic(owner) -> None:
                for role in sorted(tuple(owner.descriptors), reverse=True):
                    os.close(owner.descriptors.pop(role))
                raise runner.EvidenceError("injected cleanup diagnostic")

            with mock.patch.object(
                    runner, "validate_topology",
                    return_value=({0, 1, 2}, {2})), \
                 mock.patch.object(runner, "host_identity", return_value={}), \
                 mock.patch.object(runner, "PairLease", FakeContext), \
                 mock.patch.object(runner, "Reservation", FakeContext), \
                 mock.patch.object(runner.os, "sched_setaffinity"), \
                 mock.patch.object(runner, "input_snapshot", return_value={}), \
                 mock.patch.object(
                     runner, "capture_campaign_executables",
                     side_effect=capture), \
                 mock.patch.object(
                     runner, "cpu_stat_snapshot",
                     side_effect=KeyboardInterrupt("primary interrupt")), \
                 mock.patch.object(
                     runner.ExecutableSnapshotOwner, "close",
                     side_effect=close_with_diagnostic, autospec=True), \
                 mock.patch.object(
                     runner, "retained_file_records", return_value=[]), \
                 mock.patch.object(runner, "publish_failure_record"), \
                 self.assertRaisesRegex(
                     runner.EvidenceError,
                     "primary interrupt; sealed executable cleanup failed: "
                     "injected cleanup diagnostic"):
                runner._run_campaign_owned(options, FakeEvidence(root))
            self.assertEqual(len(captured), 2)
            for descriptor in captured:
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    def test_execution_descriptor_collision_is_fail_closed_and_high(self) -> None:
        source = runner.artifact_identity(
            Path("/usr/bin/true"), "executable")
        descriptor, _ = runner.capture_sealed_executable(
            Path("/usr/bin/true"), source, "collision benchmark")
        duplicate = -1
        try:
            duplicate = runner.duplicate_snapshot_for_execution(
                descriptor, "collision benchmark")
            self.assertGreaterEqual(
                duplicate, runner.EXECUTION_DESCRIPTOR_FLOOR)
            self.assertNotEqual(duplicate, descriptor)
            native_fcntl = runner.fcntl.fcntl

            def collide(fd, command, argument=0):
                if command == runner.LINUX_F_DUPFD_CLOEXEC:
                    return fd
                return native_fcntl(fd, command, argument)

            with mock.patch.object(
                    runner.fcntl, "fcntl", side_effect=collide), \
                 self.assertRaisesRegex(
                    runner.EvidenceError, "allocation collided"):
                runner.duplicate_snapshot_for_execution(
                    descriptor, "collision benchmark")
            runner.sealed_executable_identity(
                descriptor, "collision benchmark")
        finally:
            if duplicate >= 0:
                os.close(duplicate)
            os.close(descriptor)

    def test_v9_sealed_evidence_is_bound_and_v8_contract_is_historical(
            self) -> None:
        value = synthetic_raw()
        for path, replacement in (
            (("executable_snapshots", "baseline", "snapshot", "sha256"),
             "f" * 64),
            (("executable_snapshots", "candidate", "snapshot", "seals"), 0),
            (("invocations", 0, "execution_protocol"), "mutable-path"),
        ):
            changed = copy.deepcopy(value)
            target = changed
            for component in path[:-1]:
                target = target[component]
            target[path[-1]] = replacement
            self.assert_rejected(changed)

        historical = synthetic_raw(raw_schema=runner.RAW_SCHEMA_V8)
        historical_build = historical["identities_initial"]["candidate_build"]
        historical_compile = historical_build["validated_compile_commands"]
        historical_configuration = \
            historical_compile["effective_build_configuration"]
        self.assertIn("executable_snapshots", historical)
        self.assertEqual(
            historical_compile["schema"], runner.COMPILE_COMMANDS_SCHEMA_V4)
        self.assertEqual(
            historical_configuration["schema"],
            runner.BUILD_CONFIGURATION_RECORD_SCHEMA_V2)
        self.assertEqual(
            historical_configuration["configuration_schema"],
            runner.BUILD_CONFIGURATION_FILE_SCHEMA_V2)
        self.assertNotIn(
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
            historical_configuration["entries"])
        runner.validate_raw(
            historical, None, check_files=False, check_current_inputs=False)
        with tempfile.TemporaryDirectory() as directory:
            historical_manifest = write_complete_evidence_bundle(
                Path(directory), historical, runner.MANIFEST_SCHEMA_V8)
            verified, retained, _, _ = runner.verified_campaign_bundle(
                historical_manifest, no_current_input_check=True)
            self.assertEqual(
                (verified["schema"], retained["schema"]),
                (runner.MANIFEST_SCHEMA_V8, runner.RAW_SCHEMA_V8))
            promotion = load_plan_runner()
            _, historical_scope, _ = promotion.verify_exact_manifest(
                historical_manifest)
            self.assertEqual(
                historical_scope["schema"],
                promotion.EVIDENCE_SCOPE_SCHEMA_V4)

        current_build = value["identities_initial"]["candidate_build"]
        current_compile = current_build["validated_compile_commands"]
        current_configuration = current_compile["effective_build_configuration"]
        self.assertEqual(
            current_compile["schema"], runner.COMPILE_COMMANDS_SCHEMA)
        self.assertEqual(
            current_configuration["schema"],
            runner.BUILD_CONFIGURATION_RECORD_SCHEMA)
        self.assertEqual(
            current_configuration["configuration_schema"],
            runner.BUILD_CONFIGURATION_FILE_SCHEMA)
        self.assertEqual(
            current_configuration["entries"][
                "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"], "OFF")

        # The old body remains internally coherent, but changing only its
        # enclosing schema cannot upgrade its compile/configuration contract.
        relabeled = copy.deepcopy(historical)
        relabeled["schema"] = runner.RAW_SCHEMA
        self.assert_rejected(relabeled)

        # Conversely, current v3 selector closure cannot be presented as the
        # historical v8 contract.
        downgraded = copy.deepcopy(value)
        downgraded["schema"] = runner.RAW_SCHEMA_V8
        self.assert_rejected(downgraded)

        historical_variables = (
            *runner.BUILD_CONFIGURATION_VARIABLES_V2[:-1],
            "LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT",
            runner.BUILD_CONFIGURATION_VARIABLES_V2[-1],
        )
        extended_entries = dict(historical_configuration["entries"])
        extended_entries["LEO2_EXPERIMENT_GENERAL_ONE_LOSS_DIRECT"] = "OFF"
        extended_material = runner.build_configuration_material(
            extended_entries, historical_variables)
        extended_sidecar = (
            f"schema={runner.BUILD_CONFIGURATION_FILE_SCHEMA_V2}\n"
            f"sha256={runner.sha256_bytes(extended_material)}\n"
        ).encode("ascii") + extended_material
        with self.assertRaisesRegex(
                runner.EvidenceError,
                "effective-configuration.*(framing|line count) differs"):
            runner.parse_build_configuration_bytes(
                extended_sidecar, runner.RAW_SCHEMA_V8)

    def test_bounded_file_snapshot_rejects_fifo_without_open_block(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fifo = Path(directory) / "identity.fifo"
            os.mkfifo(fifo, 0o600)
            with self.assertRaises(runner.EvidenceError):
                runner.bounded_file_snapshot(fifo)
            with self.assertRaises(runner.EvidenceError):
                runner.bounded_file_contents_snapshot(fifo)

    def test_evidence_directory_is_owner_only_and_inode_bound(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "evidence"
            previous_umask = os.umask(0o777)
            try:
                evidence = runner.EvidenceDirectory.create_new(output)
            finally:
                os.umask(previous_umask)
            try:
                self.assertEqual(
                    stat.S_IMODE(output.stat().st_mode), 0o700)
                evidence.write_exclusive("nested/result.json", b"stable")
                result = output / "nested/result.json"
                self.assertEqual(
                    stat.S_IMODE(result.parent.stat().st_mode), 0o700)
                self.assertEqual(stat.S_IMODE(result.stat().st_mode), 0o600)
                self.assertEqual(
                    evidence.snapshot("nested/result.json", 64)[1],
                    b"stable")
                result.parent.chmod(0o755)
                with self.assertRaisesRegex(
                        runner.EvidenceError,
                        "child directory is unsafe"):
                    evidence.snapshot("nested/result.json", 64)
                result.parent.chmod(0o700)

                victim = root / "victim"
                victim.write_bytes(b"untouched")
                require_descriptor = evidence.descriptor
                self.assertIsNotNone(require_descriptor)
                os.symlink(
                    "../../victim", "symlink.json",
                    dir_fd=require_descriptor)
                with self.assertRaisesRegex(
                        runner.EvidenceError, "refusing to replace"):
                    evidence.write_exclusive("symlink.json", b"overwrite")
                self.assertEqual(victim.read_bytes(), b"untouched")
            finally:
                evidence.close()

    def test_current_evidence_requires_exact_owner_only_modes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = write_complete_evidence_bundle(
                root, synthetic_raw())
            root.chmod(0o777)
            try:
                with self.assertRaisesRegex(
                        runner.EvidenceError,
                        "evidence directory is not safely owned"):
                    runner.verified_campaign_bundle(
                        manifest, no_current_input_check=True)
            finally:
                root.chmod(0o700)

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = write_complete_evidence_bundle(
                root, synthetic_raw())
            manifest.chmod(0o666)
            with self.assertRaisesRegex(
                    runner.EvidenceError,
                    "bounded safe regular file"):
                runner.verified_campaign_bundle(
                    manifest, no_current_input_check=True)

    def test_versioned_v6_replay_retains_legacy_mode_policy(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = write_complete_evidence_bundle(
                root, synthetic_raw(raw_schema=runner.RAW_SCHEMA_V6),
                runner.MANIFEST_SCHEMA_V6)
            for path in root.iterdir():
                if path.is_file():
                    path.chmod(0o666)
            root.chmod(0o777)
            try:
                verified, raw, _, _ = runner.verified_campaign_bundle(
                    manifest, no_current_input_check=True)
                self.assertEqual(
                    verified["schema"], runner.MANIFEST_SCHEMA_V6)
                self.assertEqual(raw["schema"], runner.RAW_SCHEMA_V6)
            finally:
                root.chmod(0o700)

    def test_evidence_publication_recovers_from_partial_write(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            real_write = os.write
            writes = 0

            def partial_then_enospc(descriptor: int, value: bytes) -> int:
                nonlocal writes
                writes += 1
                if writes == 1:
                    return real_write(descriptor, value[:3])
                raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC))

            try:
                with mock.patch.object(
                        runner.os, "write",
                        side_effect=partial_then_enospc):
                    with self.assertRaisesRegex(
                            runner.EvidenceError, "cannot publish"):
                        evidence.write_exclusive(
                            "artifact.json", b"complete-payload")
                self.assertFalse((output / "artifact.json").exists())
                self.assertFalse(any(
                    name.startswith(".leopard2-evidence-pending-")
                    for name in os.listdir(output)))

                evidence.write_exclusive(
                    "artifact.json", b"complete-payload")
                self.assertEqual(
                    evidence.snapshot("artifact.json", 64)[1],
                    b"complete-payload")
            finally:
                evidence.close()

    def test_evidence_publication_rollback_always_releases_descriptors_and_lock(
            self) -> None:
        for failed_operation in ("unlink", "fstat", "fsync"):
            with self.subTest(failed_operation=failed_operation), \
                    tempfile.TemporaryDirectory() as directory:
                output = Path(directory) / "evidence"
                evidence = runner.EvidenceDirectory.create_new(output)
                descriptors_before = len(os.listdir("/proc/self/fd"))
                real_unlink = os.unlink
                real_fstat = os.fstat
                real_fsync = os.fsync
                regular_fstats = 0

                def injected_unlink(*arguments, **keywords):
                    if failed_operation == "unlink":
                        raise OSError(
                            errno.EIO, "injected rollback unlink EIO")
                    return real_unlink(*arguments, **keywords)

                def injected_fstat(descriptor):
                    nonlocal regular_fstats
                    metadata = real_fstat(descriptor)
                    if stat.S_ISREG(metadata.st_mode):
                        regular_fstats += 1
                        if failed_operation == "fstat" and regular_fstats == 2:
                            raise OSError(
                                errno.EIO, "injected rollback fstat EIO")
                    return metadata

                def injected_fsync(descriptor):
                    if failed_operation == "fsync":
                        raise OSError(
                            errno.EIO, "injected rollback fsync EIO")
                    return real_fsync(descriptor)

                try:
                    with mock.patch.object(
                            runner.os, "write",
                            side_effect=OSError(
                                errno.ENOSPC, "injected primary ENOSPC")), \
                         mock.patch.object(
                             runner.os, "unlink", new=injected_unlink), \
                         mock.patch.object(
                             runner.os, "fstat", new=injected_fstat), \
                         mock.patch.object(
                             runner.os, "fsync", new=injected_fsync), \
                         self.assertRaises(runner.EvidenceError) as captured:
                        evidence.write_exclusive("artifact.json", b"payload")
                    message = str(captured.exception)
                    self.assertIn("injected primary ENOSPC", message)
                    self.assertIn(
                        f"injected rollback {failed_operation} EIO", message)
                    self.assertIn("cleanup also failed", message)
                    self.assertEqual(
                        len(os.listdir("/proc/self/fd")), descriptors_before)

                    # Once the exception-domain and descriptor-count checks
                    # pass, retry also proves the parent flock was released.
                    evidence.write_exclusive("artifact.json", b"payload")
                    self.assertEqual(
                        evidence.snapshot("artifact.json", 64)[1], b"payload")
                finally:
                    evidence.close()

    def test_evidence_publication_resumes_after_crash_staging_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            relative = "artifact.json"
            pending = (
                ".leopard2-evidence-pending-" +
                hashlib.sha256(os.fsencode(relative)).hexdigest()
            )
            try:
                staging = output / pending
                staging.write_bytes(b"partial")
                staging.chmod(0o600)
                evidence.write_exclusive(relative, b"complete")
                self.assertFalse(staging.exists())
                self.assertEqual((output / relative).read_bytes(), b"complete")

                with self.assertRaisesRegex(
                        runner.EvidenceError, "refusing to replace"):
                    evidence.write_exclusive(relative, b"replacement")
                self.assertEqual((output / relative).read_bytes(), b"complete")
            finally:
                evidence.close()

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            relative = "artifact.json"
            pending = (
                ".leopard2-evidence-pending-" +
                hashlib.sha256(os.fsencode(relative)).hexdigest()
            )
            try:
                staging = output / pending
                staging.write_bytes(b"opened-before-fchmod")
                staging.chmod(0o000)
                evidence.write_exclusive(relative, b"complete")
                self.assertFalse(staging.exists())
                self.assertEqual((output / relative).read_bytes(), b"complete")
            finally:
                evidence.close()

    def test_failure_publication_rejects_preexisting_entries_and_validates_storage(
            self) -> None:
        for kind in ("regular", "symlink", "fifo"):
            with self.subTest(kind=kind), \
                    tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                output = root / "evidence"
                evidence = runner.EvidenceDirectory.create_new(output)
                victim = root / "victim"
                victim.write_bytes(b"untouched")
                failure_path = output / "failure.json"
                if kind == "regular":
                    failure_path.write_bytes(b"stale")
                    failure_path.chmod(0o600)
                elif kind == "symlink":
                    failure_path.symlink_to(victim)
                else:
                    os.mkfifo(failure_path, 0o600)
                try:
                    with self.assertRaisesRegex(
                            runner.EvidenceError, "refusing to replace"):
                        runner.publish_failure_record(
                            evidence,
                            synthetic_failure(runner.RAW_SCHEMA),
                            output)
                    self.assertEqual(victim.read_bytes(), b"untouched")
                finally:
                    evidence.close()

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            failure = synthetic_failure(runner.RAW_SCHEMA)
            try:
                stored = runner.publish_failure_record(
                    evidence, failure, output)
                self.assertEqual(stored, failure)
                self.assertEqual(
                    runner.strict_json_loads(
                        evidence.snapshot("failure.json")[1],
                        "stored failure"),
                    failure)
            finally:
                evidence.close()

    def test_failure_publication_rejects_postwrite_snapshot_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            real_snapshot = evidence.snapshot

            def changed_snapshot(*arguments, **keywords):
                metadata, value = real_snapshot(*arguments, **keywords)
                return metadata, value + b" "

            try:
                with mock.patch.object(
                        evidence, "snapshot",
                        side_effect=changed_snapshot):
                    with self.assertRaisesRegex(
                            runner.EvidenceError,
                            "stored failure record differs"):
                        runner.publish_failure_record(
                            evidence,
                            synthetic_failure(runner.RAW_SCHEMA),
                            output)
            finally:
                evidence.close()

    def test_evidence_snapshot_rejects_nested_parent_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            try:
                evidence.write_exclusive("nested/result.json", b"stable")

                def replace_parent() -> None:
                    (output / "nested").rename(output / "retained-nested")
                    (output / "nested").mkdir(mode=0o700)
                    (output / "nested/result.json").write_bytes(b"stable")

                with self.assertRaisesRegex(
                        runner.EvidenceError, "parent directory was replaced"):
                    evidence.snapshot(
                        "nested/result.json", 64,
                        mutation_hook=replace_parent)
            finally:
                evidence.close()

    def test_evidence_directory_validation_failures_do_not_leak_fds(
        self,
    ) -> None:
        def fd_count() -> int:
            return len(os.listdir("/proc/self/fd"))

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            child = root / "child"
            child.mkdir(mode=0o700)
            baseline = fd_count()
            native_stat = runner.os.stat

            def mismatched_child(path, *args, **kwargs):
                result = native_stat(path, *args, **kwargs)
                if path == "child" and \
                        kwargs.get("follow_symlinks") is False:
                    values = list(result)
                    values[1] += 1
                    return os.stat_result(values)
                return result

            with mock.patch.object(
                    runner.os, "stat", side_effect=mismatched_child):
                for _ in range(32):
                    with self.assertRaises(runner.EvidenceError):
                        runner.EvidenceDirectory._open_absolute_directory(child)
                    self.assertEqual(fd_count(), baseline)

            evidence = runner.EvidenceDirectory.open_existing(root)
            try:
                held_baseline = fd_count()
                with mock.patch.object(
                        runner.EvidenceDirectory,
                        "_validate_child_directory",
                        side_effect=runner.EvidenceError(
                            "injected child validation failure")):
                    for _ in range(32):
                        with self.assertRaises(runner.EvidenceError):
                            evidence._open_parent("child/result", create=False)
                        self.assertEqual(fd_count(), held_baseline)
            finally:
                evidence.close()
            self.assertEqual(fd_count(), baseline)

            for index in range(32):
                output = root / f"create-failure-{index}"
                with mock.patch.object(
                        runner.EvidenceDirectory, "validate_current",
                        side_effect=runner.EvidenceError(
                            "injected final validation failure")), \
                     self.assertRaises(runner.EvidenceError):
                    runner.EvidenceDirectory.create_new(output)
                self.assertEqual(fd_count(), baseline)

            for _ in range(32):
                with mock.patch.object(
                        runner.EvidenceDirectory, "validate_current",
                        side_effect=runner.EvidenceError(
                            "injected open validation failure")), \
                     self.assertRaises(runner.EvidenceError):
                    runner.EvidenceDirectory.open_existing(root)
                self.assertEqual(fd_count(), baseline)

    def test_evidence_directory_rejects_parent_and_lifetime_replacement(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            victim = root / "victim"
            victim.mkdir(mode=0o700)
            redirected = root / "redirected"
            redirected.symlink_to(victim.name)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "symlinked ancestor"):
                runner.EvidenceDirectory.create_new(
                    redirected / "evidence")

            output = root / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            retained = root / "retained-evidence"
            output.rename(retained)
            output.mkdir(mode=0o700)
            try:
                with self.assertRaisesRegex(
                        runner.EvidenceError, "replaced"):
                    evidence.write_exclusive("result.json", b"unsafe")
                self.assertFalse((output / "result.json").exists())
                self.assertFalse((retained / "result.json").exists())
            finally:
                evidence.close()

    def test_evidence_snapshot_rejects_in_place_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "evidence"
            evidence = runner.EvidenceDirectory.create_new(output)
            try:
                evidence.write_exclusive(
                    "result.json", b"A" * (2 * 1024 * 1024))

                def mutate() -> None:
                    descriptor = os.open(
                        output / "result.json",
                        os.O_WRONLY | getattr(os, "O_CLOEXEC", 0))
                    try:
                        os.pwrite(descriptor, b"B", 32)
                        os.fsync(descriptor)
                    finally:
                        os.close(descriptor)

                with self.assertRaisesRegex(
                        runner.EvidenceError, "changed while reading"):
                    evidence.snapshot(
                        "result.json", 4 * 1024 * 1024,
                        mutation_hook=mutate)
            finally:
                evidence.close()

    def test_process_group_reap_never_uses_unbounded_wait(self) -> None:
        class NeverReaps:
            pid = 123456
            returncode = None

            def __init__(self) -> None:
                self.calls: list[float | None] = []

            def wait(self, timeout: float | None = None) -> int:
                self.calls.append(timeout)
                raise subprocess.TimeoutExpired(("never",), timeout)

        process = NeverReaps()
        with mock.patch.object(os, "killpg"):
            reaped, returncode = runner.terminate_process_group_bounded(process)
        self.assertFalse(reaped)
        self.assertEqual(returncode, -9)
        self.assertEqual(process.calls, [5.0])
        self.assertNotIn(None, process.calls)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux retained pidfd test")
    def test_pidfd_retention_rejects_same_tick_pid_reuse_by_proc_inode(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "first-task"
            second = Path(directory) / "second-task"
            first.mkdir()
            second.mkdir()
            first_fd = os.open(first, os.O_RDONLY | os.O_DIRECTORY)
            second_fd = os.open(second, os.O_RDONLY | os.O_DIRECTORY)
            pidfd = os.open("/dev/null", os.O_RDONLY)
            record = (os.getpid(), 100, 100, 424242, "S")
            first_metadata = os.fstat(first_fd)
            expected = (
                12345, 424242, first_metadata.st_dev, first_metadata.st_ino)
            with mock.patch.object(
                    runner, "_open_proc_task_directory",
                    side_effect=(first_fd, second_fd)), \
                 mock.patch.object(
                     runner, "_proc_record_from_task_directory",
                     return_value=record), \
                 mock.patch.object(
                     runner, "_linux_pidfd_open",
                     return_value=pidfd), \
                 self.assertRaisesRegex(
                     runner.EvidenceError,
                     "changed while retaining"):
                runner._retain_linux_process(12345, expected)
            for descriptor in (first_fd, second_fd, pidfd):
                with self.assertRaises(OSError):
                    os.fstat(descriptor)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux procfs inode snapshot test")
    def test_proc_snapshots_bind_the_held_no_follow_directory_inode(self) -> None:
        record = (os.getpid(), 100, 100, 424242, "S")
        for emergency in (False, True):
            with self.subTest(emergency=emergency), \
                    tempfile.TemporaryDirectory() as directory:
                first = Path(directory) / "first-task"
                second = Path(directory) / "second-task"
                first.mkdir()
                second.mkdir()
                first_fd = os.open(first, os.O_RDONLY | os.O_DIRECTORY)
                first_metadata = os.fstat(first_fd)
                second_metadata = os.stat(second)
                open_name = (
                    "_emergency_open_proc_task_directory"
                    if emergency else "_open_proc_task_directory")
                record_name = (
                    "_emergency_proc_process_record"
                    if emergency else "_proc_process_record")
                with mock.patch.object(
                        runner, open_name, return_value=first_fd), \
                     mock.patch.object(
                         runner, "_proc_record_from_task_directory",
                         return_value=record), \
                     mock.patch.object(
                         runner.os, "stat", return_value=second_metadata):
                    self.assertIsNone(getattr(runner, record_name)(12345))
                with self.assertRaises(OSError):
                    os.fstat(first_fd)

                retained_fd = os.open(first, os.O_RDONLY | os.O_DIRECTORY)
                with mock.patch.object(
                        runner, open_name, return_value=retained_fd), \
                     mock.patch.object(
                         runner, "_proc_record_from_task_directory",
                         return_value=record), \
                     mock.patch.object(
                         runner.os, "stat", return_value=first_metadata):
                    self.assertEqual(
                        getattr(runner, record_name)(12345),
                        (*record, first_metadata.st_dev,
                         first_metadata.st_ino))
                with self.assertRaises(OSError):
                    os.fstat(retained_fd)

    def _check_post_retention_pid_reuse(
            self, new_starttime: int, emergency: bool) -> None:
        containment = runner.LinuxDescendantContainment()
        pid = 12345
        old_identity = (pid, 424242, 10, 20)
        new_identity = (pid, new_starttime, 10, 21)
        old_pidfd = os.open("/dev/null", os.O_RDONLY)
        old_procfd = os.open("/dev/null", os.O_RDONLY)
        new_pidfd = os.open("/dev/null", os.O_RDONLY)
        new_procfd = os.open("/dev/null", os.O_RDONLY)
        old_handle = runner._RetainedLinuxProcess(
            identity=old_identity, pidfd=old_pidfd, procfd=old_procfd,
            proc_identity=old_identity[2:],
            record=(containment.runner_pid, pid, pid, old_identity[1], "Z"))
        new_record = (
            containment.runner_pid, pid, pid, new_identity[1], "S")
        new_handle = runner._RetainedLinuxProcess(
            identity=new_identity, pidfd=new_pidfd, procfd=new_procfd,
            proc_identity=new_identity[2:], record=new_record)
        containment.handles[old_identity] = old_handle
        snapshot = {pid: (*new_record, *new_identity[2:])}
        retain_name = (
            "_emergency_retain_linux_process"
            if emergency else "_retain_linux_process")
        signal_name = (
            "_emergency_pidfd_signal"
            if emergency else "_linux_pidfd_signal")
        try:
            self.assertEqual(
                containment._candidate_identities(snapshot), {new_identity})
            with mock.patch.object(
                    runner, retain_name, return_value=new_handle) as retain:
                self.assertEqual(
                    containment._discover(snapshot, emergency=emergency),
                    {old_identity, new_identity})
            retain.assert_called_once_with(pid, new_identity)
            with mock.patch.object(
                    runner, "_pidfd_exited",
                    side_effect=lambda descriptor: descriptor == old_pidfd), \
                 mock.patch.object(runner, signal_name) as retained_signal:
                containment._signal_handles(emergency=emergency)
            retained_signal.assert_called_once_with(
                new_pidfd, signal.SIGKILL)
        finally:
            if new_identity not in containment.handles:
                new_handle.close()
            containment._close_handles()

    def test_containment_tracks_different_tick_pid_reuse_as_new_lifetime(
            self) -> None:
        for emergency in (False, True):
            with self.subTest(emergency=emergency):
                self._check_post_retention_pid_reuse(424243, emergency)

    def test_containment_tracks_same_tick_pid_reuse_by_proc_inode(self) -> None:
        for emergency in (False, True):
            with self.subTest(emergency=emergency):
                self._check_post_retention_pid_reuse(424242, emergency)

    def test_containment_signals_only_lifetime_retained_pidfds(self) -> None:
        containment = runner.LinuxDescendantContainment()
        identity = (12345, 67890, 1, 2)
        containment.handles[identity] = runner._RetainedLinuxProcess(
            identity=identity, pidfd=91, procfd=92,
            proc_identity=(1, 2),
            record=(containment.runner_pid, 12345, 12345, 67890, "S"))
        with mock.patch.object(
                runner, "_pidfd_exited", return_value=False), \
             mock.patch.object(
                 runner, "_linux_pidfd_signal") as retained_signal, \
             mock.patch.object(
                 runner, "_linux_pidfd_open",
                 side_effect=AssertionError(
                     "signal path reopened a numeric PID")):
            containment._signal_handles(emergency=False)
        retained_signal.assert_called_once_with(91, signal.SIGKILL)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux retained leader test")
    def test_leader_is_not_reaped_before_pidfd_tree_quiescence(self) -> None:
        real_wait = subprocess.Popen.wait
        wait_observations: list[bool] = []

        def observed_wait(
            process: subprocess.Popen[bytes],
            timeout: float | None = None,
        ) -> int:
            wait_observations.append(
                Path("/proc", str(process.pid)).exists())
            return real_wait(process, timeout=timeout)

        with mock.patch.object(
                runner.subprocess.Popen, "wait",
                new=observed_wait), \
             mock.patch.object(
                 runner.subprocess.Popen, "poll",
                 side_effect=AssertionError(
                     "leader was polled/reaped before quiescence")):
            completed = runner.run_process_bounded(
                [sys.executable, "-c", "raise SystemExit(0)"],
                timeout=2.0, max_stdout=1024, max_stderr=1024)
        self.assertEqual(completed.returncode, 0)
        self.assertEqual(wait_observations, [True])

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux emergency retained-handle test")
    def test_post_popen_attachment_fault_uses_retained_emergency_cleanup(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            marker = Path(directory) / "escaped-marker"
            child = (
                "import pathlib,sys,time\n"
                "time.sleep(.6)\n"
                "pathlib.Path(sys.argv[1]).write_text("
                "'escaped',encoding='ascii')\n"
            )
            subreaper_before = runner._get_child_subreaper()
            descriptors_before = len(os.listdir("/proc/self/fd"))
            with mock.patch.object(
                    runner.LinuxDescendantContainment, "attach",
                    side_effect=runner.EvidenceError(
                        "injected attachment fault")), \
                 self.assertRaisesRegex(
                     runner.EvidenceError, "injected attachment fault"):
                runner.run_process_bounded(
                    [sys.executable, "-c", child, str(marker)],
                    timeout=2.0, max_stdout=1024, max_stderr=1024)
            time.sleep(.7)
            self.assertFalse(marker.exists())
            self.assertEqual(
                runner._get_child_subreaper(), subreaper_before)
            self.assertEqual(
                len(os.listdir("/proc/self/fd")), descriptors_before)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_timeout_kills_setsid_double_fork_descendant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "escaped-marker"
            ready = root / "escaped-ready"
            subreaper_before = runner._get_child_subreaper()
            child = (
                "import os,pathlib,sys,time\n"
                "pid=os.fork()\n"
                "if pid == 0:\n"
                " os.setsid()\n"
                " daemon=os.fork()\n"
                " if daemon != 0: os._exit(0)\n"
                " ready=pathlib.Path(sys.argv[2])\n"
                " temporary=ready.with_name("
                "ready.name+'.tmp-'+str(os.getpid()))\n"
                " temporary.write_text("
                "str(os.getpid()),encoding='ascii')\n"
                " os.replace(temporary,ready)\n"
                " time.sleep(1.5)\n"
                " open(sys.argv[1], 'w').write('escaped')\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+5\n"
                "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
                " time.sleep(.01)\n"
                "while True: time.sleep(1)\n"
            )
            with self.assertRaisesRegex(runner.EvidenceError, "exceeded"):
                runner.run_process_bounded(
                    [sys.executable, "-c", child, str(marker), str(ready)],
                    timeout=0.8, max_stdout=1024, max_stderr=1024)
            escaped_pid = int(ready.read_text(encoding="utf-8"))
            self.assertEqual(runner._get_child_subreaper(), subreaper_before)
            time.sleep(1.0)
            self.assertFalse(marker.exists())
            self.assertFalse(Path("/proc", str(escaped_pid)).exists())

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_successful_leader_cannot_leave_daemon_descendant(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            marker = root / "daemon-marker"
            ready = root / "daemon-ready"
            child = (
                "import os,pathlib,sys,time\n"
                "pid=os.fork()\n"
                "if pid == 0:\n"
                " os.setsid()\n"
                " daemon=os.fork()\n"
                " if daemon != 0: os._exit(0)\n"
                " null=os.open('/dev/null', os.O_WRONLY)\n"
                " os.dup2(null, 1);os.dup2(null, 2);os.close(null)\n"
                " ready=pathlib.Path(sys.argv[2])\n"
                " temporary=ready.with_name("
                "ready.name+'.tmp-'+str(os.getpid()))\n"
                " temporary.write_text("
                "str(os.getpid()),encoding='ascii')\n"
                " os.replace(temporary,ready)\n"
                " time.sleep(1.0)\n"
                " open(sys.argv[1], 'w').write('escaped')\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+5\n"
                "while not os.path.exists(sys.argv[2]) and time.monotonic()<deadline:\n"
                " time.sleep(.01)\n"
            )
            completed = runner.run_process_bounded(
                [sys.executable, "-c", child, str(marker), str(ready)],
                timeout=2.0, max_stdout=1024, max_stderr=1024)
            self.assertEqual(completed.returncode, 0)
            escaped_pid = int(ready.read_text(encoding="utf-8"))
            time.sleep(1.1)
            self.assertFalse(marker.exists())
            self.assertFalse(Path("/proc", str(escaped_pid)).exists())

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux descendant containment test")
    def test_zero_exit_leader_with_retained_pipe_is_reaped_immediately(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pid_path = Path(directory) / "retained-pipe.pid"
            child = (
                "import os,pathlib,sys,time\n"
                "pid=os.fork()\n"
                "if pid==0:\n"
                " os.setsid()\n"
                " ready=pathlib.Path(sys.argv[1])\n"
                " temporary=ready.with_name("
                "ready.name+'.tmp-'+str(os.getpid()))\n"
                " temporary.write_text("
                "str(os.getpid()),encoding='ascii')\n"
                " os.replace(temporary,ready)\n"
                " time.sleep(30)\n"
                " os._exit(0)\n"
                "deadline=time.monotonic()+5\n"
                "while not os.path.exists(sys.argv[1]) and "
                "time.monotonic()<deadline: time.sleep(.01)\n"
                "sys.exit(0 if os.path.exists(sys.argv[1]) else 91)\n"
            )
            started = time.monotonic()
            completed = runner.run_process_bounded(
                [sys.executable, "-c", child, str(pid_path)],
                timeout=4.0, max_stdout=1024, max_stderr=1024)
            elapsed = time.monotonic() - started
            retained_pid = int(pid_path.read_text(encoding="ascii"))
            self.assertEqual(completed.returncode, 0)
            self.assertLess(elapsed, 2.0)
            self.assertFalse(Path("/proc", str(retained_pid)).exists())

    @unittest.skipUnless(
        sys.platform.startswith("linux") and hasattr(fcntl, "F_SETPIPE_SZ"),
        "Linux small-pipe capture test",
    )
    def test_small_capture_pipes_are_drained_concurrently(self) -> None:
        real_popen = subprocess.Popen

        def small_pipe_popen(*arguments, **keywords):
            process = real_popen(*arguments, **keywords)
            self.assertIsNotNone(process.stdout)
            self.assertIsNotNone(process.stderr)
            fcntl.fcntl(process.stdout.fileno(), fcntl.F_SETPIPE_SZ, 4096)
            fcntl.fcntl(process.stderr.fileno(), fcntl.F_SETPIPE_SZ, 4096)
            return process

        payload_size = 256 * 1024
        child = (
            "import os\n"
            f"value=b'x'*{payload_size}\n"
            "offset=0\n"
            "while offset<len(value): offset+=os.write(1,value[offset:])\n"
            "value=b'y'*len(value)\n"
            "offset=0\n"
            "while offset<len(value): offset+=os.write(2,value[offset:])\n"
        )
        with mock.patch.object(
                runner.subprocess, "Popen", side_effect=small_pipe_popen):
            completed = runner.run_process_bounded(
                [sys.executable, "-c", child], timeout=4.0,
                max_stdout=payload_size, max_stderr=payload_size)
        self.assertEqual(completed.returncode, 0)
        self.assertEqual(completed.stdout, b"x" * payload_size)
        self.assertEqual(completed.stderr, b"y" * payload_size)

    @unittest.skipUnless(sys.platform.startswith("linux"),
                         "Linux bounded capture test")
    def test_closed_capture_pipes_do_not_extend_command_timeout(self) -> None:
        child = (
            "import os,time\n"
            "os.close(1)\n"
            "os.close(2)\n"
            "time.sleep(30)\n"
        )
        started = time.monotonic()
        with self.assertRaisesRegex(runner.EvidenceError, "exceeded"):
            runner.run_process_bounded(
                [sys.executable, "-c", child], timeout=0.2,
                max_stdout=1024, max_stderr=1024)
        self.assertLess(time.monotonic() - started, 1.0)

    def test_descendant_containment_fails_closed_before_spawn(self) -> None:
        with mock.patch.object(runner.sys, "platform", "not-linux"), \
             mock.patch.object(runner.subprocess, "Popen") as popen, \
             self.assertRaisesRegex(runner.EvidenceError, "requires Linux"):
            runner.run_process_bounded(
                [sys.executable, "-c", "pass"], timeout=1.0,
                max_stdout=1024, max_stderr=1024)
        popen.assert_not_called()

        for unavailable in (
                "_validate_linux_pidfd_support", "_get_child_subreaper",
                "_proc_process_snapshot"):
            with self.subTest(unavailable=unavailable), \
                 mock.patch.object(
                     runner, unavailable,
                     side_effect=runner.EvidenceError(unavailable + " unavailable")), \
                 mock.patch.object(runner.subprocess, "Popen") as popen, \
                 self.assertRaisesRegex(runner.EvidenceError, "unavailable"):
                runner.run_process_bounded(
                    [sys.executable, "-c", "pass"], timeout=1.0,
                    max_stdout=1024, max_stderr=1024)
            popen.assert_not_called()

    def test_reservation_is_locked_and_canonical(self) -> None:
        payload = {
            "benchmark_cpu": 0,
            "nonce": "fixture-nonce",
            "owner": "unit test",
            "reserved_sibling": 1,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reservation.json"
            path.write_bytes(runner.canonical_bytes(payload))
            with runner.Reservation(path, 0, 1) as identity:
                self.assertEqual(identity["payload"], payload)
                with self.assertRaises(runner.EvidenceError):
                    with runner.Reservation(path, 0, 1):
                        pass
            path.write_bytes(runner.canonical_bytes(payload) + b"\n")
            with self.assertRaises(runner.EvidenceError):
                with runner.Reservation(path, 0, 1):
                    pass

    def test_pair_lease_serializes_different_reservation_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            with runner.PairLease(0, 1, root=root) as identity:
                self.assertEqual(identity["payload"], runner.pair_lease_payload(0, 1))
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(1, 0, root=root):
                        pass
            path = root / runner.pair_lease_name(0, 1)
            path.chmod(0o644)
            with self.assertRaises(runner.EvidenceError):
                with runner.PairLease(0, 1, root=root):
                    pass

    def test_pair_lease_detects_unlink_and_replacement(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            first = runner.PairLease(0, 1, root=root)
            with first:
                first.path.unlink()
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(0, 1, root=root):
                        pass
                with self.assertRaises(runner.EvidenceError):
                    first.validate_current()

    def test_pair_lease_failed_enter_releases_kernel_lease(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            root.mkdir(mode=0o700)
            lease = runner.PairLease(0, 1, root=root)
            lease.path.mkdir(mode=0o700)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "not a regular file"):
                lease.__enter__()
            self.assertIsNone(lease.descriptor)
            self.assertIsNone(lease.kernel_socket)
            lease.path.rmdir()
            with runner.PairLease(1, 0, root=root):
                pass

    @unittest.skipUnless(Path("/proc/self/fd").is_dir(),
                         "descriptor-count regression needs procfs")
    def test_pair_lease_enter_recovers_from_first_fstat_interrupt(
            self) -> None:
        real_fstat = os.fstat
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            descriptors_before = len(os.listdir("/proc/self/fd"))
            lease = runner.PairLease(0, 1, root=root)
            injected = False

            def interrupt_first_lease_status(descriptor: int):
                nonlocal injected
                if descriptor == lease.descriptor and \
                        lease.descriptor_identity is None and not injected:
                    injected = True
                    raise KeyboardInterrupt(
                        "injected post-open fstat interruption")
                return real_fstat(descriptor)

            with mock.patch.object(
                    runner.os, "fstat",
                    side_effect=interrupt_first_lease_status), \
                 self.assertRaisesRegex(
                     KeyboardInterrupt, "post-open fstat interruption"):
                lease.__enter__()
            self.assertTrue(injected)
            self.assertIsNone(lease.descriptor)
            self.assertIsNone(lease.descriptor_identity)
            self.assertIsNone(lease.kernel_socket)
            with runner.PairLease(1, 0, root=root):
                pass
            self.assertEqual(
                len(os.listdir("/proc/self/fd")), descriptors_before)

    @unittest.skipUnless(Path("/proc/self/fd").is_dir(),
                         "descriptor-count regression needs procfs")
    def test_pair_lease_failed_enter_retains_fd_after_close_failure(
            self) -> None:
        real_close = os.close
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            root.mkdir(mode=0o700)
            path = root / runner.pair_lease_name(0, 1)
            path.write_bytes(b"noncanonical")
            path.chmod(0o600)
            descriptors_before = len(os.listdir("/proc/self/fd"))
            lease = runner.PairLease(0, 1, root=root)
            injected = False

            def fail_before_close(descriptor: int) -> None:
                nonlocal injected
                if descriptor == lease.descriptor and not injected:
                    injected = True
                    raise OSError(errno.EIO, "injected close failure")
                real_close(descriptor)

            with mock.patch.object(
                    runner.os, "close", side_effect=fail_before_close), \
                 self.assertRaisesRegex(
                     runner.EvidenceError, "unexpected or noncanonical"):
                lease.__enter__()
            self.assertTrue(injected)
            self.assertIsNotNone(lease.descriptor)
            self.assertIsNotNone(lease.descriptor_identity)
            os.fstat(lease.descriptor)
            self.assertIsNone(lease.kernel_socket)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "already leased"):
                with runner.PairLease(1, 0, root=root):
                    pass
            lease.__exit__(None, None, None)
            self.assertIsNone(lease.descriptor)
            self.assertIsNone(lease.descriptor_identity)
            path.unlink()
            with runner.PairLease(1, 0, root=root):
                pass
            self.assertEqual(
                len(os.listdir("/proc/self/fd")), descriptors_before)

    @unittest.skipUnless(Path("/proc/self/fd").is_dir(),
                         "descriptor-count regression needs procfs")
    def test_pair_lease_exit_retains_fd_after_close_failure(self) -> None:
        real_close = os.close
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            descriptors_before = len(os.listdir("/proc/self/fd"))
            lease = runner.PairLease(0, 1, root=root)
            lease.__enter__()
            retained_descriptor = lease.descriptor
            injected = False

            def fail_before_close(descriptor: int) -> None:
                nonlocal injected
                if descriptor == retained_descriptor and not injected:
                    injected = True
                    raise OSError(errno.EIO, "injected close failure")
                real_close(descriptor)

            with mock.patch.object(
                    runner.os, "close", side_effect=fail_before_close), \
                 self.assertRaisesRegex(OSError, "injected close failure"):
                lease.__exit__(None, None, None)
            self.assertTrue(injected)
            self.assertEqual(lease.descriptor, retained_descriptor)
            self.assertIsNotNone(lease.descriptor_identity)
            os.fstat(retained_descriptor)
            self.assertIsNone(lease.kernel_socket)
            with self.assertRaisesRegex(
                    runner.EvidenceError, "already leased"):
                with runner.PairLease(1, 0, root=root):
                    pass
            lease.__exit__(None, None, None)
            with runner.PairLease(1, 0, root=root):
                pass
            self.assertEqual(
                len(os.listdir("/proc/self/fd")), descriptors_before)

    def test_pair_lease_exit_clears_already_closed_ebadf(self) -> None:
        real_close = os.close
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            lease = runner.PairLease(0, 1, root=root)
            lease.__enter__()
            retained_descriptor = lease.descriptor
            injected = False

            def close_then_report_ebadf(descriptor: int) -> None:
                nonlocal injected
                if descriptor == retained_descriptor and not injected:
                    injected = True
                    real_close(descriptor)
                    raise OSError(errno.EBADF, "injected already-closed fd")
                real_close(descriptor)

            with mock.patch.object(
                    runner.os, "close", side_effect=close_then_report_ebadf), \
                 self.assertRaisesRegex(OSError, "already-closed"):
                lease.__exit__(None, None, None)
            self.assertTrue(injected)
            self.assertIsNone(lease.descriptor)
            self.assertIsNone(lease.descriptor_identity)
            self.assertIsNone(lease.kernel_socket)
            with runner.PairLease(1, 0, root=root):
                pass

    def test_pair_lease_exit_refuses_recycled_descriptor(self) -> None:
        real_close = os.close
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            foreign = Path(directory) / "foreign"
            foreign.write_bytes(b"foreign")
            lease = runner.PairLease(0, 1, root=root)
            lease.__enter__()
            retained_descriptor = lease.descriptor
            recycled: list[int] = []

            def close_and_recycle(descriptor: int) -> None:
                if descriptor == retained_descriptor and not recycled:
                    real_close(descriptor)
                    replacement = os.open(foreign, os.O_RDONLY)
                    self.assertEqual(replacement, descriptor)
                    recycled.append(replacement)
                    raise OSError(errno.EIO, "injected post-close failure")
                real_close(descriptor)

            try:
                with mock.patch.object(
                        runner.os, "close",
                        side_effect=close_and_recycle), \
                     self.assertRaisesRegex(
                         runner.EvidenceError, "descriptor number was recycled"):
                    lease.__exit__(None, None, None)
                self.assertEqual(recycled, [retained_descriptor])
                self.assertIsNone(lease.descriptor)
                self.assertIsNone(lease.descriptor_identity)
                self.assertIsNone(lease.kernel_socket)
                os.fstat(recycled[0])
            finally:
                for descriptor in recycled:
                    real_close(descriptor)
            with runner.PairLease(1, 0, root=root):
                pass

    def test_pair_lease_exit_releases_kernel_socket_on_status_interrupt(
            self) -> None:
        real_fstat = os.fstat
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            lease = runner.PairLease(0, 1, root=root)
            lease.__enter__()
            retained_descriptor = lease.descriptor

            def interrupt_status(descriptor: int):
                if descriptor == retained_descriptor:
                    raise KeyboardInterrupt(
                        "injected descriptor-status interruption")
                return real_fstat(descriptor)

            with mock.patch.object(
                    runner.os, "fstat",
                    side_effect=interrupt_status), \
                 self.assertRaisesRegex(
                     KeyboardInterrupt, "status interruption"):
                lease.__exit__(None, None, None)
            self.assertEqual(lease.descriptor, retained_descriptor)
            self.assertIsNotNone(lease.descriptor_identity)
            self.assertIsNone(lease.kernel_socket)
            lease.__exit__(None, None, None)
            with runner.PairLease(1, 0, root=root):
                pass

    def test_pair_lease_interoperates_with_jerasure_runner(self) -> None:
        jerasure_path = MODULE_PATH.resolve().parents[3] / \
            "tools/leopard2_jerasure_compare.py"
        specification = importlib.util.spec_from_file_location(
            "jerasure_compare_pair_lease_test", jerasure_path)
        self.assertIsNotNone(specification)
        self.assertIsNotNone(specification.loader)
        jerasure = importlib.util.module_from_spec(specification)
        sys.modules[specification.name] = jerasure
        specification.loader.exec_module(jerasure)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lease_directory = root / jerasure.pair_lease_directory_name()
            with runner.PairLease(0, 1, root=lease_directory):
                with self.assertRaises(jerasure.ComparisonError):
                    with jerasure.PairLease(1, 0, root=root):
                        pass
            with jerasure.PairLease(0, 1, root=root):
                with self.assertRaises(runner.EvidenceError):
                    with runner.PairLease(1, 0, root=lease_directory):
                        pass

    def test_pair_lease_creation_ignores_restrictive_umask(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / runner.pair_lease_directory().name
            previous = os.umask(0o777)
            try:
                with runner.PairLease(0, 1, root=root) as identity:
                    self.assertEqual(
                        stat.S_IMODE(os.lstat(Path(identity["path"])).st_mode),
                        0o600)
            finally:
                os.umask(previous)
    def test_stable_anchor_blocks_recreated_reservation_inode(self) -> None:
        payload = {
            "benchmark_cpu": 0,
            "nonce": "stable-anchor-fixture",
            "owner": "unit test",
            "reserved_sibling": 1,
            "schema": runner.RESERVATION_SCHEMA,
            "status": "held",
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            runtime_root = root / "runtime"
            runtime_root.mkdir(mode=0o700)
            runtime_root.chmod(0o700)
            path = root / "reservation.json"
            old = root / "reservation.old"
            path.write_bytes(runner.canonical_bytes(payload))
            with runner.Reservation(
                    path, 0, 1, runtime_root=runtime_root):
                path.rename(old)
                path.write_bytes(runner.canonical_bytes(payload))
                try:
                    with self.assertRaises(runner.EvidenceError):
                        with runner.Reservation(
                                path, 0, 1, runtime_root=runtime_root):
                            pass
                finally:
                    path.unlink()
                    old.rename(path)

    def test_custom_cell_parser_rejects_non_wire_cases(self) -> None:
        good = runner.parse_cell("cell:240:16:65536:8:1")
        self.assertEqual((good.k, good.r, good.losses), (240, 16, 8))
        for bad in (
            "bad:8:9:64:1:1",
            "bad:8:4:65:1:1",
            "bad:8:4:64:5:1",
            "missing-fields",
        ):
            with self.subTest(cell=bad), self.assertRaises(runner.EvidenceError):
                runner.parse_cell(bad)


if __name__ == "__main__":
    unittest.main(verbosity=2)
