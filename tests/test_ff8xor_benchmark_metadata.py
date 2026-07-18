#!/usr/bin/env python3
"""Check that finite benchmark rows describe their actual sampling order."""

from __future__ import print_function

import json
import subprocess
import sys


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    require(len(sys.argv) == 2, "expected benchmark executable path")
    result = subprocess.run(
        [sys.argv[1], "--quick", "--json", "--no-pin"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    require(result.returncode == 0,
            "quick benchmark failed (%d): %s" %
            (result.returncode, result.stderr))

    records = [json.loads(line) for line in result.stdout.splitlines()
               if line.strip()]
    metadata = [record for record in records
                if record.get("record") == "metadata"]
    require(len(metadata) == 1,
            "benchmark did not emit exactly one metadata record")
    require("XOR-batch pairs ABBA" in metadata[0]["measurement_order"],
            "metadata omitted unconditional XOR-batch ABBA sampling")

    expected_operations = {
        "xor2_sequential", "xor2_batched",
        "xor4_sequential", "xor4_batched",
    }
    rows = [record for record in records
            if record.get("record") == "microbenchmark" and
            record.get("operation") in expected_operations]
    observed_operations = {record["operation"] for record in rows}
    require(observed_operations == expected_operations,
            "quick benchmark XOR-batch row set changed: %r" %
            sorted(observed_operations))
    require(len(rows) == 4 * len(expected_operations),
            "quick benchmark did not cover four XOR-batch sizes")
    for record in rows:
        require(record.get("measurement_order") == "ABBA",
                "%s mislabeled paired sampling as %r" %
                (record["operation"], record.get("measurement_order")))
        require("paired ABBA" in record.get("note", ""),
                "%s row lost paired-sampling note" % record["operation"])

    print("FF8 XOR benchmark metadata tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
