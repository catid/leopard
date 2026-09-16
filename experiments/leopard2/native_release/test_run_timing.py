#!/usr/bin/env python3
"""Full collector/replay simulation and evidence mutations; no codec executions."""
import contextlib
import copy
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from check_encode import sha256, BASELINE
from frozen_inputs import FrozenInputs
import run_timing as runner
import timing_contract as contract
from test_timing_contract import fixture_expected, timed_record


class CollectorTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.expected = fixture_expected()
        self.bundle = self.root / ".research/leopard-79h/campaign/artifacts"
        self.bundle.mkdir(parents=True)
        for name in runner.RUN_ARTIFACTS:
            value = json.dumps(self.expected) if name == "expected.json" else name
            (self.bundle / name).write_text(value)
            (self.bundle / name).chmod(0o444)
        manifest = self.bundle / "SHA256SUMS"
        manifest.write_text("".join(f"{sha256(self.bundle / name)}  {name}\n"
                                    for name in sorted(runner.RUN_ARTIFACTS)))
        manifest.chmod(0o444)
        frozen = FrozenInputs(self.bundle, runner.RUN_ARTIFACTS)
        self.plan = {
            "schema": "native-release-encode-plan/v1", "candidate_commit": contract.CANDIDATE,
            "baseline_commit": BASELINE, "cpu": 52, "sibling": 116, "controller_cpu": 0,
            "passive_seconds": 10, "rounds": 3, "samples": 9, "attempt_budget": 1,
            "groups": list(contract.GROUPS), "minimum_group_ns": 20000000,
            "control_factor": 1.02, "warmup_calls": 4, "initial_checks": 2,
            "child_cpu_seconds": 180, "child_timeout_seconds": 180,
            "child_address_space_bytes": 268435456, "child_file_bytes": 1048576,
            "promotion": False, "cells": [list(cell) for cell in contract.CELLS],
            "comparisons": [[name, list(order)] for name, order in contract.COMPARISONS],
            "source_sha256": {name: sha256(Path(runner.__file__).parent / name)
                              for name in runner.SOURCE_NAMES},
            "runtime_sha256": {name: sha256(Path(name)) for name in runner.RUNTIME_NAMES},
            "artifact_sha256": frozen.pins, "artifact_manifest_sha256": frozen.manifest_sha256,
            "artifact_directory": str(self.bundle.relative_to(self.root)),
            "attempt_directory": ".research/leopard-79h/campaign/attempt-1", "host": runner.host_identity()}
        self.path = self.root / "plan.json"
        self.path.write_text(json.dumps(self.plan))
        self.output = self.root / self.plan["attempt_directory"]
        self.resources = {"cgroup": "/fixture", "memory.max": 268435456, "memory.peak": 123456,
                          "memory.swap.max": 0, "memory.swap.current": 0,
                          "events": dict.fromkeys(("low", "high", "max", "oom", "oom_kill", "oom_group_kill"), 0)}

    def simulate(self, fail=False):
        ticks = [0]
        def counter(cpu):
            if cpu == 52: ticks[0] += 1
            return ticks[0] if cpu == 52 else 0
        def child(command, **kwargs):
            implementation = Path(command[-3]).name.split("-")[0]
            cell = int(command[-1])
            measured = command[-2] == "--measure"
            if measured and fail:
                kwargs["stderr"].write(b"injected failure\n")
                return subprocess.CompletedProcess(command, 1)
            value = timed_record(self.expected, cell, implementation) if measured else self.expected[cell][implementation]
            kwargs["stdout"].write(json.dumps(value).encode())
            return subprocess.CompletedProcess(command, 0)
        with contextlib.ExitStack() as stack:
            for target, options in (
                ("load_published_plan", {"return_value": (self.root, self.plan)}),
                ("resource_snapshot", {"return_value": self.resources}),
                ("cpu_ticks", {"side_effect": counter}),
                ("os.sched_getaffinity", {"return_value": {0, 52, 116}}),
                ("os.sched_setaffinity", {}), ("fcntl.flock", {}),
                ("subprocess.run", {"side_effect": child}), ("time.sleep", {}),
                ("time.monotonic_ns", {"side_effect": (0, 10000000000)})):
                stack.enter_context(patch("run_timing." + target, **options))
            stack.enter_context(contextlib.redirect_stdout(io.StringIO()))
            runner.run(self.path, "a" * 40)

    def replay(self):
        with contextlib.redirect_stdout(io.StringIO()):
            return runner.replay(self.path, self.output, self.bundle)

    def test_full_runtime_and_replay(self):
        runner.validate_plan(self.plan)
        self.simulate()
        self.assertEqual(self.replay()["decision"], "complete")
        with self.assertRaises(FileExistsError): self.simulate()

    def test_error_retained_and_unreplayable(self):
        with self.assertRaises(ValueError): self.simulate(fail=True)
        state = json.loads((self.output / "attempt.json").read_text())
        self.assertFalse(state["complete"])
        self.assertIsNone(state["analysis"])
        self.assertEqual(len(state["preflight"]), 16)
        self.assertIn("benchmark child failed", state["failure"]["message"])
        self.assertTrue((self.output / "cell-0-round-0-native_vs_current-slot-0.stderr").exists())
        with self.assertRaises(ValueError): self.replay()

    def test_replay_rejects_observation_reuse_missing_resources_and_drift(self):
        self.simulate()
        path = self.output / "attempt.json"
        original = json.loads(path.read_text())
        for mutation in ("reuse", "events", "command", "source", "partial"):
            state = copy.deepcopy(original)
            changed_path = None
            if mutation == "reuse":
                state["invocations"][3]["metadata"] = state["invocations"][0]["metadata"]
            elif mutation == "events": state["resources_after"]["events"] = {}
            elif mutation == "source": state["preflight"][0]["record"]["codec_commit"] = "b" * 40
            elif mutation == "partial": state["invocations"].pop()
            else:
                changed_path = self.output / state["invocations"][0]["metadata"]
                old_bytes = changed_path.read_bytes()
                value = json.loads(old_bytes)
                value["command"][2] = "53"
                changed_path.write_text(json.dumps(value))
            path.write_text(json.dumps(state))
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): self.replay()
            if changed_path: changed_path.write_bytes(old_bytes)
        path.write_text(json.dumps(original))
        self.assertEqual(self.replay()["decision"], "complete")

    def test_plan_mutations(self):
        for key, value in (("groups", [1] * 8), ("runtime_sha256", {}),
                           ("artifact_sha256", {}), ("source_sha256", {}),
                           ("attempt_directory", self.plan["artifact_directory"] + "/attempt"),
                           ("artifact_manifest_sha256", "bad"), ("attempt_budget", True)):
            bad = copy.deepcopy(self.plan)
            bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): runner.validate_plan(bad)

    def test_repository_preregistration_sources(self):
        directory = Path(runner.__file__).parent
        plan = json.loads((directory / "timing_plan_v1.json").read_text())
        runner.validate_plan(plan)
        for name, digest in plan["source_sha256"].items():
            self.assertEqual(sha256(directory / name), digest)


if __name__ == "__main__":
    unittest.main()
