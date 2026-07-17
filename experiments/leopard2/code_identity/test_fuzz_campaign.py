#!/usr/bin/env python3
"""Fail-closed tests for the Experiment W campaign runner."""

from __future__ import annotations

import json
import os
import pathlib
import stat
import sys
import tempfile
import time
import unittest
from types import SimpleNamespace
from unittest import mock

sys.dont_write_bytecode = True

import make_fuzz_corpus
import run_fuzz_campaign as campaign


FAKE_FUZZER = r'''#!/usr/bin/env python3
import pathlib
import sys
import time

values = {}
corpus = pathlib.Path(sys.argv[1])
for argument in sys.argv[2:]:
    if argument.startswith("-") and "=" in argument:
        key, value = argument[1:].split("=", 1)
        values[key] = value
mode = pathlib.Path(__file__).read_text().splitlines()[1].strip()
if mode == "# timeout":
    time.sleep(30)
    raise SystemExit(0)
if mode == "# memory":
    allocation = bytearray(64 * 1024 * 1024)
    allocation[0] = 1
    time.sleep(30)
    raise SystemExit(0)
seed = int(values["seed"])
runs = int(values["runs"])
if mode != "# malformed":
    print(f"INFO: Seed: {seed}", file=sys.stderr, flush=True)
    print("#2 cov: 3 ft: 4 corp: 2/4b", file=sys.stderr, flush=True)
    print(f"stat::number_of_executed_units: {runs}", file=sys.stderr, flush=True)
    print("stat::peak_rss_mb: 7", file=sys.stderr, flush=True)
    (corpus / f"seed-{seed}").write_bytes(seed.to_bytes(4, "little"))
if mode == "# sanitizer":
    print("ERROR: AddressSanitizer: fake finding", file=sys.stderr, flush=True)
'''


class CampaignTests(unittest.TestCase):
    def test_distinct_stable_seeds(self) -> None:
        expected = campaign.derive_seeds(1279609156, 28)
        self.assertEqual(expected, campaign.derive_seeds(1279609156, 28))
        self.assertEqual(28, len(set(expected)))
        self.assertNotIn(0, expected)

    def test_output_directories_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            fresh = root / "fresh"
            campaign.require_empty_output(fresh)
            self.assertTrue(fresh.is_dir())
            (fresh / "stale").write_bytes(b"x")
            with self.assertRaisesRegex(ValueError, "not empty"):
                campaign.require_empty_output(fresh)

            corpus = root / "corpus"
            corpus.mkdir()
            (corpus / "stale").write_bytes(b"x")
            with self.assertRaisesRegex(ValueError, "not empty"):
                make_fuzz_corpus.require_clean_output(corpus)

            target = root / "target"
            target.mkdir()
            link = root / "link"
            link.symlink_to(target, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "real directory"):
                campaign.require_empty_output(link)
            with self.assertRaisesRegex(ValueError, "real directory"):
                make_fuzz_corpus.require_clean_output(link)

    def test_content_addressed_merge_is_order_independent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            workers = []
            for index, values in enumerate(((b"a", b"shared"),
                                            (b"b", b"shared"))):
                worker_root = root / f"worker-{index}"
                corpus = worker_root / "corpus"
                corpus.mkdir(parents=True)
                for item_index, value in enumerate(values):
                    (corpus / str(item_index)).write_bytes(value)
                workers.append(SimpleNamespace(root=worker_root))
            first = campaign.merge_corpora(workers, root / "merged-first")
            second = campaign.merge_corpora(
                list(reversed(workers)), root / "merged-second"
            )
            self.assertEqual(first, second)
            self.assertEqual(3, first["files"])

    def _make_fuzzer(self, root: pathlib.Path, mode: str) -> pathlib.Path:
        fuzzer = root / "fake-fuzzer"
        source = FAKE_FUZZER.replace(
            "#!/usr/bin/env python3\n", f"#!/usr/bin/env python3\n# {mode}\n", 1
        )
        fuzzer.write_text(source, encoding="utf-8")
        fuzzer.chmod(fuzzer.stat().st_mode | stat.S_IXUSR)
        return fuzzer

    def _run(self, root: pathlib.Path, mode: str, results_name: str,
             timeout: float = 5.0,
             external_rss_limit_mb: int = 10240) -> pathlib.Path:
        fuzzer = self._make_fuzzer(root, mode)
        corpus = root / "initial-corpus"
        if not corpus.exists():
            corpus.mkdir()
            (corpus / "seed").write_bytes(b"seed")
        results = root / results_name
        arguments = [
            "run_fuzz_campaign.py", "--fuzzer", str(fuzzer),
            "--corpus", str(corpus), "--results", str(results),
            "--workers", "2", "--runs", "7",
            "--timeout-seconds", str(timeout),
            "--external-rss-limit-mb", str(external_rss_limit_mb),
        ]
        with mock.patch.object(sys, "argv", arguments), mock.patch.object(
            campaign, "git_identity", return_value=("f" * 40, False)
        ):
            campaign.main()
        return results

    def test_isolated_workers_and_deterministic_campaign(self) -> None:
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("two allowed CPUs required")
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            first = self._run(root, "success", "results-first")
            second = self._run(root, "success", "results-second")
            one = json.loads((first / "manifest.json").read_text())
            two = json.loads((second / "manifest.json").read_text())
            self.assertEqual(2, one["aggregate"]["distinct_seeds"])
            self.assertEqual(14, one["aggregate"]["total_runs"])
            self.assertEqual(
                one["aggregate"]["merged_corpus"],
                two["aggregate"]["merged_corpus"],
            )
            self.assertNotEqual(
                one["workers"][0]["corpus"]["content_sha256"],
                one["workers"][1]["corpus"]["content_sha256"],
            )
            self.assertEqual(
                one["workers"][0]["seed_input_corpus"],
                one["workers"][1]["seed_input_corpus"],
            )
            self.assertEqual(
                one["workers"][0]["seed_input_list_sha256"],
                one["workers"][1]["seed_input_list_sha256"],
            )
            for worker_index in range(2):
                worker = first / f"worker-{worker_index:03d}"
                self.assertEqual(1, len(list((worker / "seed-inputs").iterdir())))
                self.assertFalse((worker / "corpus" / "seed").exists())

    def test_malformed_log_and_sanitizer_marker_fail(self) -> None:
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("two allowed CPUs required")
        for mode in ("malformed", "sanitizer"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as tmp:
                with self.assertRaises(ValueError):
                    self._run(pathlib.Path(tmp), mode, "results")

    def test_timeout_is_bounded(self) -> None:
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("two allowed CPUs required")
        with tempfile.TemporaryDirectory() as temporary:
            started = time.monotonic()
            with self.assertRaises(ValueError):
                self._run(
                    pathlib.Path(temporary), "timeout", "results", timeout=0.05
                )
            self.assertLess(time.monotonic() - started, 2.0)

    def test_external_rss_cap_is_bounded(self) -> None:
        if len(os.sched_getaffinity(0)) < 2:
            self.skipTest("two allowed CPUs required")
        with tempfile.TemporaryDirectory() as temporary:
            started = time.monotonic()
            with self.assertRaises(ValueError):
                self._run(
                    pathlib.Path(temporary), "memory", "results",
                    external_rss_limit_mb=16,
                )
            self.assertLess(time.monotonic() - started, 2.0)


if __name__ == "__main__":
    unittest.main()
