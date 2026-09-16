#!/usr/bin/env python3
"""Check displayed geometry independently against the retained measurements."""

import copy
import importlib.util
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "native_plots", ROOT / "tools/leopard2_native_snapshot_plots.py")
PLOTS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOTS)


class NativeSnapshotPlotsTest(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(PLOTS.SUMMARY.read_text(encoding="utf-8"))
        self.rendered = PLOTS.render(self.data)

    def elements(self, name):
        root = ET.fromstring(self.rendered[name])
        return {item.attrib["id"]: item for item in root.iter() if "id" in item.attrib}

    def test_checked_in_plots_are_reproducible(self):
        for name, expected in self.rendered.items():
            self.assertEqual(expected, (PLOTS.PLOTS / name).read_text(encoding="utf-8"))

    def test_speedup_and_both_ci_endpoints_share_axis(self):
        elements = self.elements("final_native_gfni_encode_speedup.svg")
        bar = elements["encode-speedup"].attrib
        ci = elements["encode-ci"].attrib
        expected = self.data["encode"]
        self.assertAlmostEqual(float(bar["height"])/140, expected["geometric_speedup"], places=5)
        self.assertAlmostEqual((350-float(bar["y"]))/140, expected["geometric_speedup"], places=5)
        self.assertAlmostEqual((350-float(ci["y1"]))/140, expected["ci95"][1], places=5)
        self.assertAlmostEqual((350-float(ci["y2"]))/140, expected["ci95"][0], places=5)

    def test_throughput_uses_one_scale_for_all_four_bars(self):
        elements = self.elements("final_native_gfni_metrics.svg")
        for scenario, suffix in (("one_loss", "one"), ("full_loss", "full")):
            for prefix, key in (("l1", "native_leopard1_encode_gb_s"), ("l2", "leopard2_encode_gb_s")):
                expected = self.data["steady_state_measurements"][scenario][key]
                bar = elements[f"{prefix}-{suffix}"].attrib
                self.assertAlmostEqual(float(bar["height"])*12/210, expected, delta=0.00004)
                self.assertAlmostEqual(float(bar["y"])+float(bar["height"]), 300)

    def test_memory_uses_common_mib_scale_and_full_loss_values(self):
        elements = self.elements("final_native_gfni_metrics.svg")
        for identifier, key in (("l1-enc", "native_encode_work_bytes"),
                                ("l2-enc", "leopard2_encode_scratch_bytes"),
                                ("l1-dec", "native_decode_work_bytes"),
                                ("l2-dec", "leopard2_decode_scratch_bytes")):
            expected = self.data["steady_state_measurements"]["full_loss"][key] / 1048576
            bar = elements[identifier].attrib
            self.assertAlmostEqual(float(bar["height"])*128/210, expected, delta=0.0004)
            self.assertAlmostEqual(float(bar["y"])+float(bar["height"]), 300)

    def test_snapshot_identity_and_limitations_are_visible(self):
        for svg in self.rendered.values():
            self.assertIn(self.data["source_commit"][:12], svg)
            self.assertIn(self.data["baseline_commit"][:12], svg)
            self.assertIn("not a later-release claim", svg)
            self.assertNotIn("Final-source", svg)
        self.assertIn("predate", self.rendered["final_native_gfni_metrics.svg"])

    def test_bad_intervals_and_nonfinite_values_fail_closed(self):
        for interval in ([1.6, 1.3], [1.0, 2.1], [float("nan"), 1.5]):
            data = copy.deepcopy(self.data)
            data["encode"]["ci95"] = interval
            with self.assertRaises(ValueError):
                PLOTS.render(data)
        data = copy.deepcopy(self.data)
        data["steady_state_measurements"]["one_loss"]["leopard2_encode_gb_s"] = float("inf")
        with self.assertRaises(ValueError):
            PLOTS.render(data)

    def test_native_and_restricted_claims_are_distinct(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertIn("AVX2-restricted Leopard1", readme)
        self.assertIn("not final-release evidence", readme)
        self.assertIn("representative remaining losses", readme)
        atlas = ROOT / "docs/performance/leopard2_atlas"
        metadata = json.loads((atlas / "run_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual("ON", metadata["build_closures"]["leopard1"]["cache_contract"]["LEO_MAIN_PURE_AVX2"])
        for path in (atlas / "plots").glob("*speedup_vs_leopard1.svg"):
            self.assertIn("AVX2-restricted", path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
