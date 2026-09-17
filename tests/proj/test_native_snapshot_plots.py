#!/usr/bin/env python3
"""Check displayed geometry independently against the retained measurements."""

import copy
import importlib.util
import json
from pathlib import Path
import re
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
        flat = " ".join(readme.split())
        self.assertIn("AVX2-restricted Leopard1", flat)
        self.assertIn("not final-release evidence", flat)
        self.assertIn("Shipped in `master`", flat)
        self.assertIn("retained research records", flat)
        self.assertNotIn("1.52×", flat)
        self.assertNotIn("cross-process stability control was inconclusive", flat)
        atlas = ROOT / "docs/performance/leopard2_atlas"
        metadata = json.loads((atlas / "run_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual("ON", metadata["build_closures"]["leopard1"]["cache_contract"]["LEO_MAIN_PURE_AVX2"])
        for path in (atlas / "plots").glob("*speedup_vs_leopard1.svg"):
            self.assertIn("AVX2-restricted", path.read_text(encoding="utf-8"))

    def test_readme_embeds_multi_size_leopard_comparisons(self):
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        expected = (
            "docs/performance/leopard2_atlas/plots/encode_speedup_vs_leopard1.svg",
            "docs/performance/leopard2_atlas/plots/decode_one_speedup_vs_leopard1.svg",
            "docs/performance/leopard2_atlas/plots/decode_full_speedup_vs_leopard1.svg",
        )
        for relative in expected:
            self.assertIn(f"]({relative})", readme)
            svg = (ROOT / relative).read_text(encoding="utf-8")
            for size in ("64 B", "1 KiB", "4 KiB", "1 MiB"):
                self.assertIn(f">{size}</text>", svg)
        self.assertIn("not native-Leopard1 release claims", readme)


class CurrentNativeTimingPlotTest(unittest.TestCase):
    def setUp(self):
        report = ROOT / "docs/performance/native_release_encode_timing_v1.md"
        svg = report.with_suffix(".svg")
        self.report = report.read_text(encoding="utf-8")
        self.svg = svg.read_text(encoding="utf-8")
        self.root = ET.fromstring(self.svg)
        self.elements = {e.attrib["id"]: e for e in self.root.iter() if "id" in e.attrib}
        self.ratios = dict((name, float(value)) for name, value in re.findall(
            r"^\| ([a-z0-9-]+) \| ([0-9.]+) \|", self.report, re.M))

    def check_geometry(self):
        expected_names = {"copy", "small", "gf8-high", "gf8-balanced", "gf16-inflation",
                          "gf16-gfni-region", "gf16-explicit-avx2", "gf16-large"}
        self.assertEqual(set(self.ratios), expected_names)
        self.assertEqual({name[4:] for name in self.elements if name.startswith("bar-")}, expected_names)
        baseline = self.elements["native-current-baseline"].attrib
        self.assertEqual(float(baseline["y1"]), 276)
        self.assertEqual(float(baseline["y2"]), 276)
        for tick in (0.8, 1.0, 1.2, 1.4, 1.5):
            line = self.elements[f"tick-{tick:.1f}"].attrib
            self.assertAlmostEqual(float(line["y1"]), 276 - (tick-1)*370)
            self.assertEqual(line["y1"], line["y2"])
        for name, ratio in self.ratios.items():
            bar = self.elements[f"bar-{name}"].attrib
            top, height = float(bar["y"]), float(bar["height"])
            expected_y = 276 - (ratio-1)*370
            self.assertAlmostEqual(top, min(276, expected_y), delta=0.001)
            self.assertAlmostEqual(height, abs(expected_y-276), delta=0.001)
            self.assertGreaterEqual(top, 91)
            self.assertLessEqual(top+height, 350)
            self.assertEqual(bar["fill"], "#777")  # Inconclusive, not win/loss colors.
            value = self.elements[f"value-{name}"]
            self.assertEqual(value.text, f"{ratio:.3f}")
            # Keep values clear of title/subtitle and workload labels.
            self.assertGreaterEqual(float(value.attrib["y"]), 80)
            self.assertLessEqual(float(value.attrib["y"]), 350)

    def test_geometry_and_labels_match_all_eight_report_rows(self):
        self.check_geometry()

    def test_old_axis_and_title_overlap_regressions_are_detected(self):
        for identifier, field, value in (("tick-1.4", "y1", "86"),
                                          ("bar-small", "height", "205"),
                                          ("value-small", "y", "46")):
            element = self.elements[identifier]
            previous = element.attrib[field]
            element.set(field, value)
            with self.assertRaises(AssertionError):
                self.check_geometry()
            element.set(field, previous)

    def test_provenance_and_inconclusive_boundary_remain_explicit(self):
        self.assertIn("preregistration/attempt commit", self.report)
        self.assertIn("`15756b1`", self.report)
        self.assertIn("`e35b1f0c3a5ef9f241ad3b174b0fbfdb61bbccab`", self.report)
        self.assertIn("inconclusive_controls", self.report)
        self.assertIn("no qualified win/loss claim", self.svg)
        self.assertNotIn("1/1.02 gate", self.svg)  # This gate applies to controls, not these ratios.


if __name__ == "__main__":
    unittest.main()
