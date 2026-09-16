#!/usr/bin/env python3
"""Render the retained native snapshot, not a claim about later releases."""

import argparse
import html
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "docs/performance/final_native_gfni_summary.json"
PLOTS = ROOT / "docs/performance/leopard2_atlas/plots"


def number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("plot value must be numeric")
    if not math.isfinite(value) or value < 0:
        raise ValueError("plot value must be finite and nonnegative")
    return value


def text(x, y, label, size=12, anchor="middle"):
    return (f'<text x="{x}" y="{y:.3f}" text-anchor="{anchor}" '
            f'font-family="sans-serif" font-size="{size}">'
            f'{html.escape(str(label))}</text>')


def line(x1, y1, x2, y2, identifier="", color="#222"):
    tag = f' id="{identifier}"' if identifier else ""
    return (f'<line{tag} x1="{x1}" y1="{y1:.3f}" x2="{x2}" '
            f'y2="{y2:.3f}" stroke="{color}"/>')


def bar(identifier, x, value, maximum, bottom, height, color, width=45):
    value = number(value)
    if value > maximum:
        raise ValueError("plot value exceeds axis")
    scaled = height * value / maximum
    return (f'<rect id="{identifier}" x="{x}" y="{bottom-scaled:.3f}" '
            f'width="{width}" height="{scaled:.3f}" fill="{color}"/>')


def canvas(width, height, title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}">',
            f'<rect width="{width}" height="{height}" fill="white"/>',
            text(width / 2, 30, title, 19)]


def snapshot_label(data):
    return (f'Measured L2 {data["source_commit"][:12]}; '
            f'native L1 {data["baseline_commit"][:12]}; not a later-release claim')


def render(data):
    if (data["schema"] != "leopard2-final-source-native-leopard1-gfni-v20/v1"
            or data["status"] != "valid"):
        raise ValueError("expected validated native snapshot summary")
    ratio = number(data["encode"]["geometric_speedup"])
    low, high = map(number, data["encode"]["ci95"])
    if not 0 < low <= ratio <= high <= 2:
        raise ValueError("invalid confidence interval or axis overflow")
    speed = canvas(760, 450, "GF16 encode vs native Leopard1 — measured snapshot")
    speed.append(text(380, 52, "K=1000, R=200, 64 KiB; full parity; one thread"))
    for tick in (0, 1, 2):
        y = 350 - tick * 140
        speed += [line(90, y, 700, y, color="#ddd"),
                  text(78, y + 5, f"{tick}×", anchor="end")]
    speed += [line(90, 70, 90, 350),
              bar("encode-speedup", 290, ratio, 2, 350, 280, "#2b6cb0", 180),
              line(380, 350-high*140, 380, 350-low*140, "encode-ci"),
              line(365, 350-high*140, 395, 350-high*140),
              line(365, 350-low*140, 395, 350-low*140),
              text(380, 350-high*140-12, f"{ratio:.3f}×", 18),
              text(380, 380, f"95% CI {low:.3f}–{high:.3f}; "
                   f'{data["rounds"]} ABBA rounds; native L1 time / L2 time'),
              text(380, 410, snapshot_label(data), 11),
              text(380, 432, "AMD Threadripper 9980X; CPU 52, sibling 116 reserved", 11),
              "</svg>"]

    one = data["steady_state_measurements"]["one_loss"]
    full = data["steady_state_measurements"]["full_loss"]
    metrics = canvas(960, 530, "GF16 metrics vs native Leopard1 — measured snapshot")
    metrics += [text(260, 68, "Encode throughput (GB/s)", 15),
                text(730, 68, "Work / scratch allocation (MiB)", 15)]
    for left, right, maximum, ticks in (
            (75, 455, 12, (0, 4, 8, 12)),
            (550, 930, 128, (0, 32, 64, 96, 128))):
        for tick in ticks:
            y = 300 - 210 * tick / maximum
            metrics += [line(left, y, right, y, color="#ddd"),
                        text(left-8, y+4, tick, anchor="end")]
        metrics.append(line(left, 90, left, 300))
    rates = [("l1-one", one["native_leopard1_encode_gb_s"], "L1 one", "#777"),
             ("l2-one", one["leopard2_encode_gb_s"], "L2 one", "#2b6cb0"),
             ("l1-full", full["native_leopard1_encode_gb_s"], "L1 full", "#777"),
             ("l2-full", full["leopard2_encode_gb_s"], "L2 full", "#2b6cb0")]
    for x, (identifier, value, label, color) in zip((105, 185, 305, 385), rates):
        metrics += [bar(identifier, x, value, 12, 300, 210, color),
                    text(x+22.5, 292-value*210/12, f"{value:.2f}"),
                    text(x+22.5, 320, label)]
    memory = [("l1-enc", "native_encode_work_bytes", "L1 enc", "#777"),
              ("l2-enc", "leopard2_encode_scratch_bytes", "L2 enc", "#2b6cb0"),
              ("l1-dec", "native_decode_work_bytes", "L1 dec", "#777"),
              ("l2-dec", "leopard2_decode_scratch_bytes", "L2 dec", "#2b6cb0")]
    for x, (identifier, key, label, color) in zip((580, 660, 780, 860), memory):
        value = number(full[key]) / (1024 * 1024)
        metrics += [bar(identifier, x, value, 128, 300, 210, color),
                    text(x+22.5, 292-value*210/128, f"{value:.2f}"),
                    text(x+22.5, 320, label)]
    codec = [number(item["leopard2_codec_setup_us"]) for item in (one, full)]
    plan = [number(item["leopard2_decode_plan_setup_us"]) for item in (one, full)]
    metrics += [text(480, 360, "K=1000, R=200, 64 KiB; one thread; one/full = decode loss scenario"),
                text(480, 385, "Memory panel: full-loss scenario, common linear scale; not total resident memory"),
                text(480, 415, f"L2 codec setup {min(codec):.3f}–{max(codec):.3f} µs; "
                     f"decode-plan setup {min(plan):.3f}–{max(plan):.3f} µs"),
                text(480, 440, "Standalone medians, not ABBA speedups; gray = native L1, blue = L2"),
                text(480, 475, snapshot_label(data), 11),
                text(480, 500, "Setup values predate the later GF16 Walsh-locator optimization", 11),
                "</svg>"]
    return {"final_native_gfni_encode_speedup.svg": "\n".join(speed) + "\n",
            "final_native_gfni_metrics.svg": "\n".join(metrics) + "\n"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify without writing")
    args = parser.parse_args()
    plots = render(json.loads(SUMMARY.read_text(encoding="utf-8")))
    for name, svg in plots.items():
        path = PLOTS / name
        if args.check:
            if path.read_text(encoding="utf-8") != svg:
                raise SystemExit(f"stale native snapshot plot: {name}")
        else:
            path.write_text(svg, encoding="utf-8")


if __name__ == "__main__":
    main()
