#!/usr/bin/env python3
"""Evaluate run_sweep_jureca*.sh output: cost vs. ND (channel count) and NG (gas count),
per geometry (limb/nadir/zenith)."""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_series
from plot_style import PALETTE
from likwid_parsing import parse_run_dir, collect_runtime, collect

GEOMETRY_COLORS = {"limb": PALETTE[0], "nadir": PALETTE[1], "zenith": PALETTE[2]}

DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

ND_RE = re.compile(r"^\s*ND\s*=\s*(\d+)", re.MULTILINE)
NG_RE = re.compile(r"^\s*NG\s*=\s*(\d+)", re.MULTILINE)

# Sweep labels: channels_<nd>_<geom> / gases_<set>_<geom>. The geometry suffix is
# optional so run directories from before the geometry axis still evaluate.
# Gas-set names may contain underscores, hence the anchored match.
LABEL_RE = re.compile(r"^(channels|gases)_(.+?)(?:_(limb|nadir|zenith))?$")


def parse_label(label: str):
    """Return (axis, geometry) for a sweep label; axis is None for other labels."""
    m = LABEL_RE.match(label)
    if not m:
        return None, None
    return m.group(1), m.group(3)


def read_nd_ng(ctl_path: Path):
    if not ctl_path.exists():
        return None, None
    text = ctl_path.read_text()
    nd = ND_RE.search(text)
    ng = NG_RE.search(text)
    return (int(nd.group(1)) if nd else None), (int(ng.group(1)) if ng else None)


def is_bandwidth_like(metric: str) -> bool:
    m = metric.lower()
    return "bandwidth" in m or "mflop/s" in m


def safe_name(metric: str) -> str:
    return metric.split("[")[0].strip().replace(" ", "_").lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--region", default="formod", help="likwid-marker region to read (default: formod)")
    parser.add_argument("--metrics", action="append",
                         help="metric names to evaluate; default: %s" % ", ".join(DEFAULT_METRICS))
    parser.add_argument("--stream-bw", type=float, default=None,
                         help="measured STREAM bandwidth ceiling [MBytes/s] for bandwidth-like metrics")
    parser.add_argument("--out", type=Path, default=None, help="output directory (default: <run_dir>/plots)")
    args = parser.parse_args()

    metrics = args.metrics or DEFAULT_METRICS
    res_dir = args.out or (args.run_dir / "plots")
    res_dir.mkdir(parents=True, exist_ok=True)
    ctl_dir = args.run_dir / "ctl"

    configs = parse_run_dir(args.run_dir)
    if not configs:
        print(f"No parsed configs found under {args.run_dir}.")
        sys.exit(1)

    groups: dict[tuple, list] = {}
    for c in configs:
        key = (c["label"], c["threads"], c["group"], c["batch_size"])
        groups.setdefault(key, []).append(c)

    metric_keys = ["wall_time"] + metrics
    medians: dict[str, dict[str, float]] = {m: {} for m in metric_keys}
    nd_by_label: dict[str, int] = {}
    ng_by_label: dict[str, int] = {}

    max_cv, max_cv_desc = 0.0, None

    header = (f"{'label':>22} {'nd':>5} {'ng':>5} {'thr':>4} {'group':>8} {'batch':>6} | " +
              " | ".join(f"{m[:22]:>22}" for m in ["runtime/call [s]"] + metrics))
    print(header)
    print("-" * len(header))

    for (label, threads, group, batch), entries in sorted(groups.items()):
        nd, ng = read_nd_ng(ctl_dir / f"{label}.ctl")
        nd_by_label[label] = nd
        ng_by_label[label] = ng

        kept = sorted(entries, key=lambda e: e["rep"])

        row_values = []
        mean, median, stdev, cv, _, _ = collect_runtime(kept)
        row_values.append(
            f"{f'{mean:.4g}' if mean is not None else 'N/A'} "
            f"(sd={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
            f"cv={f'{cv:.1%}' if cv is not None and not math.isnan(cv) else 'N/A'})"
        )
        medians["wall_time"][label] = median
        if cv is not None and not math.isnan(cv) and cv > max_cv:
            max_cv, max_cv_desc = cv, (label, "wall_time")

        for metric in metrics:
            mean, median, stdev, cv, _, _ = collect(kept, args.region, metric)
            row_values.append(
                f"{f'{mean:.4g}' if mean is not None else 'N/A'} "
                f"(sd={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
                f"cv={f'{cv:.1%}' if cv is not None and not math.isnan(cv) else 'N/A'})"
            )
            medians[metric][label] = median
            if cv is not None and not math.isnan(cv) and cv > max_cv:
                max_cv, max_cv_desc = cv, (label, metric)

        print(f"{label:>22} {nd!s:>5} {ng!s:>5} {threads:>4} {group:>8} {batch:>6} | " +
              " | ".join(f"{v:>22}" for v in row_values))

    if max_cv_desc:
        print(f"\nLargest observed coefficient of variation: {max_cv:.1%} "
              f"(label={max_cv_desc[0]}, metric='{max_cv_desc[1]}')")
    else:
        print("\nNo metric produced a usable coefficient of variation.")

    def plot_axis(axis_vals_by_label: dict, axis: str, xlabel: str):
        # One plot per metric with one curve per geometry (legacy unsuffixed
        # labels, geom=None, form a single curve named "all")
        geoms = sorted({parse_label(l)[1] for l, v in axis_vals_by_label.items()
                        if v is not None and parse_label(l)[0] == axis},
                       key=lambda g: (g is None, g or ""))
        if not geoms:
            return

        series_by_geom = {}
        for geom in geoms:
            labels = sorted(
                (l for l, v in axis_vals_by_label.items()
                 if v is not None and parse_label(l) == (axis, geom)),
                key=lambda l: axis_vals_by_label[l],
            )
            if len(labels) < 2:
                print(f"\nFewer than 2 '{axis}' points for geometry={geom or 'n/a'}; skipping that curve.")
                continue
            series_by_geom[geom] = labels

        if not series_by_geom:
            print(f"\nFewer than 2 '{axis}' points with a resolved axis value. Skipping plots.")
            return

        print(f"\n=== {xlabel} scaling ({axis}) ===")

        for metric in ["wall_time"] + metrics:
            series = []
            for geom, labels in series_by_geom.items():
                metric_medians = [medians[metric].get(l) for l in labels]
                if any(v is None for v in metric_medians):
                    print(f"WARNING: missing '{metric}' for some {axis} points (geometry={geom or 'n/a'}); "
                          "leaving that curve out.")
                    continue
                x = [axis_vals_by_label[l] for l in labels]
                series.append((geom or "all", x, metric_medians, GEOMETRY_COLORS.get(geom, PALETTE[0])))
            if not series:
                continue

            bw_like = metric != "wall_time" and is_bandwidth_like(metric)
            ylabel = "Wall-clock time [s] per call" if metric == "wall_time" else metric
            metric_part = "wallclock" if metric == "wall_time" else safe_name(metric)
            fname = f"e3_{axis}_{metric_part}.png"
            stream_ceiling = args.stream_bw if bw_like else None

            # No "ideal" reference line: channel/gas count isn't a uniform
            # work unit here (gases differ widely in table size/line count,
            # channels in which tables they hit), so there's no growth curve
            # to assert as ideal -- only the measured values are plotted.
            plot_series(
                series, f"{ylabel} ({axis})", res_dir, fname,
                stream_ceiling=stream_ceiling, xlabel=xlabel,
                yscale="linear" if bw_like else "log",
            )

    plot_axis(nd_by_label, "channels", "Channels (ND)")
    plot_axis(ng_by_label, "gases", "Emitters (NG)")

    print(f"\nPlots written to {res_dir}/")


if __name__ == "__main__":
    main()
