#!/usr/bin/env python3
"""Evaluate run_sweep_jureca*.sh output: cost vs. ND (channel count) and NG (gas count)."""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_scaling
from likwid_parsing import parse_run_dir, collect_runtime, collect

DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

ND_RE = re.compile(r"^\s*ND\s*=\s*(\d+)", re.MULTILINE)
NG_RE = re.compile(r"^\s*NG\s*=\s*(\d+)", re.MULTILINE)


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
    parser.add_argument("--warmup", type=int, default=0,
                         help="leading repetitions to discard (default: 0; run_sweep_jureca*.sh "
                              "calls bench_run_forward once per point, not in a rep loop, so "
                              "there's normally exactly one sample per point already)")
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

        entries = sorted(entries, key=lambda e: e["rep"])
        kept = entries[args.warmup:]
        if not kept:
            print(f"{label:>22} {nd!s:>5} {ng!s:>5} {threads:>4} {group:>8} {batch:>6} | "
                  f"no data left after discarding {args.warmup} warmup repetition(s)")
            continue

        row_values = []
        mean, median, stdev, cv, _, _ = collect_runtime(kept, warmup=0)
        row_values.append(
            f"{f'{mean:.4g}' if mean is not None else 'N/A'} "
            f"(sd={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
            f"cv={f'{cv:.1%}' if cv is not None and not math.isnan(cv) else 'N/A'})"
        )
        medians["wall_time"][label] = median
        if cv is not None and not math.isnan(cv) and cv > max_cv:
            max_cv, max_cv_desc = cv, (label, "wall_time")

        for metric in metrics:
            mean, median, stdev, cv, _, _ = collect(kept, args.region, metric, warmup=0)
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

    def plot_axis(axis_vals_by_label: dict, prefix: str, xlabel: str):
        labels = sorted(
            (l for l, v in axis_vals_by_label.items() if l.startswith(prefix) and v is not None),
            key=lambda l: axis_vals_by_label[l],
        )
        if len(labels) < 2:
            print(f"\nFewer than 2 '{prefix}*' points with a resolved axis value. Skipping plots.")
            return

        print(f"\n=== {xlabel} scaling ({prefix.rstrip('_')}) ===")
        x = np.asarray([axis_vals_by_label[l] for l in labels], dtype=float)

        for metric in ["wall_time"] + metrics:
            metric_medians = [medians[metric].get(l) for l in labels]
            if any(v is None for v in metric_medians):
                print(f"WARNING: missing '{metric}' for some {prefix}* points. Skipping plot.")
                continue

            bw_like = metric != "wall_time" and is_bandwidth_like(metric)
            ylabel = "Wall-clock time [s] per call" if metric == "wall_time" else metric
            fname = f"e3_{prefix}{'wallclock' if metric == 'wall_time' else safe_name(metric)}.png"
            stream_ceiling = args.stream_bw if bw_like else None

            # No "ideal" reference line: channel/gas count isn't a uniform
            # work unit here (gases differ widely in table size/line count,
            # channels in which tables they hit), so there's no growth curve
            # to assert as ideal -- only the measured values are plotted.
            plot_scaling(
                x, np.asarray(metric_medians, dtype=float),
                f"{ylabel} ({prefix.rstrip('_')})", "#2a78d6", res_dir, fname,
                ideal=None, stream_ceiling=stream_ceiling, xlabel=xlabel,
                yscale="linear" if bw_like else "log",
            )

    plot_axis(nd_by_label, "channels_", "Channels (ND)")
    plot_axis(ng_by_label, "gases_", "Emitters (NG)")

    print(f"\nPlots written to {res_dir}/")


if __name__ == "__main__":
    main()
