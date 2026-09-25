import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_scaling, boxplot
from likwid_parsing import parse_run_dir, collect_runtime, collect

DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

def metric_plot_style(metric_name: str):
    """Return (ideal, higher_is_better, is_ceiling_metric, filename_stub, ylabel, color)."""
    lower = metric_name.lower()
    safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
    if "volume" in lower or "energy" in lower:
        return "constant", True, False, safe, f"{metric_name} (socket-wide / batch-size)", "#3ab9dc"
    if "bandwidth" in lower or "mflop/s" in lower:
        return None, True, True, safe, f"{metric_name} (socket-wide)", "#c76ce0"
    if metric_name.startswith("CAS_COUNT"):
        color = "#c0392b" if metric_name.endswith("_RD") else "#2980b9"
        return "constant", False, False, metric_name.lower(), f"{metric_name} [GBytes-equiv] (socket-wide / batch-size)", color
    return "linear", True, False, safe, f"{metric_name}/call", "#7d8f69"

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--region", default="formod", help="likwid-marker region to read (default: formod)")
    parser.add_argument("--metrics", action="append",
                         help="metric names to evaluate; default: %s" % ", ".join(DEFAULT_METRICS))
    parser.add_argument("--warmup", type=int, default=1,
                         help="number of leading repetitions to discard (default: 1)")
    parser.add_argument("--stream-bw", type=float, default=None,
                         help="measured STREAM bandwidth ceiling [MBytes/s] for this node, "
                              "e.g. from likwid-bench -t load_avx. If omitted, bandwidth-like "
                              "metrics are plotted without a reference ceiling.")
    parser.add_argument("--out", type=Path, default=None,
                         help="output directory (default: <run_dir>/plots)")
    args = parser.parse_args()

    metrics = args.metrics or DEFAULT_METRICS
    res_dir = args.out or (args.run_dir / "plots")
    res_dir.mkdir(parents=True, exist_ok=True)

    configs = parse_run_dir(args.run_dir)
    if not configs:
        print(f"No parsed configs found under {args.run_dir}.")
        sys.exit(1)

    groups: dict[tuple, list] = {}
    for c in configs:
        key = (c["label"], c["threads"], c["group"], c["batch_size"])
        groups.setdefault(key, []).append(c)

    metric_keys = ["wall_time"] + metrics
    medians: dict[str, dict[tuple, float]] = {m: {} for m in metric_keys}
    values: dict[str, dict[tuple, list]] = {m: {} for m in metric_keys}
    batch_size_dict: dict[tuple, int] = {}

    max_cv, max_cv_desc = 0.0, None
    detected_modes = set()

    header = (f"{'label':>25} {'thr':>4} {'group':>8} {'batch':>6} | " +
              " | ".join(f"{m[:22]:>22}" for m in ["runtime/call [s]"] + metrics))
    print(header)
    print("-" * len(header))

    for (label, threads, group, batch), entries in sorted(groups.items()):
        if label.endswith("_strong"):
            detected_modes.add("strong")
        if label.endswith("_weak"):
            detected_modes.add("weak")

        batch_size_dict[(label, threads)] = batch
        entries = sorted(entries, key=lambda e: e["rep"])
        kept = entries[args.warmup:]
        if len(kept) < 2:
            print(f"{label:>25} {threads:>4} {group:>8} {batch:>6} | "
                  f"insufficient repetitions after warmup discard ({len(kept)} left, need >= 2)")
            continue

        row_values = []
        mean, median, stdev, cv, n, vals = collect_runtime(kept, warmup=0)
        row_values.append(
            f"{f'{mean:.4g}' if mean is not None else 'N/A'} "
            f"(med={f'{median:.4g}' if median is not None else 'N/A'}, "
            f"sd={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
            f"cv={f'{cv:.1%}' if cv is not None and not math.isnan(cv) else 'N/A'})"
        )
        values["wall_time"][(label, threads)] = vals
        medians["wall_time"][(label, threads)] = median
        if cv is not None and not math.isnan(cv) and cv > max_cv:
            max_cv, max_cv_desc = cv, (label, threads, group, batch, "wall_time")

        for metric in metrics:
            mean, median, stdev, cv, n, vals = collect(kept, args.region, metric, warmup=0)
            row_values.append(
                f"{f'{mean:.4g}' if mean is not None else 'N/A'} "
                f"(med={f'{median:.4g}' if median is not None else 'N/A'}, "
                f"sd={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
                f"cv={f'{cv:.1%}' if cv is not None and not math.isnan(cv) else 'N/A'})"
            )
            values[metric][(label, threads)] = vals
            medians[metric][(label, threads)] = median
            if cv is not None and not math.isnan(cv) and cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, metric)

        print(f"{label:>25} {threads:>4} {group:>8} {batch:>6} | " +
              " | ".join(f"{v:>22}" for v in row_values))

    if max_cv_desc:
        print(f"\nLargest observed coefficient of variation: {max_cv:.1%} "
              f"(config={max_cv_desc[:4]}, metric='{max_cv_desc[4]}')")
    else:
        print("\nNo metric produced a usable coefficient of variation.")

    if not detected_modes:
        print("\nNo strong/weak-scaling labels detected (expected '..._strong' / '..._weak'). "
              "Skipping scaling plots.")
        return

    for mode in sorted(detected_modes):
        print(f"\n=== Scaling analysis: {mode.upper()} ===")

        lbl_intra = f"intra_socket_{mode}"
        lbl_smt = f"smt_socket_{mode}"
        lbl_compact = f"inter_compact_{mode}"
        lbl_spread = f"inter_spread_{mode}"

        intra_threads = sorted({t for (lbl, t) in medians["wall_time"] if lbl == lbl_intra})
        if not intra_threads:
            print(f"No '{lbl_intra}' data found. Skipping.")
            continue

        t1 = medians["wall_time"].get((lbl_intra, 1))
        if t1 is None:
            print(f"No 1-thread baseline for '{lbl_intra}'. Skipping speedup/efficiency plots for {mode}.")
            continue

        def speedup_and_efficiency(t_time, n_threads):
            speedup = t1 / t_time
            efficiency = speedup / n_threads if mode == "strong" else speedup
            return speedup, efficiency

        print(f"\n{'category':>15} {'threads':>8} {'wall_time [s]':>14} {'speedup':>9} {'efficiency':>11}")
        print("-" * 65)

        intra_speedups: dict[int, float] = {}
        intra_efficiency: dict[int, float] = {}
        for t in intra_threads:
            t_time = medians["wall_time"][(lbl_intra, t)]
            speedup, efficiency = speedup_and_efficiency(t_time, t)
            intra_speedups[t] = speedup
            intra_efficiency[t] = efficiency
            print(f"{'intra_socket':>15} {t:>8} {t_time:>14.4g} {speedup:>9.2f} {efficiency:>10.1%}")

        extra_speedup_points = []   # (threads, speedup, tag)
        extra_wall_points = []      # (threads, wall_time_s, tag)
        for lbl, tag in ((lbl_smt, "SMT"), (lbl_compact, "inter (compact)"), (lbl_spread, "inter (spread)")):
            t_target = next((t for (l, t) in medians["wall_time"] if l == lbl), None)
            if t_target is None:
                continue
            t_time = medians["wall_time"][(lbl, t_target)]
            speedup, efficiency = speedup_and_efficiency(t_time, t_target)
            print(f"{lbl:>15} {t_target:>8} {t_time:>14.4g} {speedup:>9.2f} {efficiency:>10.1%}")
            extra_speedup_points.append((t_target, speedup, tag))
            extra_wall_points.append((t_target, t_time, tag))

        intra_t_arr = np.asarray(intra_threads)

        if mode == "strong":
            plot_scaling(
                intra_t_arr,
                np.asarray([intra_speedups[t] for t in intra_threads]),
                "Wall-clock speedup (strong scaling)", "#efb239", res_dir,
                f"e2_{mode}_speedup.png", ideal="linear", higher_is_better=True,
                smt_points=extra_speedup_points or None,
            )
        else:
            plot_scaling(
                intra_t_arr,
                np.asarray([intra_efficiency[t] for t in intra_threads]),
                "Parallel efficiency (weak scaling, T1/Tn)", "#efb239", res_dir,
                f"e2_{mode}_efficiency.png", ideal="constant", higher_is_better=True,
                smt_points=extra_speedup_points or None,
            )

        for metric in ["wall_time"] + metrics:
            metric_medians = [medians[metric].get((lbl_intra, t)) for t in intra_threads]
            if any(v is None for v in metric_medians):
                print(f"WARNING: missing '{metric}' for some {lbl_intra} thread counts. Skipping plot.")
                continue

            if metric == "wall_time":
                b0 = batch_size_dict.get((lbl_intra, intra_threads[0]), "?")
                plot_scaling(
                    intra_t_arr, np.asarray(metric_medians, dtype=float),
                    f"Wall-clock time [s] ({b0} scenes, {mode})", "#efb239", res_dir,
                    f"e2_{mode}_wallclock_scaling.png", ideal="linear", higher_is_better=False,
                    smt_points=extra_wall_points or None,
                )
                continue

            ideal, higher_is_better, is_ceiling, safe, ylabel, color = metric_plot_style(metric)
            fname = f"e2_{mode}_{safe}_scaling.png"
            stream_ceiling = args.stream_bw if is_ceiling else None
            if is_ceiling and stream_ceiling is None:
                print(f"NOTE: --stream-bw not given, plotting '{metric}' without a reference ceiling.")

            plot_scaling(
                intra_t_arr, np.asarray(metric_medians, dtype=float),
                f"{ylabel} ({mode})", color, res_dir, fname,
                ideal=ideal, higher_is_better=higher_is_better,
                stream_ceiling=stream_ceiling,
            )

            metric_value_lists = [values[metric].get((lbl_intra, t)) or [] for t in intra_threads]
            if any(metric_value_lists):
                boxplot(intra_t_arr, metric_value_lists, f"{ylabel} ({mode})", res_dir,
                        fname.replace(".png", "_box.png"))

    print(f"\nPlots written to {res_dir}/")


if __name__ == "__main__":
    main()