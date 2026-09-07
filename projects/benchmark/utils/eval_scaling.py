import argparse
import statistics as st
import sys
from pathlib import Path
 
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import _plot_scaling
from likwid_parsing import parse_run_dir, get_metric, get_call_count, get_region_runtime
 
DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

def per_call(raw, call_count):
    if raw is None or call_count in (None, 0) or isinstance(raw, list):
        return None 
    return raw / call_count

def coefficient_of_variation(values):
    if len(values) < 2:
        return float("nan")
    mean = st.mean(values)
    if mean == 0:
        return float("nan")
    return st.pstdev(values) / mean

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--region", default="formod", help="likwid-marker region to read (default: formod)")
    parser.add_argument("--metrics", action="append", help="metric names to evaluate; "
                            "default: %s" % ", ".join(DEFAULT_METRICS))
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
    configs = parse_run_dir(args.run_dir)
    if not configs:
        print(f"No parsed configs found under {args.run_dir}/out/. ")
        sys.exit(1)

    # Group by everything except rep.
    groups: dict[tuple, list] = {}
    for c in configs:
        key = (c["label"], c["threads"], c["group"], c["batch_size"])
        groups.setdefault(key, []).append(c)

    # Coefficient of Variation
    max_cv = 0.0
    max_cv_desc = None
    warned_metrics = set()

    medians: dict[str, dict[tuple, float]] = {m: {} for m in ["wall_time", "runtime/call"] + metrics}
 
    header = f"{'label':>8} {'thr':>4} {'group':>8} {'batch':>6} | " + \
             " | ".join(f"{m[:22]:>22}" for m in ["runtime/call [s]"] + metrics)
    print(header)
    print("-" * len(header))

    for (label, threads, group, batch), entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda e: e["rep"])
        kept = entries[args.warmup:]
        if len(kept) < 2:
            print(f"{label:>8} {threads:>4} {group:>8} {batch:>6} | "
                    f"insufficient repetitions after warmup discard "
                    f"({len(kept)} left, need >= 2)")
            continue

        row_values = []

        # total wall-clock time for the batch call, from omp_get_wtime()
        # batch size is held FIXED across thread counts, so this is directly the time-to-solution for the same problem.
        wall_time = []
        for e in kept:
            b = e.get("batch")
            if b:
                wall_time.append(b["mean_s"])
        if wall_time:
            med, cv = st.median(wall_time), coefficient_of_variation(wall_time)
            row_values.append(f"{med:.4g} (cv={cv:.1%})")
            medians["wall_time"][(label, threads)] = med
            if cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, "wall_time")
        else:
            row_values.append("n/a")

        #  marker-based runtime/call, kept for reference only
        rt_per_call = []
        for e in kept:
            rt = get_region_runtime(e, args.region)
            cc = get_call_count(e, args.region)
            v = per_call(rt, cc)
            if v is not None:
                rt_per_call.append(v)
        if rt_per_call:
            med, cv = st.median(rt_per_call), coefficient_of_variation(rt_per_call)
            row_values.append(f"{med:.4g} (cv={cv:.1%})")
            medians["runtime/call"][(label, threads)] = med
            if cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, "runtime/call")
        else:
            row_values.append("n/a")

        for metric in metrics:
            vals = []
            for e in kept:
                raw = get_metric(e, args.region, metric)
                cc = get_call_count(e, args.region)
                if "bandwidth" not in  metric.lower(): 
                    v = per_call(raw, cc)
                else: 
                    v = raw
                if v is not None:
                    vals.append(v)
            if not vals:
                if metric not in warned_metrics:
                    available = set()
                    for e in kept:
                        t = e.get("regions", {}).get(args.region, {}).get("tables", {})
                        for tbl in t.values():
                            available.update(tbl.keys())
                    print(f"  WARNING: metric '{metric}' not found for group "
                          f"'{group}', region '{args.region}'. Available metrics: "
                          f"{sorted(available)}", file=sys.stderr)
                    warned_metrics.add(metric)
                row_values.append("n/a")
                continue
            med, cv = st.median(vals), coefficient_of_variation(vals)
            row_values.append(f"{med:.4g} (cv={cv:.1%})")
            medians[metric][(label, threads)] = med
            if cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, metric)
 
        print(f"{label:>8} {threads:>4} {group:>8} {batch:>6} | " +
              " | ".join(f"{v:>22}" for v in row_values))

    if max_cv_desc:
        print(f"Largest observed coefficient of variation: {max_cv:.1%}  "
            f"(config={max_cv_desc[:4]}, metric='{max_cv_desc[4]}')")
    else:
        print("No metric produced an usable coefficient of variation.")

    phys_threads = sorted({t for (lbl, t) in medians["wall_time"] if lbl == "phys"})
    smt_threads = sorted({t for (lbl, t) in medians["wall_time"] if lbl == "smt"})
    t1 = medians["wall_time"].get(("phys", 1))

    if t1 is not None:
        # Wall-clock speedup relative to the 1-thread physical-core baseline
        phys_speedups = {}
        for t in phys_threads:
            tn = medians["wall_time"].get(("phys", t))
            if tn is not None:
                phys_speedups[t] = t1 / tn

        smt_speedups = {}
        for t in smt_threads:
            tn = medians["wall_time"].get(("smt", t))
            if tn is not None:
                smt_speedups[t] = t1 / tn

        if len(phys_threads) > 1:
            print()
            print(f"{'threads':>8} {'wall_time [s]':>14} {'speedup':>9} {'efficiency':>11}")
            print("-" * 45)

            for t in phys_threads:
                tn = medians["wall_time"].get(("phys", t))
                speedup = phys_speedups.get(t)

                if tn is None or speedup is None:
                    continue

                print(
                    f"{t:>8} {tn:>14.4g} "
                    f"{speedup:>9.2f} {speedup / t:>10.1%}"
                )

            for t in smt_threads:
                tn = medians["wall_time"].get(("smt", t))
                speedup = smt_speedups.get(t)

                if tn is None or speedup is None:
                    continue

                print(
                    f"{t:>8} {tn:>14.4g} "
                    f"{speedup:>9.2f} {'(SMT)':>11}"
                )

        # Plot wall-clock speedup
        speedup_threads = sorted(phys_speedups)
        speedup_values = [phys_speedups[t] for t in speedup_threads]

        smt_speedup_points = [(t, smt_speedups[t], "SMT") for t in sorted(smt_speedups)]

        if speedup_threads:
            _plot_scaling(
                np.asarray(speedup_threads),
                np.asarray(speedup_values),
                "Wall-clock speedup",
                "#efb239",
                res_dir,
                "e2_wallclock_speedup.png",
                ideal="linear",
                higher_is_better=True,
                smt_points=smt_speedup_points,
            )

    else:
        print("\nNo 1-thread baseline found. Skipping speedup plot.")

    def _smt_point(metric_name):
        if not smt_threads:
            return None
        t = smt_threads[0]
        v = medians[metric_name].get(("smt", t))
        return [(t, v, "SMT (48 threads)")] if v is not None else None

    for metric_name in ["wall_time"] + metrics:
        vals = [medians[metric_name].get(("phys", t)) for t in phys_threads]
        if any(v is None for v in vals):
            print(f"WARNING: missing '{metric_name}' for some thread counts. Skipping plot.")
            continue

        lower = metric_name.lower()
        if metric_name == "wall_time":
            ideal, higher_is_better, ceiling = "linear", False, None
            fname, ylabel, color = "e2_wallclock_scaling.png", "Wall-clock time [s]", "#efb239"
        elif "volume" in lower:
            ideal, higher_is_better, ceiling = "constant", True, None
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", f"{metric_name}/call", "#3ab9dc"

        elif "bandwidth" in lower:
            ideal, higher_is_better, ceiling = None, True, args.stream_bw
            if ceiling is None:
                print(f"NOTE: --stream-bw not given, plotting '{metric_name}' without a reference ceiling.")
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", metric_name, "#c76ce0"
        else:
            ideal, higher_is_better, ceiling = "linear", True, None
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", f"{metric_name}/call", "#7d8f69"

        _plot_scaling(phys_threads, vals, ylabel, color, res_dir, fname,
                    ideal=ideal, higher_is_better=higher_is_better,
                    stream_ceiling=ceiling, smt_points=_smt_point(metric_name))

    print(f"\nPlots written to {res_dir}/")

if __name__ == "__main__":
    main()

