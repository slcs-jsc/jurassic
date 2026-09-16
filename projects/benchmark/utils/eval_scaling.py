import argparse
import statistics as st
import sys
from pathlib import Path
 
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_scaling, boxplot
from likwid_parsing import parse_run_dir, collect_runtime, collect
 
DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

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
    medians: dict[str, dict[tuple, float]] = {m: {} for m in ["wall_time", "runtime/call"] + metrics}
    values: dict[str, dict[tuple, list]] = {m: {} for m in ["wall_time", "runtime/call"] + metrics}
    batch_size_dict = {}
    
    header = f"{'label':>8} {'thr':>4} {'group':>8} {'batch':>6} | " + \
             " | ".join(f"{m[:22]:>22}" for m in ["runtime/call [s]"] + metrics)
    print(header)
    print("-" * len(header))

    for (label, threads, group, batch), entries in sorted(groups.items()):
        batch_size_dict[(label, threads)] = batch
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
        mean, median, stdev, cv, _, vals =  collect_runtime(kept, warmup=1)
        row_values.append(
            f"{f'{mean:.4g}' if mean is not None else 'N/A'}, "
            f"{f'{median:.4g}' if median is not None else 'N/A'}, "
            f"stdev={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
            f"cv={f'{cv:.1%}' if cv is not None else 'N/A'})"
        )
        values["wall_time"][(label, threads)] = vals
        medians["wall_time"][(label, threads)] = median

        if cv > max_cv:
            max_cv, max_cv_desc = cv, (label, threads, group, batch, "wall_time")
        else:
            row_values.append("n/a")

        for metric in metrics:
            mean, median, stdev, cv, _, vals = collect(kept, args.region, metric, warmup=1)
            row_values.append(
                f"{f'{mean:.4g}' if mean is not None else 'N/A'}, "
                f"{f'{median:.4g}' if median is not None else 'N/A'}, "
                f"stdev={f'{stdev:.4g}' if stdev is not None else 'N/A'}, "
                f"cv={f'{cv:.1%}' if cv is not None else 'N/A'})"
            )
            values[metric][(label, threads)] = vals
            medians[metric][(label, threads)] = median
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
    spread_threads = sorted({t for (lbl, t) in medians["wall_time"] if lbl == "spread"})

    t1 = medians["wall_time"].get(("phys", 1))
    if t1 is not None:
        # Wall-clock speedup relative to the 1-thread physical-core baseline
        phys_speedups = {t: t1 / medians["wall_time"][("phys", t)] for t in phys_threads
                          if medians["wall_time"].get(("phys", t)) is not None}
        smt_speedups = {t: t1 / medians["wall_time"][("smt", t)] for t in smt_threads
                         if medians["wall_time"].get(("smt", t)) is not None}
        spread_speedups = {t: t1 / medians["wall_time"][("spread", t)] for t in spread_threads
                            if medians["wall_time"].get(("spread", t)) is not None}

        if len(phys_threads) > 1:
            print()
            print(f"{'threads':>8} {'wall_time [s]':>14} {'speedup':>9} {'efficiency':>11}")
            print("-" * 55)

            for lbl, speedups in (("phys", phys_speedups), ("spread", spread_speedups)):
                for t in sorted(speedups):
                    tn = medians["wall_time"].get((lbl, t))
                    print(f"{lbl:>8} {t:>8} {tn:>14.4g} {speedups[t]:>9.2f} {speedups[t] / t:>10.1%}")

            for t in sorted(smt_speedups):
                tn = medians["wall_time"].get(("smt", t))
                print(f"{'smt':>8} {t:>8} {tn:>14.4g} {smt_speedups[t]:>9.2f} {'(SMT)':>11}")

        # Plot wall-clock speedup
        smt_speedup_points = [(t, smt_speedups[t], "SMT") for t in sorted(smt_speedups)]
        if phys_speedups:
            plot_scaling(
                np.asarray(sorted(phys_speedups)),
                np.asarray([phys_speedups[t] for t in sorted(phys_speedups)]),
                "Wall-clock speedup",
                "#efb239",
                res_dir,
                "e2_wallclock_speedup.png",
                ideal="linear",
                higher_is_better=True,
                smt_points=smt_speedup_points,
            )
        if spread_speedups:
            plot_scaling(
                np.asarray(sorted(spread_speedups)),
                np.asarray([spread_speedups[t] for t in sorted(spread_speedups)]),
                "Wall-clock speedup (spread across sockets)",
                "#89dd29",
                res_dir,
                "e2_wallclock_speedup_spread.png",
                ideal="linear",
                higher_is_better=True,
            )

        # Plot efficiency (speedup/p) of physical threads
        eff_threads = sorted(phys_speedups)
        eff_values = [phys_speedups[t] / t for t in eff_threads]
        if eff_threads:
            fig, ax = plt.subplots()
            ax.plot(eff_threads, eff_values, "o-", color="#efb239", label="phys", zorder=3)
            ax.axhline(1.0, linestyle="--", color="#898781", label="Ideal (100%)", zorder=2)
            ax.set_xlabel("threads")
            ax.set_ylabel("Parallel efficiency (speedup / threads)")
            ax.set_xscale("log")
            ax.set_xticks(eff_threads)
            ax.xaxis.set_major_formatter(plt.ScalarFormatter())
            ax.set_ylim(0, 1.15)
            ax.set_title("Parallel efficiency vs thread count")
            ax.legend(fontsize=8)
            fig.tight_layout()
            res_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(res_dir / "e2_efficiency.png", bbox_inches="tight")
            plt.close(fig)
            print(f"  wrote {res_dir / 'e2_efficiency.png'}")
    else:
        print("\nNo 1-thread baseline found. Skipping speedup plot.")

    def _smt_point(metric_name):
        if not smt_threads:
            return None
        t = smt_threads[0]
        v = medians[metric_name].get(("smt", t))
        return [(t, v, "SMT (48 threads)")] if v is not None else None

    for metric_name in ["wall_time"] + metrics:
        phys_vals_med = [medians[metric_name].get(("phys", t)) for t in phys_threads]
        phys_vals = [values[metric_name].get(("phys", t)) for t in phys_threads]
        if any(v is None for v in phys_vals_med):
            print(f"WARNING: missing '{metric_name}' for some phys thread counts. Skipping plot.")
            continue

        lower = metric_name.lower()
        if metric_name == "wall_time":
            b0 = batch_size_dict.get(("phys", phys_threads[0]), "?")
            ideal, higher_is_better, ceiling = "linear", False, None
            fname, ylabel, color = "e2_wallclock_scaling.png", f"Wall-clock time [s] ({b0} scenes)", "#efb239"

            plot_scaling(phys_threads, phys_vals_med, ylabel, color, res_dir, fname,
                                ideal=ideal, higher_is_better=higher_is_better,
                                stream_ceiling=ceiling, smt_points=_smt_point(metric_name))
            
            spread_vals = [medians[metric_name].get(("spread", t)) for t in spread_threads]
            if any(v is None for v in spread_vals):
                print(f"WARNING: missing '{metric_name}' for some phys thread counts. Skipping plot.")
                continue
    
            spread_fname = fname.replace(".png", "_spread.png")
            plot_scaling(spread_threads, spread_vals, f"{ylabel} (spread)", color, res_dir, spread_fname,
                                ideal=ideal, higher_is_better=higher_is_better,
                                stream_ceiling=ceiling)
            
        elif "volume" in lower:
            ideal, higher_is_better, ceiling = "constant", True, None
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", f"{metric_name} (socket-wide ÷ batch-size)", "#3ab9dc"
            
        elif "bandwidth" in lower:
            ideal, higher_is_better, ceiling = None, True, args.stream_bw
            if ceiling is None:
                print(f"NOTE: --stream-bw not given, plotting '{metric_name}' without a reference ceiling.")
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", f"{metric_name} (socket-wide, not normalized)", "#c76ce0"

            print(phys_vals)
            boxplot(phys_threads, phys_vals, ylabel, res_dir, fname)

        elif metric_name.startswith("CAS_COUNT"):
            ideal, higher_is_better, ceiling = "constant", False, None
            fname = f"e2_{metric_name.lower()}_scaling.png"
            ylabel = f"{metric_name} [GBytes-equiv] (socket-wide / batch-size)"
            color = "#c0392b" if metric_name.endswith("_RD") else "#2980b9"
        else:
            ideal, higher_is_better, ceiling = "linear", True, None
            safe = metric_name.split("[")[0].strip().replace(" ", "_").lower()
            fname, ylabel, color = f"e2_{safe}_scaling.png", f"{metric_name}/call", "#7d8f69"

    print(f"\nPlots written to {res_dir}/")

if __name__ == "__main__":
    main()
