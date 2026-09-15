import argparse
import statistics as st
import sys
from pathlib import Path
 
import numpy as np
import matplotlib.pyplot as plt
 
sys.path.insert(0, str(Path(__file__).resolve().parent))
from likwid_parsing import parse_run_dir, collect, collect_runtime
 
DEFAULT_METRICS = [
    "Memory write data volume [GBytes]",   
    "Memory read data volume [GBytes]",    
    "Memory data volume [GBytes]", 
    "CAS_COUNT_WR",
    "CAS_COUNT_RD",
]

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--baseline", default="base", help="baseline variant label (default: base)")
    ap.add_argument("--variant", default="nomemset", help="changed variant label (default: nomemset)")
    ap.add_argument("--region", default="formod")
    ap.add_argument("--metrics", action="append")
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    metrics = args.metrics or DEFAULT_METRICS
    res_dir = args.out or (args.run_dir / "plots")

    configs = parse_run_dir(args.run_dir)
    if not configs:
        print(f"No parsed configs found under {args.run_dir}/out/.", file=sys.stderr)
        sys.exit(1)

    labels = sorted({c["label"] for c in configs})
    for needed in (args.baseline, args.variant):
        if needed not in labels:
            print(f"Label '{needed}' not present in this run. Labels found: {labels}",
                  file=sys.stderr)
            sys.exit(1)

    # (label, threads, group) -> entries
    groups: dict[tuple, list] = {}
    for c in configs:
        groups.setdefault((c["label"], c["threads"], c["group"]), []).append(c)

    thread_counts = sorted({t for (_, t, _) in groups})
    groups_present = sorted({g for (_, _, g) in groups})

    plot_data: dict[str, dict[int, tuple]] = {}

    for threads in thread_counts:
        
        # total wall-clock time for the batch call, from omp_get_wtime()
        b_entries = groups.get((args.baseline, threads, groups_present[0]), [])
        v_entries = groups.get((args.variant, threads, groups_present[0]), [])
        if b_entries and v_entries:
            _, bm, _, bcv, bn = collect_runtime(b_entries, args.region, args.warmup)
            _, vm, _, vcv, vn = collect_runtime(v_entries, args.region, args.warmup)
            if bm and vm:
                change = (vm - bm) / bm
                noise = max(bcv, vcv) if not (np.isnan(bcv) or np.isnan(vcv)) else float("nan")
                print(f"{'wall_time [s]':>38} {bm:>13.5g} {vm:>13.5g} "
                      f"{change:>+8.1%} {noise:>7.1%}")
                plot_data.setdefault("wall_time [s]", {})[threads] = (bm, vm)

        for group in groups_present:
            b_entries = groups.get((args.baseline, threads, group), [])
            v_entries = groups.get((args.variant, threads, group), [])
            if not b_entries or not v_entries:
                continue
            for metric in metrics:
                _, bm, _, bcv, bn = collect(b_entries, args.region, metric, args.warmup)
                _, vm, _, vcv, vn = collect(v_entries, args.region, metric, args.warmup)
                if bm is None or vm is None:
                    continue
                if bn < 2 or vn < 2:
                    print(f"{metric[:38]:>38} {bm:>13.5g} {vm:>13.5g} "
                          f"{'':>9} {'':>8}  too few reps after warmup")
                    continue
                change = (vm - bm) / bm if bm != 0 else float("nan")
                noise = max(bcv, vcv)
                print(f"{metric[:38]:>38} {bm:>13.5g} {vm:>13.5g} "
                      f"{change:>+8.1%} {noise:>7.1%}")
                plot_data.setdefault(metric, {})[threads] = (bm, vm)

    _plot_ab(plot_data, args.baseline, args.variant, res_dir)
    print(f"\nPlots written to {res_dir}/")


def _plot_ab(plot_data, baseline, variant, res_dir):
    """Grouped bar chart per metric: baseline vs variant, one pair per thread count."""
    if not plot_data:
        return
    res_dir.mkdir(parents=True, exist_ok=True)
    for metric, by_threads in plot_data.items():
        threads = sorted(by_threads)
        base_vals = [by_threads[t][0] for t in threads]
        var_vals = [by_threads[t][1] for t in threads]
 
        x = np.arange(len(threads), dtype=float)
        width = 0.38
        fig, ax = plt.subplots()
        ax.bar(x - width / 2, base_vals, width, label=baseline, color="#3ab9dc")
        ax.bar(x + width / 2, var_vals, width, label=variant, color="#efb239")
        ax.set_xticks(x)
        ax.set_xticklabels([str(t) for t in threads])
        ax.set_xlabel("threads")
        ax.set_ylabel(metric)
        ax.set_title(f"{metric}: {baseline} vs {variant}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        safe = (metric.split("[")[0].strip().replace(" ", "_").replace("/", "_").lower())
        fname = f"ab_{safe}.png"
        fig.savefig(res_dir / fname, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {res_dir / fname}")

if __name__ == "__main__":
    main() 