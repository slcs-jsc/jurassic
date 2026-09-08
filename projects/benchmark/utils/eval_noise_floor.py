import argparse
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from likwid_parsing import parse_run_dir, get_metric, get_call_count, get_region_runtime 
 
DEFAULT_METRICS = [
    "Memory data volume [GBytes]",
    "Memory bandwidth [MBytes/s]",
]

def per_call(raw, call_count):
    if raw is None or call_count in (None, 0) or isinstance(raw, list):
        return None 
    return raw / call_count
 
def get_stats(values):
    mean = st.mean(values)
    median = st.median(values)
    stdev = st.pstdev(values)
    if mean == 0 : 
        cv = float("nan")
    else: 
        cv = stdev / mean
    return mean, median, stdev, cv

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--region", default="formod", help="likwid-marker region to read (default: formod)")
    parser.add_argument("--metrics", action="append", help="metric names to evaluate; "
                         "default: %s" % ", ".join(DEFAULT_METRICS))
    parser.add_argument("--warmup", type=int, default=1,
                    help="number of leading repetitions to discard (default: 1)")
    args = parser.parse_args()
    metrics = args.metrics or DEFAULT_METRICS

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
 
        # runtime per call, from the marker region info block
        rt_per_call = []
        for e in kept:
            rt = get_region_runtime(e, args.region)
            cc = get_call_count(e, args.region)
            v = per_call(rt, cc)
            if v is not None:
                rt_per_call.append(v)
        if rt_per_call:
            mean, median, stdev, cv = get_stats(rt_per_call)
            row_values.append(f"{mean:.4g}, {median:.4g}, stdev={stdev:.4g}, cv={cv:.1%}")
            if cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, "runtime/call")
        else:
            row_values.append("n/a")
 
        for metric in metrics:
            vals = []
            for e in kept:
                raw = get_metric(e, args.region, metric)
                cc = get_call_count(e, args.region)
                v = per_call(raw, cc)
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
            mean, median, stdev, cv = get_stats(vals)
            row_values.append(f"{mean:.4g}, {median:.4g}, stdev={stdev:.4g}, cv={cv:.1%}")
            if cv > max_cv:
                max_cv, max_cv_desc = cv, (label, threads, group, batch, metric)
 
        print(f"{label:>8} {threads:>4} {group:>8} {batch:>6} | " +
              " | ".join(f"{v:>22}" for v in row_values))
 
    if max_cv_desc:
        print(f"Largest observed coefficient of variation: {max_cv:.1%}  "
              f"(config={max_cv_desc[:4]}, metric='{max_cv_desc[4]}')")
        print(f"\nRecommended significance threshold for further comparisons: "
              f"differences below ~{max_cv:.1%} of the median are not "
              f"distinguishable from measurement noise at this repetition count.")
    else:
        print("No metric produced an usable coefficient of variation.")

if __name__ == "__main__":
    main()


