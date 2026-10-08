#!/usr/bin/env python3
"""
Compare two hermes profile reports (run_hermes_profile.sh) in one diagram.

Primary metric: memory data volume per scene [GBytes], read from the
log.omp<threads>.MEM_DP.csv LIKWID marker output (region `formod`). Like
eval_scaling.py, the socket-wide volume is divided by the batch size, so it
should stay flat with thread count. A second diagram shows the formod runtime
from the `RUNTIME:` line of the matching log.omp<threads>.<GROUP>.txt files.

Usage: eval_hermes_compare.py BASELINE_RUN_DIR OPTIMIZED_RUN_DIR [--out DIR]
"""
import argparse
import re
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from likwid_parsing import parse_formod_log, parse_likwid_profile, get_metric
from plot_results import plot_compare_scaling

CSV_RE = re.compile(r"^log\.omp(?P<threads>\d+)\.(?P<group>[A-Za-z0-9_]+)\.csv$")
TXT_RE = re.compile(r"^log\.omp(?P<threads>\d+)\.(?P<group>[A-Za-z0-9_]+)\.txt$")
BATCH_RE = re.compile(r"^CPU_BATCH_SIZE=(\d+)\s*$", re.MULTILINE)
VOLUME = "Memory data volume [GBytes]"
BLUE, ORANGE = "#2a78d6", "#eb6834"


def load_volume(run_dir: Path, region: str, group: str, default_batch: int) -> dict:
    """Return {threads: GBytes per scene} for one report."""
    out = {}
    for path in sorted(run_dir.glob(f"log.omp*.{group}.csv")):
        m = CSV_RE.match(path.name)
        if not m:
            continue
        regions = parse_likwid_profile(path)["regions"]
        if region not in regions:
            sys.exit(f"Region '{region}' not in {path}. Regions found: {sorted(regions)}")
        raw = get_metric({"regions": regions, "group": group}, region, VOLUME)
        if raw is None:
            print(f"warning: no '{VOLUME}' in {path}", file=sys.stderr)
            continue
        # The profile script appends CPU_BATCH_SIZE=<n> to the matching .txt.
        txt = path.with_suffix(".txt")
        found = BATCH_RE.search(txt.read_text()) if txt.exists() else None
        batch = int(found.group(1)) if found else default_batch
        out[int(m["threads"])] = raw / batch
    return dict(sorted(out.items()))


def load_runtime(run_dir: Path) -> dict:
    """Return {threads: (mean_s, stddev_s)}, averaged over the LIKWID groups."""
    per_thread: dict[int, list] = {}
    for path in sorted(run_dir.glob("log.omp*.txt")):
        m = TXT_RE.match(path.name)
        if not m:
            continue
        batch = parse_formod_log(path)["batch"]
        if batch is not None:
            per_thread.setdefault(int(m["threads"]), []).append(batch)
    return {t: (st.mean(b["mean_s"] for b in runs), st.mean(b["stddev_s"] for b in runs))
            for t, runs in sorted(per_thread.items())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("optimized", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--region", default="formod", help="LIKWID marker region (default: formod)")
    ap.add_argument("--group", default="MEM_DP", help="LIKWID group holding the volume (default: MEM_DP)")
    ap.add_argument("--batch-size", type=int, default=48,
                    help="batch size if a log carries no CPU_BATCH_SIZE line (default: 48)")
    ap.add_argument("--baseline-label", default="baseline")
    ap.add_argument("--optimized-label", default="optimized")
    args = ap.parse_args()
    out = args.out or args.optimized.parent / (args.optimized.name + "_vs_baseline")
    out.mkdir(parents=True, exist_ok=True)
    title = "Hermes profile: baseline vs optimized"
    variants = ((args.baseline_label, args.baseline, BLUE),
                (args.optimized_label, args.optimized, ORANGE))

    # Memory data volume (primary)...
    vol = {label: load_volume(p, args.region, args.group, args.batch_size)
           for label, p, _ in variants}
    for label, p, _ in variants:
        if not vol[label]:
            sys.exit(f"No '{VOLUME}' found in the {label} report {p} "
                     f"(expected log.omp<N>.{args.group}.csv)")
    plot_compare_scaling(
        [(label, list(vol[label]), list(vol[label].values()), None, color)
         for label, _, color in variants],
        "Memory data volume per scene [GBytes]",
        out, "memory_volume_compare.png", title=title, ideal="constant")

    lines = ["threads\tbaseline_GB\toptimized_GB\treduction"]
    b, o = vol[args.baseline_label], vol[args.optimized_label]
    for t in sorted(set(b) & set(o)):
        lines.append(f"{t}\t{b[t]:.4f}\t{o[t]:.4f}\t{b[t] / o[t]:.3f}")
    (out / "summary_memory_volume.tsv").write_text("\n".join(lines) + "\n")
    print("Memory data volume per scene:\n" + "\n".join(lines))

    # Runtime (secondary, skipped if the reports carry no RUNTIME lines)...
    rt = {label: load_runtime(p) for label, p, _ in variants}
    if all(rt.values()):
        plot_compare_scaling(
            [(label, list(rt[label]), [v[0] for v in rt[label].values()],
              [v[1] for v in rt[label].values()], color) for label, _, color in variants],
            "formod runtime [s]", out, "runtime_compare.png", title=title, ideal="linear")
        lines = ["threads\tbaseline_s\toptimized_s\tspeedup"]
        b, o = rt[args.baseline_label], rt[args.optimized_label]
        for t in sorted(set(b) & set(o)):
            lines.append(f"{t}\t{b[t][0]:.3f}\t{o[t][0]:.3f}\t{b[t][0] / o[t][0]:.3f}")
        (out / "summary_runtime.tsv").write_text("\n".join(lines) + "\n")
        print("Runtime:\n" + "\n".join(lines))
    else:
        print("warning: runtime plot skipped (no RUNTIME lines in one of the reports)",
              file=sys.stderr)
    print(f"Wrote plots and summaries to {out}")


if __name__ == "__main__":
    main()
