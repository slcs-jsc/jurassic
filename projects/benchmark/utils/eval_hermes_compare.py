#!/usr/bin/env python3
"""
Compare two hermes profile reports (run_hermes_profile.sh) in one diagram.

Each report directory holds log.omp<threads>.<GROUP>.txt files with the
`RUNTIME:` line printed by formod. Runtimes are averaged over the LIKWID
groups of a thread count (each group is a separate run of the same case).

Usage: eval_hermes_compare.py BASELINE_RUN_DIR OPTIMIZED_RUN_DIR [--out DIR]
"""
import argparse
import re
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from likwid_parsing import parse_formod_log
from plot_results import plot_compare_scaling

LOG_RE = re.compile(r"^log\.omp(?P<threads>\d+)\.(?P<group>[A-Za-z0-9_]+)\.txt$")
BLUE, ORANGE = "#2a78d6", "#eb6834"


def load_runtimes(run_dir: Path) -> dict:
    """Return {threads: (mean_s, stddev_s, n_groups)} for one report."""
    per_thread: dict[int, list] = {}
    for path in sorted(run_dir.glob("log.omp*.txt")):
        m = LOG_RE.match(path.name)
        if not m:
            continue
        batch = parse_formod_log(path)["batch"]
        if batch is None:
            print(f"warning: no RUNTIME line in {path}", file=sys.stderr)
            continue
        per_thread.setdefault(int(m["threads"]), []).append(batch)
    return {t: (st.mean(b["mean_s"] for b in runs),
                st.mean(b["stddev_s"] for b in runs), len(runs))
            for t, runs in sorted(per_thread.items())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("baseline", type=Path)
    ap.add_argument("optimized", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--baseline-label", default="baseline")
    ap.add_argument("--optimized-label", default="optimized")
    args = ap.parse_args()
    out = args.out or args.optimized.parent / (args.optimized.name + "_vs_baseline")

    base, opt = load_runtimes(args.baseline), load_runtimes(args.optimized)
    for name, d, p in (("baseline", base, args.baseline), ("optimized", opt, args.optimized)):
        if not d:
            sys.exit(f"No usable log.omp<N>.<GROUP>.txt runtimes in {name} report {p}")

    series = []
    for label, d, color in ((args.baseline_label, base, BLUE),
                            (args.optimized_label, opt, ORANGE)):
        series.append((label, list(d), [v[0] for v in d.values()],
                       [v[1] for v in d.values()], color))
    plot_compare_scaling(series, "formod runtime [s]", out, "runtime_compare.png",
                         title="Hermes profile: baseline vs optimized")

    out.mkdir(parents=True, exist_ok=True)
    lines = ["threads\tbaseline_s\toptimized_s\tspeedup"]
    for t in sorted(set(base) & set(opt)):
        lines.append(f"{t}\t{base[t][0]:.3f}\t{opt[t][0]:.3f}\t{base[t][0] / opt[t][0]:.3f}")
    (out / "summary.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {out}/runtime_compare.png and summary.tsv")


if __name__ == "__main__":
    main()
