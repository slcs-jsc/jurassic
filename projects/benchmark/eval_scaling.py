"""Evaluate run_scaling_*.sh output: runtime, speedup and parallel efficiency.

Reads formod's own "RUNTIME: ... mean= ..." line from every
out/<label>.t<threads>.b<batch>.rep<rep>.txt (older runs with a LIKWID group in
the name, <label>.t<threads>.<GROUP>.b<batch>.rep<rep>.txt, are read as well).

Two speedups are reported: the kernel (one formod_batch call, the RUNTIME mean)
and the whole application, i.e. the serial phases from the TIMER_* lines
(READ_*, FORMOD_REFERENCE, WRITE_OBS, FINALIZE) plus one batch. The gap between
them is the Amdahl limit set by the serial phases.
"""
import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_scaling
from likwid_parsing import parse_formod_log, get_stats

SERIAL_TIMERS = {"TIMER_FORMOD_REFERENCE", "TIMER_WRITE_OBS", "TIMER_FINALIZE"}

RUN_TXT_RE = re.compile(
    r"^(?P<label>[A-Za-z0-9_]+)\.t(?P<threads>\d+)(?:\.[A-Za-z0-9_]+)?"
    r"\.b(?P<batch>\d+)\.rep(?P<rep>\d+)\.txt$"
)


def load_runs(run_dir: Path) -> list:
    out_dir = run_dir / "out" if (run_dir / "out").is_dir() else run_dir
    runs = []
    for txt in sorted(out_dir.glob("*.txt")):
        m = RUN_TXT_RE.match(txt.name)
        if not m:
            continue
        log = parse_formod_log(txt)
        batch = log["batch"]
        if batch is None:
            print(f"WARNING: no RUNTIME line in {txt.name} (run failed?), skipping.")
            continue
        runs.append({
            "label": m.group("label"),
            "threads": int(m.group("threads")),
            "batch_size": int(m.group("batch")),
            "rep": int(m.group("rep")),
            "mean_s": batch["mean_s"],
            "serial_s": sum(v for k, v in log["timers"].items()
                            if k.startswith("TIMER_READ_") or k in SERIAL_TIMERS)
                        if log["timers"] else None,
        })
    return runs


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--out", type=Path, default=None,
                        help="output directory (default: <run_dir>/plots)")
    args = parser.parse_args()

    res_dir = args.out or (args.run_dir / "plots")
    res_dir.mkdir(parents=True, exist_ok=True)

    runs = load_runs(args.run_dir)
    if not runs:
        print(f"No runtime logs found under {args.run_dir}.")
        sys.exit(1)

    groups: dict[tuple, list] = {}
    for r in runs:
        groups.setdefault((r["label"], r["threads"]), []).append(r)

    medians: dict[tuple, float] = {}
    serial: dict[tuple, float] = {}
    batch_of: dict[tuple, int] = {}
    max_cv, max_cv_key = 0.0, None

    header = (f"{'label':>25} {'thr':>4} {'batch':>6} {'n':>3} | {'runtime/call [s]':>40} | "
              f"{'serial [s]':>10}")
    print(header)
    print("-" * len(header))

    for key, entries in sorted(groups.items()):
        label, threads = key
        vals = [e["mean_s"] for e in entries]
        batch_of[key] = entries[0]["batch_size"]
        if len(vals) < 2:
            print(f"{label:>25} {threads:>4} {batch_of[key]:>6} {len(vals):>3} | "
                  f"insufficient repetitions (need >= 2)")
            continue
        mean, median, stdev, cv = get_stats(vals)
        medians[key] = median
        serial_vals = [e["serial_s"] for e in entries if e["serial_s"] is not None]
        if serial_vals:
            serial[key] = float(np.median(serial_vals))
        if not math.isnan(cv) and cv > max_cv:
            max_cv, max_cv_key = cv, key
        print(f"{label:>25} {threads:>4} {batch_of[key]:>6} {len(vals):>3} | "
              f"{f'{mean:.4g} (med={median:.4g}, sd={stdev:.3g}, cv={cv:.1%})':>40} | "
              f"{serial.get(key, float('nan')):>10.3g}")

    if max_cv_key:
        print(f"\nLargest coefficient of variation: {max_cv:.1%} (label={max_cv_key[0]}, threads={max_cv_key[1]})")

    modes = sorted({m for (lbl, _) in medians for m in ("strong", "weak") if lbl.endswith(f"_{m}")})
    if not modes:
        print("\nNo strong/weak-scaling labels detected (expected '..._strong' / '..._weak'). Skipping plots.")
        return

    for mode in modes:
        print(f"\n=== Scaling analysis: {mode.upper()} ===")
        lbl_intra = f"intra_socket_{mode}"
        intra_threads = sorted(t for (lbl, t) in medians if lbl == lbl_intra)
        t1 = medians.get((lbl_intra, 1))
        if not intra_threads or t1 is None:
            print(f"No '{lbl_intra}' data with a 1-thread baseline. Skipping {mode} plots.")
            continue

        def speedup_and_efficiency(t_time, n, base=t1):
            speedup = base / t_time
            return speedup, (speedup / n if mode == "strong" else speedup)

        print(f"\n{'category':>15} {'threads':>8} {'wall_time [s]':>14} {'speedup':>9} {'efficiency':>11}")
        print("-" * 61)

        intra_speedup, intra_eff = {}, {}
        for t in intra_threads:
            s, e = speedup_and_efficiency(medians[(lbl_intra, t)], t)
            intra_speedup[t], intra_eff[t] = s, e
            print(f"{'intra_socket':>15} {t:>8} {medians[(lbl_intra, t)]:>14.4g} {s:>9.2f} {e:>10.1%}")

        app_threads = [t for t in intra_threads if (lbl_intra, t) in serial]
        app_speedup, app_eff = {}, {}
        if (lbl_intra, 1) in serial:
            app1 = serial[(lbl_intra, 1)] + t1
            print(f"\nWhole application = serial phases + one batch "
                  f"(kernel speedup vs. application speedup):")
            print(f"{'threads':>8} {'serial [s]':>11} {'batch [s]':>10} {'total [s]':>10} "
                  f"{'kernel spd':>11} {'app spd':>8} {'app eff':>8} {'serial share':>13}")
            print("-" * 86)
            for t in app_threads:
                total = serial[(lbl_intra, t)] + medians[(lbl_intra, t)]
                s, e = speedup_and_efficiency(total, t, base=app1)
                app_speedup[t], app_eff[t] = s, e
                print(f"{t:>8} {serial[(lbl_intra, t)]:>11.3g} {medians[(lbl_intra, t)]:>10.4g} "
                      f"{total:>10.4g} {intra_speedup[t]:>11.2f} {s:>8.2f} {e:>8.1%} "
                      f"{serial[(lbl_intra, t)] / total:>13.1%}")
        else:
            print("\nNo TIMER_* lines for the 1-thread run: skipping application speedup.")

        extra_speedup, extra_eff, extra_wall = [], [], []
        for prefix, tag in (("smt_socket", "SMT"), ("inter_compact", "inter (compact)"),
                            ("inter_spread", "inter (spread)")):
            lbl = f"{prefix}_{mode}"
            t = next((t for (l, t) in medians if l == lbl), None)
            if t is None:
                continue
            s, e = speedup_and_efficiency(medians[(lbl, t)], t)
            print(f"{prefix:>15} {t:>8} {medians[(lbl, t)]:>14.4g} {s:>9.2f} {e:>10.1%}")
            extra_speedup.append((t, s, tag))
            extra_eff.append((t, e, tag))
            extra_wall.append((t, medians[(lbl, t)], tag))

        x = np.asarray(intra_threads)
        if mode == "strong":
            plot_scaling(x, np.asarray([intra_speedup[t] for t in intra_threads]),
                         "Wall-clock speedup (strong scaling)", "#2a78d6", res_dir,
                         f"e2_{mode}_speedup.png", ideal="linear", higher_is_better=True,
                         smt_points=extra_speedup or None)
            eff_label = "Parallel efficiency (strong scaling, T1/(n*Tn))"
        else:
            eff_label = "Parallel efficiency (weak scaling, T1/Tn)"
        plot_scaling(x, np.asarray([intra_eff[t] for t in intra_threads]),
                     eff_label, "#2a78d6", res_dir, f"e2_{mode}_efficiency.png",
                     ideal="constant", higher_is_better=True,
                     smt_points=extra_eff or None,
                     yscale="linear", percent=True, ylim=(0, 1.1))

        if app_speedup:
            ax_t = np.asarray(app_threads)
            if mode == "strong":
                plot_scaling(ax_t, np.asarray([app_speedup[t] for t in app_threads]),
                             "Application speedup (serial phases + one batch, strong scaling)",
                             "#2a78d6", res_dir, f"e2_{mode}_app_speedup.png",
                             ideal="linear", higher_is_better=True)
            else:
                plot_scaling(ax_t, np.asarray([app_eff[t] for t in app_threads]),
                             "Application efficiency (serial phases + one batch, weak scaling)",
                             "#2a78d6", res_dir, f"e2_{mode}_app_efficiency.png",
                             ideal="constant", higher_is_better=True,
                             yscale="linear", percent=True, ylim=(0, 1.1))

        b0 = batch_of[(lbl_intra, intra_threads[0])]
        wall_label = f"Wall-clock time [s] ({b0} scenes, {mode})"
        plot_scaling(x, np.asarray([medians[(lbl_intra, t)] for t in intra_threads]),
                     wall_label, "#2a78d6", res_dir, f"e2_{mode}_wallclock_scaling.png",
                     ideal="linear" if mode == "strong" else "constant", higher_is_better=False,
                     smt_points=extra_wall or None)

    print(f"\nPlots written to {res_dir}/")


if __name__ == "__main__":
    main()
