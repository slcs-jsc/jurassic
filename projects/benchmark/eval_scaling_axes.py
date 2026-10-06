"""Evaluate scaling_axes.sh runs.

Efficiency = throughput(n) / (n * throughput(1)); for strong runs T1 is the batch
x the 1-thread time per scene. Per axis (geometry, channels@<case>, gases@<case>)
it prints tables and writes e4_<axis>_*.png: strong speedup/efficiency vs threads
for settings with a full curve (the baselines), strong efficiency vs ND or NG,
single-thread cost, and the weak-scaling plots if weak runs exist.
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from eval_scaling import load_runs
from plot_results import plot_series
from plot_style import PALETTE

AXIS_X = {"geometry": None, "channels": ("nd", "Channels (ND)"), "gases": ("ng", "Emitters (NG)")}


def read_kv(path: Path) -> dict:
    if not path.is_file():
        return {}
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--out", type=Path, default=None, help="output directory (default: <run_dir>/plots)")
    args = parser.parse_args()

    res_dir = args.out or (args.run_dir / "plots")
    res_dir.mkdir(parents=True, exist_ok=True)

    # a single run, or an array run with one subdirectory per case
    run_dirs = ([args.run_dir] if (args.run_dir / "settings.tsv").is_file()
                else sorted(d for d in args.run_dir.iterdir() if (d / "settings.tsv").is_file()))
    if not run_dirs:
        print(f"No settings.tsv in {args.run_dir} or its subdirectories.")
        sys.exit(1)
    settings, runs = {}, []
    for d in run_dirs:
        with open(d / "settings.tsv") as f:
            settings.update({r["setting"]: r for r in csv.DictReader(f, delimiter="\t")})
        runs += load_runs(d)
    per_socket = int(read_kv(run_dirs[0] / "topology.txt").get("phys_cores_per_socket", 0)) or None

    # (setting, placement, threads) -> lists of throughput / serial time
    thr: dict[tuple, list] = {}
    serial: dict[tuple, list] = {}
    for r in runs:
        setting, _, placement = r["label"].rpartition("_")
        if setting not in settings:
            continue
        key = (setting, placement, r["threads"])
        thr.setdefault(key, []).append(r["batch_size"] / r["mean_s"])
        if r["serial_s"] is not None:
            serial.setdefault(key, []).append(r["serial_s"])
    if not thr:
        print(f"No usable runs under {args.run_dir}.")
        sys.exit(1)
    thr = {k: float(np.median(v)) for k, v in thr.items()}
    serial = {k: float(np.median(v)) for k, v in serial.items()}

    def eff(setting, placement, n):
        base = thr.get((setting, "compact", 1))
        t = thr.get((setting, placement, n))
        return t / (n * base) if base and t else None

    def fmt(v, spec):
        return format(v, spec) if v is not None else format("n/a", ">" + spec.split(".")[0])

    all_threads = sorted({n for (_, p, n) in thr if p == "compact"})
    has_weak = len(all_threads) > 1
    full = max(n for (_, p, n) in thr if p in ("compact", "strongcompact"))
    one_socket = per_socket if any(n == per_socket for (_, p, n) in thr if p != "spread") else None

    # axis tags: "geometry", "channels@<geom>", "gases@<geom>" (older runs: "channels", "gases")
    tags = list(dict.fromkeys(t for r in settings.values() for t in r["axes"].split(",")))
    tags.sort(key=lambda t: list(AXIS_X).index(t.split("@")[0]))
    for axis in tags:
        xinfo = AXIS_X[axis.split("@")[0]]
        axis_name, stem = axis.replace("@", " around "), axis.replace("@", "_")
        members = [s for s, r in settings.items() if axis in r["axes"].split(",")]
        if not members:
            continue
        if xinfo:
            members.sort(key=lambda s: int(settings[s][xinfo[0]]))

        print(f"\n=== {axis_name} ===")
        if has_weak:
            print(f"{'setting':>26} {'geom':>7} {'ND':>4} {'NG':>3} {'rays':>4} | {'s/scene':>8} "
              f"{'ms/ray':>7} | {'eff@' + str(one_socket or '-'):>8} {'eff@' + str(full):>8} "
              f"{'spread':>7} {'2-socket':>9} | {'serial@1':>8} {f'serial@{full}':>10}")
            print("-" * 122)

        series, markers, summary = [], [], {"1 socket": [], "all sockets": [], "spread": []}
        for i, s in enumerate(members):
            r = settings[s]
            t1 = thr.get((s, "compact", 1))
            per_scene = 1 / t1 if t1 else None
            e_one = eff(s, "compact", one_socket) if one_socket else None
            e_full = eff(s, "compact", full)
            e_spread = eff(s, "spread", one_socket) if one_socket else None
            gain = (thr[(s, "spread", one_socket)] / thr[(s, "compact", one_socket)]
                    if one_socket and (s, "spread", one_socket) in thr
                    and (s, "compact", one_socket) in thr else None)
            if has_weak:
                print(f"{s:>26} {r['geometry']:>7} {r['nd']:>4} {r['ng']:>3} {r['nr']:>4} | "
                  f"{fmt(per_scene, '8.3g')} {fmt(per_scene and 1e3 * per_scene / int(r['nr']), '7.3g')} | "
                  f"{fmt(e_one, '8.1%')} {fmt(e_full, '8.1%')} {fmt(e_spread, '7.1%')} "
                  f"{fmt(gain, '8.2f') + 'x' if gain else '      n/a'} | "
                  f"{fmt(serial.get((s, 'compact', 1)), '8.3g')} {fmt(serial.get((s, 'compact', full)), '10.3g')}")

            color = PALETTE[i % len(PALETTE)]
            ns = [n for n in all_threads if eff(s, "compact", n) is not None]
            if has_weak and ns:
                name = r["geometry"] if axis == "geometry" else f"{xinfo[0].upper()}={r[xinfo[0]]}"
                series.append((name, ns, [eff(s, "compact", n) for n in ns], color))
            if e_spread is not None:
                markers.append((one_socket, e_spread, color))
            if has_weak and xinfo:
                x = int(r[xinfo[0]])
                for name, v in (("1 socket", e_one), ("all sockets", e_full), ("spread", e_spread)):
                    if v is not None:
                        summary[name].append((x, v))

        if series:
            plot_series(series, f"Parallel efficiency (weak scaling, {axis_name})", res_dir,
                        f"e4_{stem}_efficiency.png", yscale="linear", percent=True,
                        ylim=(0, 1.15),
                        vline=(one_socket, "1 socket") if one_socket and one_socket < full else None,
                        markers=markers, marker_label=f"{one_socket} threads over all sockets")

        if xinfo:
            labels = {"1 socket": f"{one_socket} threads, 1 socket",
                      "spread": f"{one_socket} threads, all sockets",
                      "all sockets": f"{full} threads, all sockets"}
            curves = [(labels[k], [p[0] for p in v], [p[1] for p in v], PALETTE[j])
                      for j, (k, v) in enumerate(summary.items()) if len(v) > 1]
            if curves:
                plot_series(curves, f"Parallel efficiency vs. {xinfo[1]} ({axis_name})", res_dir,
                            f"e4_{stem}_summary.png", yscale="linear", percent=True,
                            ylim=(0, 1.15), xlabel=xinfo[1])
            cost = [(int(settings[s][xinfo[0]]), 1 / thr[(s, "compact", 1)])
                    for s in members if (s, "compact", 1) in thr]
            if len(cost) > 1:
                plot_series([("1 thread", [c[0] for c in cost], [c[1] for c in cost], PALETTE[0])],
                            "Time per scene, single thread [s]", res_dir, f"e4_{stem}_cost.png",
                            xlabel=xinfo[1])

        strong_ns = sorted({n for (st, p, n) in thr if st in members and p == "strongcompact"})
        if strong_ns:
            batch = next(r["batch_size"] for r in runs if r["label"].endswith("_strongcompact"))
            cols = [(n, "strongcompact") for n in dict.fromkeys((one_socket, max(strong_ns))) if n in strong_ns]
            if one_socket and any((st, "strongspread", one_socket) in thr for st in members):
                cols.append((one_socket, "strongspread"))
            print(f"\n--- {axis_name}: strong scaling, batch {batch} (T1 = batch x 1-thread time per scene) ---")
            print(f"{'setting':>26} {'rays':>4} {'s/scene':>8} {'T1 [s]':>9} | " + " | ".join(
                f"{('spread ' if p == 'strongspread' else '') + str(n) + ' thr':>23}" for n, p in cols))
            print(f"{'':>26} {'':>4} {'':>8} {'':>9} | " + " | ".join(f"{'Tn [s]':>7} {'spd':>6} {'eff':>7}" for _ in cols))
            print("-" * (53 + 26 * len(cols)))
            strong_curves = {c: [] for c in cols}
            full_curves, curve_markers = [], []
            for i, s_ in enumerate(members):
                t1 = thr.get((s_, "compact", 1))
                cells = []
                for n, p in cols:
                    e = eff(s_, p, n)
                    tn = batch / thr[(s_, p, n)] if (s_, p, n) in thr else None
                    cells.append(f"{fmt(tn, '7.3g')} {fmt(e and e * n, '6.1f')} {fmt(e, '7.1%')}")
                    if e is not None and xinfo:
                        strong_curves[(n, p)].append((int(settings[s_][xinfo[0]]), e))
                print(f"{s_:>26} {settings[s_]['nr']:>4} {fmt(t1 and 1 / t1, '8.3g')} "
                      f"{fmt(t1 and batch / t1, '9.4g')} | " + " | ".join(cells))
                ns = [n for n in strong_ns if (s_, "strongcompact", n) in thr]
                if len(ns) > 1 and t1:
                    name = (settings[s_]["geometry"] if axis == "geometry"
                            else f"{xinfo[0].upper()}={settings[s_][xinfo[0]]}")
                    full_curves.append((s_, name, [1] + ns, [1.0] + [eff(s_, "strongcompact", n) for n in ns],
                                        PALETTE[i % len(PALETTE)]))
                    if (s_, "strongspread", one_socket) in thr:
                        curve_markers.append((one_socket, eff(s_, "strongspread", one_socket), PALETTE[i % len(PALETTE)]))

            for s_, _, ns, effs, _ in full_curves:
                print(f"  full curve {s_}: " + "  ".join(f"{n}: {n * e:.1f}x/{e:.0%}" for n, e in zip(ns, effs)))
            if full_curves:
                vline = (one_socket, "1 socket") if one_socket and one_socket < full else None
                plot_series([(name, ns, [n * e for n, e in zip(ns, effs)], c) for _, name, ns, effs, c in full_curves],
                            f"Speedup (strong scaling, batch {batch}, {axis_name})", res_dir,
                            f"e4_{stem}_strong_speedup.png", ideal_linear=True, vline=vline,
                            markers=[(x, x * e, c) for x, e, c in curve_markers],
                            marker_label=f"{one_socket} threads over all sockets")
                plot_series([(name, ns, effs, c) for _, name, ns, effs, c in full_curves],
                            f"Parallel efficiency (strong scaling, batch {batch}, {axis_name})", res_dir,
                            f"e4_{stem}_strong_efficiency.png", yscale="linear", percent=True,
                            ylim=(0, 1.15), vline=vline, markers=curve_markers,
                            marker_label=f"{one_socket} threads over all sockets")

            curves = [(f"{n} threads" + (", all sockets" if p == "strongspread" or n > (one_socket or n) else ", 1 socket"),
                       [q[0] for q in v], [q[1] for q in v], PALETTE[j])
                      for j, ((n, p), v) in enumerate(strong_curves.items()) if len(v) > 1]
            if curves:
                plot_series(curves, f"Parallel efficiency (strong scaling, batch {batch}) vs. {xinfo[1]}",
                            res_dir, f"e4_{stem}_strong.png", yscale="linear", percent=True,
                            ylim=(0, 1.15), xlabel=xinfo[1])

    print(f"\nPlots written to {res_dir}/")


if __name__ == "__main__":
    main()
