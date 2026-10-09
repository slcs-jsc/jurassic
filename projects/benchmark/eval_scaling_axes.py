"""Evaluate scaling_axes.sh runs.

Efficiency = throughput(n) / (n * throughput(1)), T1 = batch x 1-thread time per scene.
"formod batch" times one formod_batch call; "whole application" adds the
READ_*, FORMOD_REFERENCE, WRITE_OBS and FINALIZE timers.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from plot_results import plot_panels, plot_series
from plot_style import PALETTE

AXIS_X = {"geometry": None, "channels": ("nd", "Channels (ND)"), "gases": ("ng", "Emitters (NG)")}


MEAN_RE = re.compile(r"RUNTIME:.*?\bmean=\s*([\d.eE+-]+)\s*s")
TIMER_RE = re.compile(r"^(TIMER_\w+)\s*=\s*([\d.eE+-]+)\s*s", re.MULTILINE)
TABLE_RE = re.compile(r"Read emissivity table: tbl_")
TABLES_READ_RE = re.compile(r"^tables_read=(\d+)", re.MULTILINE)  # stripped logs
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
        text = txt.read_text()
        mean = MEAN_RE.search(text)
        if mean is None:
            print(f"WARNING: no RUNTIME line in {txt.name} (run failed?), skipping.")
            continue
        timers = {k: float(v) for k, v in TIMER_RE.findall(text)}
        runs.append({
            "label": m.group("label"),
            "threads": int(m.group("threads")),
            "batch_size": int(m.group("batch")),
            "rep": int(m.group("rep")),
            "mean_s": float(mean.group(1)),
            "tables": (int(m_t.group(1)) if (m_t := TABLES_READ_RE.search(text))
                       else len(TABLE_RE.findall(text))),
            "serial_s": sum(v for k, v in timers.items()
                            if k.startswith("TIMER_READ_") or k in SERIAL_TIMERS)
                        if timers else None,
        })
    return runs


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

    # single run or array run (one subdirectory per case)
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

    tables: dict[str, int] = {}
    # (setting, placement, threads) -> throughput, serial time
    thr: dict[tuple, list] = {}
    serial: dict[tuple, list] = {}
    # batch sweep: (setting, threads, batch)
    bthr: dict[tuple, list] = {}
    bserial: dict[tuple, list] = {}
    for r in runs:
        setting, _, placement = r["label"].rpartition("_")
        if setting not in settings:
            continue
        tables[setting] = max(tables.get(setting, 0), r["tables"])
        if placement == "batch":
            key = (setting, r["threads"], r["batch_size"])
            bthr.setdefault(key, []).append(r["batch_size"] / r["mean_s"])
            if r["serial_s"] is not None:
                bserial.setdefault(key, []).append(r["serial_s"])
            continue
        key = (setting, placement, r["threads"])
        thr.setdefault(key, []).append(r["batch_size"] / r["mean_s"])
        if r["serial_s"] is not None:
            serial.setdefault(key, []).append(r["serial_s"])
    if not thr:
        print(f"No usable runs under {args.run_dir}.")
        sys.exit(1)
    thr = {k: float(np.median(v)) for k, v in thr.items()}
    # min: first run of a setting reads the tables from a cold cache
    serial = {k: float(np.min(v)) for k, v in serial.items()}
    bthr = {k: float(np.median(v)) for k, v in bthr.items()}
    bserial = {k: float(np.min(v)) for k, v in bserial.items()}

    def eff(setting, placement, n):
        base = thr.get((setting, "compact", 1))
        t = thr.get((setting, placement, n))
        return t / (n * base) if base and t else None

    def app_eff(setting, n, batch):
        t1, tn = thr.get((setting, "compact", 1)), thr.get((setting, "strongcompact", n))
        s1, sn = serial.get((setting, "compact", 1)), serial.get((setting, "strongcompact", n))
        if None in (t1, tn, s1, sn):
            return None
        return (s1 + batch / t1) / (n * (sn + batch / tn))

    def fmt(v, spec):
        return format(v, spec) if v is not None else format("n/a", ">" + re.match(r"\d*", spec).group())

    all_threads = sorted({n for (_, p, n) in thr if p == "compact"})
    has_weak = len(all_threads) > 1
    full = max(n for (_, p, n) in thr if p in ("compact", "strongcompact"))
    one_socket = per_socket if any(n == per_socket for (_, p, n) in thr if p != "spread") else None

    # active (channel, gas) pairs: expected from the channel list, found = tables read
    print("\n=== 1-thread cost ===")
    print(f"{'setting':>22} {'ND':>4} {'NG':>3} {'rays':>4} {'pairs':>6} {'found':>6} | "
          f"{'s/scene':>8} {'ms/pair':>8}")
    print("-" * 72)
    for s, r in settings.items():
        t1 = thr.get((s, "compact", 1))
        pairs = int(r["active_pairs"]) if r.get("active_pairs") else None
        found = tables.get(s)
        print(f"{s:>22} {r['nd']:>4} {r['ng']:>3} {r['nr']:>4} {fmt(pairs, '6d')} {fmt(found, '6d')} | "
              f"{fmt(t1 and 1 / t1, '8.3g')} {fmt(t1 and pairs and 1e3 / t1 / pairs, '8.3g')}"
              + ("  <- tables missing" if pairs and found is not None and found < pairs else ""))

    # channels/gases plots per geometry, drawn side by side after the loop
    panels: dict[tuple, list] = {}

    def add_panel(axis, ref, plot, series, batch=None):
        kind, _, geom = axis.partition("@")
        fixed = f"NG={ref['ng']}" if kind == "channels" else f"ND={ref['nd']}"
        panels.setdefault((kind, plot), []).append((geom or "all", fixed, batch, series))
    # axis tags: geometry, channels@<geom>, gases@<geom>
    tags = list(dict.fromkeys(t for r in settings.values() for t in r["axes"].split(",")))
    tags.sort(key=lambda t: list(AXIS_X).index(t.split("@")[0]))
    for axis in tags:
        xinfo = AXIS_X[axis.split("@")[0]]
        stem = axis.replace("@", "_")
        members = [s for s, r in settings.items() if axis in r["axes"].split(",")]
        if not members:
            continue
        ref = settings[members[0]]
        axis_name = {"geometry": f"Geometries (ND={ref['nd']}, NG={ref['ng']})",
                     "channels": f"Channel count, {ref['geometry']} (NG={ref['ng']})",
                     "gases": f"Gas sets, {ref['geometry']} (ND={ref['nd']})"}[axis.split("@")[0]]
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
            plot_series(series, "Parallel efficiency", res_dir,
                        f"e4_{stem}_efficiency.png", yscale="linear", percent=True,
                        title=f"{axis_name} · weak scaling",
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
                plot_series(curves, "Parallel efficiency", res_dir,
                            f"e4_{stem}_summary.png", yscale="linear", percent=True,
                            title=f"{axis_name} · weak scaling",
                            ylim=(0, 1.15), xlabel=xinfo[1])
            cost = [(int(settings[s][xinfo[0]]), 1 / thr[(s, "compact", 1)])
                    for s in members if (s, "compact", 1) in thr]
            if len(cost) > 1:
                add_panel(axis, ref, "cost",
                          [("1 thread", [c[0] for c in cost], [c[1] for c in cost], PALETTE[0])])

        strong_ns = sorted({n for (st, p, n) in thr if st in members and p == "strongcompact"})
        if strong_ns:
            batch = next(r["batch_size"] for r in runs if r["label"].endswith("_strongcompact"))
            cols = [(n, "strongcompact") for n in dict.fromkeys((one_socket, max(strong_ns))) if n in strong_ns]
            if one_socket and any((st, "strongspread", one_socket) in thr for st in members):
                cols.append((one_socket, "strongspread"))
            print(f"\n--- {axis_name}: strong scaling, formod batch, batch {batch} (T1 = batch x 1-thread time per scene) ---")
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
                print(f"  formod batch {s_}: " + "  ".join(f"{n}: {n * e:.1f}x/{e:.0%}" for n, e in zip(ns, effs)))
            app_curves = []
            for s_, name, ns, _, c in full_curves:
                ae = [app_eff(s_, n, batch) if n > 1 else 1.0 for n in ns]
                if None not in ae:
                    app_curves.append((name, ns, ae, c))
                    print(f"  application  {s_}: " + "  ".join(f"{n}: {n * e:.1f}x/{e:.0%}" for n, e in zip(ns, ae)))
            per_geom = "@" in axis
            if full_curves:
                vline = (one_socket, "1 socket") if one_socket and one_socket < full else None
                spd = [(name, ns, [n * e for n, e in zip(ns, effs)], c) for _, name, ns, effs, c in full_curves]
                effc = [(name, ns, effs, c) for _, name, ns, effs, c in full_curves]
                if per_geom:
                    add_panel(axis, ref, "strong_speedup", spd, batch)
                    add_panel(axis, ref, "strong_efficiency", effc, batch)
                else:
                    plot_series(spd, "Speedup", res_dir, f"e4_{stem}_strong_speedup.png",
                                ideal_linear=True, vline=vline,
                                title=f"{axis_name} · formod batch, strong scaling, batch {batch}",
                                markers=[(x, x * e, c) for x, e, c in curve_markers],
                                marker_label=f"{one_socket} threads over all sockets")
                    plot_series(effc, "Parallel efficiency", res_dir, f"e4_{stem}_strong_efficiency.png",
                                yscale="linear", percent=True, ylim=(0, 1.15), vline=vline,
                                title=f"{axis_name} · formod batch, strong scaling, batch {batch}",
                                markers=curve_markers, marker_label=f"{one_socket} threads over all sockets")
            if app_curves:
                app_spd = [(name, ns, [n * e for n, e in zip(ns, effs)], c) for name, ns, effs, c in app_curves]
                if per_geom:
                    add_panel(axis, ref, "app_speedup", app_spd, batch)
                    add_panel(axis, ref, "app_efficiency", app_curves, batch)
                else:
                    plot_series(app_spd, "Speedup", res_dir, f"e4_{stem}_app_speedup.png",
                                ideal_linear=True, vline=vline,
                                title=f"{axis_name} · whole application, strong scaling, batch {batch}")
                    plot_series(app_curves, "Parallel efficiency", res_dir, f"e4_{stem}_app_efficiency.png",
                                yscale="linear", percent=True, ylim=(0, 1.15), vline=vline,
                                title=f"{axis_name} · whole application, strong scaling, batch {batch}")

            curves = [(f"{n} threads" + (", all sockets" if p == "strongspread" or n > (one_socket or n) else ", 1 socket"),
                       [q[0] for q in v], [q[1] for q in v], PALETTE[j])
                      for j, ((n, p), v) in enumerate(strong_curves.items()) if len(v) > 1]
            if curves:
                add_panel(axis, ref, "strong", curves, batch)

    names = {"channels": "Channel count", "gases": "Gas sets"}
    specs = {
        "strong_speedup": ("Speedup", dict(yscale="log", ideal_linear=True), "formod batch, strong scaling"),
        "strong_efficiency": ("Parallel efficiency", dict(percent=True, ylim=(0, 1.15)), "formod batch, strong scaling"),
        "app_speedup": ("Speedup", dict(yscale="log", ideal_linear=True), "whole application, strong scaling"),
        "app_efficiency": ("Parallel efficiency", dict(percent=True, ylim=(0, 1.15)), "whole application, strong scaling"),
        "strong": ("Parallel efficiency", dict(percent=True, ylim=(0, 1.15)), "formod batch, strong scaling"),
        "cost": ("Time per scene [s]", dict(yscale="log"), "1 thread"),
    }
    for (kind, plot), groups in panels.items():
        label, opts, what = specs[plot]
        if plot in ("strong", "cost"):
            opts = dict(opts, xlabel=AXIS_X[kind][1])
        groups.sort(key=lambda g: g[0])
        _, fixed, batch, _ = groups[0]
        plot_panels([(geom, series) for geom, _, _, series in groups], label, res_dir,
                    f"e4_{kind}_{plot}.png", **opts,
                    title=f"{names[kind]} ({fixed}) · {what}" + (f", batch {batch}" if batch else ""))

    checks = [(st, b) for (st, p, n) in thr if p == "t1full" and (st, "compact", 1) in thr
              for b in {r["batch_size"] for r in runs if r["label"] == f"{st}_t1full"}]
    if checks:
        print("\n=== T1 check: 1 thread on the full batch vs. batch x 1-thread time per scene ===")
        for st, b in sorted(checks):
            est, meas = b / thr[(st, "compact", 1)], b / thr[(st, "t1full", 1)]
            print(f"{st:>26} batch {b}: extrapolated {est:.4g} s, measured {meas:.4g} s ({est / meas - 1:+.1%})")

    if bthr:
        batch_sweep(settings, thr, serial, bthr, bserial, res_dir)

    print(f"\nPlots written to {res_dir}/")


def batch_sweep(settings, thr, serial, bthr, bserial, res_dir):
    """Efficiency vs batch size at a fixed thread count."""
    for n in sorted({k[1] for k in bthr}):
        print(f"\n=== Batch size sweep, {n} threads (T1 = batch x 1-thread time per scene) ===")
        print(f"{'setting':>26} {'batch':>6} {'per thr':>7} | {'Tn [s]':>8} {'spd':>6} {'eff':>7} | "
              f"{'app spd':>7} {'app eff':>7}")
        print("-" * 86)
        formod = []
        for i, s in enumerate(sorted({k[0] for k in bthr if k[1] == n})):
            t1, s1 = thr.get((s, "compact", 1)), serial.get((s, "compact", 1))
            if not t1:
                continue
            bs = sorted(b for (st, nn, b) in bthr if st == s and nn == n)
            effs, aeffs = [], []
            for b in bs:
                tn = b / bthr[(s, n, b)]
                e = bthr[(s, n, b)] / (n * t1)
                sn = bserial.get((s, n, b))
                ae = (s1 + b / t1) / (n * (sn + tn)) if s1 is not None and sn is not None else None
                effs.append(e)
                aeffs.append(ae)
                print(f"{s:>26} {b:>6} {b // n:>7} | {tn:8.3g} {n * e:6.1f} {e:7.1%} | "
                      + (f"{n * ae:7.1f} {ae:7.1%}" if ae is not None else f"{'n/a':>7} {'n/a':>7}"))
            name = settings[s]["geometry"]
            color = PALETTE[i % len(PALETTE)]
            formod.append((name, bs, effs, color))
        if formod and max(len(c[1]) for c in formod) > 1:
            plot_series(formod, "Parallel efficiency", res_dir, f"e4_batch_t{n}.png",
                        yscale="linear", percent=True, ylim=(0, 1.15), xlabel="Batch size (scenes)",
                        title=f"Batch size · formod batch, {n} threads")


if __name__ == "__main__":
    main()
