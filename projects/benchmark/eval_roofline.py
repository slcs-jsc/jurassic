"""
Roofline evaluation for JURASSIC's forward model.

Reads LIKWID FLOPS_DP + MEM_DP marker-region output produced by
run_e1_roofline.sh, loads the compute and bandwidth ceilings from
the ceilings.txt file written by likwid-bench in that same run, and
plots each case as a single point on a log-log roofline diagram.

Usage:
    python eval_roofline.py <run_dir>

Override ceilings if needed (skips ceilings.txt):
    python eval_roofline.py <run_dir> --peak-flops 546000 --stream-bw 110000

Optional extra roofs: add lines like L2_bw_mbytes=900000 (bandwidth) or
scalar_flops_mflops=80000 (compute) to ceilings.txt, or pass them on the command line
(repeatable), e.g.:
    python eval_roofline.py <run_dir> --bw-ceiling L2=900000 --bw-ceiling L3=400000 \\
        --compute-ceiling scalar=80000
"""
import argparse
import sys
from pathlib import Path

import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))

from likwid_parsing import parse_run_dir, collect, collect_intensity

def load_ceilings(run_dir: Path) -> dict:
    """Parse ceilings.txt written by likwid-bench in the run script."""
    path = run_dir / "ceilings.txt"
    if not path.is_file():
        return {}
    result = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            result[k.strip()] = v.strip()
    return result

def parse_named_value(text: str) -> tuple:
    """Parse NAME=VALUE command-line arguments."""
    name, sep, value = text.partition("=")
    try:
        if not sep or not name.strip():
            raise ValueError
        return name.strip(), float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected NAME=NUMBER, got '{text}'")

def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("run_dir", type=Path)
    parser.add_argument(
        "--region", default="formod",
        help="likwid-marker region name to read (default: formod)"
    )
    parser.add_argument(
        "--warmup", type=int, default=1,
        help="number of leading repetitions to discard (default: 1)"
    )
    parser.add_argument(
        "--peak-flops", type=float, default=None,
        help="override compute ceiling [MFLOP/s] (default: read from ceilings.txt)"
    )
    parser.add_argument(
        "--stream-bw", type=float, default=None,
        help="override bandwidth ceiling [MBytes/s] (default: read from ceilings.txt)"
    )
    parser.add_argument(
        "--bw-ceiling", action="append", default=[], metavar="NAME=MBYTES_S",
        type=parse_named_value,
        help="additional bandwidth roof, e.g. L2=800000 (repeatable). Also read from "
             "ceilings.txt keys '<name>_bw_mbytes' (except stream_bw_mbytes)"
    )
    parser.add_argument(
        "--compute-ceiling", action="append", default=[], metavar="NAME=MFLOPS",
        type=parse_named_value,
        help="additional (lower) compute ceiling, e.g. scalar=80000 (repeatable). Also read "
             "from ceilings.txt keys '<name>_flops_mflops' (except peak_flops_mflops)"
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="output directory (default: <run_dir>/plots)"
    )
    args = parser.parse_args()

    res_dir = args.out or (args.run_dir / "plots")
    res_dir.mkdir(parents=True, exist_ok=True)

    # --- Load ceilings ---
    ceilings = load_ceilings(args.run_dir)

    peak_flops = args.peak_flops or (
        float(ceilings["peak_flops_mflops"]) if "peak_flops_mflops" in ceilings else None
    )
    stream_bw = args.stream_bw or (
        float(ceilings["stream_bw_mbytes"]) if "stream_bw_mbytes" in ceilings else None
    )

    if peak_flops is None or stream_bw is None:
        print(
            "ERROR: no ceilings.txt found in run_dir and --peak-flops/--stream-bw not given.\n"
            "       Either rerun with the updated run script (which calls likwid-bench),\n"
            "       or pass --peak-flops and --stream-bw manually.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Optional extra roofs: file first, CLI overrides
    def extra_ceilings(suffix: str, main_key: str) -> dict:
        found = {}
        for k, v in ceilings.items():
            if k.endswith(suffix) and k != main_key:
                try:
                    found[k[:-len(suffix)]] = float(v)
                except ValueError:
                    pass
        return found

    bw_roofs = extra_ceilings("_bw_mbytes", "stream_bw_mbytes")
    bw_roofs.update(dict(args.bw_ceiling))
    compute_roofs = extra_ceilings("_flops_mflops", "peak_flops_mflops")
    compute_roofs.update(dict(args.compute_ceiling))

    ridge_x = peak_flops / stream_bw
    bench_threads = int(ceilings.get("threads", 0)) or "unknown"

    print("\n=== Roofline ceilings ===")
    print(f"  compute ceiling:   {peak_flops:.0f} MFLOP/s")
    print(f"  bandwidth ceiling: {stream_bw:.0f} MBytes/s")
    print(f"  ridge point:       {ridge_x:.3f} FLOP/Byte")
    for name, bw in bw_roofs.items():
        print(f"  {name + ' bandwidth:':<18} {bw:.0f} MBytes/s (ridge {peak_flops / bw:.3f} FLOP/Byte)")
    for name, cp in compute_roofs.items():
        print(f"  {name + ' compute:':<18} {cp:.0f} MFLOP/s")
    print(f"  measured at:       {bench_threads} threads\n")

    # --- Parse LIKWID profiling data ---
    configs = parse_run_dir(args.run_dir)
    if not configs:
        print(f"No parsed configs found under {args.run_dir}/out/.")
        sys.exit(1)

    # Group by (case, threads) -- no batch dimension
    groups: dict[tuple, list] = {}
    for c in configs:
        key = (c["label"], c["threads"])
        groups.setdefault(key, []).append(c)

    points = []  # (label, threads, intensity, mflops, intensity_cv, mflops_cv)

    header = (f"{'case':>20} {'threads':>8} | "
              f"{'FLOP/Byte':>12} | {'DP MFLOP/s':>12} | {'% of roof':>9} | {'cv (AI/perf)':>14}")
    print(header)
    print("-" * len(header))

    for (label, threads), entries in sorted(groups.items()):
        flops_entries = sorted(
            [e for e in entries if e["group"] == "FLOPS_DP"], key=lambda e: e["rep"]
        )
        mem_entries = sorted(
            [e for e in entries if e["group"] == "MEM_DP"], key=lambda e: e["rep"]
        )

        if not flops_entries or not mem_entries:
            print(f"{label:>20} {threads:>8} | missing FLOPS_DP and/or MEM_DP, skipping")
            continue
        if len(flops_entries) < args.warmup or len(mem_entries) <= args.warmup:
            print(f"{label:>20} {threads:>8} | insufficient reps after warmup discard, skipping")
            continue

        _, mflops_med, _, mflops_cv, _, _ = collect(
            flops_entries, args.region, "DP [MFLOP/s]", warmup=args.warmup
        )

        if mflops_med is None:
            print(f"{label:>20} {threads:>8} | missing DP [MFLOP/s], skipping")
            continue

        _, intensity_med, _, intensity_cv, _, _ = collect_intensity(
            mem_entries, args.region, warmup=args.warmup
        )
        if intensity_med is None:
            print(f"{label:>20} {threads:>8} | could not derive operational intensity, skipping")
            continue
 
        cv_str = (f"{intensity_cv:.1%}/{mflops_cv:.1%}"
                  if intensity_cv is not None and mflops_cv is not None else "n/a")
        roof_here = min(peak_flops, stream_bw * intensity_med)
        print(f"{label:>20} {threads:>8} | {intensity_med:>12.4g} | {mflops_med:>12.4g} | "
              f"{mflops_med / roof_here:>9.1%} | {cv_str:>14}")
        if bench_threads != "unknown" and threads != bench_threads:
            print(f"{'':>20} {'':>8} | WARNING: ceilings measured at {bench_threads} threads")
        points.append((label, threads, intensity_med, mflops_med, intensity_cv, mflops_cv))

    if not points:
        print("\nNo usable roofline points found.")
        sys.exit(1)

    # --- Plot (axes in GFLOP/s and GB/s) ---
    G = 1e3
    peak = peak_flops / G
    bw_main = stream_bw / G
    all_bw = {"DRAM": bw_main, **{n: b / G for n, b in bw_roofs.items()}}
    ridge = peak / bw_main
    multi_threads = len({p[1] for p in points}) > 1

    xs = [p[2] for p in points]
    ys = [p[3] / G for p in points]
    x_lo = min(min(xs), peak / max(all_bw.values())) / 2
    x_hi = max(ridge, max(xs)) * 3
    y_lo = min(ys) / 3
    y_hi = peak * 2

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)

    def slope_deg(bw: float) -> float:
        """On-screen angle of a bandwidth roof (slope 1 in log-log data space)."""
        fig.canvas.draw()  # finalize axes geometry
        box = ax.get_window_extent()
        px_per_x = box.width / math.log10(x_hi / x_lo)
        px_per_y = box.height / math.log10(y_hi / y_lo)
        return math.degrees(math.atan2(px_per_y, px_per_x))

    # bandwidth roofs (DRAM solid, cache levels dashed); each ends at the compute peak
    for name, bw in all_bw.items():
        knee = peak / bw
        main = name == "DRAM"
        color = "#4a4945" if main else "#898781"
        ax.plot([x_lo, knee], [bw * x_lo, peak], "-" if main else "--",
                color=color, linewidth=1.8 if main else 1.2, zorder=2)
        # label along the slope, kept inside the visible area
        xt = min(max(x_lo * 1.3, y_lo * 1.5 / bw), knee / 1.5)
        ax.text(xt, bw * xt * 1.12, f"{name} {bw:.0f} GB/s", fontsize=8, color=color,
                rotation=slope_deg(bw), rotation_mode="anchor", ha="left", va="bottom")
    # compute roof and optional lower compute ceilings
    ax.plot([ridge, x_hi], [peak, peak], "-", color="#4a4945", linewidth=1.8, zorder=2)
    ax.text(x_hi / 1.1, peak * 1.06, f"peak {peak:.0f} GFLOP/s", fontsize=8,
            color="#4a4945", ha="right", va="bottom")
    for name, cp in compute_roofs.items():
        cp /= G
        ax.plot([cp / bw_main, x_hi], [cp, cp], "--", color="#898781", linewidth=1.2, zorder=2)
        ax.text(x_hi / 1.1, cp * 1.06, f"{name} {cp:.0f} GFLOP/s", fontsize=8,
                color="#898781", ha="right", va="bottom")

    ax.axvline(ridge, color="#898781", linewidth=0.6, linestyle=":", zorder=1)
    ax.text(ridge * 1.04, y_lo * 1.15, f"ridge {ridge:.2f} FLOP/B", fontsize=7,
            color="#898781", rotation=90, va="bottom")

    palette = ["#2a78d6", "#c76ce0", "#efb239", "#3ab9dc", "#89dd29", "#e05c2a"]
    # labels fan out (left / below / right) in order of intensity to avoid overlaps
    offsets = [(-14, -30, "right"), (0, -48, "center"), (14, 10, "left")]
    by_intensity = sorted(range(len(points)), key=lambda k: points[k][2])
    label_slot = {k: r for r, k in enumerate(by_intensity)}
    for i, (label, threads, intensity, mflops, icv, pcv) in enumerate(points):
        color = palette[i % len(palette)]
        perf = mflops / G
        ax.errorbar([intensity], [perf], fmt="o", color=color, markersize=7, zorder=4,
                    xerr=[[intensity * icv]] if icv else None,
                    yerr=[[perf * pcv]] if pcv else None,
                    capsize=2, elinewidth=1)
        # thin guide up to the DRAM/compute roof directly above the point
        roof_here = min(peak, bw_main * intensity)
        ax.plot([intensity, intensity], [perf, roof_here], ":", color=color,
                linewidth=0.8, zorder=1)
        name = f"{label} [{threads}T]" if multi_threads else label
        dx, dy, ha = offsets[label_slot[i] % len(offsets)]
        ax.annotate(f"{name}\n{perf / roof_here:.0%} of roof", (intensity, perf),
                    textcoords="offset points", xytext=(dx, dy), fontsize=8, color=color,
                    ha=ha, arrowprops=dict(arrowstyle="-", color=color, linewidth=0.6))

    # readable ticks: plain numbers, with 2 and 5 labeled between the decades
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
        axis.set_minor_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.tick_params(axis="both", which="minor", labelsize=7)
    ax.grid(True, which="major", color="#e4e2dc", linewidth=0.6, zorder=0)
    ax.set_xlabel("Operational intensity [FLOP/Byte]")
    ax.set_ylabel("Performance [GFLOP/s]")
    suffix = f"  ({bench_threads} threads)" if bench_threads != "unknown" else ""
    ax.set_title(f"JURASSIC forward model — Roofline{suffix}")
    fig.tight_layout()

    outpath = res_dir / "e1_roofline.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {outpath}")

if __name__ == "__main__":
    main()