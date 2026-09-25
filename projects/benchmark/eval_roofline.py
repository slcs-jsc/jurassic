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
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt

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

    ridge_x = peak_flops / stream_bw
    bench_threads = int(ceilings.get("threads", 0)) or "unknown"

    print("\n=== Roofline ceilings ===")
    print(f"  compute ceiling:   {peak_flops:.0f} MFLOP/s")
    print(f"  bandwidth ceiling: {stream_bw:.0f} MBytes/s")
    print(f"  ridge point:       {ridge_x:.3f} FLOP/Byte")
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

    points = []  # (label, intensity, mflops)

    header = (f"{'case':>20} {'threads':>8} | "
              f"{'FLOP/Byte':>12} | {'DP MFLOP/s':>12} | {'cv (AI/perf)':>14}")
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
            print(f"{label:>20} {threads:>8} | could not derive arithmetic intensity, skipping")
            continue
 
        cv_str = (f"{intensity_cv:.1%}/{mflops_cv:.1%}"
                  if intensity_cv is not None and mflops_cv is not None else "n/a")
        print(f"{label:>20} {threads:>8} | {intensity_med:>12.4g} | {mflops_med:>12.4g} | {cv_str:>14}")
        points.append((label, intensity_med, mflops_med))

    if not points:
        print("\nNo usable roofline points found.")
        sys.exit(1)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
 
    x_lo = min(p[1] for p in points) / 3
    x_hi = max(ridge_x, max(p[1] for p in points)) * 3
 
    roof_x = [x_lo, ridge_x, x_hi]
    roof_y = [stream_bw * x_lo, peak_flops, peak_flops]
    ax.plot(roof_x, roof_y, "--", color="#898781", linewidth=1.5,
            label=f"Roofline ({peak_flops/1e3:.0f} GFLOP/s, {stream_bw/1e3:.0f} GB/s)",
            zorder=2)
 
    ax.axvline(ridge_x, color="#898781", linewidth=0.6, linestyle=":", zorder=1)
    ax.text(ridge_x * 1.05, peak_flops * 0.6,
            f"ridge\n{ridge_x:.2f} F/B", fontsize=7, color="#898781")
 
    palette = ["#2a78d6", "#c76ce0", "#efb239", "#3ab9dc", "#89dd29", "#e05c2a"]
    for i, (label, intensity, mflops) in enumerate(sorted(points)):
        color = palette[i % len(palette)]
        ax.scatter([intensity], [mflops], color=color, s=60, zorder=4, label=label)
        bound = "mem" if intensity < ridge_x else "compute"
        ax.annotate(f"{label}\n({bound})", (intensity, mflops),
                    textcoords="offset points", xytext=(6, 4),
                    fontsize=7, color=color)
 
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Arithmetic intensity [FLOP/Byte]  (derived: ΣDP flops / ΣMemory bandwidth)")
    ax.set_ylabel("Performance [MFLOP/s]")
    ax.set_title(f"JURASSIC forward model — roofline  ({bench_threads} threads)")
    ax.legend(fontsize=8)
    fig.tight_layout()
 
    outpath = res_dir / "e1_roofline.png"
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {outpath}")

if __name__ == "__main__":
    main()