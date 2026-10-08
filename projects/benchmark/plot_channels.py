#!/usr/bin/env python3
"""Plot a channel list (NU, NGAS, GASES): active gases per channel and which gas
is active in which channel. The channels picked for --nd are marked."""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generate_ctl import pick_channels, read_channels
from plot_style import BLUE, GRID, INK_2, ORANGE, SURFACE, apply_style, save_figure

CELL_EDGE = "#a9c8ee"  # faint blue grid between cells

apply_style()


def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("channels", type=Path, nargs="?", default=here / "configs/channels_alt3.tsv")
    p.add_argument("--nd", type=int, default=32, help="channel count to mark (default: 32, the baseline)")
    p.add_argument("--out", type=Path, default=here / "results/channels")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows = read_channels(args.channels)
    nu = np.array([n for n, _ in rows])
    nactive = np.array([len(g) for _, g in rows])
    picked = {n for n, _ in pick_channels(rows, args.nd)}
    sel = np.array([n in picked for n in nu])
    stem = args.channels.stem

    fig, ax = plt.subplots()
    ax.step(nu, nactive, where="mid", color=BLUE, linewidth=2, label="All channels")
    ax.plot(nu[sel], nactive[sel], linestyle="none", marker="o", markersize=8,
            markerfacecolor=SURFACE, markeredgecolor=ORANGE, markeredgewidth=1.8,
            label=f"ND = {args.nd} channels")
    ax.set_xlabel("Wavenumber [cm$^{-1}$]")
    ax.set_ylabel("Active gases per channel")
    ax.set_ylim(0, nactive.max() + 2)
    ax.set_title(f"{stem}: {len(rows)} channels, {nactive.mean():.1f} active gases on average")
    ax.legend(loc="lower right")
    save_figure(fig, args.out / f"{stem}_active_gases.png")

    # gases sorted by the number of channels they are active in
    gases = sorted({g for _, gs in rows for g in gs},
                   key=lambda g: (-sum(g in gs for _, gs in rows), g))
    active = np.array([[g in gs for _, gs in rows] for g in gases], dtype=int)

    # cell edges halfway between channels (the list is not evenly spaced)
    mid = (nu[:-1] + nu[1:]) / 2
    edges = np.concatenate(([nu[0] - (mid[0] - nu[0])], mid, [nu[-1] + (nu[-1] - mid[-1])]))
    cell = dict(edgecolors=CELL_EDGE, linewidth=0.4, vmin=0, vmax=1)

    fig, (ax_pick, ax) = plt.subplots(
        2, 1, sharex=True, figsize=(6.4, 0.22 * len(gases) + 1.8),
        gridspec_kw=dict(height_ratios=[1, len(gases)], hspace=0.08))
    ax_pick.pcolormesh(edges, [0, 1], sel[None, :].astype(int),
                       cmap=ListedColormap([SURFACE, ORANGE]), **cell)
    ax_pick.tick_params(axis="x", length=0)
    ax_pick.set_yticks([0.5])
    ax_pick.set_yticklabels([f"ND = {args.nd}"], fontsize=8)
    ax_pick.set_title(f"Active gases per channel (blue) and the channels picked for ND = {args.nd} (orange)",
                      color=INK_2, fontsize=10)

    ax.pcolormesh(edges, np.arange(len(gases) + 1), active,
                  cmap=ListedColormap([GRID, BLUE]), **cell)
    ax.set_ylim(len(gases), 0)
    ax.set_yticks(np.arange(len(gases)) + 0.5)
    ax.set_yticklabels([f"{g} ({active[i].sum()})" for i, g in enumerate(gases)], fontsize=8)
    ax.set_xlabel("Wavenumber [cm$^{-1}$]")
    for a in (ax_pick, ax):
        a.tick_params(axis="y", length=0)
        a.grid(False)
        for spine in a.spines.values():
            spine.set_visible(False)
    save_figure(fig, args.out / f"{stem}_gas_matrix.png")

    print(f"Plots written to {args.out}/")


if __name__ == "__main__":
    main()
