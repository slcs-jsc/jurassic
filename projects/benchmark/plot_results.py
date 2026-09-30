#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, PercentFormatter, ScalarFormatter

from plot_style import (BLUE, INK, INK_2, INK_MUTED, ORANGE, SURFACE,
                        apply_style, line_marker_kwargs, save_figure, wrap)

apply_style()


def _tint(color: str, amount: float = 0.82) -> tuple:
    """Blend `color` towards white by `amount` (0 = colour, 1 = white)."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


def _plain_number(v: float, _pos=None) -> str:
    """Tick label: plain number with thousands separator (1, 10, 17,500)."""
    return f"{v:,.10g}"


def _value_text(v: float) -> str:
    """Compact text for a direct value label."""
    return f"{v:,.0f}" if abs(v) >= 100 else f"{v:.1f}"


def _thread_axis(ax, ticks) -> None:
    ax.set_xscale("log", base=2)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    ax.set_xlabel("Threads (physical cores)")


def boxplot(
    threads: np.ndarray,
    value_arr: np.ndarray,
    label: str,
    res_dir: Path,
    filename: str,
    show_outliers: bool = False,
    color: str = BLUE,
) -> None:
    res_dir.mkdir(parents=True, exist_ok=True)
    threads = np.asarray(threads, dtype=float)

    fig, ax = plt.subplots()
    ax.boxplot(
        value_arr, positions=range(len(threads)), widths=0.5,
        patch_artist=True, showfliers=show_outliers,
        boxprops=dict(facecolor=_tint(color), edgecolor=color, linewidth=1.2),
        medianprops=dict(color=color, linewidth=2.0),
        whiskerprops=dict(color=color, linewidth=1.2),
        capprops=dict(color=color, linewidth=1.2),
        flierprops=dict(marker="o", markersize=4, markerfacecolor=color,
                        markeredgecolor=SURFACE),
    )
    ax.set_xticks(range(len(threads)))
    ax.set_xticklabels([str(int(t)) for t in threads])
    ax.set_xlabel("Threads (physical cores)")
    ax.set_ylabel(wrap(label, 40))
    ax.yaxis.set_major_formatter(FuncFormatter(_plain_number))
    save_figure(fig, res_dir / filename)


def plot_scaling(
    threads: np.ndarray,
    values: np.ndarray,
    label: str,
    color: str,
    res_dir: Path,
    filename: str,
    ideal="linear",
    higher_is_better: bool = True,
    stream_ceiling=None,
    smt_points=None,
    yscale: str | None = None,
    percent: bool = False,
    ylim: tuple | None = None,
    title: str | None = None,
) -> None:
    """
    ideal: "linear" | "constant" | None
    higher_is_better: True for throughput-like metrics (bandwidth, FLOP/s);
                      False for cost-like metrics (runtime).
    smt_points: optional list of (threads, value, label) plotted as separate markers
    yscale: "log" or "linear"; default: log for linear-ideal plots (speedup,
            runtime), linear otherwise
    percent: show the y-axis as a percentage (values are fractions)
    ylim: optional (low, high) y-axis limits
    title: optional left-aligned plot title (default: none, the axis label
           carries the quantity)
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    threads = np.asarray(threads, dtype=float)
    values = np.asarray(values, dtype=float)
    yscale = yscale or ("log" if ideal == "linear" else "linear")

    fig, ax = plt.subplots()

    # Reference lines first (behind the data), muted and dashed
    if ideal == "linear":
        ideal_vals = (
            values[0] * (threads / threads[0]) if higher_is_better
            else values[0] * (threads[0] / threads)
        )
        ax.plot(threads, ideal_vals, linestyle=(0, (5, 3)), color=INK_MUTED,
                linewidth=1.3, label="Ideal scaling", zorder=2)
    elif ideal == "constant":
        ax.axhline(values[0], linestyle=(0, (5, 3)), color=INK_MUTED,
                   linewidth=1.3, label="Ideal", zorder=2)

    if stream_ceiling is not None:
        ax.axhline(stream_ceiling, linestyle=(0, (1, 2)), color=INK_2,
                   linewidth=1.3, zorder=2,
                   label=f"STREAM ceiling ({stream_ceiling:.0f})")

    ax.plot(threads, values, label="Measured", zorder=3,
            **line_marker_kwargs(color))

    # Label the last point with its value (selective direct label)
    last_txt = f"{values[-1]:.0%}" if percent else _value_text(values[-1])
    ax.annotate(last_txt, (threads[-1], values[-1]), textcoords="offset points",
                xytext=(8, 0), ha="left", va="center", fontsize=9, color=INK)

    if smt_points:
        for t, v, l in smt_points:
            ax.plot([t], [v], linestyle="none", marker="s", markersize=7,
                    markerfacecolor=ORANGE, markeredgecolor=SURFACE,
                    zorder=4, label=l)

    ax.set_yscale(yscale)
    _thread_axis(ax, list(threads) + ([smt_points[0][0]] if smt_points else []))
    ax.set_ylabel(wrap(label, 40))
    if percent:
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(_plain_number))
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif yscale == "linear":
        ax.set_ylim(bottom=0)
    ax.margins(x=0.06)
    if title:
        ax.set_title(title)
    ax.legend(loc="best")
    save_figure(fig, res_dir / filename)
