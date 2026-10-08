#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, LogFormatter, PercentFormatter, ScalarFormatter

from plot_style import (INK_2, INK_MUTED, SURFACE, apply_style,
                        line_marker_kwargs, save_figure, wrap)

apply_style()


def _plain_number(v: float, _pos=None) -> str:
    """Tick label: plain number with thousands separator (1, 10, 17,500)."""
    return f"{v:,.10g}"


def _value_text(v: float) -> str:
    """Compact text for a direct value label."""
    return f"{v:,.0f}" if abs(v) >= 100 else f"{v:.1f}"


def _thread_axis(ax, ticks, xlabel: str = "Threads (physical cores)") -> None:
    ax.set_xscale("log", base=2)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda *_: ""))
    ax.set_xlabel(xlabel)


def plot_series(
    series: list,
    label: str,
    res_dir: Path,
    filename: str,
    stream_ceiling=None,
    yscale: str = "log",
    xlabel: str = "Threads (physical cores)",
    percent: bool = False,
    ylim: tuple | None = None,
    vline: tuple | None = None,
    markers: list | None = None,
    marker_label: str | None = None,
    ideal_linear: bool = False,
    title: str | None = None,
) -> None:
    """
    Several measured curves on one axis, e.g. one per geometry.
    series: list of (name, x, y, color) tuples
    vline: optional (x, text) vertical reference, e.g. the socket boundary
    markers: optional (x, y, color) hollow squares, labelled once as marker_label
    ideal_linear: dashed y = x reference (ideal speedup)
    """
    res_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    xs_all, labelled = [], []
    for name, x, y, color in series:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        xs_all.extend(x)
        ax.plot(x, y, label=name, zorder=3, **line_marker_kwargs(color))
        labelled.append((x[-1], y[-1], color))

    if ideal_linear and xs_all:
        lo, hi = min(xs_all), max(xs_all)
        ax.plot([lo, hi], [lo, hi], linestyle=(0, (5, 3)), color=INK_MUTED,
                linewidth=1.3, label="Ideal", zorder=2)

    for i, (mx, my, color) in enumerate(markers or []):
        ax.plot([mx], [my], linestyle="none", marker="s", markersize=8,
                markerfacecolor=SURFACE, markeredgecolor=color, markeredgewidth=1.8,
                zorder=4, label=marker_label if i == 0 else None)

    if stream_ceiling is not None:
        ax.axhline(stream_ceiling, linestyle=(0, (1, 2)), color=INK_2,
                   linewidth=1.3, zorder=2,
                   label=f"STREAM ceiling ({stream_ceiling:.0f})")

    if vline is not None:
        ax.axvline(vline[0], linestyle=(0, (1, 2)), color=INK_MUTED, linewidth=1.0, zorder=1)
        ax.annotate(vline[1], (vline[0], 0.02), xycoords=("data", "axes fraction"),
                    xytext=(-4, 0), textcoords="offset points", rotation=90,
                    ha="right", va="bottom", fontsize=8, color=INK_2)

    ax.set_yscale(yscale)
    _thread_axis(ax, sorted(set(xs_all)), xlabel=xlabel)
    ax.set_ylabel(wrap(label, 40))
    if percent:
        ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    else:
        ax.yaxis.set_major_formatter(FuncFormatter(_plain_number))
    if yscale == "log":
        ax.yaxis.set_minor_formatter(LogFormatter(base=10))
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif yscale != "log":
        ax.set_ylim(bottom=0)
    ax.margins(x=0.06)
    if title:
        ax.set_title(title)
    _end_labels(ax, labelled, percent)
    ax.legend(loc="best")
    save_figure(fig, res_dir / filename)


def _end_labels(ax, points, percent, gap=11.0):
    """Value labels at the line ends, pushed apart vertically so they don't overlap."""
    fig = ax.figure
    fig.canvas.draw()
    to_pt = 72.0 / fig.dpi
    by_x = {}
    for x, y, color in points:
        by_x.setdefault(x, []).append((ax.transData.transform((x, y))[1] * to_pt, y, color))
    for x, group in by_x.items():
        group.sort()
        placed = []
        for yp, _, _ in group:
            placed.append(max(yp, placed[-1] + gap) if placed else yp)
        shift = (sum(p - g[0] for p, g in zip(placed, group))) / len(group)
        for p, (yp, y, color) in zip(placed, group):
            ax.annotate(f"{y:.0%}" if percent else _value_text(y), (x, y),
                        textcoords="offset points", xytext=(8, p - shift - yp),
                        ha="left", va="center", fontsize=9, color=color)


def plot_panels(
    panels: list,
    label: str,
    res_dir: Path,
    filename: str,
    yscale: str = "linear",
    percent: bool = False,
    ylim: tuple | None = None,
    ideal_linear: bool = False,
    title: str | None = None,
    xlabel: str = "Threads (physical cores)",
) -> None:
    """
    One panel per group (e.g. per geometry) with a shared y-axis and one legend.
    panels: list of (panel_title, series), series as in plot_series
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(panels), sharey=True,
                             figsize=(4.2 * len(panels) + 1.0, 4.2))
    axes = np.atleast_1d(axes)
    for ax, (panel_title, series) in zip(axes, panels):
        xs_all = []
        for name, x, y, color in series:
            ax.plot(x, y, label=name, zorder=3, **line_marker_kwargs(color))
            xs_all.extend(x)
        if ideal_linear and xs_all:
            lo, hi = min(xs_all), max(xs_all)
            ax.plot([lo, hi], [lo, hi], linestyle=(0, (5, 3)), color=INK_MUTED,
                    linewidth=1.3, label="Ideal", zorder=2)
        ax.set_yscale(yscale)
        _thread_axis(ax, sorted(set(xs_all)), xlabel=xlabel)
        ax.set_title(panel_title, fontsize=10)
        ax.margins(x=0.08)
    ax0 = axes[0]
    if percent:
        ax0.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
    else:
        ax0.yaxis.set_major_formatter(FuncFormatter(_plain_number))
        if yscale == "log":
            ax0.yaxis.set_minor_formatter(LogFormatter(base=10))
    if ylim is not None:
        ax0.set_ylim(*ylim)
    elif yscale != "log":
        ax0.set_ylim(bottom=0)
    ax0.set_ylabel(wrap(label, 40))
    for ax, (_, series) in zip(axes, panels):
        _end_labels(ax, [(x[-1], y[-1], c) for _, x, y, c in series], percent)
    handles, names = ax0.get_legend_handles_labels()
    fig.legend(handles, names, loc="lower center", ncol=len(names), frameon=False)
    if title:
        fig.suptitle(title, x=0.01, ha="left", fontweight="semibold", fontsize=11)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save_figure(fig, res_dir / filename)
