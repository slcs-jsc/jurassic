#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import os
import csv
import math
import sys
import json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def boxplot(
    threads: np.ndarray, 
    value_arr: np.ndarray,
    label: str, 
    res_dir: Path,
    filename: str,
    show_outliers: bool = False
) -> None:
    res_dir.mkdir(parents=True, exist_ok=True)
    threads = np.asarray(threads, dtype=float)

    fig, ax = plt.subplots()
    if show_outliers:
        ax.boxplot(value_arr)
    else:
        ax.boxplot(value_arr, sym='')

    ax.set_xlabel("threads (physical cores)")
    ax.set_ylabel(label)
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_title(f"{label} vs thread count")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(res_dir / filename, bbox_inches="tight")
    plt.close(fig)

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
    smt_points=None
) -> None:
    """
    ideal: "linear" | "constant" | None
    higher_is_better: True for throughput-like metrics (bandwidth, FLOP/s);
                      False for cost-like metrics (runtime).
    smt_points: optional list of plotted as separate markers (e.g. result for SMT-threads)
    """
    res_dir.mkdir(parents=True, exist_ok=True)
    threads = np.asarray(threads, dtype=float)
    values = np.asarray(values, dtype=float)

    fig, ax = plt.subplots()
    ax.plot(threads, values, "o-", color=color, label="Measured (phys cores)", zorder=3)

    if ideal == "linear":
        ideal_vals = (
            values[0] * (threads / threads[0]) if higher_is_better
            else values[0] * (threads[0] / threads)
        )
        ax.plot(threads, ideal_vals, "--", color="#898781", label="Ideal linear scaling", zorder=2)

    elif ideal == "constant":
        ax.axhline(values[0], linestyle="--", color="#898781",
                   label="Ideal (constant, same work/call)", zorder=2)

    if stream_ceiling is not None:
        ax.axhline(stream_ceiling, linestyle=":", color="#c0392b",
                   label=f"STREAM bandwidth ceiling ({stream_ceiling:.0f})", zorder=2)

    if smt_points:
        for t, v, l in smt_points:
            ax.scatter([t], [v], marker="s", s=60, color="#e67e22", zorder=4, label=l)

    ax.set_xlabel("threads (physical cores)")
    ax.set_ylabel(label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    xticks = list(threads) + ([smt_points[0][0]] if smt_points else [])
    ax.set_xticks(xticks)
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_title(f"{label} vs thread count")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(res_dir / filename, bbox_inches="tight")
    plt.close(fig)







def plot_compare_scaling(
    series: list,
    label: str,
    res_dir: Path,
    filename: str,
    title: str | None = None,
    ideal: bool = True,
    log_y: bool = True
) -> None:
    """
    Plot several variants of one metric against thread count in a single diagram.

    series: list of (name, threads, values, errors, color); errors may be None.
    The first entry is the reference. ideal adds a dashed ideal-strong-scaling
    line from its first point. log_y=False uses a linear y axis from zero,
    which shows the size of a reduction directly. Every point of the other
    entries is annotated with its value as a percentage of the reference at
    the same thread count. Meant for cost-like metrics (runtime, memory
    volume), where lower is better, so below 100 % is an improvement.
    """
    res_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots()
    ref_threads = np.asarray(series[0][1], dtype=float)
    ref_values = np.asarray(series[0][2], dtype=float)

    if ideal:
        ax.plot(ref_threads, ref_values[0] * ref_threads[0] / ref_threads, "--",
                color="#898781", label="Ideal linear scaling (baseline)", zorder=2)

    for k, (name, threads, values, errors, color) in enumerate(series):
        threads = np.asarray(threads, dtype=float)
        values = np.asarray(values, dtype=float)
        if errors is None:
            ax.plot(threads, values, "o-", color=color, label=name, zorder=3)
        else:
            ax.errorbar(threads, values, yerr=np.asarray(errors, dtype=float),
                        fmt="o-", color=color, capsize=3, label=name, zorder=3)

        # Percentage of the reference at the thread counts both measured.
        if k > 0:
            for t, v in zip(threads, values):
                hit = np.where(ref_threads == t)[0]
                if len(hit) and ref_values[hit[0]] > 0:
                    ax.annotate(f"{100 * v / ref_values[hit[0]]:.0f} %", (t, v),
                                textcoords="offset points", xytext=(0, -14),
                                ha="center", fontsize=8, color=color)

    ax.set_xlabel("threads (physical cores)")
    ax.set_ylabel(label)
    ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    else:
        ax.set_ylim(0, 1.15 * max(float(np.max(s[2])) for s in series))
    ax.set_xticks(sorted({float(t) for s in series for t in s[1]}))
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_title(title or f"{label} vs thread count")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(res_dir / filename, bbox_inches="tight")
    plt.close(fig)
