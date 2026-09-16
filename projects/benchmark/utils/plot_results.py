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






