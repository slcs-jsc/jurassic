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

def plot_agent_progress(source_file: Path, res_dir: Path) -> None:

    scores, status = [], []
    with open(source_file) as f:
        for line in f:
            item = json.loads(line)
            if item["score"] is not None:  # skip entries where profiling was skipped
                scores.append(item["score"])
                status.append("accepted" if item["accepted"] else "rejected")

    scores = np.array(scores)
    status = np.array(status)

    res_dir.mkdir(parents=True, exist_ok=True)
    iterations = np.arange(1, len(scores)+1)
    accepted = status == "accepted"

    fig, ax = plt.subplots()
    ax.scatter(iterations[accepted], scores[accepted], c="#28af36", label="Accepted", zorder=3)
    ax.scatter(iterations[~accepted], scores[~accepted], c="#e62a2a", label="Rejected", zorder=3)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Score (s)")
    ax.set_title("Optimization loop progress")
    ax.legend()

    fig.tight_layout()
    fig.savefig(res_dir / "agent_progress.png", bbox_inches="tight")
    plt.close(fig)

def plot_metrics(configs: list, res_dir: Path):
    threads = sorted({c["omp_threads"] for c in configs if c["group"] == "MEM_DP"})

    bandwidth = [
        next(
            (c["metrics_stat"]["Memory bandwidth [MBytes/s]"]["sum"]
                for c in configs if c["omp_threads"] == t and c["group"] == "MEM_DP"),
            None,
        )
        for t in threads
    ]
    if any(v is None for v in bandwidth):
        raise ValueError(f"Missing MEM_DP bandwidth data for one or more of {threads} threads")
    
    flops = [
        next(
            (c["metrics_stat"]["DP [MFLOP/s]"]["sum"]
             for c in configs if c["omp_threads"] == t and c["group"] == "FLOPS_DP"),
            None,
        )
        for t in threads
    ]
    if any(v is None for v in flops):
        raise ValueError(f"Missing FLOPS_DP data for one or more of {threads} threads")

    runtime = [
        next(
            (c["timers"]["TIMER_GROUP_ANALYSIS"]
             for c in configs if c["omp_threads"] == t and c["group"] == "FLOPS_DP"),
            None,
        )
        for t in threads
    ]
    if any(v is None for v in runtime):
        raise ValueError(f"Missing TIMER_GROUP_ANALYSIS for one or more of {threads} threads")

    _plot_scaling(threads, bandwidth, "Memory bandwidth [MB/s]", "#3ab9dc", res_dir, "bandwidth_scaling.png")
    _plot_scaling(threads, flops, "DP [MFLOP/s]", "#c76ce0", res_dir, "flops_scaling.png")
    _plot_scaling(threads, runtime, "TIMER_GROUP_ANALYSIS [s]", "#efb239", res_dir, "runtime_scaling.png", higher_is_better=False)

def _plot_scaling(
    threads: np.ndarray, 
    values: np.ndarray,
    label: str,
    color: str,
    res_dir: Path,
    filename: str,
    higher_is_better: bool = True,
) -> None:
    
    res_dir.mkdir(parents=True, exist_ok=True)

    ideal = (
        values[0] * (threads / threads[0])
        if higher_is_better
        else values[0] * (threads[0] / threads)
    )

    fig, ax = plt.subplots()
    ax.scatter(threads, values, "o-", color=color, label="Measured", zorder=3)
    ax.scatter(threads, ideal, "--", color="#898781", label="Ideal linear scaling", zorder=2)
    ax.set_xlabel("OMP_NUM_THREADS")
    ax.set_ylabel(label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(threads)
    ax.xaxis.set_major_formatter(plt.ScalarFormatter())
    ax.set_title(f"{label} vs thread count")
    ax.legend()
    fig.tight_layout()
    fig.savefig(res_dir / filename, bbox_inches="tight")
    plt.close(fig)






