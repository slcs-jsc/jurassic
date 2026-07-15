#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import math
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args(argv: list[str]):
    parser = argparse.ArgumentParser(description='Plot JURASSIC benchmark summaries.')
    parser.add_argument('summary_tsv', help='TSV file created by summarize_time_logs.py')
    parser.add_argument('--output', help='Output PNG path; defaults to <summary_tsv>.png')
    parser.add_argument('--title', default='JURASSIC benchmark summary')
    return parser.parse_args(argv[1:])


def load_rows(path: Path):
    with path.open(newline='') as f:
        reader = csv.DictReader(f, delimiter='	')
        rows = list(reader)
    if not rows:
        raise ValueError(f'No rows in {path}')
    return rows


def as_float(value: str):
    if value in ('', 'None', None):
        return None
    return float(value)


def as_int(value: str):
    if value in ('', 'None', None):
        return None
    return int(float(value))


def plot_omp(ax_runtime, ax_speedup, rows):
    x = [as_int(row['label']) for row in rows]
    batch_size = as_int(rows[0]['batch_size']) if rows else None
    benchmark = [as_float(row['benchmark_s']) for row in rows]

    per_case = [as_float(row['time_per_case_s']) for row in rows]
    throughput = [as_float(row['throughput_cases_per_s']) for row in rows]
    baseline = per_case[0] if per_case and per_case[0] not in (None, 0.0) else None

    # Use the first valid run as the scaling reference to keep the summary TSV
    # self-contained and avoid hard-coding a thread count here.
    speedup = [baseline / y if baseline is not None and y not in (None, 0.0) else math.nan for y in per_case]
    efficiency = [s / xi if not math.isnan(s) and xi not in (None, 0) else math.nan for s, xi in zip(speedup, x)]

    ax_runtime.plot(x, benchmark, marker='o', linewidth=2, label='Batch runtime')
    ax_runtime.set_xlabel('OMP_NUM_THREADS')
    ax_runtime.set_ylabel('Time [s]')
    title = 'CPU scaling'
    if batch_size is not None:
        title += f' (batch={batch_size})'
    ax_runtime.set_title(title)
    ax_runtime.grid(True, alpha=0.3)
    ax_runtime.legend()

    ax_speedup.plot(x, per_case, marker='o', linewidth=2, label='Time / case')
    twin = None
    if any(v is not None for v in throughput):
        # Plot throughput on a second axis so the time-per-case trend remains
        # readable even when the absolute scales differ strongly.
        twin = ax_speedup.twinx()
        twin.plot(x, throughput, marker='s', linewidth=2, color='tab:red', label='Throughput')
        twin.set_ylabel('Cases / s', color='tab:red')
        twin.tick_params(axis='y', labelcolor='tab:red')
    ax_speedup.plot(x, speedup, marker='^', linewidth=2, label='Speedup', color='tab:green')
    ax_speedup.plot(x, efficiency, marker='d', linewidth=2, label='Efficiency', color='tab:purple')
    lines = list(ax_speedup.get_lines())
    if twin is not None:
        lines += list(twin.get_lines())
    labels = [line.get_label() for line in lines]
    ax_speedup.legend(lines, labels, loc='best')
    ax_speedup.set_xlabel('OMP_NUM_THREADS')
    ax_speedup.set_ylabel('Per-case / relative metric')
    ax_speedup.set_title('CPU throughput / scaling')
    ax_speedup.grid(True, alpha=0.3)


def plot_batch(ax_runtime, ax_throughput, rows):
    x = [as_int(row['label']) for row in rows]
    benchmark = [as_float(row['benchmark_s']) for row in rows]
    per_case = [as_float(row['time_per_case_s']) for row in rows]
    throughput = [as_float(row['throughput_cases_per_s']) for row in rows]

    ax_runtime.plot(x, benchmark, marker='o', linewidth=2, label='Benchmark')

    # Batch sizes usually follow powers of two, so a log2 axis makes the GPU
    # scaling trend easier to compare across small and large batches.
    ax_runtime.set_xscale('log', base=2)
    ax_runtime.set_xticks(x)
    ax_runtime.set_xticklabels([str(value) for value in x])
    ax_runtime.set_xlabel('BATCH_SIZE')
    ax_runtime.set_ylabel('Time [s]')
    ax_runtime.set_title('GPU batch runtime')
    ax_runtime.grid(True, alpha=0.3)
    ax_runtime.legend()

    ax_throughput.plot(x, per_case, marker='o', linewidth=2, label='Time / case')
    ax_throughput.set_xscale('log', base=2)
    ax_throughput.set_xticks(x)
    ax_throughput.set_xticklabels([str(value) for value in x])
    if any(v is not None for v in throughput):
        twin = ax_throughput.twinx()
        twin.plot(x, throughput, marker='^', linewidth=2, color='tab:red', label='Throughput')
        twin.set_ylabel('Cases / s', color='tab:red')
        twin.tick_params(axis='y', labelcolor='tab:red')
        lines = ax_throughput.get_lines() + twin.get_lines()
        labels = [line.get_label() for line in lines]
        ax_throughput.legend(lines, labels, loc='best')
    else:
        ax_throughput.legend()
    ax_throughput.set_xlabel('BATCH_SIZE')
    ax_throughput.set_ylabel('Time / case [s]')
    ax_throughput.set_title('GPU throughput')
    ax_throughput.grid(True, alpha=0.3)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    path = Path(args.summary_tsv)
    rows = load_rows(path)
    mode = rows[0]['mode']
    output = Path(args.output) if args.output else path.with_suffix('.png')

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    fig.suptitle(args.title)

    if mode == 'omp':
        plot_omp(axes[0], axes[1], rows)
    elif mode == 'batch':
        plot_batch(axes[0], axes[1], rows)
    else:
        raise ValueError(f'Unsupported mode in TSV: {mode}')

    fig.savefig(output, dpi=150)
    plt.close(fig)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
