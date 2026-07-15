#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import re
import sys


def extract_float(pattern: str, text: str):
    match = re.search(pattern, text, re.S)
    return float(match.group(1)) if match else None


def extract_int(pattern: str, text: str):
    match = re.search(pattern, text, re.S)
    return int(match.group(1).replace(',', '')) if match else None


def label_value(path: Path, text: str, mode: str):
    if mode == 'omp':
        match = re.search(r'OMP_NUM_THREADS=(\d+)', text)
        if match:
            return int(match.group(1))
        match = re.search(r'log\.omp(\d+)', path.name)
        return int(match.group(1)) if match else None
    if mode == 'batch':
        match = re.search(r'BATCH_SIZE = (\d+)', text)
        if match:
            return int(match.group(1))
        match = re.search(r'log\.batch(\d+)', path.name)
        return int(match.group(1)) if match else None
    raise ValueError(mode)


def parse_batch_size(path: Path, text: str):
    match = re.search(r'BATCH_SIZE = (\d+)', text)
    if match:
        return int(match.group(1))
    match = re.search(r'CPU_BATCH_SIZE=(\d+)', text)
    if match:
        return int(match.group(1))
    match = re.search(r'log\.batch(\d+)', path.name)
    if match:
        return int(match.group(1))
    return 1


def parse_args(argv: list[str]):
    parser = argparse.ArgumentParser(
        description='Summarize JURASSIC TASK=time benchmark logs.'
    )
    parser.add_argument('mode', choices=['omp', 'batch'])
    parser.add_argument('pattern', help='Glob pattern for log files, e.g. log.omp*')
    parser.add_argument('--tsv-out', help='Optional TSV output path')
    return parser.parse_args(argv[1:])


def rows_from_logs(mode: str, pattern: str):
    paths = sorted(Path('.').glob(pattern))
    if not paths:
        raise FileNotFoundError(f'No logs matched pattern: {pattern}')

    rows = []
    for path in paths:
        text = path.read_text(errors='ignore')
        label = label_value(path, text, mode)
        formod_ref = extract_float(r'TIMER_FORMOD_REFERENCE = ([0-9.]+) s', text)
        formod = extract_float(r'TIMER_FORMOD = ([0-9.]+) s', text)
        total = extract_float(r'TIMER_TOTAL = ([0-9.]+) s', text)
        sample = extract_float(r'TIMER_BENCHMARK_SAMPLE = ([0-9.]+) s', text)
        runtime = extract_float(r'RUNTIME: .*?mean= ([0-9.]+) s', text)
        kernel = extract_int(r'formod_batch\s+NVIDIA.*?device time\(us\): total=([0-9,]+)', text)

        # Prefer the dedicated benchmark timer, but fall back to the reported
        # runtime summary for older logs.
        benchmark = runtime if runtime is not None else sample
        kernel_s = (kernel / 1e6) if kernel is not None else None
        batch_size = parse_batch_size(path, text)
        per_case = None
        throughput = None
        if benchmark is not None and batch_size > 0:
            per_case = benchmark / batch_size
            throughput = batch_size / benchmark if benchmark > 0 else None
        rows.append({
            'mode': mode,
            'label': label,
            'batch_size': batch_size,
            'formod_reference_s': formod_ref,
            'timer_formod_s': formod,
            'benchmark_s': benchmark,
            'timer_total_s': total,
            'kernel_s': kernel_s,
            'time_per_case_s': per_case,
            'throughput_cases_per_s': throughput,
            'log': path.name,
        })

    rows.sort(key=lambda row: (row['label'] is None, row['label']))
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]):
    fieldnames = [
        'mode',
        'label',
        'batch_size',
        'formod_reference_s',
        'timer_formod_s',
        'benchmark_s',
        'timer_total_s',
        'kernel_s',
        'time_per_case_s',
        'throughput_cases_per_s',
        'log',
    ]
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt_seconds(value: object, digits: int = 3) -> str:
    if value is None:
        return ''
    return f'{float(value):.{digits}f} s'


def print_table(mode: str, rows: list[dict[str, object]]):
    first_col = 'OMP' if mode == 'omp' else 'Batch'
    headers = [first_col, 'Batch', 'Benchmark', 'Time/Case', 'Throughput', 'Log']
    table_rows = []
    for row in rows:
        label = '' if row['label'] is None else str(row['label'])
        per_case = '' if row['time_per_case_s'] is None else f"{float(row['time_per_case_s']):.6f} s"
        batch = '' if row['batch_size'] is None else str(row['batch_size'])
        throughput = '' if row['throughput_cases_per_s'] is None else f"{float(row['throughput_cases_per_s']):.3f} /s"
        table_rows.append([
            label,
            batch,
            fmt_seconds(row['benchmark_s']),
            per_case,
            throughput,
            str(row['log']),
        ])

    widths = [len(header) for header in headers]
    for row in table_rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    sep = '+' + '+'.join('-' * (width + 2) for width in widths) + '+'
    print(sep)
    print('| ' + ' | '.join(header.ljust(widths[idx]) for idx, header in enumerate(headers)) + ' |')
    print(sep)
    for row in table_rows:
        print('| ' + ' | '.join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + ' |')
    print(sep)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    try:
        rows = rows_from_logs(args.mode, args.pattern)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    if args.tsv_out:
        write_tsv(Path(args.tsv_out), rows)

    print_table(args.mode, rows)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
