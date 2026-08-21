#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import argparse
import csv
import re
import sys
import statistics

AGG_FIELDNAMES = [
    'mode',
    'label',
    'batch_size',
    'n_reps',
    'benchmark_s_min',
    'benchmark_s_mean',
    'benchmark_s_stddev',
    'benchmark_s_relstddev',
    'time_per_case_s',
    'time_per_case_s_mean',
    'throughput_cases_per_s',
    'formod_reference_s',
    'timer_formod_s',
    'timer_total_s',
    'kernel_s',
    'logs',
]
 
PER_REP_FIELDNAMES = [
    'mode',
    'label',
    'rep',
    'batch_size',
    'formod_reference_s',
    'timer_formod_s',
    'benchmark_s',
    'timer_total_s',
    'kernel_s',
    'log',
]

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

def rep_value(path:Path, text: str):
    match = re.search(r'REP=(\d+)', text)
    if match:
        return int(match.group(1))
    match = re.search(r'\.rep(\d+)', path.name)
    if match: 
        return int(match.group(1))
    return None

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
    parser.add_argument('--tsv-out', help='Optional TSV output path for aggregated result over all samples')
    parser.add_argument('--per-rep-tsv-out', help='Optional TSV output path for per-repetition raw samples')
    return parser.parse_args(argv[1:])


def samples_from_logs(mode: str, pattern: str):
    paths = sorted(Path('.').glob(pattern))
    if not paths:
        raise FileNotFoundError(f'No logs matched pattern: {pattern}')

    samples = []
    for path in paths:
        text = path.read_text(errors='ignore')
        label = label_value(path, text, mode)
        rep = rep_value(path, text)
        formod_ref = extract_float(r'TIMER_FORMOD_REFERENCE = ([0-9.]+) s', text)
        formod = extract_float(r'TIMER_FORMOD = ([0-9.]+) s', text)
        total = extract_float(r'TIMER_TOTAL = ([0-9.]+) s', text)
        sample_timer = extract_float(r'TIMER_BENCHMARK_SAMPLE = ([0-9.]+) s', text)
        runtime = extract_float(r'RUNTIME: .*?mean= ([0-9.]+) s', text)
        kernel = extract_int(r'formod_batch\s+NVIDIA.*?device time\(us\): total=([0-9,]+)', text)

        # Prefer the dedicated benchmark timer, but fall back to the reported
        # runtime summary for older logs.
        benchmark = runtime if runtime is not None else sample_timer
        kernel_s = (kernel / 1e6) if kernel is not None else None
        batch_size = parse_batch_size(path, text)

        samples.append({
            'mode': mode,
            'label': label,
            'rep': rep,
            'batch_size': batch_size,
            'formod_reference_s': formod_ref,
            'timer_formod_s': formod,
            'benchmark_s': benchmark,
            'timer_total_s': total,
            'kernel_s': kernel_s,
            'log': path.name,
        })

    samples.sort(key=lambda row: (row['label'] is None, row['label'], row['rep'] or 0))
    return samples

def _mean(values):
    values = [v for v in values if v is not None]
    return statistics.fmean(values) if values else None

def _min(values):
    values = [v for v in values if v is not None]
    return min(values) if values else None

def _stddev(values):
    values = [v for v in values if v is not None]
    if len(values) < 2:
        return 0.0 if values else None                
    return statistics.stdev(values)

def aggregate_rows(samples: list[dict[str, object]]):
    groups: dict[object, list[dict[str, object]]] = {}
    
    for sample in samples:
        groups.setdefault(sample['label'], []).append(sample)

    rows = []

    for label in sorted(groups, key=lambda x: (x is None, x)):
        group = groups[label]
        benchmark_values = [s['benchmark_s'] for s in group]
        batch_size = group[0]['batch_size']

        benchmark_min = _min(benchmark_values)
        benchmark_mean = _mean(benchmark_values)
        benchmark_stddev = _stddev(benchmark_values)
        n_reps = sum(1 for v in benchmark_values if v is not None)

        per_case_min = benchmark_min / batch_size if benchmark_min is not None and batch_size > 0 else None
        per_case_mean = benchmark_mean / batch_size if benchmark_mean is not None and batch_size > 0 else None
        throughput_min = batch_size / benchmark_min if benchmark_min not in (None, 0) else None

        rows.append({
            'mode': group[0]['mode'],
            'label': label,
            'batch_size': batch_size,
            'n_reps': n_reps,
            'benchmark_s_min': benchmark_min,
            'benchmark_s_mean': benchmark_mean,
            'benchmark_s_stddev': benchmark_stddev,
            'benchmark_s_relstddev': (
                100.0 * benchmark_stddev / benchmark_mean
                if benchmark_stddev is not None and benchmark_mean not in (None, 0)
                else None
            ),
            'time_per_case_s': per_case_min,
            'time_per_case_s_mean': per_case_mean,
            'throughput_cases_per_s': throughput_min,
            'formod_reference_s': _mean([s['formod_reference_s'] for s in group]),
            'timer_formod_s': _mean([s['timer_formod_s'] for s in group]),
            'timer_total_s': _mean([s['timer_total_s'] for s in group]),
            'kernel_s': _mean([s['kernel_s'] for s in group]),
            'logs': ','.join(s['log'] for s in group),
        })
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]):
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
    headers = [first_col, 'Batch', 'N', 'Min', 'Mean', 'StDev', 'RelStDev', 'Time/Case', 'Throughput']
    table_rows = []
    for row in rows:
        label = '' if row['label'] is None else str(row['label'])
        n_reps = str(row['n_reps'])
        batch = '' if row['batch_size'] is None else str(row['batch_size'])
        per_case = '' if row['time_per_case_s'] is None else f"{float(row['time_per_case_s']):.6f} s"
        throughput = '' if row['throughput_cases_per_s'] is None else f"{float(row['throughput_cases_per_s']):.3f} /s"
        relstdev = '' if row['benchmark_s_relstddev'] is None else f"{float(row['benchmark_s_relstddev']):.1f}"
        table_rows.append([
            label,
            batch,
            n_reps,
            fmt_seconds(row['benchmark_s_min']),
            fmt_seconds(row['benchmark_s_mean']),
            fmt_seconds(row['benchmark_s_stddev']),
            relstdev,
            per_case,
            throughput,
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

    noise = [r for r in rows if r['benchmark_s_relstddev'] is not None and r['benchmark_s_relstddev'] > 5.0]
    if noise: 
        print('\nWARNING: relative stddev > 5% for:', file=sys.stderr)
        for row in noise:
            print(
                f"  label={row['label']} n={row['n_reps']} relstddev={row['benchmark_s_relstddev']:.1f}%",
                file=sys.stderr,
            )

def main(argv: list[str]) -> int:
    args = parse_args(argv)
    try:
        samples = samples_from_logs(args.mode, args.pattern)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    rows = aggregate_rows(samples)

    if args.tsv_out:
        write_tsv(Path(args.tsv_out), rows, AGG_FIELDNAMES)
    if args.per_rep_tsv_out: 
        write_tsv(Path(args.args.per_rep_tsv_out), samples, PER_REP_FIELDNAMES)

    print_table(args.mode, rows)
    return 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv))
