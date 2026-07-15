# JURASSIC Benchmark Setup

This directory contains the repository-side setup for systematic CPU and GPU
benchmarks of JURASSIC forward-model workloads.

The focus is on representative `formod` performance across:

- viewing geometry: `zenith`, `nadir`, `limb`
- number of lines of sight / rays
- number of spectral channels
- gas-set complexity
- CPU OpenMP scaling
- GPU batch-throughput scaling

This is intentionally separate from regression tests in `tests/`. Benchmarks are
hardware-dependent, often scheduler-dependent, and are not intended to be part of
`make check`.

## Design Principles

- Treat `zenith`, `nadir`, and `limb` as equally important benchmark classes.
- Use realistic lookup-table inventories instead of only the minimal example cases.
- Keep CPU and GPU benchmark definitions aligned wherever possible.
- Measure throughput with metrics that match the actual execution mode.
- Start from stable, conservative OpenACC parallelization and iterate only with data.

## LUT Strategy

Benchmark runs assume NetCDF `tria` lookup tables and do not use the shipped
small example LUT directories. The baseline CTLs are configured for `TBLFMT = 3`,
and the benchmark runners resolve `TBLBASE` from `BENCH_TBLBASE` and write a
temporary active CTL into the run directory.

Current local default:

```text
$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6
```

On HPC systems, set `BENCH_TBLBASE` explicitly to the mounted or local benchmark LUT
path for that machine.

## External LUT Inventory

The current benchmark planning assumes the accurate 1/cm `tria` table family,
with a local default such as:

```text
$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6
```

and a site-specific HPC path supplied via `BENCH_TBLBASE`.

Important subsets of the `tria_1cm` inventory include:

- `nc_1e-6/`
- `tria_500/` ... `tria_2900/`

The `nc_1e-6` set currently exposes 36 gas tables, for example `H2O`, `CO2`, `O3`,
`CH4`, `N2O`, `CO`, `HNO3`, `NO2`, `O2`, and `N2`.

Not every gas is relevant in every spectral band. Some tables may be missing, may
contain only trivial values, or may be intentionally omitted. Benchmark cases should
therefore resolve the final gas list per spectral band instead of assuming that every
configured gas is active everywhere.

## Benchmark Axes

### Geometries

- `zenith`
- `nadir`
- `limb`

### Standard Channel Counts

The standard channel-scaling axis for version 1 is:

```text
1, 2, 4, 8, 16, 32, 64, 128, 256
```

This explicitly includes the `1`-channel case, which is important for later nadir GPU
throughput studies with little or no spectral parallelism.

Optional stretch cases such as `512`, `1024`, or full 1/cm spectra are intentionally
kept outside the version-1 default matrix because they may exceed practical memory
limits.

### Ray Sets

The benchmark ray-count classes are geometry-specific:

| Geometry | `small` | `medium` | `large` |
|---|---:|---:|---:|
| `zenith` | 16 | 64 | 256 |
| `nadir`  | 1  | 8  | 64  |
| `limb`   | 16 | 64 | 128 |

### Gas Sets

Gas sets are defined in HITRAN-like priority order and are intended to be filtered by
actual LUT availability in the selected band.

- `core`
- `priority_mid`
- `priority_full`

The exact ordered candidate lists are stored in `configs/gas_sets/`.

### Spectral Bands

The standard version-1 spectral bands are:

- `0500_0900`
- `0900_1500`
- `1500_2200`
- `2200_3000`

Definitions are stored in `configs/spectral_sets/`.

### Baseline Cases

The benchmark defaults are tied to three explicit `tria`-based baseline cases
listed in `configs/baseline_cases.tsv`:

| Case | Geometry | Control file | Rays | Gases | Channels | Table base |
|---|---|---|---:|---:|---:|---|
| `zenith_baseline` | `zenith` | `projects/benchmark/cases/zenith_baseline.ctl` | 64 | 7 | 32 | `BENCH_TBLBASE` |
| `nadir_baseline` | `nadir` | `projects/benchmark/cases/nadir_baseline.ctl` | 8 | 7 | 32 | `BENCH_TBLBASE` |
| `limb_baseline` | `limb` | `projects/benchmark/cases/limb_baseline.ctl` | 64 | 7 | 32 | `BENCH_TBLBASE` |

The baseline CTLs are separate from the shipped example projects. The example
projects under `projects/zenith`, `projects/nadir`, and `projects/limb` stay as-is
with their small test LUTs; the benchmark runners only use the dedicated baseline
CTLs and inject `TBLBASE` from `BENCH_TBLBASE` at runtime. These baseline cases
are the concrete comparison anchor points for CPU/GPU, compiler, and system results.

### Targets

#### CPU

OpenMP scaling axis:

```text
OMP_NUM_THREADS = 1, 2, 4, 8, 12
```

#### GPU

Batch-throughput scaling axis:

```text
BATCH_SIZE = 1, 8, 64, 256
```

## Version-1 Benchmark Families

### `geometry_baseline`

Purpose: stable cross-geometry baseline case.

- geometries: all three
- ray set: `medium`
- band: `0900_1500`
- channels: `32`
- gas set: `core`
- CPU: `OMP=1,12`
- GPU: `BATCH=1,8,64`

### `channel_scaling`

Purpose: quantify the effect of channel count.

- geometries: all three
- ray set: `medium`
- band: `0900_1500`
- gas set: `core`
- channels: `1,2,4,8,16,32,64,128,256`
- CPU: `OMP=1,12`
- GPU: `BATCH=1,8,64,256`

### `gas_scaling`

Purpose: quantify the effect of gas-set size and composition.

- geometries: all three
- ray set: `medium`
- band: `0900_1500`
- channels: `64`, `256`
- gas sets: `core`, `priority_mid`, `priority_full`
- CPU: `OMP=1,12`
- GPU: `BATCH=1,8,64`

### `spectral_band_scaling`

Purpose: compare physically different spectral regions.

- geometries: all three
- ray set: `medium`
- bands: all four standard bands
- channels: `64`, `256`
- gas sets: `core`, `priority_full`
- CPU: `OMP=1`
- GPU: `BATCH=1,8,64`

### `cpu_strong_scaling`

Purpose: OpenMP scaling of representative single-case workloads.

- geometries: all three
- ray set: `medium`
- band: `0900_1500`
- channels: `64`
- gas set: `core`
- CPU only: `OMP=1,2,4,8,12`

### `gpu_throughput_scaling`

Purpose: throughput scaling of the current stable OpenACC path.

- geometries: all three
- ray sets:
  - `zenith=medium`
  - `nadir=small`
  - `limb=medium`
- band: `0900_1500`
- channels: `1,16,64,256`
- gas sets: `core`, `priority_mid`
- GPU only: `BATCH=1,8,64,256`

### `ray_scaling`

Purpose: sensitivity to line-of-sight count.

- geometries: all three
- ray sets: `small`, `medium`, `large`
- band: `0900_1500`
- channels: `64`
- gas set: `core`
- CPU: `OMP=1,12`
- GPU: `BATCH=1,8,64`

## Reported Metrics

Every benchmark run should, as far as available, record:

- geometry
- ray set
- number of rays (`NR`)
- spectral band
- number of channels (`ND`)
- gas set
- number of active gases (`NG`)
- target (`cpu` or `gpu`)
- OpenMP thread count
- batch size
- `RUNTIME: execution= ... mean`
- time per case
- throughput
- derived speedup / efficiency

For `TASK=time`, the primary benchmark metric is the batch runtime. Derived metrics
such as time per case and throughput are computed from that runtime. `TIMER_TOTAL`
and other auxiliary timers may still be logged, but they are not the main comparison
metric in the current benchmark workflow.

## Current Repository Contents

- `configs/gas_sets/`:
  ordered gas candidate lists
- `configs/spectral_sets/`:
  standard benchmark band definitions
- `configs/rays.tsv`:
  geometry-dependent ray-count classes
- `configs/channel_counts.txt`:
  standard version-1 channel counts
- `configs/baseline_cases.tsv`:
  canonical `tria`-based baseline cases for `zenith`, `nadir`, and `limb`
- `cases/*.ctl`:
  dedicated baseline CTLs with runtime `TBLBASE` injection
- `scripts/plan_matrix.py`:
  prints the planned version-1 benchmark matrix as TSV
- `scripts/summarize_time_logs.py`:
  summarizes `TASK=time` benchmark logs into markdown and TSV tables
- `scripts/plot_benchmark_results.py`:
  creates PNG plots from summary TSV files using `matplotlib`
- `scripts/run_local_cpu.sh`:
  local notebook/workstation CPU runner
- `scripts/run_juwels_booster.sh`:
  JUWELS Booster runner for `sbatch` execution

## Next Steps

This setup is intentionally the first iteration. Likely follow-up work includes:

- resolving `priority_full` dynamically from actual LUT availability
- generating control files and observation geometries automatically
- extending the runner set for JURECA and JUPITER
- refining which benchmark families are required or optional
- adding guarded stretch cases for very large channel counts or broad spectra


## Runner Scripts

The repository now provides two initial system-specific runners based on explicit `tria` baseline cases:

- `scripts/run_local_cpu.sh`:
  local CPU benchmark runner for notebook or workstation use
- `scripts/run_juwels_booster.sh`:
  JUWELS Booster runner for batch submission via `sbatch`

Both runners currently execute `formod ... TASK time` benchmarks for one selected
baseline case from `configs/baseline_cases.tsv`, generate an active CTL with the
resolved `BENCH_TBLBASE`, and store all outputs under `projects/benchmark/runs/`.
They are meant
for stable system bring-up and reproducible timing runs before the full benchmark
matrix is automated.

### Local Notebook CPU Runner

Default behavior:

- case: `zenith_baseline`
- geometry: `zenith`
- control file: `projects/benchmark/cases/zenith_baseline.ctl`
- `BENCH_TBLBASE`: `$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6`
- threads: `1 2 4 8 12`
- CPU batch size: `64`
- compiler: `gcc`
- GPU: disabled

Default invocation:

```bash
cd projects/benchmark/scripts
./run_local_cpu.sh
```

Explicit invocation:

```bash
CASE_NAME=limb_baseline THREADS="1 6 12" CPU_BATCH_SIZE=64 COMPILER=clang ./run_local_cpu.sh
```

`CASE_NAME` currently accepts `zenith_baseline`, `nadir_baseline`, or `limb_baseline`.
Each baseline uses `core` gases, 32 channels in the `0900_1500` band, and a medium
ray count for its geometry. The CPU runners benchmark batch throughput by default
with `CPU_BATCH_SIZE=64`, so that `OMP_NUM_THREADS` measures the currently relevant
formod-batch parallelism instead of single-case latency. In the current clean local
notebook runs, `CPU_BATCH_SIZE=64` gives a useful compromise across geometries:
`zenith_baseline` largely saturates by `OMP=8-12`, while `limb_baseline` and
`nadir_baseline` still improve up to `OMP=12`. Override `BENCH_TBLBASE` when the
`tria` tables live elsewhere.

### JUWELS Booster Runner

Default behavior:

- case: `zenith_baseline`
- geometry: `zenith`
- control file: `projects/benchmark/cases/zenith_baseline.ctl`
- `BENCH_TBLBASE`: `$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6`
- target: `gpu`
- batches: `1 8 64 256`
- GPU compiler: `nvc`
- CPU compiler: `gcc`
- MPI: enabled

The script loads the JUWELS module stack internally and writes results to a run
directory named from the Slurm job ID. `CASE_NAME` selects the baseline geometry
and CTL template; `BENCH_TBLBASE` selects the actual `tria` directory; `CTLFILE` can
still override the template explicitly.

GPU invocation:

```bash
cd projects/benchmark/scripts
sbatch --export=ALL,BENCH_TBLBASE=/path/to/tria_1cm/nc_1e-6 run_juwels_booster.sh
```

CPU invocation on Booster:

```bash
cd projects/benchmark/scripts
sbatch --export=ALL,CASE_NAME=nadir_baseline,BENCH_TBLBASE=/path/to/tria_1cm/nc_1e-6,TARGET=cpu,THREADS="1 2 4 8 12" run_juwels_booster.sh
```

CPU+GPU invocation:

```bash
cd projects/benchmark/scripts
sbatch --export=ALL,CASE_NAME=limb_baseline,BENCH_TBLBASE=/path/to/tria_1cm/nc_1e-6,TARGET=both run_juwels_booster.sh
```

### Current Local CPU Observations

With the current local notebook setup, `CPU_BATCH_SIZE=64`, and the clean baseline
runs stored under `projects/benchmark/runs/`, the three baseline geometries show the
following qualitative behavior:

- `zenith_baseline`: strong gain from `OMP=1` to `OMP=8`, then near-saturation at `OMP=12`
- `limb_baseline`: continues to improve through `OMP=12`
- `nadir_baseline`: also continues to improve through `OMP=12`, but is much lighter per case than `zenith` or `limb`

This means that the local CPU optimum is geometry-dependent. For general notebook
benchmarking, `CPU_BATCH_SIZE=64` remains a reasonable default because it avoids the
clear underutilization seen with smaller batches while keeping runtime moderate.

### Optional Instrumentation

A useful follow-up extension is optional low-level instrumentation for selected
benchmark runs. This is not part of the current default workflow, but it would be
valuable for understanding scaling limits and runtime bottlenecks.

Candidate mechanisms include:

- CPU profilers and counters, for example `perf` or `likwid`
- OpenMP runtime environment reporting where useful
- NVIDIA OpenACC diagnostics such as `NVCOMPILER_ACC_TIME` and related `NV_ACC_*`
  or compiler/runtime diagnostic flags
- system-specific profiler hooks on HPC platforms

Potential uses:

- separating compute time from launch and data-movement overhead
- checking whether poor scaling comes from memory bandwidth, synchronization, or
  under-filled parallel regions
- comparing CPU and GPU runs with richer evidence than wall-clock time alone
- collecting supplementary data for compiler and system comparisons

If added later, this should remain opt-in and should write separate raw diagnostic
artifacts into the run directory, so that the standard benchmark workflow stays
lightweight and reproducible.

### Summary Files and Plots

The runners now generate both human-readable and machine-readable summaries:

- markdown: `summary*.md`
- TSV: `summary*.tsv`
- plots: `plot_*.png`

Current standard plots are:

- CPU runs: batch runtime, time per case, throughput, speedup, and efficiency versus `OMP_NUM_THREADS`; interpretation depends on `CPU_BATCH_SIZE`
- GPU runs: batch runtime, time per case, and throughput versus `BATCH_SIZE`

### Run Directory Layout

Each benchmark run writes into:

```text
projects/benchmark/runs/<run_id>/
```

This contains at least:

- `config.txt`
- `summary.md` or `summary.cpu.md` / `summary.gpu.md`
- raw `log.*` files
- generated `data/` snapshots

These run directories are intentionally ignored by git.
