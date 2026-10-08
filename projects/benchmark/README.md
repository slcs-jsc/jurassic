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
temporary active CTL into the run directory. The H2O table must have been
generated with RFM `H2O(sub)` for consistency with JURASSIC's MT_CKD 4.1
continuum.

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

Not every gas is relevant in every channel. Some tables may be missing, may contain
only trivial values, or may be intentionally omitted. `configs/channels_alt3.tsv`
therefore lists for each channel the gases that have a table there.

## Benchmark Axes

### Geometries

- `zenith`
- `nadir`
- `limb`

### Channels and gases

JURASSIC only computes a gas in a channel if a lookup table exists for that gas at the
channel's wavenumber. The runtime therefore depends on the number of active
(channel, gas) pairs, i.e. the pairs for which a lookup table exists, and not simply
on ND × NG.

The file `configs/channels_alt3.tsv` contains a list of 128 channels between 587 and
739 cm⁻¹, each listed with the gases for which a lookup table exists. For ND channels,
the ND channels are picked evenly spaced in the list, including the first and the last
one (`experiments/generate_ctl.py`), so they always cover the whole range.

The gas sets are defined in `configs/gas_sets/`. The name gives the number of gases,
and each set contains the previous one plus some more gases:

| set | gases |
|---|---|
| `ng04` | CO2, H2O, O3, HNO3 |
| `ng08` | `ng04` + CH4, N2O, NH3, SO2 |
| `ng13` | `ng08` + C2H2, H2O2, HCN, HF, NO2 |
| `ng18` | `ng13` + C2H6, COF2, N2O5, HCl, ClO |

`ng13` contains exactly the 13 gases that are active in all 128 channels. `ng18` adds
the 5 gases that are only active in parts of the range and is the full gas list.

### Reference Cases

The benchmark defaults are tied to three reference cases listed in
`configs/baseline_cases.tsv`, all with 32 channels from the channel list and all 18
gases (`ng18`, 458 active pairs):

| Case | Geometry | Control file | Rays | Gases | Channels |
|---|---|---|---:|---:|---:|
| `zenith_baseline` | `zenith` | `projects/benchmark/cases/zenith_baseline.ctl` | 64 | 18 | 32 |
| `nadir_baseline` | `nadir` | `projects/benchmark/cases/nadir_baseline.ctl` | 8 | 18 | 32 |
| `limb_baseline` | `limb` | `projects/benchmark/cases/limb_baseline.ctl` | 64 | 18 | 32 |

JURASSIC is compiled for at most 8 gases by default (`NG` in `src/jurassic.h`), so all
benchmark scripts build with `DEFINES="-DNG=18"`.

The reference CTLs are separate from the shipped example projects. The example
projects under `projects/examples/zenith`, `projects/examples/nadir`, and
`projects/examples/limb` stay as-is with their small test LUTs; the benchmark runners
only use the dedicated reference CTLs and inject `TBLBASE` from `BENCH_TBLBASE` at
runtime.

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

## Reported Metrics

Every benchmark run should, as far as available, record:

- geometry
- number of rays (`NR`)
- number of channels (`ND`)
- gas set
- number of gases (`NG`) and active (channel, gas) pairs
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

- `configs/channels_alt3.tsv`:
  channel list with the active gases per channel
- `configs/gas_sets/`:
  nested gas sets `ng04`, `ng08`, `ng13`, `ng18`
- `configs/baseline_cases.tsv`:
  reference cases for `zenith`, `nadir`, and `limb`
- `cases/*.ctl`:
  reference CTLs with runtime `TBLBASE` injection
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
Each case uses all 18 gases (`ng18`), 32 channels from `configs/channels_alt3.tsv`,
and the ray count listed in the reference case table above. The CPU runners benchmark batch throughput by default
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

Current JUWELS Booster reference results with the corrected 1/cm `tria`
benchmark channels on Friday, July 17, 2026 are:

| Case | Batch | Time/Case [s] | Throughput [cases/s] |
|---|---:|---:|---:|
| `zenith_baseline` | 256 | 1.0612 | 0.942 |
| `zenith_baseline` | 512 | 0.5370 | 1.862 |
| `zenith_baseline` | 1024 | 0.3084 | 3.243 |
| `zenith_baseline` | 2048 | 0.1723 | 5.805 |
| `nadir_baseline` | 256 | 0.1202 | 8.316 |
| `nadir_baseline` | 512 | 0.06535 | 15.302 |
| `nadir_baseline` | 1024 | 0.03869 | 25.846 |
| `nadir_baseline` | 2048 | 0.02675 | 37.382 |
| `limb_baseline` | 256 | 0.8211 | 1.218 |
| `limb_baseline` | 512 | 0.4234 | 2.362 |
| `limb_baseline` | 1024 | 0.2536 | 3.944 |
| `limb_baseline` | 2048 | 0.1352 | 7.399 |

Within this setup, `2048` is the largest successful JUWELS Booster batch size
observed so far for `zenith_baseline`; a targeted `4096` run failed with
`CUDA_ERROR_OUT_OF_MEMORY`, so `2048` is the current practical upper bound for
that case on a single A100.

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
- NVIDIA OpenACC diagnostics such as `NVCOMPILER_ACC_TIME`, `NVCOMPILER_ACC_NOTIFY`,
  and related `NV_ACC_*` or compiler/runtime diagnostic flags
- system-specific profiler hooks on HPC platforms

Potential uses:

- separating compute time from launch and data-movement overhead
- checking whether poor scaling comes from memory bandwidth, synchronization, or
  under-filled parallel regions
- comparing CPU and GPU runs with richer evidence than wall-clock time alone
- collecting supplementary data for compiler and system comparisons

The JUWELS Booster runner already exposes `ACC_TIME` and `ACC_NOTIFY` as opt-in
diagnostic controls. Both are disabled by default so standard benchmark runs stay
lightweight and reproducible; enable them explicitly when raw compiler/runtime
diagnostics are needed in the batch logs.

### Summary Files and Plots

The runners now generate both human-readable and machine-readable summaries:

- plain-text ASCII table: `summary*.txt`
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


# LIKWID profiling

https://github.com/rrze-hpc/likwid

### Measurement Regions
base: LIKWID_MARKER_START/STOP("formod") located in formod_batch() (jurassic.c:3619/3622)

### Runner Scripts 
| Script  | Config   | Metrics    | Purpose |
| :---:   | :---: | :---: | :---: |
| run_noise_floor.sh | single thread, fixed batch size -> N identical runs | Runtime, MEM_DP | Determine measurement noise |
| run_roofline.sh | varies problem size, single thread | FLOPS_DP + MEM_DP → operational intensity | Memory- or compute-bound? |
| run_scaling.sh | varies thread count (1,2,4,8,12,24) + 48 SMT separately | MEM_DP, Runtime | Analyse scaling and saturation |
| run_tma.sh | compare single thread vs max physical thread count (24) | TMA, Cache volume + miss ratio | Perform Top-down Microarchitecture Analysis |
| run_compare_ab.sh | varies code, single thread, same job + same node | Write-/Call-/Read- volume, Runtime | Compare efficiency of two code versions |

### Results - Forward Model (CPU version)

#### Noise 

   label  thr    group  batch | runtime/call [s] | Memory data volume | Memory bandwidth
--------------------------------------------------------------------------------------------------------
   noise    1   MEM_DP     48 | 4.696, 4.686, stdev=0.03716, cv=0.8% | 1.465, 1.467, stdev=0.01842, cv=1.3% | 6.368, 6.357, stdev=0.06141, cv=1.0%
Largest observed coefficient of variation: 1.3%  (config=('noise', 1, 'MEM_DP', 48), metric='Memory data volume [GBytes]')

Recommended significance threshold for further comparisons: differences below ~1.3% of the median are not distinguishable from measurement noise at this repetition count.

#### Parallel efficiency
How does performance scale with computational resources?
(Case: zenith, tria_1cm/nc_1e-6/tria, batch size: 48, strong scaling)
[Thread Count vs. Runtime](results/baseline/e2_wallclock_scaling.png)
[Thread Count vs. Speedup](results/baseline/e2_wallclock_speedup.png)
[Thread Count vs. Memory bandwith](results/baseline/e2_memory_bandwidth_scaling.png)

SMT gives no Speedup and produces 3x higher DRAM read traffic. 

##### Resource contention
Does one core's TMA profile change when 23 other cores are active and sharing the same memory, versus running alone?

