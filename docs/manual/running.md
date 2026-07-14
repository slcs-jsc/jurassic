# Running JURASSIC

This page describes how to run JURASSIC applications from the command
line, both for local executions and for parallel runs on HPC systems. It
complements the [Quickstart](quickstart.md) by providing a systematic
overview of runtime usage.

---

## General execution model

JURASSIC is provided as a set of small command-line applications (e.g.
`formod`, `kernel`, `retrieval`). Many of the scientific applications take a
control file plus additional positional input/output files on the command
line. Smaller utilities such as `day2doy`, `time2jsec`, `planck`, and
`brightness` accept only direct numeric arguments. For control-file-driven
applications, parameters may optionally be overridden via `KEY VALUE` pairs.

### Common pattern

```bash
./application <positional-args...> [KEY VALUE ...]
```

- `positional-args` are **application-specific** (see below).
- `KEY VALUE` pairs apply only to applications that read a control file.
- Place overrides **after** the positional arguments.

---

## Application-specific command lines

| Application | Command line | Meaning of positional arguments |
|------------|--------------|---------------------------------|
| Forward model | `./formod <ctl> <obs_in> <atm_in> <rad_out>` | control file, observation geometry input, atmosphere input, radiance/observation output |
| Kernels (Jacobians) | `./kernel <ctl> <obs_in> <atm_in> <kernel_out>` | control file, observation geometry input, atmosphere input, kernel matrix output |
| Retrieval (optimal estimation) | `./retrieval <ctl> <dirlist-or->` | control file, either a text file listing working directories or `-` for shared-file/profile-list mode |

!!! note
    For `formod` and `kernel`, the filenames are usually `*.tab` files,
    but the actual names are user-defined.

    For `formod` and `kernel`, the control parameter `DIRLIST` can
    switch between a single-directory and a multi-directory workflow.

---

## Running a forward simulation (`formod`)

### Single run

```bash
./formod run.ctl obs.tab atm.tab rad.tab
```

This will:

- read control parameters from `run.ctl`,
- read observation geometry from `obs.tab`,
- read the atmospheric state from `atm.tab`,
- run the forward model,
- write simulated radiances or brightness temperatures to `rad.tab`.

### Benchmark mode with `TASK=time`

Use `TASK=time` to run the built-in runtime benchmark:

```bash
./formod run.ctl obs.tab atm.tab rad.tab TASK time
```

Optional benchmark-specific parameters are:

- `EXECUTION batch|scalar` to select the batch path or the scalar debug path.
- `BATCH_SIZE <n>` to set the batch size for `EXECUTION=batch`.

Notes for interpreting the benchmark output:

- `TASK=time` first performs one regular forward calculation and writes its result to
  `<rad>`.
- This initial reference run is timed as `FORMOD_REFERENCE`.
- The actual benchmark loop is timed separately as `BENCHMARK_SAMPLE` and printed as
  `RUNTIME: execution= ...`.
- For throughput comparisons, use the reported benchmark runtime, not
  `TIMER_FORMOD` or `FORMOD_REFERENCE`.
- When `BATCH_SIZE > 1`, the benchmark internally evaluates perturbed atmospheric
  cases and returns the first batch element as the written reference output.

---

### Multi-directory run with `DIRLIST`

If the control file sets `DIRLIST` to a filename (instead of `-`),
`formod` loops over the directories listed there (one directory per
line). In each working directory it uses the same positional filenames,
interpreted relative to that directory.

Example:

```bash
# in run.ctl:  DIRLIST dirlist.txt
./formod run.ctl obs.tab atm.tab rad.tab
```

A minimal directory layout for this workflow is:

```text
run.ctl
dirlist.txt
case_001/
  obs.tab
  atm.tab
case_002/
  obs.tab
  atm.tab
```

with `dirlist.txt` containing:

```text
case_001
case_002
```

For this example, JURASSIC reads `case_001/obs.tab` and
`case_001/atm.tab`, writes `case_001/rad.tab`, and then repeats the same
file-name pattern for `case_002`.

!!! note
    Lookup tables are loaded once before the directory loop starts and
    then reused for all cases in the `DIRLIST` run. This avoids the
    expensive table-loading step for every individual case and is one of
    the main performance advantages of this workflow.

    `TBLBASE` is therefore resolved relative to the directory where
    `formod` is launched, or as an absolute path, not relative to each
    case directory.

---

## Running kernel calculations (`kernel`)

Kernel calculations compute Jacobians (sensitivities) of the simulated
output with respect to atmospheric/state variables.

```bash
./kernel run.ctl obs.tab atm.tab kernel.tab
```

As with `formod`, `kernel` honors `DIRLIST` for multi-directory
workflows.

---

## Running retrievals (`retrieval`)

The retrieval application performs inverse modelling using optimal
estimation.

```bash
./retrieval run.ctl dirlist.txt
```

Here, `dirlist.txt` is a plain text file containing one working directory
per line. Each directory is processed independently.

For shared netCDF workflows, retrieval can also read and write selected
profile records from common files. In that mode, pass `-` instead of a
directory list and provide:

- `SHARED_IO_PROFLIST` — a text file with one profile index per line
- shared input files such as `SHARED_IO_ATM_APR_FILE` and
  `SHARED_IO_OBS_MEAS_FILE`
- optional shared output files such as `SHARED_IO_ATM_FINAL_FILE`,
  `SHARED_IO_OBS_FINAL_FILE`, and `SHARED_IO_MATRIX_KERNEL_FILE`

Example:

```bash
./retrieval run.ctl - \
  SHARED_IO_PROFLIST proflist.txt \
  ATMFMT 3 OBSFMT 3 MATRIXFMT 3 \
  SHARED_IO_ATM_APR_FILE atm_apr.nc \
  SHARED_IO_OBS_MEAS_FILE obs_meas.nc \
  SHARED_IO_ATM_FINAL_FILE atm_final.nc \
  SHARED_IO_OBS_FINAL_FILE obs_final.nc \
  SHARED_IO_MATRIX_KERNEL_FILE matrix_kernel.nc
```

In shared-file mode, the profile indices listed in `SHARED_IO_PROFLIST`
select the records read from shared atmospheric and observation files and
the records written to shared output files.

At least one of the following must be provided:

- a directory list as the second positional argument, or
- `SHARED_IO_PROFLIST` for shared-file profile processing.

### Inputs expected in each working directory

- `atm_apr.tab` — a priori atmospheric state
- `obs_meas.tab` — measured (or synthetic) observations

### Outputs written in each working directory

- `atm_final.tab` — retrieved atmospheric state
- `obs_final.tab` — modelled observations
- `matrix_kernel.tab` — final Jacobian matrix

Additional matrices (covariance, averaging kernels, etc.) are written if
enabled via control parameters.

---

## Parallel execution

Parallel execution in JURASSIC is achieved through workflow-level
parallelism, optional MPI in retrievals, and OpenMP threading.
There is no general hybrid MPI–OpenMP execution model across the full
tool suite. MPI is limited to retrieval workflows.

---

### MPI execution

MPI parallelization is implemented only in the retrieval code and is
used to distribute independent retrieval cases across MPI ranks.

```bash
mpirun -np 8 ./retrieval run.ctl dirlist.txt
```

Each MPI rank processes a subset of directories from `dirlist.txt`. MPI
ranks do not communicate during execution.

!!! note
    Running non-retrieval executables (`formod`, `kernel`, tools) under
    `mpirun` provides no performance benefit.

---

### OpenMP threading

OpenMP is used within a single process to accelerate computationally
intensive loops. It is available for both retrieval and non-retrieval
executables.

The number of OpenMP threads is controlled via:

```bash
export OMP_NUM_THREADS=4
```

---

### Combining MPI and OpenMP

For retrievals, MPI and OpenMP can be combined:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./retrieval run.ctl dirlist.txt
```

This uses up to 32 CPU cores in total.

---

## HPC batch systems

On HPC systems, JURASSIC is typically run inside a batch job script.
Example (Slurm):

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./retrieval run.ctl dirlist.txt
```

Adjust resource requests according to problem size and architecture.

---

## Error handling and logging

- Runtime errors are reported to standard output.
- Fatal configuration errors typically cause immediate program exit.
- Numerical warnings may indicate configuration or lookup-table issues.

Always inspect log output, especially when developing new workflows.

---

## Reproducibility

To ensure reproducible runs:

- record the Git version string printed by the executables,
- archive control files and input data,
- document compiler and MPI/OpenMP settings.

---

## Related pages

- [Quickstart](quickstart.md)
- [Configuration](configuration.md)
- [Parallelization](parallelization.md)
