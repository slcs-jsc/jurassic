# Running JURASSIC

This page describes how to run JURASSIC applications from the command line, both
for local executions and for parallel runs on HPC systems. It complements the
[Quickstart](quickstart.md) by providing a systematic overview of runtime usage.

---

## General execution model

JURASSIC is provided as a set of small command-line applications (e.g. `formod`,
`kernel`, `retrieval`). **All applications take a control file and, depending on
the application, additional positional input/output files on the command line.**
Control parameters may optionally be overridden via `KEY VALUE` pairs.

### Common pattern

```bash
./application <positional-args...> [KEY VALUE ...]
```

- `positional-args` are **application-specific** (see below).
- `KEY VALUE` pairs override control-file settings (e.g. `RAYDS 5 REFRAC 0`).
  Place overrides **after** the positional arguments.

### Application-specific command lines

| Application | Command line | Meaning of positional arguments |
|---|---|---|
| Forward model | `./formod <ctl> <obs_in> <atm_in> <rad_out>` | control file, observation geometry input, atmosphere input, radiance/observation output |
| Kernels (Jacobians) | `./kernel <ctl> <obs_in> <atm_in> <kernel_out>` | control file, observation geometry input, atmosphere input, kernel matrix output |
| Retrieval (optimal estimation) | `./retrieval <ctl> <dirlist>` | control file, text file listing working directories (one per line) |

> Notes  
> - For `formod` and `kernel`, the filenames (`<obs_in>`, `<atm_in>`, `<rad_out>`,
>   `<kernel_out>`) are usually `*.tab` files (e.g. `obs.tab`, `atm.tab`,
>   `rad.tab`, `kernel.tab`), but the actual names are up to you.
> - For `formod` and `kernel`, the control parameter `DIRLIST` can switch between
>   a *single run* and a *multi-directory run* (details below).

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
- run the forward model, and
- write simulated radiances to `rad.tab`.

### Multi-directory run with `DIRLIST`

If the control file sets `DIRLIST` to a filename (instead of `-`), `formod` will
loop over the directories listed there (one directory per line). In each working
directory it will use the same positional filenames, interpreted relative to the
working directory.

Example:

```bash
# in run.ctl:  DIRLIST dirlist.txt
./formod run.ctl obs.tab atm.tab rad.tab
```

If `dirlist.txt` contains e.g.

```text
case_0001
case_0002
case_0003
```

then `formod` will work in `case_0001/`, `case_0002/`, … and read/write:

- `case_0001/obs.tab`, `case_0001/atm.tab`, `case_0001/rad.tab`
- `case_0002/obs.tab`, `case_0002/atm.tab`, `case_0002/rad.tab`
- …

### Optional control parameters used by `formod`

`formod` additionally checks (via the control file or command-line overrides):

- `TASK` (e.g. to choose which forward-model task to run),
- `OBSREF` (optional reference observation data, if used by the selected task),
- `DIRLIST` (see above).

### Command-line overrides

```bash
./formod run.ctl obs.tab atm.tab rad.tab RAYDS 5 REFRAC 0
```

---

## Running kernel calculations (`kernel`)

Kernel calculations compute Jacobians (sensitivities) of radiances with respect
to atmospheric/state variables.

### Single run

```bash
./kernel run.ctl obs.tab atm.tab kernel.tab
```

This will write the kernel (Jacobian) matrix to `kernel.tab` (matrix format as
implemented by JURASSIC).

### Multi-directory run with `DIRLIST`

As with `formod`, `kernel` honors `DIRLIST` in the control file. If `DIRLIST` is
set to a list file, `kernel` loops over the listed working directories and reads
`<obs_in>` and `<atm_in>` and writes `<kernel_out>` inside each directory.

---

## Running retrievals (`retrieval`)

The retrieval application performs inverse modelling using optimal estimation.

```bash
./retrieval run.ctl dirlist.txt
```

Here, `dirlist.txt` is a plain text file containing one working directory per
line. For each directory, `retrieval` reads fixed-named inputs and writes
fixed-named outputs (relative to that directory).

### Inputs expected in each working directory

- `atm_apr.tab` — a priori atmospheric state  
- `obs_meas.tab` — measured (or synthetic) observations/radiances

### Outputs written in each working directory

At the end of the retrieval, JURASSIC writes at least:

- `atm_final.tab` — retrieved atmospheric state
- `obs_final.tab` — modelled observations for the retrieved state
- `matrix_kernel.tab` — final Jacobian (kernel) matrix

If analytical error analysis is enabled, additional matrices are written:

- `matrix_cov_apr.tab` — a priori covariance
- `matrix_cov_ret.tab` — retrieval covariance
- `matrix_corr.tab` — correlation matrix
- `matrix_gain.tab` — gain matrix
- `matrix_avk.tab` — averaging kernel matrix

(Exact outputs depend on the control parameters that enable/disable retrieval
diagnostics.)

---

## Standard file types

Depending on the workflow, typical `*.tab` files include:

- atmosphere profiles (e.g. `atm.tab`, `atm_apr.tab`, `atm_final.tab`)
- observation geometry and radiances (e.g. `obs.tab`, `obs_meas.tab`, `rad.tab`,
  `obs_final.tab`)
- matrices (e.g. `kernel.tab`, `matrix_kernel.tab`, `matrix_avk.tab`, …)
- spectroscopic/emissivity lookup tables (configured via `TBLBASE`, `TBLFMT`, …)

---

## Parallel execution

JURASSIC supports **hybrid MPI–OpenMP parallelization**.

### MPI execution

MPI parallelization is typically used to distribute independent observations
(rays) across processes.

```bash
mpirun -np 8 ./formod run.ctl obs.tab atm.tab rad.tab
```

### OpenMP threading

Within each MPI process, OpenMP is used to exploit shared-memory parallelism.
The number of OpenMP threads is controlled via `OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./formod run.ctl obs.tab atm.tab rad.tab
```

---

## HPC batch systems

On HPC systems, JURASSIC is typically run inside a batch job script. A minimal
example (Slurm-style) is shown below:

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./formod run.ctl obs.tab atm.tab rad.tab
```

Adjust the resource settings according to problem size and system architecture.

---

## Error handling and logging

- Most runtime errors are reported to standard output.
- Fatal configuration errors typically cause immediate program exit.
- Numerical warnings may indicate configuration or lookup-table issues and
  should be investigated.

Always check log output, especially when developing new configurations.

---

## Reproducibility

To ensure reproducible runs:

- record the Git version string embedded in the executable output,
- archive control files and input tab files,
- document compiler and MPI/OpenMP settings.

---

## Related pages

- [Quickstart](quickstart.md)
- [Configuration](configuration.md)
- [Parallelization](parallelization.md)
