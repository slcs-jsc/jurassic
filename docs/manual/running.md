# Running JURASSIC

This page describes how to run JURASSIC applications from the command
line, both for local executions and for parallel runs on HPC systems.
It complements the [Quickstart](quickstart.md) by providing a more
systematic overview of runtime usage.

---

## General execution model

JURASSIC is provided as a set of small command-line applications (e.g.
`formod`, `kernel`, `retrieval`) that all follow the same basic runtime
pattern:

```bash
./application control.ctl [KEY VALUE ...]
```

Where:

- `application` is the executable to run,
- `control.ctl` is the control file,
- optional `KEY VALUE` pairs override control-file settings.

All applications:

1. read the control file,
2. read input profiles and geometry,
3. execute the forward model (and optionally kernels/retrieval),
4. write output files.

---

## Running a forward simulation

A typical forward-model run uses the `formod` application.

### Example

```bash
./formod run.ctl
```

This will:

- read `run.ctl`,
- read the atmosphere and observation files referenced by the control,
- compute radiances and transmittances,
- write results to `obs.tab` and any configured output files.

### Command-line overrides

Control parameters can be overridden directly on the command line:

```bash
./formod run.ctl RAYDS 5 REFRAC 0
```

This is useful for quick sensitivity tests without modifying the
control file.

---

## Running kernel calculations

Kernel calculations compute Jacobians (sensitivities) of radiances with
respect to atmospheric state variables.

```bash
./kernel run.ctl
```

Kernel output is written to diagnostic files and/or matrix output files,
depending on the configuration (see [Output diagnostics](output_diagnostics.md)).

Kernel calculations are typically more expensive than pure forward
simulations.

---

## Running retrievals

Retrieval applications perform inverse modelling using optimal
estimation.

```bash
./retrieval run.ctl
```

In a retrieval run:

- `obs.tab` usually contains measured (or synthetic) radiances,
- the forward model is evaluated iteratively,
- retrieved state vectors and diagnostics are written to output files.

Successful retrievals depend critically on consistent configuration of
the forward model, state vector, and a priori information.

---

## Standard input and output

### Input files

Most applications expect the following files to be present:

- control file (`*.ctl`)
- atmospheric profile file (e.g. `atm.tab`)
- observation geometry file (e.g. `obs.tab`)
- spectroscopic lookup tables (via `TBLBASE`)

### Output files

Depending on the application, output may include:

- updated `obs.tab` with simulated or fitted radiances,
- kernel and diagnostic files,
- matrix output (Jacobians, averaging kernels, error covariances),
- log output to standard output.

---

## Parallel execution

JURASSIC supports **hybrid MPI–OpenMP parallelization**.

### MPI execution

MPI parallelization is typically used to distribute independent
observations (rays) across processes.

Example:

```bash
mpirun -np 8 ./formod run.ctl
```

This runs the forward model on 8 MPI processes.

### OpenMP threading

Within each MPI process, OpenMP is used to exploit shared-memory
parallelism.

The number of OpenMP threads is controlled via the environment variable
`OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./formod run.ctl
```

This configuration uses up to 32 CPU cores in total.

---

## HPC batch systems

On HPC systems, JURASSIC is typically run inside a batch job script.
A minimal example (Slurm-style) is shown below:

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./formod run.ctl
```

Adjust the resource settings according to problem size and system
architecture.

---

## Performance considerations

- Use MPI to parallelize over many observations.
- Use OpenMP to accelerate ray tracing and radiative transfer within
  each process.
- Ensure lookup tables are stored on fast file systems.
- Avoid enabling expensive diagnostics unless needed.

Performance tuning is discussed further in the
[Parallelization](parallelization.md) section.

---

## Error handling and logging

- Most runtime errors are reported to standard output.
- Fatal configuration errors typically cause immediate program exit.
- Numerical warnings may indicate configuration or table issues and
  should be investigated.

Always check log output, especially when developing new configurations.

---

## Reproducibility

To ensure reproducible runs:

- record the Git version string embedded in the executable,
- archive control files and input profiles,
- document compiler and MPI/OpenMP settings.

---

## Summary

Running JURASSIC consists of executing small, purpose-specific
applications with a shared configuration and input-file model.
The design supports flexible workflows ranging from small test cases
to large-scale parallel production runs on HPC systems.

---

## Related pages

- [Quickstart](quickstart.md)
- [Configuration](configuration.md)
- [Parallelization](parallelization.md)
