# HPC workflows

This page describes recommended workflows for running JURASSIC on
high-performance computing (HPC) systems. It focuses on practical
aspects such as job organization, batch scheduling, data management,
and performance-oriented execution strategies.

The goal is to help users run **large production simulations and
retrievals** efficiently, reproducibly, and robustly.

---

## Typical HPC use cases

JURASSIC is commonly used on HPC systems for:

- processing large satellite datasets,
- global or regional forward simulations,
- ensemble sensitivity studies,
- large-scale retrieval campaigns,
- benchmarking and intercomparison experiments.

These use cases typically involve many independent observations and
benefit strongly from parallel execution.

---

## Directory and workflow structure

A recommended directory layout for directory-based HPC workflows is:

```text
project/
├── bin/            # JURASSIC executables
├── ctl/            # control files
├── input/
│   ├── atm/        # atmospheric profiles
│   ├── obs/        # observation geometry files
│   └── luts/       # spectroscopic lookup tables
├── output/
│   ├── radiance/
│   ├── diagnostics/
│   └── logs/
└── scripts/        # batch and helper scripts
```

Separating input, output, and scripts simplifies automation and
post-processing.

!!! note
    Large retrieval campaigns may instead use shared netCDF atmospheric,
    observation, and matrix files together with `SHARED_IO_PROFLIST`.
    This reduces the number of directories and small ASCII files that
    need to be managed by the file system.

---

## Batch job submission

### Example: Slurm batch script

```bash
#!/bin/bash
#SBATCH --job-name=jurassic_formod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun ./bin/formod ctl/run.ctl input/obs/obs.tab input/atm/atm.tab output/radiance/rad.tab
```

Key points:

- For non-MPI executables such as `formod`, use a single task and scale with
  `cpus-per-task` / `OMP_NUM_THREADS`.
- Set `OMP_NUM_THREADS` explicitly.
- Redirect output to per-job log files.

---

## Choosing MPI and OpenMP settings

### General guidelines

- Use **workflow-level parallelism** or job arrays for independent forward-model runs.
- Use **MPI** only for retrieval workloads that process independent directory
  lists or shared-file profile lists.
- Use **OpenMP** to accelerate selected computational loops, especially
  source-function table initialization and kernel/Jacobian calculations.
- Avoid oversubscription (total threads > physical cores).

### Example configurations

| Nodes | MPI tasks | OpenMP threads | Total cores |
|------:|----------:|---------------:|------------:|
| 1     | 8         | 4              | 32          |
| 2     | 16        | 4              | 64          |
| 4     | 32        | 4              | 128         |

Always benchmark representative workloads to find the optimal balance.

### Example: MPI retrieval

For retrieval workloads, MPI distributes independent retrieval cases
across ranks:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./bin/retrieval ctl/retrieval.ctl dirlist.txt
```

For shared-file retrieval workflows, pass `-` instead of a directory
list and provide a profile-index list:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./bin/retrieval ctl/retrieval.ctl - \
  SHARED_IO_PROFLIST input/proflist.txt
```

---

## Splitting large workloads

Very large datasets are often split into **independent chunks** that
are processed in separate jobs.

Common splitting strategies:

- by time (e.g. day, orbit, month),
- by latitude or longitude bands,
- by instrument scan segments.

Each job runs JURASSIC on a subset of the data, and results are merged
in a post-processing step.

---

## Lookup table management

Lookup tables are read frequently and can dominate I/O costs if not
handled carefully.

Recommendations:

- store LUTs on a **fast parallel file system**,
- prefer netCDF lookup tables for new workflows, because they use fast
  packed I/O and group multiple channels into one file per gas,
- avoid duplicating large LUTs across directories,
- reuse LUTs across many jobs,
- verify LUT availability before job submission.

On some systems, staging LUTs to node-local storage can improve
performance.

---

## I/O and output handling

- Write output to dedicated directories to avoid file-name collisions.
- Minimize diagnostic output for large production runs.
- Use binary or netCDF matrix formats (`MATRIXFMT = 2` or `3`) for large
  matrix outputs.
- Prefer shared netCDF workflows for large retrieval campaigns when
  practical, to avoid thousands of small ASCII files and directories.
- Clean up temporary files in long-running workflows.

---

## Job arrays and ensembles

Many HPC systems support **job arrays**, which are ideal for ensemble
runs or parameter sweeps.

Example (Slurm):

```bash
#SBATCH --array=0-31

./bin/formod ctl/run_${SLURM_ARRAY_TASK_ID}.ctl \
  input/obs/obs_${SLURM_ARRAY_TASK_ID}.tab \
  input/atm/atm_${SLURM_ARRAY_TASK_ID}.tab \
  output/radiance/rad_${SLURM_ARRAY_TASK_ID}.tab
```

Each array element runs an independent JURASSIC configuration.

---

## Fault tolerance and restarts

For long runs:

- keep individual jobs reasonably short,
- split workloads so failed jobs can be rerun easily,
- archive control files and logs for each job.

JURASSIC itself does not implement checkpointing, so job-level
fault tolerance is handled at the workflow level.

---

## Reproducibility and provenance

For reproducible HPC workflows:

- record the JURASSIC version string embedded in the executable,
- archive control files, input profiles, and LUT metadata,
- document compiler, MPI, and OpenMP settings,
- keep batch scripts under version control.

This is particularly important for large retrieval campaigns and
published results.

---

## Summary

HPC workflows with JURASSIC typically consist of many independent,
parallel jobs orchestrated by batch scripts and workflow logic outside
the core model.

Careful organization of inputs, outputs, and job configuration enables
efficient, scalable, and reproducible use of JURASSIC on modern
supercomputing systems.

---

## Related pages

- [Running JURASSIC](running.md)
- [Parallelization](parallelization.md)
- [Building from source](building.md)
