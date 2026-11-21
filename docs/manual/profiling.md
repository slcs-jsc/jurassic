# Profiling and benchmarking

This page describes a practical approach to investigating JURASSIC
performance. For user-facing configuration choices and scaling
considerations, see [Performance tuning](performance.md). For build
options, see [Building from source](building.md).

## Choose a representative workload

Use a reproducible input configuration and record the code revision,
compiler, optimization flags, hardware, thread/rank counts, lookup
tables, and input size. Compare runs only when these conditions are
understood. For retrieval workloads, record iteration counts and
convergence behavior as well as elapsed time.

## Identify the expensive work

The relevant cost centers depend on the workflow:

- Forward calculations: radiative-transfer evaluations and lookup-table
  access.
- Kernel calculations: repeated evaluations used to construct Jacobians.
- Retrievals: forward-model and Jacobian calls over multiple iterations.
- Parallel runs: work distribution, I/O, and synchronization overhead.

Measure before optimizing. A smaller elapsed time is not sufficient
evidence of an equivalent scientific result; compare the corresponding
outputs and relevant regression tests.

## Compiler and profiling options

The `src/Makefile` documents the build options `OPT`, `INFO`, and
`PROF`. For example, from `src/`:

```bash
make clean
make COMPILER=gcc OPT="-O2 -g"
```

The `PROF=1` option adds `-pg` for profiling builds, and `INFO=1`
enables compiler optimization reports where supported. The available
profiling tools and their output depend on the selected compiler and
platform. Do not assume that one profiler works with every CPU, GPU,
or MPI configuration.

## CPU and parallel benchmarking

Start with a serial or single-process baseline, then vary OpenMP
threads or MPI ranks only for executables and modes that support them.
Keep input data and scientific settings fixed. Record wall-clock time,
memory usage where available, and numerical differences.

For HPC batch runs, see [HPC workflows](hpc_workflows.md) and
[Parallelization](parallelization.md). The current MPI retrieval
implementation and workflow-level parallelism have different scaling
characteristics; do not interpret them as interchangeable.

## Reporting results

A useful benchmark record includes:

- Git commit or release and whether the checkout was modified.
- Compiler, libraries, build flags, and execution environment.
- Workload and relevant dimensions (rays, channels, retrieval cases).
- Runtime, rank/thread configuration, and repeatability.
- Regression-test status and any numerical differences.

For optimization changes, report both performance and correctness
relative to a documented baseline.
