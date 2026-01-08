# Parallelization and performance

JURASSIC is designed for efficient execution on modern multicore CPUs
and high-performance computing (HPC) systems. To achieve this, the code
implements a **hybrid MPI–OpenMP parallelization model** that allows it
to scale from laptops to large clusters.

This page describes the parallelization strategy, typical execution
modes, and best practices for performance tuning.

---

## Parallelization model overview

JURASSIC uses two complementary levels of parallelism:

- **MPI (Message Passing Interface)** for distributed-memory
  parallelism across nodes and processes
- **OpenMP** for shared-memory parallelism within a single process

This hybrid approach is well suited for radiative transfer applications,
where many calculations (e.g. individual rays or observations) are
largely independent but each involves non-trivial computational work.

---

## MPI parallelization

### What is parallelized?

At the MPI level, JURASSIC typically parallelizes over:

- independent **observations / rays**,
- sets of profiles or time steps,
- ensembles of forward or retrieval runs.

Each MPI process works on a subset of the total workload, minimizing
communication overhead.

### Characteristics

- MPI communication is relatively light-weight.
- Most computation occurs locally within each process.
- Scaling is close to linear as long as enough independent work items
  are available.

### Typical usage

```bash
mpirun -np 16 ./formod run.ctl
```

This launches 16 MPI processes, each handling a portion of the
observation set.

---

## OpenMP parallelization

### What is parallelized?

Within each MPI process, OpenMP is used to parallelize:

- ray tracing and radiative transfer along different rays,
- loops over spectral channels or windows,
- internal numerical kernels.

### Controlling OpenMP threads

The number of OpenMP threads is controlled via the environment variable:

```bash
export OMP_NUM_THREADS=4
```

OpenMP threads share memory within a process and therefore have low
overhead compared to MPI communication.

---

## Hybrid MPI–OpenMP execution

In a hybrid setup, MPI distributes work across processes, while OpenMP
accelerates computations within each process.

Example:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./formod run.ctl
```

This configuration uses up to **32 CPU cores** in total.

Hybrid execution is usually more efficient than pure MPI or pure
OpenMP, especially on multi-socket nodes.

---

## Load balancing considerations

For good scalability, it is important that each MPI process receives a
similar amount of work.

Potential load-imbalance sources include:

- uneven distribution of observation geometries,
- variable ray-path lengths (e.g. limb vs. nadir),
- retrieval runs with different convergence behavior.

**Best practices:**

- group similar observation types together,
- avoid mixing very small and very large workloads in one run,
- test scaling behavior for representative cases.

---

## Parallel retrievals

Retrieval applications benefit strongly from parallelization because:

- each observation can often be retrieved independently,
- forward-model and Jacobian calculations are computationally intensive.

MPI parallelization is usually applied across retrieval cases, while
OpenMP accelerates the internal linear algebra and radiative transfer.

---

## I/O considerations

Parallel performance can be limited by file-system access if not
handled carefully.

Recommendations:

- store lookup tables on fast local or parallel file systems,
- avoid excessive per-process file I/O,
- minimize diagnostic output in large production runs,
- reuse lookup tables across many runs to amortize I/O cost.

---

## Numerical reproducibility

Due to parallel execution:

- floating-point round-off differences may occur,
- bitwise-identical results across different MPI/OpenMP layouts are not
  guaranteed.

These differences are typically negligible for scientific applications
but should be considered when comparing results across platforms.

---

## Performance tuning tips

- Prefer **hybrid MPI–OpenMP** over pure MPI for multicore nodes.
- Choose MPI process counts that align with node topology.
- Tune `OMP_NUM_THREADS` to avoid oversubscription.
- Validate performance with small test cases before large runs.
- Disable expensive diagnostics unless explicitly needed.

---

## Scaling expectations

JURASSIC has demonstrated good scalability for:

- large numbers of independent observations,
- global simulations,
- long time series,
- ensemble studies.

Actual scaling depends on the balance between computation and I/O and
on the characteristics of the HPC system.

---

## Summary

The hybrid MPI–OpenMP parallelization in JURASSIC enables efficient use
of modern HPC systems while remaining flexible for smaller-scale runs.
With appropriate configuration, JURASSIC can process large atmospheric
datasets and retrieval problems within practical time constraints.

---

## Related pages

- [Running JURASSIC](running.md)
- [Building from source](building.md)
- [Configuration](configuration.md)
