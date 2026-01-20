# Parallelization and performance

JURASSIC is designed for efficient execution on modern multicore CPUs
and high-performance computing (HPC) systems. Parallel execution is
achieved through a combination of workflow-level parallelism,
optional MPI task distribution in retrievals, and OpenMP
threading.

This page describes the actual parallelization mechanisms implemented
in JURASSIC, typical execution modes, and best practices for performance
tuning.

---

## Parallelization model overview

JURASSIC employs three complementary levels of parallelism:

- **Workflow-level parallelism**  
  Independent jobs (e.g. different days, orbits, ensembles) submitted
  via batch systems or job arrays.

- **MPI (Message Passing Interface)**  
  MPI is used to distribute independent retrieval tasks across processes.
  Optional and limited to retrieval executables only.

- **OpenMP**  
  Shared-memory parallelism within a single process, used by several
  computational kernels.

There is no global hybrid MPI–OpenMP model across the entire code
base. MPI is currently implemented only in the retrieval code.

---

## MPI parallelization

### Scope of MPI support

MPI parallelization is implemented exclusively in the retrieval
executable. No forward models, tools, or utilities use MPI internally.

MPI support is optional and enabled at compile time.

---

### What is parallelized?

At the MPI level, retrieval runs are parallelized over:

- independent retrieval cases,
- observation directories or datasets.

Each MPI rank processes a disjoint subset of retrieval tasks.

---

### Characteristics

- MPI is used only for task distribution.
- There is no communication between MPI ranks during the retrieval.
- No domain decomposition or collective algorithms are employed.
- Each rank performs its own I/O.

This design is robust and scales well for large numbers of independent
retrievals.

---

### Typical usage (retrieval only)

```bash
mpirun -np 16 ./retrieval ctl/retrieval.ctl
```

Running non-retrieval executables under `mpirun` provides no performance
benefit.

---

## OpenMP parallelization

### What is parallelized?

OpenMP is used to accelerate computationally intensive loops within a
single process, such as:

- radiative transfer calculations,
- spectral loops,
- internal numerical kernels.

OpenMP is available in both retrieval and non-retrieval executables.

---

### Controlling OpenMP threads

The number of OpenMP threads is controlled via:

```bash
export OMP_NUM_THREADS=4
```

OpenMP threads share memory and therefore incur low overhead.

---

## Combining MPI and OpenMP (retrieval only)

For retrievals, MPI and OpenMP can be combined:

- MPI distributes independent retrieval tasks across processes.
- OpenMP accelerates computations within each retrieval.

Example:

```bash
export OMP_NUM_THREADS=4
mpirun -np 8 ./retrieval ctl/retrieval.ctl
```

This configuration uses up to 32 CPU cores in total.

---

## Load balancing considerations

Since MPI distributes tasks statically, load balance depends on the
uniformity of retrieval cases.

Potential imbalance sources include:

- varying convergence behavior,
- different observation geometries,
- heterogeneous input datasets.

**Best practices:**

- group similar retrieval cases together,
- avoid mixing very small and very large problems in one run,
- benchmark representative workloads.

---

## I/O considerations

Parallel performance can be limited by file-system access.

Recommendations:

- store lookup tables on fast parallel file systems,
- ensure unique output paths for each MPI rank,
- minimize diagnostic output in large runs,
- reuse lookup tables across many retrievals.

---

## Numerical reproducibility

Due to parallel execution:

- floating-point round-off differences may occur,
- bitwise-identical results across different OpenMP thread counts or MPI
  layouts are not guaranteed.

These effects are typically small and scientifically insignificant but
should be considered in regression testing.

---

## Performance tuning tips

- Use MPI only for retrieval workloads.
- Prefer OpenMP for accelerating single-case computations.
- Avoid oversubscription (MPI ranks × OpenMP threads > available cores).
- Validate scaling with small test cases before large campaigns.
- Disable expensive diagnostics unless needed.

---

## Scaling expectations

JURASSIC scales well for:

- large numbers of independent retrieval cases,
- ensemble and campaign-style workflows,
- OpenMP-accelerated single-node execution.

Scaling efficiency depends primarily on I/O performance and workload
uniformity.

---

## Summary

JURASSIC does not implement a general hybrid MPI–OpenMP model across
the entire code base. MPI parallelization is currently limited to the
retrieval executables and is used solely for distributing independent
retrieval tasks. OpenMP provides shared-memory acceleration within a
single process, and large-scale parallelism is typically achieved at
the workflow level.

---

## Related pages

- Running JURASSIC
- Building from source
- Configuration
