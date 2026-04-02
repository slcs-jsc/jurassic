# Performance considerations

This page discusses performance aspects of JURASSIC, including expected
scaling behavior, computational cost drivers, and practical strategies
for achieving efficient runtimes on both local systems and HPC
platforms.

The focus is on **understanding where time is spent** and **how users
can influence performance** through configuration and workflow design,
given the *actual parallelization mechanisms implemented in JURASSIC*.

---

## Main performance drivers

The computational cost of a JURASSIC run is primarily determined by:

- number of **observations / rays**,
- number of **detector channels** (`ND`),
- number of **spectral windows** (`NW`),
- complexity of the **ray geometry** (limb vs. nadir),
- use of **Jacobians or retrievals**,
- choice of **lookup tables** and interpolation density.

Forward simulations are generally much cheaper than kernel or retrieval
runs.

---

## Forward model performance

### Radiative transfer cost

For forward modelling, runtime scales approximately with:

```text
O(N_obs × N_ray_segments × ND × NG)
```

where:

- `N_obs` is the number of observations,
- `N_ray_segments` depends on ray length and discretization (`RAYDS`, `RAYDZ`),
- `ND` is the number of detector channels,
- `NG` is the number of emitting gases.

Limb geometries typically require more ray segments than nadir
geometries and are therefore more expensive.

Forward-model executables are built with OpenMP support, but the most
visible OpenMP-parallel regions in the current code base are in lookup
table initialization and kernel finite-difference evaluation rather than
in the main `formod` ray loop itself.

---

## Kernel and retrieval performance

### Jacobians

Kernel calculations require evaluating sensitivities with respect to
each retrieved parameter. In JURASSIC this is done by finite differences,
so runtime typically increases by a factor of:

- **O(N_state)** relative to a pure forward run.

Because the Jacobians are evaluated by repeated forward-model calls,
kernel runs remain substantially more expensive than pure forward
simulations.

---

### Retrieval iterations

Retrieval runs require multiple forward-model and Jacobian evaluations.
Total runtime scales with:

- number of iterations until convergence,
- size of the state vector,
- cost of matrix operations.

Poorly chosen a priori constraints or noisy input data can significantly
increase iteration counts.

---

## Lookup table access

Lookup table interpolation is a frequent operation in JURASSIC and can
become a bottleneck if not handled efficiently.

Recommendations:

- prefer **binary table formats** over ASCII,
- place LUTs on fast local or parallel file systems,
- avoid unnecessary table reloads between runs.

For large campaigns, LUT I/O costs are typically amortized over many
simulations.

---

## Parallel scaling behavior

### Workflow-level scaling

The most common and scalable form of parallelism in JURASSIC is
**workflow-level parallelization**, where independent simulations or
retrievals are executed as separate jobs (e.g. job arrays, campaign
splitting).

This approach applies to *all* executables and scales trivially as long
as sufficient resources are available.

---

### MPI scaling

MPI parallelization is implemented only in the retrieval executables
and is used exclusively to distribute independent retrieval tasks
across MPI ranks.

- Each MPI rank processes a subset of retrieval cases.
- There is no communication between ranks during execution.
- Scaling is close to linear as long as enough retrieval cases are
  available.

Scaling efficiency decreases if the number of retrieval cases per MPI
rank becomes too small or if I/O dominates runtime.

---

### OpenMP scaling

OpenMP is used within a single process to accelerate selected
computationally intensive loops, such as lookup-table initialization and
kernel finite-difference calculations.

OpenMP scaling is typically limited by:

- memory bandwidth,
- cache behavior,
- load imbalance in inner loops.

Best performance is usually achieved with a moderate number of threads
per process.

---

## MPI vs. OpenMP considerations

There is no general hybrid MPI–OpenMP execution model across the full
tool suite.

General guidance:

- Use MPI only for retrieval workloads with many independent cases.
- Use OpenMP to accelerate single-case computations.
- Avoid oversubscription (MPI ranks × OpenMP threads > physical cores).
- Running non-retrieval executables under `mpirun` provides no benefit.

Optimal configurations depend on hardware and problem size and should
be determined empirically.

---

## Configuration parameters affecting performance

Several control parameters have a direct impact on performance:

- `RAYDS`, `RAYDZ`  
  Smaller step sizes increase accuracy but also increase runtime.

- `ND`, `NW`, `NG`  
  Increasing spectral or chemical complexity increases cost.

- `WRITE_MATRIX`, `WRITE_BBT`  
  Diagnostic output can significantly increase runtime and I/O.

Users should balance accuracy requirements against computational cost.

---

## Memory usage

Memory consumption in JURASSIC is generally modest compared to many
large-scale models but increases with:

- number of detector channels,
- size of lookup tables,
- enabled diagnostic matrices.

Memory usage is usually dominated by lookup tables rather than by
per-observation data.

---

## Benchmarking and validation

Performance tuning should always be accompanied by validation:

- compare results before and after performance-related changes,
- ensure numerical accuracy remains acceptable,
- benchmark representative workloads rather than minimal test cases.

Small configuration changes can have a large impact on both performance
and accuracy.

---

## Practical performance tips

- Start with example configurations and modify incrementally.
- Disable diagnostics for production runs.
- Use OpenMP primarily to accelerate kernel-heavy workloads and other
  OpenMP-enabled code paths.
- Use MPI only for retrieval campaigns with many independent cases.
- Split very large workloads into multiple jobs when appropriate.
- Monitor runtime and scaling behavior during pilot runs.

---

## Summary

JURASSIC performance is driven primarily by problem size, numerical
configuration, and workflow design.

MPI parallelization is limited to retrieval executables and is used
solely for distributing independent retrieval tasks. OpenMP provides
shared-memory acceleration within a single process, while large-scale
throughput is typically achieved via workflow-level parallelization.

Understanding these distinctions allows users to choose efficient and
robust execution strategies for both small experiments and large HPC
campaigns.

---

## Related pages

- [Parallelization](parallelization.md)
- [HPC workflows](hpc_workflows.md)
- [Configuration](configuration.md)
