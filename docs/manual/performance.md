# Performance considerations

This page discusses performance aspects of JURASSIC, including expected
scaling behavior, computational cost drivers, and practical strategies
for achieving efficient runtimes on both local systems and HPC
platforms.

The focus is on **understanding where time is spent** and **how users
can influence performance** through configuration and workflow design.

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

---

## Kernel and retrieval performance

### Jacobians

Kernel calculations require evaluating sensitivities with respect to
each retrieved parameter. This typically increases runtime by a factor
of:

- **O(N_state)** relative to a pure forward run.

Analytic Jacobians significantly reduce this overhead compared to
finite-difference approaches but kernel runs are still substantially
more expensive than forward simulations.

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

### Strong scaling

JURASSIC shows good strong scaling as long as:

- there are enough independent observations to distribute,
- I/O does not dominate runtime.

Scaling efficiency decreases when the number of observations per MPI
process becomes too small.

### Weak scaling

Weak scaling is generally favorable because each MPI process performs a
similar amount of work when observations are added proportionally to
resources.

---

## MPI vs. OpenMP balance

Choosing an appropriate balance between MPI tasks and OpenMP threads is
crucial for performance.

General guidance:

- use MPI to distribute observations across nodes,
- use OpenMP to exploit shared-memory parallelism within a node,
- avoid oversubscription of CPU cores.

Optimal settings depend on node architecture and problem size and
should be determined empirically.

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

Small configuration changes can have a large impact on performance and
accuracy.

---

## Practical performance tips

- Start with example configurations and modify incrementally.
- Disable diagnostics for production runs.
- Use hybrid MPI–OpenMP execution on multicore nodes.
- Split very large workloads into multiple jobs.
- Monitor runtime and scaling behavior during pilot runs.

---

## Summary

JURASSIC is designed to deliver high performance for infrared radiative
transfer and retrieval applications by combining efficient spectral
approximations with scalable parallelization.

By understanding the main performance drivers and tuning configuration
and workflow parameters accordingly, users can achieve efficient and
robust runtimes for both small experiments and large HPC campaigns.

---

## Related pages

- [Parallelization](parallelization.md)
- [HPC workflows](hpc_workflows.md)
- [Configuration](configuration.md)
