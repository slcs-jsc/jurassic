# CPU/OpenMP Benchmark

This directory contains scripts to profile and evaluate the performance of JURASSIC's forward model ('formod').

The framework is divided into two distinct analysis runs:
1. **`run_roofline`**: Determines if the algorithm is compute-bound or memory-bound.
2. **`run_scaling`**: Measures multi-core scalability, SMT overhead, and cross-socket NUMA impact.

## Prerequisites

The framework relies on **LIKWID** for performance counter collection and thread/memory pinning.

* **Affinity Control**: Standard OpenMP variables (`OMP_PLACES`, `OMP_PROC_BIND`) are explicitly unset in `base.sh`. LIKWID handles explicit core placement via hardware expressions (`E:S<socket>:<count>`).
* **Memory Placement**: Control memory allocations via LIKWID-native flags (`-m` for local NUMA binding, `-i` for striated interleave).
* **Compiler Flags**: Ensure JURASSIC is compiled with vectorization enabled (`-O3 -march=native` or `-march=znver2`) to allow the application to utilize the AVX2 and FMA3 hardware pipelines monitored by the `FLOPS_DP` counters.

## Scaling ('run_scaling_*')

### Evaluation Axes
1. **Intra-Socket Scaling**: Scale across the physical cores of a single socket
    * **Strong Scaling (`strong`)**: `BATCH_SIZE` is fixed. Runtime plots directly scale as "Time-to-Solution".
    * **Weak Scaling (`weak`)**: `BATCH_SIZE` scales linearly with the thread count (\(BASE\_BATCH \times t\))
2. **SMT (Hyperthreading)**: Runs using logical hardware threads on Socket 0 to measure pipeline occupancy stalls and cache behavior
3. **Inter-Socket / Numa Layout**: 
   * **`inter_compact`**: All 64 threads loaded onto Socket 0 (`E:S0:64`, memory allocation restricted locally via `-m`)
   * **`inter_spread`**: Threads split evenly across both physical processors (`E:S0:32@E:S1:32`, memory allocation striped via `-i`

### Execution
Submit the job to Slurm, specifying the scaling mode:
```bash
sbatch --export=SCALING_MODE=strong run_roofline_jureca.sh
sbatch --export=SCALING_MODE=weak run_roofline_jureca.sh
```

### Parsing & Visualization
The evaluation script parses the structured output labels (`intra_socket_strong`, `inter_spread_weak`, etc.) 
and plots normalized speedup, efficiency (\(E = \frac{T_1}{T_n}\) for weak scaling), and memory bandwidth behavior:
```python
python evaluate_scaling.py /path/to/run_dir --stream-bw <measured_stream_ceiling>
```

This workflow times the current `master` forward model for the checked-in limb,
nadir, and zenith example inputs. It uses `formod TASK t`, whose reported time
measures repeated forward-model calls on perturbed atmospheres. The initial
lookup-table read and output write are outside that timer. Each invocation runs
for at least ten seconds of accumulated model time.

Build the CPU executable. Supply an external TRIA directory containing
`tria_CO2.nc`, `tria_H2O.nc`, `tria_O3.nc`, `tria_F11.nc`, and
`tria_CCl4.nc` (the nadir case uses only CO2). The scripts pass its
`tria` file prefix as `TBLBASE` with `TBLFMT 3`; no tables are copied
into the repository. The H2O table must have been generated with RFM
`H2O(sub)` for consistency with JURASSIC's MT_CKD 4.1 continuum. From the
repository root:

```sh
cd src && make -j && cd ..
JURASSIC_TBL_DIR=/path/to/tria-directory python3 projects/benchmark/run.py --threads 1 2 4 8
```

`--tbl-dir /path/to/tria-directory` can replace the environment variable.
The runner checks required files before starting and reports missing tables.
`--cases limb nadir zenith` selects cases; `--bin` selects an existing CPU
`formod` executable; `--output` selects the TSV result path. The default result
is `projects/benchmark/runs/cpu.tsv`. Logs and generated radiances are saved
beside it. The runner sets `OMP_NUM_THREADS` for each run and `LC_ALL=C` for
stable number parsing. Other OpenMP settings, such as `OMP_PROC_BIND` and
`OMP_PLACES`, are inherited from the environment; record their values and the
compiler, host, and commit when publishing measurements. Use the same machine
and build for scaling comparisons. The benchmark reports runtime and OpenMP scaling only; it makes no
accuracy claim. Use the same TRIA inventory as validation when comparing
workflow configurations.
