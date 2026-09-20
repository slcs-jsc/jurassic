# CPU/OpenMP benchmark

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
