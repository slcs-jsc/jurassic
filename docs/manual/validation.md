# Validation and verification

This page summarizes how JURASSIC has been validated and how users can
verify correct installation and numerical behavior for their own
applications. Validation is a critical aspect of radiative transfer and
retrieval modelling and underpins the scientific credibility of the
results produced with JURASSIC.

---

## Validation philosophy

JURASSIC validation follows three complementary approaches:

1. **Intercomparison with reference models**
2. **Regression testing and example workflows**
3. **Scientific application and peer-reviewed publications**

Together, these approaches ensure that the model is both numerically
correct and scientifically reliable across a wide range of use cases.

---

## Intercomparison with reference models

The core radiative transfer algorithms implemented in JURASSIC have been
extensively benchmarked against established, high-accuracy reference
models, including:

- the Karlsruhe Optimized and Precise Radiative Transfer Algorithm (KOPRA),
- the Reference Forward Model (RFM),
- the Stand-alone AIRS Radiative Transfer Algorithm (SARTA).

These intercomparisons cover:

- limb and nadir viewing geometries,
- temperature and trace gas sensitivities,
- clear-sky and simplified aerosol conditions,
- a broad range of atmospheric states.

Results demonstrate that JURASSIC reproduces reference-model radiances,
Jacobians, and retrieval-relevant quantities within the accuracy
expected from its spectral approximations.

Published validation and application studies are listed on the
[References](references.md) page, including comparisons with KOPRA and
retrieval studies using satellite observations.

---

## Validation of spectral approximations

Particular attention has been given to validating the:

- Emissivity Growth Approximation (EGA),
- Curtis–Godson Approximation (CGA),
- band-averaged emissivity lookup table approach.

Comparisons against line-by-line calculations show that these
approximations provide a good balance between accuracy and performance
for typical infrared remote sensing applications, especially at
moderate spectral resolution.

Residual differences relative to line-by-line models are well
characterized and documented in the literature.

---

## Retrieval validation

The optimal estimation retrieval framework in JURASSIC has been
validated through:

- synthetic retrieval experiments,
- comparison with independent retrieval systems,
- application to real satellite measurements.

Key aspects of retrieval validation include:

- convergence behavior,
- consistency of Jacobians and averaging kernels,
- realism of retrieved error estimates,
- physical plausibility of retrieved atmospheric states.

Retrieval results obtained with JURASSIC have been published in numerous
peer-reviewed studies.

---

## Example projects and regression tests

The JURASSIC distribution includes example projects (e.g. limb, nadir,
and zenith configurations) that serve as both tutorials and regression
tests.

These examples:

- generate forward-model output for known configurations,
- compare results against reference data,
- produce diagnostic plots for visual inspection.

Running these examples after installation is the recommended way to
verify that JURASSIC is functioning correctly on a given system.

From the repository root, example projects are typically executed via:

```bash
cd projects/limb
./run.sh
```

and similarly for other configurations.

---

## Numerical validation workflow

Performance benchmarks and numerical validation are deliberately separate.
The runners in `projects/benchmark/` measure performance, whereas
`projects/validation/` compares candidate calculations with frozen reference data.
This prevents an optimized candidate from regenerating the data used to judge its
own accuracy.

The validation matrix is defined in
`projects/validation/configs/validation_cases.tsv`. The smoke profile covers limb,
nadir, and zenith geometries at 125 channels from 500 to 2980 cm⁻¹ in 20 cm⁻¹
steps. The full profile covers the complete 500 to 2999 cm⁻¹ range in 1 cm⁻¹ steps
and uses the full trace-gas set. Full validation assumes a build with sufficient gas
capacity, for example:

```bash
cd src
make clean
make DEFINES=-DNG=40 -j
```

Reference calculations are generated separately from a trusted scalar build and
stored as double-precision netCDF files. They should only be replaced deliberately;
normal candidate runs never modify them. From `projects/validation/`, run the quick
profile with:

```bash
VALIDATION_TBLBASE=/path/to/tria scripts/run_validation.py --profile smoke
```

Use `--profile full` for the complete matrix. The default forward-model method is
EGA. Other approximations can be evaluated against the same frozen references, for
example:

```bash
VALIDATION_TBLBASE=/path/to/tria \
  scripts/run_validation.py --profile smoke --formod-method cga
```

The pass/fail gate compares every brightness-temperature or radiance value and every
transmittance value directly. Paired non-finite values are accepted, one-sided
non-finite values fail, and each case uses the relative tolerance recorded in the
manifest. Per-case logs, netCDF candidate results, accuracy reports, and a combined
`summary.tsv` are written below `projects/validation/runs/`.

Generate channel-wise diagnostic plots for a completed run with:

```bash
scripts/plot_summary.py runs/validation_<timestamp>
```

This creates three complementary summaries:

- `summary_plot.png` is the main diagnostic: limb radiance uses relative error,
  while nadir/zenith brightness temperature and all transmittances use absolute
  error;
- `summary_plot_re.png` shows relative-error statistics for every panel;
- `summary_plot_ae.png` shows absolute-error statistics for every panel.

The plotted channel statistics are `MeanRE`, `SDRE`, `MinRE`, and `MaxRE`, or their
absolute-error counterparts `MeanAE`, `SDAE`, `MinAE`, and `MaxAE`. Here, “absolute”
means an error in the native output units rather than a relative percentage; the
values retain the sign of `simulation - reference`. Consequently, `MeanAE` measures
bias, and `MinAE` and `MaxAE` show the signed extrema. See the project-level
`projects/validation/README.md` for reference generation, provenance checks, custom
output paths, and individual-case plots.

---

## Automated test suite

In addition to example projects, the build system provides a test
target. From the `src/` directory, run:

```bash
make check
```

This target runs a predefined set of tests and reports pass/fail
status. These tests are intended to catch regressions caused by code
changes or build-system issues.

The default `make check` target currently runs the following suites from
`tests/`:

- `cli_test` - command-line behavior and help output for all executables
- `tools_test` - small utility programs such as time conversion and
  Planck/brightness-temperature tools
- `atm_test` - atmospheric-profile generation, interpolation, hydrostatic
  adjustment, and format conversion
- `mat_test` - matrix input/output and format conversion
- `obs_test` - observation geometry, ray tracing, and observation format
  conversion
- `tbl_test` - filter functions, lookup tables, and table format
  conversion
- `formod_test` - limb, nadir, and zenith forward-model calculations and
  kernels
- `ret_test` - retrieval workflow checks, including shared netCDF input
  and output

Most regression tests generate fresh output files and compare them with
checked-in reference files under the corresponding `data.ref/`
directory. This makes the test suite useful for detecting unintended
changes in numerical output, file-format behavior, and command-line
interfaces.

### `OBSREF` diagnostics

Supplying `OBSREF <file>` to `formod` prints two `EVAL` records per spectral
channel. The relative record reports the signed relative error statistics
`MeanRE`, `SDRE`, `MinRE`, and `MaxRE` in percent. The absolute record reports
`MeanAE`, `SDAE`, `MinAE`, and `MaxAE` in the active observation output units
(brightness temperature or radiance). Despite the `AE` label, these values use
the signed difference `simulation - reference`; `MeanAE` therefore measures
bias rather than mean absolute magnitude. `SDAE` is non-negative, while
`MinAE` and `MaxAE` retain the negative and positive extrema. Non-finite
simulation/reference pairs are excluded from both relative and absolute statistics.

An additional `mpi_test` suite is available for MPI-enabled retrieval
checks, but it is not part of the default `make check` target.

---

## Numerical reproducibility

Due to floating-point arithmetic and parallel execution, JURASSIC does
not guarantee bitwise-identical results across:

- different compilers,
- different MPI/OpenMP configurations,
- different hardware architectures.

The regression tests use exact file comparisons for controlled test
setups and supported reference outputs. For broader user workflows,
small numerical differences may occur across compilers, hardware, and
parallel configurations. Validation should therefore focus on physical
and statistical consistency rather than bitwise identity alone.

---

## User-level validation recommendations

Users are encouraged to perform application-specific validation by:

- comparing selected results against reference models or datasets,
- testing sensitivity to configuration parameters,
- inspecting Jacobians, averaging kernels, and residuals,
- documenting configuration choices and assumptions.

Such validation is especially important when introducing new lookup
tables, instrument configurations, or retrieval setups.

---

## Summary

JURASSIC has undergone extensive validation through model
intercomparisons, regression testing, and scientific application. These
efforts demonstrate that the model provides reliable and accurate
results within the scope of its documented assumptions and
approximations.

Users are encouraged to make use of the provided example projects and
tests to verify correct installation and to perform additional
application-specific validation as needed.

---

## Related pages

- [Performance considerations](performance.md)
- [Quickstart examples](quickstart.md)
- [Model assumptions & limitations](limitations.md)
