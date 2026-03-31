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

Example projects are typically executed via:

```bash
cd projects/limb
./run.sh
```

and similarly for other configurations.

---

## Automated test suite

In addition to example projects, the build system provides a test
target:

```bash
make check
```

This target runs a predefined set of tests and reports pass/fail
status. These tests are intended to catch regressions caused by code
changes or build-system issues.

---

## Numerical reproducibility

Due to floating-point arithmetic and parallel execution, JURASSIC does
not guarantee bitwise-identical results across:

- different compilers,
- different MPI/OpenMP configurations,
- different hardware architectures.

However, numerical differences are typically small and do not affect
scientific conclusions. Validation should therefore focus on physical
and statistical consistency rather than exact numerical identity.

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
- [Testing and examples](quickstart.md)
- [Model assumptions & limitations](limitations.md)
