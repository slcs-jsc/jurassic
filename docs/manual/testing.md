# Testing and code quality

This page summarizes the development checks for JURASSIC. For the full
test catalog, input provenance, and reference-data rules, see
[tests/README.md](https://github.com/slcs-jsc/jurassic/blob/master/tests/README.md).
For scientific model assessment, see [Validation and verification](validation.md).

## Regression tests

From the repository root:

```bash
cd src
make
make check
```

The default `make check` executes nine suites: `cli_test`, `tools_test`,
`clim_test`, `atm_test`, `mat_test`, `obs_test`, `tbl_test`,
`formod_test`, and `ret_test`. For targeted checks, use the corresponding
Makefile target, for example:

```bash
make formod_test
make ret_test
```

The optional `make mpi_test` checks MPI retrievals against serial results
and requires a suitable MPI build and runtime. The optional `make arm_test`
uses ARM-specific forward-model references on a compatible ARM system.
Neither is included in the default `make check`.

After changing compiler or build options, use `make clean` before rebuilding:
the Makefile does not automatically detect changes to command-line flags.

## Reference results

Most test scripts write generated results to their own `data/` directory
and compare them exactly with checked-in `data.ref/` files. Some tests
compare files generated within the same run instead. Exact regression
comparisons help detect changes to numerical outputs and file formats;
they do not prove scientific accuracy across platforms or applications.

If a regression result changes:

1. Inspect the generated output and the corresponding reference.
2. Determine whether the change is expected and scientifically or
   technically justified.
3. Update only the affected reference files, after verification.
4. Explain the change and its impact in the commit or pull request.

Do not update reference data simply to make a failing test pass.
See [tests/README.md](https://github.com/slcs-jsc/jurassic/blob/master/tests/README.md)
for provenance and limitations of the reference inputs.

## Static analysis and coverage

The build system provides optional developer targets:

```bash
cd src
make cppcheck
make lizard
make COV=1
make coverage
```

These require the corresponding external tools. Coverage instrumentation
may require a clean rebuild when switching build flags. Use coverage
results to identify untested code paths, not as a substitute for
scientific validation.

## Automated checks

Review [GitHub Actions](https://github.com/slcs-jsc/jurassic/actions)
after pushing changes. Additional information is available from
[Codecov](https://codecov.io/gh/slcs-jsc/jurassic),
[Codacy](https://app.codacy.com/gh/slcs-jsc/jurassic), and the
[JURASSIC nightly builds](https://datapub.fz-juelich.de/slcs/jurassic/nightly_builds).
Nightly results are particularly useful for configurations not exercised
on a developer workstation.

## Reproducibility and scientific validation

Exact comparisons are meaningful for the controlled reference
configurations, but floating-point results may differ across compilers,
hardware, and MPI/OpenMP configurations. Investigate changes in
physical or statistical consistency as well as test failures.
Independent reference-model comparisons and application-specific
validation are described in [Validation and verification](validation.md)
and the repository's `projects/validation/` directory.
