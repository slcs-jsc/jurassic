# Contributing to JURASSIC (Juelich Rapid Spectral Simulation Code)

Thank you for your interest in contributing to JURASSIC, an infrared radiative transfer and retrieval model for atmospheric remote sensing. Bug reports, documentation improvements, tests, and code contributions are welcome.

## Code of Conduct

All contributors are expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## How to Contribute

### Reporting Issues

Please open a [GitHub issue](https://github.com/slcs-jsc/jurassic/issues) with:

- A clear description and, for bugs, steps to reproduce the problem.
- Relevant control files, commands, logs, error messages, or output excerpts.
- The JURASSIC version or commit, operating system, compiler, and relevant libraries.

### Submitting Changes

1. Fork the repository and clone your fork:

   ```bash
   git clone https://github.com/your-username/jurassic.git
   cd jurassic
   ```

2. Create a feature branch:

   ```bash
   git checkout -b your-branch-name
   ```

3. Make focused changes and update documentation and tests where appropriate.
4. Run the relevant tests described below. Inspect the changes before committing:

   ```bash
   git diff --check
   git status
   git add <changed-files>
   git commit -m "Describe the change"
   ```

5. Push your branch and open a pull request against the JURASSIC repository:

   ```bash
   git push origin your-branch-name
   ```

   Explain the motivation, affected workflows, tests performed, and any intentional changes to numerical results or reference data. Link related issues when applicable. The [pull request template](.github/pull_request_template.md) provides a checklist.

## Building and Testing

Install the required libraries as described in [DEPENDENCIES.md](DEPENDENCIES.md). If they are unavailable on your system, the bundled source archives can be built locally using `cd libs && bash build.sh`. For additional build options, see the [building guide](docs/manual/building.md).

From the repository root, build JURASSIC and run the default regression suite:

```bash
cd src
make
make check
```

The default `make check` runs nine test suites, including forward-model and retrieval tests. The same suites are exercised in GitHub Actions. You can also run a single suite, for example:

```bash
make formod_test
make ret_test
```

Optional, environment-specific tests are intentionally separate from `make check`:

- `make mpi_test` checks MPI retrieval results against serial results; it requires an appropriate MPI build and runtime (for example, configure the build with `MPI=1`).
- `make arm_test` uses ARM-specific forward-model reference files and is intended for a compatible ARM environment.

After changing compiler or build options, use `make clean` before rebuilding: the Makefile does not automatically detect changes to command-line flags. See [tests/README.md](tests/README.md) for the full test catalog and instructions.

### Regression Tests and Reference Data

Most suites generate outputs under their local `data/` directory and compare them **exactly** with checked-in `data.ref/` files using `diff`. This strict comparison is intentional: it detects unexpected changes to numerical results and file formats in controlled configurations. It is not a substitute for scientific accuracy validation, which is addressed separately in [the validation documentation](docs/manual/validation.md) and `projects/validation/`.

When a change affects a test result:

1. Run the relevant test and inspect the differences against `data.ref/`.
2. Determine whether the change is expected and scientifically or technically justified; do not update references merely to make a test pass.
3. Update **only** the affected reference files after verification.
4. Include the reference updates with the corresponding code change and explain their rationale in the commit or pull request.

Reference data provenance and further details are documented in [tests/README.md](tests/README.md).

## Coding Standards

- Follow the existing C style and naming conventions; avoid unrelated formatting changes.
- Keep builds free of compiler warnings under the project's strict warning settings, and document new functions and interfaces in `src/jurassic.h` as appropriate.
- Add or adjust the most relevant regression tests when changing behavior, and run `make check` before requesting review.
- For changes to shared data structures or retrieval workflows, consider their use in MPI paths.
- Keep changes focused, and update user documentation for changed functionality or configuration.

## License

By contributing to JURASSIC, you agree that your contributions will be licensed under the GNU General Public License v3.0 (see [COPYING](COPYING)).
