# Contributing to JURASSIC (Juelich Rapid Spectral Simulation Code)

Thank you for contributing to JURASSIC, an infrared radiative transfer and retrieval model for atmospheric remote sensing. Bug reports, documentation improvements, tests, and code contributions are welcome.

## Code of Conduct

All contributors are expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/slcs-jsc/jurassic/issues) with a clear description and, for bugs, steps to reproduce. Include relevant control files, commands, logs or output excerpts, plus the version or commit, operating system, compiler, and libraries.

## Submitting changes

1. Fork the repository, clone your fork, and create a focused branch.
2. Implement the change and update documentation and tests as appropriate.
3. Inspect the diff, run the relevant tests, and commit with a descriptive message.
4. Push your branch and open a pull request explaining the motivation, affected workflows, tests performed, and any intentional numerical or reference-data changes. Link related issues when applicable. The [pull request template](.github/pull_request_template.md) provides a checklist.

For example:

```bash
git diff --check
git status
cd src
make
make check
```

## Building and testing

Install the required libraries as described in [DEPENDENCIES.md](DEPENDENCIES.md). See [Building from source](docs/manual/building.md) for compilation and optional build modes.

The default regression suite runs with `make check`. Optional MPI and ARM tests require suitable environments. For the test catalog, reference-data provenance, comparison rules, optional checks, and CI, see [Testing and code quality](docs/manual/testing.md) and [tests/README.md](tests/README.md).

Do not update reference results solely to make a test pass. Document and justify intentional changes to numerical outputs. Regression tests are not a substitute for [scientific validation](docs/manual/validation.md).

## Coding standards

Follow the existing C style and naming conventions; avoid unrelated formatting. Keep builds free of compiler warnings under the project's strict warning settings. Document changed public interfaces in `src/jurassic.h` as appropriate, and consider MPI retrieval paths when changing shared data structures or retrieval workflows.

See [Code development](docs/manual/code-development.md) and [Code structure](docs/manual/code_structure.md) for further developer guidance.

## License

Contributions are licensed under the [GNU General Public License v3.0](COPYING).
