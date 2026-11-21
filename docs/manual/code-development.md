# Code development

This page provides an entry point for developers contributing to JURASSIC.

## Development workflow

Follow the existing C style and naming conventions. Keep changes focused
and avoid unrelated formatting changes. New or modified public functions
and data structures should be documented in `src/jurassic.h` as appropriate;
the [Doxygen reference](https://slcs-jsc.github.io/jurassic/doxygen/)
provides source-level API documentation.

Before committing, inspect the diff and run the relevant tests:

```bash
git diff --check
git status
```

Use the project's [contribution guidelines](https://github.com/slcs-jsc/jurassic/blob/master/CONTRIBUTING.md)
for issue reports, feature branches, pull requests, and review expectations.
Small maintainer documentation changes may be committed directly to the
default branch; substantial or externally contributed changes generally
benefit from review.

## Formatting, analysis, and documentation

The `src/Makefile` provides developer utilities, including:

```bash
cd src
make cppcheck
make lizard
make indent
make mkdocs
make doxygen
```

The `indent` target changes source formatting; review its diff carefully
and avoid broad unrelated reformatting. Static analysis and documentation
targets may require optional tools to be installed.

Update the relevant MkDocs page when behavior, configuration, or user
workflows change. For tests, reference data, CI, and coverage, see
[Testing and code quality](testing.md).

## Developer manual

- [Code structure](code_structure.md): C library, public API, and tools.
- [Building from source](building.md): compilation and build options.
- [Testing and code quality](testing.md): regression tests, reference
  data, automated checks, and static analysis.
- [Profiling and benchmarking](profiling.md): performance investigations.
- [Release management](releases.md): coordinated release policy and checklist.

For scientific validation beyond regression tests, see
[Validation and verification](validation.md).
