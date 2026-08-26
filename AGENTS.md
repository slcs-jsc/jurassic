# Repository Guidelines

## Project Structure & Module Organization
Core source code lives in `src/`. The main library and command-line tools are built from C sources there, with outputs installed to `bin/` by default via `make install`. Regression tests live under `tests/`, with one directory per suite such as `tests/atm_test`, `tests/formod_test`, and `tests/tools_test`; each suite is driven by a `run.sh` script plus control or reference data files. Example workflows are in `projects/examples/limb`, `projects/examples/nadir`, and `projects/examples/zenith`. Documentation sources are under `docs/`, and bundled third-party dependency build scripts live in `libs/`.

## Build, Test, and Development Commands
From the repository root, build bundled libraries first when needed: `cd libs && bash build.sh`. Build the code with `cd src && make -j`; this compiles the executables listed in `src/Makefile`. Run the regression suite with `cd src && make check`, which executes each `../tests/<name>/run.sh` script. Use `cd src && make clean` to remove binaries and coverage files. For analysis and docs, use `cd src && make cppcheck`, `make lizard`, `make doxygen`, or `make mkdocs`. Generate coverage with `cd src && make COV=1 && make check && make coverage`.

## Coding Style & Naming Conventions
This is a C codebase compiled with strict warnings enabled (`-Werror -Wall -Wconversion ...`). Keep changes warning-free with GCC. Match the existing formatting style: two-space indentation, opening braces on the same line for functions and control statements, and concise block comments such as `/* Allocate... */`. Use `snake_case` for functions and variables, and keep new executable or test names aligned with existing patterns like `obsfmt`, `tblgen`, and `atm_test`. If you reformat, use `cd src && make indent`.

## Testing Guidelines
Add or update regression coverage in the nearest existing suite under `tests/`, or create a new `<feature>_test` directory with a `run.sh` entry and reference outputs. Keep test scripts POSIX shell or Bash, fail fast with `set -euo pipefail`, and write deterministic outputs to local `data/` subdirectories before diffing against checked-in references. Run `cd src && make check` before opening a PR; use the coverage workflow when touching numerical kernels or I/O-heavy paths.
When public I/O APIs or user-visible file-format behavior change, update both the Doxygen comments in `src/` and the MkDocs manual in `docs/manual/` as part of the same change.

## Commit & Pull Request Guidelines
Recent history uses short, imperative commit subjects such as `Update workflow.` and `Fixed doxygen comments.` Keep subjects concise and specific. Pull requests should describe the scientific or tooling impact, note any changed commands or dependencies, and link related issues. Include sample output or plots when modifying example projects, docs rendering, or user-facing diagnostics. Confirm that GitHub Actions `tests` passes before requesting review.
