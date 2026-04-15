# Building from source

This page describes how to build the JURASSIC command-line tools from
source using the provided `Makefile`. It is intended for developers and
advanced users who want to compile JURASSIC themselves, change compiler
flags, or enable optional build modes.

Unless stated otherwise, the `make` commands on this page are run from
the `src/` directory of the JURASSIC source tree.

!!! tip "First run"
    If you are only looking for a first run, see the
    [Quickstart](quickstart.md). For general setup and prerequisites, see
    [Installation](installation.md).

---

## What the Makefile builds

The build in `src/` produces a set of small executables (command-line
tools) that are linked with the JURASSIC object code and the dependency
libraries in `../libs/build/`. The Makefile defines the executable list
in the `EXC` variable.

Typical tools include drivers for:

- forward modelling (e.g. `formod`)
- kernels / Jacobians (e.g. `kernel`)
- retrieval (e.g. `retrieval`)
- format conversion utilities (e.g. `atmfmt`, `obsfmt`, `matfmt`)
- table utilities (`tblgen`, `tblfmt`, ...)

The exact list can be inspected in the Makefile via the `EXC = ...` line.

---

## Default build

From the directory that contains the Makefile, run:

```bash
make
```

This builds the default target `all` (i.e., all executables).

For a faster build on multi-core systems:

```bash
make -j
```

---

## Build outputs

By default, the resulting executables are placed in the current build
directory (the same directory as the Makefile).

A version string is compiled into the code via:

- `VERSION ?= $(shell git describe --abbrev=6 --dirty --always --tags)`

This means builds from a Git checkout will embed a descriptive version
(e.g. tag + commit hash), and "dirty" working trees are marked.

---

## Running the test suite

The Makefile defines a test list in `TESTS` and provides a `check` target
that runs them:

```bash
make check
```

Each test is executed by running the corresponding `../tests/<name>/run.sh`
script (the Makefile prints a pass/fail banner for each test).

---

## Installing and uninstalling

The Makefile supports a simple "copy install" into a binary directory:

```bash
make install
```

The default install destination is:

- `DESTDIR ?= ../bin`

To install into a custom location:

```bash
make install DESTDIR=/path/to/bin
```

To uninstall from the same destination:

```bash
make uninstall DESTDIR=/path/to/bin
```

---

## Important build variables

You can override Makefile variables on the command line. For most local
builds, prefer the high-level variables below over replacing `CFLAGS`
directly.

### Compiler and flags

- `CC` - C compiler (default: `gcc`)
- `CFLAGS` - compile flags; overriding this variable replaces the
  Makefile defaults, including include paths, preprocessor definitions,
  warning flags, and OpenMP flags
- `LDFLAGS` - link flags

Examples:

```bash
make CC=clang
make OPT="-O2 -march=native"
```

### Optimization and diagnostics flags

- `OPT` (default: `-O3`) - optimization flags passed into `CFLAGS`
- `INFO` (default: `0`) - if set to `1`, adds `-fopt-info` for optimization reports
- `PROF` (default: `0`) - if set to `1`, adds `-pg` for profiling builds
- `COV` (default: `0`) - if set to `1`, enables coverage flags (`--coverage`)

Examples:

```bash
make OPT=-O0
make INFO=1
make PROF=1
make COV=1
```

### MPI support

- `MPI` (default: `0`) - if set to `1`, builds with MPI support enabled
  and defines the `MPI` preprocessor macro

```bash
make MPI=1
```

When `MPI=1`, the Makefile uses `mpicc` unless `CC` is set explicitly.
MPI-specific runtime behavior is implemented only in the retrieval
executable, but the full tool suite is rebuilt with the selected MPI
compiler wrapper.

### Static linking

- `STATIC` (default: `0`) - if set to `1`, attempts to add `-static`

```bash
make STATIC=1
```

Note: the Makefile explicitly prevents static builds when `MPI=1` or
`UNIFIED=1`.

### Library paths

The build expects headers and libraries under:

- includes: `../libs/build/include` (added via `INCDIR`)
- libs: `../libs/build/lib` (added via `LIBDIR`)

If these directories are missing, ensure the dependency libraries have been
built/installed as described in [Installation](installation.md).

---

## Optional: build with the JURASSIC-UNIFIED library

The Makefile contains an option to link against a "JURASSIC-UNIFIED"
library (including CUDA runtime linkage). This is controlled by:

- `UNIFIED ?= 0` (default off)
- `UNIDIR = ../../jurassic-unified` (path to the unified repository)

To enable:

```bash
make UNIFIED=1
```

When `UNIFIED=1`, the Makefile:

- adds include paths from `$(UNIDIR)/include`
- adds library paths from `$(UNIDIR)/unified_library` and `$(CUDA_PATH)/lib64`
- links against `-lcudart -ljurassic_unified` (and `-lstdc++`)

You may need to define `CUDA_PATH` in your environment, for example:

```bash
export CUDA_PATH=/usr/local/cuda
make UNIFIED=1
```

---

## Developer utility targets

The Makefile provides several targets that are helpful for development:

- `make clean` - remove build products (executables, objects, coverage files)
- `make doxygen` - build local Doxygen documentation; the hosted
  [Doxygen manual](https://slcs-jsc.github.io/jurassic/doxygen/) is also
  available online
- `make mkdocs` - build MkDocs documentation (runs `mkdocs build` in `../docs`)
- `make cppcheck` - run static analysis (requires `cppcheck`)
- `make indent` - apply code formatting rules (requires `indent`)
- `make lizard` - run cyclomatic complexity analysis (requires `lizard`)
- `make coverage` - generate coverage output (requires `lcov`, `genhtml`, `gcov`)
- `make strip` - strip symbols from executables
- `make dist` - create a source distribution artifact (if used in your workflow)

Not all developer tools are required for normal builds. If a target fails
because a tool is missing, either install the tool or skip that target.

---

## Troubleshooting

### "Cannot find headers/libraries under ../libs/build/..."

This typically means dependency libraries have not been built or the build
directory layout differs from what the Makefile expects. Ensure your
dependency build step populates:

- `../libs/build/include`
- `../libs/build/lib`

### MPI / OpenMP issues

The Makefile enables OpenMP via `-fopenmp` in `CFLAGS` by default (for GCC).
If you use a different compiler, you may need compiler-specific OpenMP flags.

### CUDA / unified build issues

If `UNIFIED=1` fails:

- confirm `UNIDIR` points to the correct checkout of `jurassic-unified`
- ensure `CUDA_PATH` is set and points to a valid CUDA installation
- confirm the unified library has been built before building this repository

---

## Next steps

After a successful build:

- run `make check` to validate basic functionality
- proceed to the [Quickstart](quickstart.md) examples
- consult [Code structure](code_structure.md) if you plan to modify or extend the code
