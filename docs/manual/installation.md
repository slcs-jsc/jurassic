# Installation

This section describes how to build and install JURASSIC on a Linux
system or high-performance computing (HPC) environment. JURASSIC is
designed for batch-style execution and is typically compiled from
source.

------------------------------------------------------------------------

## System requirements

JURASSIC is primarily developed and tested on 64-bit Linux systems,
including high-performance computing environments. It is built from
source using a Makefile-based workflow.

------------------------------------------------------------------------

## Required software

The following software components are required to build JURASSIC:

- **C compiler:** required for auxiliary components and libraries, for
  example `gcc`.
- **GNU Make:** required to build the bundled libraries and the JURASSIC
  executables.
- **GNU Scientific Library (GSL):** required for numerical kernels and
  linear algebra support.
- **netCDF-C library:** required for netCDF input/output support used
  throughout the toolchain.
- **HDF5 library:** required for netCDF-4/HDF5 input/output workflows.
- **MPI library (optional, retrieval only):** required only for
  MPI-enabled retrieval builds. Supported MPI implementations include
  OpenMPI and MPICH-derived implementations such as ParaStation MPI. MPI
  is used exclusively by the retrieval code to distribute independent
  retrieval tasks across processes.
- **Gnuplot (optional):** used by example projects to generate
  diagnostic plots.

------------------------------------------------------------------------

## Obtaining the source code

The JURASSIC source code is hosted in the
[GitHub repository](https://github.com/slcs-jsc/jurassic). To obtain the
most recent development version from the default branch, clone the
repository using:

```bash
git clone https://github.com/slcs-jsc/jurassic.git
cd jurassic
```

Alternatively, download a release archive from the
[GitHub releases page](https://github.com/slcs-jsc/jurassic/releases).

!!! note "Recommended version"
    We generally recommend using the current development version from the
    default branch. Bug fixes are usually applied there and are not
    routinely backported to older releases.

------------------------------------------------------------------------

## Configuring the build

JURASSIC uses a Makefile-based build system in `src/Makefile`. Prior to
compilation, you may need to edit that file or set make variables to match
your local compiler and MPI setup.

Typical configuration options include:

- Selection of the compiler preset via `COMPILER`
- Selection of the MPI wrapper via `MPICC` when building with `MPI=1`
- Compiler optimization and debugging flags
- Enabling or disabling MPI support (retrieval only)
- Enabling OpenMP parallelization
- Enabling optional GPU/OpenACC builds via `GPU=1`

On HPC systems, it is recommended to load the appropriate compiler and
MPI modules before configuring the build.

------------------------------------------------------------------------

## Building JURASSIC

The dependency build command is run from the `libs/` directory. The
JURASSIC build, install, and test commands are run from the `src/`
directory.

### Build bundled dependencies

From `libs/`, compile the bundled third-party libraries:

```bash
./build.sh
```

This populates `libs/build/`, which is the default include/library location
used by `src/Makefile`.

!!! note "Using system libraries"
    If you prefer to use system-provided GSL, HDF5, and netCDF libraries,
    this bundled-library step can be skipped, but the include and library
    paths in `src/Makefile` must point to the corresponding system
    locations.

### Default build (no MPI)

To build JURASSIC without MPI support:

```bash
make -j
```

This builds all executables in serial/OpenMP mode. MPI is not required
for this configuration.

---

### Building with MPI-enabled retrieval

To enable MPI support for the retrieval executable, build with:

```bash
make -j MPI=1
```

If your system provides multiple MPI wrappers, select the desired one explicitly:

```bash
make -j COMPILER=gcc MPI=1 MPICC=mpicc.openmpi
```

This will:

- compile the full tool suite with the MPI compiler wrapper,
- use `MPICC` (default: `mpicc`) as the MPI compiler wrapper,
- define the `MPI` preprocessor macro used by the retrieval source code.

MPI-specific runtime behavior is implemented only in the retrieval
executable. The other binaries are still rebuilt under the selected MPI
toolchain when `MPI=1` is used.

### Building with GPU support

To enable the GPU/OpenACC build path, select a supported compiler and build with:

```bash
make -j COMPILER=gcc GPU=1
make -j COMPILER=nvc GPU=1
```

For MPI-enabled GPU builds with the NVIDIA HPC SDK:

```bash
make -j COMPILER=nvc MPI=1 MPICC=mpicc GPU=1
```

Optional pinned-memory support for suitable `nvc` targets can be enabled with:

```bash
make -j COMPILER=nvc GPU=1 GPU_PIN=1
```

`clang` is supported for CPU builds, but the Makefile intentionally rejects
`COMPILER=clang GPU=1` because this configuration does not provide the required
OpenACC support.

---

### Clean rebuild

To perform a clean rebuild:

```bash
make clean
make -j
```

---

### Install binaries

The `src/Makefile` provides a simple copy-install target. By default,
executables are installed to `../bin` relative to the `src` directory:

```bash
make install
```

To install into a custom binary directory:

```bash
make install DESTDIR=/path/to/bin
```

To remove binaries from the same destination:

```bash
make uninstall DESTDIR=/path/to/bin
```

------------------------------------------------------------------------

## Verifying the installation

After compilation, verify the installation by running the test suite or the
example projects described in the [Quickstart](quickstart.md).

```bash
make check
```

A successful test run or example simulation indicates that JURASSIC has been
built correctly and that all required dependencies are working as expected.

------------------------------------------------------------------------

## Installation on HPC systems

On shared HPC systems, JURASSIC is typically installed in a user
workspace rather than system-wide. Recommended practices include:

- Building JURASSIC with the same compiler and MPI library used for
  production retrieval runs
- Enabling MPI only when running MPI-enabled retrievals
- Using environment modules to manage compiler and MPI versions
- Testing retrieval scalability with a small number of MPI ranks before
  large-scale production runs

Further details on MPI execution and performance considerations are
provided in [Running JURASSIC](running.md) and
[HPC workflows](hpc_workflows.md).

------------------------------------------------------------------------

## Troubleshooting

Common issues during installation include:

- Missing or incompatible compiler versions
- Using an MPI-enabled build without a working `MPICC` wrapper or MPI runtime
- Mismatches between compile-time and runtime MPI environments
- Incorrect OpenMP settings for the selected `COMPILER` preset
- Unsupported GPU/compiler combinations such as `COMPILER=clang GPU=1`
- Runtime linker errors caused by missing shared-library paths

If you encounter problems, consult the build output carefully and verify
that your compiler, library paths, and MPI environment are correctly
configured.
Additional help may be available through the project maintainers.

!!! warning "Bundled shared libraries"
    When using the bundled dynamic libraries, executables may need the
    library directory in `LD_LIBRARY_PATH` at runtime:

    ```bash
    export LD_LIBRARY_PATH=/path/to/jurassic/libs/build/lib:$LD_LIBRARY_PATH
    ```

    The example scripts in `projects/` set this path automatically for
    the default bundled-library layout.

------------------------------------------------------------------------

## Next steps

Once JURASSIC is installed, proceed to the [Quickstart](quickstart.md)
to run your first simulation, or consult the User Manual for detailed
information on configuration and usage.
