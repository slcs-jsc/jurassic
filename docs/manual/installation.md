# Installation

This section describes how to build and install JURASSIC on a Linux
system or high-performance computing (HPC) environment. JURASSIC is
designed for batch-style execution and is typically compiled from
source.

------------------------------------------------------------------------

## System requirements

JURASSIC is primarily developed and tested on Linux systems. The
following requirements apply:

- 64-bit Linux operating system
- C compiler with OpenMP support
- GNU Scientific Library (GSL)
- netCDF-C library
- MPI library (optional, retrieval only)
- GNU Make or a compatible build system

MPI is required only if MPI-enabled retrieval executables are built.
All other components can be built and run without MPI.

------------------------------------------------------------------------

## Required software

The following software components are required to build JURASSIC:

- **C compiler**  
  Required for auxiliary components and libraries (e.g. `gcc`, `icc`).

- **GNU Make**  
  Required to build the bundled libraries and the JURASSIC executables.

- **GNU Scientific Library (GSL)**  
  Required for numerical kernels and linear algebra support.

- **netCDF-C library**  
  Required for netCDF input/output support used throughout the toolchain.

- **MPI library** (optional, retrieval only)  
  For example:

  - OpenMPI
  - MPICH
  - Intel MPI

  MPI is used exclusively by the retrieval code to distribute
  independent retrieval tasks across processes. No other JURASSIC
  executables use MPI internally.

- **GNU Plot** (optional)  
  Used by example projects to generate diagnostic plots.

- **Optional Fortran compiler**  
  A Fortran compiler is not required for the default repository build. It may
  be useful on HPC systems or when integrating external tooling that depends on
  Fortran-enabled libraries.

------------------------------------------------------------------------

## Obtaining the source code

The JURASSIC source code is hosted on GitHub. Clone the repository
using:

```bash
git clone https://github.com/slcs-jsc/jurassic.git
cd jurassic
```

Alternatively, you may download a source archive from the GitHub
repository.

------------------------------------------------------------------------

## Configuring the build

JURASSIC uses a Makefile-based build system in `src/Makefile`. Prior to
compilation, you may need to edit that file or set make variables to match
your local compiler and MPI setup.

Typical configuration options include:

- Selection of the C compiler and any optional external toolchains
- Compiler optimization and debugging flags
- Enabling or disabling MPI support (retrieval only)
- Enabling OpenMP parallelization

On HPC systems, it is recommended to load the appropriate compiler and
MPI modules before configuring the build.

------------------------------------------------------------------------

## Building JURASSIC

### Build bundled dependencies

For the default repository build, first compile the bundled third-party
libraries:

```bash
cd libs
bash build.sh
```

This populates `libs/build/`, which is the default include/library location
used by `src/Makefile`.

### Default build (no MPI)

To build JURASSIC without MPI support:

```bash
cd src
make
```

This builds all executables in serial/OpenMP mode. MPI is not required
for this configuration.

---

### Building with MPI-enabled retrieval

To enable MPI support for the retrieval executable, build with:

```bash
cd src
make MPI=1
```

This will:

- compile the full tool suite with the MPI compiler wrapper,
- automatically select `mpicc` (unless `CC` is set explicitly),
- define the `MPI` preprocessor macro used by the retrieval source code.

MPI-specific runtime behavior is implemented only in the retrieval
executable. The other binaries are still rebuilt under the selected MPI
toolchain when `MPI=1` is used.

---

### Clean rebuild

To perform a clean rebuild:

```bash
cd src
make clean
make
```

------------------------------------------------------------------------

## Verifying the installation

After compilation, verify the installation by running the test suite or the
example projects described in the [Quickstart](quickstart.md).

```bash
cd src
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
provided in the HPC workflows documentation.

------------------------------------------------------------------------

## Troubleshooting

Common issues during installation include:

- Missing or incompatible compiler versions
- Using an MPI-enabled build without `mpicc` or an MPI runtime
- Mismatches between compile-time and runtime MPI environments
- Incorrect OpenMP settings

If you encounter problems, consult the build output carefully and verify
that your compiler and MPI environment are correctly configured.
Additional help may be available through the project maintainers.

------------------------------------------------------------------------

## Next steps

Once JURASSIC is installed, proceed to the [Quickstart](quickstart.md)
to run your first simulation, or consult the User Manual for detailed
information on configuration and usage.
