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
- C and Fortran compilers with OpenMP support
- MPI library (optional, retrieval only)
- GNU Make or a compatible build system

MPI is required only if MPI-enabled retrieval executables are built.
All other components can be built and run without MPI.

------------------------------------------------------------------------

## Required software

The following software components are required to build JURASSIC:

- **Fortran compiler**  
  A modern Fortran compiler such as:

  - GNU Fortran (`gfortran`)
  - Intel oneAPI Fortran (`ifx` / `ifort`)
  - NVHPC Fortran (`nvfortran`)

- **C compiler**  
  Required for auxiliary components and libraries (e.g. `gcc`, `icc`).

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

JURASSIC uses a Makefile-based build system. Prior to compilation, you
may need to edit the Makefile or set make variables to match your local
compiler and MPI setup.

Typical configuration options include:

- Selection of the Fortran and C compilers
- Compiler optimization and debugging flags
- Enabling or disabling MPI support (retrieval only)
- Enabling OpenMP parallelization

On HPC systems, it is recommended to load the appropriate compiler and
MPI modules before configuring the build.

------------------------------------------------------------------------

## Building JURASSIC

### Default build (no MPI)

To build JURASSIC without MPI support:

```bash
make
```

This builds all executables in serial/OpenMP mode. MPI is not required
for this configuration.

---

### Building with MPI-enabled retrieval

To enable MPI support for the retrieval executable, build with:

```bash
make MPI=1
```

This will:

- compile the retrieval code with MPI support,
- automatically select `mpicc` (unless `CC` is set explicitly),
- define the `MPI` preprocessor macro used by the retrieval source code.

All other executables remain non-MPI and are unaffected by this option.

---

### Clean rebuild

To perform a clean rebuild:

```bash
make clean
make
```

------------------------------------------------------------------------

## Verifying the installation

After compilation, verify the installation by running the example
projects described in the [Quickstart](quickstart.md).

A successful run of the example simulations indicates that JURASSIC has
been built correctly and that all required dependencies are working as
expected.

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
