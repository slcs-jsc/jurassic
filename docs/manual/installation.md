# Installation

This section describes how to build and install JURASSIC on a Linux
system or high-performance computing (HPC) environment. JURASSIC is
designed for batch-style execution and is typically compiled from
source.

------------------------------------------------------------------------

## System requirements

JURASSIC is primarily developed and tested on Linux systems. The
following requirements apply:

-   64-bit Linux operating system
-   C and Fortran compilers with OpenMP support
-   MPI library for distributed-memory parallel execution (optional but
    recommended)
-   GNU Make or a compatible build system

The model is commonly deployed on HPC clusters, but it can also be run
on modern workstations.

------------------------------------------------------------------------

## Required software

The following software components are required to build JURASSIC:

-   **Fortran compiler**\
    A modern Fortran compiler such as:

    -   GNU Fortran (`gfortran`)
    -   Intel oneAPI Fortran (`ifx` / `ifort`)
    -   NVHPC Fortran (`nvfortran`)

-   **C compiler**\
    Required for auxiliary components and libraries (e.g. `gcc`, `icc`).

-   **MPI library** (optional)\
    For example:

    -   OpenMPI
    -   MPICH
    -   Intel MPI

    If MPI is not available, JURASSIC can still be built and run in
    shared-memory (OpenMP-only) mode.

-   **GNU Plot** (optional)\
    Used by the example projects to generate diagnostic plots.

------------------------------------------------------------------------

## Obtaining the source code

The JURASSIC source code is hosted on GitHub. Clone the repository
using:

``` bash
git clone https://github.com/slcs-jsc/jurassic.git
cd jurassic
```

Alternatively, you may download a source archive from the GitHub
repository.

------------------------------------------------------------------------

## Configuring the build

JURASSIC uses a Makefile-based build system. Prior to compilation, you
may need to edit the configuration files to match your local compiler
and MPI setup.

Typical configuration options include:

-   Selection of the Fortran compiler
-   Compiler optimization and debugging flags
-   Enabling or disabling MPI support
-   Enabling OpenMP parallelization

On HPC systems, it is recommended to load the appropriate compiler and
MPI modules before configuring the build.

------------------------------------------------------------------------

## Building JURASSIC

To compile JURASSIC, run:

``` bash
make
```

This will build the JURASSIC executable and supporting libraries.

On successful completion, the main executable will be available in the
project directory.

To perform a clean rebuild, use:

``` bash
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

-   Building JURASSIC with the same compiler and MPI library used for
    production runs
-   Using environment modules to manage compiler and MPI versions
-   Testing scalability with a small number of MPI ranks before
    large-scale production runs

Further details on parallel execution and performance tuning are
provided in the Advanced Usage section of the User Manual.

------------------------------------------------------------------------

## Troubleshooting

Common issues during installation include:

-   Missing or incompatible compiler versions
-   MPI library mismatches between compile time and runtime
-   Incorrect OpenMP settings

If you encounter problems, consult the build output carefully and verify
that your compiler and MPI environment are correctly configured.
Additional help may be available through the project maintainers.

------------------------------------------------------------------------

## Next steps

Once JURASSIC is installed, proceed to the [Quickstart](quickstart.md)
to run your first simulation, or consult the User Manual for detailed
information on configuration and usage.
