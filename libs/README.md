# JURASSIC third-party libraries

This directory contains bundled source archives of the third-party
libraries used by JURASSIC and a [`build.sh`](build.sh) script that
builds and installs them locally.

## Bundled sources

| Archive | Version | Purpose |
| --- | --- | --- |
| `gsl-2.7.1.tar.bz2` | 2.7.1 | GNU Scientific Library for numerical calculations |
| `zlib-1.3.1.tar.bz2` | 1.3.1 | Compression dependency of the HDF5/netCDF stack |
| `szip-2.1.1.tar.bz2` | 2.1.1 | Additional HDF5 compression support |
| `hdf5-1.14.4-3.tar.bz2` | 1.14.4-3 | Hierarchical scientific data storage |
| `netcdf-c-4.9.2.tar.bz2` | 4.9.2 | netCDF input/output |

The archives are included in the repository, so these sources do not
have to be downloaded separately. The build still requires suitable
compilers, `make`, and the usual tools needed by the bundled source
packages.

## Building

Run the script **from the `libs/` directory**:

```bash
cd libs
./build.sh
```

The script removes any previous `libs/build/` directory, extracts
the source archives, and builds all five libraries in dependency order:
GSL, zlib, SZIP, HDF5, and netCDF-C. There are no command-line flags
for selecting individual libraries.

The local installation prefix is `libs/build/`, with generated
`bin/`, `lib/`, `src/`, and other installation directories.
The script sets `LD_LIBRARY_PATH` for the build process. Depending on
the host environment, the installed shared libraries may also need to
be made available when running JURASSIC.

HDF5 is configured with its high-level library, without its Fortran
or C++ interfaces. netCDF-C uses the locally built HDF5 installation.
The build script runs upstream test suites; if the HDF5 or netCDF test
suite reports failures, it prints a warning and continues with
installation, whereas a failed build or installation stops the script.

**Note:** Each invocation starts by deleting `libs/build/`. Do not
use this script expecting it to preserve an existing local installation.

For building the JURASSIC programs after installing the dependencies,
see the repository's main build instructions and
[`src/Makefile`](../src/Makefile).
