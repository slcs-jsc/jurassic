# JURASSIC Dependencies

This file lists the software dependencies required to build and run JURASSIC.
Versions correspond to the packages provided by Ubuntu 24.04 LTS (Noble Numbat),
which is the reference platform used for development and testing.

## Mandatory Dependencies

| Dependency | Version  | Description |
|------------|----------|-------------|
| Git | 2.43.0 | Clone the repository and derive version info |
| GNU Make | 4.3 | Build system |
| GCC | 13.3.0 | C compiler |
| GNU Scientific Library | 2.7.1 | Numerical routines |
| netCDF-C | 4.9.2 | File I/O (netCDF) |

## Optional Dependencies

| Dependency | Version | Description |
|------------|---------|-------------|
| HDF5 | 1.14.4 | Required for netCDF-4 support |
| zlib | 1.3.1 | Required for netCDF-4 compression |
| Zstandard (zstd) | 1.5.5 | Optional netCDF/HDF5 filter-based compression |
| gnuplot | 6.0 | Used by example scripts |
| Cppcheck | 2.13 | Static code analysis  |
| gprofng | 2.42 | Performance profiling tool |
| gcov | 13.2.0 | Coverage testing tool |
| lcov | 2.0 | Code coverage visualization tool |
| Doxygen | 1.9.8 | Automatic source code documentation |
| MkDocs | 1.5.3 | Build documentation locally |

## Notes

- **Using bundled libraries:** JURASSIC includes a `libs/` directory that can build GSL and netCDF if system libraries are not available or if a self-contained build is preferred.

- **Bundled build toolchain:** `libs/build.sh` builds HDF5 with the C++ interface disabled (`--disable-cxx`). A GNU C compiler is required, but a C++ compiler is not required for the bundled-library build.

- **netCDF version:** netCDF classic only requires `libnetcdf-dev`. netCDF-4 additionally requires HDF5 and zlib. Using filter-based compression (e.g. zstd) requires the HDF5 plugin package and possibly setting `HDF5_PLUGIN_PATH`.

- **Static linking**: Fully static builds require static versions of all libraries and usually do not support HDF5 filter plugins. Dynamic linking is recommended for most users.

## Installing dependencies

### Ubuntu 24.04 LTS

Use the following commands to install the dependencies on Ubuntu 24.04 LTS.

```bash
sudo apt update
sudo apt install \
  git make gcc \
  libgsl-dev \
  libnetcdf-dev
```

Optional netCDF-4 and compression support:

```bash
sudo apt install \
  libhdf5-dev zlib1g-dev libzstd-dev
```

Optional tools:

```bash
sudo apt install \
  binutils cppcheck lcov \
  doxygen mkdocs \
  gnuplot
```

### Fedora

Use the following commands to install the dependencies on Fedora.

```bash
sudo dnf install \
  git make gcc \
  gsl-devel \
  netcdf-devel
```

Optional netCDF-4 and compression support:

```bash
sudo dnf install \
  hdf5-devel zlib-devel libzstd-devel
```

Optional tools:

```bash
sudo dnf install \
  binutils cppcheck lcov \
  doxygen mkdocs \
  gnuplot
```
