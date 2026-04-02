# Code structure

This page provides a developer-oriented overview of the **C library layer**
of JURASSIC and the small command-line applications built on top of it.

The intent is to help new contributors understand:

- where the *public API* lives,
- how *data structures* (control, atmosphere, observation) are represented,
- how the forward model, Jacobians, and retrieval are wired together, and
- how to add new small tools that reuse the JURASSIC library.

---

## Core library files

### `jurassic.h` — public API and core data types

`jurassic.h` is the central header that defines the main data structures
and declares the public functions used throughout the code base.

Key data structures (typedefs):

- `ctl_t` — **control parameters** (model switches, spectral settings, runtime options)
- `atm_t` — **atmospheric state** (profile data, optional cloud/surface fields, etc.)
- `obs_t` — **observation geometry + radiance container** (viewing geometry and per-channel radiance/Jacobian outputs)
- `los_t` — **line-of-sight representation** used internally during ray tracing / integration
- `ret_t` — **retrieval control parameters** (optimal estimation settings, convergence, regularization)
- `tbl_t` — **lookup-table container** (precomputed spectroscopy / emissivity/transmittance tables)

The header is Doxygen-friendly and many functions/types have docstrings.
If you are looking for “what is the intended contract of this function?”,
start with the corresponding Doxygen block in `jurassic.h`.

### `jurassic.c` — implementation

`jurassic.c` implements the functionality declared in `jurassic.h` and
contains most of the numerical machinery used by the example tools:

- input/output helpers (`read_*`, `write_*`)
- forward modelling (`formod*`)
- Jacobians / kernels (`kernel`, plus operator helpers)
- retrieval / optimal estimation (`optimal_estimation`, minimization utilities)
- LUT handling (`read_tbl_*`, interpolation helpers)

---

## High-level workflow

At a high level, most JURASSIC workflows follow this pattern:

1. **Read control parameters** into `ctl_t` (e.g. from a `.ctl` file).
2. **Read observation geometry** into `obs_t` (e.g. `obs.tab`).
3. **Read atmospheric state** into `atm_t` (e.g. `atm.tab`).
4. **Load lookup tables** into `tbl_t` (often driven by `ctl_t` settings).
5. Run one of:
   - **forward model** → compute radiances (`formod`)
   - **kernel/Jacobian** → sensitivity matrices (`kernel`)
   - **retrieval** → inverse problem (`optimal_estimation`)
6. **Write outputs** (tables, matrices, diagnostics) and optionally plot.

---

## Public API by category (developer map)

Below is a *conceptual* grouping of commonly used entry points from
`jurassic.h`. Exact details are documented in the header.

### Input / output

Typical project-style tools rely on table-based I/O helpers:

- Control / configuration:
  - `read_ctl`, `scan_ctl`
- Atmosphere and observation:
  - `read_atm(dirname, filename, ctl, atm, profile)`
  - `write_atm(dirname, filename, ctl, atm, profile)`
  - `read_obs(dirname, filename, ctl, obs, profile)`
  - `write_obs(dirname, filename, ctl, obs, profile)`
- Retrieval configuration:
  - `read_ret`
- Generic helpers:
  - `read_matrix(dirname, filename, ctl, matrix, dataset)`
  - `write_matrix(dirname, filename, ctl, matrix, atm, obs, rowspace, colspace, sort, dataset)`

For ASCII and binary formats, the `profile` / `dataset` arguments are ignored.
For netCDF formats, they select the record along the unlimited dimension.

> Tip: the small example tools (`formod.c`, `kernel.c`, `retrieval.c`)
> are good minimal references for correct I/O sequences.

### Lookup tables (spectroscopy / transmittance)

Lookup tables are handled via `tbl_t` and helpers such as:

- Read/write:
  - `read_tbl`, `read_tbl_asc`, `read_tbl_bin`, `read_tbl_nc_channel`
  - `write_tbl`, `write_tbl_asc`, `write_tbl_bin`, `write_tbl_nc`
- Interpolation utilities:
  - `intpol_tbl_ega`, `intpol_tbl_cga`, `intpol_tbl_eps`, `intpol_tbl_u`
  - `locate_tbl`

### Forward model (radiative transfer)

The main forward-model entry point is:

- `formod(const ctl_t *ctl, const tbl_t *tbl, atm_t *atm, obs_t *obs)`

Internally, `formod` can select different forward-model backends based on
the control settings (for example pencil-beam vs. line-by-line coupling).
The header comments document the selection logic.

### Operators and Jacobians

To compute Jacobians (kernels), the typical operator decomposition is:

- **state vector mapping**:
  - `atm2x` (atmosphere → state vector)
  - `x2atm` (state vector → atmosphere; if used)
- **measurement mapping**:
  - `obs2y` (observation container → measurement vector)
  - `y2obs` (measurement vector → observation container; if used)
- **kernel computation**:
  - `kernel(...)` (computes Jacobians / weighting functions)

### Retrieval / optimal estimation

The retrieval layer implements an optimal estimation workflow:

- `optimal_estimation(...)` — high-level retrieval driver
- minimization + diagnostics helpers (commonly used internally):
  - `cost_function(...)`
  - `analyze_avk(...)`, `analyze_avk_quantity(...)` (averaging kernels / diagnostics)
- shared-output helpers:
  - `shared_io_output_target(...)` (select legacy vs shared retrieval output target)
  - `shared_io_output_file(...)`, `shared_io_lock(...)`, `shared_io_unlock(...)`
    (coordinate shared netCDF retrieval output)

---

## Small command-line applications

In addition to the core library, JURASSIC includes small C programs that
use the library to perform common tasks. These are intentionally small,
and therefore serve as excellent “how do I call the library correctly?”
references.

### `formod.c` — forward model driver

`formod.c` demonstrates the minimal forward-modelling workflow:

- read `ctl_t` via `read_ctl`
- read `obs_t` via `read_obs`
- read `atm_t` via `read_atm`
- run `formod(...)`
- write/compare output (project scripts typically handle plotting/verification)

It also supports batching over a directory list (project-style runs).

### `kernel.c` — Jacobian / kernel driver

`kernel.c` demonstrates:

- reading control + inputs (`read_ctl`, `read_obs`, `read_atm`)
- mapping into vector form (`atm2x`, `obs2y`)
- calling `kernel(...)`
- writing matrices (e.g. via `write_matrix`)

### `retrieval.c` — optimal estimation driver

`retrieval.c` demonstrates:

- reading control parameters (`read_ctl`)
- reading retrieval parameters (`read_ret`)
- reading atmosphere + observation containers
- running `optimal_estimation(...)`

---

## Typical call sequences

### Forward model (radiances)

```c
ctl_t ctl;
atm_t atm;
obs_t obs;
tbl_t *tbl = ...; /* load tables as required */

read_ctl(argc, argv, &ctl);
read_obs(NULL, "obs.tab", &ctl, &obs, 0);
read_atm(NULL, "atm.tab", &ctl, &atm, 0);

formod(&ctl, tbl, &atm, &obs);

/* write/inspect obs.rad / diagnostics */
```

### Kernel (Jacobian) calculation

```c
read_ctl(argc, argv, &ctl);
read_obs(NULL, "obs.tab", &ctl, &obs, 0);
read_atm(NULL, "atm.tab", &ctl, &atm, 0);

atm2x(&ctl, &atm, x);
obs2y(&ctl, &obs, y, NULL, NULL);

kernel(&ctl, tbl, &atm, &obs, K);

/* write_matrix(NULL, "matrix.nc", &ctl, K, &atm, &obs, "y", "x", "r", 0) */
```

### Retrieval (optimal estimation)

```c
read_ctl(argc, argv, &ctl);
read_ret(argc, argv, &ctl, &ret);
read_obs(NULL, "obs_meas.tab", &ctl, &obs_meas, 0);
read_obs(NULL, "obs_apr.tab", &ctl, &obs_apr, 0);
read_atm(NULL, "atm_apr.tab", &ctl, &atm_apr, 0);
read_atm(NULL, "atm_ret.tab", &ctl, &atm_ret, 0);

optimal_estimation(&ret, &ctl, tbl,
                   &obs_meas, &obs_apr,
                   &atm_apr, &atm_ret, &chisq);

/* results are written into obs_apr/atm_ret and/or retrieval output files */
```

The exact argument lists and memory ownership rules are documented in
`jurassic.h`.

---

## Adding a new small tool

A typical new tool (e.g. `src/mytool.c`) should:

1. include `jurassic.h`
2. parse command-line arguments consistently with other tools
3. use the existing `read_*` helpers rather than duplicating parsers
4. keep the tool thin: put the science/numerics in the library, not in the CLI

Recommended starting point:

- copy the skeleton of `formod.c` (for forward runs) or `kernel.c`
  (for matrix-style outputs)
- replace only the “core call” (e.g. `formod` → your new library routine)

---

## Where to look next

If you are trying to understand a specific capability:

- **Forward modelling**: start at `formod(...)` in `jurassic.c`, then follow
  the called helpers for ray tracing / transmittance evaluation.
- **Lookup tables**: search for `read_tbl_*` and `intpol_tbl_*`.
- **Retrieval**: start at `optimal_estimation(...)` and follow the minimizer
  logic (Levenberg–Marquardt and diagnostics).

For API-level documentation, prefer the Doxygen comments in `jurassic.h`
as the primary reference.
