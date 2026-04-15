# Configuration

JURASSIC is configured through a **control file** (often called `*.ctl`)
and optional **command-line overrides**. Most applications in the
repository, such as the forward model, kernel, and retrieval tools, use
the same control-file syntax and keyword names.

This page documents the configuration **syntax** and the most important **control keywords**.

---

## Control file syntax

### Basic format

JURASSIC reads a control file from the **first command-line argument** if it does **not** start with `-`.
Lines starting with `#` are treated as comments and ignored. Each other
non-empty line is interpreted as:

```text
KEY = VALUE
```

The separator is required but otherwise only serves to separate the key
from the value.

### Arrays and wildcards

Many parameters are arrays (e.g. `NU[0]`, `EMITTER[1]`). Array indices
are zero-based, so the first entry is index `0`. You can use:

- `KEY[i] = VALUE` for a specific index `i`
- `KEY[*] = VALUE` as a **default for all indices**

When reading a control file, JURASSIC uses the first matching entry it
finds. If you want to combine a wildcard default with index-specific
values, put the index-specific entries before the wildcard entry.
Command-line overrides are evaluated after the control file and can
therefore override both exact and wildcard entries.

Examples:

```text
WINDOW[3] = 1
WINDOW[*] = 0

RETQ_ZMAX[1] = 60
RETQ_ZMIN[*] = 0
RETQ_ZMAX[*] = 80
```

### Command-line overrides

All control keywords can also be provided on the command line as pairs:

```bash
./formod run.ctl obs.tab atm.tab rad.tab RAYDS 5 REFRAC 0
```

Command-line values override values from the control file.

---

## Key concepts

Most JURASSIC applications follow the same high-level flow:

1. Read control parameters
2. Read observation geometry (from application inputs interpreted using `OBSFMT`)
3. Read atmospheric profiles (from application inputs interpreted using `ATMFMT`)
4. Configure spectral definition (channels + spectral windows + lookup tables)
5. Run the forward model (and optionally kernels / retrieval)
6. Write output products

---

## Control keywords

Below is a practical overview of commonly used control keywords and
their default values.
“Array” means that the keyword is indexed, e.g. `NU[0]`.

### Emitters and gases

These control which gases/emitters are considered in radiative transfer.

- `NG` *(default: `0`)*  
  Number of emitters.

- `EMITTER[i]` *(array, default: required if `NG>0`)*  
  Name of emitter *i* (e.g. `CO2`, `H2O`, ...).

!!! tip
    You typically choose emitters that have lookup tables available under
    your `TBLBASE`.

### Spectral channels and windows

Channels are detector wavenumbers (or band centers) and each channel belongs to a spectral window.

- `ND` *(default: `0`)*  
  Number of detector channels.

- `NU[i]` *(array, default: required if `ND>0`)*  
  Wavenumber of channel *i*.

- `NW` *(default: `1`)*  
  Number of spectral windows.

- `WINDOW[i]` *(array, default: `0`)*  
  Assign window index for channel *i*. (Typical use: group channels into bands for lookup tables and continua.)

### Continua

Enable/disable built-in continuum contributions (interpreted as on/off switches).

- `CTM_CO2` *(default: `1`)*
- `CTM_H2O` *(default: `1`)*
- `CTM_N2`  *(default: `1`)*
- `CTM_O2`  *(default: `1`)*

### Refraction and ray tracing

- `REFRAC` *(default: `1`)*  
  Toggle atmospheric refraction in ray tracing (`1` = on, `0` = off).
  When enabled, JURASSIC uses a dry-air Edlen-type refractivity
  approximation based on pressure and temperature, evaluated at a fixed
  reference wavelength of 10 um.

- `RAYDS` *(default: `10`)*  
  Maximum step length along the ray path (typically in km; see the theory/limitations pages for guidance).

- `RAYDZ` *(default: `0.1`)*  
  Additional vertical step-size constraint (used to keep vertical sampling sufficiently fine, especially for nadir).

### File formats and lookup tables

- `TBLBASE` *(default: `-`)*  
  Base path/prefix for emissivity lookup tables.

- `TBLFMT` *(default: `1`)*  
  Lookup table format selector:
  `1` = ASCII, `2` = binary, `3` = netCDF.

- `ATMFMT` *(default: `1`)*  
  Atmosphere file format selector.

- `OBSFMT` *(default: `1`)*  
  Observation geometry file format selector.

!!! note
    The format selectors must match the files you provide. For example,
    `ATMFMT` must match the atmosphere file, `OBSFMT` must match the
    observation file, and `TBLFMT` must match the lookup-table files
    under `TBLBASE`.

    With `TBLFMT = 3`, lookup tables are read from netCDF files. In this
    layout, JURASSIC expects one netCDF file per emitter and one packed
    lookup-table variable per channel inside that file.

    When the same files are used by several tools, use the same format
    selectors in each control file or command-line override.

### Hydrostasy

- `HYDZ` *(default: `-999`)*  
  Optional hydrostatic reconstruction setting / control parameter.
  If unused in your workflow, keep the default.

### Field of view

- `FOV` *(default: `-`)*  
  Field-of-view configuration string (application-specific). If unused, keep `-`.

### Optional: retrieval and inversion configuration

These keywords are used by applications that perform optimal-estimation retrievals (or compute retrieval-related products).

**Vertical limits (z-min/z-max):**

- `RETP_ZMIN`, `RETP_ZMAX` *(defaults: `-999`)*  
  Pressure-related retrieval limits (application meaning depends on the retrieval setup).

- `RETT_ZMIN`, `RETT_ZMAX` *(defaults: `-999`)*  
  Temperature retrieval limits.

- `RETQ_ZMIN[i]`, `RETQ_ZMAX[i]` *(array, defaults: `-999`)*  
  VMR retrieval limits for emitter *i*.

- `RETK_ZMIN[iw]`, `RETK_ZMAX[iw]` *(array, defaults: `-999`)*  
  Retrieval limits for aerosol/extinction or window-dependent quantities (indexed by spectral window).

**Cloud / surface fit toggles:**

- `RET_CLZ` *(default: `0`)*
- `RET_CLDZ` *(default: `0`)*
- `RET_CLK` *(default: `0`)*
- `RET_SFT` *(default: `0`)*
- `RET_SFEPS` *(default: `0`)*

(These switches enable/disable optional retrieval degrees of freedom. Exact interpretation depends on the chosen
retrieval application and state-vector definition.)

**Shared retrieval file mode:**

- `SHARED_IO_PROFLIST` *(default: `-`)*  
  Optional text file containing one retrieval profile index per line.

- `SHARED_IO_ATM_APR_FILE`, `SHARED_IO_OBS_MEAS_FILE` *(default: `-`)*  
  Optional shared retrieval input files. When set together with
  `SHARED_IO_PROFLIST`, the retrieval application reads atmospheric and
  observation records from these files instead of from per-directory
  `atm_apr.*` / `obs_meas.*` files.

- `SHARED_IO_ATM_FINAL_FILE`, `SHARED_IO_OBS_FINAL_FILE` *(default: `-`)*  
  Optional shared retrieval output files.

- `SHARED_IO_MATRIX_COV_APR_FILE`, `SHARED_IO_MATRIX_KERNEL_FILE`,
  `SHARED_IO_MATRIX_COV_RET_FILE`, `SHARED_IO_MATRIX_CORR_FILE`,
  `SHARED_IO_MATRIX_GAIN_FILE`, `SHARED_IO_MATRIX_AVK_FILE`
  *(default: `-`)*  
  Optional shared matrix output files.

- `SHARED_IO_ATM_ERR_TOTAL_FILE`, `SHARED_IO_ATM_ERR_NOISE_FILE`,
  `SHARED_IO_ATM_ERR_FORMOD_FILE`, `SHARED_IO_ATM_CONT_FILE`,
  `SHARED_IO_ATM_RES_FILE` *(default: `-`)*  
  Optional shared atmospheric diagnostic output files.

### Output control

- `WRITE_BBT` *(default: `0`)*  
  Write brightness temperatures instead of radiances, where applicable.

- `WRITE_MATRIX` *(default: `0`)*  
  Write retrieval matrices (e.g. Jacobian/averaging kernel/error covariance), where supported.

### Forward model backend selection

- `FORMOD` *(default: `1`)*  
  Select the forward model implementation / backend (application-specific).
  Keep the default unless you explicitly built and intend to use an alternative backend.

### Line-by-line / Reference Forward Model (RFM) integration hooks

Some workflows use external line-by-line assets or tools, for example
the [Reference Forward Model (RFM)](https://eodg.atm.ox.ac.uk/RFM/) for
table generation. These parameters provide paths.

- `RFMBIN` *(default: `-`)*  
  Path to the RFM binary (if used).

- `RFMHIT` *(default: `-`)*  
  Path to the HITRAN file used by RFM (if used).

- `RFMXSC[i]` *(array, default: `-`)*  
  Per-emitter cross-section file (if used).

---

## Minimal examples

### Forward simulation (single window)

```text
# minimal_forward.ctl

# Lookup-table prefix.
TBLBASE = ../tests/data/airs

# Active emitters.
NG = 1
EMITTER[0] = CO2

# Spectral channels.
ND = 3
NU[0] = 667.7820
NU[1] = 668.5410
NU[2] = 669.8110

# Write brightness temperatures instead of radiances.
WRITE_BBT = 1
```

Run from `src/` with compatible `obs.tab` and `atm.tab` files:

```bash
./formod minimal_forward.ctl obs.tab atm.tab rad.tab
```

For a complete runnable setup, start from the example projects described
in the [Quickstart](quickstart.md).

### Overriding parameters on the command line

```bash
./formod minimal_forward.ctl obs.tab atm.tab rad.tab RAYDS 5 REFRAC 0
```

---

## Tips for robust configurations

- Start from the **example projects** (`projects/limb`, `projects/nadir`, `projects/zenith`) and modify incrementally.
- Keep `ND`, `NW`, and `WINDOW[i]` consistent; mismatches are a common source of confusing results.
- When changing `NG`/`EMITTER`, ensure the corresponding lookup tables exist under `TBLBASE`.
- If you tune `RAYDS`/`RAYDZ`, validate against reference output to ensure accuracy remains acceptable.
