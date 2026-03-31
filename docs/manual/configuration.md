# Configuration

JURASSIC is configured through a **control file** (often called `*.ctl`) and optional **command-line overrides**.
Most applications in the repository (e.g. forward model, kernel, retrieval tools) share the same control-parameter
infrastructure and populate a `ctl_t` structure internally.

This page documents the configuration **syntax** and the most important **control keywords**.

---

## Control file syntax

### Basic format

JURASSIC reads a control file from the **first command-line argument** if it does **not** start with `-`.
Each non-empty line is interpreted as:

```text
KEY = VALUE
```

Internally, this is parsed by reading three tokens (`KEY`, a dummy token such as `=`, and `VALUE`).

### Arrays and wildcards

Many parameters are arrays (e.g. `NU[0]`, `EMITTER[1]`). The parser supports:

- `KEY[i] = VALUE` for a specific index `i`
- `KEY[*] = VALUE` as a **default for all indices** (can be overridden by specific indices)

Examples:

```text
NU[*] = 1000.0
NU[0] = 925.0

EMITTER[*] = CO2
EMITTER[1] = H2O
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

1. Read control parameters (`read_ctl()` → `ctl_t`)
2. Read observation geometry (from `OBS...` / `OBSFMT`)
3. Read atmospheric profiles (from `ATM...` / `ATMFMT`)
4. Configure spectral definition (channels + spectral windows + lookup tables)
5. Run the forward model (and optionally kernels / retrieval)
6. Write output products

---

## Control keywords

Below is a practical overview of the control keywords evaluated by `read_ctl()` and their default values.
“Array” means that the keyword is indexed, e.g. `NU[0]`.

### Emitters and gases

These control which gases/emitters are considered in radiative transfer.

- `NG` *(default: `0`)*  
  Number of emitters.

- `EMITTER[i]` *(array, default: required if `NG>0`)*  
  Name of emitter *i* (e.g. `CO2`, `H2O`, ...).

> Tip: You typically choose emitters that have lookup tables available under your `TBLBASE`.

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

- `RAYDS` *(default: `10`)*  
  Maximum step length along the ray path (typically in km; see the theory/limitations pages for guidance).

- `RAYDZ` *(default: `0.1`)*  
  Additional vertical step-size constraint (used to keep vertical sampling sufficiently fine, especially for nadir).

### File formats and lookup tables

- `TBLBASE` *(default: `-`)*  
  Base path/prefix for emissivity lookup tables.

- `TBLFMT` *(default: `1`)*  
  Lookup table format selector.

- `ATMFMT` *(default: `1`)*  
  Atmosphere file format selector.

- `OBSFMT` *(default: `1`)*  
  Observation geometry file format selector.

> Notes:
> - In the example projects, the formats and table layout are chosen to match the provided `atm.tab`, `obs.tab`,
>   `rad.tab`, and associated lookup tables.
> - If you introduce new formats, keep the “format selectors” consistent across applications.

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

### Output control

- `WRITE_BBT` *(default: `0`)*  
  Write brightness temperature output (if applicable).

- `WRITE_MATRIX` *(default: `0`)*  
  Write retrieval matrices (e.g. Jacobian/averaging kernel/error covariance), where supported.

### Forward model backend selection

- `FORMOD` *(default: `1`)*  
  Select the forward model implementation / backend (application-specific).
  Keep the default unless you explicitly built and intend to use an alternative backend.

### Line-by-line / RFM integration hooks

Some workflows use external line-by-line assets or tools (e.g. for table generation). These parameters provide paths.

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
NG = 2
EMITTER[0] = CO2
EMITTER[1] = H2O

ND = 3
NU[0] = 930.0
NU[1] = 950.0
NU[2] = 970.0
WINDOW[*] = 0
NW = 1

TBLBASE = tables/midlat_example
TBLFMT  = 1

ATMFMT = 1
OBSFMT = 1

REFRAC = 1
RAYDS  = 10
RAYDZ  = 0.1

CTM_CO2 = 1
CTM_H2O = 1
CTM_N2  = 1
CTM_O2  = 1
```

Run (example):

```bash
./formod minimal_forward.ctl obs.tab atm.tab rad.tab
```

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
