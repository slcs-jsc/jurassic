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

### Cloud and surface spectral grids

- `NCL` *(default: `0`)* and `CLNU[icl]` *(cm⁻¹)*
  Number and wavenumbers of cloud extinction grid points. When `NCL > 0`,
  the atmospheric input supplies one `clk[icl]` coefficient per point.
- `NSF` *(default: `0`)* and `SFNU[isf]` *(cm⁻¹)*
  Number and wavenumbers of surface emissivity grid points. When `NSF > 0`,
  the atmospheric input supplies surface temperature `sft` and one
  emissivity `sfeps[isf]` per point. Emissivity is linearly interpolated
  to detector channel wavenumbers. With `NSF = 0`, channel emissivity is 1
  and the surface temperature comes from the lower atmospheric boundary.
  With `NSF > 0`, a positive `sft` overrides that boundary temperature.

Each grid may be empty (`0` points) or contain at least two points;
`NCL = 1` and `NSF = 1` are rejected by the control parser.

### Surface treatment

- `SFTYPE` *(default: `2`)*
  `0`: no surface radiance; `1`: surface thermal emission;
  `2`: thermal emission plus simplified reflection of downwelling
  atmospheric thermal radiance; `3`: the preceding terms plus a
  simplified surface-reflected solar term. Reflection uses
  `1 - emissivity`. Surface terms apply only when the ray reaches the
  lowest altitude of the atmospheric profile. The solar option does not
  include atmospheric solar scattering.
- `SFSZA` *(default: `-999`, degrees)*
  Solar zenith angle at the surface for `SFTYPE = 3`. A negative value
  calculates the angle from observation time and surface location;
  a nonnegative value supplies the angle directly.

### Continua

Enable/disable built-in continuum contributions (interpreted as on/off switches).

- `CTM_CO2` *(default: `1`)*
- `CTM_H2O` *(default: `1`)*
  Uses the MT_CKD 4.1 water vapor continuum. H2O line lookup tables must be generated with RFM `H2O(sub)` so that the 25 cm⁻¹ line-wing contribution is not counted twice.
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

### Climatology

- `CLIM_FILE` *(default: `-`)*  
  Climatological data source for the `climatology` tool. `-` uses the
  built-in JURASSIC climatology. A JuMACS netCDF file provides temperature
  and mapped trace gases from JuMACS, while quantities not supplied by
  JuMACS retain built-in behavior.

!!! note
    JuMACS products are monthly zonal-mean climatologies, and the
    `climatology` tool interpolates them in climatological time, latitude,
    and altitude. Harmonized JuMACS products provide 12 monthly means on
    36 area-weighted 5-degree latitude bands; their upper altitude is
    model-dependent. Pressure always remains based on the built-in
    JURASSIC climatology, because the harmonized JuMACS products do not
    provide a pressure field.

    JuMACS values are used unchanged up to and including the model top of
    the product. Above the model top, temperature and most mapped gases
    smoothly transition back to the built-in climatology over a fixed
    20 km layer using a cosine taper: the JuMACS value is retained at the
    model top, the built-in climatology is recovered exactly at the end of
    the transition layer, and above the transition layer the built-in
    climatology is used unchanged. For selected gases whose built-in
    upper-atmosphere extrapolation is unsuitable, the top-level JuMACS
    value is instead kept constant with altitude above the model top:
    C2H6, CCl4, F11, F12, F22, H2CO, HNO3, HNO4, N2O, N2O5, NO2,
    and OCS. ClONO2 transitions to the built-in climatology
    logarithmically (geometrically) with altitude within the transition
    region.

    The mapped trace gases are C2H2, C2H6, CCl4, CH4, ClO, ClONO2, CO,
    F11, F12, F22, H2CO, H2O, H2O2, HCl, HNO3, HNO4, HOCl, N2O, N2O5, NO,
    NO2, O3, OCS, and SO2. Species not mapped to the product, including
    the main gases CO2, O2, and N2, continue to use the built-in
    climatology. Not every JuMACS product provides all mapped species;
    if a mapped species is missing from a specific product, JURASSIC
    keeps the built-in climatological value for that species. This is
    product-dependent and intentional.

    Above 120 km, the built-in temperature and trace-gas profiles are
    held constant at their 120 km values, while the built-in pressure
    continues to decrease with its existing logarithmic high-altitude
    extrapolation. Where a JuMACS product provides data above 120 km
    (e.g. WACCM-X up to 150 km), its values are used there according to
    the rules above.

    JURASSIC distributes three harmonized JuMACS climatology files in
    `data/`:

    | File | Model | Approximate vertical coverage |
    |------|-------|-------------------------------|
    | `jumacs_emac_refd1_1985-2014_5deg_1km.nc` | EMAC | 0–76 km |
    | `jumacs_geosccm_refd1_1985-2014_5deg_1km.nc` | GEOSCCM | 0–75 km |
    | `jumacs_waccmx_1985-2014_5deg_1km.nc` | WACCM-X | 0–150 km |

    The products differ in their covered species and upper altitude.

### Hydrostasy

- `HYDZ` *(default: `-999`)*  
  Optional hydrostatic reconstruction setting / control parameter.
  If unused in your workflow, keep the default.

### Field of view

- `FOV` *(default: `-`)*  
  Field-of-view configuration string (application-specific). If unused, keep `-`.

### Optional: retrieval and inversion configuration

These keywords select elements of the core optimal-estimation state
vector. The `ZMIN`/`ZMAX` pairs are inclusive altitude bounds in km:
each profile level within the pair enters the state vector. Both bounds
default to `-999`, which selects no ordinary profile levels. Cloud and
surface `RET_*` switches default to `0` (disabled); set one to `1` to
retrieve the corresponding quantity.

| Quantity | Atmospheric field | Control keyword(s) | State-vector elements |
| --- | --- | --- | --- |
| Pressure (hPa) | `p[ip]` | `RETP_ZMIN`, `RETP_ZMAX` | One per selected altitude level; pressure itself is retrieved |
| Temperature (K) | `t[ip]` | `RETT_ZMIN`, `RETT_ZMAX` | One per selected altitude level |
| Gas volume mixing ratio (ppv) | `q[ig][ip]` | `RETQ_ZMIN[ig]`, `RETQ_ZMAX[ig]` | One per selected level and emitter `ig` |
| Atmospheric extinction (km⁻¹) | `k[iw][ip]` | `RETK_ZMIN[iw]`, `RETK_ZMAX[iw]` | One per selected level and spectral window `iw` |
| Cloud center altitude (km) | `clz` | `RET_CLZ` | One per profile |
| Gaussian cloud vertical width (km) | `cldz` | `RET_CLDZ` | One per profile |
| Cloud extinction (km⁻¹) | `clk[icl]` | `RET_CLK` | One per `CLNU[icl]` point when `NCL > 0` |
| Surface temperature (K) | `sft` | `RET_SFT` | One per profile |
| Surface emissivity (dimensionless) | `sfeps[isf]` | `RET_SFEPS` | One per `SFNU[isf]` point when `NSF > 0` |

The cloud extinction is added to the atmospheric extinction with a
Gaussian altitude dependence centered at `clz` and scaled by `cldz`.
The surface emissivity elements are interpolated to detector channels.
Instrument parameters, including spectral shift, gain, pointing, field
of view, channel response, and offset, are not part of the current core
retrieval state vector.

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
RFM ASCII spectra read by JURASSIC must contain at least two spectral points.

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

- Start from the **example projects** (`projects/examples/limb`, `projects/examples/nadir`, `projects/examples/zenith`) and modify incrementally.
- Keep `ND`, `NW`, and `WINDOW[i]` consistent; mismatches are a common source of confusing results.
- When changing `NG`/`EMITTER`, ensure the corresponding lookup tables exist under `TBLBASE`.
- If you tune `RAYDS`/`RAYDZ`, validate against reference output to ensure accuracy remains acceptable.
