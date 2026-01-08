# Input profiles and geometry files

JURASSIC simulations typically use two primary tabular input files:

- an **atmospheric profile file** (commonly `atm.tab`)
- an **observation geometry file** (commonly `obs.tab`)

Both files are plain-text tables in the default configuration and are
read by the JURASSIC I/O routines (`read_atm()` and `read_obs()`), which
populate the internal `atm_t` and `obs_t` structures.

This page documents the **default ASCII formats** used by the example
projects and by many user workflows.

> The exact file formats are selected via the control parameters
> `ATMFMT` and `OBSFMT`. The formats documented here correspond to the
> default ASCII readers (`read_atm_asc()` and `read_obs_asc()`), typically
> used with `ATMFMT = 1` and `OBSFMT = 1`.

---

## General format rules

### Columns and delimiters

- Columns are separated by **spaces** or **tabs**.
- There is **no header requirement**. If you include a header line, it
  must not start with a valid numeric value in the first column, so it
  can be skipped by the reader.
- Extra columns are not supported; provide exactly the required number
  of fields per row for the chosen configuration.

### Optional columns depend on your configuration

The number of required columns depends on your control parameters:

- `NG` / `EMITTER[i]` controls how many gas VMR columns are read.
- `NW` controls how many extinction columns are read (one per spectral window).
- Retrieval and surface/cloud options can add additional fields (see below).

---

## Atmospheric profile file (`atm.tab`)

### Purpose

The atmospheric profile file defines the thermodynamic and compositional
state of the atmosphere along the simulated ray paths. In the default
ASCII format, each row describes one vertical grid point.

### Units (conventions used in the code)

The file reader expects the following unit conventions:

- time: **seconds**
- altitude `z`: **km**
- longitude/latitude: **degrees**
- pressure `p`: **hPa**
- temperature `t`: **K**
- volume mixing ratios `q`: **ppv** (typically dimensionless, e.g. mol/mol)
- extinction `k`: **km⁻¹**

> Note: These conventions are consistent with the logging messages in the
> JURASSIC library and the example projects.

### Column layout (default ASCII format)

For each profile level (row), the expected columns are:

1. `time`
2. `z`
3. `lon`
4. `lat`
5. `p`
6. `t`
7. `q[0] ... q[NG-1]`  (one column per gas/emitter)
8. `k[0] ... k[NW-1]`  (one column per spectral window)

In other words, the number of columns is:

- **6 + NG + NW** (plus optional cloud/surface fields as described below)

### Optional: cloud layer parameters (single set)

If cloud-layer parameters are enabled (`NCL > 0` in the control), then
the **first profile row only** (`np == 0`) is expected to include
additional columns:

- `clz`   (cloud layer altitude, km)
- `cldz`  (cloud layer thickness, km)
- `clk[0] ... clk[NCL-1]`  (cloud extinction per cloud channel/bin, km⁻¹)

These values describe a simplified cloud representation that can be used
by certain workflows.

### Optional: surface parameters (single set)

If surface parameters are enabled (`NSF > 0` in the control), then
the **first profile row only** (`np == 0`) is expected to include:

- `sft` (surface temperature, K)
- `sfeps[0] ... sfeps[NSF-1]` (surface emissivity parameters)

### Example (schematic)

For `NG=2` gases and `NW=1` window, a single profile row would look like:

```text
time   z    lon   lat    p      t      q_CO2      q_H2O      k_win0
0.0    0.0  0.0   50.0   1013.0  288.0  4.2e-4    7.0e-3     0.0
```

---

## Observation geometry file (`obs.tab`)

### Purpose

The observation geometry file defines the viewing geometry for each
simulated ray (line of sight). In the default ASCII format, each row
describes one ray and can also contain measured radiances and
transmittances (used, for example, in retrieval workflows).

### Units (conventions used in the code)

- time: **seconds**
- altitudes: **km**
- longitude/latitude: **degrees**
- radiance values: model-specific units (must be consistent across the workflow)
- transmittance `tau`: **dimensionless** (typically 0–1)

### Column layout (default ASCII format)

For each ray (row), the expected columns are:

1. `time`
2. `obsz`    (observer altitude, km)
3. `obslon`  (observer longitude, deg)
4. `obslat`  (observer latitude, deg)
5. `vpz`     (view point altitude, km)
6. `vplon`   (view point longitude, deg)
7. `vplat`   (view point latitude, deg)
8. `tpz`     (tangent point altitude, km)
9. `tplon`   (tangent point longitude, deg)
10. `tplat`  (tangent point latitude, deg)
11. `rad[0] ... rad[ND-1]` (radiances per detector channel)
12. `tau[0] ... tau[ND-1]` (transmittances per detector channel)

So the number of columns is:

- **10 + 2·ND**

### Notes on `rad[]` and `tau[]`

Depending on the application:

- In **retrieval mode**, `rad[]` often represents the measured radiances.
- In **forward model mode**, these fields may be initialized by the
  project scripts or can be overwritten by the model outputs.

The example projects generate an `obs.tab` file that is consistent with
the selected detector configuration (`ND` and `NU[i]`) in the control file.

### Example (schematic)

For `ND=3` channels, a single ray row would look like:

```text
time  obsz  obslon  obslat  vpz  vplon  vplat  tpz  tplon  tplat  rad0 rad1 rad2  tau0 tau1 tau2
0.0   700.0 0.0     50.0    0.0  0.0    50.0   20.0 0.0   50.0   0.0  0.0  0.0   1.0  1.0  1.0
```

---

## Practical tips

- Keep `NG`, `NW`, and `ND` in your control file consistent with the
  number of columns in `atm.tab` and `obs.tab`.
- Start from the provided **example projects** and adjust incrementally.
- If you change the list of emitters (`EMITTER[i]`), update the gas
  VMR columns in `atm.tab` accordingly.
- If you use retrieval tools, ensure that `obs.tab` contains the
  radiance values you intend to fit (measured or synthetic).

---

## Related pages

- [Configuration](configuration.md) — control file keywords and overrides
- [Theory](theory.md) — radiative transfer background
- [Retrieval theory](retrieval_theory.md) — optimal estimation concepts
