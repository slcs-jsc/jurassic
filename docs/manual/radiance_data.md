# Radiance and observation data

This page documents the tabular data files used by JURASSIC to describe
**observation geometry** and to store **radiances** (measured or
simulated). In typical workflows you will encounter two closely related
files:

- `obs.tab` — **observation definition** (geometry input for radiative transfer)
- `rad.tab` — **radiance data** (geometry + computed tangent point + radiances)

Both files share the same leading geometry columns so they can be
compared, merged, or processed with the same downstream tools.

> Recommended convention: treat `obs.tab` as *input* (requested geometry)
> and `rad.tab` as *output* (evaluated geometry + radiances).  
> Some tools/workflows may overwrite radiance columns in `obs.tab`, but
> using a distinct `rad.tab` is usually clearer.

---

## Role in typical workflows

### Forward simulation

1. Prepare `obs.tab` with the requested observation geometry (radiances
   may be dummy values).
2. Run the forward model (`formod`).
3. JURASSIC writes `rad.tab` containing:
   - the same geometry block (copied through),
   - the **computed** tangent point (after refraction-aware ray tracing),
   - **simulated radiances** (per detector channel),
   - optionally transmittances and other per-channel quantities.

### Retrieval (inverse modelling)

1. Provide `rad.tab` where the radiance columns contain **measured**
   radiances (and geometry is known).
2. The retrieval tool iteratively runs the forward model internally.
3. Outputs may include:
   - fitted radiances,
   - residuals,
   - retrieval diagnostics (Jacobians, averaging kernels, covariances).

---

## File format selection

The observation/radiance table format is selected via the control
parameter:

- `OBSFMT`

The description below refers to the **default ASCII format**, typically
used with:

```text
OBSFMT = 1
```

This is the format used by the example projects distributed with
JURASSIC.

---

## General format rules

- Plain-text ASCII table.
- Columns separated by spaces or tabs.
- Header lines are optional; if present, they must not start with a
  valid numeric value in the first column.
- Each row corresponds to one observation (one line of sight / one ray).

---

## Geometry block (shared by `obs.tab` and `rad.tab`)

In the default ASCII format, the **first 10 columns** define the
observation geometry:

1. `time` — time stamp (s)
2. `obsz` — observer altitude (km)
3. `obslon` — observer longitude (deg)
4. `obslat` — observer latitude (deg)
5. `vpz` — view point altitude (km)
6. `vplon` — view point longitude (deg)
7. `vplat` — view point latitude (deg)
8. `tpz` — tangent point altitude (km)
9. `tplon` — tangent point longitude (deg)
10. `tplat` — tangent point latitude (deg)

### Interpretation

- The observer position (`obs*`) is typically the instrument location
  (satellite altitude for limb/nadir, ground-based for zenith, etc.).
- The view point (`vp*`) is a point used to define the line of sight.
  Depending on geometry, this may be a ground intercept, a target point,
  or another convenient reference.
- The tangent point (`tp*`) in `obs.tab` is typically an **initial** or
  **nominal** value (or can be set to a placeholder).  
  In `rad.tab` it is typically replaced/filled with the **computed**
  tangent point derived from the actual refracted ray path.

> Note: The exact interpretation of `vp*` and the use of tangent point
> fields depends on geometry (limb, nadir, occultation). The key point
> is that `rad.tab` should contain the *evaluated* geometry consistent
> with the ray tracing used by the forward model.

---

## Radiance and transmittance blocks

After the 10 geometry columns, the table contains per-channel data.

The number of detector channels is defined in the control file by:

- `ND`
- `NU[i]` (channel wavenumbers)

### Radiances

- `rad[i]` — radiance for detector channel *i*

There are **ND** radiance columns:

```text
rad[0] ... rad[ND-1]
```

#### Meaning in `obs.tab`
Radiance values in `obs.tab` are often **dummy** (e.g. 0.0) when the file
is used purely to define geometry for forward simulations.

#### Meaning in `rad.tab`
Radiance values in `rad.tab` are typically **physical radiances**:
- simulated radiances (forward model output), or
- measured radiances (retrieval input).

### Transmittances

Some workflows also store band transmittances:

- `tau[i]` — total atmospheric transmittance for channel *i*

There are **ND** transmittance columns:

```text
tau[0] ... tau[ND-1]
```

Transmittances are dimensionless (typically 0–1) and are useful for:

- diagnostics and validation,
- sanity checks of spectral windows,
- sensitivity studies.

### Total column count

For `OBSFMT = 1`, the default layout is:

```text
10 + ND + ND  =  10 + 2·ND
```

(i.e., geometry + radiances + transmittances).

---

## Units and conventions

JURASSIC does not enforce a single radiance unit. Radiance units must be
**consistent across your workflow**, especially between:

- measured radiances (if used),
- lookup tables / Planck normalization,
- any post-processing scripts.

Typical radiance conventions include:
- spectral radiance (e.g. W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹),
- band-integrated radiance,
- brightness-temperature-equivalent radiance (after conversion).

Geometry units are as listed above (km, deg, seconds).

---

## Practical examples (schematic)

### Example `obs.tab` row (geometry defined, radiances dummy)

For `ND = 3`:

```text
time  obsz  obslon  obslat  vpz  vplon  vplat  tpz  tplon  tplat  rad0 rad1 rad2  tau0 tau1 tau2
0.0   700.0 0.0     50.0    0.0  0.0    50.0   20.0 0.0   50.0   0.0  0.0  0.0   1.0  1.0  1.0
```

### Example `rad.tab` row (computed tangent point + radiances)

The first 7 geometry columns may remain identical, while the tangent
point and radiances are filled with evaluated results:

```text
time  obsz  obslon  obslat  vpz  vplon  vplat  tpz    tplon  tplat  rad0    rad1    rad2    tau0  tau1  tau2
0.0   700.0 0.0     50.0    0.0  0.0    50.0   19.83  0.02   49.98  1.23e-7 1.18e-7 9.95e-8 0.12  0.15  0.22
```

(Values shown are illustrative.)

---

## Consistency requirements and common pitfalls

- Ensure `ND` in your control file matches the number of `rad[]` and
  `tau[]` columns in your tables.
- Keep the geometry columns consistent between `obs.tab` and `rad.tab`
  if you rely on line-by-line comparisons.
- If you import measured radiances, verify unit consistency and channel
  ordering (`NU[i]`).
- When using refraction (`REFRAC = 1`), prefer using the tangent point
  stored in `rad.tab` for scientific interpretation and plotting.

---

## Summary

`obs.tab` and `rad.tab` encode the same observation list:

- `obs.tab` defines the *requested* observation geometry used as input
  for radiative transfer calculations.
- `rad.tab` contains the *evaluated* observation geometry (including the
  computed tangent point) and the corresponding radiances (simulated or
  measured), with optional transmittances.

This separation provides a clear and reproducible workflow for forward
simulations and retrieval applications.

---

## Related pages

- [Configuration](configuration.md)
- [Running JURASSIC](running.md)
- [Atmospheric data](atmospheric_data.md)
- [Diagnostics & logs](diagnostics.md)
