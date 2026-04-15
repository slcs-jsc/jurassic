# Atmospheric data

This page describes the atmospheric input data used by JURASSIC. In the
example projects this is typically an ASCII file named `atm.tab`, but
JURASSIC can also read and write atmospheric profiles in netCDF format.
The atmospheric data file defines the thermodynamic and compositional
state of the atmosphere used for radiative transfer and retrieval
calculations.

In JURASSIC, the atmospheric data represent the **model state**, while
observation geometry and radiances are handled separately (see
[Observation and radiance files](radiance_data.md)).

---

## Purpose of the atmospheric data file

The atmospheric data file provides vertical profiles of:

- pressure
- temperature
- trace gas volume mixing ratios
- optional aerosol or cloud extinction parameters

These profiles are assumed to be horizontally homogeneous along each
ray path. They are interpolated internally to the ray-tracing grid
during radiative transfer calculations.

---

## File format selection

The format of the atmospheric data file is selected via the control
parameter:

- `ATMFMT`

The **default ASCII format** is selected with:

```text
ATMFMT = 1
```

This is the format used in the example projects distributed with
JURASSIC and is described first below.

JURASSIC also supports:

- `ATMFMT = 2` for the compact binary format used by JURASSIC,
- `ATMFMT = 3` for netCDF files.

The netCDF format is useful for workflows that manage many atmospheric
profiles in one file, for example shared-input or retrieval workflows
that select profiles by index.

!!! note
    `ATMFMT = 2` selects the JURASSIC C binary format. It was implemented
    for workflows where file-I/O performance is critical, because reading
    and writing ASCII tables can be much slower than binary I/O. For new
    workflows, users are generally encouraged to use the netCDF format
    instead.

---

## ASCII format rules

- The file is a plain-text ASCII table.
- Columns are separated by spaces or tabs.
- Header lines are optional; if present, they must not start with a
  valid numeric value in the first column.
- Each row corresponds to one vertical grid point.

---

## Units and conventions

The atmospheric readers expect the following units:

| Quantity            | Unit |
|---------------------|------|
| time                | s    |
| altitude (`z`)      | km   |
| longitude, latitude | deg  |
| pressure (`p`)      | hPa  |
| temperature (`t`)   | K    |
| volume mixing ratio | ppv |
| extinction (`k`)    | km⁻¹ |

These conventions are consistent with the internal JURASSIC data
structures and example workflows.

---

## Column layout (default ASCII format)

For each vertical level, the atmospheric file contains the following
columns:

1. `time`  
2. `z` – altitude  
3. `lon` – longitude  
4. `lat` – latitude  
5. `p` – pressure  
6. `t` – temperature  
7. `q[0] ... q[NG-1]` – volume mixing ratios for each emitter  
8. `k[0] ... k[NW-1]` – extinction coefficients for each spectral window  

The total number of columns is therefore:

```text
6 + NG + NW
```

where:

- `NG` is the number of emitters defined via `EMITTER[i]`,
- `NW` is the number of spectral windows.

---

## Trace gas volume mixing ratios

The number and order of trace gas columns must match the emitters
defined in the control file:

```text
NG = 2
EMITTER[0] = CO2
EMITTER[1] = H2O
```

In this case, the atmospheric file must contain exactly two VMR columns,
in the same order.

---

## Aerosol and cloud extinction

Extinction coefficients are provided per spectral window. If aerosol or
cloud extinction is not used, the corresponding columns should still be
present and set to zero.

Extinction values are assumed to be vertically distributed and are
applied along the ray path during radiative transfer calculations.

---

## Optional: cloud layer parameters

Some workflows support simplified cloud-layer representations defined
by a small set of parameters. If enabled via the control file, the
ASCII format used by JURASSIC appends additional columns such as:

- cloud-layer altitude
- cloud-layer thickness
- cloud extinction coefficients

These cloud-layer values are written on every row in the standard ASCII
output, even though they represent profile-level parameters. When
reading ASCII input, JURASSIC uses the values from the first row.

If these options are not used, the corresponding parameters should be
disabled in the control file and omitted from the atmospheric data.

---

## Optional: surface parameters

Surface-related parameters (e.g. surface temperature or emissivity) may
also be specified in the atmospheric data file for certain applications.

As with the cloud-layer parameters, the standard ASCII writer appends
these values on every row, but the ASCII reader consumes them from the
first row only.

If surface effects are not required, the default configuration should
be used and these fields omitted.

---

## netCDF format

For `ATMFMT = 3`, JURASSIC expects a netCDF file with a profile
dimension and a level dimension. A single file may contain several
profiles; applications select the desired zero-based profile index when
reading the file.

!!! note
    A major advantage of the netCDF format is that many atmospheric
    profiles can be stored in a single file. This avoids many of the
    file-system problems that can occur when large runs are organized as
    thousands of directories containing small ASCII files.

The standard JURASSIC writer creates the dimensions:

| Dimension | Meaning |
|-----------|---------|
| `profile` | atmospheric profile index |
| `level` | vertical grid level |

The variable `nlev` gives the number of valid vertical levels for each
profile. Variables that depend on height use the dimensions
`profile, level`; profile-level cloud and surface parameters use only
the `profile` dimension.

The required core variables are:

| Variable | Unit | Meaning |
|----------|------|---------|
| `nlev` | 1 | number of vertical levels |
| `time` | s | time in seconds since 2000-01-01, 00:00 UTC |
| `z` | km | altitude |
| `lon` | degrees_east | longitude |
| `lat` | degrees_north | latitude |
| `p` | hPa | pressure |
| `t` | K | temperature |

Trace-gas variables are named exactly as the emitters in the control
file. For example:

```text
NG = 2
EMITTER[0] = CO2
EMITTER[1] = H2O
```

requires netCDF variables named `CO2` and `H2O`, both in units of
volume mixing ratio (`ppv`). The order is still defined by the control
file, but the netCDF reader finds the data by variable name.

Extinction variables are named by spectral window:

```text
ext_win_0
ext_win_1
...
ext_win_<NW-1>
```

The number of variables must match `NW`.

If cloud-layer parameters are enabled with `NCL > 0`, the netCDF file
must also contain:

- `cld_z` for cloud-layer height in km,
- `cld_dz` for cloud-layer depth in km,
- `cld_k_<wavenumber>` for cloud-layer extinction at each `CLNU[i]`.

For example, `CLNU[0] = 800.0000` corresponds to a variable named
`cld_k_800.0000`.

If surface parameters are enabled with `NSF > 0`, the netCDF file must
also contain:

- `srf_t` for surface temperature in K,
- `srf_eps_<wavenumber>` for surface emissivity at each `SFNU[i]`.

For example, `SFNU[0] = 925.0000` corresponds to a variable named
`srf_eps_925.0000`.

!!! note
    The control file and the netCDF file must describe the same
    quantities. In particular, `NG`, `EMITTER[i]`, `NW`, `NCL`, `CLNU[i]`,
    `NSF`, and `SFNU[i]` determine which netCDF variables JURASSIC will
    look for.

A compact `ncdump -h`-style skeleton for `NG = 2`, `NW = 1`, `NCL = 1`,
and `NSF = 1` could look like:

```text
netcdf atm {
dimensions:
  profile = UNLIMITED ;
  level = 121 ;

variables:
  int nlev(profile) ;
    nlev:long_name = "number of vertical levels" ;
    nlev:units = "1" ;

  double time(profile, level) ;
    time:units = "s" ;
  double z(profile, level) ;
    z:units = "km" ;
  double lon(profile, level) ;
    lon:units = "degrees_east" ;
  double lat(profile, level) ;
    lat:units = "degrees_north" ;
  double p(profile, level) ;
    p:units = "hPa" ;
  double t(profile, level) ;
    t:units = "K" ;

  double CO2(profile, level) ;
    CO2:units = "ppv" ;
  double H2O(profile, level) ;
    H2O:units = "ppv" ;

  double ext_win_0(profile, level) ;
    ext_win_0:units = "km**-1" ;

  double cld_z(profile) ;
    cld_z:units = "km" ;
  double cld_dz(profile) ;
    cld_dz:units = "km" ;
  double cld_k_800.0000(profile) ;
    cld_k_800.0000:units = "km**-1" ;

  double srf_t(profile) ;
    srf_t:units = "K" ;
  double srf_eps_925.0000(profile) ;
    srf_eps_925.0000:units = "1" ;
}
```

Cloud and surface variables are only required when the corresponding
options are enabled in the control file.

---

## Example (schematic)

For `NG = 2` emitters and `NW = 1` spectral window, a single profile row
might look like:

```text
time   z    lon   lat    p      t      q_CO2      q_H2O      k_win0
0.0    0.0  0.0   50.0   1013.0  288.0  4.2e-4    7.0e-3     0.0
```

---

## Consistency requirements

When preparing atmospheric data files, ensure that:

- the number of gas columns matches `NG`,
- the number of extinction columns matches `NW`,
- units are consistent with the lookup tables,
- the vertical range covers all ray paths used in the simulation.

Inconsistent atmospheric input is a common source of runtime errors or
unexpected results.

---

## Summary

The atmospheric data file defines the thermodynamic and compositional
state of the atmosphere used by JURASSIC. It is independent of the
observation geometry and radiance data and should be prepared carefully
to ensure physically meaningful and numerically stable results. The
ASCII format is convenient for small examples and inspection, while the
netCDF format is better suited to multi-profile workflows.

---

## Related pages

- [Configuration](configuration.md)
- [Observation and radiance files](radiance_data.md)
- [Spectroscopic data & LUTs](lookup_tables.md)
