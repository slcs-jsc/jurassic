# Atmospheric data

This page describes the atmospheric input data used by JURASSIC,
typically provided in a file named `atm.tab`. The atmospheric data file
defines the thermodynamic and compositional state of the atmosphere used
for radiative transfer and retrieval calculations.

In JURASSIC, the atmospheric data represent the **model state**, while
observation geometry and radiances are handled separately (see
*Observation and radiance files*).

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

The description below refers to the **default ASCII format**, typically
used with:

```text
ATMFMT = 1
```

This is the format used in the example projects distributed with
JURASSIC.

---

## General format rules

- The file is a plain-text ASCII table.
- Columns are separated by spaces or tabs.
- Header lines are optional; if present, they must not start with a
  valid numeric value in the first column.
- Each row corresponds to one vertical grid point.

---

## Units and conventions

The default atmospheric reader expects the following units:

| Quantity            | Unit |
|---------------------|------|
| time                | s    |
| altitude (`z`)      | km   |
| longitude, latitude | deg  |
| pressure (`p`)      | hPa  |
| temperature (`t`)   | K    |
| volume mixing ratio | mol/mol |
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
**first atmospheric level only** may include additional columns such as:

- cloud-layer altitude
- cloud-layer thickness
- cloud extinction coefficients

If these options are not used, the corresponding parameters should be
disabled in the control file and omitted from the atmospheric data.

---

## Optional: surface parameters

Surface-related parameters (e.g. surface temperature or emissivity) may
also be specified in the atmospheric data file for certain applications.

When enabled, these parameters are expected in the **first row only**.
If surface effects are not required, the default configuration should
be used and these fields omitted.

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

The atmospheric data file (`atm.tab`) defines the thermodynamic and
compositional state of the atmosphere used by JURASSIC. It is independent
of the observation geometry and radiance data and should be prepared
carefully to ensure physically meaningful and numerically stable
results.

---

## Related pages

- [Configuration](configuration.md)
- Observation and radiance files (radiance_data.md)
- [Spectroscopic data & LUTs](lookup_tables.md)
