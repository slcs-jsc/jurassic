# Spectroscopic lookup tables (LUTs)

JURASSIC uses **precomputed spectroscopic lookup tables (LUTs)** to
represent infrared absorption and emission efficiently. These tables
are a central component of the model and enable fast radiative transfer
calculations without performing line-by-line spectroscopy at runtime.

This page describes the role, structure, and configuration of LUTs in
JURASSIC.

---

## Purpose of lookup tables

Infrared radiative transfer depends on highly structured molecular
spectra with millions of individual absorption lines. Evaluating these
lines explicitly is computationally expensive and impractical for
large-scale applications.

JURASSIC avoids this cost by using lookup tables that:

- are generated **offline** using detailed line-by-line models,
- store **band-averaged emissivities or transmittances**,
- are interpolated efficiently during runtime.

This approach preserves high spectroscopic fidelity while enabling
orders-of-magnitude speedups compared to line-by-line calculations.

---

## Spectral approximations

The lookup tables are designed to support the spectral approximations
used by JURASSIC, in particular:

- the **Emissivity Growth Approximation (EGA)**,
- the **Curtis–Godson Approximation (CGA)**.

The tables provide emissivity values for homogeneous atmospheric layers
as functions of pressure, temperature, and absorber amount. These
values are then combined along an inhomogeneous ray path using EGA.

Details of the underlying theory are provided in the
[Theory](theory.md) section.

---

## Table dimensions and variables

Each lookup table is defined on a multi-dimensional grid. Typical
dimensions include:

- **pressure**  
- **temperature**  
- **absorber column density** (or equivalent quantity)

Separate tables are usually provided for each:

- molecular species (emitter),
- spectral window or band.

The exact grid definition depends on how the tables were generated and
must be compatible with the atmospheric conditions encountered during
runtime.

---

## Table organization and naming

Lookup tables are located using the control parameter:

- `TBLBASE`

`TBLBASE` defines a base path or prefix from which the model derives the
full table file names.

The current table I/O implementation expects flat filename patterns
derived directly from `TBLBASE`; the naming is not arbitrary.

---

## Table formats

JURASSIC supports multiple lookup table formats. The format is selected
via the control keyword:

- `TBLFMT`

Common formats include:

- **ASCII tables** (human-readable, larger, slower to read),
- **binary tables** (compact, fast I/O, preferred for production use).
- **netCDF tables** (portable container format; tables for one emitter are
  grouped into a single file with one packed variable per channel).

For the current table I/O implementation, the naming conventions are:

- ASCII: one file per channel and emitter, for example
  `TBLBASE_792.0000_CO2.tab`
- binary: one file per channel and emitter, for example
  `TBLBASE_792.0000_CO2.bin`
- netCDF: one file per emitter, for example `TBLBASE_CO2.nc`, containing
  one packed variable per channel such as `tbl_792.0000` and `tbl_832.0000`

The example projects included with JURASSIC demonstrate working table
layouts and formats and are the recommended starting point for new
users.

For ASCII tables, the corresponding filter-function files are also
looked up via fixed names of the form `TBLBASE_<nu>.filt`.

---

## Runtime interpolation

During execution, JURASSIC performs fast interpolation within the lookup
tables to obtain emissivities for the local atmospheric state.

Interpolation is typically:

- linear in pressure, temperature, and absorber amount,
- performed independently for each gas and spectral window.

The interpolation routines are optimized for performance and are called
frequently during ray-path integration.

---

## Validity and extrapolation

Lookup tables are valid only within the parameter ranges for which they
were generated. If runtime atmospheric conditions fall outside these
ranges, interpolation may result in:

- extrapolated values,
- reduced accuracy,
- or, in extreme cases, invalid results.

**Best practice:**  
Ensure that lookup tables cover the full range of pressure, temperature,
and gas concentrations expected in your application.

---

## Generating lookup tables

The underlying spectroscopic content for lookup tables is typically
prepared **outside** the core JURASSIC runtime, for example with a
line-by-line radiative transfer model and a chosen spectroscopic
database.

While JURASSIC is agnostic to the specific external spectroscopy tool,
the repository does include utilities that support table workflows:

- `tblgen` for generating table files from prepared spectroscopy inputs,
- `tblfmt` for converting between ASCII, binary, and netCDF formats,
- `tblrdc` for reducing or repacking table files.

These tools are part of the JURASSIC distribution, but the upstream
line-by-line spectroscopy setup and application-specific table design
remain external responsibilities.

In practice, JURASSIC lookup tables are often generated with the
[Reference Forward Model (RFM) from the University of Oxford](https://eodg.atm.ox.ac.uk/RFM/).
The control parameters `RFMBIN`, `RFMHIT`, and `RFMXSC[i]` support
RFM-based workflows where applicable.

---

## Performance considerations

- Binary lookup tables significantly reduce I/O overhead.
- netCDF lookup tables can be convenient when multiple channels for one
  emitter should be bundled into a single file.
- Keeping tables on fast local or parallel file systems improves
  performance for large runs.
- Table reuse across many simulations amortizes the offline generation
  cost.

For HPC workflows, careful table management is essential for achieving
good scalability.

---

## Summary

Spectroscopic lookup tables are a fundamental building block of
JURASSIC. They enable fast infrared radiative transfer simulations by
decoupling expensive spectroscopic calculations from runtime
execution.

Correct configuration and careful generation of lookup tables are
crucial for both accuracy and performance.

---

## Related pages

- [Configuration](configuration.md)
- [Theory](theory.md)
- [Model assumptions & limitations](limitations.md)
