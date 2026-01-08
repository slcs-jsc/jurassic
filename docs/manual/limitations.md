# Model assumptions and limitations

This page summarizes the main physical and technical assumptions made
by JURASSIC and discusses the resulting limitations. Understanding
these assumptions is important for correctly interpreting simulation
and retrieval results and for assessing the suitability of JURASSIC
for a given application.

---

## Local thermodynamic equilibrium (LTE)

JURASSIC assumes **local thermodynamic equilibrium (LTE)** throughout
the atmosphere. Under this assumption, molecular level populations are
determined solely by the local atmospheric temperature and follow the
Boltzmann distribution. Infrared emission is therefore described by
the Planck function at the local temperature.

This assumption is valid for most infrared remote sensing applications
in the troposphere and stratosphere. In the upper mesosphere and
thermosphere, non-LTE effects can become important and are **not**
represented in the core JURASSIC model.

**Implication:**  
JURASSIC should not be used for applications where non-LTE emission or
absorption plays a significant role unless dedicated extensions are
employed.

---

## Spectral approximations

JURASSIC relies on fast spectral approximations to achieve high
computational efficiency. In particular, the model uses:

- the Emissivity Growth Approximation (EGA)
- the Curtis–Godson Approximation (CGA)

These approximations replace explicit line-by-line calculations during
runtime with band-averaged representations based on precomputed
lookup tables.

While these methods have been extensively validated and shown to be
accurate for many infrared applications, they are approximations by
construction.

**Implication:**  
Small deviations relative to line-by-line models may occur, especially
for very narrow spectral features or highly non-uniform atmospheric
conditions. JURASSIC is therefore best suited for broadband or
moderate-resolution infrared simulations rather than ultra-high
spectral resolution studies.

---

## Lookup table dependence

Spectroscopic information in JURASSIC is provided via precomputed
lookup tables derived from detailed line-by-line radiative transfer
calculations.

The accuracy of JURASSIC simulations depends on:

- the quality of the underlying spectroscopic database,
- the line-by-line model used to generate the tables,
- the spectral resolution and parameter coverage of the tables.

Lookup tables are generated offline and are assumed to be valid within
the parameter ranges for which they were created.

**Implication:**  
Simulations outside the pressure, temperature, or gas concentration
ranges covered by the lookup tables may lead to extrapolation errors
or invalid results. Users should ensure that the lookup tables are
appropriate for their application.

---

## Scattering and cloud effects

The core JURASSIC model is designed for **clear-air conditions**, where
scattering of infrared radiation can be neglected.

Clouds and aerosols can be represented using simplified parameterizations,
such as grey-body or extinction-based approaches, which are suitable for
many stratospheric and upper-tropospheric applications.

More sophisticated treatments of infrared scattering, including
single- and multiple-scattering by clouds and aerosols, are **not**
part of the core distribution but have been implemented in dedicated
extensions.

**Implication:**  
Applications strongly affected by cloud or aerosol scattering require
careful validation or the use of extended model versions. The core
JURASSIC model is not intended as a full multiple-scattering infrared
radiative transfer solver.

---

## Atmospheric representation

The atmosphere in JURASSIC is represented as a vertically stratified
medium with horizontally homogeneous layers along each ray path.

Refraction is accounted for, but horizontal gradients within a single
ray path are not explicitly resolved.

**Implication:**  
Strong three-dimensional atmospheric variability (e.g. sharp horizontal
gradients or small-scale structures) may not be fully captured unless
advanced techniques such as tomographic retrievals are employed.

---

## Instrument and geometry assumptions

JURASSIC supports a wide range of observation geometries, including
limb, nadir, zenith, and occultation configurations. Instrument
characteristics are represented through configurable spectral bands
and detector definitions.

Instrument effects such as spectral response functions are represented
in a simplified manner and must be carefully configured by the user.

**Implication:**  
Users are responsible for ensuring that the instrument configuration
used in a simulation accurately reflects the characteristics of the
real instrument being modeled.

---

## Numerical and technical limitations

As a fast radiative transfer model optimized for performance, JURASSIC
makes design choices that favor efficiency and scalability:

- finite vertical and spectral resolution
- numerical tolerances chosen for stability and speed
- reliance on precomputed data products

Parallel execution using MPI and OpenMP introduces additional
considerations related to numerical reproducibility across different
hardware and process layouts.

**Implication:**  
Bitwise-identical results across different platforms or parallel
configurations are not guaranteed, although numerical differences are
typically negligible for scientific applications.

---

## Summary

JURASSIC is a fast, flexible, and well-validated infrared radiative
transfer model designed for large-scale atmospheric remote sensing
applications. Its assumptions and approximations are appropriate for a
wide range of scientific and operational use cases, but they must be
kept in mind when interpreting results.

Users requiring non-LTE physics, full multiple-scattering treatments,
or ultra-high spectral resolution should consider complementary
specialized models or validated JURASSIC extensions.
