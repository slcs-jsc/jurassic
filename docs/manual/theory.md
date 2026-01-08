# Physical background and theory

This page provides an overview of the physical background and numerical
methods implemented in JURASSIC. The description focuses on the
radiative transfer formulation, spectral approximations, and numerical
integration schemes that enable fast and accurate infrared simulations.
A detailed and formal description of the algorithms can be found in
Baumeister and Hoffmann (2022).

---

## Atmospheric representation

JURASSIC represents the atmosphere as a vertically stratified medium.
Atmospheric state variables are defined on discrete altitude levels and
include:

- pressure  
- temperature  
- volume mixing ratios of trace gases  
- optional aerosol or cloud extinction coefficients  

Between the specified levels, atmospheric quantities are interpolated
linearly (pressure in logarithmic space). For geometric calculations,
the Earth is approximated as a sphere with a fixed mean radius.

---

## Ray tracing and refraction

Radiative transfer calculations are performed along curved ray paths
that account for atmospheric refraction. Refraction causes bending of
the ray paths toward the Earth’s surface and is particularly important
for limb-sounding geometries.

Ray paths are calculated numerically by solving the Eikonal equation,
which describes the propagation of electromagnetic waves in a medium
with spatially varying refractive index. In the mid-infrared spectral
range, the refractive index depends primarily on pressure and
temperature; wavelength dependence and water vapour contributions are
small and neglected in the default formulation.

The ray tracing algorithm is generic and supports limb, nadir, zenith,
and occultation geometries for instruments located inside or outside
the atmosphere. The ray path is discretized into segments, which form
the basis for the subsequent radiative transfer calculations.

---

## Monochromatic radiative transfer equation

At the most fundamental level, radiative transfer in JURASSIC is based
on the monochromatic radiative transfer equation along a ray path:

- absorption attenuates radiation along the path,
- thermal emission contributes according to the local temperature.

Under the assumption of local thermodynamic equilibrium (LTE) and in
the absence of scattering, the source function reduces to the Planck
function evaluated at the local temperature. This formulation is valid
for most infrared applications in the troposphere and stratosphere.

A full monochromatic (line-by-line) evaluation of the radiative transfer
equation would require explicit treatment of millions of spectral lines
and is therefore computationally prohibitive for large-scale
applications.

---

## Band transmittance approximation

Satellite instruments measure radiances integrated over finite spectral
bands rather than at individual wavenumbers. JURASSIC exploits this by
using a **band transmittance approximation**, in which spectrally
averaged quantities are used instead of monochromatic ones.

In this approach:

- the Planck function is averaged over the instrument spectral response,
- the atmospheric transmissivity is represented by band-averaged
  emissivities,
- radiative transfer is performed using these spectrally averaged
  quantities.

This approximation significantly reduces computational cost but
introduces errors due to neglected spectral correlations. These errors
are minimized through the use of advanced emissivity modelling
techniques.

---

## Emissivity lookup tables

Spectral absorption and emission are represented in JURASSIC using
precomputed lookup tables of band-averaged emissivities. These tables
are generated offline using detailed line-by-line radiative transfer
models and spectroscopic databases.

Each lookup table provides emissivity as a function of:

- pressure  
- temperature  
- absorber column density  

During runtime, emissivities are obtained through fast linear
interpolation within the tables. This approach preserves much of the
accuracy of line-by-line spectroscopy while enabling orders-of-magnitude
faster calculations.

The validity of the radiative transfer results depends on the coverage
and resolution of the lookup tables.

---

## Emissivity Growth Approximation (EGA)

The **Emissivity Growth Approximation (EGA)** is a key component of
JURASSIC and is used to compute the effective emissivity of an
inhomogeneous atmospheric path.

Instead of treating each atmospheric segment independently, EGA
accounts for the cumulative growth of emissivity along the ray path:

1. The current total path emissivity is mapped to an equivalent
   “pseudo-column density” on the emissivity curve of the next segment.
2. The emissivity increase due to the additional segment is then
   computed by advancing along this curve.

This procedure allows the model to approximate the emissivity of an
inhomogeneous atmosphere using lookup tables derived for homogeneous
layers. The EGA method substantially reduces errors associated with
simple band transmittance approaches and has been shown to provide good
accuracy for a wide range of infrared remote sensing applications.

---

## Numerical integration along the line of sight

Radiative transfer along the ray path is evaluated by dividing the path
into discrete segments defined by the ray tracing algorithm. For each
segment, atmospheric properties are assumed to be homogeneous.

The radiance at the instrument is obtained by iteratively accumulating
the emission from each segment while accounting for absorption by the
overlying atmosphere. This procedure follows the Beer–Lambert law and
ensures numerical stability and efficiency.

The integration scheme is designed to be compatible with both forward
modelling and the calculation of derivatives required for retrieval
applications.

---

## Multiple absorbers and continua

JURASSIC supports multiple molecular absorbers and includes continuum
absorption and emission processes for selected gases. The total
emissivity is constructed from the individual gas contributions under a
continuum approximation, which neglects higher-order spectral
correlation terms between different species.

This approximation is accurate for many broadband and moderate-resolution
applications but represents a known limitation for very narrow spectral
channels or strongly overlapping absorption features.

---

## Relation to retrieval calculations

The radiative transfer formulation described above forms the basis for
both forward simulations and inverse modelling in JURASSIC. The same
numerical framework is used to compute radiances and their sensitivities
with respect to atmospheric state variables.

This consistency between forward model and Jacobian calculations is a
key requirement for optimal estimation retrievals and ensures numerical
stability and physical coherence of the retrieval results.

---

## Summary

JURASSIC combines physically sound radiative transfer theory with
carefully chosen approximations to achieve high computational
performance. By using curved ray tracing, emissivity lookup tables, and
the emissivity growth approximation, the model enables fast and accurate
infrared simulations suitable for large-scale atmospheric remote
sensing and retrieval applications.

For a comprehensive and formal derivation of the algorithms, users are
referred to Baumeister and Hoffmann (2022) and the references listed in
the [References](references.md) section.
