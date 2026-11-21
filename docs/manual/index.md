# Welcome to JURASSIC!

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source
infrared radiative transfer and retrieval model for atmospheric remote
sensing. It simulates radiances and transmittances in the mid-infrared
(approximately 3–20 µm) and supports applications ranging from satellite
data analysis to retrieval development and sensitivity studies.

![Infrared limb spectra for mid-latitude atmospheric conditions simulated with JURASSIC](img/spectra.png)

## Key features

- **Radiative transfer:** Fast simulations for limb, nadir, zenith, and
  occultation viewing geometries using spectral approximations and
  precomputed emissivity lookup tables.
- **Atmospheric retrievals:** Integrated optimal estimation framework
  with finite-difference Jacobians for atmospheric state variables.
- **Flexible applications:** Configurable atmospheric inputs, spectral
  channels, and instrument geometries.
- **Parallel computing:** OpenMP acceleration and optional MPI task
  distribution for retrieval workflows.
- **Open source:** C implementation released under the GNU GPL v3.

For the model's assumptions, accuracy, and applicable spectral
configurations, see [Model assumptions and limitations](limitations.md)
and [Scientific validation](validation.md).

## Getting started

Begin with the [Quickstart](quickstart.md) for a first simulation, or
consult [Installation](installation.md) and
[Configuration files](configuration.md) for setup details.

## Explore the documentation

- **Scientific background:** [Overview](overview.md),
  [Physical background](theory.md), and [Retrieval theory](retrieval_theory.md).
- **Using the model:** [Running JURASSIC](running.md),
  [Atmospheric data](atmospheric_data.md), and
  [Spectroscopic data and lookup tables](lookup_tables.md).
- **HPC and development:** [HPC workflows](hpc_workflows.md),
  [Code development](code-development.md), and
  [Testing and code quality](testing.md).

## Citation and contributions

Please consult [References](references.md) for scientific publications
and the [contribution guidelines](https://github.com/slcs-jsc/jurassic/blob/master/CONTRIBUTING.md)
if you would like to report an issue or contribute.

## Contact

Dr. Lars Hoffmann, <l.hoffmann@fz-juelich.de>  
Jülich Supercomputing Centre, Forschungszentrum Jülich, Germany
