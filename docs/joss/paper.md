---
title: "JURASSIC: A fast infrared radiative transfer model for atmospheric remote sensing"

tags:
  - radiative transfer
  - infrared remote sensing
  - atmospheric modelling
  - high-performance computing
  - retrieval algorithms

authors:
  - name: Lars Hoffmann
    affiliation: 1
    orcid: 0000-0003-3773-4377
  - name: Sabine Griessbach
    affiliation: 1
    orcid: 0000-0003-3792-3573
  - name: Paul F. Baumeister
    orcid: 0000-0002-2005-4474
    affiliation: 1

affiliations:
  - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
    index: 1

date: 2026-01-07

bibliography: paper.bib
---

# Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source infrared radiative transfer model for simulating atmospheric radiances and transmittances in the mid-infrared spectral region. The model is designed to balance computational efficiency and physical accuracy, enabling the processing of large datasets and ensembles typical of modern atmospheric remote sensing applications. JURASSIC supports both forward simulations and inverse modelling through an integrated optimal estimation framework, making it suitable for satellite radiance simulation, sensitivity studies, and retrieval algorithm development.

# Statement of Need

Infrared remote sensing instruments provide critical observations of atmospheric temperature, trace gases, aerosols, and clouds. The increasing data volume from modern satellite missions places strong demands on radiative transfer models, which must be both computationally efficient and scientifically reliable. While line-by-line radiative transfer models offer high accuracy, their computational cost limits their applicability for large-scale or near-real-time applications.

Fast radiative transfer models address this challenge by using spectral approximations and precomputed spectroscopic information, but many existing tools are tailored to specific instruments or operational contexts. JURASSIC addresses this gap by providing a general-purpose, modular infrared radiative transfer framework that supports multiple observation geometries, customizable spectral configurations, and high-performance computing environments, while being openly available to promote transparency and reproducibility.

# Software Description

## Model architecture and design

JURASSIC is implemented in the C programming language with a modular design that separates ray tracing, radiative transfer, spectroscopy, and retrieval components. The atmosphere is represented as a vertically stratified medium, and radiative transfer calculations are performed along curved ray paths that account for atmospheric refraction. The model supports limb, nadir, zenith, and occultation geometries for instruments located inside or outside the atmosphere. Spectral absorption and emission are represented using precomputed lookup tables derived from detailed line-by-line calculations, allowing rapid evaluation of band-averaged emissivities during runtime.

## Radiative transfer and retrieval capabilities

JURASSIC implements established spectral approximations, including the Emissivity Growth Approximation (EGA) and the Curtis–Godson Approximation (CGA), to efficiently model infrared radiative transfer [@Godson1953; @Gordley1981; @Marshall1994]. These methods allow accurate simulations of atmospheric radiances and transmittances across a broad spectral range while avoiding the computational expense of full line-by-line calculations.

In addition to forward modelling, JURASSIC includes an optimal estimation retrieval module for inverse modelling of atmospheric state variables. This enables the retrieval of geophysical parameters such as temperature and trace gas volume mixing ratios directly from measured radiances within the same software framework.

A comprehensive description of the underlying algorithms, numerical methods, and GPU-accelerated implementation is provided by @Baumeister2022.

\autoref{fig:spectra} illustrates example mid-infrared limb and nadir radiance spectra simulated with JURASSIC, including the contributions of selected trace gases for a mid-latitude reference atmosphere.

![Simulated mid-infrared radiance spectra for limb (a) and nadir (b) viewing geometries calculated with JURASSIC for mid-latitude atmospheric conditions at 1 cm$^{-1}$ spectral resolution.\label{fig:spectra}](fig_spectra.pdf)

## Assumptions and treatment of clouds and aerosols

In its standard configuration, JURASSIC is designed primarily for clear-air conditions, where scattering of infrared radiation can be neglected and local thermodynamic equilibrium is assumed. Cloud and aerosol effects can be represented using simplified grey-body or extinction-based parameterizations suitable for many stratospheric and upper-tropospheric applications. More advanced treatments of infrared scattering are not part of the core distribution but have been implemented in dedicated extensions of JURASSIC that account for single and multiple scattering by aerosols and clouds and have been applied, for example, to volcanic aerosol observations and cloud-affected infrared limb measurements [@Griessbach2013; @Griessbach2016].

## Applications

JURASSIC has been applied in a wide range of atmospheric remote sensing studies, particularly for infrared limb and nadir observations. Early applications focused on trace gas retrievals from satellite measurements, including Envisat Michelson Interferometer for Passive Atmospheric Sounding (MIPAS) observations used to derive stratospheric distributions and climatologies of chlorofluorocarbons and other species [@Hoffmann2008]. The model has also been used for temperature retrievals from nadir-viewing infrared sounders, such as Atmospheric Infrared Sounder (AIRS) radiance measurements applied to gravity wave studies [@Hoffmann2009]. JURASSIC has also been extended for tomographic retrievals of mesoscale atmospheric structures from airborne and satellite measurements, supporting studies of gravity waves and stratospheric dynamics [@Ungermann2010; @Ungermann2012].

## Performance and scalability

JURASSIC supports hybrid MPI–OpenMP parallelization and is designed for efficient execution on multicore CPUs and high-performance computing clusters. The model has also been extended to GPU-accelerated architectures, demonstrating substantial speedups and improved energy efficiency compared to CPU-only implementations [@Baumeister2022]. These capabilities enable the processing of large observational datasets, global simulations, and long time series within practical time constraints.

# Quality Control

The JURASSIC codebase includes automated tests and example configurations to verify correct installation and numerical behaviour. Example projects for limb and nadir geometries are provided and can be executed as part of the test workflow. Model outputs are validated against reference data, and extensive benchmarking and intercomparison studies with established radiative transfer models have been documented in the literature. Continuous integration practices, versioned releases, and persistent digital object identifiers are used to support reproducibility and long-term usability.

# Availability

- **Source code:** https://github.com/slcs-jsc/jurassic  
- **License:** GNU General Public License v3.0  
- **Persistent identifier:** https://doi.org/10.5281/zenodo.4572889  

# Acknowledgements

We acknowledge contributions from numerous collaborators within the atmospheric remote sensing community for testing, validation, and scientific feedback. The Jülich Supercomputing Centre provided computational and storage resources.

# References
