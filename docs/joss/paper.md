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

affiliations:
  - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
    index: 1

date: 2026-01-07

bibliography: paper.bib
---

# Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source infrared radiative transfer model for simulating atmospheric radiances and transmittances in the mid-infrared spectral region (wavelengths of approximately 3–20 µm). The model is designed to balance computational efficiency and physical accuracy, enabling the processing of large datasets and ensembles typical of modern atmospheric remote sensing applications. JURASSIC supports both forward simulations and inverse modelling through an integrated optimal estimation framework, making it suitable for satellite, airborne, balloon-borne, and ground-based infrared remote sensing use cases, as well as for sensitivity studies and retrieval algorithm development.

# Statement of Need

Infrared remote sensing instruments provide comprehensive observations of atmospheric temperature, trace gases, aerosols, and clouds. The increasing data volume from modern satellite missions places strong demands on radiative transfer models, which must be both computationally efficient and scientifically reliable. While line-by-line radiative transfer models offer high accuracy, their computational cost limits their applicability for large-scale or near-real-time applications.

Fast radiative transfer models address this challenge by using spectral approximations and precomputed spectroscopic information, but many existing tools are tailored to specific instruments or operational contexts. JURASSIC addresses this gap by providing a general-purpose, modular infrared radiative transfer framework with support for multiple observation geometries, customizable spectral configurations, and high-performance computing environments. Its open-source development model promotes transparent and reproducible research workflows.

# Radiative Transfer and Retrieval Capabilities

In JURASSIC, the atmosphere is represented as a vertically stratified medium, and radiative transfer calculations are performed along curved ray paths that account for atmospheric refraction. The model supports limb, nadir, zenith, and occultation geometries for instruments located inside or outside the atmosphere.

JURASSIC employs two spectral approximations, the Emissivity Growth Approximation (EGA) and the Curtis–Godson Approximation (CGA), to efficiently model infrared radiative transfer [@Godson1953; @Gordley1981; @Marshall1994]. These methods enable accurate simulations of atmospheric radiances and transmittances across a broad spectral range while avoiding the computational expense of full line-by-line calculations.

Spectral absorption and emission in JURASSIC are represented using precomputed emissivity lookup tables derived from detailed line-by-line radiative transfer calculations. These tables can be generated offline using established line-by-line models, such as the Reference Forward Model (RFM) [@Dudhia2017], but are not tied to a specific implementation. During runtime, band-averaged emissivities are obtained through fast linear interpolation within the lookup tables, enabling computationally efficient radiative transfer calculations while preserving high spectroscopic fidelity.

In addition to forward modelling, JURASSIC features an optimal estimation retrieval module for inverse modelling of atmospheric state variables [@Rodgers2000]. This enables the retrieval of geophysical parameters such as temperature and trace gas volume mixing ratios directly from measured radiances within the same software framework.

# Applications

JURASSIC has been applied in a wide range of atmospheric remote sensing studies for infrared limb and nadir measurements. The first limb-sounding application of JURASSIC was demonstrated using Envisat Michelson Interferometer for Passive Atmospheric Sounding (MIPAS) observations, which were used to derive stratospheric distributions and climatologies of chlorofluorocarbons and other trace gases [@Hoffmann2008]. The first nadir application was subsequently demonstrated using Atmospheric Infrared Sounder (AIRS) radiance measurements for temperature retrievals and gravity wave studies [@Hoffmann2009].

\autoref{fig:spectra} illustrates example mid-infrared limb and nadir radiance spectra simulated with JURASSIC, including the contributions of selected trace gases for a mid-latitude reference atmosphere.

![Simulated mid-infrared limb (a) and nadir (b) spectra from JURASSIC for a mid-latitude atmosphere at 1 cm$^{-1}$ resolution. Coloured curves show single-gas contributions, and the black curve shows the total spectrum.\label{fig:spectra}](fig_spectra.pdf){ width=80% }

# Model Assumptions and Extensions

JURASSIC assumes local thermodynamic equilibrium, so molecular level populations are determined by the local atmospheric temperature and follow the Boltzmann distribution. Under this assumption, infrared emission is governed by the Planck function at the local temperature, which is valid for most infrared remote sensing applications in the troposphere and stratosphere.

The core JURASSIC model is designed for clear-air conditions, where scattering of infrared radiation can be neglected. Cloud and aerosol effects can be represented using simplified grey-body or extinction-based parameterizations suitable for many stratospheric and upper-tropospheric applications. More advanced treatments of infrared scattering are not part of the core distribution but have been implemented in dedicated extensions that account for single and multiple scattering by aerosols and clouds and have been applied, for example, to volcanic aerosol observations and cloud-affected infrared limb measurements [@Griessbach2014; @Griessbach2016].

Further extensions of JURASSIC include tomographic retrieval capabilities for mesoscale atmospheric structures from airborne and satellite measurements [@Ungermann2010; @@Ungermann2012b; Ungermann2012] as well as a GPU-accelerated implementation for efficient execution on heterogeneous computing architectures [@Baumeister2022]. More recently, an integration effort has focused on combining the core version, scattering-enabled, and GPU-accelerated variants into a single unified code base, referred to as [JURASSIC-UNIFIED](https://github.com/slcs-jsc/jurassic-unified), with the aim of improving maintainability, extensibility, and long-term sustainability through a modular, library-oriented software design [@Pozgaj2025].

# Software Description

JURASSIC is implemented in the C programming language with a modular design that separates ray tracing, radiative transfer, spectroscopy, and retrieval components. JURASSIC supports hybrid MPI–OpenMP parallelization and is designed for efficient execution on multicore CPUs and high-performance computing clusters. These capabilities enable the processing of large observational datasets, global simulations, and long time series within practical time constraints.

The JURASSIC codebase includes automated tests and example configurations to verify correct installation and numerical behaviour. Example projects for limb and nadir geometries are provided and can be executed as part of the test workflow. Model outputs are validated against reference data, and benchmarking and intercomparison studies with established radiative transfer models have been documented in the literature [@Griessbach2013; @Baumeister2022]. Versioned releases and persistent digital object identifiers support reproducibility and long-term usability.

User and developer documentation is provided through a comprehensive README, example projects, and online manuals hosted with the [source code repository](https://github.com/slcs-jsc/jurassic), facilitating installation, configuration, and extension of the model. JURASSIC is distributed as open-source software under the GNU General Public License v3.0, and a persistent, citable archive of the source code is available via [Zenodo](https://doi.org/10.5281/zenodo.4572889).

# Acknowledgements

We acknowledge contributions from the atmospheric remote sensing community for testing, validation, and scientific feedback. Computational and storage resources for development and applications of JURASSIC were provided by the Jülich Supercomputing Centre.

# References
