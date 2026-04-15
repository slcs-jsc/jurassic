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
    affiliation: 1
    orcid: 0000-0002-2005-4474
  - name: Stjepan Pozgaj
    affiliation: 2
    orcid: 0009-0002-4799-3911
  - name: Yiran Zhang
    affiliation: 1
    orcid: 0009-0009-9539-9598
  - name: Florian Rahlmann
    affiliation: 3
    orcid: 0009-0006-4785-6406
  - name: Amirhossein Nikfal
    affiliation: 1
    orcid: 0000-0002-6699-9473
  - name: Catrin Meyer
    affiliation: 1
    orcid: 0000-0002-9271-6174

affiliations:
  - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
    index: 1
  - name: University of Zagreb, Faculty of Electrical Engineering and Computing, Zagreb, Croatia
    index: 2
  - name: Technische Universität Hamburg, Hamburg, Germany
    index: 3

date: 2026-01-07

bibliography: paper.bib
---

# Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source software package for simulating atmospheric radiances and transmittances in the mid-infrared spectral region (approximately 3--20 µm in wavelength). It is designed for atmospheric remote sensing applications that require both physical realism and high computational throughput, such as satellite data analysis, retrieval development and processing, and large-scale sensitivity and ensemble studies. JURASSIC supports forward simulations for limb, nadir, zenith, and occultation viewing geometries and includes an integrated optimal estimation framework for inverse modelling of atmospheric state variables.

JURASSIC addresses a common need in atmospheric remote sensing. Line-by-line radiative transfer models provide high spectroscopic accuracy but are often too expensive for large satellite datasets, iterative retrieval workflows, or comprehensive sensitivity studies. JURASSIC instead combines efficient spectral approximations with precomputed spectroscopic lookup tables, enabling fast calculations while preserving the accuracy required for many research applications. The software is implemented in C, supports OpenMP acceleration together with optional MPI-based task distribution for retrieval workflows, and is intended for efficient use on multicore workstations and high-performance computing systems.

# Statement of need

Infrared remote sensing observations are widely used to derive atmospheric temperature, trace-gas abundances, aerosols, and cloud-related quantities from satellite, airborne, balloon-borne, and ground-based instruments. These applications require radiative transfer software that is accurate enough to represent the relevant spectroscopy and viewing geometry, while also being fast enough for large observational datasets and iterative inverse methods.

JURASSIC is designed for researchers working in atmospheric infrared remote sensing, especially those processing large ensembles of atmospheric scenarios, performing sensitivity and uncertainty analyses, or developing retrieval schemes. It provides a general-purpose, research-oriented framework for fast mid-infrared radiative transfer and retrieval that is not restricted to a single instrument or observing configuration.

# State of the field

Atmospheric infrared radiative transfer is supported by a mature ecosystem of software, including highly accurate line-by-line reference models [@Clough2005; @Dudhia2017; @Buehler2018] and faster approximate forward models [@Moncet2008; @Hocking2021]. Line-by-line tools are essential for spectroscopic reference calculations and benchmarking, but their runtime cost can become prohibitive when thousands to millions of forward calculations are required. JURASSIC occupies a distinct position within this landscape as a fast, research-oriented model for multiple atmospheric remote sensing geometries rather than a single observation mode or instrument family.

In addition to forward modelling, JURASSIC includes an optimal estimation retrieval framework, allowing radiative transfer and inversion to be performed within the same software environment. The core version of JURASSIC described here has been used for limb and nadir radiative transfer applications, including retrievals of atmospheric temperature and trace gas abundances. Related extensions, not covered in this paper, have also supported cloud- and aerosol-affected cases as well as tomographic and GPU-accelerated applications. Together, these developments reflect a broader research scope than that of instrument-specific forward models tightly coupled to a single observing system or operational product chain.

# Software design

JURASSIC is implemented in modular components for ray tracing, spectroscopy, radiative transfer, and retrieval. This modular design supports code reuse across applications, simplifies maintenance, and makes it easier to adapt individual components to specific research needs. The software is written in C and supports OpenMP acceleration together with optional MPI-based task distribution for retrieval workflows, enabling efficient use of multicore systems and distributed execution of independent retrieval cases. JURASSIC supports NetCDF-based input and output for atmospheric data, observation data, and lookup tables, improving interoperability with common scientific data workflows.

The radiative transfer formulation in the core version of JURASSIC described here assumes a vertically stratified atmosphere and curved ray paths that account for atmospheric refraction. The software supports limb, nadir, zenith, and occultation viewing geometries for sensors located inside or outside the atmosphere, corresponding to tangent-path, downward-looking, upward-looking, and transmission-style observation configurations, respectively.

To approximate infrared absorption and emission efficiently, JURASSIC uses the Emissivity Growth Approximation and the Curtis--Godson Approximation [@Godson1953; @Gordley1981; @Marshall1994]. In this framework, spectroscopic information is represented by precomputed emissivity lookup tables derived from detailed line-by-line radiative transfer calculations. During runtime, these tables are accessed through fast interpolation, allowing band-averaged emissivities to be evaluated efficiently. Together, these methods provide a practical balance between accuracy and speed across a broad range of atmospheric conditions relevant to infrared remote sensing.

The retrieval component implements optimal estimation methods for deriving atmospheric quantities such as temperature and trace-gas volume mixing ratios directly from radiance measurements [@Rodgers2000]. Integrating retrieval capabilities with the forward model reduces duplication of model interfaces, supports reproducible inversion workflows, and facilitates sensitivity studies of forward-model assumptions within retrieval applications.

The version of JURASSIC described in this paper is intentionally focused on infrared remote sensing applications in which local thermodynamic equilibrium is an adequate assumption and clear-air radiative transfer is the dominant use case. Under these conditions, molecular populations are determined by local temperature and emission is governed by the Planck function. The core workflow also supports simplified treatments of cloud and aerosol effects, including grey-body and extinction-based parameterizations.

Related developments outside the scope of this paper include more advanced scattering extensions for cloud- and aerosol-affected cases [@Griessbach2014; @Griessbach2016], a GPU-accelerated implementation for heterogeneous architectures [@Baumeister2022], and recent integration efforts aimed at connecting core, scattering-enabled, and accelerator-oriented variants within a more maintainable library-based framework referred to as JURASSIC-UNIFIED [@Pozgaj2022].

# Research impact statement

JURASSIC has been used in a broad range of atmospheric remote sensing studies spanning multiple instrument generations and scientific applications. Its early use in limb sounding was demonstrated with measurements from the Michelson Interferometer for Passive Atmospheric Sounding (MIPAS), where it supported retrievals and climatological analyses of chlorofluorocarbons and other trace gases [@Hoffmann2008]. It was later applied to nadir observations from the Atmospheric Infrared Sounder (AIRS) for temperature retrievals and gravity wave studies [@Hoffmann2009]. These applications established JURASSIC as a practical tool for both limb and nadir infrared remote sensing.

Subsequent developments building on JURASSIC expanded its scientific scope. Extensions for aerosols and clouds enabled applications to volcanic aerosol observations and cloud-affected infrared limb measurements [@Griessbach2014; @Griessbach2016]. Tomographic variants supported retrievals of mesoscale atmospheric structures from airborne and satellite observations [@Ungermann2010; @Ungermann2010b; @Ungermann2012]. More recent work demonstrated GPU acceleration for efficient execution on heterogeneous architectures [@Baumeister2022]. Together, these studies show that JURASSIC has served as a reusable platform for radiative transfer and retrieval research across multiple observing geometries and application settings.

The software also shows practical evidence of technical maturity and reproducibility. The code base includes automated tests and example configurations for representative limb and nadir applications. Reference cases are provided to verify installation and numerical behaviour, and model results have been compared with established radiative transfer approaches in benchmarking and intercomparison studies [@Griessbach2013; @Baumeister2022]. Public source-code hosting, versioned releases, and a persistent Zenodo archive further support citation, reproducibility, and long-term reuse. Together with the documented scientific applications above, these features indicate realized impact rather than only anticipated utility.

\autoref{fig:spectra} shows example mid-infrared limb and nadir radiance spectra simulated with JURASSIC for a mid-latitude reference atmosphere, including the contributions of selected trace gases.

![Simulated mid-infrared limb (a) and nadir (b) spectra from JURASSIC for a mid-latitude atmosphere at 1 cm$^{-1}$ resolution. Coloured curves show single-gas contributions, and the black curve shows the total spectrum.\label{fig:spectra}](fig_spectra.pdf){ width=80% }

# AI usage disclosure

The authors used generative AI tools for limited assistance in software development, testing, and documentation as well as drafting and revising this manuscript. All AI-assisted content was reviewed, corrected where necessary, and validated by the authors for scientific accuracy and software correctness.

# Acknowledgements

We acknowledge contributions from the atmospheric remote sensing community for testing, validation, and scientific feedback. Computational and storage resources for development and applications of JURASSIC were provided by the Jülich Supercomputing Centre.

# References
