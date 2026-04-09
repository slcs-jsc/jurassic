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
  - name: Yiran Zhang
    affiliation: 1
  - name: Florian Rahlmann
    affiliation: 2
    orcid: 0009-0006-4785-6406
  - name: Amir Nikfal
    affiliation: 1
    orcid: 0000-0002-6699-9473

affiliations:
  - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Jülich, Germany
    index: 1
  - name: Technische Universität Hamburg, Hamburg, Germany
    index: 2

date: 2026-01-07

bibliography: paper.bib
---

# Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source software package for simulating atmospheric radiances and transmittances in the mid-infrared spectral region (approximately 3--20 µm). It is designed for atmospheric remote sensing applications that require both physical realism and high computational throughput, such as satellite data analysis, sensitivity studies, retrieval development, and ensemble calculations. JURASSIC supports forward simulations for limb, nadir, zenith, and occultation viewing geometries and includes an integrated optimal estimation framework for inverse modelling of atmospheric state variables.

JURASSIC addresses a common need in atmospheric remote sensing. Line-by-line radiative transfer models provide high spectroscopic accuracy but are often too expensive for large datasets, repeated retrieval iterations, or operationally sized sensitivity studies. JURASSIC instead combines efficient spectral approximations with precomputed spectroscopic lookup tables, enabling fast calculations while preserving the accuracy required for many research applications. The software is implemented in C, supports OpenMP acceleration together with optional MPI-based task distribution for retrieval workflows, and is intended for efficient use on multicore workstations and high-performance computing systems.

# Statement of need

Infrared remote sensing observations are widely used to derive atmospheric temperature, trace-gas abundances, aerosols, and cloud-related quantities from satellite, airborne, balloon-borne, and ground-based instruments. These applications require radiative transfer software that is accurate enough to represent the relevant spectroscopy and viewing geometry, while also being fast enough for large observational datasets and iterative inverse methods.

JURASSIC is designed for researchers working in atmospheric infrared remote sensing, especially those developing retrieval algorithms, performing uncertainty and sensitivity analyses, or processing large ensembles of atmospheric scenarios. It provides a general-purpose, research-oriented framework for fast mid-infrared radiative transfer and retrieval that is not restricted to a single instrument or observing configuration.

# State of the field

Atmospheric infrared radiative transfer is supported by a mature ecosystem of software, including highly accurate line-by-line reference models and faster approximate forward models. Line-by-line tools are essential for spectroscopic reference calculations and benchmarking, but their runtime cost can become prohibitive when thousands to millions of forward calculations are required. JURASSIC occupies a distinct position within this landscape as a fast, research-oriented model for multiple atmospheric remote sensing geometries rather than a single observation mode or instrument family.

In addition to forward modelling, JURASSIC includes an optimal estimation retrieval framework, allowing radiative transfer and inversion to be performed within the same software environment. The core software in this repository has been used for limb and nadir sounding and for retrievals of temperature and trace gases, while the broader JURASSIC software family has also supported cloud- and aerosol-affected cases through dedicated extensions as well as tomographic and GPU-accelerated applications. These developments reflect a broader research scope than instrument-specific forward models tightly coupled to one observing system or operational product chain.

# Software design

JURASSIC is designed around a central trade-off between spectroscopic fidelity and computational efficiency. Instead of carrying out line-by-line calculations during runtime, it uses precomputed emissivity lookup tables derived from detailed line-by-line radiative transfer calculations. At runtime, these tables are accessed through fast interpolation, allowing band-averaged emissivities to be evaluated efficiently.

The radiative transfer formulation is based on a vertically stratified atmosphere and curved ray paths that account for atmospheric refraction. The software supports limb, nadir, zenith, and occultation viewing geometries for instruments located inside or outside the atmosphere. To approximate infrared absorption and emission efficiently, JURASSIC uses the Emissivity Growth Approximation and the Curtis--Godson Approximation [@Godson1953; @Gordley1981; @Marshall1994]. These approximations were chosen because they provide a practical balance between accuracy and speed across a broad range of atmospheric conditions relevant to remote sensing.

JURASSIC is implemented in modular components for ray tracing, spectroscopy, radiative transfer, and retrieval. This structure makes it easier to adapt individual parts of the workflow, reuse common functionality across applications, and maintain the code over time.

The software also reflects design choices motivated by large-scale computation. JURASSIC is written in C and supports OpenMP acceleration together with optional MPI-based task distribution for retrieval workflows, allowing efficient use of shared-memory systems and distributed execution of independent inversion cases. This matters for applications involving large observational archives, long time series, or computationally intensive inversion workflows. Beyond the core repository described here, a GPU-accelerated implementation has also been developed for heterogeneous architectures [@Baumeister2022], and recent integration work has aimed to connect the core, scattering-enabled, and accelerator-oriented variants within a more maintainable and extensible library-based code base, referred to as JURASSIC-UNIFIED [@Pozgaj2025].

JURASSIC is intentionally focused on infrared remote sensing applications in which local thermodynamic equilibrium is an adequate assumption and clear-air radiative transfer is the dominant use case. Under these conditions, molecular populations are determined by local temperature and emission is governed by the Planck function. The software has also been extensible enough to support dedicated treatments of cloud and aerosol effects, including simplified grey-body and extinction-based parameterizations in the core workflow and more advanced scattering extensions in related developments [@Griessbach2014; @Griessbach2016].

Finally, JURASSIC includes an optimal estimation framework for retrieving atmospheric quantities such as temperature and trace-gas volume mixing ratios directly from radiance measurements [@Rodgers2000]. Including retrieval capabilities in the same framework as the forward model is a deliberate design choice: it reduces duplication of model interfaces, simplifies reproducible inversion workflows, and supports direct sensitivity studies of forward-model assumptions within retrieval applications.

# Research impact statement

JURASSIC has been used in a broad range of atmospheric remote sensing studies over multiple generations of instruments and research questions. Its first limb-sounding application was demonstrated with measurements from the Michelson Interferometer for Passive Atmospheric Sounding (MIPAS), where it supported retrievals and climatological analyses of chlorofluorocarbons and other trace gases [@Hoffmann2008]. It was later applied to nadir observations from the Atmospheric Infrared Sounder (AIRS) for temperature retrievals and gravity-wave studies [@Hoffmann2009]. These early applications established the software as a practical tool for both limb and nadir infrared remote sensing.

Subsequent work expanded the scientific reach of the software. Dedicated developments added treatments for aerosols and clouds, including applications to volcanic aerosol observations and cloud-affected infrared limb measurements [@Griessbach2014; @Griessbach2016]. Tomographic extensions enabled retrievals of mesoscale atmospheric structures from airborne and satellite observations [@Ungermann2010; @Ungermann2012b; @Ungermann2012]. More recent work has demonstrated GPU acceleration for efficient execution on heterogeneous architectures [@Baumeister2022]. Together, these studies show that JURASSIC is not limited to a single measurement geometry or algorithmic niche, but has served as a reusable platform for radiative transfer and retrieval research.

The software also provides practical signals of community readiness and reproducibility. The code base includes automated tests and example configurations for representative limb and nadir applications. Reference cases are used to verify installation and numerical behaviour, and model results have been compared with established radiative transfer approaches in benchmarking and intercomparison studies [@Griessbach2013; @Baumeister2022]. Public source-code hosting, versioned releases, and a persistent Zenodo archive support citation and long-term reuse. These factors, together with the documented scientific applications above, provide evidence of realized impact rather than only anticipated utility.

\autoref{fig:spectra} illustrates example mid-infrared limb and nadir radiance spectra simulated with JURASSIC, including the contributions of selected trace gases for a mid-latitude reference atmosphere.

![Simulated mid-infrared limb (a) and nadir (b) spectra from JURASSIC for a mid-latitude atmosphere at 1 cm$^{-1}$ resolution. Coloured curves show single-gas contributions, and the black curve shows the total spectrum.\label{fig:spectra}](fig_spectra.pdf){ width=80% }

# AI usage disclosure

The authors used generative AI tools for limited assistance in software development, testing, and documentation as well as drafting and revising this manuscript. All AI-assisted content was reviewed, corrected where necessary, and validated by the authors for scientific accuracy and software correctness.

# Acknowledgements

We acknowledge contributions from the atmospheric remote sensing community for testing, validation, and scientific feedback. Computational and storage resources for development and applications of JURASSIC were provided by the Jülich Supercomputing Centre.

# References
