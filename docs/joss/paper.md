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
  - name: Stjepan Požgaj
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

Atmospheric infrared radiative transfer is supported by a mature ecosystem of software, including line-by-line reference models [@Clough2005; @Dudhia2017; @Buehler2018] and fast or parameterized forward models [@Moncet2008; @Hocking2021]. Comparative studies have examined both classes for satellite sounding and retrieval applications, including radiances, transmittances, and Jacobians [@Saunders2007; @Schreier2018; @vonClarmann2002]. JURASSIC occupies a complementary role within this landscape as a research-oriented framework for fast thermal-infrared radiative transfer and retrieval across multiple observation geometries.

In addition to forward modelling, JURASSIC includes an optimal estimation retrieval framework, allowing radiative transfer and inversion to be performed within the same software environment. The core version of JURASSIC described here has been used for limb and nadir radiative transfer applications, including retrievals of atmospheric temperature and trace gas abundances. Related extensions, not covered in this paper, have also supported cloud- and aerosol-affected cases as well as tomographic and GPU-accelerated applications. Together, these developments illustrate the use of JURASSIC as a reusable research framework across different infrared remote sensing configurations.

# Software design

JURASSIC is implemented in modular components for ray tracing, spectroscopy, radiative transfer, and retrieval. This modular design supports code reuse across applications, simplifies maintenance, and makes it easier to adapt individual components to specific research needs. The software is written in C and supports OpenMP acceleration together with optional MPI-based task distribution for retrieval workflows, enabling efficient use of multicore systems and distributed execution of independent retrieval cases. JURASSIC supports NetCDF-based input and output for atmospheric data, observation data, and lookup tables, improving interoperability with common scientific data workflows.

The radiative transfer formulation in the core version of JURASSIC described here assumes a vertically stratified atmosphere and curved ray paths that account for atmospheric refraction. The software supports limb, nadir, zenith, and occultation viewing geometries for sensors located inside or outside the atmosphere, corresponding to tangent-path, downward-looking, upward-looking, and transmission-style observation configurations, respectively. The core model assumes local thermodynamic equilibrium and focuses primarily on clear-air thermal-infrared radiative transfer, while also supporting simplified grey-body and extinction-based treatments of clouds and aerosols.

To approximate infrared absorption and emission efficiently, JURASSIC uses the Emissivity Growth Approximation (EGA) and the Curtis--Godson Approximation (CGA) [@Godson1953; @Gordley1981; @Marshall1994]. JURASSIC applies a band transmittance approach in which radiative transfer is evaluated using quantities averaged over the spectral response of each instrument channel rather than monochromatically. The corresponding precomputed emissivity lookup tables contain these band-averaged emissivities and are derived from detailed line-by-line calculations. During runtime, the tables are evaluated by interpolation, enabling inexpensive repeated forward model evaluations. The main trade-offs are the memory required for the lookup tables and the restriction of the calculations to the atmospheric and spectral ranges represented by those tables.

The accuracy of EGA, CGA, and related band model approximations depends on spectral interval, channel response, atmospheric state, and viewing geometry, with errors typically at the percent level or below in different applications [@Gordley1981; @Marshall1994; @Francis2006; @Baumeister2022]. The JURASSIC repository contains a validation test case with 2500 channels and 36 gases for a mid-latitude atmosphere over 500--2999 cm$^{-1}$. Against the Reference Forward Model (RFM), median absolute relative limb radiance differences are 0.13--1.15% for EGA and 0.25--1.11% for CGA, while RMS brightness temperature differences across the nadir and zenith cases range from 0.41 to 0.56 K. On a single-thread Intel Core i7-1365U CPU, JURASSIC is 100--790 times faster than RFM for EGA and 170--1290 times faster for CGA.

The retrieval component implements optimal estimation methods for deriving atmospheric quantities such as temperature and trace-gas volume mixing ratios directly from radiance measurements [@Rodgers2000]. Integrating retrieval capabilities with the forward model reduces duplication of model interfaces, supports reproducible inversion workflows, and facilitates sensitivity studies of forward model assumptions within retrieval applications.

Jacobians required by the retrieval are calculated using finite differences. This approach keeps the derivative calculation general and closely coupled to the forward model, at the cost of additional forward model evaluations for perturbed state vector elements. These calculations are independent and are parallelized with OpenMP. For retrieval workloads containing many independent cases, JURASSIC additionally supports MPI-based task distribution. Computational resources can therefore be used either for shared-memory acceleration of an individual retrieval or for distributing independent retrieval cases, depending on the workload.

Related developments outside the scope of this paper include more advanced scattering extensions for cloud- and aerosol-affected cases [@Griessbach2014; @Griessbach2016], a GPU-accelerated implementation for heterogeneous architectures [@Baumeister2022], and recent integration efforts aimed at connecting core, scattering-enabled, and accelerator-oriented variants within a more maintainable library-based framework referred to as JURASSIC-UNIFIED [@Pozgaj2022].

# Research impact statement

JURASSIC has been used in a broad range of atmospheric remote sensing studies spanning multiple instrument generations and scientific applications. Its early use in limb sounding was demonstrated with measurements from the Michelson Interferometer for Passive Atmospheric Sounding (MIPAS), where it supported retrievals and climatological analyses of chlorofluorocarbons and other trace gases [@Hoffmann2008]. It was later applied to nadir observations from the Atmospheric Infrared Sounder (AIRS) for temperature retrievals and gravity wave studies [@Hoffmann2009]. These applications established JURASSIC as a practical tool for both limb and nadir infrared remote sensing.

Subsequent developments building on JURASSIC expanded its scientific scope. Extensions for aerosols and clouds enabled applications to volcanic aerosol observations and cloud-affected infrared limb measurements [@Griessbach2014; @Griessbach2016]. Tomographic variants supported retrievals of mesoscale atmospheric structures from airborne and satellite observations [@Ungermann2010; @Ungermann2010b; @Ungermann2012]. More recent work demonstrated GPU acceleration for efficient execution on heterogeneous architectures [@Baumeister2022]. Together, these studies show that JURASSIC has served as a reusable framework for radiative transfer and retrieval research across multiple observing geometries and application settings.

The software also shows practical evidence of technical maturity and reproducibility. The code base includes automated tests and example configurations for representative limb and nadir applications. Reference cases are provided to verify installation and numerical behaviour, and model results have been compared with established radiative transfer approaches in benchmarking and intercomparison studies [@Griessbach2013; @Baumeister2022]. Public source code hosting, versioned releases, and a persistent Zenodo archive further support citation, reproducibility, and long-term reuse. Together with the documented scientific applications above, these features indicate realized impact rather than only anticipated utility.

\autoref{fig:spectra} shows example mid-infrared limb and nadir radiance spectra simulated with JURASSIC for a mid-latitude reference atmosphere, including the contributions of selected trace gases.

![Simulated mid-infrared limb (a) and nadir (b) spectra from JURASSIC for a mid-latitude atmosphere at 1 cm$^{-1}$ resolution. Coloured curves show single-gas contributions, and the black curve shows the total spectrum.\label{fig:spectra}](fig_spectra.pdf){ width=80% }

# AI usage disclosure

The authors used generative AI tools for limited assistance in software development, testing, and documentation as well as drafting and revising this manuscript. All AI-assisted content was reviewed, corrected where necessary, and validated by the authors for scientific accuracy and software correctness.

# Acknowledgements

We acknowledge contributions from the atmospheric remote sensing community for testing, validation, and scientific feedback. Computational and storage resources for development and applications of JURASSIC were provided by the Jülich Supercomputing Centre.

# References
