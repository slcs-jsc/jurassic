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
  - name: Luis Niewöhner
    affiliation: 3
    orcid: 0009-0004-2244-9842
  - name: Florian Rahlmann
    affiliation: 4
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
  - name: Leibniz University Hannover, Hannover, Germany
    index: 3
  - name: Technische Universität Hamburg, Hamburg, Germany
    index: 4

date: 2026-01-07

bibliography: paper.bib
---

# Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open source software package for simulating atmospheric radiances and transmittances in the mid-infrared spectral region (approximately 3--20 µm). It is designed for atmospheric remote sensing applications that require both physical realism and high computational throughput, such as satellite data analysis, retrieval development and processing, and sensitivity and ensemble studies. JURASSIC supports forward simulations for limb, nadir, zenith, and occultation viewing geometries and includes an integrated optimal estimation framework for retrieving atmospheric state variables.

JURASSIC addresses a common need in atmospheric remote sensing. Line-by-line radiative transfer models provide high spectroscopic accuracy but are often too expensive for large satellite datasets, iterative retrieval workflows, or extensive sensitivity studies. JURASSIC instead combines efficient radiative transfer approximations with precomputed emissivity lookup tables, enabling fast calculations while preserving the accuracy required for many research applications. The software is implemented in C, supports OpenMP acceleration and optional MPI task distribution for retrieval workflows, and is intended for efficient use on multicore workstations and high-performance computing systems.

# Statement of need

Infrared remote sensing observations are widely used to derive atmospheric temperature, trace gas abundances, aerosols, and cloud properties from satellite, airborne, balloon-borne, and ground-based instruments. These applications require radiative transfer software that is accurate enough to represent the relevant spectroscopy and viewing geometry, while also being fast enough for large observational datasets and iterative inverse methods.

JURASSIC is designed for researchers working in atmospheric infrared remote sensing, especially those processing large ensembles of atmospheric scenarios, performing sensitivity and uncertainty analyses, or developing retrieval schemes. It provides a research framework for fast mid-infrared radiative transfer and retrieval that is not restricted to a single instrument or observing configuration.

# State of the field

Atmospheric infrared radiative transfer is supported by a wide range of models, including line-by-line reference models [@Clough2005; @Dudhia2017; @Buehler2018] and fast or parameterized forward models [@Moncet2008; @Hocking2021]. Comparative studies have examined both classes for satellite sounding and retrieval applications, including radiances, transmittances, and Jacobians [@Saunders2007; @Schreier2018; @vonClarmann2002]. These models offer different trade-offs between spectroscopic generality, computational cost, supported observing configurations, and adaptation to new applications. JURASSIC occupies a complementary role as a research framework for fast thermal-infrared radiative transfer and retrieval across multiple observation geometries.

In addition to forward modelling, JURASSIC includes finite-difference Jacobians and an optimal estimation retrieval framework, allowing radiative transfer and inversion within the same software environment. JURASSIC is particularly useful when high computational throughput is required for a well-defined instrument or spectral configuration while retaining flexibility in viewing geometry, atmospheric state, and retrieval setup. Its computational efficiency relies on precomputed emissivity lookup tables and radiative transfer approximations. For a new instrument or channel configuration, suitable lookup tables need to be generated and the approximation accuracy validated against a line-by-line reference model. Established fast model systems may be preferable for supported operational sensor configurations, whereas line-by-line models remain the natural choice when maximum spectroscopic generality or reference model accuracy is required.

The core version of JURASSIC has been used for limb and nadir radiative transfer applications, including retrievals of atmospheric temperature and trace gas abundances [@Hoffmann2008; @Hoffmann2009]. Related extensions, not covered in this paper, have supported tomographic retrievals for novel limb imaging concepts [@Ungermann2010; @Ungermann2010b], cases affected by clouds and aerosols [@Griessbach2013; @Griessbach2014; @Griessbach2016], and applications using GPUs [@Baumeister2022]. More recent work on JURASSIC-UNIFIED aims to integrate the core model with scattering and GPU variants into a common library framework [@Pozgaj2022].

# Software design

JURASSIC is organized into modular components for ray tracing, spectroscopy, radiative transfer, and retrieval. This modular design facilitates code reuse, maintenance, and adaptation to different applications. The software is written in C and supports OpenMP parallelization and optional MPI task distribution for retrieval workflows. JURASSIC supports NetCDF input and output for atmospheric data, observation data, and lookup tables.

The radiative transfer formulation assumes a vertically stratified atmosphere and curved ray paths accounting for atmospheric refraction. The software supports limb, nadir, zenith, and occultation viewing geometries for sensors located inside or outside the atmosphere. \autoref{fig:spectra} shows example mid-infrared limb and nadir radiance spectra for a mid-latitude reference atmosphere, including contributions from selected trace gases.

![Simulated mid-infrared limb (a) and nadir (b) spectra from JURASSIC for a mid-latitude atmosphere at 1 cm$^{-1}$ resolution. Coloured curves show contributions from individual gases, and the black curve shows the total spectrum.\label{fig:spectra}](fig_spectra.pdf){ width=80% }

The core model assumes local thermodynamic equilibrium and focuses primarily on clear-air thermal infrared radiative transfer, while also supporting simplified grey-body and extinction-based treatments of clouds and aerosols. At the lower atmospheric boundary, surface temperature and, optionally, spectrally varying emissivity determine an isotropic thermal surface source; optional treatments add reflected downwelling thermal radiation and a simplified solar contribution. The core model does not include solar scattering by cloud and aerosol particles.

To approximate infrared absorption and emission efficiently, JURASSIC uses the Emissivity Growth Approximation (EGA) and the Curtis--Godson Approximation (CGA) [@Godson1953; @Gordley1981; @Marshall1994]. Instead of monochromatic calculations, JURASSIC uses channel-averaged transmittances based on precomputed emissivity lookup tables derived from line-by-line spectra integrated over the instrument spectral response functions. Interpolation of these tables enables efficient repeated forward model evaluations. The main trade-offs are memory requirements and limitations imposed by the tabulated atmospheric and spectral ranges.

Comparisons with line-by-line calculations for different satellite applications demonstrated relative radiance differences below one percent for JURASSIC, together with reductions in computational cost of about three orders of magnitude [@Hoffmann2008; @Hoffmann2009]. Further studies evaluated the spectral averaging approach against KOPRA [@Griessbach2013] and investigated GPU performance [@Baumeister2022].

For the present work, reproducible comparisons with RFM over the 500–3000 cm⁻¹ spectral range show median absolute relative limb radiance differences of about 0.1–1.2% and nadir/zenith RMS brightness temperature differences of about 0.4–0.6 K. Depending on viewing geometry and approximation, JURASSIC is roughly two to three orders of magnitude faster than RFM in single-threaded benchmarks. OpenMP benchmarks on one 64-core AMD EPYC 7742 socket achieve about 50× speedup at 64 threads, corresponding to about 80% parallel efficiency. Detailed benchmark configurations, results, and reproduction instructions are provided in the repository.

The retrieval component implements optimal estimation methods for deriving atmospheric quantities from radiance measurements [@Rodgers2000]. Its state vector can include pressure, temperature, trace gas volume mixing ratios, spectrally varying extinction, cloud parameters, surface temperature, and spectrally varying surface emissivity. Jacobians are calculated using finite differences, a simple and flexible approach that requires additional forward model evaluations. These calculations are independent and parallelized with OpenMP. For retrieval workloads containing many independent cases, JURASSIC additionally supports MPI task distribution. The framework also supports reproducible retrievals, sensitivity studies, and observing system simulation experiments (OSSEs) using synthetic radiances.

A reproducible end-to-end retrieval validation provided in the repository exercises the finite-difference Jacobians and the optimal estimation inversion for a synthetic limb retrieval of temperature and trace gases. The reference case converges in three iterations, substantially reduces profile and radiance residual errors, and provides averaging kernels, posterior uncertainties, and runtime diagnostics.

The code base includes automated tests and documented example configurations for representative limb, nadir, and zenith applications. Reference cases allow users to verify installation and numerical behaviour. Public source code hosting, versioned releases, and a persistent Zenodo archive support citation, reproducibility, and long-term reuse.

# Research impact statement

JURASSIC has been used in atmospheric remote sensing studies spanning multiple instrument generations and scientific applications. Its early use in limb sounding was demonstrated with measurements from the Michelson Interferometer for Passive Atmospheric Sounding (MIPAS), where it supported retrievals and climatological analyses of chlorofluorocarbons and other trace gases [@Hoffmann2008]. It was later applied to nadir observations from the Atmospheric Infrared Sounder (AIRS) for temperature retrievals and gravity wave studies [@Hoffmann2009].

Further applications used JURASSIC to generate synthetic AIRS radiances from atmospheric model simulations, enabling direct comparisons of observed and simulated convective and orographic gravity waves [@Grimsdell2010; @Orr2015; @Kruse2022]. JURASSIC radiative transfer calculations also helped quantify gravity wave temperature amplitudes associated with the formation of polar stratospheric clouds over Antarctica [@Eckermann2009]. High-resolution AIRS temperature retrievals were subsequently used to compare the characteristics of stratospheric gravity waves observed by AIRS and the High Resolution Dynamics Limb Sounder (HIRDLS) [@Meyer2018].

Further developments expanded its scientific scope. Extensions for aerosols and clouds enabled applications to volcanic aerosol observations and infrared limb measurements affected by clouds [@Griessbach2014; @Griessbach2016]. Radiative transfer simulations also supported the analysis of volcanic emissions in AIRS observations [@Hoffmann2014]. These capabilities were also applied to investigate aerosol and cloud top heights derived from MIPAS observations [@Griessbach2020]. Tomographic variants supported retrievals of mesoscale atmospheric structures from airborne and satellite observations [@Ungermann2010; @Ungermann2010b; @Ungermann2012]. More recent work demonstrated GPU acceleration [@Baumeister2022].

# AI usage disclosure

We used OpenAI's ChatGPT and Codex, as well as Anthropic's Claude, to assist with software development, code review, testing, documentation, and manuscript revision. All AI-assisted outputs were reviewed, edited, and validated by the authors. The authors made all core scientific and software design decisions and remain fully responsible for the submitted software and manuscript.

# Acknowledgements

We acknowledge contributions from the atmospheric remote sensing community for testing, validation, and scientific feedback. Computational and storage resources for development and applications of JURASSIC were provided by the Jülich Supercomputing Centre.

# References
