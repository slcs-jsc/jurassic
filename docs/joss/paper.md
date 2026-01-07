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
    orcid: 0000-0002-1529-1075
    affiliation: 1
  - name: Paul F. Baumeister
    orcid: 0000-0003-0123-4567
    affiliation: 1
affiliations:
  - name: Jülich Supercomputing Centre, Forschungszentrum Jülich, Germany
    index: 1
date: 2026-01-07
bibliography: paper.bib
---

## Summary

The Jülich Rapid Spectral Simulation Code (JURASSIC) is an open-source infrared radiative transfer model for simulating atmospheric radiances and transmittances in the mid-infrared spectral region. The model is designed to balance computational efficiency and physical accuracy, enabling the processing of large datasets and ensembles typical of modern atmospheric remote sensing applications. JURASSIC supports both forward simulations and inverse modelling through an integrated optimal estimation framework, making it suitable for satellite radiance simulation, sensitivity studies, and retrieval algorithm development.

## Statement of Need

Infrared remote sensing instruments provide critical observations of atmospheric temperature, trace gases, aerosols, and clouds. The increasing data volume from modern and upcoming satellite missions places strong demands on radiative transfer models, which must be both computationally efficient and scientifically reliable. While line-by-line radiative transfer models offer high accuracy, their computational cost limits their applicability for large-scale or near-real-time applications.

Fast radiative transfer models address this challenge by using spectral approximations and precomputed spectroscopic information. However, many existing tools are tailored to specific instruments, spectral ranges, or operational contexts, limiting their flexibility for research-oriented workflows. JURASSIC addresses this gap by providing a general-purpose, modular infrared radiative transfer framework that supports multiple observation geometries, customizable spectral configurations, and high-performance computing environments. Its open-source distribution further promotes transparency, reproducibility, and community-driven development.

## Software Description

### Model architecture and design

JURASSIC is implemented in the C programming language with a modular design that separates ray tracing, radiative transfer, spectroscopy, and retrieval components. The atmosphere is represented as a vertically stratified medium, and radiative transfer calculations are performed along curved ray paths that account for atmospheric refraction. The model supports limb, nadir, zenith, and occultation geometries for instruments located inside or outside the atmosphere.

Spectral absorption and emission are represented using precomputed lookup tables derived from detailed line-by-line calculations. During runtime, band-averaged emissivities are obtained through fast interpolation, enabling rapid evaluation of the radiative transfer equation without sacrificing spectroscopic fidelity.

### Radiative transfer and retrieval capabilities

JURASSIC implements established spectral approximations, including the Emissivity Growth Approximation and the Curtis–Godson Approximation, to efficiently model infrared radiative transfer. These methods allow accurate simulations of atmospheric radiances and transmittances across a broad spectral range while avoiding the computational expense of full line-by-line calculations.

In addition to forward modelling, JURASSIC includes an optimal estimation retrieval module for inverse modelling of atmospheric state variables. This enables the retrieval of geophysical parameters such as temperature and trace gas volume mixing ratios directly from measured radiances within the same software framework.

A comprehensive description of the underlying algorithms, numerical methods, and GPU-accelerated implementation is provided by Baumeister and Hoffmann (2022). :contentReference[oaicite:0]{index=0}

### Performance and scalability

JURASSIC supports hybrid MPI–OpenMP parallelization and is designed for efficient execution on multicore CPUs and high-performance computing clusters. The model has also been extended to GPU-accelerated architectures, demonstrating substantial speedups and improved energy efficiency compared to CPU-only implementations. These capabilities enable the processing of large observational datasets, global simulations, and long time series within practical time constraints.

## Quality Control

The JURASSIC codebase includes automated tests and example configurations to verify correct installation and numerical behaviour. Example projects for limb and nadir geometries are provided and can be executed as part of the test workflow. Model outputs are validated against reference data, and extensive benchmarking and intercomparison studies with established radiative transfer models have been documented in the literature.

Continuous integration practices, versioned releases, and persistent digital object identifiers are used to support reproducibility and long-term usability.

## Availability

- **Source code:** https://github.com/slcs-jsc/jurassic  
- **License:** GNU General Public License v3.0  
- **Persistent identifier:** https://doi.org/10.5281/zenodo.4572889  

## Acknowledgements

The development of JURASSIC has benefited from contributions by numerous collaborators and from computational resources provided by the Jülich Supercomputing Centre. We acknowledge support from Forschungszentrum Jülich and the broader atmospheric remote sensing community for testing, validation, and scientific feedback.

## References

