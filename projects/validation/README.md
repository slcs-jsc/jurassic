# JURASSIC–RFM spectral validation

This directory contains a reproducible comparison of the current CPU JURASSIC
implementation with the RFM line-by-line model. It covers 500–2999 cm⁻¹ at
1 cm⁻¹ sampling with the mid-latitude climatology and 36 gases. Limb rays use
geometric tangent heights of 5, 10, 20, and 50 km; refraction is enabled.
Nadir and zenith each use one representative vertical ray.

The repository contains compact reference results:

- `rfm_reference/`: channel-matched RFM spectra, common inputs, timings, and provenance;
- `test_ega/`: JURASSIC EGA spectra and timings;
- `test_cga/`: JURASSIC CGA spectra and timings;
- `analysis/`: numerical accuracy metrics and publication-ready PNG figures.

Temporary spectral chunks, NetCDF files, and model logs are written below
`work/` and are not tracked by Git.

## Inspect and reproduce the analysis

The supplied spectra are sufficient to reproduce every accuracy metric and
figure without installing RFM or obtaining the external lookup tables:

```bash
python3 projects/validation/analyze.py
```

The analysis reports relative radiance differences for all limb channels with
nonzero RFM radiance. A true zero reference would be reported as undefined.
For nadir and zenith, the channel-averaged radiances are converted with
JURASSIC's inverse Planck definition and differences are reported in kelvin.

Main outputs are:

```text
analysis/REPORT.md
analysis/accuracy_metrics.csv
analysis/limb_radiance_spectra.png
analysis/limb_relative_errors.png
analysis/nadir_brightness_temperature.png
analysis/zenith_brightness_temperature.png
analysis/limb_accuracy_summary.png
analysis/brightness_temperature_accuracy_summary.png
analysis/runtime_summary.png
```

The limb spectra use a logarithmic radiance axis. Nadir and zenith are shown
separately as brightness temperatures with their absolute differences in
kelvin. Accuracy and runtime summaries are separate figures. `analysis/REPORT.md`
embeds all figures and tabulates the principal accuracy and timing results for
review.

## Repeat the calculations

Build the CPU executables with capacity for all 36 gases:

```bash
cd src
make clean
make DEFINES=-DNG=40 -j
```

Supply the external JURASSIC lookup-table directory with `JURASSIC_TBL_DIR` or
`--tbl-dir`.
It must contain `tria_<gas>.nc` for all gases recorded in the manifests. The
lookup tables are not part of this repository. All generated control files use
`TBLFMT 3`.

Generate the RFM reference first:

```bash
python3 projects/validation/run_rfm.py --force
```

Defaults can be overridden with `RFM_BIN`, `RFM_HIT`, and `RFM_XSC_DIR`, or
the corresponding command-line options. RFM, HITRAN data, cross sections, and
the instrumented RFM timing copy remain external dependencies.

Then repeat either JURASSIC calculation:

```bash
python3 projects/validation/run_ega.py --force
python3 projects/validation/run_cga.py --force
python3 projects/validation/analyze.py
```

Each calculation uses one model thread. By default two independent spectral
chunks run concurrently and are restricted to logical CPUs 0 and 2 on the
reference notebook. Set `VALIDATION_JOBS` and `VALIDATION_CPUSET` for another
machine; an empty `VALIDATION_CPUSET` disables affinity. Existing compact
results are replaced only with `--force`, and only after the new calculation
has completed successfully.

EGA and CGA regenerate the deterministic atmosphere and geometries and compare
them against `rfm_reference/input/` before accepting their results. Thus all
three methods use the same atmospheric state and viewing geometry. In RFM
mode, JURASSIC applies the spectral response stored in each TRIA channel to
the high-resolution RFM output; arbitrary monochromatic RFM samples are not
used as references. With these 1 cm^-1 tables, adjacent channel responses overlap
and each of the 20 chunks forms one contiguous RFM block. RFM is therefore
started 20 times per geometry.

## Timing interpretation

`timings.csv` records the summed single-process times across spectral chunks.
JURASSIC model time is `TIMER_FORMOD`. Instrumented RFM model time is the
reported path plus spectral phases minus output time. For limb, the reported
per-spectrum time is the four-ray calculation divided by four. Parallel jobs
reduce elapsed execution time but do not alter these summed single-process
model times.

`READ_TBL` includes reading and preparing the lookup tables. RFM timing keeps
its HITRAN initialization and binary-read measurements separately. Hardware,
thread count, external-data provenance, and the timing definition must remain
with any reported speed-up.
