# JURASSIC Numerical Validation

This directory validates numerical results independently of the performance
benchmarks in `projects/benchmark/`. Benchmark runners measure runtime only. A code
change should pass validation before its benchmark numbers are accepted.

## Coverage

`configs/validation_cases.tsv` is the single source of truth. Every row defines an
executable validation case and its tolerance. The current matrix provides:

- `smoke`: zenith, nadir, and limb with 125 channels from 500 to 2980 cm^-1
- `full`: zenith, nadir, and limb for all integer 1-cm channels in the four ranges
  500–900, 900–1500, 1500–2200, and 2200–2999 cm^-1; bands are split into
  contiguous chunks of at most 128 channels to respect the compiled `ND` limit
- the seven-gas `core` set for smoke cases and the complete 36-gas LUT inventory
  for full cases
- medium ray sets matching the established geometry definitions

Adjacent full bands include their shared boundary channel deliberately. Chunks within
a band do not overlap or leave gaps. Add new
gas sets or atmospheric scenarios to this manifest rather than hiding them inside a
runner script.

## Frozen References

References must be generated deliberately from a trusted scalar CPU build. The
generator refuses a dirty worktree by default and never overwrites a reference
unless both an explicit `--case` and `--replace` are supplied.

```bash
cd projects/validation
VALIDATION_TBLBASE=/path/to/tria scripts/generate_references.py --profile smoke
VALIDATION_TBLBASE=/path/to/tria scripts/generate_references.py --profile full
# Run independent scalar reference cases concurrently:
VALIDATION_TBLBASE=/path/to/tria scripts/generate_references.py --profile full --jobs 6
# Resume an interrupted matrix without replacing completed cases:
VALIDATION_TBLBASE=/path/to/tria scripts/generate_references.py --profile full --missing
```

`--jobs` defaults to one. Each worker still runs `formod` with
`EXECUTION scalar`; only independent validation cases are generated concurrently.
Allow about 3 GiB of memory per worker for the 128-channel full cases and select the
worker count with sufficient memory and CPU headroom.

Each reference directory contains the canonical control file, frozen atmospheric
inputs, NetCDF observation geometry (`obs.nc`), double-precision NetCDF scalar
radiances (`rad.nc`), logs, and `metadata.json`. Metadata records
the reference commit, the resolved LUT base path, and SHA-256 checksums for the
inputs, output, and control file. The three smoke-reference directories and their reference-spectrum PNGs are
versioned with the repository, so a fresh clone can run and inspect smoke validation
immediately. Their large generation logs are excluded. Full-reference directories
remain ignored because the complete matrix is substantially larger; archive and
distribute the full reference tree as a separate versioned validation artifact.

Audit completeness and provenance without running the forward model:

```bash
scripts/check_references.py --profile smoke
VALIDATION_TBLBASE=/path/to/tria scripts/check_references.py --profile full
```

Adding or changing a manifest row makes the audit fail until its trusted reference
is explicitly generated again.

## Candidate Validation

The default run is the quick smoke profile using the current `src/formod` binary:

```bash
VALIDATION_TBLBASE=/path/to/tria scripts/run_validation.py
```

Run the complete spectral and geometry matrix before accepting an optimization:

```bash
cd src
make clean
make DEFINES=-DNG=40 -j
cd ../projects/validation
VALIDATION_TBLBASE=/path/to/tria scripts/run_validation.py --profile full
```

Full validation assumes that JURASSIC was compiled with capacity for all configured
emitters, for example with `DEFINES=-DNG=40`. The validation scripts intentionally
do not manage or inspect this compile-time setting. Smoke validation continues to
use the seven-gas `core` set and works with the default `NG=8` build.

Use `--formod` to test a separately built executable, `--execution` to select its
scalar or batch path, and `--formod-method cga` to compare the Curtis-Godson
approximation against the frozen default EGA references. `OBSFMT = 3` keeps reference and candidate results as NetCDF
doubles. The pass/fail gate compares all science variables (`bt_*` for
zenith/nadir and `rad_*` for limb) plus all `tau_*` variables directly. Paired NaNs
are accepted, while one-sided NaNs fail. `OBSREF` remains in the log as supplementary
diagnostics. The run fails for missing/stale references, changed frozen inputs,
missing science variables, incompatible NaN masks, or a tolerance violation. Results
are written below `runs/` as per-case logs, radiances, JSON accuracy reports, and a
combined `summary.tsv`.

## Plots

Generate readable channel-wise error spectra for an explicit run, a configured
default, or the latest available validation run. Separate Brightness Temperature
(or limb Radiance) and transmittance panels show MeanRE, SDRE, MinRE, and MaxRE over all
rays at every channel, together with the positive and negative tolerance. The
`summary_plot_ae.png` file shows MeanAE, SDAE, MinAE, and MaxAE spectra in radiance,
kelvin, or dimensionless transmittance units. The main `summary_plot.png` combines
limb-radiance RE with brightness-temperature and transmittance AE; the all-relative
plot is retained as `summary_plot_re.png`:

```bash
scripts/plot_summary.py runs/validation_<timestamp>
VALIDATION_RUN_DIR=runs/validation_<timestamp> scripts/plot_summary.py
scripts/plot_summary.py
```

Override the three output paths independently when needed:

```bash
scripts/plot_summary.py runs/validation_<timestamp> \
  --combined-output summary_plot.png \
  --output summary_plot_re.png \
  --absolute-output summary_plot_ae.png
```

The AE statistics retain the sign of `candidate - reference`; “absolute” distinguishes
native-unit errors from relative percentages. Thus `MeanAE` is the bias rather than
the mean absolute magnitude.

Generate full reference spectra and transmittance plots for all geometries:

```bash
scripts/plot_full_reference_spectrum.py
```

This creates separate Brightness Temperature plots for zenith/nadir, a Radiance
plot for limb, and one additional `*_tau.png` transmittance plot per geometry.

Generate the corresponding plots for the quick smoke profile separately:

```bash
scripts/plot_full_reference_spectrum.py --profile smoke
```

Smoke outputs use the `reference_smoke_*` prefix, while full outputs retain the
`reference_full_*` prefix. This keeps both plot sets unambiguous even when they are
generated into the same `plots/` directory.

Plot science and transmittance errors for one validation case:

```bash
scripts/plot_channels.py runs/validation_<timestamp>/<case_name>
```
