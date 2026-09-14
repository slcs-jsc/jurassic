# JURASSIC examples

These examples provide small, self-contained introductions showing how to run
radiative transfer calculations with JURASSIC for three common observation
geometries:

- `limb/` observes atmospheric tangent paths and writes radiances.
- `nadir/` observes the atmosphere from above and writes brightness
  temperatures.
- `zenith/` observes upward from the surface and writes brightness
  temperatures.

## Installation

Build JURASSIC and its bundled libraries from the repository root before running
an example:

```bash
cd libs
bash build.sh
cd ../src
make -j
cd ..
```

Then run an example from its own directory:

```bash
cd projects/examples/nadir
./run.sh
```

## Using precompiled binaries

The examples normally use the locally built executables in `src/`. To use an
extracted precompiled development package instead, set `JURASSIC_BIN` to its
`bin` directory:

```bash
git clone https://github.com/slcs-jsc/jurassic.git
cd jurassic

export JURASSIC_BIN=/path/to/jurassic-linux-x86_64-<commit>/bin

cd projects/examples/nadir
./run.sh
```

The repository is still needed for the example configurations, lookup tables,
scripts, and reference data, but JURASSIC itself does not need to be compiled
when `JURASSIC_BIN` points to the packaged executables. Without `JURASSIC_BIN`,
the scripts use the normal source build in `src/`. Gnuplot is still required for
the diagnostic plots.

Each script generates an atmosphere (`atm.tab`), observation geometry
(`obs.tab`), radiative-transfer result (`rad.tab`), kernel functions
(`kernel.tab`), and diagnostic PNG plots. It finishes by comparing `rad.tab`
exactly with the checked-in `rad.org` reference and returns a nonzero status if
they differ.

## Limb example

The limb case views long, nearly horizontal paths through atmospheric tangent
points at different altitudes. This geometry provides strong vertical
sensitivity to temperature and trace gases by resolving how their contributions
to the measured radiance vary with altitude.

```bash
cd projects/examples/limb
./run.sh
```

![Limb radiance profiles](limb/plot_rad.png)

*Simulated limb radiances at 792 and 832 cm<sup>-1</sup> as a function of
tangent altitude.*

The radiance profiles show the simulated signal at 792 and 832 cm<sup>-1</sup>
as a function of tangent height. Their channel-dependent shapes are physically
plausible: lowering the tangent point increases the absorbing and emitting
column, while temperature and gas abundance also vary along the path. The
result is not expected to be a simple monotonic curve.

![Limb temperature kernel at 792 cm-1](limb/plot_kernel_temperature_792.png)

*Temperature kernel profiles for the limb geometry at 792 cm<sup>-1</sup>.*

The temperature kernel is the change in radiance caused by a small temperature
change at each altitude. Each coloured profile belongs to a different tangent
height. Sensitivity is concentrated near and above the corresponding tangent
region because the ray does not sample lower altitudes, while absorption limits
how deeply an optically thick ray can sense.

## Nadir example

The nadir case looks down from a satellite across a latitude scan. It
demonstrates how thermal emission emerging from the atmosphere is represented
as brightness temperature in three nearby channels.

```bash
cd projects/examples/nadir
./run.sh
```

![Nadir brightness temperatures](nadir/plot_rad.png)

*Simulated nadir brightness temperatures in three nearby spectral channels
across the latitude scan.*

Brightness temperature is the temperature a blackbody would need to reproduce
the simulated channel radiance. Differences among the three curves are
reasonable because each channel has different CO2 absorption and therefore
samples a different range of atmospheric levels. Variation across the scan is
also expected as the viewing path becomes more slanted away from nadir.

![Nadir temperature kernel at 668.5410 cm-1](nadir/plot_kernel_temperature_668.5410.png)

*Temperature kernel profiles for the nadir scan at 668.5410 cm<sup>-1</sup>.*

This kernel shows how brightness temperature responds to temperature changes
with altitude for the different scan positions. Peaks mark the layers that
contribute most strongly; their displacement with viewing geometry reflects
the changing path length and optical depth.

## Zenith example

The zenith case observes upward from the surface through the atmosphere at a
range of viewing angles. It illustrates downwelling atmospheric emission and
the transition from a short vertical path to longer slant paths.

```bash
cd projects/examples/zenith
./run.sh
```

![Zenith brightness temperatures](zenith/plot_rad.png)

*Simulated zenith-viewing brightness temperatures at 792 and 832 cm<sup>-1</sup>
across the angular scan.*

The brightness-temperature curves are approximately symmetric around the
vertical view because opposite viewing directions traverse equivalent model
atmospheres. Toward more oblique angles, the longer atmospheric path increases
absorption and emission, so changes in brightness temperature and decreasing
transmittance are physically expected.

![Zenith temperature kernel at 792 cm-1](zenith/plot_kernel_temperature_792.0000.png)

*Temperature kernel profiles for the zenith geometry at 792 cm<sup>-1</sup>.*

The temperature kernel identifies the atmospheric layers controlling the
downwelling signal. More oblique rays generally place greater weight on the
lower, denser atmosphere, while the vertical ray samples a shorter column. The
kernel structure is a sensitivity diagnostic, not an independent accuracy
test.

## Scope

The examples use the small set of lookup tables under `tests/data/` and are
intended for installation checks, tutorials, and configuration experiments.
Numerical validation across broader spectral and gas configurations is provided
separately
under `projects/validation/`; performance measurements are under
`projects/benchmark/`.

The scripts overwrite their generated TAB and PNG outputs in place. Run
`make check` from `src/` for the regression suite.
