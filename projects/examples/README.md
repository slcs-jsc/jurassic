# JURASSIC examples

These examples provide small, self-contained introductions to three common
observation geometries:

- `limb/` observes atmospheric tangent paths and writes radiances.
- `nadir/` observes the atmosphere from above and writes brightness
  temperatures.
- `zenith/` observes upward from the surface and writes brightness
  temperatures.

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

Each script generates an atmosphere (`atm.tab`), observation geometry
(`obs.tab`), radiative-transfer result (`rad.tab`), kernel functions
(`kernel.tab`), and diagnostic PNG plots. It finishes by comparing `rad.tab`
exactly with the checked-in `rad.org` reference and returns a nonzero status if
they differ.

The examples use the small lookup tables under `tests/data/` and are intended
for installation checks, tutorials, and configuration experiments. Numerical
validation across broader spectral and gas configurations is provided separately
under `projects/validation/`; performance measurements are under
`projects/benchmark/`.

The scripts overwrite their generated TAB and PNG outputs in place. Run
`make check` from `src/` for the regression suite.
