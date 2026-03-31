# Quickstart

This Quickstart guides you through running the provided example
simulations that come with JURASSIC. The goal is to verify that the
model is correctly installed and to familiarize yourself with the basic
workflow of a JURASSIC radiative transfer simulation.

## Prerequisites

This guide assumes that:

-   JURASSIC has been successfully compiled
-   All required runtime dependencies are available
-   You are working in a shell environment on a Linux or HPC system

Detailed installation instructions are provided in the
[Installation](installation.md) section of the user manual.

------------------------------------------------------------------------

## Running the example simulations

JURASSIC includes a set of example projects that demonstrate typical use
cases and serve as regression tests. These projects are located in the
`projects` directory.

``` bash
cd jurassic/projects
```

Three example observation geometries are provided:

-   **Nadir sounding**, representative of satellite instruments viewing
    the atmosphere from above
-   **Limb sounding**, representative of instruments observing along the
    atmospheric limb
-   **Zenith sounding**, representative of upward-looking observations
    from the ground

### Nadir sounder example

To run the nadir sounder example:

``` bash
cd nadir
./run.sh
```

### Limb sounder example

To run the limb sounder example:

``` bash
cd ../limb
./run.sh
```

### Zenith sounder example

To run the zenith sounder example:

``` bash
cd ../zenith
./run.sh
```

Each example is executed via a simple wrapper script that prepares the
input data, runs the radiative transfer calculations, and performs basic
verification steps.

------------------------------------------------------------------------

## Input files

In these examples, the following files are generated or used:

-   **Observation geometry**\
    The file `obs.tab` defines the viewing geometry of the instrument,
    including sensor position and line-of-sight information.

    ``` bash
    cat obs.tab
    ```

-   **Atmospheric state**\
    The file `atm.tab` contains a standard mid-latitude atmospheric
    profile, including pressure, temperature, and trace gas
    concentrations.

    ``` bash
    cat atm.tab
    ```

-   **Radiative transfer output**\
    The file `rad.tab` contains the simulated radiances written by the
    forward model for the configured detector channels.

    ``` bash
    cat rad.tab
    ```

These files illustrate the basic structure of JURASSIC input data and
are discussed in more detail in the User Manual.

------------------------------------------------------------------------

## Output and verification

After the simulation completes, JURASSIC:

-   Computes radiances for two or three detector channels
-   Compares the results against reference data to verify correctness
-   Generates diagnostic plots of radiances and Jacobians using gnuplot

If the run completes without errors and the verification checks pass,
your JURASSIC installation is functioning correctly.

------------------------------------------------------------------------

## Next steps

After completing the Quickstart, you may want to explore:

-   Customizing atmospheric profiles and observation geometries
-   Defining your own spectral bands and instrument configurations
-   Running JURASSIC on parallel HPC systems

These topics are covered in detail in the subsequent sections of the
User Manual.

------------------------------------------------------------------------

## Further information

More detailed documentation for users and developers is available in the
main [JURASSIC documentation](https://slcs-jsc.github.io/jurassic) and
in the [GitHub wiki](https://github.com/slcs-jsc/jurassic/wiki).

A detailed description of the JURASSIC model and its applications can be
found in the publications listed in the [References](references.md).

If you are interested in using JURASSIC for operational or research
applications, please do not hesitate to contact us for support.
