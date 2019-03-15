# Juelich Rapid Spectral Simulation Code

The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast infrared radiative transfer model for the analysis of atmospheric remote sensing measurements.

![GitHub tag (latest SemVer)](https://img.shields.io/github/tag/slcs-jsc/jurassic.svg)
![GitHub top language](https://img.shields.io/github/languages/top/slcs-jsc/jurassic.svg)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/slcs-jsc/jurassic.svg)
![GitHub last commit](https://img.shields.io/github/last-commit/slcs-jsc/jurassic.svg)
![GitHub](https://img.shields.io/github/license/slcs-jsc/jurassic.svg)

## Features

* JURASSIC uses the emissivity growth approximation to conduct the radiative transfer calculations.
* Band transmittances are obtained from pre-calculated look-up tables from line-by-line calculations.
* The model was carefully tested in intercomparisons with the Karlsruhe Optimized and Precise Radiative Transfer Algorithm (KOPRA), the Reference Forward Model (RFM), and the Stand-alone AIRS Radiative Transfer Algorithm (SARTA).
* JURASSIC features an MPI/OpenMP hybrid parallelization for efficient use on supercomputers.

## Getting started

### Prerequisites

This documentation describes the installation of JURASSIC on a Linux system. A number of standard tools (gcc, git, make) and software libraries are needed to install JURASSIC. The [GNU Scientific Library](https://www.gnu.org/software/gsl) is required for numerical calculations. A copy of this library can be found in the git repository.

Start by downloading the source code from the github repository:

    git clone https://github.com/slcs-jsc/jurassic

To update an existing installation use:

    git pull https://github.com/slcs-jsc/jurassic

### Installation

First, compile the GSL library needed for JURASSIC by using the build script:

    cd jurassic/lib
    ./build.sh

Next, change to the source directory, edit the Makefile according to your needs, and try to compile the code:

    cd jurassic/src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied and run on other machines. Sometimes this causes problems. In this case remove the '-static' flag from the CFLAGS in the Makefile and compile again.

By default we use rather strict compiler warnings. All warning messages will be turned into errors and no binaries will be produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the jurassic/src/ directory.

### Try the example

It is recommended that you create a project directory for testing the example and also to store other experiments:

    mkdir -p jurassic/projects

This shows how to run the example:

    cp -a jurassic/example jurassic/projects
    cd jurassic/projects/example
    ./run.sh

In this example we generate an observation geometry for an infrared nadir sounder,

    cat obs.tab

a standard atmosphere for mid-latitudes,

    cat atm.tab

and conduct radiative transfer calculations for two detector channels in the 15 micron CO2 waveband:

    cat rad.tab

The output of the simulation is verified by comparing it to reference data.

## Further information

These are the main references for citing the JURASSIC model in scientific publications:

* Hoffmann, L., and M. J. Alexander, Retrieval of stratospheric temperatures from Atmospheric Infrared Sounder radiance measurements for gravity wave studies, J. Geophys. Res., 114, D07105, https://doi.org/10.1029/2008JD011241, 2009.

* Hoffmann, L., Kaufmann, M., Spang, R., Müller, R., Remedios, J. J., Moore, D. P., Volk, C. M., von Clarmann, T., and Riese, M.: Envisat MIPAS measurements of CFC-11: retrieval, validation, and climatology, Atmos. Chem. Phys., 8, 3671-3688, https://doi.org/10.5194/acp-8-3671-2008, 2008.

More details on the data structures and algorithms can be found in the [JURASSIC reference manual](doc/refman.pdf)

## Contributing

We are interested in sharing JURASSIC for operational or research applications.

Please do not hesitate to contact us, if you have any further questions or need support.

## License

JURASSIC is distributed under the GNU GPL v3.

## Contact

Dr. Lars Hoffmann  

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: l.hoffmann@fz-juelich.de

website: https://www.fz-juelich.de/ias/jsc/slcs
