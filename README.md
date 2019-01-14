# Juelich Rapid Spectral Simulation Code

The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast
radiative transfer model for the mid-infrared spectral region.
JURASSIC has been developed for the analysis of atmospheric
remote sensing measurements.

JURASSIC uses the emissivity growth approximation to conduct the
radiative transfer calculations. Band transmittances are obtained
from pre-calculated look-up tables from line-by-line calculations.

The model was carefully tested in intercomparisons with
the Karlsruhe Optimized and Precise Radiative Transfer Algorithm
(KOPRA), the Reference Forward Model (RFM), and the
Stand-alone AIRS Radiative Transfer Algorithm (SARTA).

Further information can be found at:
http://www.fz-juelich.de/ias/jsc/jurassic

## Installation

This documentation describes the installation of JURASSIC on a Linux system.
A number of standard tools such as the GNU Compiler Collection (gcc)
and 'make' are required to compile JURASSIC.

Start by downloading the source code from the github repository:

    git clone https://github.com/slcs-jsc/jurassic

Change to the directory jurassic/ which holds source codes,
libraries, documentation, etc:

    cd jurassic

The GNU Scientific Library
(https://www.gnu.org/software/gsl) is required for numerical operations.
A copy of the GSL can be found in the repository, if it is not available
on your system. A script is provided to build the GSL:

    cd lib
    ./build.sh

Next, change to the source directory and edit the Makefile according to
your needs. In particular, check the paths to the libraries
(INCDIR and LIBDIR). Then try to compile the code:

    cd ../src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied and run
on other machines. Sometimes this causes problems. In this case remove
the '-static' flag from the CFLAGS in the Makefile and compile again.

By default the code will be compiled with OpenMP parallelization.
If this is not desired remove '-fopenmp' from the CFLAGS.

By default we use rather strict compiler warnings.
All warning messages will be turned into errors and no binaries will be
produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the src/ directory.

## Getting started

This script illustrates how to use JURASSIC:

    cd ../project/example
    ./run.sh

In this example we generate an observation geometry for an infrared nadir
sounder,

    cat obs.tab

a standard atmosphere for mid-latitudes,

    cat atm.tab

and conduct radiative transfer calculations for two detector channels
in the 15 micron CO2 waveband:

    cat rad.tab

The output of the simulation is verified by comparing it to reference data.

More details on the control parameters, data structures, and algorithms
can be found in the JURASSIC reference manual:

     evince ../doc/refman.pdf

## Contact

We are interested in sharing JURASSIC for research applications.

Please do not hesitate to contact us if you have any further questions:

Dr. Lars Hoffmann  
Forschungszentrum Jülich  
Jülich Supercomputing Centre  
52425 Jülich  
Germany  

e-mail: l.hoffmann@fz-juelich.de

## License

JURASSIC is distributed under the GNU GPL v3.
