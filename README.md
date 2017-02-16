# Juelich Rapid Spectral Simulation Code

The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast
radiative transfer model for the mid-infrared spectral region.
JURASSIC has been developed for the analysis of atmospheric
remote sensing measurements.

## Algorithms

JURASSIC uses the emissivity growth approximation (EGA) to conduct the
radiative transfer calculations. Band transmittances are obtained
from pre-calculated look-up tables from line-by-line calculations.

JURASSIC features a 3-D raytracer to realize different observation
geometries and to enable tomographic studies. JURASSIC has also been
extended to perform Mie-scattering calculations for aerosol and ice
particles.

The model was carefully tested in intercomparisons with
the Karlsruhe Optimized and Precise Radiative Transfer Algorithm
(KOPRA), the Reference Forward Model (RFM), and the
Stand-alone AIRS Radiative Transfer Algorithm (SARTA).

## Implementation

JURASSIC is written in the C programming language. It uses the
GNU Scientific Library (GSL; https://www.gnu.org/software/gsl/)
for numerical operations. JURASSIC can be used on a desktop PC,
but it also features an MPI/OpenMP hybrid parallelization for
efficient use on supercomputers.

## Contact

We are interested in sharing JURASSIC for research applications.

Please do not hesitate to contact us with further questions:

Dr. Lars Hoffmann  
Forschungszentrum Jülich  
Jülich Supercomputing Centre  
52425 Jülich  
Germany  

e-mail: l.hoffmann@fz-juelich.de

## License

JURASSIC is distributed under the GNU GPL v3.
