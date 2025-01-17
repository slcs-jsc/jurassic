# Quickstart

## Example of a JURASSIC simulation

JURASSIC provides a project directory for testing the examples and
also to store other experiments:

    cd jurassic/projects

This shows how to run the example for the nadir sounder:

    cd nadir ./run.sh

This shows how to run the example for the limb sounder:

    cd ../limb ./run.sh

In both examples, we generate an observation geometry file,

    cat obs.tab

a standard atmosphere for mid-latitudes,

    cat atm.tab

and conduct radiative transfer calculations for two or three detector
channels:

    cat rad.tab

The output of the simulation is verified by comparing it to reference
data.  Additionally, gnuplot is used to create plots of the radiance
data:

<p align="center">
  <img src="projects/limb/plot_rad.png" alt="limb radiance data" width="45%"/>
  &emsp;
  <img src="projects/nadir/plot_rad.png" alt="nadir radiance data" width="45%"/>
</p>

Kernel functions are calculated using a finite difference method:

<p align="center">
  <img src="projects/limb/plot_kernel_temperature_792.png" alt="limb temperature kernel function" width="45%"/>
  &emsp;
  <img src="projects/nadir/plot_kernel_temperature_668.5410.png" alt="nadir temperature kernel function" width="45%"/>
</p>

<p align="center">
  <img src="projects/limb/plot_kernel_H2O_792.png" alt="limb water vapor kernel function" width="45%"/>
  &emsp;
  <img src="projects/nadir/plot_kernel_CO2_668.5410.png" alt="nadir water vapor kernel function" width="45%"/>
</p>

## Further information

More detailed information for new users and developers is provided in
the [JURASSIC user manual](https://slcs-jsc.github.io/jurassic) and
collected in the [GitHub wiki](https://github.com/slcs-jsc/jurassic/wiki).

A detailed description of the JURASSIC model is provided in this paper:

* Baumeister, P. F. and Hoffmann, L.: Fast infrared radiative transfer
  calculations using graphics processing units: JURASSIC-GPU v2.0,
  Geosci. Model Dev., 15, 1855–1874,
  <https://doi.org/10.5194/gmd-15-1855-2022>, 2022.

* Hoffmann, L., and M. J. Alexander, Retrieval of stratospheric
  temperatures from Atmospheric Infrared Sounder radiance measurements
  for gravity wave studies, J. Geophys. Res., 114, D07105,
  <https://doi.org/10.1029/2008JD011241>, 2009.

* Hoffmann, L., Kaufmann, M., Spang, R., Müller, R., Remedios, J. J.,
  Moore, D. P., Volk, C. M., von Clarmann, T., and Riese, M.: Envisat
  MIPAS measurements of CFC-11: retrieval, validation, and
  climatology, Atmos. Chem. Phys., 8, 3671-3688,
  <https://doi.org/10.5194/acp-8-3671-2008>, 2008.

We are interested in sharing JURASSIC for operational and research
applications. Please do not hesitate to contact us, if you have any
further questions or need support.
