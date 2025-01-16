/*
  This file is part of JURASSIC.
  
  JURASSIC is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  JURASSIC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with JURASSIC. If not, see <http://www.gnu.org/licenses/>.
  
  Copyright (C) 2003-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert brightness temperature to radiance.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  /* Check arguments... */
  if (argc != 3 && argc != 7)
    ERRMSG
      ("Give parameters: [ <t> <nu> |  <t0> <t1> <dt> <nu0> <nu1> <dnu> ]");

  /* Calculate single value... */
  if (argc == 3) {

    /* Read arguments... */
    double t = atof(argv[1]);
    double nu = atof(argv[2]);

    /* Compute Planck function... */
    printf("%.10g\n", PLANCK(t, nu));
  }

  /* Calculate table... */
  else if (argc == 7) {

    /* Read arguments... */
    double t0 = atof(argv[1]);
    double t1 = atof(argv[2]);
    double dt = atof(argv[3]);
    double nu0 = atof(argv[4]);
    double nu1 = atof(argv[5]);
    double dnu = atof(argv[6]);

    /* Write header... */
    printf("# $1 = brightness temperature [K]\n"
	   "# $2 = wavenumber [cm^-1]\n"
	   "# $3 = radiance [W/(m^2 sr cm^-1)]\n");

    /* Compute Planck function... */
    for (double t = t0; t <= t1; t += dt) {
      printf("\n");
      for (double nu = nu0; nu <= nu1; nu += dnu)
	printf("%.10g %.4f %.10g\n", t, nu, PLANCK(t, nu));
    }
  }

  return EXIT_SUCCESS;
}
