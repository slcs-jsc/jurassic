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
  
  Copyright (C) 2003-2026 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert brightness temperature to radiance.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
static void usage(void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc != 3 && argc != 7)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: planck <t> <nu>\n"
	   "   or: planck <t0> <t1> <dt> <nu0> <nu1> <dnu>\n\n"
	   "Use -h for full help.");

  /* Calculate single value... */
  if (argc == 3) {

    /* Read arguments... */
    const double t = atof(argv[1]);
    const double nu = atof(argv[2]);

    /* Compute Planck function... */
    printf("%.10g\n", PLANCK(t, nu));
  }

  /* Calculate table... */
  else if (argc == 7) {

    /* Read arguments... */
    const double t0 = atof(argv[1]);
    const double t1 = atof(argv[2]);
    const double dt = atof(argv[3]);
    const double nu0 = atof(argv[4]);
    const double nu1 = atof(argv[5]);
    const double dnu = atof(argv[6]);

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

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC Planck-function converter.\n\n");
  printf("Convert brightness temperature to radiance for a single value\n");
  printf("or a tabulated temperature and wavenumber range.\n\n");
  printf("Usage:\n");
  printf("  planck <t> <nu>\n");
  printf("  planck <t0> <t1> <dt> <nu0> <nu1> <dnu>\n\n");
  printf("Arguments:\n");
  printf("  <t>    Brightness temperature [K].\n");
  printf("  <nu>   Wavenumber [cm^-1].\n");
  printf("  <t0>   First temperature value for table output.\n");
  printf("  <t1>   Last temperature value for table output.\n");
  printf("  <dt>   Temperature increment for table output.\n");
  printf("  <nu0>  First wavenumber for table output.\n");
  printf("  <nu1>  Last wavenumber for table output.\n");
  printf("  <dnu>  Wavenumber increment for table output.\n\n");
  printf("Output:\n");
  printf("  Writes results to standard output.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
