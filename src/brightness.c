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
  Convert radiance to brightness temperature.
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
    ERRMSG
      ("Missing or invalid command-line arguments.\n\n"
       "Usage: brightness <rad> <nu>\n"
       "   or: brightness <rad0> <rad1> <drad> <nu0> <nu1> <dnu>\n\n"
       "Use -h for full help.");

  /* Calculate single value... */
  if (argc == 3) {

    /* Read arguments... */
    double rad = atof(argv[1]);
    double nu = atof(argv[2]);

    /* Compute brightness temperature... */
    printf("%.10g\n", BRIGHT(rad, nu));

  }

  /* Calculate table... */
  else if (argc == 7) {

    /* Read arguments... */
    double rad0 = atof(argv[1]);
    double rad1 = atof(argv[2]);
    double drad = atof(argv[3]);
    double nu0 = atof(argv[4]);
    double nu1 = atof(argv[5]);
    double dnu = atof(argv[6]);

    /* Write header... */
    printf("# $1 = radiance [W/(m^2 sr cm^-1)]\n"
	   "# $2 = wavenumber [cm^-1]\n"
	   "# $3 = brightness temperature [K]\n");

    /* Compute brightness temperature... */
    for (double rad = rad0; rad <= rad1; rad += drad) {
      printf("\n");
      for (double nu = nu0; nu <= nu1; nu += dnu)
	printf("%.10g %.4f %.10g\n", rad, nu, BRIGHT(rad, nu));
    }
  }

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC brightness-temperature converter.\n\n");
  printf("Convert radiance to brightness temperature for a single value\n");
  printf("or a tabulated radiance and wavenumber range.\n\n");
  printf("Usage:\n");
  printf("  brightness <rad> <nu>\n");
  printf("  brightness <rad0> <rad1> <drad> <nu0> <nu1> <dnu>\n\n");
  printf("Arguments:\n");
  printf("  <rad>   Radiance [W/(m^2 sr cm^-1)].\n");
  printf("  <nu>    Wavenumber [cm^-1].\n");
  printf("  <rad0>  First radiance value for table output.\n");
  printf("  <rad1>  Last radiance value for table output.\n");
  printf("  <drad>  Radiance increment for table output.\n");
  printf("  <nu0>   First wavenumber for table output.\n");
  printf("  <nu1>   Last wavenumber for table output.\n");
  printf("  <dnu>   Wavenumber increment for table output.\n\n");
  printf("Output:\n");
  printf("  Writes results to standard output.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
