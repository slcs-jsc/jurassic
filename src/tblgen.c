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
  Prepapre look-up tables from monochromatic absorption spectra.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Maximum number of grid points for filter files: */
#define MAXNF 20000

/* Maximum number of grid points for spectra: */
#define MAXNPTS 10000000

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

  FILE *in;

  static char *line = NULL;

  static double dnu, abs[MAXNPTS], epsold, f, filt[MAXNF],
    nu, nu0, nu1, nuf[MAXNF], press, temp, u;

  static int i, idx, nf, npts;

  size_t line_buf_size = 0;

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Read command line arguments... */
  if (argc != 5)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: tblgen <press> <temp> <spec> <filter>\n\n"
	   "Use -h for full help.");
  sscanf(argv[1], "%lg", &press);
  sscanf(argv[2], "%lg", &temp);

  /* Compute column density [molec/cm^2] (1 km path length, 1 ppmv)... */
  u = 1e-6 * press * 100 / (1.380658e-23 * temp) * 1000 / 1e4;

  /* Read filter function... */
  if (!(in = fopen(argv[4], "r")))
    ERRMSG("Cannot open filter file!");
  while (getline(&line, &line_buf_size, in) != -1)
    if (sscanf(line, "%lg %lg", &nuf[nf], &filt[nf]) == 2)
      if (++nf >= MAXNF)
	ERRMSG("Too many points in filter function");
  fclose(in);

  /* Read spectrum... */
  if (!(in = fopen(argv[3], "r")))
    ERRMSG("Cannot open spectrum!");
  for (i = 0; i < 4; ++i)
    if (getline(&line, &line_buf_size, in) == -1)
      ERRMSG("Error while reading spectrum!");
  sscanf(line, "%d %lg %lg %lg", &npts, &nu0, &dnu, &nu1);
  if (npts > MAXNPTS)
    ERRMSG("Too many points in optical depth spectrum!");
  i = 0;
  while (getline(&line, &line_buf_size, in) != -1) {
    char *tok = strtok(line, " \t\n");
    while (tok != NULL) {
      if (i >= MAXNPTS)
	ERRMSG("Too many data points in spectrum file!");
      sscanf(tok, "%lg", &abs[i]);
      abs[i] /= u;
      ++i;
      tok = strtok(NULL, " \t\n");
    }
  }
  fclose(in);

  /* Set grid spacing... */
  dnu = (nu1 - nu0) / ((double) npts - 1.0);
  const double r0 = (nuf[0] - nu0) / (nu1 - nu0) * (double) npts;
  const int i0 = (int) r0;

  /* Loop over column densities... */
  for (u = 1.0; u <= 1e30; u *= 1.122) {

    /* Integrate... */
    double epssum = 0, fsum = 0;
    for (i = i0; i < npts; i++) {
      nu = nu0 + dnu * (double) i;
      if (nu < nuf[0])
	continue;
      else if (nu > nuf[nf - 1])
	break;
      else {
	if (nu < nuf[idx] || nu > nuf[idx + 1])
	  idx = locate_irr(nuf, nf, nu);
	f = LIN(nuf[idx], filt[idx], nuf[idx + 1], filt[idx + 1], nu);
	fsum += f;
	epssum += f * exp(-abs[i] * u);
      }
    }
    epssum = 1 - epssum / fsum;

    /* Write output... */
    if (epssum >= 1e-6 && epssum <= 0.999999 && epssum > epsold)
      printf("%g %g %g %g\n", press, temp, u, epssum);
    epsold = epssum;

    /* Check for termination... */
    if (epssum > 0.999999)
      return EXIT_SUCCESS;
  }

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC lookup-table generator.\n\n");
  printf("Prepare transmission look-up table entries from monochromatic\n");
  printf("absorption spectra and a filter function.\n\n");
  printf("Usage:\n");
  printf("  tblgen <press> <temp> <spec> <filter>\n\n");
  printf("Arguments:\n");
  printf("  <press>   Pressure [hPa].\n");
  printf("  <temp>    Temperature [K].\n");
  printf("  <spec>    Monochromatic absorption spectrum file.\n");
  printf("  <filter>  Filter-function file.\n\n");
  printf("Output:\n");
  printf("  Writes results to standard output.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
