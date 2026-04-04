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
  Interpolate atmospheric data to another spatial grid.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
static void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static atm_t atm_in, atm_pts;
  static ctl_t ctl;

  double k[NW], q[NG];

  /* Print usage information... */
  USAGE;

  /* Interpolate atmospheric data... */

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: interpolate <ctl> <atm_in> <atm_pts> <atm_out> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm_in, 0);
  read_atm(NULL, argv[3], &ctl, &atm_pts, 0);

  /* Interpolate atmospheric data... */
  for (int ip = 0; ip < atm_pts.np; ip++) {
    intpol_atm(&ctl, &atm_in, atm_pts.z[ip],
	       &atm_pts.p[ip], &atm_pts.t[ip], q, k);
    for (int ig = 0; ig < ctl.ng; ig++)
      atm_pts.q[ig][ip] = q[ig];
    for (int iw = 0; iw < ctl.nw; iw++)
      atm_pts.k[iw][ip] = k[iw];
  }

  /* Save interpolated data... */
  write_atm(NULL, argv[4], &ctl, &atm_pts, 0);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC interpolation tool.\n\n");
  printf("Interpolate atmospheric state variables from one atmospheric\n");
  printf("profile to the grid defined by another profile.\n\n");
  printf("Usage:\n");
  printf
    ("  interpolate <ctl> <atm_in> <atm_pts> <atm_out> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>      Control file.\n");
  printf("  <atm_in>   Input atmospheric data file.\n");
  printf("  <atm_pts>  Atmospheric file defining interpolation points.\n");
  printf("  <atm_out>  Output atmospheric data file.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
