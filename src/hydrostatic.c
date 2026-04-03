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
  Recalculate pressure based on hydrostatic equilibrium.
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

  static atm_t atm;
  static ctl_t ctl;

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: hydrostatic <ctl> <atm_in> <atm_hyd> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Check reference height... */
  if (ctl.hydz < 0)
    ERRMSG("Set HYDZ>=0!");

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm, 0);

  /* Build atmosphere based on hydrostatic equilibrium... */
  hydrostatic(&ctl, &atm);

  /* Write atmospheric data... */
  write_atm(NULL, argv[3], &ctl, &atm, 0);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC hydrostatic tool.\n\n");
  printf("Recalculate pressure from hydrostatic equilibrium for an\n");
  printf("existing atmospheric profile.\n\n");
  printf("Usage:\n");
  printf("  hydrostatic <ctl> <atm_in> <atm_hyd> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>      Control file.\n");
  printf("  <atm_in>   Input atmospheric data file.\n");
  printf("  <atm_hyd>  Output atmospheric data file with hydrostatic pressure.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Notes:\n");
  printf("  HYDZ must be set to a non-negative reference height.\n");
  printf("\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
