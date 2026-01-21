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
  Calculate kernel functions.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Perform kernel calculations in a single directory. */
void call_kernel(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  char dirlist[LEN];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <kernel>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Initialize look-up tables... */
  tbl_t *tbl = read_tbl(&ctl);

  /* Get dirlist... */
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);

  /* Set flags... */
  ctl.write_matrix = 1;

  /* Single kernel calculation... */
  if (dirlist[0] == '-')
    call_kernel(&ctl, tbl, NULL, argv[2], argv[3], argv[4]);

  /* Work on directory list... */
  else {

    /* Open directory list... */
    FILE *in;
    if (!(in = fopen(dirlist, "r")))
      ERRMSG("Cannot open directory list!");

    /* Loop over directories... */
    char wrkdir[LEN];
    while (fscanf(in, "%4999s", wrkdir) != EOF) {

      /* Write info... */
      LOG(1, "\nWorking directory: %s", wrkdir);

      /* Call forward model... */
      call_kernel(&ctl, tbl, wrkdir, argv[2], argv[3], argv[4]);
    }

    /* Close dirlist... */
    fclose(in);
  }

  /* Free... */
  tbl_free(&ctl, tbl);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void call_kernel(
  const ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile) {

  static atm_t atm;
  static obs_t obs;

  /* Read observation geometry... */
  read_obs(wrkdir, obsfile, ctl, &obs);

  /* Read atmospheric data... */
  read_atm(wrkdir, atmfile, ctl, &atm);

  /* Get sizes... */
  const size_t n = atm2x(ctl, &atm, NULL, NULL, NULL);
  const size_t m = obs2y(ctl, &obs, NULL, NULL, NULL);

  /* Check sizes... */
  if (n == 0)
    ERRMSG("No state vector elements!");
  if (m == 0)
    ERRMSG("No measurement vector elements!");

  /* Allocate... */
  gsl_matrix *k = gsl_matrix_alloc(m, n);

  /* Compute kernel matrix... */
  kernel(ctl, tbl, &atm, &obs, k);

  /* Write matrix to file... */
  write_matrix(wrkdir, kernelfile, ctl, k, &atm, &obs, "y", "x", "r");

  /* Free... */
  gsl_matrix_free(k);
}
