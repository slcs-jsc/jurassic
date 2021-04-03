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
  
  Copyright (C) 2003-2021 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Calculate kernel functions.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm;
  static ctl_t ctl;
  static obs_t obs;

  gsl_matrix *k;

  size_t m, n;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <kernel>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Set flags... */
  ctl.write_matrix = 1;

  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs);

  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm);

  /* Get sizes... */
  n = atm2x(&ctl, &atm, NULL, NULL, NULL);
  m = obs2y(&ctl, &obs, NULL, NULL, NULL);

  /* Check sizes... */
  if (n <= 0)
    ERRMSG("No state vector elements!");
  if (m <= 0)
    ERRMSG("No measurement vector elements!");

  /* Allocate... */
  k = gsl_matrix_alloc(m, n);

  /* Compute kernel matrix... */
  kernel(&ctl, &atm, &obs, k);

  /* Write matrix to file... */
  write_matrix(NULL, argv[4], &ctl, k, &atm, &obs, "y", "x", "r");

  /* Free... */
  gsl_matrix_free(k);

  return EXIT_SUCCESS;
}
