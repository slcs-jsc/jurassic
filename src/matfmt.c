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
  
  Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH
*/

/*!
  \file
  Convert matrix data files.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  static atm_t atm;
  static obs_t obs;

  gsl_matrix *matrix;

  size_t nr, nc;

  /* Check arguments... */
  if (argc < 11)
    ERRMSG("Give parameters: <ctl> <atm> <obs> <rowspace> <colspace> <sort>"
	   " <matrix_in> <matrixfmt_in> <matrix_out> <matrixfmt_out>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Check matrix metadata... */
  if (argv[4][0] != 'x' && argv[4][0] != 'y')
    ERRMSG("Unknown row space, use x or y!");
  if (argv[5][0] != 'x' && argv[5][0] != 'y')
    ERRMSG("Unknown column space, use x or y!");
  if (argv[6][0] != 'r' && argv[6][0] != 'c')
    ERRMSG("Unknown sort order, use r or c!");

  /* Read atmospheric data if needed... */
  if (argv[4][0] == 'x' || argv[5][0] == 'x')
    read_atm(NULL, argv[2], &ctl, &atm, 0);

  /* Read observation data if needed... */
  if (argv[4][0] == 'y' || argv[5][0] == 'y')
    read_obs(NULL, argv[3], &ctl, &obs, 0);

  /* Determine matrix size... */
  if (argv[4][0] == 'y')
    nr = obs2y(&ctl, &obs, NULL, NULL, NULL);
  else
    nr = atm2x(&ctl, &atm, NULL, NULL, NULL);

  if (argv[5][0] == 'y')
    nc = obs2y(&ctl, &obs, NULL, NULL, NULL);
  else
    nc = atm2x(&ctl, &atm, NULL, NULL, NULL);

  /* Allocate matrix... */
  matrix = gsl_matrix_alloc(nr, nc);

  /* Read matrix data... */
  ctl.matrixfmt = atoi(argv[8]);
  read_matrix(NULL, argv[7], &ctl, matrix);

  /* Write matrix data... */
  ctl.write_matrix = 1;
  ctl.matrixfmt = atoi(argv[10]);
  write_matrix(NULL, argv[9], &ctl, matrix, &atm, &obs, argv[4], argv[5],
	       argv[6]);

  /* Free... */
  gsl_matrix_free(matrix);

  return EXIT_SUCCESS;
}
