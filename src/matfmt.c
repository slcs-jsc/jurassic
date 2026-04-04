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

  ctl_t ctl;
  int prof_in, prof_out;

  static atm_t atm;
  static obs_t obs;

  gsl_matrix *matrix;

  size_t nr, nc;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 11)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: matfmt <ctl> <atm> <obs> <rowspace> <colspace> <sort> <matrix_in> <matrixfmt_in> <matrix_out> <matrixfmt_out> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read dataset indices... */
  prof_in = (int) scan_ctl(argc, argv, "PROF_IN", -1, "0", NULL);
  prof_out = (int) scan_ctl(argc, argv, "PROF_OUT", -1, "0", NULL);

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
  read_matrix(NULL, argv[7], &ctl, matrix, prof_in);

  /* Write matrix data... */
  ctl.write_matrix = 1;
  ctl.matrixfmt = atoi(argv[10]);
  write_matrix(NULL, argv[9], &ctl, matrix, &atm, &obs, argv[4], argv[5],
	       argv[6], prof_out);

  /* Free... */
  gsl_matrix_free(matrix);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC matrix-format converter.\n\n");
  printf("Convert matrix files between supported MATRIXFMT formats and\n");
  printf("annotate them using atmospheric or observation metadata.\n\n");
  printf("Usage:\n");
  printf
    ("  matfmt <ctl> <atm> <obs> <rowspace> <colspace> <sort> <matrix_in> <matrixfmt_in> <matrix_out> <matrixfmt_out> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>           Control file.\n");
  printf("  <atm>           Atmospheric data file.\n");
  printf("  <obs>           Observation data file.\n");
  printf("  <rowspace>      Row space: x or y.\n");
  printf("  <colspace>      Column space: x or y.\n");
  printf("  <sort>          Sort order: r or c.\n");
  printf("  <matrix_in>     Input matrix file.\n");
  printf("  <matrixfmt_in>  Input matrix format identifier.\n");
  printf("  <matrix_out>    Output matrix file.\n");
  printf("  <matrixfmt_out> Output matrix format identifier.\n");
  printf("  [KEY VALUE]     Optional control parameters.\n\n");
  printf("Optional control overrides:\n");
  printf("  PROF_IN <n>     Read profile index <n> from the input file.\n");
  printf("  PROF_OUT <n>    Write output as profile index <n>.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
