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
  Convert atmospheric data files.
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

  static atm_t atm;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: atmfmt <ctl> <atm_in> <atmfmt_in> <atm_out> <atmfmt_out> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read profile indices... */
  int prof_in = (int) scan_ctl(argc, argv, "PROF_IN", -1, "0", NULL);
  int prof_out = (int) scan_ctl(argc, argv, "PROF_OUT", -1, "0", NULL);

  /* Read atmospheric data... */
  ctl.atmfmt = atoi(argv[3]);
  read_atm(NULL, argv[2], &ctl, &atm, prof_in);

  /* Write atmospheric data... */
  ctl.atmfmt = atoi(argv[5]);
  write_atm(NULL, argv[4], &ctl, &atm, prof_out);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC atmospheric-format converter.\n\n");
  printf
    ("Convert atmospheric input files between supported ATMFMT formats.\n\n");
  printf("Usage:\n");
  printf
    ("  atmfmt <ctl> <atm_in> <atmfmt_in> <atm_out> <atmfmt_out> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>         Control file.\n");
  printf("  <atm_in>      Input atmospheric data file.\n");
  printf("  <atmfmt_in>   Input atmospheric file format identifier.\n");
  printf("  <atm_out>     Output atmospheric data file.\n");
  printf("  <atmfmt_out>  Output atmospheric file format identifier.\n");
  printf("  [KEY VALUE]   Optional control parameters.\n\n");
  printf("Optional control overrides:\n");
  printf("  PROF_IN <n>   Read profile index <n> from the input file.\n");
  printf("  PROF_OUT <n>  Write output as profile index <n>.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
