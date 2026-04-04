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
  Convert look-up table file format.
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

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: tblfmt <ctl> <tblbase_in> <tblfmt_in> <tblbase_out> <tblfmt_out> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read tables... */
  sprintf(ctl.tblbase, "%s", argv[2]);
  ctl.tblfmt = atoi(argv[3]);
  tbl_t *tbl = read_tbl(&ctl);

  /* Write tables... */
  sprintf(ctl.tblbase, "%s", argv[4]);
  ctl.tblfmt = atoi(argv[5]);
  write_tbl(&ctl, tbl);

  /* Free... */
  tbl_free(&ctl, tbl);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC lookup-table format converter.\n\n");
  printf("Convert look-up tables between supported TBLFMT formats.\n\n");
  printf("Usage:\n");
  printf
    ("  tblfmt <ctl> <tblbase_in> <tblfmt_in> <tblbase_out> <tblfmt_out> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>         Control file.\n");
  printf("  <tblbase_in>  Input table base name.\n");
  printf("  <tblfmt_in>   Input table format identifier.\n");
  printf("  <tblbase_out> Output table base name.\n");
  printf("  <tblfmt_out>  Output table format identifier.\n");
  printf("  [KEY VALUE]   Optional control parameters.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
