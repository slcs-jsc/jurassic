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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Convert look-up table file format.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  tbl_t *tbl;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <tblbase_in> <tblfmt_in>"
	   " <tblbase_out> <tblfmt_out>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Allocate... */
  ALLOC(tbl, tbl_t, 1);

  /* Read tables... */
  sprintf(ctl.tblbase, "%s", argv[2]);
  ctl.tblfmt = atoi(argv[3]);
  read_tbl(&ctl, tbl);

  /* Write tables... */
  sprintf(ctl.tblbase, "%s", argv[4]);
  ctl.tblfmt = atoi(argv[5]);
  write_tbl(&ctl, tbl);

  /* Free... */
  free(tbl);

  return EXIT_SUCCESS;
}
