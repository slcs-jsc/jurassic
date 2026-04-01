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
  Convert atmospheric data files.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  static atm_t atm;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <atm_in> <atmfmt_in>"
	   " <atm_out> <atmfmt_out> [PROF_IN n] [PROF_OUT n]");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read profile indices... */
  int prof_in = (int) scan_ctl(argc, argv, "PROF_IN", -1, "0", NULL);
  int prof_out = (int) scan_ctl(argc, argv, "PROF_OUT", -1, "0", NULL);

  /* Read atmospheric data... */
  ctl.atmfmt = atoi(argv[3]);
  if (ctl.atmfmt == 3)
    read_atm_nc(argv[2], &ctl, &atm, prof_in);
  else
    read_atm(NULL, argv[2], &ctl, &atm);

  /* Write atmospheric data... */
  ctl.atmfmt = atoi(argv[5]);
  if (ctl.atmfmt == 3)
    write_atm_nc(argv[4], &ctl, &atm, prof_out);
  else
    write_atm(NULL, argv[4], &ctl, &atm);

  return EXIT_SUCCESS;
}
