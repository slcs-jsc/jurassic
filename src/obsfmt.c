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
  Convert observation data files.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  ctl_t ctl;

  static obs_t obs;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <obs_in> <obsfmt_in>"
	   " <obs_out> <obsfmt_out> [PROF_IN n] [PROF_OUT n]");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read profile indices... */
  const int prof_in = (int) scan_ctl(argc, argv, "PROF_IN", -1, "0", NULL);
  const int prof_out = (int) scan_ctl(argc, argv, "PROF_OUT", -1, "0", NULL);

  /* Read observation data... */
  ctl.obsfmt = atoi(argv[3]);
  read_obs(NULL, argv[2], &ctl, &obs, prof_in);

  /* Write observation data... */
  ctl.obsfmt = atoi(argv[5]);
  write_obs(NULL, argv[4], &ctl, &obs, prof_out);

  return EXIT_SUCCESS;
}
