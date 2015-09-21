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
  
  Copright (C) 2003-2015 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Interpolate atmospheric data to another spatial grid.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm_in, atm_pts;
  static ctl_t ctl;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_pts> <atm_out>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm_in);
  read_atm(NULL, argv[3], &ctl, &atm_pts);

  /* Interpolate atmospheric data... */
  intpol_atm(&ctl, &atm_pts, &atm_in);

  /* Save interpolated data... */
  write_atm(NULL, argv[4], &ctl, &atm_pts);

  return EXIT_SUCCESS;
}
