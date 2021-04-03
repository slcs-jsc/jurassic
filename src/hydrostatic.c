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
  Recalculate pressure based on hydrostatic equilibrium.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm;
  static ctl_t ctl;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <atm_in> <atm_hyd>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Check reference height... */
  if (ctl.hydz < 0)
    ERRMSG("Set HYDZ>=0!");

  /* Read atmospheric data... */
  read_atm(NULL, argv[2], &ctl, &atm);

  /* Build atmosphere based on hydrostatic equilibrium... */
  hydrostatic(&ctl, &atm);

  /* Write atmospheric data... */
  write_atm(NULL, argv[3], &ctl, &atm);

  return EXIT_SUCCESS;
}
