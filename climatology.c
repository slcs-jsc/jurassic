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
  Prepare atmospheric data file from climatological data.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm;
  static ctl_t ctl;

  double dz, t0, z, z0, z1;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  z0 = scan_ctl(argc, argv, "Z0", -1, "0", NULL);
  z1 = scan_ctl(argc, argv, "Z1", -1, "90", NULL);
  dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);

  /* Set atmospheric grid... */
  for (z = z0; z <= z1; z += dz) {
    atm.time[atm.np] = t0;
    atm.z[atm.np] = z;
    if ((++atm.np) >= NP)
      ERRMSG("Too many atmospheric grid points!");
  }

  /* Interpolate climatological data... */
  climatology(&ctl, &atm);

  /* Write data to disk... */
  write_atm(NULL, argv[2], &ctl, &atm);

  return EXIT_SUCCESS;
}
