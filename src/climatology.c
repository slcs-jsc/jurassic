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
  
  Copyright (C) 2003-2026 Forschungszentrum Juelich GmbH
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

  double clk[NCL], sfeps[NSF];

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <atm>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  const double t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  const double t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  const double dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  const double z0 = scan_ctl(argc, argv, "Z0", -1, "0", NULL);
  const double z1 = scan_ctl(argc, argv, "Z1", -1, "90", NULL);
  const double dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);
  const int zsurf = (int) scan_ctl(argc, argv, "ZSURF", -1, "0", NULL);
  const double clz = scan_ctl(argc, argv, "CLZ", -1, "0", NULL);
  const double cldz = scan_ctl(argc, argv, "CLDZ", -1, "0", NULL);
  for (int icl = 0; icl < ctl.ncl; icl++)
    clk[icl] = scan_ctl(argc, argv, "CLK", icl, "0", NULL);
  const double sft = scan_ctl(argc, argv, "SFT", -1, "0", NULL);
  for (int isf = 0; isf < ctl.nsf; isf++)
    sfeps[isf] = scan_ctl(argc, argv, "SFEPS", isf, "1", NULL);

  /* Loop over time steps... */
  for (double t = t0; t <= t1 + 0.5 * dt; t += dt) {

    /* Add near surface layer... */
    if (zsurf) {
      atm.np = 6;
      atm.z[0] = 0;
      atm.z[1] = 0.01;
      atm.z[2] = 0.02;
      atm.z[3] = 0.05;
      atm.z[4] = 0.1;
      atm.z[5] = 0.2;
    } else
      atm.np = 0;

    /* Add heights... */
    for (double z = z0; z <= z1; z += dz)
      if (atm.np == 0 || z > atm.z[atm.np - 1]) {
	atm.z[atm.np] = z;
	if ((++atm.np) >= NP)
	  ERRMSG("Too many atmospheric grid points!");
      }

    /* Set time... */
    for (int ip = 0; ip < atm.np; ip++)
      atm.time[ip] = t;
  }

  /* Interpolate climatological data... */
  climatology(&ctl, &atm);

  /* Set cloud layer... */
  atm.clz = clz;
  atm.cldz = cldz;
  for (int icl = 0; icl < ctl.ncl; icl++)
    atm.clk[icl] = clk[icl];

  /* Set surface layer... */
  atm.sft = sft;
  for (int isf = 0; isf < ctl.nsf; isf++)
    atm.sfeps[isf] = sfeps[isf];

  /* Write data to disk... */
  write_atm(NULL, argv[2], &ctl, &atm);

  return EXIT_SUCCESS;
}
