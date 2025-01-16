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
  
  Copyright (C) 2003-2025 Forschungszentrum Juelich GmbH
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
  double t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  double t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  double dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  double z0 = scan_ctl(argc, argv, "Z0", -1, "0", NULL);
  double z1 = scan_ctl(argc, argv, "Z1", -1, "90", NULL);
  double dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);
  double clz = scan_ctl(argc, argv, "CLZ", -1, "0", NULL);
  double cldz = scan_ctl(argc, argv, "CLDZ", -1, "0", NULL);
  for (int icl = 0; icl < ctl.ncl; icl++)
    clk[icl] = scan_ctl(argc, argv, "CLK", icl, "0", NULL);
  double sfz = scan_ctl(argc, argv, "SFZ", -1, "0", NULL);
  double sfp = scan_ctl(argc, argv, "SFP", -1, "0", NULL);
  double sft = scan_ctl(argc, argv, "SFT", -1, "0", NULL);
  for (int isf = 0; isf < ctl.nsf; isf++)
    sfeps[isf] = scan_ctl(argc, argv, "SFEPS", isf, "1", NULL);

  /* Set atmospheric grid... */
  for (double t = t0; t <= t1; t += dt)
    for (double z = z0; z <= z1; z += dz) {
      atm.time[atm.np] = t;
      atm.z[atm.np] = z;
      if ((++atm.np) >= NP)
	ERRMSG("Too many atmospheric grid points!");
    }

  /* Interpolate climatological data... */
  climatology(&ctl, &atm);

  /* Set cloud layer... */
  atm.clz = clz;
  atm.cldz = cldz;
  for (int icl = 0; icl < ctl.ncl; icl++)
    atm.clk[icl] = clk[icl];

  /* Set surface layer... */
  atm.sfz = sfz;
  atm.sfp = sfp;
  atm.sft = sft;
  for (int isf = 0; isf < ctl.nsf; isf++)
    atm.sfeps[isf] = sfeps[isf];

  /* Write data to disk... */
  write_atm(NULL, argv[2], &ctl, &atm);

  return EXIT_SUCCESS;
}
