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
  Create observation geometry for a limb sounder.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static obs_t obs;

  double dt, dz, obsz, t, t0, t1, z, z0, z1;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <obs>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  obsz = scan_ctl(argc, argv, "OBSZ", -1, "780", NULL);
  t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  z0 = scan_ctl(argc, argv, "Z0", -1, "3", NULL);
  z1 = scan_ctl(argc, argv, "Z1", -1, "68", NULL);
  dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);

  /* Create measurement geometry... */
  for (t = t0; t <= t1; t += dt)
    for (z = z0; z <= z1; z += dz) {
      obs.time[obs.nr] = t;
      obs.obsz[obs.nr] = obsz;
      obs.vpz[obs.nr] = z;
      obs.vplat[obs.nr] = 180 / M_PI * acos((RE + z) / (RE + obsz));
      if ((++obs.nr) >= NR)
	ERRMSG("Too many rays!");
    }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs);

  return EXIT_SUCCESS;
}
