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
  Create observation geometry for zenith observations.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static obs_t obs;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <obs>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  const double t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  const double t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  const double dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  const double obsz = scan_ctl(argc, argv, "OBSZ", -1, "0", NULL);
  const double vpz = scan_ctl(argc, argv, "VPZ", -1, "700", NULL);
  const double theta0 = scan_ctl(argc, argv, "THETA0", -1, "0.0", NULL);
  const double theta1 = scan_ctl(argc, argv, "THETA1", -1, "0.0", NULL);
  const double dtheta = scan_ctl(argc, argv, "DTHETA", -1, "1.0", NULL);

  /* Create measurement geometry... */
  for (double t = t0; t <= t1; t += dt)
    for (double theta = theta0; theta <= theta1; theta += dtheta) {
      obs.time[obs.nr] = t;
      obs.obsz[obs.nr] = obsz;
      obs.vpz[obs.nr] = vpz;
      obs.vplat[obs.nr] =
	theta - RAD2DEG(asin((RE + obsz) / (RE + vpz) * sin(DEG2RAD(theta))));
      if ((++obs.nr) >= NR)
	ERRMSG("Too many rays!");
    }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs);

  return EXIT_SUCCESS;
}
