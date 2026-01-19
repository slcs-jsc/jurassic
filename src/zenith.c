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
  const double vpz = scan_ctl(argc, argv, "VPZ", -1, "700", NULL);
  const double theta0 = scan_ctl(argc, argv, "THETA0", -1, "0.0", NULL);
  const double theta1 = scan_ctl(argc, argv, "THETA1", -1, "0.0", NULL);
  const double dtheta = scan_ctl(argc, argv, "DTHETA", -1, "1.0", NULL);
  const double obsz = scan_ctl(argc, argv, "OBSZ", -1, "0", NULL);
  
  /* Set distances... */
  const double ro = RE + obsz;	/* observer radius */
  const double rv = RE + vpz;	/* shell radius */

  /* Create measurement geometry... */
  for (double t = t0; t <= t1; t += dt) {
    for (double theta = theta0; theta <= theta1 + 1e-12; theta += dtheta) {

      /* Check number of ray paths... */
      if ((obs.nr) >= NR)
	ERRMSG("Too many rays!");

      /* Auxilary variables... */
      const double th = DEG2RAD(theta);
      const double sth = sin(th);
      const double cth = cos(th);

      /*
         Ray in local x-z plane:
         A = (0, 0, ro)
         u = (sin(th), 0, cos(th))  (|u|=1)
         P(s) = A + s u, s >= 0

         Intersect with sphere |P| = rv:
         s^2 + 2 ro cos(th) s + (ro^2 - rv^2) = 0
         Discriminant:
         disc = rv^2 - ro^2 sin^2(th)
       */
      double disc = rv * rv - ro * ro * sth * sth;
      if (disc < 0.0) {
	if (disc > -1e-12 * rv * rv)
	  disc = 0.0;		// tolerate tiny negatives
	else
	  continue;
      }

      /* If disc < 0, ray misses the shell: skip this ray */
      if (disc < 0.0)
	continue;

      /* Two intersections along the ray... */
      const double sq = sqrt(disc);
      const double s1 = -ro * cth - sq;
      const double s2 = -ro * cth + sq;

      /* Choose nearest forward intersection (smallest s >= 0)... */
      double s = DBL_MAX;
      if (s1 >= 0.0)
	s = s1;
      if (s2 >= 0.0 && s2 < s)
	s = s2;

      /* If both are behind observer, skip... */
      if (s == DBL_MAX)
	continue;

      /* Intersection point B in local x-z plane... */
      const double bx = s * sth;
      const double bz = ro + s * cth;

      /*
         Central angle beta = angle between OA (z-axis) and OB.
         In this local frame, OA is +z. For points in x-z plane:
         beta = atan2(bx, bz)
         This preserves sign (north/south) consistent with theta sign.
       */
      const double beta = atan2(bx, bz);

      /* Set observer and view point... */
      obs.time[obs.nr] = t;
      obs.obsz[obs.nr] = obsz;
      obs.vpz[obs.nr] = vpz;

      /* Stored as "view point latitude-like" angle for this 2-D meridional setup... */
      obs.vplat[obs.nr] = RAD2DEG(beta);

      /* Increment ray path counter... */
      obs.nr++;
    }
  }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs);

  return EXIT_SUCCESS;
}
