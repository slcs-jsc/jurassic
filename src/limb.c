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
  Create observation geometry for a limb sounder.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
static void usage(void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static obs_t obs;

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: limb <ctl> <obs> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  const double obsz = scan_ctl(argc, argv, "OBSZ", -1, "780", NULL);
  const double t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  const double t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  const double dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  const double z0 = scan_ctl(argc, argv, "Z0", -1, "3", NULL);
  const double z1 = scan_ctl(argc, argv, "Z1", -1, "68", NULL);
  const double dz = scan_ctl(argc, argv, "DZ", -1, "1", NULL);

  /* Create measurement geometry... */
  for (double t = t0; t <= t1; t += dt)
    for (double z = z0; z <= z1; z += dz) {
      obs.time[obs.nr] = t;
      obs.obsz[obs.nr] = obsz;
      obs.vpz[obs.nr] = z;
      obs.vplat[obs.nr] = 180 / M_PI * acos((RE + z) / (RE + obsz));
      if ((++obs.nr) >= NR)
	ERRMSG("Too many rays!");
    }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs, 0);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC limb-geometry tool.\n\n");
  printf("Create observation geometry for a limb sounding instrument.\n\n");
  printf("Usage:\n");
  printf("  limb <ctl> <obs> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>  Control file.\n");
  printf("  <obs>  Output observation geometry file.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Common control overrides:\n");
  printf("  OBSZ         Observer altitude [km].\n");
  printf("  T0, T1, DT   Time range and spacing.\n");
  printf("  Z0, Z1, DZ   Tangent-altitude range and spacing.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
