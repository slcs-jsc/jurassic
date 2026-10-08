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
  Create observation geometry for a nadir sounder.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
static void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static obs_t obs;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: nadir <ctl> <obs> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  const double t0 = scan_ctl(argc, argv, "T0", -1, "0", NULL);
  const double t1 = scan_ctl(argc, argv, "T1", -1, "0", NULL);
  const double dt = scan_ctl(argc, argv, "DT", -1, "1", NULL);
  const double obsz = scan_ctl(argc, argv, "OBSZ", -1, "700", NULL);
  const double lat0 = scan_ctl(argc, argv, "LAT0", -1, "-8.01", NULL);
  const double lat1 = scan_ctl(argc, argv, "LAT1", -1, "8.01", NULL);
  const double dlat = scan_ctl(argc, argv, "DLAT", -1, "0.18", NULL);

  /* Create measurement geometry... */
  for (double t = t0; t <= t1; t += dt)
    for (double lat = lat0; lat <= lat1; lat += dlat) {
      obs.time[obs.nr] = t;
      obs.obsz[obs.nr] = obsz;
      obs.vplat[obs.nr] = lat;
      if ((++obs.nr) >= NR)
	ERRMSG("Too many rays!");
    }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs, 0);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC nadir-geometry tool.\n\n");
  printf("Create observation geometry for a nadir sounding instrument.\n\n");
  printf("Usage:\n");
  printf("  nadir <ctl> <obs> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>  Control file.\n");
  printf("  <obs>  Output observation geometry file.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Tool-specific control parameters:\n");
  printf("  OBSZ           Observer altitude [km].\n");
  printf("  T0, T1, DT     Time range and spacing.\n");
  printf("  LAT0, LAT1     Latitude range.\n");
  printf("  DLAT           Latitude spacing.\n\n");
  printf("Common control parameters:\n");
  printf("  OBSFMT                 Output observation file format.\n");
  printf
    ("  ND, NU[i]              Spectral channels in observation files.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
