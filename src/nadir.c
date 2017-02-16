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
  Create observation geometry for a nadir sounder.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;
  static obs_t obs;

  double dlat, lat, lat0, lat1, obsz;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <obs>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  obsz = scan_ctl(argc, argv, "OBSZ", -1, "700", NULL);
  lat0 = scan_ctl(argc, argv, "LAT0", -1, "-8.01", NULL);
  lat1 = scan_ctl(argc, argv, "LAT1", -1, "8.01", NULL);
  dlat = scan_ctl(argc, argv, "DLAT", -1, "0.18", NULL);

  /* Create measurement geometry... */
  for (lat = lat0; lat <= lat1; lat += dlat) {
    obs.obsz[obs.nr] = obsz;
    obs.vplat[obs.nr] = lat;
    if ((++obs.nr) >= NR)
      ERRMSG("Too many rays!");
  }

  /* Write observation data... */
  write_obs(NULL, argv[2], &ctl, &obs);

  return EXIT_SUCCESS;
}
