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
  Converter for spectra.
*/

#include <omp.h>
#include "jurassic.h"		// ctl_t, obs_t, atm_t, read_*, write_*, formod_*PU

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  obs_t *obs;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <obs> <spec.tab>");

  /* Allocate... */
  ALLOC(obs, obs_t, 1);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read observation geometry... */
  read_obs(".", argv[2], &ctl, obs);

  /* Write info... */
  printf("Write spectra: %s\n", argv[3]);

  /* Create file... */
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2 = observer altitude [km]\n"
	  "# $3 = observer longitude [deg]\n"
	  "# $4 = observer latitude [deg]\n"
	  "# $5 = view point altitude [km]\n"
	  "# $6 = view point longitude [deg]\n"
	  "# $7 = view point latitude [deg]\n"
	  "# $8 = tangent point altitude [km]\n"
	  "# $9 = tangent point longitude [deg]\n"
	  "# $10 = tangent point latitude [deg]\n"
	  "# $11 = channel frequency [cm^-1]\n"
	  "# $12 = channel radiance [W/(m^2 sr cm^-1)]\n"
	  "# $13 = channel transmittance [1]\n");

  /* Write data... */
  for (int ir = 0; ir < obs->nr; ir++) {
    fprintf(out, "\n");
    for (int id = 0; id < ctl.nd; id++)
      fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %.4f %g %g\n",
	      obs->time[ir], obs->obsz[ir], obs->obslon[ir], obs->obslat[ir],
	      obs->vpz[ir], obs->vplon[ir], obs->vplat[ir], obs->tpz[ir],
	      obs->tplon[ir], obs->tplat[ir], ctl.nu[id], obs->rad[id][ir],
	      obs->tau[id][ir]);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(obs);

  return EXIT_SUCCESS;
}
