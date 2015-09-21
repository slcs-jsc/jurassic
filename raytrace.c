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
  Determine atmospheric ray paths.
*/

#include "jurassic.h"

int main(
  int argc,
  char *argv[]) {

  static atm_t atm, atm2;
  static ctl_t ctl;
  static los_t los;
  static obs_t obs;

  FILE *out;

  char filename[LEN], losbase[LEN], massbase[LEN];

  double s;

  int ig, ip, ir, iw;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <obs> <atm>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Get basenames... */
  scan_ctl(argc, argv, "LOSBASE", -1, "los", losbase);
  scan_ctl(argc, argv, "MASSBASE", -1, "-", massbase);

  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs);

  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm);

  /* Create output file... */
  if (!(out = fopen("raytrace.tab", "w")))
    ERRMSG("Cannot create raytrace.tab!");

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
	  "# $11 = ray path index\n" "# $12 = ray path length [km]\n");
  for (ig = 0; ig < ctl.ng; ig++)
    fprintf(out, "# $%d = %s column density [molec/cm^2]\n",
	    13 + ig, ctl.emitter[ig]);
  fprintf(out, "\n");

  /* Loop over rays... */
  for (ir = 0; ir < obs.nr; ir++) {

    /* Raytracing... */
    raytrace(&ctl, &atm, &obs, &los, ir);

    /* Get ray path length... */
    s = 0;
    for (ip = 0; ip < los.np; ip++)
      s += los.ds[ip];

    /* Write to file... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %d %g",
	    obs.time[ir], obs.obsz[ir], obs.obslon[ir], obs.obslat[ir],
	    obs.vpz[ir], obs.vplon[ir], obs.vplat[ir], obs.tpz[ir],
	    obs.tplon[ir], obs.tplat[ir], ir, s);
    for (ig = 0; ig < ctl.ng; ig++)
      fprintf(out, " %g", los.cgu[ig][los.np - 1]);
    fprintf(out, "\n");

    /* Atmospheric data... */
    if (losbase[0] != '-') {

      /* Copy data... */
      atm2.np = los.np;
      for (ip = 0; ip < los.np; ip++) {
	atm2.time[ip] = obs.time[ir];
	atm2.z[ip] = los.z[ip];
	atm2.lon[ip] = los.lon[ip];
	atm2.lat[ip] = los.lat[ip];
	atm2.p[ip] = los.p[ip];
	atm2.t[ip] = los.t[ip];
	for (ig = 0; ig < ctl.ng; ig++)
	  atm2.q[ig][ip] = los.q[ig][ip];
	for (iw = 0; iw < ctl.nw; iw++)
	  atm2.k[iw][ip] = los.k[iw][ip];
      }

      /* Save data... */
      sprintf(filename, "%s.%d", losbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);
    }

    /* Number and column densities... */
    if (massbase[0] != '-') {

      /* Copy data... */
      for (ip = 0; ip < los.np; ip++) {
	atm2.p[ip] = atm2.t[ip] = 0;
	for (ig = 0; ig < ctl.ng; ig++)
	  atm2.q[ig][ip] = los.u[ig][ip] / (los.ds[ip] * 1e5);
	for (iw = 0; iw < ctl.nw; iw++)
	  atm2.k[iw][ip] = 0;
      }

      /* Save data... */
      sprintf(filename, "%s.numdens.%d", massbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);

      /* Copy data... */
      for (ip = 0; ip < los.np; ip++)
	for (ig = 0; ig < ctl.ng; ig++)
	  atm2.q[ig][ip] = los.cgu[ig][ip];

      /* Save data... */
      sprintf(filename, "%s.coldens.%d", massbase, ir);
      write_atm(NULL, filename, &ctl, &atm2);
    }
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
