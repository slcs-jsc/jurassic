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
  Determine atmospheric ray paths.
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

  static atm_t atm;
  static ctl_t ctl;
  static los_t los;
  static obs_t obs;

  FILE *out;

  char filename[2 * LEN], losbase[LEN];

  double u[NG];

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: raytrace <ctl> <obs> <atm> <raytrace.tab> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Get basenames... */
  scan_ctl(argc, argv, "LOSBASE", -1, "los", losbase);

  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs, 0);

  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm, 0);

  /* Write info... */
  LOG(1, "Write raytrace data: %s", argv[4]);

  /* Create file... */
  if (!(out = fopen(argv[4], "w")))
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
	  "# $11 = ray path index\n" "# $12 = ray path length [km]\n");
  for (int ig = 0; ig < ctl.ng; ig++)
    fprintf(out, "# $%d = %s column density [molec/cm^2]\n",
	    13 + ig, ctl.emitter[ig]);
  fprintf(out, "\n");

  /* Loop over rays... */
  for (int ir = 0; ir < obs.nr; ir++) {

    /* Raytracing... */
    raytrace(&ctl, &atm, &obs, &los, ir);

    /* Set filename data... */
    sprintf(filename, "%s.%d.tab", losbase, ir);

    /* Write info... */
    LOG(1, "Write LOS data: %s", filename);

    /* Create file... */
    FILE *out2;
    if (!(out2 = fopen(filename, "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out2,
	    "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = pressure [hPa]\n" "# $6 = temperature [K]\n");
    for (int ig = 0; ig < ctl.ng; ig++)
      fprintf(out2, "# $%d = %s volume mixing ratio [ppv]\n",
	      7 + ig, ctl.emitter[ig]);
    for (int iw = 0; iw < ctl.nw; iw++)
      fprintf(out2, "# $%d = extinction (window %d) [km^-1]\n",
	      7 + ctl.ng + iw, iw);
    fprintf(out2, "\n");

    /* Write data... */
    for (int ip = 0; ip < los.np; ip++) {
      fprintf(out2, "%.2f %g %g %g %g %g", obs.time[ir], los.z[ip],
	      los.lon[ip], los.lat[ip], los.p[ip], los.t[ip]);
      for (int ig = 0; ig < ctl.ng; ig++)
	fprintf(out2, " %g", los.q[ip][ig]);
      for (int iw = 0; iw < ctl.nw; iw++)
	fprintf(out2, " %g", los.k[ip][iw]);
      fprintf(out2, "\n");
    }

    /* Close file... */
    fclose(out2);

    /* Get column densities... */
    double s = 0;
    for (int ig = 0; ig < ctl.ng; ig++)
      u[ig] = 0;
    for (int ip = 0; ip < los.np; ip++) {
      s += los.ds[ip];
      for (int ig = 0; ig < ctl.ng; ig++)
	u[ig] += los.u[ip][ig];
    }

    /* Write summary data... */
    fprintf(out, "%.2f %g %g %g %g %g %g %g %g %g %d %g",
	    obs.time[ir], obs.obsz[ir], obs.obslon[ir], obs.obslat[ir],
	    obs.vpz[ir], obs.vplon[ir], obs.vplat[ir],
	    obs.tpz[ir], obs.tplon[ir], obs.tplat[ir], ir, s);
    for (int ig = 0; ig < ctl.ng; ig++)
      fprintf(out, " %g", u[ig]);
    fprintf(out, "\n");
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC ray-tracing tool.\n\n");
  printf("Determine atmospheric line-of-sight paths and write a summary\n");
  printf("table plus one LOS profile file per ray.\n\n");
  printf("Usage:\n");
  printf("  raytrace <ctl> <obs> <atm> <raytrace.tab> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>           Control file.\n");
  printf("  <obs>           Observation geometry input file.\n");
  printf("  <atm>           Atmospheric state input file.\n");
  printf("  <raytrace.tab>  Output summary table for traced rays.\n");
  printf("  [KEY VALUE]     Optional control parameters.\n\n");
  printf("Optional control overrides:\n");
  printf("  LOSBASE <name>  Basename for per-ray LOS output files.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
