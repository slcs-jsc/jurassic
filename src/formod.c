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
  JURASSIC forward model.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Perform forward model calculations in a single directory. */
void call_formod(
  ctl_t * ctl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  const char *task);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  FILE *in;

  char dirlist[LEN], task[LEN], wrkdir[LEN];

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <rad>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Get task... */
  scan_ctl(argc, argv, "TASK", -1, "-", task);

  /* Get dirlist... */
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);

  /* Single forward calculation... */
  if (dirlist[0] == '-')
    call_formod(&ctl, NULL, argv[2], argv[3], argv[4], task);

  /* Work on directory list... */
  else {

    /* Open directory list... */
    if (!(in = fopen(dirlist, "r")))
      ERRMSG("Cannot open directory list!");

    /* Loop over directories... */
    while (fscanf(in, "%s", wrkdir) != EOF) {

      /* Write info... */
      printf("\nWorking directory: %s\n", wrkdir);

      /* Call forward model... */
      call_formod(&ctl, wrkdir, argv[2], argv[3], argv[4], task);
    }

    /* Close dirlist... */
    fclose(in);
  }

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void call_formod(
  ctl_t * ctl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  const char *task) {

  static atm_t atm, atm2;
  static obs_t obs, obs2;

  char filename[LEN];

  int id, ig, ig2, ip, ir, iw;

  /* Read observation geometry... */
  read_obs(wrkdir, obsfile, ctl, &obs);

  /* Read atmospheric data... */
  read_atm(wrkdir, atmfile, ctl, &atm);

  /* Compute multiple profiles... */
  if (task[0] == 'p' || task[0] == 'P') {

    /* Loop over ray paths... */
    for (ir = 0; ir < obs.nr; ir++) {

      /* Get atmospheric data... */
      atm2.np = 0;
      for (ip = 0; ip < atm.np; ip++)
	if (atm.time[ip] == obs.time[ir]) {
	  atm2.time[atm2.np] = atm.time[ip];
	  atm2.z[atm2.np] = atm.z[ip];
	  atm2.lon[atm2.np] = atm.lon[ip];
	  atm2.lat[atm2.np] = atm.lat[ip];
	  atm2.p[atm2.np] = atm.p[ip];
	  atm2.t[atm2.np] = atm.t[ip];
	  for (ig = 0; ig < ctl->ng; ig++)
	    atm2.q[ig][atm2.np] = atm.q[ig][ip];
	  for (iw = 0; iw < ctl->nw; iw++)
	    atm2.k[iw][atm2.np] = atm.k[iw][ip];
	  atm2.np++;
	}

      /* Get observation data... */
      obs2.nr = 1;
      obs2.time[0] = obs.time[ir];
      obs2.vpz[0] = obs.vpz[ir];
      obs2.vplon[0] = obs.vplon[ir];
      obs2.vplat[0] = obs.vplat[ir];
      obs2.obsz[0] = obs.obsz[ir];
      obs2.obslon[0] = obs.obslon[ir];
      obs2.obslat[0] = obs.obslat[ir];

      /* Check number of data points... */
      if (atm2.np > 0) {

	/* Call forward model... */
	formod(ctl, &atm2, &obs2);

	/* Save radiance data... */
	for (id = 0; id < ctl->nd; id++) {
	  obs.rad[id][ir] = obs2.rad[id][0];
	  obs.tau[id][ir] = obs2.tau[id][0];
	}
      }
    }

    /* Write radiance data... */
    write_obs(wrkdir, radfile, ctl, &obs);
  }

  /* Compute single profile... */
  else {

    /* Call forward model... */
    formod(ctl, &atm, &obs);

    /* Save radiance data... */
    write_obs(wrkdir, radfile, ctl, &obs);

    /* Compute contributions... */
    if (task[0] == 'c' || task[0] == 'C') {

      /* Switch off continua... */
      ctl->ctm_co2 = 0;
      ctl->ctm_h2o = 0;
      ctl->ctm_n2 = 0;
      ctl->ctm_o2 = 0;

      /* Loop over emitters... */
      for (ig = 0; ig < ctl->ng; ig++) {

	/* Copy atmospheric data... */
	copy_atm(ctl, &atm2, &atm, 0);

	/* Set extinction to zero... */
	for (iw = 0; iw < ctl->nw; iw++)
	  for (ip = 0; ip < atm2.np; ip++)
	    atm2.k[iw][ip] = 0;

	/* Set volume mixing ratios to zero... */
	for (ig2 = 0; ig2 < ctl->ng; ig2++)
	  if (ig2 != ig)
	    for (ip = 0; ip < atm2.np; ip++)
	      atm2.q[ig2][ip] = 0;

	/* Call forward model... */
	formod(ctl, &atm2, &obs);

	/* Save radiance data... */
	sprintf(filename, "%s.%s", radfile, ctl->emitter[ig]);
	write_obs(wrkdir, filename, ctl, &obs);
      }

      /* Copy atmospheric data... */
      copy_atm(ctl, &atm2, &atm, 0);

      /* Set volume mixing ratios to zero... */
      for (ig = 0; ig < ctl->ng; ig++)
	for (ip = 0; ip < atm2.np; ip++)
	  atm2.q[ig][ip] = 0;

      /* Call forward model... */
      formod(ctl, &atm2, &obs);

      /* Save radiance data... */
      sprintf(filename, "%s.EXTINCT", radfile);
      write_obs(wrkdir, filename, ctl, &obs);
    }

    /* Measure CPU-time... */
    if (task[0] == 't' || task[0] == 'T') {
      TIMER("formod", 1);
      formod(ctl, &atm, &obs);
      TIMER("formod", 3);
    }
  }
}
