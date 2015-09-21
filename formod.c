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
  static obs_t obs;

  char filename[LEN];

  int ig, ig2, ip, iw;

  /* Read observation geometry... */
  read_obs(wrkdir, obsfile, ctl, &obs);

  /* Read atmospheric data... */
  read_atm(wrkdir, atmfile, ctl, &atm);

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
