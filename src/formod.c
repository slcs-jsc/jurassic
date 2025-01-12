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
#ifdef UNIFIED
#include "jurassic_unified_library.h"
#endif

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

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <obs> <atm> <rad>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

#ifdef UNIFIED

  static atm_t atm;
  static obs_t obs;

  /* Read observation geometry... */
  read_obs(NULL, argv[2], &ctl, &obs);

  /* Read atmospheric data... */
  read_atm(NULL, argv[3], &ctl, &atm);

  /* Call forward model... */
  jur_unified_init(argc, argv);
  jur_unified_formod_multiple_packages(&atm, &obs, 1, NULL);

  /* Save radiance data... */
  write_obs(NULL, argv[4], &ctl, &obs);

#else

  char dirlist[LEN], task[LEN];

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
    FILE *in;
    if (!(in = fopen(dirlist, "r")))
      ERRMSG("Cannot open directory list!");

    /* Loop over directories... */
    char wrkdir[LEN];
    while (fscanf(in, "%s", wrkdir) != EOF) {

      /* Write info... */
      LOG(1, "\nWorking directory: %s", wrkdir);

      /* Call forward model... */
      call_formod(&ctl, wrkdir, argv[2], argv[3], argv[4], task);
    }

    /* Close dirlist... */
    fclose(in);
  }

#endif

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void call_formod(
  ctl_t *ctl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  const char *task) {

  static atm_t atm, atm2;
  static obs_t obs, obs2;

  /* Read observation geometry... */
  read_obs(wrkdir, obsfile, ctl, &obs);

  /* Read atmospheric data... */
  read_atm(wrkdir, atmfile, ctl, &atm);

  /* Compute multiple profiles... */
  if (task[0] == 'p' || task[0] == 'P') {

    /* Loop over ray paths... */
    for (int ir = 0; ir < obs.nr; ir++) {

      /* Get atmospheric data... */
      atm2.np = 0;
      for (int ip = 0; ip < atm.np; ip++)
	if (atm.time[ip] == obs.time[ir]) {
	  atm2.time[atm2.np] = atm.time[ip];
	  atm2.z[atm2.np] = atm.z[ip];
	  atm2.lon[atm2.np] = atm.lon[ip];
	  atm2.lat[atm2.np] = atm.lat[ip];
	  atm2.p[atm2.np] = atm.p[ip];
	  atm2.t[atm2.np] = atm.t[ip];
	  for (int ig = 0; ig < ctl->ng; ig++)
	    atm2.q[ig][atm2.np] = atm.q[ig][ip];
	  for (int iw = 0; iw < ctl->nw; iw++)
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
	for (int id = 0; id < ctl->nd; id++) {
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

      char filename[LEN];

      /* Switch off continua... */
      ctl->ctm_co2 = 0;
      ctl->ctm_h2o = 0;
      ctl->ctm_n2 = 0;
      ctl->ctm_o2 = 0;

      /* Loop over emitters... */
      for (int ig = 0; ig < ctl->ng; ig++) {

	/* Copy atmospheric data... */
	copy_atm(ctl, &atm2, &atm, 0);

	/* Set extinction to zero... */
	for (int iw = 0; iw < ctl->nw; iw++)
	  for (int ip = 0; ip < atm2.np; ip++)
	    atm2.k[iw][ip] = 0;

	/* Set volume mixing ratios to zero... */
	for (int ig2 = 0; ig2 < ctl->ng; ig2++)
	  if (ig2 != ig)
	    for (int ip = 0; ip < atm2.np; ip++)
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
      for (int ig = 0; ig < ctl->ng; ig++)
	for (int ip = 0; ip < atm2.np; ip++)
	  atm2.q[ig][ip] = 0;

      /* Call forward model... */
      formod(ctl, &atm2, &obs);

      /* Save radiance data... */
      sprintf(filename, "%s.EXTINCT", radfile);
      write_obs(wrkdir, filename, ctl, &obs);
    }

    /* Measure CPU-time... */
    if (task[0] == 't' || task[0] == 'T') {

      /* Init... */
      double t_min, t_max, t_mean = 0, t_sigma = 0;
      int n = 0;

      /* Initialize random number generator... */
      gsl_rng_env_setup();
      gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

      /* Loop over profiles... */
      do {

	/* Create random atmosphere... */
	copy_atm(ctl, &atm2, &atm, 0);
	double dtemp = 40. * (gsl_rng_uniform(rng) - 0.5);
	double dpress = 1. - 0.1 * gsl_rng_uniform(rng);
	double dq[NG];
	for (int ig = 0; ig < ctl->ng; ig++)
	  dq[ig] = 0.8 + 0.4 * gsl_rng_uniform(rng);
	for (int ip = 0; ip < atm2.np; ip++) {
	  atm.t[ip] += dtemp;
	  atm.p[ip] *= dpress;
	  for (int ig = 0; ig < ctl->ng; ig++)
	    atm.q[ig][ip] *= dq[ig];
	}

	/* Measure runtime... */
	double t0 = omp_get_wtime();
	formod(ctl, &atm2, &obs);
	double dt = omp_get_wtime() - t0;

	/* Get runtime statistics... */
	t_mean += (dt);
	t_sigma += POW2(dt);
	if (n == 0 || dt < t_min)
	  t_min = dt;
	if (n == 0 || dt > t_max)
	  t_max = dt;
	n++;

      } while (t_mean < 10.0);

      /* Write results... */
      t_mean /= (double) n;
      t_sigma = sqrt(t_sigma / (double) n - POW2(t_mean));
      printf("RUNTIME_MEAN = %g s\n", t_mean);
      printf("RUNTIME_SIGMA = %g s\n", t_sigma);
      printf("RUNTIME_MIN = %g s\n", t_min);
      printf("RUNTIME_MAX = %g s\n", t_max);
      printf("RAYS_PER_SECOND = %g", (double) obs.nr / t_mean);
    }

    /* Analyze effect of step size... */
    if (task[0] == 's' || task[0] == 'S') {

      /* Reference run... */
      ctl->rayds = 0.1;
      ctl->raydz = 0.01;
      formod(ctl, &atm, &obs);
      copy_obs(ctl, &obs2, &obs, 0);

      /* Loop over step size... */
      for (double dz = 0.01; dz <= 2; dz *= 1.1) {
	printf("STEPSIZE: \n");
	for (double ds = 0.1; ds <= 50; ds *= 1.1) {

	  /* Set step size... */
	  ctl->rayds = ds;
	  ctl->raydz = dz;

	  /* Measure runtime... */
	  double t0 = omp_get_wtime();
	  formod(ctl, &atm, &obs);
	  double dt = omp_get_wtime() - t0;

	  /* Get differences... */
	  double mean[ND], sigma[ND];
	  for (int id = 0; id < ctl->nd; id++) {
	    mean[id] = sigma[id] = 0;
	    int n = 0;
	    for (int ir = 0; ir < obs.nr; ir++) {
	      double err = 200. * (obs.rad[id][ir] - obs2.rad[id][ir])
		/ (obs.rad[id][ir] + obs2.rad[id][ir]);
	      mean[id] += err;
	      sigma[id] += POW2(err);
	      n++;
	    }
	    mean[id] /= n;
	    sigma[id] = sqrt(sigma[id] / n - POW2(mean[id]));
	  }

	  /* Write results... */
	  printf("STEPSIZE: %g %g %g", ds, dz, dt);
	  for (int id = 0; id < ctl->nd; id++)
	    printf(" %g %g", mean[id], sigma[id]);
	  printf("\n");
	}
      }
    }
  }
}
