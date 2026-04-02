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
  const tbl_t * tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  const char *task,
  const char *obsref);

/*! Calculate relative errors. */
void compute_rel_errors(
  const ctl_t * ctl,
  const obs_t * obs_test,
  const obs_t * obs_ref,
  double *mre,
  double *sdre,
  double *minre,
  double *maxre);

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
  SELECT_TIMER("READ_CTL", "INPUT");
  read_ctl(argc, argv, &ctl);

#ifdef UNIFIED

  static atm_t atm;
  static obs_t obs;

  /* Read observation geometry... */
  SELECT_TIMER("READ_OBS", "INPUT");
  read_obs(NULL, argv[2], &ctl, &obs, 0);

  /* Read atmospheric data... */
  SELECT_TIMER("READ_ATM", "INPUT");
  read_atm(NULL, argv[3], &ctl, &atm, 0);

  /* Call forward model... */
  SELECT_TIMER("FORMOD_UNIFIED", "FORWARD");
  jur_unified_init(argc, argv);
  jur_unified_formod_multiple_packages(&atm, &obs, 1, NULL);

  /* Save radiance data... */
  SELECT_TIMER("WRITE_OBS", "OUTPUT");
  write_obs(NULL, argv[4], &ctl, &obs, 0);

#else

  char dirlist[LEN], obsref[LEN], task[LEN];

  /* Initialize look-up tables... */
  SELECT_TIMER("READ_TBL", "INPUT");
  tbl_t *tbl = read_tbl(&ctl);

  /* Get task... */
  SELECT_TIMER("READ_TASK", "INPUT");
  scan_ctl(argc, argv, "TASK", -1, "-", task);

  /* Get dirlist... */
  SELECT_TIMER("READ_DIRLIST", "INPUT");
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);

  /* Get reference data... */
  SELECT_TIMER("READ_OBSREF", "INPUT");
  scan_ctl(argc, argv, "OBSREF", -1, "-", obsref);

  /* Single forward calculation... */
  if (dirlist[0] == '-')
    call_formod(&ctl, tbl, NULL, argv[2], argv[3], argv[4], task, obsref);

  /* Work on directory list... */
  else {

    /* Open directory list... */
    FILE *in;
    if (!(in = fopen(dirlist, "r")))
      ERRMSG("Cannot open directory list!");

    /* Loop over directories... */
    char wrkdir[LEN];
    while (fscanf(in, "%4999s", wrkdir) != EOF) {

      /* Write info... */
      LOG(1, "\nWorking directory: %s", wrkdir);

      /* Call forward model... */
      call_formod(&ctl, tbl, wrkdir, argv[2], argv[3], argv[4], task, obsref);
    }

    /* Close dirlist... */
    fclose(in);
  }

#endif

  /* Free... */
  SELECT_TIMER("FINALIZE", "OVERHEAD");
  tbl_free(&ctl, tbl);
  PRINT_TIMERS;

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void call_formod(
  ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  const char *task,
  const char *obsref) {

  static atm_t atm, atm2;
  static obs_t obs, obs2;

  /* Read atmospheric data... */
  SELECT_TIMER("READ_ATM", "INPUT");
  read_atm(wrkdir, atmfile, ctl, &atm, 0);

  /* Read observation geometry... */
  SELECT_TIMER("READ_OBS", "INPUT");
  read_obs(wrkdir, obsfile, ctl, &obs, 0);

  /* Compute multiple profiles... */
  if (task[0] == 'p' || task[0] == 'P') {

    SELECT_TIMER("FORMOD_PROFILE_LOOP", "FORWARD");

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
	formod(ctl, tbl, &atm2, &obs2);

	/* Save radiance data... */
	for (int id = 0; id < ctl->nd; id++) {
	  obs.rad[id][ir] = obs2.rad[id][0];
	  obs.tau[id][ir] = obs2.tau[id][0];
	}
      }
    }

    /* Write radiance data... */
    SELECT_TIMER("WRITE_OBS", "OUTPUT");
    write_obs(wrkdir, radfile, ctl, &obs, 0);
  }

  /* Compute single profile... */
  else {

    /* Call forward model... */
    SELECT_TIMER("FORMOD", "FORWARD");
    formod(ctl, tbl, &atm, &obs);

    /* Save radiance data... */
    SELECT_TIMER("WRITE_OBS", "OUTPUT");
    write_obs(wrkdir, radfile, ctl, &obs, 0);

    /* Evaluate results... */
    if (obsref[0] != '-') {

      SELECT_TIMER("EVAL", "OUTPUT");

      /* Read reference data... */
      read_obs(wrkdir, obsref, ctl, &obs2, 0);

      /* Calculate relative errors... */
      double mre[ND], sdre[ND], minre[ND], maxre[ND];
      compute_rel_errors(ctl, &obs, &obs2, mre, sdre, minre, maxre);

      /* Write results... */
      for (int id = 0; id < ctl->nd; id++)
	printf
	  ("EVAL: nu= %.4f cm^-1 | MRE= %g %% | SDRE= %g %% | MinRE= %g %% | MaxRE= %g %%\n",
	   ctl->nu[id], mre[id], sdre[id], minre[id], maxre[id]);
    }

    /* Compute contributions of emitters... */
    if (task[0] == 'c' || task[0] == 'C') {

      SELECT_TIMER("CONTRIBUTIONS", "FORWARD");

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

	/* Select emitter... */
	for (int ig2 = 0; ig2 < ctl->ng; ig2++)
	  if (ig2 != ig)
	    for (int ip = 0; ip < atm2.np; ip++)
	      atm2.q[ig2][ip] = 0;

	/* Call forward model... */
	formod(ctl, tbl, &atm2, &obs);

	/* Save radiance data... */
	sprintf(filename, "%s.%s", radfile, ctl->emitter[ig]);
	write_obs(wrkdir, filename, ctl, &obs, 0);
      }

      /* Copy atmospheric data... */
      copy_atm(ctl, &atm2, &atm, 0);

      /* Set volume mixing ratios to zero, keep extinction... */
      for (int ig = 0; ig < ctl->ng; ig++)
	for (int ip = 0; ip < atm2.np; ip++)
	  atm2.q[ig][ip] = 0;

      /* Call forward model... */
      formod(ctl, tbl, &atm2, &obs);

      /* Save radiance data... */
      sprintf(filename, "%s.EXTINCT", radfile);
      write_obs(wrkdir, filename, ctl, &obs, 0);
    }

    /* Measure CPU-time... */
    if (task[0] == 't' || task[0] == 'T') {

      /* Init... */
      double t_min = 0, t_max = 0, t_mean = 0, t_sd = 0;
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
	  atm2.t[ip] += dtemp;
	  atm2.p[ip] *= dpress;
	  for (int ig = 0; ig < ctl->ng; ig++)
	    atm2.q[ig][ip] *= dq[ig];
	}

	/* Measure runtime... */
	SELECT_TIMER("BENCHMARK_SAMPLE", "ANALYSIS");
	double t0 = omp_get_wtime();
	formod(ctl, tbl, &atm2, &obs);
	double dt = omp_get_wtime() - t0;

	/* Get runtime statistics... */
	t_mean += dt;
	t_sd += POW2(dt);
	if (n == 0 || dt < t_min)
	  t_min = dt;
	if (n == 0 || dt > t_max)
	  t_max = dt;
	n++;

      } while (t_mean < 10.0);

      /* Keep the legacy benchmark line and the new timer-style overlap metric. */
      t_mean /= (double) n;
      t_sd = sqrt(t_sd / (double) n - POW2(t_mean));
      printf("RUNTIME: mean= %g s | stddev= %g s | min= %g s | max= %g s\n",
	     t_mean, t_sd, t_min, t_max);

      /* Free... */
      gsl_rng_free(rng);
    }

    /* Analyze impact of step size... */
    if (task[0] == 's' || task[0] == 'S') {

      SELECT_TIMER("STEP_ANALYSIS", "ANALYSIS");

      /* Reference run... */
      ctl->rayds = 0.1;
      ctl->raydz = 0.01;
      formod(ctl, tbl, &atm, &obs);
      copy_obs(ctl, &obs2, &obs, 0);

      /* Loop over step size... */
      for (double dz = 0.01; dz <= 2; dz *= 1.1)
	for (double ds = 0.1; ds <= 50; ds *= 1.1) {

	  /* Set step size... */
	  ctl->rayds = ds;
	  ctl->raydz = dz;

	  /* Measure runtime... */
	  double t0 = omp_get_wtime();
	  formod(ctl, tbl, &atm, &obs);
	  double dt = omp_get_wtime() - t0;

	  /* Calculate relative errors... */
	  double mre[ND], sdre[ND], minre[ND], maxre[ND];
	  compute_rel_errors(ctl, &obs, &obs2, mre, sdre, minre, maxre);

	  /* Write results... */
	  for (int id = 0; id < ctl->nd; id++)
	    printf
	      ("STEPSIZE: ds= %.4f km | dz= %g km | t= %g s | nu= %.4f cm^-1"
	       " | MRE= %g %% | SDRE= %g %% | MinRE= %g %% | MaxRE= %g %%\n",
	       ds, dz, dt, ctl->nu[id], mre[id], sdre[id], minre[id],
	       maxre[id]);
	}
    }
  }
}

/*****************************************************************************/

void compute_rel_errors(
  const ctl_t *ctl,
  const obs_t *obs_test,
  const obs_t *obs_ref,
  double *mre,
  double *sdre,
  double *minre,
  double *maxre) {

  /* Loop over channels... */
  for (int id = 0; id < ctl->nd; id++) {

    double sum = 0, sum2 = 0;
    minre[id] = +1e9;
    maxre[id] = -1e9;
    int n = 0;

    /* Loop over ray paths... */
    for (int ir = 0; ir < obs_test->nr; ir++) {

      /* Check for zero... */
      if (obs_ref->rad[id][ir] == 0)
	continue;

      /* Calculate relative error... */
      double err = 100.0 * (obs_test->rad[id][ir] - obs_ref->rad[id][ir])
	/ obs_ref->rad[id][ir];

      /* Get statistics... */
      sum += err;
      sum2 += err * err;
      if (err > maxre[id])
	maxre[id] = err;
      if (err < minre[id])
	minre[id] = err;
      n++;
    }

    /* Get mean and standard deviaton... */
    mre[id] = sum / n;
    sdre[id] = sqrt(sum2 / n - mre[id] * mre[id]);
  }
}
