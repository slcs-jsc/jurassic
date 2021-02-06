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
  
  Copyright (C) 2019-2021 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Inversion tool for MPTRAC.
*/

#include "jurassic.h"
#include <gsl/gsl_fit.h>

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Maximum number of data points... */
#define NMAX 1000

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  static atm_t atm, atm2;

  static obs_t obs;

  static gsl_matrix *k;

  static FILE *in, *out;

  static char line[LEN];

  static double rtime, rz, rlon, rlat, rp, rt, rso2, rh2o, ro3, robs,
    obs_min, obs_meas, obs_sim, scl = 1.0, scl_err, scl_old, c0, c1,
    cov00, cov01, cov11, sumsq, x[NMAX], x2[NMAX], y[NMAX], y2[NMAX],
    t0 = GSL_NAN, dt, tol;
  
  static int i, ig, ip, it, itmax, method, n, nmean[NMAX], nprof;

  static size_t mk, nk;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Give parameters: <ctl> <prof> <inv> <kernel>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  dt = scan_ctl(argc, argv, "INVERT_DT", -1, "86400", NULL);
  itmax = (int) scan_ctl(argc, argv, "INVERT_ITMAX", -1, "10", NULL);
  tol = scan_ctl(argc, argv, "INVERT_TOL", -1, "1e-4", NULL);
  method = (int) scan_ctl(argc, argv, "INVERT_METHOD", -1, "1", NULL);
  obs_min = scan_ctl(argc, argv, "INVERT_OBSMIN", -1, "-1e99", NULL);

  /* Check control parameters... */
  if (ctl.ng != 4)
    ERRMSG("Set NG = 4!");
  if (strcmp(ctl.emitter[0], "SO2") != 0)
    ERRMSG("Set EMITTER[0] = SO2!");
  if (strcmp(ctl.emitter[1], "H2O") != 0)
    ERRMSG("Set EMITTER[1] = H2O!");
  if (strcmp(ctl.emitter[2], "O3") != 0)
    ERRMSG("Set EMITTER[2] = O3!");
  if (strcmp(ctl.emitter[3], "CO2") != 0)
    ERRMSG("Set EMITTER[3] = CO2!");
  if (ctl.nd != 2)
    ERRMSG("Set ND = 2!");

  /* Set control parameters... */
  ctl.write_bbt = 1;
  ctl.write_matrix = 1;

  /* Set observation data... */
  obs.nr = 1;
  obs.obsz[0] = 705;

  /* ------------------------------------------------------------
     Fit scaling factor for total mass...
     ------------------------------------------------------------ */

  /* Iterations... */
  for (it = 0; it < itmax; it++) {

    /* Initialize... */
    atm.np = n = 0;
    if (method == 1 || method == 4)
      for (i = 0; i < NMAX; i++)
	x[i] = y[i] = GSL_NAN;
    else if (method == 2 || method == 3 || method == 5 || method == 6)
      for (i = 0; i < NMAX; i++) {
	x[i] = y[i] = 0;
	nmean[i] = 0;
      }
    else
      ERRMSG("Check INVERT_METHOD!");
    
    /* Read profile data... */
    printf("Read profile data: %s\n", argv[2]);

    /* Open file... */
    if (!(in = fopen(argv[2], "r")))
      ERRMSG("Cannot open file!");

    /* Write inversion data... */
    printf("Write inversion data: %s\n", argv[3]);

    /* Create file... */
    if (!(out = fopen(argv[3], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	    "# $2 = altitude [km]\n"
	    "# $3 = longitude [deg]\n"
	    "# $4 = latitude [deg]\n"
	    "# $5 = simulated SO2 index [K]\n"
	    "# $6 = measured SO2 index [K]\n\n");

    /* Read line... */
    while (fgets(line, LEN, in)) {

      /* Read data... */
      if (sscanf
	  (line, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &rtime, &rz,
	   &rlon, &rlat, &rp, &rt, &rso2, &rh2o, &ro3, &robs) != 10)
	continue;

      /* Save initial time... */
      if (!gsl_finite(t0))
	t0 = rtime;

      /* Check for new profile... */
      if ((rtime != atm.time[0] || rlon != atm.lon[0] || rlat != atm.lat[0])
	  && atm.np > 0) {

	/* Call forward model... */
	formod(&ctl, &atm, &obs);
	obs_sim = obs.rad[0][0] - obs.rad[1][0];

	/* Write output... */
	fprintf(out, "%.2f %g %g %g %g %g\n", atm.time[0], 0.0, atm.lon[0],
		atm.lat[0], obs_sim, obs_meas);

	/* Get time index... */
	i = (int) ((atm.time[0] - t0) / dt);
	if (i < 0 && i >= NMAX)
	  ERRMSG("Index out of range!");

	/* Get maxima... */
	if (method == 1 || method == 4) {
	  if (gsl_finite(x[i]))
	    x[i] = GSL_MAX(x[i], obs_sim);
	  else
	    x[i] = obs_sim;
	  if (gsl_finite(y[i]))
	    y[i] = GSL_MAX(y[i], obs_meas);
	  else
	    y[i] = obs_meas;
	}
	
	/* Get means... */
	else if ((method == 2 || method == 3 || method == 5 || method == 6)
		 && (obs_meas >= obs_min || obs_sim >= obs_min)) {
	  x[i] += obs_sim;
	  y[i] += obs_meas;
	  nmean[i]++;
	}

	/* Calculate mean atmospheric profile... */
	nprof++;
	atm2.np = atm.np;
	for (ip = 0; ip < atm.np; ip++) {
	  atm2.time[ip] += atm.time[ip];
	  atm2.z[ip] += atm.z[ip];
	  atm2.lon[ip] += atm.lon[ip];
	  atm2.lat[ip] += atm.lat[ip];
	  atm2.p[ip] += atm.p[ip];
	  atm2.t[ip] += atm.t[ip];
	  for (ig = 1; ig < ctl.ng; ig++)
	    atm2.q[ig][ip] += atm.q[ig][ip];
	}

	/* Reset counter... */
	atm.np = 0;
      }

      /* Save data... */
      obs_meas = robs;
      atm.time[atm.np] = rtime;
      atm.z[atm.np] = rz;
      atm.lon[atm.np] = rlon;
      atm.lat[atm.np] = rlat;
      atm.p[atm.np] = rp;
      atm.t[atm.np] = rt;
      atm.q[0][atm.np] = rso2 * scl;
      atm.q[1][atm.np] = rh2o;
      atm.q[2][atm.np] = ro3;
      atm.q[3][atm.np] = 371.789948e-6 + 2.026214e-6
	* (atm.time[atm.np] - 63158400.) / 31557600.;
      if ((++atm.np) > NP)
	ERRMSG("Too many data points!");
    }

    /* Calculate means... */
    if(method == 2 || method == 5)
      for (i = 0; i < NMAX; i++) {
	x[i] /= nmean[i];
	y[i] /= nmean[i];
      }
    
    /* Filter data... */
    n = 0;
    for (i = 0; i < NMAX; i++)
      if (gsl_finite(x[i]) && gsl_finite(y[i])) {
	x2[n] = x[i];
	y2[n] = y[i];
	n++;
      }

    /* Report data... */
    fprintf(out, "\n");
    for (i = 0; i < n; i++)
      fprintf(out, "# time= %.2f | si_sim= %g | si_obs= %g\n",
	      t0 + (i + 0.5) * dt, x2[i], y2[i]);

    /* Report statistics... */
    fprintf(out, "\n");
    fprintf(out, "# scl= %g +/- %g\n", scl, scl_err);
    fprintf(out, "# RMSE= %g\n", sqrt(sumsq / n));
    fprintf(out, "# n= %d\n", n);

    /* Close files... */
    fclose(out);
    fclose(in);

    /* Write info... */
    printf("  it= %d | scl= %g +/- %g | RMSE= %g\n", it, scl, scl_err,
	   sqrt(sumsq / n));

    /* Get new scaling factor... */
    if(method >=1 && method <= 3)
      gsl_fit_mul(x2, 1, y2, 1, (size_t) n, &c1, &cov11, &sumsq);
    else if(method >= 4 && method <= 6)
      gsl_fit_linear(x2, 1, y2, 1, (size_t) n, &c0, &c1,
		     &cov00, &cov01, &cov11, &sumsq);
    scl *= c1;
    scl_err = scl * sqrt(cov11 / POW2(c1));
    
    /* Convergence test... */
    if (fabs(2.0 * (scl - scl_old) / (scl + scl_old)) < tol)
      break;
    scl_old = scl;
  }

  /* ------------------------------------------------------------
     Calculate kernel...
     ------------------------------------------------------------ */

  /* Set atmospheric data... */
  for (ip = 0; ip < atm2.np; ip++) {
    atm2.time[ip] /= nprof;
    atm2.z[ip] /= nprof;
    atm2.lon[ip] /= nprof;
    atm2.lat[ip] /= nprof;
    atm2.p[ip] /= nprof;
    atm2.t[ip] /= nprof;
    for (ig = 0; ig < ctl.ng; ig++)
      atm2.q[ig][ip] /= nprof;
  }

  /* Get sizes... */
  nk = atm2x(&ctl, &atm2, NULL, NULL, NULL);
  mk = obs2y(&ctl, &obs, NULL, NULL, NULL);

  /* Allocate... */
  k = gsl_matrix_alloc(mk, nk);

  /* Compute kernel matrix... */
  kernel(&ctl, &atm2, &obs, k);

  /* Write matrix to file... */
  write_matrix(NULL, argv[4], &ctl, k, &atm2, &obs, "y", "x", "r");

  /* Free... */
  gsl_matrix_free(k);

  return EXIT_SUCCESS;
}
