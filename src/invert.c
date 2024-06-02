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

/*! Maximum number of data lines... */
#define NLMAX 30000000

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

  static double rtime[NLMAX], rz[NLMAX], rlon[NLMAX], rlat[NLMAX], obs_meas,
    obs_sim, scl = 1.0, scl_err, c0, c1, cov00, cov01, cov11, sumsq,
    x[NMAX], x2[NMAX], y[NMAX], y_err[NMAX], y2[NMAX], y2_err[NMAX],
    y2_sim[NMAX], y2_sim_err[NMAX], w2[NMAX], dt, tol, obs_err;

  static float rp[NLMAX], rt[NLMAX], rso2[NLMAX], rh2o[NLMAX],
    ro3[NLMAX], robs[NLMAX];

  static int data, fit, i, ig, il, ip, it, itmax, n, nl, ndata[NMAX], nprof;

  static size_t mk, nk;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Give parameters: <ctl> <prof> <inv> <atm> <kernel>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  dt = scan_ctl(argc, argv, "INVERT_DT", -1, "86400", NULL);
  obs_err = scan_ctl(argc, argv, "INVERT_OBS_ERR", -1, "1.0", NULL);
  data = (int) scan_ctl(argc, argv, "INVERT_DATA", -1, "2", NULL);
  fit = (int) scan_ctl(argc, argv, "INVERT_FIT", -1, "3", NULL);
  itmax = (int) scan_ctl(argc, argv, "INVERT_ITMAX", -1, "10", NULL);
  tol = scan_ctl(argc, argv, "INVERT_TOL", -1, "1e-4", NULL);

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
     Read profiles...
     ------------------------------------------------------------ */

  /* Read profile data... */
  LOG(1, "Read profile data: %s", argv[2]);

  /* Open file... */
  if (!(in = fopen(argv[2], "r")))
    ERRMSG("Cannot open file!");

  /* Read file... */
  while (fgets(line, LEN, in)) {

    /* Read data... */
    if (sscanf(line, "%lg %lg %lg %lg %g %g %g %g %g %g", &rtime[nl],
	       &rz[nl], &rlon[nl], &rlat[nl], &rp[nl], &rt[nl], &rso2[nl],
	       &rh2o[nl], &ro3[nl], &robs[nl]) != 10)
      continue;
    if ((++nl) > NLMAX)
      ERRMSG("Too many profile data points!");
  }

  /* Close files... */
  fclose(in);

  /* ------------------------------------------------------------
     Fit scaling factor for total mass...
     ------------------------------------------------------------ */

  /* Iterations... */
  for (it = 0; it < itmax; it++) {

    /* Init... */
    atm.np = n = 0;
    for (i = 0; i < NMAX; i++) {
      ndata[i] = 0;
      x[i] = y[i] = NAN;
    }

    /* Loop over lines... */
    for (il = 0; il < nl; il++) {

      /* Check for new profile... */
      if ((rtime[il] != atm.time[0]
	   || rlon[il] != atm.lon[0]
	   || rlat[il] != atm.lat[0])
	  && atm.np > 0) {

	/* Call forward model... */
	formod(&ctl, &atm, &obs);
	obs_sim = obs.rad[0][0] - obs.rad[1][0];

	/* Get time index... */
	i = (int) ((atm.time[0] - rtime[0]) / dt);
	if (i < 0 && i >= NMAX)
	  ERRMSG("Time index out of range!");

	/* Get maxima... */
	if (data == 1) {
	  x[i] = (isfinite(x[i]) ? MAX(x[i], obs_sim) : obs_sim);
	  y[i] = (isfinite(y[i]) ? MAX(y[i], obs_meas) : obs_meas);
	  y_err[i] = obs_err;
	  if (isfinite(x[i]) && isfinite(y[i]))
	    ndata[i] = 1;
	}

	/* Get means... */
	else if (data == 2) {
	  if (ndata[i] == 0) {
	    x[i] = obs_sim;
	    y[i] = obs_meas;
	    y_err[i] = POW2(obs_meas);
	  } else {
	    x[i] += obs_sim;
	    y[i] += obs_meas;
	    y_err[i] += POW2(obs_meas);
	  }
	  ndata[i]++;
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
      obs_meas = robs[il];
      atm.time[atm.np] = rtime[il];
      atm.z[atm.np] = rz[il];
      atm.lon[atm.np] = rlon[il];
      atm.lat[atm.np] = rlat[il];
      atm.p[atm.np] = rp[il];
      atm.t[atm.np] = rt[il];
      atm.q[0][atm.np] = rso2[il] * scl;
      atm.q[1][atm.np] = rh2o[il];
      atm.q[2][atm.np] = ro3[il];
      atm.q[3][atm.np] = 371.789948e-6 + 2.026214e-6
	* (atm.time[atm.np] - 63158400.) / 31557600.;
      if ((++atm.np) > NP)
	ERRMSG("Too many data points!");
    }

    /* Calculate means... */
    if (data == 2)
      for (i = 0; i < NMAX; i++)
	if (ndata[i] > 0) {
	  x[i] /= ndata[i];
	  y[i] /= ndata[i];
	  y_err[i] = sqrt(MAX(y_err[i] / ndata[i] - POW2(y[i]), 0.0))
	    / sqrt(ndata[i]);	/* standard error! */
	}

    /* Filter data... */
    n = 0;
    for (i = 0; i < NMAX; i++)
      if (ndata[i] > 0 && isfinite(x[i]) && isfinite(y[i])
	  && isfinite(y_err[i])) {
	x2[n] = x[i];
	y2[n] = y[i];
	y2_err[n] = y_err[i];
	w2[n] = 1. / POW2(y_err[i]);
	n++;
      }

    /* Fit radiance data... */
    if (fit == 1)
      gsl_fit_mul(x2, 1, y2, 1, (size_t) n, &c1, &cov11, &sumsq);
    else if (fit == 2)
      gsl_fit_wmul(x2, 1, w2, 1, y2, 1, (size_t) n, &c1, &cov11, &sumsq);
    else if (fit == 3)
      gsl_fit_linear(x2, 1, y2, 1, (size_t) n, &c0, &c1, &cov00, &cov01,
		     &cov11, &sumsq);
    else if (fit == 4)
      gsl_fit_wlinear(x2, 1, w2, 1, y2, 1, (size_t) n, &c0, &c1, &cov00,
		      &cov01, &cov11, &sumsq);
    else
      ERRMSG("Check INVERT_FIT!");

    /* Get new scaling factor... */
    double scl_old = scl;
    scl_err = scl * sqrt(cov11);
    scl *= c1;

    /* Write info... */
    LOG(1, "  it= %d | scl= %g +/- %g | RMSE= %g",
	it, scl, scl_err, sqrt(sumsq / n));

    /* Convergence test... */
    if (fabs(2.0 * (scl - scl_old) / (scl + scl_old)) < tol)
      break;
  }

  /* ------------------------------------------------------------
     Write inversion data...
     ------------------------------------------------------------ */

  /* Write info... */
  LOG(1, "Write inversion data: %s", argv[3]);

  /* Create file... */
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 2000-01-01T00:00Z)\n"
	  "# $2 = simulated SO2 index [K]\n"
	  "# $3 = scaled simulated SO2 index [K]\n"
	  "# $4 = error of scaled simulated SO2 index [K]\n"
	  "# $5 = observed SO2 index [K]\n"
	  "# $6 = error of observed SO2 index [K]\n\n");

  /* Write data... */
  for (i = 0; i < n; i++) {

    /* Calculate scaled SO2 index... */
    if (fit == 1 || fit == 2)
      gsl_fit_mul_est(x2[i], c1, cov11, &y2_sim[i], &y2_sim_err[i]);
    else if (fit == 3 || fit == 4)
      gsl_fit_linear_est(x2[i], c0, c1, cov00, cov01, cov11, &y2_sim[i],
			 &y2_sim_err[i]);

    /* Write output... */
    fprintf(out, "%.2f %g %g %g %g %g\n", rtime[0] + (i + 0.5) * dt,
	    x2[i], y2_sim[i], y2_sim_err[i], y2[i], y2_err[i]);
  }

  /* Report scaling factor for total mass... */
  fprintf(out, "\n");
  fprintf(out, "#    scl= %g +/- %g\n", scl, scl_err);
  fprintf(out, "#     c1= %g +/- %g\n", c1, sqrt(cov11));
  if (fit == 3 || fit == 4) {
    fprintf(out, "#     c0= %g +/- %g\n", c0, sqrt(cov00));
    fprintf(out, "#   corr= %g\n", cov01 / (sqrt(cov00) * sqrt(cov11)));
  }
  fprintf(out, "#   RMSE= %g\n", sqrt(sumsq / n));
  fprintf(out, "#      n= %d\n", n);

  /* Close files... */
  fclose(out);

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

  /* Write atmospheric data... */
  write_atm(NULL, argv[4], &ctl, &atm);

  /* Write matrix to file... */
  write_matrix(NULL, argv[5], &ctl, k, &atm2, &obs, "y", "x", "r");

  /* Free... */
  gsl_matrix_free(k);

  return EXIT_SUCCESS;
}
