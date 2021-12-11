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
  JURASSIC retrieval processor.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Retrieval control parameters. */
typedef struct {

  /*! Working directory. */
  char dir[LEN];

  /*! Re-computation of kernel matrix (number of iterations). */
  int kernel_recomp;

  /*! Maximum number of iterations. */
  int conv_itmax;

  /*! Minimum normalized step size in state space. */
  double conv_dmin;

  /*! Carry out error analysis (0=no, 1=yes). */
  int err_ana;

  /*! Forward model error [%]. */
  double err_formod[ND];

  /*! Noise error [W/(m^2 sr cm^-1)]. */
  double err_noise[ND];

  /*! Pressure error [%]. */
  double err_press;

  /*! Vertical correlation length for pressure error [km]. */
  double err_press_cz;

  /*! Horizontal correlation length for pressure error [km]. */
  double err_press_ch;

  /*! Temperature error [K]. */
  double err_temp;

  /*! Vertical correlation length for temperature error [km]. */
  double err_temp_cz;

  /*! Horizontal correlation length for temperature error [km]. */
  double err_temp_ch;

  /*! Volume mixing ratio error [%]. */
  double err_q[NG];

  /*! Vertical correlation length for volume mixing ratio error [km]. */
  double err_q_cz[NG];

  /*! Horizontal correlation length for volume mixing ratio error [km]. */
  double err_q_ch[NG];

  /*! Extinction error [1/km]. */
  double err_k[NW];

  /*! Vertical correlation length for extinction error [km]. */
  double err_k_cz[NW];

  /*! Horizontal correlation length for extinction error [km]. */
  double err_k_ch[NW];

  /*! Cloud height error [km]. */
  double err_clz;

  /*! Cloud depth error [km]. */
  double err_cldz;

  /*! Cloud extinction error [1/km]. */
  double err_clk[NCL];

  /*! Surface height error [km]. */
  double err_sfz;

  /*! Surface pressure error [hPa]. */
  double err_sfp;

  /*! Surface temperature error [K]. */
  double err_sft;

  /*! Surface emissivity error. */
  double err_sfeps[NSF];

} ret_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Compute information content and resolution. */
void analyze_avk(
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  int *iqa,
  int *ipa,
  gsl_matrix * avk);

/*! Analyze averaging kernels for individual retrieval target. */
void analyze_avk_quantity(
  gsl_matrix * avk,
  int iq,
  int *ipa,
  size_t *n0,
  size_t *n1,
  double *cont,
  double *res);

/*! Compute cost function. */
double cost_function(
  gsl_vector * dx,
  gsl_vector * dy,
  gsl_matrix * s_a_inv,
  gsl_vector * sig_eps_inv);

/*! Invert symmetric matrix. */
void matrix_invert(
  gsl_matrix * a);

/*! Compute matrix product A^TBA or ABA^T for diagonal matrix B. */
void matrix_product(
  gsl_matrix * a,
  gsl_vector * b,
  int transpose,
  gsl_matrix * c);

/*! Carry out optimal estimation retrieval. */
void optimal_estimation(
  ret_t * ret,
  ctl_t * ctl,
  obs_t * obs_meas,
  obs_t * obs_i,
  atm_t * atm_apr,
  atm_t * atm_i);

/*! Read retrieval control parameters. */
void read_ret(
  int argc,
  char *argv[],
  ctl_t * ctl,
  ret_t * ret);

/*! Set a priori covariance. */
void set_cov_apr(
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  int *iqa,
  int *ipa,
  gsl_matrix * s_a);

/*! Set measurement errors. */
void set_cov_meas(
  ret_t * ret,
  ctl_t * ctl,
  obs_t * obs,
  gsl_vector * sig_noise,
  gsl_vector * sig_formod,
  gsl_vector * sig_eps_inv);

/*! Write retrieval error to file. */
void write_stddev(
  const char *quantity,
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  gsl_matrix * s);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static atm_t atm_i, atm_apr;
  static ctl_t ctl;
  static obs_t obs_i, obs_meas;
  static ret_t ret;

  FILE *dirlist;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <dirlist>");

  /* Measure CPU-time... */
  TIMER("total", 1);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  read_ret(argc, argv, &ctl, &ret);

  /* Open directory list... */
  if (!(dirlist = fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%s", ret.dir) != EOF) {

    /* Write info... */
    LOG(1, "\nRetrieve in directory %s...\n", ret.dir);

    /* Read atmospheric data... */
    read_atm(ret.dir, "atm_apr.tab", &ctl, &atm_apr);

    /* Read observation data... */
    read_obs(ret.dir, "obs_meas.tab", &ctl, &obs_meas);

    /* Run retrieval... */
    optimal_estimation(&ret, &ctl, &obs_meas, &obs_i, &atm_apr, &atm_i);

    /* Measure CPU-time... */
    TIMER("total", 2);
  }

  /* Write info... */
  LOG(1, "\nRetrieval done...");

  /* Measure CPU-time... */
  TIMER("total", 3);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void analyze_avk(
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  int *iqa,
  int *ipa,
  gsl_matrix * avk) {

  static atm_t atm_cont, atm_res;

  int icl, ig, iq, isf, iw;

  size_t i, n, n0[NQ], n1[NQ];

  /* Get sizes... */
  n = avk->size1;

  /* Find sub-matrices for different quantities... */
  for (iq = 0; iq < NQ; iq++) {
    n0[iq] = N;
    for (i = 0; i < n; i++) {
      if (iqa[i] == iq && n0[iq] == N)
	n0[iq] = i;
      if (iqa[i] == iq)
	n1[iq] = i - n0[iq] + 1;
    }
  }

  /* Initialize... */
  copy_atm(ctl, &atm_cont, atm, 1);
  copy_atm(ctl, &atm_res, atm, 1);

  /* Analyze quantities... */
  analyze_avk_quantity(avk, IDXP, ipa, n0, n1, atm_cont.p, atm_res.p);
  analyze_avk_quantity(avk, IDXT, ipa, n0, n1, atm_cont.t, atm_res.t);
  for (ig = 0; ig < ctl->ng; ig++)
    analyze_avk_quantity(avk, IDXQ(ig), ipa, n0, n1,
			 atm_cont.q[ig], atm_res.q[ig]);
  for (iw = 0; iw < ctl->nw; iw++)
    analyze_avk_quantity(avk, IDXK(iw), ipa, n0, n1,
			 atm_cont.k[iw], atm_res.k[iw]);
  analyze_avk_quantity(avk, IDXCLZ, ipa, n0, n1, &atm_cont.clz, &atm_res.clz);
  analyze_avk_quantity(avk, IDXCLDZ, ipa, n0, n1, &atm_cont.cldz,
		       &atm_res.cldz);
  for (icl = 0; icl < ctl->ncl; icl++)
    analyze_avk_quantity(avk, IDXCLK(icl), ipa, n0, n1,
			 &atm_cont.clk[icl], &atm_res.clk[icl]);
  analyze_avk_quantity(avk, IDXSFZ, ipa, n0, n1, &atm_cont.sfz, &atm_res.sfz);
  analyze_avk_quantity(avk, IDXSFP, ipa, n0, n1, &atm_cont.sfp, &atm_res.sfp);
  analyze_avk_quantity(avk, IDXSFT, ipa, n0, n1, &atm_cont.sft, &atm_res.sft);
  for (isf = 0; isf < ctl->nsf; isf++)
    analyze_avk_quantity(avk, IDXSFEPS(isf), ipa, n0, n1,
			 &atm_cont.sfeps[isf], &atm_res.sfeps[isf]);

  /* Write results to disk... */
  write_atm(ret->dir, "atm_cont.tab", ctl, &atm_cont);
  write_atm(ret->dir, "atm_res.tab", ctl, &atm_res);
}

/*****************************************************************************/

void analyze_avk_quantity(
  gsl_matrix * avk,
  int iq,
  int *ipa,
  size_t *n0,
  size_t *n1,
  double *cont,
  double *res) {

  /* Loop over state vector elements... */
  if (n0[iq] < N)
    for (size_t i = 0; i < n1[iq]; i++) {

      /* Get area of averaging kernel... */
      for (size_t j = 0; j < n1[iq]; j++)
	cont[ipa[n0[iq] + i]] += gsl_matrix_get(avk, n0[iq] + i, n0[iq] + j);

      /* Get information density... */
      res[ipa[n0[iq] + i]] = 1 / gsl_matrix_get(avk, n0[iq] + i, n0[iq] + i);
    }
}

/*****************************************************************************/

double cost_function(
  gsl_vector * dx,
  gsl_vector * dy,
  gsl_matrix * s_a_inv,
  gsl_vector * sig_eps_inv) {

  gsl_vector *x_aux, *y_aux;

  double chisq_a, chisq_m = 0;

  /* Get sizes... */
  size_t m = dy->size;
  size_t n = dx->size;

  /* Allocate... */
  x_aux = gsl_vector_alloc(n);
  y_aux = gsl_vector_alloc(m);

  /* Determine normalized cost function...
     (chi^2 = 1/m * [dy^T * S_eps^{-1} * dy + dx^T * S_a^{-1} * dx]) */
  for (size_t i = 0; i < m; i++)
    chisq_m += POW2(gsl_vector_get(dy, i) * gsl_vector_get(sig_eps_inv, i));
  gsl_blas_dgemv(CblasNoTrans, 1.0, s_a_inv, dx, 0.0, x_aux);
  gsl_blas_ddot(dx, x_aux, &chisq_a);

  /* Free... */
  gsl_vector_free(x_aux);
  gsl_vector_free(y_aux);

  /* Return cost function value... */
  return (chisq_m + chisq_a) / (double) m;
}

/*****************************************************************************/

void matrix_invert(
  gsl_matrix * a) {

  size_t diag = 1;

  /* Get size... */
  size_t n = a->size1;

  /* Check if matrix is diagonal... */
  for (size_t i = 0; i < n && diag; i++)
    for (size_t j = i + 1; j < n; j++)
      if (gsl_matrix_get(a, i, j) != 0) {
	diag = 0;
	break;
      }

  /* Quick inversion of diagonal matrix... */
  if (diag)
    for (size_t i = 0; i < n; i++)
      gsl_matrix_set(a, i, i, 1 / gsl_matrix_get(a, i, i));

  /* Matrix inversion by means of Cholesky decomposition... */
  else {
    gsl_linalg_cholesky_decomp(a);
    gsl_linalg_cholesky_invert(a);
  }
}

/*****************************************************************************/

void matrix_product(
  gsl_matrix * a,
  gsl_vector * b,
  int transpose,
  gsl_matrix * c) {

  gsl_matrix *aux;

  /* Set sizes... */
  size_t m = a->size1;
  size_t n = a->size2;

  /* Allocate... */
  aux = gsl_matrix_alloc(m, n);

  /* Compute A^T B A... */
  if (transpose == 1) {

    /* Compute B^1/2 A... */
    for (size_t i = 0; i < m; i++)
      for (size_t j = 0; j < n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_vector_get(b, i) * gsl_matrix_get(a, i, j));

    /* Compute A^T B A = (B^1/2 A)^T (B^1/2 A)... */
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, aux, aux, 0.0, c);
  }

  /* Compute A B A^T... */
  else if (transpose == 2) {

    /* Compute A B^1/2... */
    for (size_t i = 0; i < m; i++)
      for (size_t j = 0; j < n; j++)
	gsl_matrix_set(aux, i, j,
		       gsl_matrix_get(a, i, j) * gsl_vector_get(b, j));

    /* Compute A B A^T = (A B^1/2) (A B^1/2)^T... */
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, aux, aux, 0.0, c);
  }

  /* Free... */
  gsl_matrix_free(aux);
}

/*****************************************************************************/

void optimal_estimation(
  ret_t * ret,
  ctl_t * ctl,
  obs_t * obs_meas,
  obs_t * obs_i,
  atm_t * atm_apr,
  atm_t * atm_i) {

  static int ipa[N], iqa[N];

  gsl_matrix *a, *auxnm, *corr, *cov, *gain, *k_i, *s_a_inv;
  gsl_vector *b, *dx, *dy, *sig_eps_inv, *sig_formod, *sig_noise,
    *x_a, *x_i, *x_step, *y_aux, *y_i, *y_m;

  FILE *out;

  char filename[2 * LEN];

  double chisq, chisq_old, disq = 0, lmpar = 0.001;

  int icl, ig, ip, isf, it = 0, it2, iw;

  size_t i, j, m, n;

  /* ------------------------------------------------------------
     Initialize...
     ------------------------------------------------------------ */

  /* Get sizes... */
  m = obs2y(ctl, obs_meas, NULL, NULL, NULL);
  n = atm2x(ctl, atm_apr, NULL, iqa, ipa);
  if (m <= 0 || n <= 0)
    ERRMSG("Check problem definition!");

  /* Write info... */
  LOG(1, "Problem size: m= %d / n= %d "
      "(alloc= %.4g MB / stat= %.4g MB)",
      (int) m, (int) n,
      (double) (3 * m * n + 4 * n * n + 8 * m +
		8 * n) * sizeof(double) / 1024. / 1024.,
      (double) (5 * sizeof(atm_t) + 3 * sizeof(obs_t)
		+ 2 * N * sizeof(int)) / 1024. / 1024.);

  /* Allocate... */
  a = gsl_matrix_alloc(n, n);
  cov = gsl_matrix_alloc(n, n);
  k_i = gsl_matrix_alloc(m, n);
  s_a_inv = gsl_matrix_alloc(n, n);

  b = gsl_vector_alloc(n);
  dx = gsl_vector_alloc(n);
  dy = gsl_vector_alloc(m);
  sig_eps_inv = gsl_vector_alloc(m);
  sig_formod = gsl_vector_alloc(m);
  sig_noise = gsl_vector_alloc(m);
  x_a = gsl_vector_alloc(n);
  x_i = gsl_vector_alloc(n);
  x_step = gsl_vector_alloc(n);
  y_aux = gsl_vector_alloc(m);
  y_i = gsl_vector_alloc(m);
  y_m = gsl_vector_alloc(m);

  /* Set initial state... */
  copy_atm(ctl, atm_i, atm_apr, 0);
  copy_obs(ctl, obs_i, obs_meas, 0);
  formod(ctl, atm_i, obs_i);

  /* Set state vectors and observation vectors... */
  atm2x(ctl, atm_apr, x_a, NULL, NULL);
  atm2x(ctl, atm_i, x_i, NULL, NULL);
  obs2y(ctl, obs_meas, y_m, NULL, NULL);
  obs2y(ctl, obs_i, y_i, NULL, NULL);

  /* Set inverse a priori covariance S_a^-1... */
  set_cov_apr(ret, ctl, atm_apr, iqa, ipa, s_a_inv);
  write_matrix(ret->dir, "matrix_cov_apr.tab", ctl, s_a_inv,
	       atm_i, obs_i, "x", "x", "r");
  matrix_invert(s_a_inv);

  /* Get measurement errors... */
  set_cov_meas(ret, ctl, obs_meas, sig_noise, sig_formod, sig_eps_inv);

  /* Create cost function file... */
  sprintf(filename, "%s/costs.tab", ret->dir);
  if (!(out = fopen(filename, "w")))
    ERRMSG("Cannot create cost function file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = iteration number\n"
	  "# $2 = normalized cost function\n"
	  "# $3 = number of measurements\n"
	  "# $4 = number of state vector elements\n\n");

  /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
  gsl_vector_memcpy(dx, x_i);
  gsl_vector_sub(dx, x_a);
  gsl_vector_memcpy(dy, y_m);
  gsl_vector_sub(dy, y_i);

  /* Compute cost function... */
  chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

  /* Write info... */
  LOG(1, "it= %d / chi^2/m= %g", it, chisq);

  /* Write to cost function file... */
  fprintf(out, "%d %g %d %d\n", it, chisq, (int) m, (int) n);

  /* Compute initial kernel... */
  kernel(ctl, atm_i, obs_i, k_i);

  /* ------------------------------------------------------------
     Levenberg-Marquardt minimization...
     ------------------------------------------------------------ */

  /* Outer loop... */
  for (it = 1; it <= ret->conv_itmax; it++) {

    /* Store current cost function value... */
    chisq_old = chisq;

    /* Compute kernel matrix K_i... */
    if (it > 1 && it % ret->kernel_recomp == 0)
      kernel(ctl, atm_i, obs_i, k_i);

    /* Compute K_i^T * S_eps^{-1} * K_i ... */
    if (it == 1 || it % ret->kernel_recomp == 0)
      matrix_product(k_i, sig_eps_inv, 1, cov);

    /* Determine b = K_i^T * S_eps^{-1} * dy - S_a^{-1} * dx ... */
    for (i = 0; i < m; i++)
      gsl_vector_set(y_aux, i, gsl_vector_get(dy, i)
		     * POW2(gsl_vector_get(sig_eps_inv, i)));
    gsl_blas_dgemv(CblasTrans, 1.0, k_i, y_aux, 0.0, b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, s_a_inv, dx, 1.0, b);

    /* Inner loop... */
    for (it2 = 0; it2 < 20; it2++) {

      /* Compute A = (1 + lmpar) * S_a^{-1} + K_i^T * S_eps^{-1} * K_i ... */
      gsl_matrix_memcpy(a, s_a_inv);
      gsl_matrix_scale(a, 1 + lmpar);
      gsl_matrix_add(a, cov);

      /* Solve A * x_step = b by means of Cholesky decomposition... */
      gsl_linalg_cholesky_decomp(a);
      gsl_linalg_cholesky_solve(a, b, x_step);

      /* Update atmospheric state... */
      gsl_vector_add(x_i, x_step);
      copy_atm(ctl, atm_i, atm_apr, 0);
      copy_obs(ctl, obs_i, obs_meas, 0);
      x2atm(ctl, x_i, atm_i);

      /* Check atmospheric state... */
      for (ip = 0; ip < atm_i->np; ip++) {
	atm_i->p[ip] = GSL_MIN(GSL_MAX(atm_i->p[ip], 5e-7), 5e4);
	atm_i->t[ip] = GSL_MIN(GSL_MAX(atm_i->t[ip], 100), 400);
	for (ig = 0; ig < ctl->ng; ig++)
	  atm_i->q[ig][ip] = GSL_MIN(GSL_MAX(atm_i->q[ig][ip], 0), 1);
	for (iw = 0; iw < ctl->nw; iw++)
	  atm_i->k[iw][ip] = GSL_MAX(atm_i->k[iw][ip], 0);
      }
      atm_i->clz = GSL_MAX(atm_i->clz, 0);
      atm_i->cldz = GSL_MAX(atm_i->cldz, 0.1);
      for (icl = 0; icl < ctl->ncl; icl++)
	atm_i->clk[icl] = GSL_MAX(atm_i->clk[icl], 0);
      atm_i->sfz = GSL_MAX(atm_i->sfz, 0);
      atm_i->sfp = GSL_MAX(atm_i->sfp, 0);
      atm_i->sft = GSL_MIN(GSL_MAX(atm_i->sft, 100), 400);
      for (isf = 0; isf < ctl->nsf; isf++)
	atm_i->sfeps[isf] = GSL_MIN(GSL_MAX(atm_i->sfeps[isf], 0), 1);

      /* Forward calculation... */
      formod(ctl, atm_i, obs_i);
      obs2y(ctl, obs_i, y_i, NULL, NULL);

      /* Determine dx = x_i - x_a and dy = y - F(x_i) ... */
      gsl_vector_memcpy(dx, x_i);
      gsl_vector_sub(dx, x_a);
      gsl_vector_memcpy(dy, y_m);
      gsl_vector_sub(dy, y_i);

      /* Compute cost function... */
      chisq = cost_function(dx, dy, s_a_inv, sig_eps_inv);

      /* Modify Levenberg-Marquardt parameter... */
      if (chisq > chisq_old) {
	lmpar *= 10;
	gsl_vector_sub(x_i, x_step);
      } else {
	lmpar /= 10;
	break;
      }
    }

    /* Write info... */
    LOG(1, "it= %d / chi^2/m= %g", it, chisq);

    /* Write to cost function file... */
    fprintf(out, "%d %g %d %d\n", it, chisq, (int) m, (int) n);

    /* Get normalized step size in state space... */
    gsl_blas_ddot(x_step, b, &disq);
    disq /= (double) n;

    /* Convergence test... */
    if ((it == 1 || it % ret->kernel_recomp == 0) && disq < ret->conv_dmin)
      break;
  }

  /* Close cost function file... */
  fclose(out);

  /* Store results... */
  write_atm(ret->dir, "atm_final.tab", ctl, atm_i);
  write_obs(ret->dir, "obs_final.tab", ctl, obs_i);
  write_matrix(ret->dir, "matrix_kernel.tab", ctl, k_i,
	       atm_i, obs_i, "y", "x", "r");

  /* ------------------------------------------------------------
     Analysis of retrieval results...
     ------------------------------------------------------------ */

  /* Check if error analysis is requested... */
  if (ret->err_ana) {

    /* Allocate... */
    auxnm = gsl_matrix_alloc(n, m);
    corr = gsl_matrix_alloc(n, n);
    gain = gsl_matrix_alloc(n, m);

    /* Compute inverse retrieval covariance...
       cov^{-1} = S_a^{-1} + K_i^T * S_eps^{-1} * K_i */
    matrix_product(k_i, sig_eps_inv, 1, cov);
    gsl_matrix_add(cov, s_a_inv);

    /* Compute retrieval covariance... */
    matrix_invert(cov);
    write_matrix(ret->dir, "matrix_cov_ret.tab", ctl, cov,
		 atm_i, obs_i, "x", "x", "r");
    write_stddev("total", ret, ctl, atm_i, cov);

    /* Compute correlation matrix... */
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
	gsl_matrix_set(corr, i, j, gsl_matrix_get(cov, i, j)
		       / sqrt(gsl_matrix_get(cov, i, i))
		       / sqrt(gsl_matrix_get(cov, j, j)));
    write_matrix(ret->dir, "matrix_corr.tab", ctl, corr,
		 atm_i, obs_i, "x", "x", "r");

    /* Compute gain matrix...
       G = cov * K^T * S_eps^{-1} */
    for (i = 0; i < n; i++)
      for (j = 0; j < m; j++)
	gsl_matrix_set(auxnm, i, j, gsl_matrix_get(k_i, j, i)
		       * POW2(gsl_vector_get(sig_eps_inv, j)));
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, cov, auxnm, 0.0, gain);
    write_matrix(ret->dir, "matrix_gain.tab", ctl, gain,
		 atm_i, obs_i, "x", "y", "c");

    /* Compute retrieval error due to noise... */
    matrix_product(gain, sig_noise, 2, a);
    write_stddev("noise", ret, ctl, atm_i, a);

    /* Compute retrieval error  due to forward model errors... */
    matrix_product(gain, sig_formod, 2, a);
    write_stddev("formod", ret, ctl, atm_i, a);

    /* Compute averaging kernel matrix
       A = G * K ... */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, gain, k_i, 0.0, a);
    write_matrix(ret->dir, "matrix_avk.tab", ctl, a,
		 atm_i, obs_i, "x", "x", "r");

    /* Analyze averaging kernel matrix... */
    analyze_avk(ret, ctl, atm_i, iqa, ipa, a);

    /* Free... */
    gsl_matrix_free(auxnm);
    gsl_matrix_free(corr);
    gsl_matrix_free(gain);
  }

  /* ------------------------------------------------------------
     Finalize...
     ------------------------------------------------------------ */

  gsl_matrix_free(a);
  gsl_matrix_free(cov);
  gsl_matrix_free(k_i);
  gsl_matrix_free(s_a_inv);

  gsl_vector_free(b);
  gsl_vector_free(dx);
  gsl_vector_free(dy);
  gsl_vector_free(sig_eps_inv);
  gsl_vector_free(sig_formod);
  gsl_vector_free(sig_noise);
  gsl_vector_free(x_a);
  gsl_vector_free(x_i);
  gsl_vector_free(x_step);
  gsl_vector_free(y_aux);
  gsl_vector_free(y_i);
  gsl_vector_free(y_m);
}

/*****************************************************************************/

void read_ret(
  int argc,
  char *argv[],
  ctl_t * ctl,
  ret_t * ret) {

  /* Iteration control... */
  ret->kernel_recomp =
    (int) scan_ctl(argc, argv, "KERNEL_RECOMP", -1, "3", NULL);
  ret->conv_itmax = (int) scan_ctl(argc, argv, "CONV_ITMAX", -1, "30", NULL);
  ret->conv_dmin = scan_ctl(argc, argv, "CONV_DMIN", -1, "0.1", NULL);

  /* Error analysis... */
  ret->err_ana = (int) scan_ctl(argc, argv, "ERR_ANA", -1, "1", NULL);

  for (int id = 0; id < ctl->nd; id++)
    ret->err_formod[id] = scan_ctl(argc, argv, "ERR_FORMOD", id, "0", NULL);

  for (int id = 0; id < ctl->nd; id++)
    ret->err_noise[id] = scan_ctl(argc, argv, "ERR_NOISE", id, "0", NULL);

  ret->err_press = scan_ctl(argc, argv, "ERR_PRESS", -1, "0", NULL);
  ret->err_press_cz = scan_ctl(argc, argv, "ERR_PRESS_CZ", -1, "-999", NULL);
  ret->err_press_ch = scan_ctl(argc, argv, "ERR_PRESS_CH", -1, "-999", NULL);

  ret->err_temp = scan_ctl(argc, argv, "ERR_TEMP", -1, "0", NULL);
  ret->err_temp_cz = scan_ctl(argc, argv, "ERR_TEMP_CZ", -1, "-999", NULL);
  ret->err_temp_ch = scan_ctl(argc, argv, "ERR_TEMP_CH", -1, "-999", NULL);

  for (int ig = 0; ig < ctl->ng; ig++) {
    ret->err_q[ig] = scan_ctl(argc, argv, "ERR_Q", ig, "0", NULL);
    ret->err_q_cz[ig] = scan_ctl(argc, argv, "ERR_Q_CZ", ig, "-999", NULL);
    ret->err_q_ch[ig] = scan_ctl(argc, argv, "ERR_Q_CH", ig, "-999", NULL);
  }

  for (int iw = 0; iw < ctl->nw; iw++) {
    ret->err_k[iw] = scan_ctl(argc, argv, "ERR_K", iw, "0", NULL);
    ret->err_k_cz[iw] = scan_ctl(argc, argv, "ERR_K_CZ", iw, "-999", NULL);
    ret->err_k_ch[iw] = scan_ctl(argc, argv, "ERR_K_CH", iw, "-999", NULL);
  }

  ret->err_clz = scan_ctl(argc, argv, "ERR_CLZ", -1, "0", NULL);
  ret->err_cldz = scan_ctl(argc, argv, "ERR_CLDZ", -1, "0", NULL);
  for (int icl = 0; icl < ctl->ncl; icl++)
    ret->err_clk[icl] = scan_ctl(argc, argv, "ERR_CLK", icl, "0", NULL);

  ret->err_sfz = scan_ctl(argc, argv, "ERR_SFZ", -1, "0", NULL);
  ret->err_sfp = scan_ctl(argc, argv, "ERR_SFP", -1, "0", NULL);
  ret->err_sft = scan_ctl(argc, argv, "ERR_SFT", -1, "0", NULL);
  for (int isf = 0; isf < ctl->nsf; isf++)
    ret->err_sfeps[isf] = scan_ctl(argc, argv, "ERR_SFEPS", isf, "0", NULL);
}

/*****************************************************************************/

void set_cov_apr(
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  int *iqa,
  int *ipa,
  gsl_matrix * s_a) {

  gsl_vector *x_a;

  double x0[3], x1[3];

  /* Get sizes... */
  size_t n = s_a->size1;

  /* Allocate... */
  x_a = gsl_vector_alloc(n);

  /* Get sigma vector... */
  atm2x(ctl, atm, x_a, NULL, NULL);
  for (size_t i = 0; i < n; i++) {
    if (iqa[i] == IDXP)
      gsl_vector_set(x_a, i, ret->err_press / 100 * gsl_vector_get(x_a, i));
    if (iqa[i] == IDXT)
      gsl_vector_set(x_a, i, ret->err_temp);
    for (int ig = 0; ig < ctl->ng; ig++)
      if (iqa[i] == IDXQ(ig))
	gsl_vector_set(x_a, i, ret->err_q[ig] / 100 * gsl_vector_get(x_a, i));
    for (int iw = 0; iw < ctl->nw; iw++)
      if (iqa[i] == IDXK(iw))
	gsl_vector_set(x_a, i, ret->err_k[iw]);
    if (iqa[i] == IDXCLZ)
      gsl_vector_set(x_a, i, ret->err_clz);
    if (iqa[i] == IDXCLDZ)
      gsl_vector_set(x_a, i, ret->err_cldz);
    for (int icl = 0; icl < ctl->ncl; icl++)
      if (iqa[i] == IDXCLK(icl))
	gsl_vector_set(x_a, i, ret->err_clk[icl]);
    if (iqa[i] == IDXSFZ)
      gsl_vector_set(x_a, i, ret->err_sfz);
    if (iqa[i] == IDXSFP)
      gsl_vector_set(x_a, i, ret->err_sfp);
    if (iqa[i] == IDXSFT)
      gsl_vector_set(x_a, i, ret->err_sft);
    for (int isf = 0; isf < ctl->nsf; isf++)
      if (iqa[i] == IDXSFEPS(isf))
	gsl_vector_set(x_a, i, ret->err_sfeps[isf]);
  }

  /* Check standard deviations... */
  for (size_t i = 0; i < n; i++)
    if (POW2(gsl_vector_get(x_a, i)) <= 0)
      ERRMSG("Check a priori data (zero standard deviation)!");

  /* Initialize diagonal covariance... */
  gsl_matrix_set_zero(s_a);
  for (size_t i = 0; i < n; i++)
    gsl_matrix_set(s_a, i, i, POW2(gsl_vector_get(x_a, i)));

  /* Loop over matrix elements... */
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      if (i != j && iqa[i] == iqa[j]) {

	/* Initialize... */
	double cz = 0;
	double ch = 0;

	/* Set correlation lengths for pressure... */
	if (iqa[i] == IDXP) {
	  cz = ret->err_press_cz;
	  ch = ret->err_press_ch;
	}

	/* Set correlation lengths for temperature... */
	if (iqa[i] == IDXT) {
	  cz = ret->err_temp_cz;
	  ch = ret->err_temp_ch;
	}

	/* Set correlation lengths for volume mixing ratios... */
	for (int ig = 0; ig < ctl->ng; ig++)
	  if (iqa[i] == IDXQ(ig)) {
	    cz = ret->err_q_cz[ig];
	    ch = ret->err_q_ch[ig];
	  }

	/* Set correlation lengths for extinction... */
	for (int iw = 0; iw < ctl->nw; iw++)
	  if (iqa[i] == IDXK(iw)) {
	    cz = ret->err_k_cz[iw];
	    ch = ret->err_k_ch[iw];
	  }

	/* Compute correlations... */
	if (cz > 0 && ch > 0) {

	  /* Get Cartesian coordinates... */
	  geo2cart(0, atm->lon[ipa[i]], atm->lat[ipa[i]], x0);
	  geo2cart(0, atm->lon[ipa[j]], atm->lat[ipa[j]], x1);

	  /* Compute correlations... */
	  double rho =
	    exp(-DIST(x0, x1) / ch -
		fabs(atm->z[ipa[i]] - atm->z[ipa[j]]) / cz);

	  /* Set covariance... */
	  gsl_matrix_set(s_a, i, j, gsl_vector_get(x_a, i)
			 * gsl_vector_get(x_a, j) * rho);
	}
      }

  /* Free... */
  gsl_vector_free(x_a);
}

/*****************************************************************************/

void set_cov_meas(
  ret_t * ret,
  ctl_t * ctl,
  obs_t * obs,
  gsl_vector * sig_noise,
  gsl_vector * sig_formod,
  gsl_vector * sig_eps_inv) {

  static obs_t obs_err;

  /* Get size... */
  size_t m = sig_eps_inv->size;

  /* Noise error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for (int ir = 0; ir < obs_err.nr; ir++)
    for (int id = 0; id < ctl->nd; id++)
      obs_err.rad[id][ir]
	= (gsl_finite(obs->rad[id][ir]) ? ret->err_noise[id] : GSL_NAN);
  obs2y(ctl, &obs_err, sig_noise, NULL, NULL);

  /* Forward model error (always considered in retrieval fit)... */
  copy_obs(ctl, &obs_err, obs, 1);
  for (int ir = 0; ir < obs_err.nr; ir++)
    for (int id = 0; id < ctl->nd; id++)
      obs_err.rad[id][ir]
	= fabs(ret->err_formod[id] / 100 * obs->rad[id][ir]);
  obs2y(ctl, &obs_err, sig_formod, NULL, NULL);

  /* Total error... */
  for (size_t i = 0; i < m; i++)
    gsl_vector_set(sig_eps_inv, i, 1 / sqrt(POW2(gsl_vector_get(sig_noise, i))
					    +
					    POW2(gsl_vector_get
						 (sig_formod, i))));

  /* Check standard deviations... */
  for (size_t i = 0; i < m; i++)
    if (gsl_vector_get(sig_eps_inv, i) <= 0)
      ERRMSG("Check measurement errors (zero standard deviation)!");
}

/*****************************************************************************/

void write_stddev(
  const char *quantity,
  ret_t * ret,
  ctl_t * ctl,
  atm_t * atm,
  gsl_matrix * s) {

  static atm_t atm_aux;

  gsl_vector *x_aux;

  char filename[LEN];

  /* Get sizes... */
  size_t n = s->size1;

  /* Allocate... */
  x_aux = gsl_vector_alloc(n);

  /* Compute standard deviation... */
  for (size_t i = 0; i < n; i++)
    gsl_vector_set(x_aux, i, sqrt(gsl_matrix_get(s, i, i)));

  /* Write to disk... */
  copy_atm(ctl, &atm_aux, atm, 1);
  x2atm(ctl, x_aux, &atm_aux);
  sprintf(filename, "atm_err_%s.tab", quantity);
  write_atm(ret->dir, filename, ctl, &atm_aux);

  /* Free... */
  gsl_vector_free(x_aux);
}
