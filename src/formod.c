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
  char task_mode,
  const char *obsref,
  int formod_scalar,
  int batch_size);

/*! Parse the formod task selector and map aliases to an internal task code. */
char parse_task_mode(
  const char * task);

/*! Parse the formod execution selector. */
int parse_formod_scalar(
  int argc,
  char *argv[]);

/*! Parse the benchmark batch size. */
int parse_formod_batch_size(
  int argc,
  char *argv[]);

/*! Execute a single forward model call via the default batch path or the
  scalar debug path. */
void exec_formod_single(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar);

/*! Execute a batch throughput benchmark with perturbed atmospheric cases and
  return the first result in @p obs. */
void exec_formod_batch_repeat(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const atm_t * atm,
  obs_t * obs,
  int batch_size);

/*! Calculate relative errors. */
void compute_rel_errors(
  const ctl_t * ctl,
  const obs_t * obs_test,
  const obs_t * obs_ref,
  double *mre,
  double *sdre,
  double *minre,
  double *maxre);

/*! Run the default forward-model path. */
void exec_formod_default(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar);

/*! Process observation data profile by profile using time-based matching. */
void exec_formod_profiles(
  ctl_t * ctl,
  const tbl_t * tbl,
  const char *wrkdir,
  const char *radfile,
  const atm_t * atm,
  obs_t * obs,
  atm_t * atm_scratch,
  obs_t * obs_single,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar);

/*! Print relative-error diagnostics against reference observations. */
void eval_formod_results(
  const ctl_t * ctl,
  const char *wrkdir,
  const char *obsref,
  const obs_t * obs,
  obs_t * obs_ref);

/*! Write emitter and extinction contribution outputs. */
void exec_formod_contributions(
  ctl_t * ctl,
  const tbl_t * tbl,
  const char *wrkdir,
  const char *radfile,
  const atm_t * atm,
  obs_t * obs,
  atm_t * atm_scratch,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar);

/*! Run the built-in runtime benchmark. */
void exec_formod_benchmark(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const atm_t * atm,
  obs_t * obs,
  atm_t * atm_scratch,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar,
  int batch_size);

/*! Analyze sensitivity to ray-tracing step sizes. */
void exec_formod_stepsize(
  ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs,
  obs_t * obs_ref,
  los_t * los_scratch,
  obs_t * obs_scratch,
  int formod_scalar);

/*! Print command-line help. */
void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: formod <ctl> <obs> <atm> <rad> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

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
  const int formod_scalar = parse_formod_scalar(argc, argv);
  const int batch_size = parse_formod_batch_size(argc, argv);

  /* Initialize look-up tables... */
  SELECT_TIMER("READ_TBL", "INPUT");
  tbl_t *tbl = read_tbl(&ctl);

  /* Get task... */
  SELECT_TIMER("READ_TASK", "INPUT");
  scan_ctl(argc, argv, "TASK", -1, "forward", task);
  const char task_mode = parse_task_mode(task);

  /* Get dirlist... */
  SELECT_TIMER("READ_DIRLIST", "INPUT");
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);

  /* Get reference data... */
  SELECT_TIMER("READ_OBSREF", "INPUT");
  scan_ctl(argc, argv, "OBSREF", -1, "-", obsref);

  /* Single forward calculation... */
  if (dirlist[0] == '-')
    call_formod(&ctl, tbl, NULL, argv[2], argv[3], argv[4], task_mode, obsref,
		formod_scalar, batch_size);

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
      call_formod(&ctl, tbl, wrkdir, argv[2], argv[3], argv[4], task_mode,
		  obsref,
		  formod_scalar, batch_size);
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

void usage(
  void) {

  printf("\nJURASSIC forward model.\n\n");
  printf
    ("Compute simulated radiances or transmittances for the configured\n");
  printf
    ("channels from observation geometry, atmospheric state, and control\n");
  printf("settings.\n\n");
  printf("Usage:\n");
  printf("  formod <ctl> <obs> <atm> <rad> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>  Control file.\n");
  printf("  <obs>  Observation geometry input file.\n");
  printf("  <atm>  Atmospheric state input file.\n");
  printf("  <rad>  Output file for simulated radiances or transmittances.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Tool-specific control parameters:\n");
  printf("  TASK <name>      Select the operation mode.\n");
  printf
    ("                   forward or omitted: regular forward-model simulation.\n");
  printf
    ("                   contrib: write separate outputs for individual emitter\n");
  printf("                      extinction contributions.\n");
  printf
    ("                   profiles: process observation data profile by profile,\n");
  printf
    ("                      matching each profile to the corresponding\n");
  printf("                      atmospheric profile by time.\n");
  printf
    ("                   stepsize: analyze sensitivity to ray-tracing step\n");
  printf("                      sizes.\n");
  printf
    ("                   time: run the built-in runtime benchmark using\n");
  printf
    ("                      perturbed versions of the input atmosphere.\n");
  printf
    ("  DIRLIST <file>   Read working directories from <file> and run one case\n");
  printf
    ("                   per directory using the same <obs>, <atm>, and <rad>\n");
  printf("                   filenames relative to each listed directory.\n");
  printf
    ("  OBSREF <file>    Read reference observations from <file> and print\n");
  printf
    ("                   relative-error diagnostics against the simulated output.\n");
  printf
    ("  EXECUTION <mode> Select execution backend: batch (default) or scalar.\n");
  printf
    ("  BATCH_SIZE <n>   Batch size for TASK=time with EXECUTION=batch.\n");
  printf("\n");
  printf("Common control parameters:\n");
  printf
    ("  TBLBASE, TBLFMT               Lookup-table base name and format.\n");
  printf
    ("  ATMFMT, OBSFMT                Atmospheric and observation file formats.\n");
  printf("  NG, EMITTER[i]                Active emitters.\n");
  printf("  ND, NU[i], NW, WINDOW[i]      Spectral channels and windows.\n");
  printf
    ("  NCL, CLNU[i], NSF, SFNU[i]    Cloud and surface spectral grids.\n");
  printf
    ("  WRITE_BBT, FORMOD             Output units and forward-model selection.\n");
  printf("  CTM_*, REFRAC                 Continua and refractivity.\n");
  printf
    ("  RAYDS, RAYDZ, FOV             Ray tracing and field of view.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}

/*****************************************************************************/

char parse_task_mode(
  const char *task) {

  if (task[0] == '-' || strcmp(task, "forward") == 0)
    return 'f';
  if (task[0] == 'c' || strcmp(task, "contrib") == 0
      || strcmp(task, "contributions") == 0)
    return 'c';
  if (task[0] == 'p' || strcmp(task, "profile") == 0
      || strcmp(task, "profiles") == 0)
    return 'p';
  if (task[0] == 's' || strcmp(task, "stepsize") == 0)
    return 's';
  if (task[0] == 't' || strcmp(task, "time") == 0
      || strcmp(task, "benchmark") == 0)
    return 't';
  if (task[0] == 'b' || strcmp(task, "batch") == 0
      || strcmp(task, "batch_benchmark") == 0)
    return 't';

  ERRMSG("Unknown TASK '%s'!", task);
  return 'f';
}

/*****************************************************************************/

int parse_formod_scalar(
  int argc,
  char *argv[]) {

  char execution[LEN];
  scan_ctl(argc, argv, "EXECUTION", -1, "batch", execution);

  if (strcmp(execution, "batch") == 0)
    return 0;
  if (strcmp(execution, "scalar") == 0)
    return 1;

  ERRMSG("Unknown EXECUTION '%s'!", execution);
  return 0;
}

/*****************************************************************************/

int parse_formod_batch_size(
  int argc,
  char *argv[]) {

  const int batch_size =
    (int) scan_ctl(argc, argv, "BATCH_SIZE", -1, "1", NULL);

  if (batch_size < 1)
    ERRMSG("BATCH_SIZE must be positive!");

  return batch_size;
}

/*****************************************************************************/

void exec_formod_single(
  const ctl_t *ctl,
  const tbl_t *tbl,
  atm_t *atm,
  obs_t *obs,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar) {

  if (formod_scalar) {
    const int status = formod(ctl, tbl, atm, obs, los_scratch, obs_scratch);
    if (status != FORMOD_STATUS_OK)
      ERRMSG("Forward model failed with status code %d!", status);
    return;
  }

  int status[1];
#if defined(_OPENACC)
#pragma acc enter data create(atm[0:1],obs[0:1],status[0:1],los_scratch[0:1],   obs_scratch[0:1])
#endif
  formod_batch(ctl, tbl, atm, obs, 1, status, los_scratch, obs_scratch);
#if defined(_OPENACC)
#pragma acc exit data delete(atm[0:1],obs[0:1],status[0:1],los_scratch[0:1],   obs_scratch[0:1])
#endif
  if (status[0] != FORMOD_STATUS_OK)
    ERRMSG("Forward model failed with status code %d!", status[0]);
}

/*****************************************************************************/

void exec_formod_batch_repeat(
  const ctl_t *ctl,
  const tbl_t *tbl,
  const atm_t *atm,
  obs_t *obs,
  int batch_size) {

  atm_t *atm_batch;
  obs_t *obs_batch;
  los_t *los_batch;
  obs_t *obs_scratch_batch;
  int *status;
  gsl_rng *rng;

  if (batch_size < 1)
    ERRMSG("BATCH_SIZE must be positive!");

  ALLOC(atm_batch, atm_t, batch_size);
  ALLOC(obs_batch, obs_t, batch_size);
  ALLOC(los_batch, los_t, batch_size);
  ALLOC(obs_scratch_batch, obs_t, batch_size);
  ALLOC(status, int,
	batch_size);

  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, 0UL);

  for (int ib = 0; ib < batch_size; ib++) {
    atm_batch[ib] = *atm;
    obs_batch[ib] = *obs;

    if (ib == 0)
      continue;

    const double dtemp = 40.0 * (gsl_rng_uniform(rng) - 0.5);
    const double dpress = 1.0 - 0.1 * gsl_rng_uniform(rng);
    double dq[NG];
    for (int ig = 0; ig < ctl->ng; ig++)
      dq[ig] = 0.8 + 0.4 * gsl_rng_uniform(rng);
    for (int ip = 0; ip < atm_batch[ib].np; ip++) {
      atm_batch[ib].t[ip] += dtemp;
      atm_batch[ib].p[ip] *= dpress;
      for (int ig = 0; ig < ctl->ng; ig++)
	atm_batch[ib].q[ig][ip] *= dq[ig];
    }
  }

#if defined(_OPENACC)
#pragma acc enter data create(atm_batch[0:batch_size],obs_batch[0:batch_size],status[0:batch_size],los_batch[0:batch_size],obs_scratch_batch[0:batch_size])
#endif
  formod_batch(ctl, tbl, atm_batch, obs_batch, batch_size, status,
	       los_batch, obs_scratch_batch);
#if defined(_OPENACC)
#pragma acc exit data delete(atm_batch[0:batch_size],obs_batch[0:batch_size],status[0:batch_size],los_batch[0:batch_size],obs_scratch_batch[0:batch_size])
#endif

  if (status[0] != FORMOD_STATUS_OK)
    ERRMSG("Forward model failed with status code %d!", status[0]);
  for (int ib = 1; ib < batch_size; ib++)
    if (status[ib] != FORMOD_STATUS_OK)
      ERRMSG("Batch benchmark failed with status code %d at element %d!",
	     status[ib], ib);

  *obs = obs_batch[0];

  gsl_rng_free(rng);
  free(status);
  free(obs_scratch_batch);
  free(los_batch);
  free(obs_batch);
  free(atm_batch);
}

/*****************************************************************************/

void exec_formod_default(
  const ctl_t *ctl,
  const tbl_t *tbl,
  atm_t *atm,
  obs_t *obs,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar) {

  exec_formod_single(ctl, tbl, atm, obs, los_scratch, obs_scratch,
		     formod_scalar);
}

/*****************************************************************************/

void exec_formod_profiles(
  ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *radfile,
  const atm_t *atm,
  obs_t *obs,
  atm_t *atm_scratch,
  obs_t *obs_single,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar) {

  SELECT_TIMER("FORMOD_PROFILE_LOOP", "FORWARD");

  for (int ir = 0; ir < obs->nr; ir++) {

    atm_scratch->np = 0;
    for (int ip = 0; ip < atm->np; ip++)
      if (atm->time[ip] == obs->time[ir]) {
	atm_scratch->time[atm_scratch->np] = atm->time[ip];
	atm_scratch->z[atm_scratch->np] = atm->z[ip];
	atm_scratch->lon[atm_scratch->np] = atm->lon[ip];
	atm_scratch->lat[atm_scratch->np] = atm->lat[ip];
	atm_scratch->p[atm_scratch->np] = atm->p[ip];
	atm_scratch->t[atm_scratch->np] = atm->t[ip];
	for (int ig = 0; ig < ctl->ng; ig++)
	  atm_scratch->q[ig][atm_scratch->np] = atm->q[ig][ip];
	for (int iw = 0; iw < ctl->nw; iw++)
	  atm_scratch->k[iw][atm_scratch->np] = atm->k[iw][ip];
	atm_scratch->np++;
      }

    obs_single->nr = 1;
    obs_single->time[0] = obs->time[ir];
    obs_single->vpz[0] = obs->vpz[ir];
    obs_single->vplon[0] = obs->vplon[ir];
    obs_single->vplat[0] = obs->vplat[ir];
    obs_single->obsz[0] = obs->obsz[ir];
    obs_single->obslon[0] = obs->obslon[ir];
    obs_single->obslat[0] = obs->obslat[ir];

    if (atm_scratch->np > 0) {
      exec_formod_single(ctl, tbl, atm_scratch, obs_single, los_scratch,
			 obs_scratch, formod_scalar);

      for (int id = 0; id < ctl->nd; id++) {
	obs->rad[id][ir] = obs_single->rad[id][0];
	obs->tau[id][ir] = obs_single->tau[id][0];
      }
    }
  }

  SELECT_TIMER("WRITE_OBS", "OUTPUT");
  write_obs(wrkdir, radfile, ctl, obs, 0);
}

/*****************************************************************************/

void eval_formod_results(
  const ctl_t *ctl,
  const char *wrkdir,
  const char *obsref,
  const obs_t *obs,
  obs_t *obs_ref) {

  SELECT_TIMER("EVAL", "OUTPUT");

  read_obs(wrkdir, obsref, ctl, obs_ref, 0);

  double mre[ND], sdre[ND], minre[ND], maxre[ND];
  compute_rel_errors(ctl, obs, obs_ref, mre, sdre, minre, maxre);

  for (int id = 0; id < ctl->nd; id++)
    printf
      ("EVAL: nu= %.4f cm^-1 | MRE= %g %% | SDRE= %g %% | MinRE= %g %% | MaxRE= %g %%\n",
	 ctl->nu[id], mre[id], sdre[id], minre[id], maxre[id]);
}

/*****************************************************************************/

void exec_formod_contributions(
  ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *radfile,
  const atm_t *atm,
  obs_t *obs,
  atm_t *atm_scratch,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar) {

  char filename[LEN];

  SELECT_TIMER("CONTRIBUTIONS", "FORWARD");

  ctl->ctm_co2 = 0;
  ctl->ctm_h2o = 0;
  ctl->ctm_n2 = 0;
  ctl->ctm_o2 = 0;

  for (int ig = 0; ig < ctl->ng; ig++) {
    copy_atm(ctl, atm_scratch, atm, 0);

    for (int iw = 0; iw < ctl->nw; iw++)
      for (int ip = 0; ip < atm_scratch->np; ip++)
	atm_scratch->k[iw][ip] = 0;

    for (int ig2 = 0; ig2 < ctl->ng; ig2++)
      if (ig2 != ig)
	for (int ip = 0; ip < atm_scratch->np; ip++)
	  atm_scratch->q[ig2][ip] = 0;

    exec_formod_single(ctl, tbl, atm_scratch, obs, los_scratch,
		       obs_scratch, formod_scalar);

    sprintf(filename, "%s.%s", radfile, ctl->emitter[ig]);
    write_obs(wrkdir, filename, ctl, obs, 0);
  }

  copy_atm(ctl, atm_scratch, atm, 0);

  for (int ig = 0; ig < ctl->ng; ig++)
    for (int ip = 0; ip < atm_scratch->np; ip++)
      atm_scratch->q[ig][ip] = 0;

  exec_formod_single(ctl, tbl, atm_scratch, obs, los_scratch, obs_scratch,
		     formod_scalar);

  sprintf(filename, "%s.EXTINCT", radfile);
  write_obs(wrkdir, filename, ctl, obs, 0);
}

/*****************************************************************************/

void exec_formod_benchmark(
  const ctl_t *ctl,
  const tbl_t *tbl,
  const atm_t *atm,
  obs_t *obs,
  atm_t *atm_scratch,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar,
  int batch_size) {

  double t_min = 0, t_max = 0, t_mean = 0, t_sd = 0;
  int n = 0;

  gsl_rng_env_setup();
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  if (formod_scalar && batch_size > 1)
    ERRMSG("BATCH_SIZE > 1 cannot be combined with EXECUTION scalar!");

  do {
    SELECT_TIMER("BENCHMARK_SAMPLE", "ANALYSIS");
    double t0 = omp_get_wtime();

    if (batch_size > 1)
      exec_formod_batch_repeat(ctl, tbl, atm, obs, batch_size);
    else {
      copy_atm(ctl, atm_scratch, atm, 0);
      double dtemp = 40. * (gsl_rng_uniform(rng) - 0.5);
      double dpress = 1. - 0.1 * gsl_rng_uniform(rng);
      double dq[NG];
      for (int ig = 0; ig < ctl->ng; ig++)
	dq[ig] = 0.8 + 0.4 * gsl_rng_uniform(rng);
      for (int ip = 0; ip < atm_scratch->np; ip++) {
	atm_scratch->t[ip] += dtemp;
	atm_scratch->p[ip] *= dpress;
	for (int ig = 0; ig < ctl->ng; ig++)
	  atm_scratch->q[ig][ip] *= dq[ig];
      }

      exec_formod_single(ctl, tbl, atm_scratch, obs, los_scratch,
			 obs_scratch, formod_scalar);
    }

    double dt = omp_get_wtime() - t0;

    t_mean += dt;
    t_sd += POW2(dt);
    if (n == 0 || dt < t_min)
      t_min = dt;
    if (n == 0 || dt > t_max)
      t_max = dt;
    n++;

  } while (t_mean < 10.0);

  t_mean /= (double) n;
  t_sd = sqrt(t_sd / (double) n - POW2(t_mean));
  printf("RUNTIME: execution= %s | batch_size= %d | mean= %g s"
	 " | stddev= %g s | min= %g s | max= %g s\n",
	 formod_scalar ? "scalar" : "batch", batch_size, t_mean, t_sd, t_min,
	 t_max);

  gsl_rng_free(rng);
}

/*****************************************************************************/

void exec_formod_stepsize(
  ctl_t *ctl,
  const tbl_t *tbl,
  atm_t *atm,
  obs_t *obs,
  obs_t *obs_ref,
  los_t *los_scratch,
  obs_t *obs_scratch,
  int formod_scalar) {

  SELECT_TIMER("STEP_ANALYSIS", "ANALYSIS");

  ctl->rayds = 0.1;
  ctl->raydz = 0.01;
  exec_formod_single(ctl, tbl, atm, obs, los_scratch, obs_scratch,
		     formod_scalar);
  copy_obs(ctl, obs_ref, obs, 0);

  for (double dz = 0.01; dz <= 2; dz *= 1.1)
    for (double ds = 0.1; ds <= 50; ds *= 1.1) {
      ctl->rayds = ds;
      ctl->raydz = dz;

      double t0 = omp_get_wtime();
      exec_formod_single(ctl, tbl, atm, obs, los_scratch,
			 obs_scratch, formod_scalar);
      double dt = omp_get_wtime() - t0;

      double mre[ND], sdre[ND], minre[ND], maxre[ND];
      compute_rel_errors(ctl, obs, obs_ref, mre, sdre, minre, maxre);

      for (int id = 0; id < ctl->nd; id++)
	printf
	  ("STEPSIZE: ds= %.4f km | dz= %g km | t= %g s | nu= %.4f cm^-1"
	   " | MRE= %g %% | SDRE= %g %% | MinRE= %g %% | MaxRE= %g %%\n",
	   ds, dz, dt, ctl->nu[id], mre[id], sdre[id], minre[id],
	   maxre[id]);
    }
}

/*****************************************************************************/

void call_formod(
  ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *radfile,
  char task_mode,
  const char *obsref,
  int formod_scalar,
  int batch_size) {

  static atm_t atm, atm2;
  static obs_t obs, obs2, obs_ref, obs_scratch;
  static los_t los_scratch;

  /* Read atmospheric data... */
  SELECT_TIMER("READ_ATM", "INPUT");
  read_atm(wrkdir, atmfile, ctl, &atm, 0);

  /* Read observation geometry... */
  SELECT_TIMER("READ_OBS", "INPUT");
  read_obs(wrkdir, obsfile, ctl, &obs, 0);

  /* Compute multiple profiles... */
  if (task_mode == 'p') {
    exec_formod_profiles(ctl, tbl, wrkdir, radfile, &atm, &obs, &atm2, &obs2,
			 &los_scratch, &obs_scratch, formod_scalar);
  }

  /* Compute single profile... */
  else {

    if (task_mode != 't' && batch_size > 1)
      ERRMSG("BATCH_SIZE > 1 is only supported for TASK=time!");

    SELECT_TIMER("FORMOD", "FORWARD");
    exec_formod_default(ctl, tbl, &atm, &obs, &los_scratch, &obs_scratch,
			formod_scalar);

    SELECT_TIMER("WRITE_OBS", "OUTPUT");
    write_obs(wrkdir, radfile, ctl, &obs, 0);

    if (obsref[0] != '-') {
      eval_formod_results(ctl, wrkdir, obsref, &obs, &obs_ref);
    }

    if (task_mode == 'c') {
      exec_formod_contributions(ctl, tbl, wrkdir, radfile, &atm, &obs, &atm2,
				   &los_scratch, &obs_scratch, formod_scalar);
    }

    if (task_mode == 't') {
      exec_formod_benchmark(ctl, tbl, &atm, &obs, &atm2, &los_scratch,
			    &obs_scratch, formod_scalar, batch_size);
    }

    if (task_mode == 's') {
      exec_formod_stepsize(ctl, tbl, &atm, &obs, &obs_ref, &los_scratch,
			   &obs_scratch, formod_scalar);
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
