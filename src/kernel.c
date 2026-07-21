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
  Calculate kernel functions.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Perform kernel calculations in a single directory. */
void call_kernel(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile);

/*! Read DIRLIST cases in bounded tiles and write one kernel matrix per case. */
void call_kernel_batch(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const char *dirlist,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile,
  int batch_size);

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

  char dirlist[LEN];

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 5)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: kernel <ctl> <obs> <atm> <kernel> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  SELECT_TIMER("READ_CTL", "INPUT");
  read_ctl(argc, argv, &ctl);

  /* Initialize look-up tables... */
  SELECT_TIMER("READ_TBL", "INPUT");
  tbl_t *tbl = read_tbl(&ctl);

  /* Get dirlist... */
  SELECT_TIMER("READ_DIRLIST", "INPUT");
  scan_ctl(argc, argv, "DIRLIST", -1, "-", dirlist);
  const int kernel_batch_size =
    (int) scan_ctl(argc, argv, "KERNEL_BATCH", -1, "0", NULL);
  if (kernel_batch_size < 0)
    ERRMSG("KERNEL_BATCH must not be negative!");
  /* Keep the established scalar and sequential DIRLIST behavior by default. */
  if (dirlist[0] == '-' && kernel_batch_size > 0)
    ERRMSG("KERNEL_BATCH requires DIRLIST!");

  /* Set flags... */
  ctl.write_matrix = 1;

  /* Single kernel calculation... */
  if (dirlist[0] == '-')
    call_kernel(&ctl, tbl, NULL, argv[2], argv[3], argv[4]);

  /* Work on directory list... */
  /* The batch path is meaningful only when independent cases come from DIRLIST. */
  else if (kernel_batch_size > 0)
    call_kernel_batch(&ctl, tbl, dirlist, argv[2], argv[3], argv[4],
                      kernel_batch_size);
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
      call_kernel(&ctl, tbl, wrkdir, argv[2], argv[3], argv[4]);
    }

    /* Close dirlist... */
    fclose(in);
  }

  /* Free... */
  SELECT_TIMER("FINALIZE", "OVERHEAD");
  tbl_free(&ctl, tbl);
  PRINT_TIMERS;

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void usage(
  void) {

  printf("\nJURASSIC kernel tool.\n\n");
  printf("Compute Jacobian or kernel matrices for the configured channels\n");
  printf
    ("from observation geometry, atmospheric state, and control settings.\n\n");
  printf("Usage:\n");
  printf("  kernel <ctl> <obs> <atm> <kernel> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>     Control file.\n");
  printf("  <obs>     Observation geometry input file.\n");
  printf("  <atm>     Atmospheric state input file.\n");
  printf("  <kernel>  Output file for the kernel matrix.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Tool-specific control parameters:\n");
  printf
    ("  DIRLIST <file>   Read working directories from <file> and run one case\n");
  printf
    ("                   per directory using the same <obs>, <atm>, and\n");
  printf
    ("                   <kernel> filenames relative to each listed directory.\n");
  printf("\n");
  printf("  KERNEL_BATCH <n>  Process DIRLIST cases in batches of <n> retrievals.\n");
  printf("Common control parameters:\n");
  printf
    ("  TBLBASE, TBLFMT               Lookup-table base name and format.\n");
  printf("  ATMFMT, OBSFMT, MATRIXFMT     Input/output file formats.\n");
  printf("  NG, EMITTER[i]                Active emitters.\n");
  printf("  ND, NU[i], NW, WINDOW[i]      Spectral channels and windows.\n");
  printf
    ("  NCL, CLNU[i], NSF, SFNU[i]    Cloud and surface spectral grids.\n");
  printf("  RET*_ZMIN, RET*_ZMAX          State-vector altitude limits.\n");
  printf
    ("  WRITE_BBT, FORMOD             Output units and forward-model selection.\n");
  printf("  CTM_*, REFRAC                 Continua and refractivity.\n");
  printf
    ("  RAYDS, RAYDZ, FOV             Ray tracing and field of view.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}

/*****************************************************************************/

void call_kernel(
  const ctl_t *ctl,
  const tbl_t *tbl,
  const char *wrkdir,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile) {

  static atm_t atm;
  static obs_t obs;

  /* Read observation geometry... */
  SELECT_TIMER("READ_OBS", "INPUT");
  read_obs(wrkdir, obsfile, ctl, &obs, 0);

  /* Read atmospheric data... */
  SELECT_TIMER("READ_ATM", "INPUT");
  read_atm(wrkdir, atmfile, ctl, &atm, 0);

  /* Get sizes... */
  const size_t n = atm2x(ctl, &atm, NULL, NULL, NULL);
  const size_t m = obs2y(ctl, &obs, NULL, NULL, NULL);

  /* Check sizes... */
  if (n == 0)
    ERRMSG("No state vector elements!");
  if (m == 0)
    ERRMSG("No measurement vector elements!");

  /* Allocate... */
  gsl_matrix *k = gsl_matrix_alloc(m, n);

  /* Compute kernel matrix... */
    /* Each tile is independent, so its matrices can be released immediately. */
  SELECT_TIMER("KERNEL", "FORWARD");
  kernel(ctl, tbl, &atm, &obs, k);

  /* Write matrix to file... */
  SELECT_TIMER("WRITE_MATRIX", "OUTPUT");
  write_matrix(wrkdir, kernelfile, ctl, k, &atm, &obs, "y", "x", "r", 0);

  /* Free... */
  gsl_matrix_free(k);
}

/*******************************************************************************/

void call_kernel_batch(
  const ctl_t *ctl,
  const tbl_t *tbl,
  const char *dirlist,
  const char *obsfile,
  const char *atmfile,
  const char *kernelfile,
  int batch_size) {

  FILE *in;
  if (!(in = fopen(dirlist, "r")))
    ERRMSG("Cannot open directory list!");

  char *dirs;
  atm_t *atm;
  obs_t *obs;
  gsl_matrix **k;
  ALLOC(dirs, char, (size_t) batch_size * LEN);
  ALLOC(atm, atm_t, batch_size);
  ALLOC(obs, obs_t, batch_size);
  /* Keep only one CLI tile in host memory, regardless of DIRLIST length. */
  ALLOC(k, gsl_matrix *, batch_size);

  /* The final tile may contain fewer cases than KERNEL_BATCH. */
  while (1) {
    int nret = 0;
    size_t n = 0;
    size_t m = 0;

    /* Read one directory-list tile. */
    while (nret < batch_size) {
      char *dir = &dirs[(size_t) nret * LEN];
      if (fscanf(in, "%4999s", dir) != 1)
        break;
      LOG(1, "\nWorking directory: %s", dir);

      SELECT_TIMER("READ_OBS", "INPUT");
      read_obs(dir, obsfile, ctl, &obs[nret], 0);
      SELECT_TIMER("READ_ATM", "INPUT");
      read_atm(dir, atmfile, ctl, &atm[nret], 0);

      /* kernel_batch requires common vector sizes within the current tile. */
      const size_t n0 = atm2x(ctl, &atm[nret], NULL, NULL, NULL);
      const size_t m0 = obs2y(ctl, &obs[nret], NULL, NULL, NULL);
      if (n0 == 0)
        ERRMSG("No state vector elements!");
      if (m0 == 0)
        ERRMSG("No measurement vector elements!");
      if (nret == 0) {
        n = n0;
        m = m0;
      } else if (n0 != n || m0 != m)
        ERRMSG("Kernel batch requires matching vector sizes!");

      k[nret] = gsl_matrix_alloc(m, n);
      nret++;
    }
    if (nret == 0)
      break;

    /* Each tile is independent, so its matrices can be released immediately. */
    SELECT_TIMER("KERNEL", "FORWARD");
    kernel_batch(ctl, tbl, atm, obs, k, nret, nret);

    for (int ir = 0; ir < nret; ir++) {
      SELECT_TIMER("WRITE_MATRIX", "OUTPUT");
      /* Preserve the normal DIRLIST output convention for every case. */
      write_matrix(&dirs[(size_t) ir * LEN], kernelfile, ctl, k[ir],
                   &atm[ir], &obs[ir], "y", "x", "r", 0);
      gsl_matrix_free(k[ir]);
    }
  }

  free(k);
  free(obs);
  free(atm);
  free(dirs);
  fclose(in);
}
