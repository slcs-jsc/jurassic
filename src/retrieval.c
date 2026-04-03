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
  JURASSIC retrieval processor.
*/

#include "jurassic.h"

#ifdef MPI
#include <mpi.h>
#endif

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
void usage(void);

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
  FILE *proflist = NULL;

  /* MPI task distribution (optional)... */
  int ntask = -1;
  int rank = 0;
  int size = 1;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  if (argc == 2
      && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: retrieval <ctl> <dirlist> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  SELECT_TIMER("READ_CTL", "INPUT");
  read_ctl(argc, argv, &ctl);
  SELECT_TIMER("READ_RET", "INPUT");
  read_ret(argc, argv, &ctl, &ret);

  /* Initialize look-up tables... */
  SELECT_TIMER("READ_TBL", "INPUT");
  tbl_t *tbl = read_tbl(&ctl);

  /* Open directory list... */
  SELECT_TIMER("OPEN_CASELISTS", "INPUT");
  dirlist = NULL;
  if (argv[2][0] != '-') {
    if (!(dirlist = fopen(argv[2], "r")))
      ERRMSG("Cannot open directory list!");
  }
  if (ret.shared_io_proflist[0] != '-')
    if (!(proflist = fopen(ret.shared_io_proflist, "r")))
      ERRMSG("Cannot open profile list!");
  if (dirlist == NULL && proflist == NULL)
    ERRMSG("Give a directory list or set SHARED_IO_PROFLIST!");

  /* Loop over retrieval cases... */
  while (1) {

    int have_dir = 0, have_profile = 0;

    if (dirlist != NULL)
      have_dir = fscanf(dirlist, "%4999s", ret.dir) == 1;
    else {
      sprintf(ret.dir, ".");
      have_dir = 1;
    }

    if (proflist != NULL)
      have_profile = fscanf(proflist, "%d", &ret.shared_io_profile) == 1;
    else {
      ret.shared_io_profile = 0;
      have_profile = 1;
    }

    if ((dirlist != NULL && !have_dir) || (proflist != NULL && !have_profile)) {
      if ((dirlist != NULL && have_dir) || (proflist != NULL && have_profile))
	ERRMSG("DIRLIST and SHARED_IO_PROFLIST have different lengths!");
      break;
    }

    /* Distribute directories with MPI (optional)... */
    if ((++ntask) % size != rank)
      continue;

    /* Write info... */
    if (size > 1 && proflist != NULL && dirlist != NULL) {
      LOG(1,
	  "\nRetrieve profile %d in directory %s on rank %d of %d...\n",
	  ret.shared_io_profile, ret.dir, rank + 1, size);
    } else if (size > 1 && proflist != NULL) {
      LOG(1, "\nRetrieve profile %d on rank %d of %d...\n",
	  ret.shared_io_profile, rank + 1, size);
    } else if (size > 1) {
      LOG(1, "\nRetrieve in directory %s on rank %d of %d...\n",
	  ret.dir, rank + 1, size);
    } else if (proflist != NULL && dirlist != NULL) {
      LOG(1, "\nRetrieve profile %d in directory %s...\n",
	  ret.shared_io_profile, ret.dir);
    } else if (proflist != NULL) {
      LOG(1, "\nRetrieve profile %d...\n", ret.shared_io_profile);
    } else
      LOG(1, "\nRetrieve in directory %s...\n", ret.dir);

    /* Read atmospheric data... */
    SELECT_TIMER("READ_ATM", "INPUT");
    const char *dirname;
    int profile;
    const char *filename =
      shared_io_input_target(&ret, ret.shared_io_atm_apr_file, "atm_apr.tab",
			     &dirname, &profile);
    read_atm(dirname, filename, &ctl, &atm_apr, profile);

    /* Read observation data... */
    SELECT_TIMER("READ_OBS", "INPUT");
    filename =
      shared_io_input_target(&ret, ret.shared_io_obs_meas_file, "obs_meas.tab",
			     &dirname, &profile);
    read_obs(dirname, filename, &ctl, &obs_meas, profile);

    /* Run retrieval... */
    double chisq;
    optimal_estimation(&ret, &ctl, tbl, &obs_meas, &obs_i, &atm_apr, &atm_i,
		       &chisq);
  }

  /* Write info... */
  SELECT_TIMER("FINALIZE", "OVERHEAD");
  LOG(1, "\nRetrieval done...");

  /* Free... */
  if (dirlist != NULL)
    fclose(dirlist);
  if (proflist != NULL)
    fclose(proflist);
  tbl_free(&ctl, tbl);

#ifdef MPI
  MPI_Finalize();
#endif

  PRINT_TIMERS;

  return EXIT_SUCCESS;
}

/*****************************************************************************/

void usage(void) {

  printf("\nJURASSIC retrieval tool.\n\n");
  printf("Run optimal-estimation retrievals from measured observations,\n");
  printf("a priori atmospheric input, and control settings.\n\n");
  printf("Usage:\n");
  printf("  retrieval <ctl> <dirlist> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>      Control file.\n");
  printf("  <dirlist>  File containing working directories to process.\n");
  printf("             Use '-' together with shared I/O settings such as\n");
  printf("             SHARED_IO_PROFLIST for shared-file workflows.\n");
  printf("  [KEY VALUE]  Optional control parameters.\n\n");
  printf("Optional control overrides:\n");
  printf("  SHARED_IO_PROFLIST <file>\n");
  printf("             Read profile indices from <file> for shared-file retrievals.\n");
  printf("  SHARED_IO_* <file>\n");
  printf("             Override shared input or output filenames on the command line.\n");
  printf("\nFurther information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
