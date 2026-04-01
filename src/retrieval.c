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

static const char *ret_input_target(
  const ret_t *ret,
  const char *shared_file,
  const char *legacy_file,
  const char **dirname,
  int *profile) {

  if (shared_file[0] != '-') {
    *dirname = NULL;
    *profile = ret->profile;
    return shared_file;
  }

  *dirname = ret->dir;
  *profile = 0;
  return legacy_file;
}

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

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <dirlist>");

  /* Measure CPU-time... */
  TIMER("total", 1);

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  read_ret(argc, argv, &ctl, &ret);

  /* Initialize look-up tables... */
  tbl_t *tbl = read_tbl(&ctl);

  /* Open directory list... */
  dirlist = NULL;
  if (argv[2][0] != '-') {
    if (!(dirlist = fopen(argv[2], "r")))
      ERRMSG("Cannot open directory list!");
  }
  if (ret.proflist[0] != '-')
    if (!(proflist = fopen(ret.proflist, "r")))
      ERRMSG("Cannot open profile list!");
  if (dirlist == NULL && proflist == NULL)
    ERRMSG("Give a directory list or set PROFLIST!");

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
      have_profile = fscanf(proflist, "%d", &ret.profile) == 1;
    else {
      ret.profile = 0;
      have_profile = 1;
    }

    if ((dirlist != NULL && !have_dir) || (proflist != NULL && !have_profile)) {
      if ((dirlist != NULL && have_dir) || (proflist != NULL && have_profile))
	ERRMSG("DIRLIST and PROFLIST have different lengths!");
      break;
    }

    /* Distribute directories with MPI (optional)... */
    if ((++ntask) % size != rank)
      continue;

    /* Write info... */
    if (size > 1 && proflist != NULL && dirlist != NULL) {
      LOG(1,
	  "\nRetrieve profile %d in directory %s on rank %d of %d...\n",
	  ret.profile, ret.dir, rank + 1, size);
    } else if (size > 1 && proflist != NULL) {
      LOG(1, "\nRetrieve profile %d on rank %d of %d...\n",
	  ret.profile, rank + 1, size);
    } else if (size > 1) {
      LOG(1, "\nRetrieve in directory %s on rank %d of %d...\n",
	  ret.dir, rank + 1, size);
    } else if (proflist != NULL && dirlist != NULL) {
      LOG(1, "\nRetrieve profile %d in directory %s...\n",
	  ret.profile, ret.dir);
    } else if (proflist != NULL) {
      LOG(1, "\nRetrieve profile %d...\n", ret.profile);
    } else
      LOG(1, "\nRetrieve in directory %s...\n", ret.dir);

    /* Read atmospheric data... */
    const char *dirname;
    int profile;
    const char *filename =
      ret_input_target(&ret, ret.atm_apr_file, "atm_apr.tab",
		       &dirname, &profile);
    read_atm(dirname, filename, &ctl, &atm_apr, profile);

    /* Read observation data... */
    filename = ret_input_target(&ret, ret.obs_meas_file, "obs_meas.tab",
				&dirname, &profile);
    read_obs(dirname, filename, &ctl, &obs_meas, profile);

    /* Run retrieval... */
    double chisq;
    optimal_estimation(&ret, &ctl, tbl, &obs_meas, &obs_i, &atm_apr, &atm_i,
		       &chisq);

    /* Measure CPU-time... */
    TIMER("total", 2);
  }

  /* Write info... */
  LOG(1, "\nRetrieval done...");

  /* Measure CPU-time... */
  TIMER("total", 3);

  /* Free... */
  if (dirlist != NULL)
    fclose(dirlist);
  if (proflist != NULL)
    fclose(proflist);
  tbl_free(&ctl, tbl);

#ifdef MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
