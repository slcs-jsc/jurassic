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
  if (!(dirlist = fopen(argv[2], "r")))
    ERRMSG("Cannot open directory list!");

  /* Loop over directories... */
  while (fscanf(dirlist, "%4999s", ret.dir) != EOF) {

    /* Distribute directories with MPI (optional)... */
    if ((++ntask) % size != rank)
      continue;

    /* Write info... */
    if (size > 1) {
      LOG(1, "\nRetrieve in directory %s on rank %d of %d...\n",
	  ret.dir, rank + 1, size);
    } else
      LOG(1, "\nRetrieve in directory %s...\n", ret.dir);

    /* Read atmospheric data... */
    read_atm(ret.dir, "atm_apr.tab", &ctl, &atm_apr);

    /* Read observation data... */
    read_obs(ret.dir, "obs_meas.tab", &ctl, &obs_meas);

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
  tbl_free(&ctl, tbl);

#ifdef MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;
}
