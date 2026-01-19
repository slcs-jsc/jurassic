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
  
  Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Reduce look-up table size using recursive, error-controlled
  interpolation (Ramer-Douglas-Peucker-style).
*/


#include "jurassic.h"

#define ABS_ERROR_RATE (1e-2)
#define REL_ERROR_RATE (1e-2)

#define REL_ERROR_VIOLATION(diff, val) \
    ((diff) > (REL_ERROR_RATE * fabs(val)))
#define ABS_ERROR_VIOLATION(diff, val) \
    ((diff) > ABS_ERROR_RATE)

/* Copy density entry for given indices... */
static void copy_density(
  tbl_t *tbl,
  int id,
  int ig,
  int ip,
  int it,
  int dst,
  int src) {
  tbl->logu[id][ig][ip][it][dst] = tbl->logu[id][ig][ip][it][src];
  tbl->logeps[id][ig][ip][it][dst] = tbl->logeps[id][ig][ip][it][src];
}

/**
 * Adaptive reduction of a (u, eps) table segment.
 *
 * Endpoints are preserved; interior points are retained only if their
 * linear interpolation error exceeds the specified tolerance.
 * Reduction is performed in-place and preserves ordering in u.
 */
static void reduction(
  tbl_t *tbl,
  int id,
  int ig,
  int ip,
  int it,
  int left,
  int right,
  int *write_idx) {

  int worstindex = -1;
  double max_err = 0.0;

  /* Scan interior points and find the largest interpolation error... */
  for (int iu = left + 1; iu < right; iu++) {

    /* Linear interpolation in log-log space: interpolate log(eps) linearly over log(u)... */
    const double interp_logeps =
      LIN(tbl->logu[id][ig][ip][it][left], tbl->logeps[id][ig][ip][it][left],
	  tbl->logu[id][ig][ip][it][right],
	  tbl->logeps[id][ig][ip][it][right],
	  tbl->logu[id][ig][ip][it][iu]);

    /* Compute relative error in linear eps (eps = exp(logeps)) */
    const double eps_i = exp((double) tbl->logeps[id][ig][ip][it][iu]);
    const double interp_eps = exp(interp_logeps);
    const double diff = fabs(eps_i - interp_eps);

    /* Check interpolation errors and find maximum error... */
    if (REL_ERROR_VIOLATION(diff, eps_i) && diff > max_err) {
      max_err = diff;
      worstindex = iu;
    }
  }

  /* Check whether any point violated the error threshold... */
  if (worstindex != -1) {

    /* Recursion on left segment... */
    if (worstindex - left > 1)
      reduction(tbl, id, ig, ip, it, left, worstindex, write_idx);

    /* Keep worst entry... */
    copy_density(tbl, id, ig, ip, it, *write_idx, worstindex);
    (*write_idx)++;

    /* Recursion on right segment... */
    if (right - worstindex > 1)
      reduction(tbl, id, ig, ip, it, worstindex, right, write_idx);
  }
}


int main(
  int argc,
  char *argv[]) {
  ctl_t ctl;

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("tblrdc: Give parameters: <ctl> <tblbase_in> <tblfmt_in>"
	   " <tblbase_out> <tblfmt_out>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);

  /* Read tables... */
  sprintf(ctl.tblbase, "%s", argv[2]);
  ctl.tblfmt = atoi(argv[3]);
  tbl_t *tbl = read_tbl(&ctl);

  /* Loop over trace gases and channels... */
  for (int id = 0; id < ctl.nd; id++)
    for (int ig = 0; ig < ctl.ng; ig++) {

      /* Init counters... */
      int total_input = 0;
      int total_output = 0;

      /* Loop over (p, T) combinations... */
      for (int ip = 0; ip < tbl->np[id][ig]; ip++) {
	for (int it = 0; it < tbl->nt[id][ig][ip]; it++) {

	  const int num_points = tbl->nu[id][ig][ip][it];

	  /* First entry is implicitly kept (in-place reduction)... */
	  int write_idx = 1;

	  /* Reduce Middle-Entries if more than 2 entries... */
	  if (num_points > 2)
	    reduction(tbl, id, ig, ip, it, 0, num_points - 1, &write_idx);

	  /* Keep last entry if more than 1 entry... */
	  if (num_points > 1) {
	    copy_density(tbl, id, ig, ip, it, write_idx, num_points - 1);
	    write_idx++;
	  }

	  /* Update number of kept points... */
	  printf("Block of %d Entries has been reduced to: %d \n",
		 num_points, write_idx);
	  tbl->nu[id][ig][ip][it] = write_idx;
	  total_input += num_points;
	  total_output += write_idx;
	}
      }

      /* Write info... */
      if (total_input != 0) {
	printf
	  ("\n In total the table has been reduced from %d entries down to %d entries, this represents a %.1f%% reduction.\n\n",
	   total_input, total_output,
	   100.0 * (total_input - total_output) / total_input);
      }
    }

  /* Write tables... */
  sprintf(ctl.tblbase, "%s", argv[4]);
  ctl.tblfmt = atoi(argv[5]);
  write_tbl(&ctl, tbl);

  /* Free... */
  free(tbl);

  return EXIT_SUCCESS;
}
