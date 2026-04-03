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
  
  Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  Reduce look-up table size using recursive, error-controlled interpolation.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Copy density entry for given indices. */
static void copy_density(
  tbl_t * tbl,
  int id,
  int ig,
  int ip,
  int it,
  int dst,
  int src);

/*! Adaptive reduction of a (u, eps) table segment. */
static void reduction(
  tbl_t * tbl,
  int id,
  int ig,
  int ip,
  int it,
  int left,
  int right,
  int *write_idx,
  double atol,
  double rtol);

/*! Print command-line help. */
static void usage(void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {
  ctl_t ctl;

  if (argc == 2
      && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
    usage();
    return EXIT_SUCCESS;
  }

  /* Check arguments... */
  if (argc < 6)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: tblrdc <ctl> <tblbase_in> <tblfmt_in> <tblbase_out> <tblfmt_out> [KEY VALUE ...]\n\n"
	   "Use -h for full help.");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  const double atol = scan_ctl(argc, argv, "ATOL", -1, "1e-9", NULL);
  const double rtol = scan_ctl(argc, argv, "RTOL", -1, "0.01", NULL);

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
	    reduction(tbl, id, ig, ip, it, 0, num_points - 1, &write_idx,
		      atol, rtol);

	  /* Keep last entry if more than 1 entry... */
	  if (num_points > 1) {
	    copy_density(tbl, id, ig, ip, it, write_idx, num_points - 1);
	    write_idx++;
	  }

	  /* Update number of kept points... */
	  tbl->nu[id][ig][ip][it] = write_idx;
	  total_input += num_points;
	  total_output += write_idx;
	}
      }

      /* Write info... */
      LOG(1, "Reduction: %s | %.4f cm^-1 | n= %d | CR= %g", ctl.emitter[ig],
	  ctl.nu[id], total_input,
	  total_output > 0 ? (double) total_input / total_output : NAN);

    }

  /* Write tables... */
  sprintf(ctl.tblbase, "%s", argv[4]);
  ctl.tblfmt = atoi(argv[5]);
  write_tbl(&ctl, tbl);

  /* Free... */
  tbl_free(&ctl, tbl);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(void) {
  printf("\nJURASSIC lookup-table reduction tool.\n\n");
  printf("Reduce look-up table size by adaptive interpolation with\n");
  printf("absolute and relative error thresholds.\n\n");
  printf("Usage:\n");
  printf("  tblrdc <ctl> <tblbase_in> <tblfmt_in> <tblbase_out> <tblfmt_out> [KEY VALUE ...]\n\n");
  printf("Arguments:\n");
  printf("  <ctl>         Control file.\n");
  printf("  <tblbase_in>  Input table base name.\n");
  printf("  <tblfmt_in>   Input table format identifier.\n");
  printf("  <tblbase_out> Output table base name.\n");
  printf("  <tblfmt_out>  Output table format identifier.\n");
  printf("  [KEY VALUE]   Optional control parameters.\n\n");
  printf("Optional control overrides:\n");
  printf("  ATOL <x>      Absolute interpolation-error threshold.\n");
  printf("  RTOL <x>      Relative interpolation-error threshold.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}

/*****************************************************************************/

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

/*****************************************************************************/

static void reduction(
  tbl_t *tbl,
  int id,
  int ig,
  int ip,
  int it,
  int left,
  int right,
  int *write_idx,
  double atol,
  double rtol) {

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
    if (diff > atol + rtol * fabs(eps_i) && diff > max_err) {
      max_err = diff;
      worstindex = iu;
    }
  }

  /* Check whether any point violated the error threshold... */
  if (worstindex != -1) {

    /* Recursion on left segment... */
    if (worstindex - left > 1)
      reduction(tbl, id, ig, ip, it, left, worstindex, write_idx, atol, rtol);

    /* Keep worst entry... */
    copy_density(tbl, id, ig, ip, it, *write_idx, worstindex);
    (*write_idx)++;

    /* Recursion on right segment... */
    if (right - worstindex > 1)
      reduction(tbl, id, ig, ip, it, worstindex, right, write_idx, atol,
		rtol);
  }
}
