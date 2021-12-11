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
  Create radiometric filter functions.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Compute apodized instrument line shape. */
double ails(
  int apo,
  double opl,
  double dnu);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static ctl_t ctl;

  static double ff[NSHAPE], fnu[NSHAPE];

  double fsum = 0.0;

  int fn = 0;

  /* Write info... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <filter>");

  /* Read control parameters... */
  read_ctl(argc, argv, &ctl);
  int type = (int) scan_ctl(argc, argv, "FILTER_TYPE", -1, "1", NULL);
  double opd = scan_ctl(argc, argv, "FILTER_OPD", -1, "10.0", NULL);
  double fwhm = scan_ctl(argc, argv, "FILTER_FWHM", -1, "1.0", NULL);
  double center = scan_ctl(argc, argv, "FILTER_CENTER", -1, "1000.0", NULL);
  double width = scan_ctl(argc, argv, "FILTER_WIDTH", -1, "2.1", NULL);
  double samp = scan_ctl(argc, argv, "FILTER_SAMP", -1, "0.0005", NULL);

  /* Compute filter function... */
  for (double nu = center - 0.5 * width;
       nu <= center + 0.5 * width; nu += samp) {

    /* Set frequency... */
    fnu[fn] = nu;

    /* Boxcar... */
    if (type == 0)
      ff[fn] = (fabs(nu - center) <= 0.5 * fwhm ? 1.0 : 0.0);

    /* Triangle... */
    else if (type == 1) {
      ff[fn] = 1.0 - fabs(nu - center) / fwhm;
      ff[fn] = GSL_MAX(ff[fn], 0.0);
    }

    /* Gaussian... */
    else if (type == 2) {
      double sigma = fwhm / 2.355;
      ff[fn] = exp(-0.5 * POW2((nu - center) / sigma));
    }

    /* Sinc function... */
    else if (type == 3)
      ff[fn] = ails(0, opd, nu - center);

    /* Norton-Beer strong apodization... */
    else if (type == 4)
      ff[fn] = ails(1, opd, nu - center);

    /* Error message... */
    else
      ERRMSG("Filter function type unknown!");

    /* Count spectral grid points... */
    if ((++fn) >= NSHAPE)
      ERRMSG("Too many filter function data points!");
  }

  /* Normalize filter function... */
  for (int i = 0; i < fn; i++)
    fsum += ff[i];
  for (int i = 0; i < fn; i++)
    ff[i] /= (fsum * samp);

  /* Write to file... */
  write_shape(argv[2], fnu, ff, fn);

  return (EXIT_SUCCESS);
}

/*****************************************************************************/

double ails(
  int apo,
  double opl,
  double dnu) {

  double a, a2, a4, a6, a8, cosa, q0, q2, q4, sinca;

  /* Auxiliary quantities... */
  a = 2 * M_PI * dnu * opl;
  a2 = a * a;
  a4 = a2 * a2;
  a6 = a4 * a2;
  a8 = a4 * a4;

  /* Sinc function... */
  if (apo == 0) {
    if (fabs(a) < 0.7)
      return 1 - a2 / 6 + a4 / 120 - a6 / 5040 + a8 / 362880;
    else
      return sin(a) / a;
  }

  /* Norton-Beer strong apodization... */
  else if (apo == 1) {
    if (fabs(a) < 0.7) {
      q0 = 1 - a2 / 6 + a4 / 120 - a6 / 5040 + a8 / 362880;
      q2 = 1 - a2 / 14 + a4 / 504 - a6 / 33264 + a8 / 3459456;
      q4 = 1 - a2 / 22 + a4 / 1144 - a6 / 102960 + a8 / 14002560;
    } else {
      sinca = sin(a) / a;
      cosa = cos(a);
      q0 = sinca;
      q2 = -15 / a2 * ((1 - 3 / a2) * sinca + (3 / a2) * cosa);
      q4 =
	945 / a4 * ((1 - 45 / a2 + 105 / a4) * sinca +
		    5 / a2 * (2 - 21 / a2) * cosa);
    }
    return 0.045335 * q0 + 0.554883 * q2 * 8. / 15. +
      0.399782 * q4 * 384. / 945.;
  }

  /* Error message.... */
  else
    ERRMSG("Unknown apodization!");
}
