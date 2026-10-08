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
  Convert Julian seconds to date.
*/

#include "jurassic.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Print command-line help. */
static void usage(
  void);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  double remain;

  int day, hour, min, mon, sec, year;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 2)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: jsec2time <jsec>\n\n" "Use -h for full help.");

  /* Read arguments... */
  const double jsec = atof(argv[1]);

  /* Convert time... */
  jsec2time(jsec, &year, &mon, &day, &hour, &min, &sec, &remain);
  printf("%d %d %d %d %d %d %g\n", year, mon, day, hour, min, sec, remain);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC time converter.\n\n");
  printf
    ("Convert Julian seconds since 2000-01-01T00:00Z to calendar time.\n\n");
  printf("Usage:\n");
  printf("  jsec2time <jsec>\n\n");
  printf("Arguments:\n");
  printf("  <jsec>  Seconds since 2000-01-01T00:00Z.\n\n");
  printf("Output:\n");
  printf("  Writes results to standard output as:\n");
  printf("  year month day hour minute second remainder\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
