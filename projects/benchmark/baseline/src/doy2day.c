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
  Convert day of year to date.
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

  int day, mon;

  /* Print usage information... */
  USAGE;

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Missing or invalid command-line arguments.\n\n"
	   "Usage: doy2day <year> <doy>\n\n" "Use -h for full help.");

  /* Read arguments... */
  const int year = atoi(argv[1]);
  const int doy = atoi(argv[2]);

  /* Convert... */
  doy2day(year, doy, &mon, &day);
  printf("%d %d %d\n", year, mon, day);

  return EXIT_SUCCESS;
}

/*****************************************************************************/

static void usage(
  void) {
  printf("\nJURASSIC calendar converter.\n\n");
  printf("Convert year and day-of-year to a calendar date.\n\n");
  printf("Usage:\n");
  printf("  doy2day <year> <doy>\n\n");
  printf("Arguments:\n");
  printf("  <year>  Calendar year.\n");
  printf("  <doy>   Day of year.\n\n");
  printf("Output:\n");
  printf("  Writes results to standard output.\n\n");
  printf("Further information:\n");
  printf("  Manual: https://slcs-jsc.github.io/jurassic/\n");
}
