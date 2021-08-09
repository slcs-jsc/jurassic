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
  JURASSIC library declarations.
*/

/*!
  \mainpage
  
  The JUelich RApid Spectral SImulation Code (JURASSIC) is a fast radiative
  transfer model for the mid-infrared spectral region.

  This reference manual provides information on the algorithms
  and data structures used in the code.

  Further information can be found at: https://github.com/slcs-jsc/jurassic
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Allocate memory. */
#define ALLOC(ptr, type, n)				\
  if((ptr=malloc((size_t)(n)*sizeof(type)))==NULL)	\
    ERRMSG("Out of memory!");

/*! Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)							\
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/*! Compute exponential interpolation. */
#define EXP(x0, y0, x1, y1, x)					\
  (((y0)>0 && (y1)>0)						\
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))          \
   : LIN(x0, y0, x1, y1, x))

/*! Read binary data. */
#define FREAD(ptr, type, size, out) {					\
    if(fread(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while reading!");					\
  }

/*! Write binary data. */
#define FWRITE(ptr, type, size, out) {					\
    if(fwrite(ptr, sizeof(type), size, out)!=size)			\
      ERRMSG("Error while writing!");					\
  }

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/*! Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/*! Compute x^2. */
#define POW2(x) ((x)*(x))

/*! Compute x^3. */
#define POW3(x) ((x)*(x)*(x))

/*! Start or stop a timer. */
#define TIMER(name, mode)				\
  {timer(name, __FILE__, __func__, __LINE__, mode);}

/*! Read string tokens. */
#define TOK(line, tok, format, var) {			\
    if(((tok)=strtok((line), " \t"))) {			\
      if(sscanf(tok, format, &(var))!=1) continue;	\
    } else ERRMSG("Error while reading!");		\
  }

/* ------------------------------------------------------------
   Log messages...
   ------------------------------------------------------------ */

/*! Level of log messages (0=none, 1=basic, 2=detailed, 3=debug). */
#ifndef LOGLEV
#define LOGLEV 2
#endif

/*! Print log message. */
#define LOG(level, ...) {						\
    if(level >= 2)							\
      printf("  ");							\
    if(level <= LOGLEV) {						\
      printf(__VA_ARGS__);						\
      printf("\n");							\
    }									\
  }

/*! Print warning message. */
#define WARN(...) {							\
    printf("\nWarning (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
  }

/*! Print error message and quit program. */
#define ERRMSG(...) {							\
    printf("\nError (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
    exit(EXIT_FAILURE);							\
  }

/*! Print macro for debugging. */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Minimum temperature for source function [K]. */
#define TMIN 100.

/*! Maximum temperature for source function [K]. */
#define TMAX 400.

/*! First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#define C1 1.19104259e-8

/*! Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#define C2 1.43877506

/*! Standard gravity [m/s^2]. */
#define G0 9.80665

/*! Boltzmann constant [kg m^2/(K s^2)]. */
#define KB 1.3806504e-23

/*! Avogadro's number. */
#define NA 6.02214199e23

/*! Standard scale height [km]. */
#define H0 7.0

/*! Standard pressure [hPa]. */
#define P0 1013.25

/*! Standard temperature [K]. */
#define T0 273.15

/*! Mean radius of Earth [km]. */
#define RE 6367.421

/*! Ideal gas constant [J/(mol K)]. */
#define RI 8.3144598

/*! Mass of Earth [kg]. */
#define ME 5.976e24

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of cloud layer spectral grid points. */
#ifndef NCL
#define NCL 8
#endif

/*! Maximum number of radiance channels. */
#ifndef ND
#define ND 64
#endif

/*! Maximum number of emitters. */
#ifndef NG
#define NG 8
#endif

/*! Maximum number of atmospheric data points. */
#ifndef NP
#define NP 256
#endif

/*! Maximum number of ray paths. */
#ifndef NR
#define NR 256
#endif

/*! Maximum number of surface layer spectral grid points. */
#ifndef NSF
#define NSF 8
#endif

/*! Maximum number of spectral windows. */
#ifndef NW
#define NW 4
#endif

/*! Maximum length of ASCII data lines. */
#ifndef LEN
#define LEN 10000
#endif

/*! Maximum size of measurement vector. */
#ifndef M
#define M (NR*ND)
#endif

/*! Maximum size of state vector. */
#ifndef N
#define N ((2+NG+NW)*NP+NCL+NSF+5)
#endif

/*! Maximum number of quantities. */
#ifndef NQ
#define NQ (7+NG+NW+NCL+NSF)
#endif

/*! Maximum number of LOS points. */
#ifndef NLOS
#define NLOS 2048
#endif

/*! Maximum number of shape function grid points. */
#ifndef NSHAPE
#define NSHAPE 10000
#endif

/*! Number of ray paths used for FOV calculations. */
#ifndef NFOV
#define NFOV 5
#endif

/*! Maximum number of pressure levels in emissivity tables. */
#ifndef TBLNP
#define TBLNP 41
#endif

/*! Maximum number of temperatures in emissivity tables. */
#ifndef TBLNT
#define TBLNT 30
#endif

/*! Maximum number of column densities in emissivity tables. */
#ifndef TBLNU
#define TBLNU 320
#endif

/*! Maximum number of source function temperature levels. */
#ifndef TBLNS
#define TBLNS 1200
#endif

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/*! Index for pressure. */
#define IDXP 0

/*! Index for temperature. */
#define IDXT 1

/*! Indices for volume mixing ratios. */
#define IDXQ(ig) (2+ig)

/*! Indices for extinction. */
#define IDXK(iw) (2+ctl->ng+iw)

/*! Index for cloud layer height. */
#define IDXCLZ (2+ctl->ng+ctl->nw)

/*! Index for cloud layer depth. */
#define IDXCLDZ (3+ctl->ng+ctl->nw)

/*! Indices for cloud layer extinction. */
#define IDXCLK(icl) (4+ctl->ng+ctl->nw+icl)

/*! Index for surface layer height. */
#define IDXSFZ (4+ctl->ng+ctl->nw+ctl->ncl)

/*! Index for surface layer pressure. */
#define IDXSFP (5+ctl->ng+ctl->nw+ctl->ncl)

/*! Index for surface layer temperature. */
#define IDXSFT (6+ctl->ng+ctl->nw+ctl->ncl)

/*! Indices for surface layer emissivity. */
#define IDXSFEPS(isf) (7+ctl->ng+ctl->nw+ctl->ncl+isf)

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! Atmospheric data. */
typedef struct {

  /*! Number of data points. */
  int np;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NP];

  /*! Altitude [km]. */
  double z[NP];

  /*! Longitude [deg]. */
  double lon[NP];

  /*! Latitude [deg]. */
  double lat[NP];

  /*! Pressure [hPa]. */
  double p[NP];

  /*! Temperature [K]. */
  double t[NP];

  /*! Volume mixing ratio [ppv]. */
  double q[NG][NP];

  /*! Extinction [1/km]. */
  double k[NW][NP];

  /*! Cloud layer height [km]. */
  double clz;

  /*! Cloud layer depth [km]. */
  double cldz;

  /*! Cloud layer extinction [1/km]. */
  double clk[NCL];

  /*! Surface height [km]. */
  double sfz;

  /*! Surface pressure [hPa]. */
  double sfp;

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[NCL];

} atm_t;

/*! Forward model control parameters. */
typedef struct {

  /*! Number of emitters. */
  int ng;

  /*! Name of each emitter. */
  char emitter[NG][LEN];

  /*! Number of radiance channels. */
  int nd;

  /*! Centroid wavenumber of each channel [cm^-1]. */
  double nu[ND];

  /*! Number of spectral windows. */
  int nw;

  /*! Window index of each channel. */
  int window[ND];

  /*! Number of cloud layer spectral grid points. */
  int ncl;

  /*! Cloud layer wavenumber [cm^-1]. */
  double clnu[NCL];

  /*! Number of surface layer spectral grid points. */
  int nsf;

  /*! Surface layer wavenumber [cm^-1]. */
  double sfnu[NSF];

  /*! Basename for table files and filter function files. */
  char tblbase[LEN];

  /*! Look-up table file format (1=ASCII, 2=binary). */
  int tblfmt;

  /*! Reference height for hydrostatic pressure profile (-999 to skip) [km]. */
  double hydz;

  /*! Compute CO2 continuum (0=no, 1=yes). */
  int ctm_co2;

  /*! Compute H2O continuum (0=no, 1=yes). */
  int ctm_h2o;

  /*! Compute N2 continuum (0=no, 1=yes). */
  int ctm_n2;

  /*! Compute O2 continuum (0=no, 1=yes). */
  int ctm_o2;

  /*! Take into account refractivity (0=no, 1=yes). */
  int refrac;

  /*! Maximum step length for raytracing [km]. */
  double rayds;

  /*! Vertical step length for raytracing [km]. */
  double raydz;

  /*! Field-of-view data file. */
  char fov[LEN];

  /*! Minimum altitude for pressure retrieval [km]. */
  double retp_zmin;

  /*! Maximum altitude for pressure retrieval [km]. */
  double retp_zmax;

  /*! Minimum altitude for temperature retrieval [km]. */
  double rett_zmin;

  /*! Maximum altitude for temperature retrieval [km]. */
  double rett_zmax;

  /*! Minimum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmin[NG];

  /*! Maximum altitude for volume mixing ratio retrieval [km]. */
  double retq_zmax[NG];

  /*! Minimum altitude for extinction retrieval [km]. */
  double retk_zmin[NW];

  /*! Maximum altitude for extinction retrieval [km]. */
  double retk_zmax[NW];

  /*! Retrieve cloud layer height (0=no, 1=yes). */
  int ret_clz;

  /*! Retrieve cloud layer depth (0=no, 1=yes). */
  int ret_cldz;

  /*! Retrieve cloud layer extinction (0=no, 1=yes). */
  int ret_clk;

  /*! Retrieve surface layer height (0=no, 1=yes). */
  int ret_sfz;

  /*! Retrieve surface layer pressure (0=no, 1=yes). */
  int ret_sfp;

  /*! Retrieve surface layer temperature (0=no, 1=yes). */
  int ret_sft;

  /*! Retrieve surface layer emissivity (0=no, 1=yes). */
  int ret_sfeps;

  /*! Use brightness temperature instead of radiance (0=no, 1=yes). */
  int write_bbt;

  /*! Write matrix file (0=no, 1=yes). */
  int write_matrix;

} ctl_t;

/*! Line-of-sight data. */
typedef struct {

  /*! Number of LOS points. */
  int np;

  /*! Altitude [km]. */
  double z[NLOS];

  /*! Longitude [deg]. */
  double lon[NLOS];

  /*! Latitude [deg]. */
  double lat[NLOS];

  /*! Pressure [hPa]. */
  double p[NLOS];

  /*! Temperature [K]. */
  double t[NLOS];

  /*! Volume mixing ratio [ppv]. */
  double q[NG][NLOS];

  /*! Extinction [1/km]. */
  double k[ND][NLOS];

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[NSF];

  /*! Segment length [km]. */
  double ds[NLOS];

  /*! Column density [molecules/cm^2]. */
  double u[NG][NLOS];

} los_t;

/*! Observation geometry and radiance data. */
typedef struct {

  /*! Number of ray paths. */
  int nr;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NR];

  /*! Observer altitude [km]. */
  double obsz[NR];

  /*! Observer longitude [deg]. */
  double obslon[NR];

  /*! Observer latitude [deg]. */
  double obslat[NR];

  /*! View point altitude [km]. */
  double vpz[NR];

  /*! View point longitude [deg]. */
  double vplon[NR];

  /*! View point latitude [deg]. */
  double vplat[NR];

  /*! Tangent point altitude [km]. */
  double tpz[NR];

  /*! Tangent point longitude [deg]. */
  double tplon[NR];

  /*! Tangent point latitude [deg]. */
  double tplat[NR];

  /*! Transmittance of ray path. */
  double tau[ND][NR];

  /*! Radiance [W/(m^2 sr cm^-1)]. */
  double rad[ND][NR];

} obs_t;

/*! Emissivity look-up tables. */
typedef struct {

  /*! Number of pressure levels. */
  int np[NG][ND];

  /*! Number of temperatures. */
  int nt[NG][ND][TBLNP];

  /*! Number of column densities. */
  int nu[NG][ND][TBLNP][TBLNT];

  /*! Pressure [hPa]. */
  double p[NG][ND][TBLNP];

  /*! Temperature [K]. */
  double t[NG][ND][TBLNP][TBLNT];

  /*! Column density [molecules/cm^2]. */
  float u[NG][ND][TBLNP][TBLNT][TBLNU];

  /*! Emissivity. */
  float eps[NG][ND][TBLNP][TBLNT][TBLNU];

  /*! Source function temperature [K]. */
  double st[TBLNS];

  /*! Source function radiance [W/(m^2 sr cm^-1)]. */
  double sr[ND][TBLNS];

} tbl_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Compose state vector or parameter vector. */
size_t atm2x(
  ctl_t * ctl,
  atm_t * atm,
  gsl_vector * x,
  int *iqa,
  int *ipa);

/*! Add element to state vector. */
void atm2x_help(
  double value,
  int value_iqa,
  int value_ip,
  gsl_vector * x,
  int *iqa,
  int *ipa,
  size_t *n);

/*! Compute brightness temperature. */
double brightness(
  double rad,
  double nu);

/*! Convert Cartesian coordinates to geolocation. */
void cart2geo(
  double *x,
  double *z,
  double *lon,
  double *lat);

/*! Interpolate climatological data. */
void climatology(
  ctl_t * ctl,
  atm_t * atm_mean);

/*! Compute carbon dioxide continuum (optical depth). */
double ctmco2(
  double nu,
  double p,
  double t,
  double u);

/*! Compute water vapor continuum (optical depth). */
double ctmh2o(
  double nu,
  double p,
  double t,
  double q,
  double u);

/*! Compute nitrogen continuum (absorption coefficient). */
double ctmn2(
  double nu,
  double p,
  double t);

/*! Compute oxygen continuum (absorption coefficient). */
double ctmo2(
  double nu,
  double p,
  double t);

/*! Copy and initialize atmospheric data. */
void copy_atm(
  ctl_t * ctl,
  atm_t * atm_dest,
  atm_t * atm_src,
  int init);

/*! Copy and initialize observation data. */
void copy_obs(
  ctl_t * ctl,
  obs_t * obs_dest,
  obs_t * obs_src,
  int init);

/*! Find index of an emitter. */
int find_emitter(
  ctl_t * ctl,
  const char *emitter);

/*! Determine ray paths and compute radiative transfer. */
void formod(
  ctl_t * ctl,
  atm_t * atm,
  obs_t * obs);

/*! Compute absorption coefficient of continua. */
void formod_continua(
  ctl_t * ctl,
  los_t * los,
  int ip,
  double *beta);

/*! Apply field of view convolution. */
void formod_fov(
  ctl_t * ctl,
  obs_t * obs);

/*! Compute radiative transfer for a pencil beam. */
void formod_pencil(
  ctl_t * ctl,
  atm_t * atm,
  obs_t * obs,
  int ir);

/*! Compute Planck source function. */
void formod_srcfunc(
  ctl_t * ctl,
  tbl_t * tbl,
  double t,
  double *src);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  double z,
  double lon,
  double lat,
  double *x);

/*! Set hydrostatic equilibrium. */
void hydrostatic(
  ctl_t * ctl,
  atm_t * atm);

/*! Determine name of state vector quantity for given index. */
void idx2name(
  ctl_t * ctl,
  int idx,
  char *quantity);

/*! Initialize source function table. */
void init_srcfunc(
  ctl_t * ctl,
  tbl_t * tbl);

/*! Interpolate atmospheric data. */
void intpol_atm(
  ctl_t * ctl,
  atm_t * atm,
  double z,
  double *p,
  double *t,
  double *q,
  double *k);

/*! Get transmittance from look-up tables. */
void intpol_tbl(
  ctl_t * ctl,
  tbl_t * tbl,
  los_t * los,
  int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/*! Interpolate emissivity from look-up tables. */
double intpol_tbl_eps(
  tbl_t * tbl,
  int ig,
  int id,
  int ip,
  int it,
  double u);

/*! Interpolate column density from look-up tables. */
double intpol_tbl_u(
  tbl_t * tbl,
  int ig,
  int id,
  int ip,
  int it,
  double eps);

/*! Convert seconds to date. */
void jsec2time(
  double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/*! Compute Jacobians. */
void kernel(
  ctl_t * ctl,
  atm_t * atm,
  obs_t * obs,
  gsl_matrix * k);

/*! Find array index for irregular grid. */
int locate_irr(
  double *xx,
  int n,
  double x);

/*! Find array index for regular grid. */
int locate_reg(
  double *xx,
  int n,
  double x);

/*! Find array index in float array. */
int locate_tbl(
  float *xx,
  int n,
  double x);

/*! Compose measurement vector. */
size_t obs2y(
  ctl_t * ctl,
  obs_t * obs,
  gsl_vector * y,
  int *ida,
  int *ira);

/*! Compute Planck function. */
double planck(
  double t,
  double nu);

/*! Do ray-tracing to determine LOS. */
void raytrace(
  ctl_t * ctl,
  atm_t * atm,
  obs_t * obs,
  los_t * los,
  int ir);

/*! Read atmospheric data. */
void read_atm(
  const char *dirname,
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Read forward model control parameters. */
void read_ctl(
  int argc,
  char *argv[],
  ctl_t * ctl);

/*! Read matrix. */
void read_matrix(
  const char *dirname,
  const char *filename,
  gsl_matrix * matrix);

/*! Read observation data. */
void read_obs(
  const char *dirname,
  const char *filename,
  ctl_t * ctl,
  obs_t * obs);

/*! Read shape function. */
void read_shape(
  const char *filename,
  double *x,
  double *y,
  int *n);

/*! Read look-up table data. */
void read_tbl(
  ctl_t * ctl,
  tbl_t * tbl);

/*! Compute refractivity (return value is n - 1). */
double refractivity(
  double p,
  double t);

/*! Search control parameter file for variable entry. */
double scan_ctl(
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value);

/*! Find tangent point of a given LOS. */
void tangent_point(
  los_t * los,
  double *tpz,
  double *tplon,
  double *tplat);

/*! Convert date to seconds. */
void time2jsec(
  int year,
  int mon,
  int day,
  int hour,
  int min,
  int sec,
  double remain,
  double *jsec);

/*! Measure wall-clock time. */
void timer(
  const char *name,
  const char *file,
  const char *func,
  int line,
  int mode);

/*! Write atmospheric data. */
void write_atm(
  const char *dirname,
  const char *filename,
  ctl_t * ctl,
  atm_t * atm);

/*! Write matrix. */
void write_matrix(
  const char *dirname,
  const char *filename,
  ctl_t * ctl,
  gsl_matrix * matrix,
  atm_t * atm,
  obs_t * obs,
  const char *rowspace,
  const char *colspace,
  const char *sort);

/*! Write observation data. */
void write_obs(
  const char *dirname,
  const char *filename,
  ctl_t * ctl,
  obs_t * obs);

/*! Write shape function. */
void write_shape(
  const char *filename,
  double *x,
  double *y,
  int n);

/*! Write look-up table data. */
void write_tbl(
  ctl_t * ctl,
  tbl_t * tbl);

/*! Decompose parameter vector or state vector. */
void x2atm(
  ctl_t * ctl,
  gsl_vector * x,
  atm_t * atm);

/*! Get element from state vector. */
void x2atm_help(
  double *value,
  gsl_vector * x,
  size_t *n);

/*! Decompose measurement vector. */
void y2obs(
  ctl_t * ctl,
  gsl_vector * y,
  obs_t * obs);
