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
  
  Copyright (C) 2003-2025 Forschungszentrum Juelich GmbH
*/

/*! 
  \file
  JURASSIC library declarations.
*/

/*! 
  \mainpage
  
  The Juelich Rapid Spectral Simulation Code (JURASSIC) is a fast
  infrared radiative transfer model for the analysis of atmospheric
  remote sensing measurements.

  \section Introduction

  The source code of JURASSIC is available from the
  [git repository](https://github.com/slcs-jsc/jurassic). Please see the
  [README.md](https://github.com/slcs-jsc/jurassic/blob/master/README.md)
  in the git repository for introductory information. More information
  can be found in the [user manual](https://slcs-jsc.github.io/jurassic).
  
  This doxygen manual contains information about the algorithms and
  data structures used in the code. Please refer to the `jurassic.h'
  documentation for a first overview.
  
  \section References
  
  For citing the model in scientific publications, please see
  [CITATION.cff](https://github.com/slcs-jsc/jurassic/blob/master/CITATION.cff)
  and refer to the following papers:
  
  _Baumeister, P. F. and Hoffmann, L.: Fast infrared radiative
  transfer calculations using graphics processing units: JURASSIC-GPU
  v2.0, Geosci. Model Dev., 15, 1855–1874,
  https://doi.org/10.5194/gmd-15-1855-2022, 2022._
  
  _Hoffmann, L., and M. J. Alexander, Retrieval of stratospheric
  temperatures from Atmospheric Infrared Sounder radiance measurements
  for gravity wave studies, J. Geophys. Res., 114, D07105,
  https://doi.org/10.1029/2008JD011241, 2009._
  
  _Hoffmann, L., Kaufmann, M., Spang, R., Müller, R., Remedios, J. J.,
  Moore, D. P., Volk, C. M., von Clarmann, T., and Riese, M.: Envisat
  MIPAS measurements of CFC-11: retrieval, validation, and
  climatology, Atmos. Chem. Phys., 8, 3671-3688,
  https://doi.org/10.5194/acp-8-3671-2008, 2008._
  
  Additional references are collected here:
  https://slcs-jsc.github.io/jurassic/references
  
  \section License
  
  JURASSIC is being develop at the Jülich Supercomputing Centre,
  Forschungszentrum Jülich, Germany.
  
  JURASSIC is distributed under the terms of the
  [GNU General Public License v3.0](https://github.com/slcs-jsc/jurassic/blob/master/COPYING).
  
  \section Contributing
  
  We are interested in supporting operational and research
  applications with JURASSIC.
  
  You can submit bug reports or feature requests on the
  [issue tracker](https://github.com/slcs-jsc/jurassic/issues).
  
  Proposed code changes and fixes can be submitted as
  [pull requests](https://github.com/slcs-jsc/jurassic/pulls).
  
  Please do not hesitate to contact us if you have any questions or
  need assistance.
  
  \section Contact
  
  Dr. Lars Hoffmann
  
  Jülich Supercomputing Centre, Forschungszentrum Jülich
  
  e-mail: <l.hoffmann@fz-juelich.de>
*/

#ifndef JURASSIC_H
#define JURASSIC_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
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

/*! Compute brightness temperature. */
#define BRIGHT(rad, nu)					\
  (C2 * (nu) / gsl_log1p(C1 * POW3(nu) / (rad)))

/*! Convert degrees to radians. */
#define DEG2RAD(deg)				\
  ((deg) * (M_PI / 180.0))

/*! Compute Cartesian distance between two vectors. */
#define DIST(a, b) sqrt(DIST2(a, b))

/*! Compute squared distance between two vectors. */
#define DIST2(a, b)							\
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/*! Compute dot product of two vectors. */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

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

/*! @brief Macro to determine the maximum of two values. */
#define MAX(a,b)				\
  (((a)>(b))?(a):(b))

/*! Macro to determine the minimum of two values. */
#define MIN(a,b)				\
  (((a)<(b))?(a):(b))

/*! Compute linear interpolation. */
#define LIN(x0, y0, x1, y1, x)			\
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/*! Compute logarithmic interpolation in x. */
#define LOGX(x0, y0, x1, y1, x)				\
  (((x)/(x0)>0 && (x1)/(x0)>0)				\
   ? ((y0)+((y1)-(y0))*log((x)/(x0))/log((x1)/(x0)))	\
   : LIN(x0, y0, x1, y1, x))

/*! Compute logarithmic interpolation in y. */
#define LOGY(x0, y0, x1, y1, x)					\
  (((y1)/(y0)>0)						\
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0))))          \
   : LIN(x0, y0, x1, y1, x))

/*! Compute norm of a vector. */
#define NORM(a) sqrt(DOTP(a, a))

/*! Compute Planck function. */
#define PLANCK(T, nu)				\
  (C1 * POW3(nu) / gsl_expm1(C2 * (nu) / (T)))

/*! Compute x^2. */
#define POW2(x) ((x)*(x))

/*! Compute x^3. */
#define POW3(x) ((x)*(x)*(x))

/*! Convert radians to degrees. */
#define RAD2DEG(rad)				\
  ((rad) * (180.0 / M_PI))

/*! Compute refractivity (return value is n - 1). */
#define REFRAC(p, T)				\
  (7.753e-05 * (p) / (T))

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

/*! First spectroscopic constant (c_1 = 2 h c^2) [W/(m^2 sr cm^-4)]. */
#ifndef C1
#define C1 1.19104259e-8
#endif

/*! Second spectroscopic constant (c_2 = h c / k) [K/cm^-1]. */
#ifndef C2
#define C2 1.43877506
#endif

/*! Minimum emissivity. */
#ifndef EPSMIN
#define EPSMIN 0
#endif

/*! Maximum emissivity. */
#ifndef EPSMAX
#define EPSMAX 1
#endif

/*! Standard gravity [m/s^2]. */
#ifndef G0
#define G0 9.80665
#endif

/*! Standard scale height [km]. */
#ifndef H0
#define H0 7.0
#endif

/*! Boltzmann constant [kg m^2/(K s^2)]. */
#ifndef KB
#define KB 1.3806504e-23
#endif

/*! Mass of Earth [kg]. */
#ifndef ME
#define ME 5.976e24
#endif

/*! Avogadro's number. */
#ifndef NA
#define NA 6.02214199e23
#endif

/*! Mean radius of Earth [km]. */
#ifndef RE
#define RE 6367.421
#endif

/*! Ideal gas constant [J/(mol K)]. */
#ifndef RI
#define RI 8.3144598
#endif

/*! Standard pressure [hPa]. */
#ifndef P0
#define P0 1013.25
#endif

/*! Standard temperature [K]. */
#ifndef T0
#define T0 273.15
#endif

/*! Minimum temperature for source function [K]. */
#ifndef TMIN
#define TMIN 100.
#endif

/*! Maximum temperature for source function [K]. */
#ifndef TMAX
#define TMAX 400.
#endif

/*! Effective temperature of the sun [K]. */
#ifndef TSUN
#define TSUN 5780.
#endif

/*! Minimum column density [molecules/cm^2]. */
#ifndef UMIN
#define UMIN 0
#endif

/*! Maximum column density [molecules/cm^2]. */
#ifndef UMAX
#define UMAX 1e30
#endif

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of cloud layer spectral grid points. */
#ifndef NCL
#define NCL 8
#endif

/*! Maximum number of radiance channels. */
#ifndef ND
#define ND 128
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
#define NLOS 4096
#endif

/*! Maximum number of shape function grid points. */
#ifndef NSHAPE
#define NSHAPE 20000
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

/*! Maximum number of RFM spectral grid points. */
#ifndef RFMNPTS
#define RFMNPTS 10000000
#endif

/*! Maximum length of RFM data lines. */
#ifndef RFMLINE
#define RFMLINE 100000
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

  /*! Extinction [km^-1]. */
  double k[NW][NP];

  /*! Cloud layer height [km]. */
  double clz;

  /*! Cloud layer depth [km]. */
  double cldz;

  /*! Cloud layer extinction [km^-1]. */
  double clk[NCL];

  /*! Surface height [km]. */
  double sfz;

  /*! Surface pressure [hPa]. */
  double sfp;

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[NSF];

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

  /*! Surface treatment (0=none, 1=emissions, 2=downward, 3=solar). */
  int sftype;

  /*! Solar zenith angle at the surface [deg] (-999=auto). */
  double sfsza;

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

  /*! Forward model (0=CGA, 1=EGA, 2=RFM). */
  int formod;

  /*! Path to RFM binary. */
  char rfmbin[LEN];

  /*! HITRAN file for RFM. */
  char rfmhit[LEN];

  /*! Emitter cross-section files for RFM. */
  char rfmxsc[NG][LEN];

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
  double q[NLOS][NG];

  /*! Extinction [km^-1]. */
  double k[NLOS][ND];

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[ND];

  /*! Segment length [km]. */
  double ds[NLOS];

  /*! Column density [molecules/cm^2]. */
  double u[NLOS][NG];

  /*! Curtis-Godson pressure [hPa]. */
  double cgp[NLOS][NG];

  /*! Curtis-Godson temperature [K]. */
  double cgt[NLOS][NG];

  /*! Curtis-Godson column density [molecules/cm^2]. */
  double cgu[NLOS][NG];

  /*! Segment emissivity. */
  double eps[NLOS][ND];

  /*! Segment source function [W/(m^2 sr cm^-1)]. */
  double src[NLOS][ND];

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
  int np[ND][NG];

  /*! Number of temperatures. */
  int nt[ND][NG][TBLNP];

  /*! Number of column densities. */
  int nu[ND][NG][TBLNP][TBLNT];

  /*! Pressure [hPa]. */
  double p[ND][NG][TBLNP];

  /*! Temperature [K]. */
  double t[ND][NG][TBLNP][TBLNT];

  /*! Column density [molecules/cm^2]. */
  float u[ND][NG][TBLNP][TBLNT][TBLNU];

  /*! Emissivity. */
  float eps[ND][NG][TBLNP][TBLNT][TBLNU];

  /*! Source function temperature [K]. */
  double st[TBLNS];

  /*! Source function radiance [W/(m^2 sr cm^-1)]. */
  double sr[TBLNS][ND];

} tbl_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Compose state vector or parameter vector. */
size_t atm2x(
  const ctl_t * ctl,
  const atm_t * atm,
  gsl_vector * x,
  int *iqa,
  int *ipa);

/*! Add element to state vector. */
void atm2x_help(
  const double value,
  const int value_iqa,
  const int value_ip,
  gsl_vector * x,
  int *iqa,
  int *ipa,
  size_t *n);

/*! Convert Cartesian coordinates to geolocation. */
void cart2geo(
  const double *x,
  double *z,
  double *lon,
  double *lat);

/*! Interpolate climatological data. */
void climatology(
  const ctl_t * ctl,
  atm_t * atm_mean);

/*! Compute carbon dioxide continuum (optical depth). */
double ctmco2(
  const double nu,
  const double p,
  const double t,
  const double u);

/*! Compute water vapor continuum (optical depth). */
double ctmh2o(
  const double nu,
  const double p,
  const double t,
  const double q,
  const double u);

/*! Compute nitrogen continuum (absorption coefficient). */
double ctmn2(
  const double nu,
  const double p,
  const double t);

/*! Compute oxygen continuum (absorption coefficient). */
double ctmo2(
  const double nu,
  const double p,
  const double t);

/*! Copy and initialize atmospheric data. */
void copy_atm(
  const ctl_t * ctl,
  atm_t * atm_dest,
  const atm_t * atm_src,
  const int init);

/*! Copy and initialize observation data. */
void copy_obs(
  const ctl_t * ctl,
  obs_t * obs_dest,
  const obs_t * obs_src,
  const int init);

/*! Find index of an emitter. */
int find_emitter(
  const ctl_t * ctl,
  const char *emitter);

/*! Determine ray paths and compute radiative transfer. */
void formod(
  const ctl_t * ctl,
  atm_t * atm,
  obs_t * obs);

/*! Compute absorption coefficient of continua. */
void formod_continua(
  const ctl_t * ctl,
  const los_t * los,
  const int ip,
  double *beta);

/*! Apply field of view convolution. */
void formod_fov(
  const ctl_t * ctl,
  obs_t * obs);

/*! Compute radiative transfer for a pencil beam. */
void formod_pencil(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs,
  const int ir);

/*! Apply RFM for radiative transfer calculations. */
void formod_rfm(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs);

/*! Compute Planck source function. */
void formod_srcfunc(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const double t,
  double *src);

/*! Convert geolocation to Cartesian coordinates. */
void geo2cart(
  const double z,
  const double lon,
  const double lat,
  double *x);

/*! Set hydrostatic equilibrium. */
void hydrostatic(
  const ctl_t * ctl,
  atm_t * atm);

/*! Determine name of state vector quantity for given index. */
void idx2name(
  const ctl_t * ctl,
  const int idx,
  char *quantity);

/*! Initialize source function table. */
void init_srcfunc(
  const ctl_t * ctl,
  tbl_t * tbl);

/*! Interpolate atmospheric data. */
void intpol_atm(
  const ctl_t * ctl,
  const atm_t * atm,
  const double z,
  double *p,
  double *t,
  double *q,
  double *k);

/*! Get transmittance from look-up tables (CGA method). */
void intpol_tbl_cga(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/*! Get transmittance from look-up tables (EGA method). */
void intpol_tbl_ega(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/*! Interpolate emissivity from look-up tables. */
double intpol_tbl_eps(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double u);

/*! Interpolate column density from look-up tables. */
double intpol_tbl_u(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double eps);

/*! Convert seconds to date. */
void jsec2time(
  const double jsec,
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
  const double *xx,
  const int n,
  const double x);

/*! Find array index for regular grid. */
int locate_reg(
  const double *xx,
  const int n,
  const double x);

/*! Find array index in float array. */
int locate_tbl(
  const float *xx,
  const int n,
  const double x);

/*! Compose measurement vector. */
size_t obs2y(
  const ctl_t * ctl,
  const obs_t * obs,
  gsl_vector * y,
  int *ida,
  int *ira);

/*! Do ray-tracing to determine LOS. */
void raytrace(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs,
  los_t * los,
  const int ir);

/*! Read atmospheric data. */
void read_atm(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
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
  const ctl_t * ctl,
  obs_t * obs);

/*! Read observation data in RFM format. */
double read_obs_rfm(
  const char *basename,
  const double z,
  double *nu,
  double *f,
  int n);

/*! Read RFM spectrum. */
void read_rfm_spec(
  const char *filename,
  double *nu,
  double *rad,
  int *npts);

/*! Read shape function. */
void read_shape(
  const char *filename,
  double *x,
  double *y,
  int *n);

/*! Read look-up table data. */
void read_tbl(
  const ctl_t * ctl,
  tbl_t * tbl);

/*! Search control parameter file for variable entry. */
double scan_ctl(
  int argc,
  char *argv[],
  const char *varname,
  int arridx,
  const char *defvalue,
  char *value);

/*! Calculate solar zenith angle. */
double sza(
  double sec,
  double lon,
  double lat);

/*! Find tangent point of a given LOS. */
void tangent_point(
  const los_t * los,
  double *tpz,
  double *tplon,
  double *tplat);

/*! Convert date to seconds. */
void time2jsec(
  const int year,
  const int mon,
  const int day,
  const int hour,
  const int min,
  const int sec,
  const double remain,
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
  const ctl_t * ctl,
  const atm_t * atm);

/*! Write atmospheric data in RFM format. */
void write_atm_rfm(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/*! Write matrix. */
void write_matrix(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const gsl_matrix * matrix,
  const atm_t * atm,
  const obs_t * obs,
  const char *rowspace,
  const char *colspace,
  const char *sort);

/*! Write observation data. */
void write_obs(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs);

/*! Write shape function. */
void write_shape(
  const char *filename,
  const double *x,
  const double *y,
  const int n);

/*! Write look-up table data. */
void write_tbl(
  const ctl_t * ctl,
  const tbl_t * tbl);

/*! Decompose parameter vector or state vector. */
void x2atm(
  const ctl_t * ctl,
  const gsl_vector * x,
  atm_t * atm);

/*! Get element from state vector. */
void x2atm_help(
  double *value,
  const gsl_vector * x,
  size_t *n);

/*! Decompose measurement vector. */
void y2obs(
  const ctl_t * ctl,
  const gsl_vector * y,
  obs_t * obs);

#endif
