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

/* ------------------------------------------------------------
   Includes...
   ------------------------------------------------------------ */

#include <errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include <netcdf.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

/*! Nitrogen concentration. */
#ifndef N2
#define N2 0.78084
#endif

/*! Oxygen concentration. */
#ifndef O2
#define O2 0.20946
#endif

/*! Standard pressure [hPa]. */
#ifndef P0
#define P0 1013.25
#endif

/*! Mean radius of Earth [km]. */
#ifndef RE
#define RE 6367.421
#endif

/*! Ideal gas constant [J/(mol K)]. */
#ifndef RI
#define RI 8.3144598
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
#define N ((2 + NG + NW) * NP + NCL + NSF + 3)
#endif

/*! Maximum number of quantities. */
#ifndef NQ
#define NQ (5 + NG + NW + NCL + NSF)
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

/*! Maximum number of column densities per emissivity curve. */
#ifndef TBLNU
#define TBLNU 512
#endif

/*! Maximum number of source function temperature levels. */
#ifndef TBLNS
#define TBLNS 1200
#endif

/*! Maximum number of RFM spectral grid points. */
#ifndef RFMNPTS
#define RFMNPTS 10000000
#endif

/* ------------------------------------------------------------
   Quantity indices...
   ------------------------------------------------------------ */

/*! Index for pressure. */
#define IDXP 0

/*! Index for temperature. */
#define IDXT 1

/*! Indices for volume mixing ratios. */
#define IDXQ(ig) (2 + (ig))

/*! Indices for extinction. */
#define IDXK(iw) (2 + (ctl->ng) + (iw))

/*! Index for cloud layer height. */
#define IDXCLZ (2 + (ctl->ng) + (ctl->nw))

/*! Index for cloud layer depth. */
#define IDXCLDZ (3 + (ctl->ng) + (ctl->nw))

/*! Indices for cloud layer extinction. */
#define IDXCLK(icl) (4 + (ctl->ng) + (ctl->nw) + (icl))

/*! Index for surface layer temperature. */
#define IDXSFT (4 + (ctl->ng) + (ctl->nw) + (ctl->ncl))

/*! Indices for surface layer emissivity. */
#define IDXSFEPS(isf) (5 + (ctl->ng) + (ctl->nw) + (ctl->ncl) + (isf))

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/**
 * @brief Allocate memory for an array.
 *
 * Allocates a contiguous block of memory for an array of `n` elements of
 * type `type` and assigns the pointer to `ptr`. If allocation fails, the
 * program prints an error message and terminates.
 *
 * @param[out] ptr Pointer to be allocated.
 * @param[in] type Data type of the elements to allocate.
 * @param[in] n Number of elements to allocate.
 *
 * @note Wraps `malloc()` with built-in error handling.
 *
 * @see FREAD, FWRITE
 *
 * @author Lars Hoffmann
 */
#define ALLOC(ptr, type, n)				 \
  if((ptr=calloc((size_t)(n), sizeof(type)))==NULL)      \
    ERRMSG("Out of memory!");

/**
 * @brief Compute brightness temperature from radiance.
 *
 * Computes the equivalent blackbody (brightness) temperature corresponding to
 * a given spectral radiance and wavenumber, using the inverse Planck function.
 * This form assumes the spectroscopic constants `C1` and `C2` are defined for
 * wavenumber units (cm⁻¹).
 *
 * @param[in] rad Spectral radiance [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 * @param[in] nu  Wavenumber [cm⁻¹].
 *
 * @return Brightness temperature [K].
 *
 * @see PLANCK, C1, C2
 *
 * @note Based on Planck’s law in wavenumber form:
 *       \f$ T_b = \frac{c_2 \nu}{\ln\!\left(1 + \frac{c_1 \nu^3}{L_\nu}\right)} \f$
 *       where \f$L_\nu\f$ is the radiance. The ν³ dependence reflects the
 *       definition of `C1` for radiance per unit wavenumber in cm⁻¹.
 *
 * @warning The radiance MUST be in W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹. Using radiance per meter
 *          will produce brightness temperatures incorrect by six orders of
 *          magnitude.
 *
 * @author Lars Hoffmann
 */
#define BRIGHT(rad, nu)					\
  (C2 * (nu) / gsl_log1p(C1 * POW3(nu) / (rad)))

/**
 * @brief Clamp a value to a specified range.
 *
 * Ensures that @p v lies between @p lo and @p hi.
 * If @p v < @p lo, returns @p lo.
 * If @p v > @p hi, returns @p hi.
 * Otherwise, returns @p v unchanged.
 *
 * This macro works with any numeric type (e.g., int, float, double).
 * All arguments are evaluated exactly once — avoid passing expressions
 * with side effects (e.g., ++ operators or function calls).
 *
 * @param v   Input value to clamp.
 * @param lo  Lower bound.
 * @param hi  Upper bound.
 * @return The clamped value between @p lo and @p hi.
 *
 * @author Lars Hoffmann
 */
#define CLAMP(v, lo, hi)				\
  (((v) < (lo)) ? (lo) : (((v) > (hi)) ? (hi) : (v)))

/**
 * @brief Convert degrees to radians.
 *
 * Converts an angle measured in degrees to radians using:
 * \f$ \text{radians} = \text{degrees} \times \frac{\pi}{180} \f$.
 *
 * @param[in] deg Angle in degrees.
 *
 * @return Angle in radians.
 *
 * @see RAD2DEG
 *
 * @author Lars Hoffmann
 */
#define DEG2RAD(deg) ((deg) * (M_PI / 180.0))

/**
 * @brief Compute Cartesian distance between two 3D vectors.
 *
 * Computes the Euclidean distance between two 3D vectors or points.
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Euclidean distance between a and b.
 *
 * @see DIST2
 *
 * @note Equivalent to \f$ \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2 + (z_1-z_2)^2} \f$.
 *
 * @author Lars Hoffmann
 */
#define DIST(a, b) sqrt(DIST2(a, b))

/**
 * @brief Compute squared distance between two 3D vectors.
 *
 * Computes the square of the Euclidean distance between two 3D vectors.
 * Useful when only relative distances are needed (avoids square root).
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Squared distance between a and b.
 *
 * @see DIST
 *
 * @author Lars Hoffmann
 */
#define DIST2(a, b) \
  ((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])+(a[2]-b[2])*(a[2]-b[2]))

/**
 * @brief Compute dot product of two 3D vectors.
 *
 * Computes the scalar (dot) product between two 3-element vectors:
 * \f$ a \cdot b = a_x b_x + a_y b_y + a_z b_z \f$.
 *
 * @param[in] a First vector (array of length 3).
 * @param[in] b Second vector (array of length 3).
 *
 * @return Scalar dot product of a and b.
 *
 * @see NORM
 *
 * @author Lars Hoffmann
 */
#define DOTP(a, b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/**
 * @brief Read binary data from a file.
 *
 * Reads `size` elements of type `type` from the file stream `out` into
 * the memory pointed to by `ptr`. If the number of elements read differs
 * from `size`, an error message is printed and the program terminates.
 *
 * @param[out] ptr Pointer to destination memory buffer.
 * @param[in] type Type of each element to read.
 * @param[in] size Number of elements to read.
 * @param[in,out] out File stream opened for reading.
 *
 * @see FWRITE, ALLOC
 *
 * @note Wraps `fread()` with error handling.
 *
 * @author Lars Hoffmann
 */
#define FREAD(ptr, type, size, out) { \
    if(fread(ptr, sizeof(type), size, out)!=size) \
      ERRMSG("Error while reading!"); \
  }

/**
 * @brief Write binary data to a file.
 *
 * Writes `size` elements of type `type` from the memory pointed to by `ptr`
 * to the file stream `out`. If the number of elements written differs
 * from `size`, an error message is printed and the program terminates.
 *
 * @param[in] ptr Pointer to memory buffer containing data to write.
 * @param[in] type Type of each element to write.
 * @param[in] size Number of elements to write.
 * @param[in,out] out File stream opened for writing.
 *
 * @see FREAD
 *
 * @note Wraps `fwrite()` with error handling.
 *
 * @author Lars Hoffmann
 */
#define FWRITE(ptr, type, size, out) { \
    if(fwrite(ptr, sizeof(type), size, out)!=size) \
      ERRMSG("Error while writing!"); \
  }

/**
 * @brief Compute linear interpolation.
 *
 * Performs simple linear interpolation between two known points
 * (x₀, y₀) and (x₁, y₁) to estimate the value of y at a given x.
 *
 * @param[in] x0 Lower x-value.
 * @param[in] y0 Function value at x₀.
 * @param[in] x1 Upper x-value.
 * @param[in] y1 Function value at x₁.
 * @param[in] x Interpolation point.
 *
 * @return Linearly interpolated y-value at x.
 *
 * @see LOGX, LOGY
 *
 * @author Lars Hoffmann
 */
#define LIN(x0, y0, x1, y1, x) \
  ((y0)+((y1)-(y0))/((x1)-(x0))*((x)-(x0)))

/**
 * @brief Compute logarithmic interpolation in y.
 *
 * Performs interpolation assuming exponential variation in the y-axis
 * (logarithmic in y). If y₁/y₀ is nonpositive, reverts to linear interpolation.
 *
 * @param[in] x0 Lower x-value.
 * @param[in] y0 Function value at x₀.
 * @param[in] x1 Upper x-value.
 * @param[in] y1 Function value at x₁.
 * @param[in] x Interpolation point.
 *
 * @return Interpolated y-value at x.
 *
 * @see LIN, LOGX
 *
 * @author Lars Hoffmann
 */
#define LOGY(x0, y0, x1, y1, x) \
  (((y1)/(y0)>0) \
   ? ((y0)*exp(log((y1)/(y0))/((x1)-(x0))*((x)-(x0)))) \
   : LIN(x0, y0, x1, y1, x))

/**
 * @brief Determine the maximum of two values.
 *
 * Returns the greater of two scalar values.
 *
 * @param[in] a First value.
 * @param[in] b Second value.
 *
 * @return Maximum of a and b.
 *
 * @see MIN
 *
 * @note Both arguments are evaluated multiple times; avoid side effects.
 *
 * @author Lars Hoffmann
 */
#define MAX(a,b) (((a)>(b))?(a):(b))

/**
 * @brief Determine the minimum of two values.
 *
 * Returns the smaller of two scalar values.
 *
 * @param[in] a First value.
 * @param[in] b Second value.
 *
 * @return Minimum of a and b.
 *
 * @see MAX
 *
 * @note Both arguments are evaluated multiple times; avoid side effects.
 *
 * @author Lars Hoffmann
 */
#define MIN(a,b) (((a)<(b))?(a):(b))

/**
 * @brief Execute a NetCDF command and check for errors.
 *
 * This macro executes a NetCDF command and checks the result. If the
 * result indicates an error, it prints the error message using
 * ERRMSG.
 *
 * @param cmd NetCDF command to execute.
 *
 * @author Lars Hoffmann
 */
#define NC(cmd) {				     \
    int nc_result=(cmd);			     \
    if(nc_result!=NC_NOERR)			     \
      ERRMSG("%s", nc_strerror(nc_result));	     \
  }

/**
 * @brief Define a NetCDF variable with attributes.
 *
 * This macro defines a NetCDF variable with the specified name, data
 * type, dimensions, long name, and units. It also sets the
 * `long_name` and `units` attributes for the variable.
 * It enables compression and quantizatio of the data.
 *
 * @param varname Name of the variable.
 * @param type Data type of the variable.
 * @param ndims Number of dimensions for the variable.
 * @param dims Array of dimension IDs.
 * @param long_name Long name of the variable.
 * @param units Units of the variable.
 * @param level zlib compression level (0 = off).
 * @param quant Number of digits for quantization (0 = off).
 *
 * @note To enable ZSTD compression, replace `nc_def_var_deflate()` by
 * `nc_def_var_filter()` below. Use dynamic linking, static linking does not work.
 * Set environment variable `HDF5_PLUGIN_PATH` to `./libs/build/share/netcdf-plugins/`.
 *
 * @author Lars Hoffmann
 */
#define NC_DEF_VAR(varname, type, ndims, dims, long_name, units, level, quant) { \
    NC(nc_def_var(ncid, varname, type, ndims, dims, &varid));		\
    NC(nc_put_att_text(ncid, varid, "long_name", strnlen(long_name, LEN), long_name)); \
    NC(nc_put_att_text(ncid, varid, "units", strnlen(units, LEN), units)); \
    if((quant) > 0)							\
      NC(nc_def_var_quantize(ncid, varid, NC_QUANTIZE_GRANULARBR, quant)); \
    if((level) != 0) {							\
      NC(nc_def_var_deflate(ncid, varid, 1, 1, level));			\
      /* unsigned int ulevel = (unsigned int)level; */			\
      /* NC(nc_def_var_filter(ncid, varid, 32015, 1, (unsigned int[]){ulevel})); */ \
    }									\
  }

/**
 * @brief Retrieve a double-precision variable from a NetCDF file.
 *
 * This macro retrieves a double-precision variable from a NetCDF
 * file. It first checks if the variable exists in the file and then
 * reads its data into the specified pointer. If the `force` parameter
 * is set to true, it forces the retrieval of the variable, raising an
 * error if the variable does not exist.  If `force` is false, it
 * retrieves the variable if it exists and issues a warning if it does
 * not.
 *
 * @param varname Name of the variable to retrieve.
 * @param ptr Pointer to the memory location where the data will be stored.
 * @param force Boolean flag indicating whether to force retrieval (true) or not (false).
 *
 * @author Lars Hoffmann
 */
#define NC_GET_DOUBLE(varname, ptr, force) {			\
    if(force) {							\
      NC(nc_inq_varid(ncid, varname, &varid));			\
      NC(nc_get_var_double(ncid, varid, ptr));			\
    } else {							\
      if(nc_inq_varid(ncid, varname, &varid) == NC_NOERR) {	\
	NC(nc_get_var_double(ncid, varid, ptr));		\
      } else							\
	WARN("netCDF variable %s is missing!", varname);	\
    }								\
  }

/**
 * @brief Inquire the length of a dimension in a NetCDF file.
 *
 * This macro retrieves the length of a specified dimension from a
 * NetCDF file.  It checks if the length of the dimension is within a
 * specified range and assigns the length to the provided pointer. If
 * the length is outside the specified range, an error message is
 * raised.
 *
 * @param dimname Name of the dimension to inquire.
 * @param ptr Pointer to an integer where the dimension length will be stored.
 * @param min Minimum acceptable length for the dimension.
 * @param max Maximum acceptable length for the dimension.
 * @param check Flag to check bounds. Set to 1 for bounds check.
 *
 * @author Lars Hoffmann
 * @author Jan Clemens
 */
#define NC_INQ_DIM(dimname, ptr, min, max, check) {       \
    int dimid; size_t naux;				  \
    NC(nc_inq_dimid(ncid, dimname, &dimid));		  \
    NC(nc_inq_dimlen(ncid, dimid, &naux));		  \
    *ptr = (int)naux;                                     \
    if (check)		                                  \
      if ((*ptr) < (min) || (*ptr) > (max))		  \
        ERRMSG("Dimension %s is out of range!", dimname); \
  }

/**
 * @brief Write double precision data to a NetCDF variable.
 *
 * This macro writes data to a specified NetCDF variable. It can
 * handle both full variable writes and hyperslab writes depending on
 * the `hyperslab` parameter. If `hyperslab` is true, the data is
 * written as a hyperslab; otherwise, the entire variable is written.
 *
 * @param varname Name of the NetCDF variable to write to.
 * @param ptr Pointer to the data to be written.
 * @param hyperslab Boolean indicating whether to write the data as a hyperslab.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_DOUBLE(varname, ptr, hyperslab) {		\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_double(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_double(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Write a float array to a NetCDF file.
 *
 * This macro writes a float array to a specified variable in a NetCDF
 * file.  Depending on the value of the hyperslab parameter, the data
 * can be written as a hyperslab or as a whole variable.
 *
 * @param varname Name of the variable to which the float array will be written.
 * @param ptr Pointer to the float array to be written.
 * @param hyperslab Boolean flag indicating if the data should be written as a hyperslab. 
 *        - If true, the data will be written as a hyperslab using the start and count arrays.
 *        - If false, the data will be written to the entire variable.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_FLOAT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_float(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_float(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Write integer data to a NetCDF variable.
 *
 * This macro writes data to a specified NetCDF variable. It can
 * handle both full variable writes and hyperslab writes depending on
 * the `hyperslab` parameter. If `hyperslab` is true, the data is
 * written as a hyperslab; otherwise, the entire variable is written.
 *
 * @param varname Name of the NetCDF variable to write to.
 * @param ptr Pointer to the data to be written.
 * @param hyperslab Boolean indicating whether to write the data as a hyperslab.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_INT(varname, ptr, hyperslab) {			\
    NC(nc_inq_varid(ncid, varname, &varid));			\
    if(hyperslab) {						\
      NC(nc_put_vara_int(ncid, varid, start, count, ptr));	\
    } else {							\
      NC(nc_put_var_int(ncid, varid, ptr));			\
    }								\
  }

/**
 * @brief Add a text attribute to a NetCDF variable.
 *
 * This macro adds a text attribute to a specified NetCDF variable. It
 * first retrieves the variable ID using its name, then it attaches
 * the text attribute to the variable.
 *
 * @param varname Name of the NetCDF variable to which the attribute will be added.
 * @param attname Name of the attribute to be added.
 * @param text Text of the attribute to be added.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_ATT(varname, attname, text) {				\
    NC(nc_inq_varid(ncid, varname, &varid));				\
    NC(nc_put_att_text(ncid, varid, attname, strnlen(text, LEN), text)); \
  }

/**
 * @brief Add a global text attribute to a NetCDF file.
 *
 * This macro adds a text attribute to the global attributes of a
 * NetCDF file.  It directly attaches the attribute to the file,
 * rather than to a specific variable.
 *
 * @param attname Name of the global attribute to be added.
 * @param text Text of the attribute to be added.
 *
 * @author Lars Hoffmann
 */
#define NC_PUT_ATT_GLOBAL(attname, text)				\
  NC(nc_put_att_text(ncid, NC_GLOBAL, attname, strnlen(text, LEN), text));

/**
 * @brief Convert noise-equivalent spectral radiance (NESR) to
 *        noise-equivalent delta temperature (NEDT).
 *
 * Computes the temperature noise level corresponding to a given NESR at a
 * background radiance, using the inverse Planck function. The NEDT expresses
 * radiometric sensitivity in temperature units, assuming blackbody emission.
 *
 * @param[in] t_bg  Background (scene) brightness temperature [K].
 * @param[in] nesr  Noise-equivalent spectral radiance
 *                  [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 * @param[in] nu    Wavenumber [cm⁻¹].
 *
 * @return Noise-equivalent delta temperature [K].
 *
 * @note This macro computes
 *       \f[
 *          \mathrm{NEDT} = T_b(L_\nu(T_\mathrm{bg}) + \mathrm{NESR})
 *                        - T_\mathrm{bg}
 *       \f]
 *       where \f$T_b\f$ is the brightness temperature from the inverse
 *       Planck function.
 *
 * @warning All units must be consistent with `PLANCK` and `BRIGHT`. Mixing
 *          radiance units per meter instead of per cm⁻¹ will yield incorrect
 *          results by orders of magnitude.
 *
 * @see NESR, PLANCK, BRIGHT, C1, C2
 */
#define NEDT(t_bg, nesr, nu)					\
  (BRIGHT(PLANCK((t_bg), (nu)) + (nesr), (nu)) - (t_bg))

/**
 * @brief Convert noise-equivalent delta temperature (NEDT) to
 *        noise-equivalent spectral radiance (NESR).
 *
 * Computes the radiance difference corresponding to a temperature noise level
 * at a given background temperature and wavenumber, using Planck’s law.
 * The NESR quantifies the radiometric sensitivity in spectral radiance units.
 *
 * @param[in] t_bg  Background (scene) brightness temperature [K].
 * @param[in] nedt  Noise-equivalent delta temperature [K].
 * @param[in] nu    Wavenumber [cm⁻¹].
 *
 * @return Noise-equivalent spectral radiance [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 *
 * @note This macro computes
 *       \f[
 *          \mathrm{NESR} = L_\nu(T_\mathrm{bg} + \mathrm{NEDT})
 *                        - L_\nu(T_\mathrm{bg})
 *       \f]
 *       where \f$L_\nu\f$ is spectral radiance from Planck’s law.
 *
 * @warning All units must be consistent with `PLANCK`. In particular,
 *          `nu` must be given in cm⁻¹.
 *
 * @see PLANCK, NEDT, BRIGHT
 */
#define NESR(t_bg, nedt, nu)					\
  (PLANCK((t_bg) + (nedt), (nu)) - PLANCK((t_bg), (nu)))

/**
 * @brief Compute the norm (magnitude) of a 3D vector.
 *
 * Computes the Euclidean norm using the dot product:
 * \f$ |a| = \sqrt{a \cdot a} \f$.
 *
 * @param[in] a Input vector (array of length 3).
 *
 * @return Magnitude (norm) of vector a.
 *
 * @see DOTP
 *
 * @author Lars Hoffmann
 */
#define NORM(a) sqrt(DOTP(a, a))

/**
 * @brief Compute spectral radiance using Planck’s law.
 *
 * Evaluates Planck’s blackbody spectral radiance for a given temperature and
 * wavenumber. The expression is formulated in wavenumber space (cm⁻¹) and uses
 * the spectroscopic constants `C1` and `C2` in wavenumber units.
 *
 * @param[in] nu  Wavenumber [cm⁻¹].
 * @param[in] T   Temperature [K].
 *
 * @return Spectral radiance \f$L_\nu\f$ in
 *         [W·m⁻²·sr⁻¹·(cm⁻¹)⁻¹].
 *
 * @note The Planck law implemented here is
 *       \f[
 *           L_\nu = \frac{C_{1}\,\nu^{3}}
 *                        {\exp\!\left(\frac{C_{2}\,\nu}{T}\right) - 1}
 *       \f]
 *       with \f$\nu\f$ in cm⁻¹.
 *       The \f$\nu^{3}\f$ dependence follows from the definition of `C1`
 *       for radiance per wavenumber in cm⁻¹.
 *
 * @warning The units of all inputs MUST be consistent:
 *          `nu` in cm⁻¹ and `T` in kelvin. Mixing wavenumber and
 *          frequency formulations, or using radiance per meter instead of
 *          per cm⁻¹, will produce results off by orders of magnitude.
 *
 * @see BRIGHT, C1, C2
 */
#define PLANCK(T, nu) \
  (C1 * POW3(nu) / gsl_expm1(C2 * (nu) / (T)))

/**
 * @brief Compute the square of a value.
 *
 * Returns x². Inline alternative to `pow(x,2)`.
 *
 * @param[in] x Input value.
 *
 * @return x squared.
 *
 * @author Lars Hoffmann
 */
#define POW2(x) ((x)*(x))

/**
 * @brief Compute the cube of a value.
 *
 * Returns x³. Inline alternative to `pow(x,3)`.
 *
 * @param[in] x Input value.
 *
 * @return x cubed.
 *
 * @author Lars Hoffmann
 */
#define POW3(x) ((x)*(x)*(x))

/**
 * @brief Convert radians to degrees.
 *
 * Converts an angle measured in radians to degrees using:
 * \f$ \text{degrees} = \text{radians} \times \frac{180}{\pi} \f$.
 *
 * @param[in] rad Angle in radians.
 *
 * @return Angle in degrees.
 *
 * @see DEG2RAD
 *
 * @author Lars Hoffmann
 */
#define RAD2DEG(rad) ((rad) * (180.0 / M_PI))

/**
 * @brief Compute air refractivity (n - 1).
 *
 * Approximates the refractivity of air under standard conditions using:
 * \f$ n - 1 = 7.753\times10^{-5} \frac{p}{T} \f$,
 * where p is pressure (hPa) and T is temperature (K).
 *
 * @param[in] p Pressure [hPa].
 * @param[in] T Temperature [K].
 *
 * @return Refractivity (dimensionless, n - 1).
 *
 * @author Lars Hoffmann
 */
#define REFRAC(p, T) (7.753e-05 * (p) / (T))

/**
 * @brief Log detailed statistics of an emissivity look-up table.
 *
 * This macro logs pressure, temperature, absorber column, and emissivity
 * ranges for all pressure levels of a given (channel, emitter) table.
 * Output is written at verbosity level 2 and is intended for diagnostics
 * and debugging.
 *
 * @param tbl Pointer to the look-up table structure.
 * @param id  Channel (detector) index.
 * @param ig  Emitter (trace gas) index.
 *
 * @author Lars Hoffmann
 */
#define TBL_LOG(tbl, id, ig)						\
  do {									\
    for (int ip = 0; ip < (tbl)->np[(id)][(ig)]; ip++)			\
      LOG(2,								\
          "p[%2d]= %.5e hPa | T[0:%2d]= %.2f ... %.2f K | "		\
          "u[0:%3d]= %.5e ... %.5e molec/cm^2 | "			\
          "eps[0:%3d]= %.5e ... %.5e",					\
          (ip),								\
          (tbl)->p[(id)][(ig)][(ip)],					\
          (tbl)->nt[(id)][(ig)][(ip)] - 1,				\
          (tbl)->t[(id)][(ig)][(ip)][0],				\
          (tbl)->t[(id)][(ig)][(ip)]					\
	  [(tbl)->nt[(id)][(ig)][(ip)] - 1],                            \
          (tbl)->nu[(id)][(ig)][(ip)][0] - 1,				\
          exp((tbl)->logu[(id)][(ig)][(ip)][0][0]),			\
          exp((tbl)->logu[(id)][(ig)][(ip)][0]				\
	      [(tbl)->nu[(id)][(ig)][(ip)][0] - 1]),			\
          (tbl)->nu[(id)][(ig)][(ip)][0] - 1,				\
          exp((tbl)->logeps[(id)][(ig)][(ip)][0][0]),			\
          exp((tbl)->logeps[(id)][(ig)][(ip)][0]			\
	      [(tbl)->nu[(id)][(ig)][(ip)][0] - 1]));			\
  } while (0)

/**
 * @brief Start or stop a named timer.
 *
 * Calls the `timer()` function with contextual information (file name, function
 * name, and line number) to start or stop timing. Useful for performance profiling.
 *
 * @param[in] name Name or label of the timer.
 * @param[in] mode Operation mode (e.g., start or stop).
 *
 * @note Relies on a user-defined `timer()` function.
 *
 * @author Lars Hoffmann
 */
#define TIMER(name, mode) \
  {timer(name, __FILE__, __func__, __LINE__, mode);}

/**
 * @brief Tokenize a string and parse a variable.
 *
 * Splits a text line into tokens separated by spaces or tabs, and reads
 * a value from the first token using `sscanf()` and the provided format.
 * If tokenization or parsing fails, an error message is printed.
 *
 * @param[in,out] line Input string buffer to tokenize (modified by strtok()).
 * @param[out] tok Pointer to token string.
 * @param[in] format Format string for `sscanf()`.
 * @param[out] var Variable to store parsed value.
 *
 * @note Uses `strtok()` internally and modifies the input buffer.
 *
 * @see FREAD, FWRITE
 *
 * @author Lars Hoffmann
 */
#define TOK(line, tok, format, var) { \
    if(((tok)=strtok((line), " \t"))) { \
      if(sscanf(tok, format, &(var))!=1) continue; \
    } else ERRMSG("Error while reading!"); \
  }

/* ------------------------------------------------------------
   Log messages...
   ------------------------------------------------------------ */

/*! Level of log messages (0=none, 1=basic, 2=detailed, 3=debug). */
#ifndef LOGLEV
#define LOGLEV 2
#endif

/*!
 * \brief Print a log message with a specified logging level.
 *
 * This macro prints a formatted log message to the standard output if
 * the specified logging level meets certain conditions. The message
 * will be indented if the logging level is greater than or equal to
 * 2.
 * 
 * \param level The logging level of the message. This should be an integer value.
 * \param ... The formatted message string and its arguments, similar to printf.
 *
 * \details
 * The `LOG` macro provides a simple way to log messages with
 * different levels of importance. The message is only printed if the
 * specified `level` is less than or equal to the pre-defined `LOGLEV`
 * macro. If the `level` is greater than or equal to 2, the message is
 * preceded by two spaces for indentation.
 *
 * The macro expands to a block of code that:
 * - Checks if the `level` is greater than or equal to 2, and if so, prints two spaces.
 * - Checks if the `level` is less than or equal to `LOGLEV`, and if so, prints the
 *   formatted message followed by a newline.
 *
 * \note
 * The `LOGLEV` macro must be defined with an appropriate logging level
 * before using the `LOG` macro.
 * 
 * @author Lars Hoffmann
 */
#define LOG(level, ...) {						\
    if(level >= 2)							\
      printf("  ");							\
    if(level <= LOGLEV) {						\
      printf(__VA_ARGS__);						\
      printf("\n");							\
    }									\
  }

/*!
 * \brief Print a warning message with contextual information.
 *
 * This macro prints a formatted warning message to the standard
 * output, including the file name, function name, and line number
 * where the warning occurred. The message is then passed to the `LOG`
 * macro with a logging level of 0.
 * 
 * \param ... The formatted warning message string and its arguments, similar to printf.
 *
 * \details
 * The `WARN` macro is used to print warning messages with additional context
 * about where the warning was triggered. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual warning
 * message. This ensures that warning messages are always logged,
 * regardless of the value of `LOGLEV`.
 *
 * \note
 * The `LOG` macro must be defined before using the `WARN` macro.
 * 
 * @author Lars Hoffmann
 */
#define WARN(...) {							\
    printf("\nWarning (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
  }

/*!
 * \brief Print an error message with contextual information and terminate the program.
 *
 * This macro prints a formatted error message to the standard output,
 * including the file name, function name, and line number where the
 * error occurred. After printing the message, the program is
 * terminated with an exit status indicating failure.
 * 
 * \param ... The formatted error message string and its arguments, similar to printf.
 *
 * \details
 * The `ERRMSG` macro is used to report critical errors that require the
 * program to terminate immediately. The message includes the following
 * contextual information:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 *
 * After printing this contextual information, the macro uses the
 * `LOG` macro with a logging level of 0 to print the actual error
 * message. Finally, the program exits with a failure status
 * (`EXIT_FAILURE`).
 *
 * \note
 * The `LOG` macro must be defined before using the `ERRMSG` macro.
 * 
 * @author Lars Hoffmann
 */
#define ERRMSG(...) {							\
    printf("\nError (%s, %s, l%d): ", __FILE__, __func__, __LINE__);	\
    LOG(0, __VA_ARGS__);						\
    exit(EXIT_FAILURE);							\
  }

/*!
 * \brief Print the value of a variable with contextual information.
 *
 * This macro prints the value of a variable to the standard output,
 * including the file name, function name, and line number where the
 * macro is called. The output also includes the variable's name and
 * value in a formatted string.
 * 
 * \param format The format string used to print the variable's value, similar to printf.
 * \param var The variable to be printed.
 *
 * \details
 * The `PRINT` macro is used to output the value of a variable along with
 * additional context about where the macro is called. The message includes:
 * - The name of the source file where the macro is called (`__FILE__`).
 * - The name of the function where the macro is called (`__func__`).
 * - The line number in the source file where the macro is called (`__LINE__`).
 * - The name of the variable being printed (`#var`).
 * - The value of the variable, formatted according to the provided format string (`format`).
 *
 * This macro is particularly useful for debugging purposes, providing
 * a convenient way to trace variable values and their locations in
 * the code.
 *
 * \note
 * The format string must be compatible with the type of the variable being printed.
 * 
 * @author Lars Hoffmann
 */
#define PRINT(format, var)						\
  printf("Print (%s, %s, l%d): %s= "format"\n",				\
	 __FILE__, __func__, __LINE__, #var, var);

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/**
 * @brief Atmospheric profile data.
 *
 * Holds one vertical atmospheric column including geolocation,
 * thermodynamic, cloud, and surface properties for radiative-transfer
 * calculations.
 */
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

  /*! Surface temperature [K]. */
  double sft;

  /*! Surface emissivity. */
  double sfeps[NSF];

} atm_t;

/**
 * @brief Control parameters.
 * 
 * This structure contains all control parameters used by the JURASSIC
 * model. The struct is used to collect and to easily pass the control
 * parameters on to the various functions.
 */
typedef struct {

  /*! Number of emitters. */
  int ng;

  /*! Name of each emitter. */
  char emitter[NG][LEN];

  /*! Emitter index of CO2. */
  int ig_co2;

  /*! Emitter index of H2O. */
  int ig_h2o;

  /*! Emitter index of N2. */
  int ig_n2;

  /*! Emitter index of O2. */
  int ig_o2;

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

  /*! Look-up table file format (1=ASCII, 2=binary, 3=netCDF). */
  int tblfmt;

  /*! Atmospheric data file format (1=ASCII, 2=binary, 3=netCDF). */
  int atmfmt;

  /*! Observation data file format (1=ASCII, 2=binary, 3=netCDF). */
  int obsfmt;

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

  /*! Field-of-view vertical distance [km]. */
  double fov_dz[NSHAPE];

  /*! Field-of-view weighting factor. */
  double fov_w[NSHAPE];

  /*! Field-of-view number of data points. */
  int fov_n;

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

/**
 * @brief Line-of-sight data.
 *
 * Contains all quantities along a ray path used for radiative-transfer
 * calculations, including geometry, thermodynamic state, gas and
 * extinction profiles, and precomputed optical parameters.
 */
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

/**
 * @brief Observation geometry and radiance data.
 *
 * Stores viewing geometry and radiative quantities for multiple ray paths.
 * Each path represents a line of sight between observer and tangent point,
 * including associated time and location data.
 */
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

/**
 * @brief Retrieval control parameters.
 *
 * The `ret_t` structure holds all parameters controlling the iterative
 * inversion process (retrieval), including convergence criteria,
 * kernel matrix updates, and uncertainty specifications for each
 * retrieved quantity (pressure, temperature, trace gas concentrations,
 * extinction, cloud, and surface parameters).
 */
typedef struct {

  /*! Working directory. */
  char dir[LEN];

  /*! Re-computation of kernel matrix (number of iterations). */
  int kernel_recomp;

  /*! Maximum number of iterations. */
  int conv_itmax;

  /*! Minimum normalized step size in state space. */
  double conv_dmin;

  /*! Carry out error analysis (0=no, 1=yes). */
  int err_ana;

  /*! Forward model error [%]. */
  double err_formod[ND];

  /*! Noise error [W/(m^2 sr cm^-1)]. */
  double err_noise[ND];

  /*! Pressure error [%]. */
  double err_press;

  /*! Vertical correlation length for pressure error [km]. */
  double err_press_cz;

  /*! Horizontal correlation length for pressure error [km]. */
  double err_press_ch;

  /*! Temperature error [K]. */
  double err_temp;

  /*! Vertical correlation length for temperature error [km]. */
  double err_temp_cz;

  /*! Horizontal correlation length for temperature error [km]. */
  double err_temp_ch;

  /*! Volume mixing ratio error [%]. */
  double err_q[NG];

  /*! Vertical correlation length for volume mixing ratio error [km]. */
  double err_q_cz[NG];

  /*! Horizontal correlation length for volume mixing ratio error [km]. */
  double err_q_ch[NG];

  /*! Extinction error [km^-1]. */
  double err_k[NW];

  /*! Vertical correlation length for extinction error [km]. */
  double err_k_cz[NW];

  /*! Horizontal correlation length for extinction error [km]. */
  double err_k_ch[NW];

  /*! Cloud height error [km]. */
  double err_clz;

  /*! Cloud depth error [km]. */
  double err_cldz;

  /*! Cloud extinction error [km^-1]. */
  double err_clk[NCL];

  /*! Surface temperature error [K]. */
  double err_sft;

  /*! Surface emissivity error. */
  double err_sfeps[NSF];

} ret_t;

/**
 * @brief Emissivity look-up tables.
 *
 * Stores precomputed emissivity and source-function data for
 * different gases, spectral channels, and emitter column densities.
 */
typedef struct {

  /*! Number of pressure levels. */
  int np[ND][NG];

  /*! Number of temperatures. */
  int nt[ND][NG][TBLNP];

  /*! Number of column densities. */
  int nu[ND][NG][TBLNP][TBLNT];

  /*! Pressure [hPa]. */
  double p[ND][NG][TBLNP];

  /*! Log-pressure [hPa]. */
  double lnp[ND][NG][TBLNP];

  /*! Temperature [K]. */
  double t[ND][NG][TBLNP][TBLNT];

  /*! Logarithm of column density [molecules/cm^2]. */
  float *logu[ND][NG][TBLNP][TBLNT];

  /*! Logarithm of emissivity. */
  float *logeps[ND][NG][TBLNP][TBLNT];

  /*! Filter function number of spectral grid points. */
  int filt_n[ND];

  /*! Filter function spectral grid points [cm^-1]. */
  double filt_nu[ND][NSHAPE];

  /*! Filter function values. */
  double filt_f[ND][NSHAPE];

  /*! Source function temperature [K]. */
  double st[TBLNS];

  /*! Source function radiance [W/(m^2 sr cm^-1)]. */
  double sr[TBLNS][ND];

} tbl_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/**
 * @brief Analyze averaging kernel (AVK) matrix for retrieval diagnostics.
 *
 * Decomposes and evaluates the averaging kernel matrix \f$\mathbf{A}\f$
 * to quantify the retrieval sensitivity and vertical resolution for each
 * retrieved quantity (pressure, temperature, trace gases, extinction, etc.).
 * The results are written to diagnostic atmospheric files for visualization.
 *
 * @param[in]  ret   Retrieval control and configuration structure (`ret_t`).
 * @param[in]  ctl   Global control structure (`ctl_t`) defining retrieval setup and quantities.
 * @param[in]  atm   Retrieved atmospheric state structure (`atm_t`).
 * @param[in]  iqa   Array mapping state vector indices to physical quantities (e.g., IDXP, IDXT, IDXQ(...)).
 * @param[in]  ipa   Array mapping state vector indices to atmospheric grid points.
 * @param[in]  avk   Averaging kernel matrix \f$\mathbf{A}\f$ (`gsl_matrix`, n×n).
 *
 * @details
 * The averaging kernel matrix \f$\mathbf{A}\f$ describes the sensitivity of the
 * retrieved state \f$\hat{x}\f$ to the true state \f$x\f$:
 * \f[
 * \hat{x} - x_a = \mathbf{A} (x - x_a)
 * \f]
 * where \f$x_a\f$ is the a priori state.  
 * This function separates and evaluates:
 * - **Contribution profiles:** diagonal dominance or total response (sensitivity).
 * - **Resolution profiles:** width or spread of the averaging kernels.
 *
 * Internally, this function:
 * 1. Identifies submatrices of \f$\mathbf{A}\f$ corresponding to each
 *    physical quantity (pressure, temperature, VMRs, extinction, clouds, surface).
 * 2. Calls `analyze_avk_quantity()` for each retrieved parameter type to compute
 *    quantitative measures of contribution and resolution.
 * 3. Writes the results as atmospheric profiles:
 *    - `atm_cont.tab` — contribution functions (sensitivity).
 *    - `atm_res.tab` — resolution functions (vertical response width).
 *
 * @see analyze_avk_quantity, write_atm, set_cov_apr, set_cov_meas
 *
 * @note
 * - Submatrices are identified by matching the quantity index `iqa[i]` to the
 *   quantity constants (e.g., `IDXT`, `IDXQ(ig)`, etc.).
 * - Cloud and surface quantities are treated as scalar elements.
 * - Output files are written to the retrieval working directory (`ret->dir`).
 *
 * @warning
 * - The averaging kernel must be fully computed and dimensionally consistent
 *   with the state vector before calling this function.
 * - File writing will overwrite existing diagnostics in `ret->dir`.
 *
 * @author
 * Lars Hoffmann
 */
void analyze_avk(
  const ret_t * ret,
  const ctl_t * ctl,
  const atm_t * atm,
  const int *iqa,
  const int *ipa,
  const gsl_matrix * avk);

/**
 * @brief Analyze averaging kernel submatrix for a specific retrieved quantity.
 *
 * Computes the contribution (sensitivity) and resolution (information density)
 * profiles for a given physical quantity from its corresponding submatrix
 * of the full averaging kernel matrix \f$\mathbf{A}\f$.
 *
 * @param[in]  avk   Averaging kernel matrix \f$\mathbf{A}\f$ (`gsl_matrix`, n×n),
 *                   describing the sensitivity of the retrieved state to the true state.
 * @param[in]  iq    Quantity index identifier (e.g., `IDXP`, `IDXT`, `IDXQ(ig)`,
 *                   `IDXK(iw)`, `IDXCLZ`, etc.).
 * @param[in]  ipa   Array mapping state vector indices to atmospheric grid indices.
 * @param[in]  n0    Array of starting indices for each quantity sub-block in the state vector.
 * @param[in]  n1    Array of lengths (number of elements) for each quantity sub-block.
 * @param[out] cont  Array of contribution values, representing the total
 *                   sensitivity or response of the retrieval at each grid point.
 * @param[out] res   Array of resolution measures, representing the vertical
 *                   information density or averaging kernel width at each grid point.
 *
 * @details
 * For a given quantity \f$q\f$, the function extracts its corresponding
 * block \f$\mathbf{A}_q\f$ from the full averaging kernel matrix and computes:
 *
 * - **Contribution function**
 *   \f[
 *   \text{cont}_i = \sum_j A_{ij}^{(q)}
 *   \f]
 *   giving the total response of retrieved element \f$i\f$ to perturbations
 *   in all true elements \f$j\f$ of the same quantity.
 *
 * - **Resolution (information density)**
 *   \f[
 *   \text{res}_i = \frac{1}{A_{ii}^{(q)}}
 *   \f]
 *   approximating the degree of vertical resolution or smoothing at level \f$i\f$.
 *
 * The results are stored in the provided arrays `cont[]` and `res[]`,
 * indexed according to the atmospheric profile index via `ipa[]`.
 *
 * @see analyze_avk, write_atm, set_cov_apr
 *
 * @note
 * - Only elements associated with the specified quantity index (`iq`) are processed.
 * - The function assumes the averaging kernel matrix is square and symmetric.
 * - The contribution and resolution arrays must be preallocated to at least
 *   the number of atmospheric grid points (`atm->np`).
 *
 * @warning
 * - If the diagonal element `A_ii` is zero or near-zero, the corresponding
 *   resolution value may be undefined or numerically unstable.
 * - Input indices `n0[iq]` and `n1[iq]` must have been computed by `analyze_avk()`.
 *
 * @author
 * Lars Hoffmann
 */
void analyze_avk_quantity(
  const gsl_matrix * avk,
  const int iq,
  const int *ipa,
  const size_t *n0,
  const size_t *n1,
  double *cont,
  double *res);

/**
 * @brief Convert atmospheric data to state vector elements.
 *
 * Extracts selected quantities from an atmospheric profile (@ref atm_t)
 * according to retrieval settings in @ref ctl_t, and appends them to
 * the state vector @p x. For each included quantity, the function also
 * stores its quantity index (@p iqa) and profile index (@p ipa).
 *
 * The function respects retrieval altitude limits defined in @p ctl
 * (e.g., `retp_zmin/zmax`, `rett_zmin/zmax`, etc.) and includes only
 * variables flagged for retrieval (e.g., `ret_clz`, `ret_sft`, etc.).
 *
 * @param[in]  ctl  Control settings defining retrieval configuration and limits.
 * @param[in]  atm  Atmospheric profile data to extract from.
 * @param[out] x    GSL vector to store state-vector elements.
 * @param[out] iqa  Quantity index array corresponding to elements in @p x.
 * @param[out] ipa  Profile index array corresponding to elements in @p x.
 *
 * @return Number of elements written to the state vector.
 *
 * @note Internally calls @ref atm2x_help() to append individual values.
 *
 * @see atm_t, ctl_t, atm2x_help
 *
 * @author Lars Hoffmann
 */
size_t atm2x(
  const ctl_t * ctl,
  const atm_t * atm,
  gsl_vector * x,
  int *iqa,
  int *ipa);

/**
 * @brief Append a single atmospheric value to the state vector.
 *
 * Helper routine for @ref atm2x(). Inserts one scalar value and its
 * corresponding quantity and profile indices into the state vector and
 * tracking arrays, then increments the element counter.
 *
 * @param[in]  value     Value to add to the state vector.
 * @param[in]  value_iqa Quantity index (e.g., @ref IDXP, @ref IDXT, etc.).
 * @param[in]  value_ip  Profile index within the atmospheric profile.
 * @param[out] x         GSL vector containing state-vector elements (may be NULL).
 * @param[out] iqa       Quantity index array corresponding to @p x (may be NULL).
 * @param[out] ipa       Profile index array corresponding to @p x (may be NULL).
 * @param[in,out] n      Current number of elements in the state vector; incremented on return.
 *
 * @note This function performs no range checking and assumes valid array bounds.
 *
 * @see atm2x
 *
 * @author Lars Hoffmann
 */
void atm2x_help(
  const double value,
  const int value_iqa,
  const int value_ip,
  gsl_vector * x,
  int *iqa,
  int *ipa,
  size_t *n);

/**
 * @brief Converts Cartesian coordinates to geographic coordinates.
 *
 * This function converts a point from Cartesian coordinates (x, y, z)
 * to geographic coordinates (longitude, latitude, and altitude).  It
 * uses the spherical Earth approximation for the conversion.
 *
 * @param x Pointer to an array containing the Cartesian coordinates (x, y, z) in kilometers.
 * @param z Pointer to a double where the computed altitude (above the reference ellipsoid) will be stored, in kilometers.
 * @param lon Pointer to a double where the computed longitude (in degrees) will be stored.
 * @param lat Pointer to a double where the computed latitude (in degrees) will be stored.
 *
 * @author Lars Hoffmann
 */
void cart2geo(
  const double *x,
  double *z,
  double *lon,
  double *lat);

/**
 * @brief Initializes atmospheric climatology profiles.
 *
 * This function populates the atmospheric state (`atm`) with standard
 * climatological profiles of pressure, temperature, and trace gas
 * concentrations (e.g., H2O, CH4, CO, O3, etc.) as a function of altitude.
 * The profiles are based on reference climatological datasets and are used
 * for atmospheric modeling and radiative transfer calculations.
 *
 * @param[in]  ctl  Control parameters structure.
 * @param[out] atm  Atmospheric state structure to be populated with
 *                  climatological data.
 *
 * @author Lars Hoffmann
 */
void climatology(
  const ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Calculates the cosine of the solar zenith angle.
 *
 * This function computes the cosine of the solar zenith angle (SZA), which describes
 * the angle between the local zenith (straight up) and the line connecting the
 * observer to the center of the Sun. The cosine of the SZA is often used directly
 * in radiative transfer and photochemical calculations to avoid unnecessary use
 * of trigonometric inverse functions.
 *
 * @param sec Seconds elapsed since 2000-01-01T12:00Z.
 * @param lon Observer's longitude in degrees.
 * @param lat Observer's latitude in degrees.
 * @return The cosine of the solar zenith angle (dimensionless, range [-1, 1]).
 *
 * The cosine of the solar zenith angle is computed based on the observer's position
 * (longitude and latitude) and the specified time in seconds elapsed since
 * 2000-01-01T12:00Z.
 *
 * @note The input longitude and latitude must be specified in degrees.
 *
 * @see acos() — can be used to convert the returned value to the solar zenith angle in radians if needed.
 *
 * @author Lars Hoffmann
 */
double cos_sza(
  const double sec,
  const double lon,
  const double lat);

/**
 * @brief Compute the normalized quadratic cost function for optimal estimation.
 *
 * Evaluates the cost function
 * \f[
 * J = \frac{1}{m} \left[ (\mathbf{y} - \mathbf{y_a})^T
 * \mathbf{S_\epsilon^{-1}} (\mathbf{y} - \mathbf{y_a})
 * + (\mathbf{x} - \mathbf{x_a})^T \mathbf{S_a^{-1}} (\mathbf{x} - \mathbf{x_a}) \right]
 * \f]
 * where \f$\mathbf{dx} = \mathbf{x} - \mathbf{x_a}\f$ and
 * \f$\mathbf{dy} = \mathbf{y} - \mathbf{y_a}\f$ represent the deviations
 * from a priori state and measurement vectors, respectively.
 *
 * @param[in] dx           State deviation vector (`x - x_a`).
 * @param[in] dy           Measurement deviation vector (`y - y_a`).
 * @param[in] s_a_inv      Inverse of the a priori covariance matrix (\f$\mathbf{S_a^{-1}}\f$).
 * @param[in] sig_eps_inv  Vector of inverse measurement uncertainties
 *                         (\f$\mathbf{S_\epsilon^{-1/2}}\f$ diagonal elements).
 *
 * @return Normalized cost function value (dimensionless).
 *
 * @details
 * - Implements the standard *Optimal Estimation Method* (Rodgers, 2000) cost function.
 * - The first term (\f$\chi^2_m\f$) measures the weighted measurement residuals:
 *   \f$\chi^2_m = \sum_i (dy_i / \sigma_{\epsilon,i})^2\f$.
 * - The second term (\f$\chi^2_a\f$) measures the deviation from the a priori state:
 *   \f$\chi^2_a = dx^T \mathbf{S_a^{-1}} dx\f$.
 * - The result is normalized by the number of measurements \f$m\f$ for scale consistency.
 *
 * @see retrieval, s_a_inv, sig_eps_inv
 *
 * @note
 * - Assumes diagonal measurement error covariance (independent measurements).
 * - The covariance matrices must be positive-definite and properly dimensioned.
 * - The function does **not** modify the input vectors.
 *
 * @author
 * Lars Hoffmann
 */
double cost_function(
  const gsl_vector * dx,
  const gsl_vector * dy,
  const gsl_matrix * s_a_inv,
  const gsl_vector * sig_eps_inv);

/**
 * @brief Compute carbon dioxide continuum (optical depth).
 *
 * @author Lars Hoffmann
 */
double ctmco2(
  const double nu,
  const double p,
  const double t,
  const double u);

/**
 * @brief Compute water vapor continuum (optical depth).
 *
 * @author Lars Hoffmann
 */
double ctmh2o(
  const double nu,
  const double p,
  const double t,
  const double q,
  const double u);

/**
 * @brief Compute N₂ collision-induced absorption coefficient.
 *
 * Calculates the nitrogen (N₂) absorption coefficient due to
 * collision-induced absorption (CIA) near the 4.3 µm CO₂ band using
 * tabulated laboratory data for the absorption strength (B) and
 * temperature exponent (β).
 *
 * The function linearly interpolates B and β as a function of
 * wavenumber, then applies a temperature- and pressure-dependent
 * scaling relation to compute the absorption coefficient.
 *
 * @param[in] nu  Wavenumber [cm⁻¹].
 * @param[in] p   Pressure [hPa].
 * @param[in] t   Temperature [K].
 *
 * @return N₂ absorption coefficient [km⁻¹].
 *
 * @note Valid for approximately 2120–2600 cm⁻¹ (4.0–4.7 µm) where
 *       N₂–N₂ and N₂–O₂ CIA dominates. Returns zero outside the tabulated
 *       range.
 *
 * @see locate_reg, LIN, P0, N2
 *
 * @par Reference
 *   Lafferty et al., J. Quant. Spectrosc. Radiat. Transf., 68, 473–479 (2001)
 *
 * @author Lars Hoffmann
 */
double ctmn2(
  const double nu,
  const double p,
  const double t);

/**
 * @brief Compute O₂ collision-induced absorption coefficient.
 *
 * Calculates the molecular oxygen (O₂) absorption coefficient due to
 * collision-induced absorption (CIA) using tabulated laboratory data
 * for the absorption strength (B) and temperature exponent (β).
 *
 * The function linearly interpolates B and β as a function of wavenumber
 * and applies a pressure- and temperature-dependent scaling relation to
 * compute the absorption coefficient.
 *
 * @param[in] nu  Wavenumber [cm⁻¹].
 * @param[in] p   Pressure [hPa].
 * @param[in] t   Temperature [K].
 *
 * @return O₂ absorption coefficient [km⁻¹].
 *
 * @note Valid for approximately 1360–1800 cm⁻¹ (∼ 7.4–5.5 µm),
 *       corresponding to the O₂ CIA band. Returns zero outside the
 *       tabulated range.
 *
 * @see locate_reg, LIN, P0, O2
 *
 * @par References
 *   Greenblatt et al., J. Quant. Spectrosc. Radiat. Transf., 33, 127–140 (1985)
 *   Smith and Newnham, Appl. Opt., 39, 318–326 (2000)
 *
 * @author Lars Hoffmann
 */
double ctmo2(
  const double nu,
  const double p,
  const double t);

/**
 * @brief Copy or initialize atmospheric profile data.
 *
 * Copies all fields from one atmospheric structure (@p atm_src) to another
 * (@p atm_dest), including geolocation, thermodynamic, gas, extinction,
 * cloud, and surface parameters. If @p init is nonzero, the destination
 * fields are instead initialized to default values (zeros for most fields,
 * unity for surface emissivity).
 *
 * @param[in]  ctl       Control structure defining array dimensions.
 * @param[out] atm_dest  Destination atmospheric structure.
 * @param[in]  atm_src   Source atmospheric structure.
 * @param[in]  init      Initialization flag:
 *                       - 0 → copy data from @p atm_src
 *                       - 1 → initialize @p atm_dest to default values
 *
 * @note The number of vertical levels (`atm_src->np`) determines the
 *       amount of data copied. The function performs shallow copies
 *       using @c memcpy for efficiency.
 *
 * @see atm_t, ctl_t
 *
 * @author Lars Hoffmann
 */
void copy_atm(
  const ctl_t * ctl,
  atm_t * atm_dest,
  const atm_t * atm_src,
  const int init);

/**
 * @brief Copy or initialize observation geometry and radiance data.
 *
 * Copies all observation fields from a source structure (@p obs_src)
 * to a destination structure (@p obs_dest), including observer, view,
 * and tangent point geometry, as well as radiance and transmittance data.
 * If @p init is nonzero, radiance and transmittance values are reset
 * to zero wherever finite values are found.
 *
 * @param[in]  ctl       Control structure defining the number of channels
 *                       (@ref ctl_t::nd) and other dimensions.
 * @param[out] obs_dest  Destination observation structure.
 * @param[in]  obs_src   Source observation structure.
 * @param[in]  init      Initialization flag:
 *                       - 0 → copy data from @p obs_src  
 *                       - 1 → initialize @p obs_dest (set @p rad and @p tau to zero)
 *
 * @note The number of ray paths (`obs_src->nr`) defines the copied data size.
 *       Shallow copies are performed via @c memcpy for efficiency.
 *
 * @see obs_t, ctl_t
 * 
 * @author Lars Hoffmann
 */
void copy_obs(
  const ctl_t * ctl,
  obs_t * obs_dest,
  const obs_t * obs_src,
  const int init);

/**
 * @brief Convert a calendar date to day-of-year.
 *
 * This function computes the day-of-year (1–365 or 1–366 for leap years)
 * corresponding to a given Gregorian date specified by year, month, and day.
 *
 * Leap years are determined using the standard Gregorian rules:
 *   - Every year divisible by 4 is a leap year,
 *   - except years divisible by 100, unless also divisible by 400.
 *
 * @param year  The full year (e.g., 2025).
 * @param mon   The month of the year (1–12).
 * @param day   The day of the month (1–31).
 * @param[out] doy Pointer to an integer where the resulting day-of-year
 *                 will be stored (1–365 or 1–366).
 *
 * @note The caller must ensure that the input date is valid.
 * 
 * @author Lars Hoffmann
 */
void day2doy(
  int year,
  int mon,
  int day,
  int *doy);

/**
 * @brief Convert a day-of-year value to a calendar date.
 *
 * This function computes the month and day corresponding to a given
 * day-of-year (1–365 or 1–366 for leap years) in the Gregorian calendar.
 *
 * Leap years are determined using the standard Gregorian rules:
 *   - Every year divisible by 4 is a leap year,
 *   - except years divisible by 100, unless also divisible by 400.
 *
 * @param year  The full year (e.g., 2025).
 * @param doy   The day-of-year (1–365 or 1–366).
 * @param[out] mon Pointer to an integer where the calculated month (1–12)
 *                 will be stored.
 * @param[out] day Pointer to an integer where the calculated day-of-month
 *                 will be stored.
 *
 * @note The caller must ensure that @p doy is within the valid range for
 *       the specified year.
 * 
 * @author Lars Hoffmann
 */
void doy2day(
  int year,
  int doy,
  int *mon,
  int *day);

/**
 * @brief Find gas species index by name.
 *
 * Searches the list of emitter (gas) names defined in the control
 * structure for a case-insensitive match to the given string.
 * Returns the corresponding gas index if found, or -1 otherwise.
 *
 * @param[in] ctl      Control structure containing the list of gas emitters.
 * @param[in] emitter  Name of the gas species to search for (e.g. "H2O", "CO2").
 *
 * @return Index of the matching emitter (0 ≤ index < ctl->ng),
 *         or -1 if no match is found.
 *
 * @note Comparison is case-insensitive using @c strcasecmp().
 *
 * @see ctl_t
 * 
 * @author Lars Hoffmann
 */
int find_emitter(
  const ctl_t * ctl,
  const char *emitter);

/**
 * @brief Execute the selected forward model.
 *
 * Computes synthetic radiances or brightness temperatures for the
 * given atmospheric state and observation geometry using the selected
 * forward-model method (CGA, EGA, or RFM). The function also applies
 * hydrostatic adjustment, field-of-view convolution, and optional
 * brightness-temperature conversion.
 *
 * @param[in]  ctl  Control structure defining model settings and options.
 * @param[in]  tbl  Emissivity and source-function lookup tables.
 * @param[in,out] atm  Atmospheric profile; may be adjusted for hydrostatic balance.
 * @param[in,out] obs  Observation geometry and radiance data; populated with model output.
 *
 * @note The model type is selected via @ref ctl_t::formod:
 *       - 0 or 1 → pencil-beam models (@ref formod_pencil)  
 *       - 2 → RFM line-by-line model (@ref formod_rfm)
 *
 * @note The function preserves @p obs->rad elements marked as invalid
 *       (NaN) by applying an internal observation mask.
 *
 * @see ctl_t, atm_t, obs_t, tbl_t, formod_pencil, formod_rfm, formod_fov, hydrostatic
 * 
 * @author Lars Hoffmann
 */
void formod(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs);

/**
 * @brief Compute total extinction including gaseous continua.
 *
 * Calculates the total extinction coefficient for a given line-of-sight
 * point by summing spectrally dependent extinction and optional gaseous
 * continuum contributions (CO₂, H₂O, N₂, O₂).
 *
 * The function updates the extinction array @p beta for each spectral
 * channel based on line-by-line extinction and enabled continua flags
 * specified in @p ctl.
 *
 * @param[in]  ctl  Control structure defining model setup and continuum options.
 * @param[in]  los  Line-of-sight data containing pressure, temperature,
 *                  gas concentrations, and extinction coefficients.
 * @param[in]  ip   Index of the line-of-sight point to process.
 * @param[out] beta Array of total extinction coefficients [km⁻¹] per channel.
 *
 * @note Each continuum component is added only if its corresponding
 *       control flag (e.g., @ref ctl_t::ctm_co2, @ref ctl_t::ctm_h2o)
 *       is enabled and the gas index is valid.
 *
 * @see ctl_t, los_t, ctmco2, ctmh2o, ctmn2, ctmo2
 * 
 * @author Lars Hoffmann
 */
void formod_continua(
  const ctl_t * ctl,
  const los_t * los,
  const int ip,
  double *beta);

/**
 * @brief Apply field-of-view (FOV) convolution to modeled radiances.
 *
 * Convolves pencil-beam radiances and transmittances along each ray path
 * with the instrument field-of-view weighting function defined in the
 * control structure. This simulates finite FOV effects on measured
 * radiances and transmittances.
 *
 * @param[in]  ctl  Control structure containing FOV parameters
 *                  (offsets @ref ctl_t::fov_dz and weights @ref ctl_t::fov_w).
 * @param[in,out] obs  Observation structure; input pencil-beam data are replaced
 *                     with FOV-convolved radiances and transmittances.
 *
 * @note The convolution is skipped if @ref ctl_t::fov starts with '-'
 *       (indicating no FOV correction). Requires at least two valid
 *       altitude samples per time step.
 *
 * @throws ERRMSG if insufficient data are available for convolution.
 *
 * @author Lars Hoffmann
 */
void formod_fov(
  const ctl_t * ctl,
  obs_t * obs);

/**
 * @brief Compute line-of-sight radiances using the pencil-beam forward model.
 *
 * Simulates monochromatic radiances and transmittances along a single
 * line of sight using a layer-by-layer pencil-beam approximation.
 * The model includes gaseous absorption, continuum extinction, surface
 * emission, reflection, and optional solar illumination.
 *
 * @param[in]  ctl  Control structure defining model configuration, gas setup,
 *                  surface type, and spectral parameters.
 * @param[in]  tbl  Emissivity and source-function lookup tables.
 * @param[in]  atm  Atmospheric state containing pressure, temperature,
 *                  and gas profiles.
 * @param[in,out] obs  Observation data; updated with modeled radiances and
 *                     transmittances for the specified ray path.
 * @param[in]  ir   Index of the current ray path in @p obs.
 *
 * @note Depending on @ref ctl_t::formod, this function calls either
 *       @ref intpol_tbl_cga() (CGA) or @ref intpol_tbl_ega() (EGA)
 *       for gas absorption interpolation.  
 *       Surface effects include emission, reflection, and—if enabled—
 *       solar illumination based on the solar zenith angle.
 *
 * @see ctl_t, atm_t, obs_t, tbl_t, los_t,
 *      raytrace, formod_continua, formod_srcfunc,
 *      intpol_tbl_cga, intpol_tbl_ega, PLANCK
 *
 * @author Lars Hoffmann
 */
void formod_pencil(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const atm_t * atm,
  obs_t * obs,
  const int ir);

/**
 * @brief Forward-model radiance and transmittance with the Reference Forward Model (RFM).
 *
 * This routine provides the JURASSIC interface to the external RFM executable. It writes an
 * RFM atmospheric profile file and a per-channel RFM driver file, runs RFM, and reads back
 * the resulting radiance and transmittance spectra. The spectra are then convolved with the
 * instrument channel filter functions and stored in the observation structure.
 *
 * Filter functions are taken exclusively from the lookup-table container @p tbl
 * (i.e., @c tbl->filt_n, @c tbl->filt_nu, @c tbl->filt_f). This makes the behavior independent
 * of the lookup-table format (ASCII, binary, or netCDF) and avoids reliance on external
 * per-channel ASCII filter files.
 *
 * Viewing geometry is classified as one of:
 * - limb: uses tangent altitude
 * - nadir: path intersects the surface (secant/air-mass-factor geometry)
 * - zenith: upward-looking path exits the atmosphere at the top boundary
 *
 * Mixed geometries (e.g., limb and nadir simultaneously) are not allowed and will trigger an error.
 *
 * Limitations:
 * - Requires identical observer positions for all rays.
 * - Does not support extinction input data (atm->k must be zero everywhere).
 * - Uses temporary files in the current working directory (e.g., @c rfm.atm, @c rfm.drv,
 *   @c rad_*.asc, @c tra_*.asc) and removes them on completion.
 *
 * @param[in]  ctl  Control and configuration parameters (RFM binary path, HITRAN/XSC settings,
 *                  channel definitions, flags such as refraction and continua).
 * @param[in]  tbl  Lookup-table container providing per-channel filter functions.
 * @param[in]  atm  Atmospheric state (altitude grid, temperature, and gas profiles).
 * @param[out] obs  Observation geometry and output arrays; on return, @c obs->rad[id][ir]
 *                  and @c obs->tau[id][ir] are filled for all channels and rays.
 *
 * @pre @p tbl contains valid filter functions for all channels used (@c tbl->filt_n[id] > 0).
 * @pre All rays share identical observer position (@c obsz/obslon/obslat).
 * @pre No extinction data is present (@c atm->k == 0 for all wavelengths/levels).
 *
 * @post @c obs->rad and @c obs->tau are updated with channel-integrated radiance and
 *       transmittance computed by RFM.
 *
 * @note The surface temperature is taken as the temperature at the lowest altitude level
 *       of the atmospheric grid.
 *
 * @warning This function executes external commands via @c system(), and creates/removes files
 *          in the current working directory.
 *
 * @par Reference
 *   Dudhia, A., "The Reference Forward Model (RFM)", JQSRT 186, 243–253 (2017)
 *
 * @author Lars Hoffmann
 */
void formod_rfm(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const atm_t * atm,
  obs_t * obs);

/**
 * @brief Interpolate the source function (Planck radiance) at a given temperature.
 *
 * Computes source function values by linearly interpolating
 * between precomputed Planck radiances in the lookup table. The resulting
 * radiance spectrum corresponds to the input temperature and is stored
 * in @p src for all spectral channels.
 *
 * @param[in]  ctl  Control structure defining the number of spectral channels.
 * @param[in]  tbl  Emissivity and source-function lookup table containing
 *                  Planck radiances (@ref tbl_t::sr) and corresponding
 *                  temperatures (@ref tbl_t::st).
 * @param[in]  t    Temperature [K] for which the source function is evaluated.
 * @param[out] src  Output array of interpolated source-function values
 *                  [W·m⁻²·sr⁻¹·cm⁻¹] per spectral channel.
 *
 * @note Linear interpolation is used between the two nearest temperature
 *       grid points. The function does not extrapolate beyond the tabulated
 *       temperature range.
 *
 * @see ctl_t, tbl_t, LIN, locate_reg
 * 
 * @author Lars Hoffmann
 */
void formod_srcfunc(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const double t,
  double *src);

/**
 * @brief Converts geographic coordinates (longitude, latitude, altitude) to Cartesian coordinates.
 *
 * This function converts geographic coordinates specified by
 * longitude, latitude, and altitude into Cartesian coordinates. The
 * Earth is approximated as a sphere with radius defined by the
 * constant `RE`.
 *
 * @param z The altitude above the Earth's surface in kilometers.
 * @param lon The longitude in degrees.
 * @param lat The latitude in degrees.
 * @param x Pointer to an array of three doubles where the computed Cartesian coordinates (x, y, z) will be stored.
 *
 * The function computes the Cartesian coordinates using the given altitude, longitude, and latitude.
 * It assumes the Earth is a perfect sphere and uses the following formulas:
 * - \f$ x = (\textrm{radius}) \cos(\textrm{lat in radians}) \cos(\textrm{lon in radians}) \f$
 * - \f$ y = (\textrm{radius}) \cos(\textrm{lat in radians}) \sin(\textrm{lon in radians}) \f$
 * - \f$ z = (\textrm{radius}) \sin(\textrm{lat in radians}) \f$
 *
 * @note The constant `RE` is defined as the Earth's radius in kilometers.
 * @note Longitude and latitude should be in degrees.
 *
 * @see https://en.wikipedia.org/wiki/Geographic_coordinate_conversion
 *
 * @author Lars Hoffmann
 */
void geo2cart(
  const double z,
  const double lon,
  const double lat,
  double *x);

/**
 * @brief Adjust pressure profile using the hydrostatic equation.
 *
 * Recomputes the atmospheric pressure field to ensure hydrostatic
 * equilibrium with respect to altitude and temperature. Starting from
 * a reference altitude (@ref ctl_t::hydz), the routine integrates the
 * hydrostatic balance equation both upward and downward through the
 * profile using small interpolation steps.
 *
 * The air density is corrected for humidity using the mean molecular
 * mass of dry air and water vapor.
 *
 * @param[in]  ctl  Control structure providing model constants and
 *                  reference altitude (@ref ctl_t::hydz) and H₂O index
 *                  (@ref ctl_t::ig_h2o).
 * @param[in,out] atm  Atmospheric state; input temperatures, heights,
 *                     and humidities are used to update pressure [hPa].
 *
 * @note The integration is performed in log-pressure space assuming
 *       hydrostatic balance:
 *       \f$ \frac{dp}{dz} = -\rho g \f$,  
 *       using 20 linear substeps between adjacent levels.
 *
 * @see ctl_t, atm_t, LIN, G0, RI
 * 
 * @author Lars Hoffmann
 */
void hydrostatic(
  const ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Convert a quantity index to a descriptive name string.
 *
 * Translates a model quantity index (e.g., pressure, temperature,
 * gas mixing ratio, extinction, or surface parameter) into a
 * human-readable name. The function uses the index mapping defined
 * in the control structure to assign appropriate labels.
 *
 * @param[in]  ctl       Control structure containing gas, window,
 *                       cloud, and surface setup information.
 * @param[in]  idx       Quantity index (see @ref IDXP, @ref IDXT,
 *                       @ref IDXQ, @ref IDXK, etc.).
 * @param[out] quantity  Character buffer to receive the descriptive
 *                       quantity name (e.g., "PRESSURE",
 *                       "H2O", "CLOUD_HEIGHT", "SURFACE_EMISSIVITY_1000.0").
 *
 * @note The function writes directly to @p quantity using @c sprintf().
 *       The caller must ensure sufficient buffer size (≥ LEN).
 *
 * @see ctl_t, IDXP, IDXT, IDXQ, IDXK, IDXCLZ, IDXCLDZ, IDXCLK, IDXSFT, IDXSFEPS
 * 
 * @author Lars Hoffmann
 */
void idx2name(
  const ctl_t * ctl,
  const int idx,
  char *quantity);

/**
 * @brief Initialize source function lookup tables from emissivity data.
 *
 * This function computes channel-dependent source function lookup tables
 * based on the spectral filter functions and the Planck function. For each
 * spectral channel, the Planck radiance is integrated over the corresponding
 * filter function and stored as a function of temperature.
 *
 * The resulting source function table represents the **band-integrated
 * thermal emission** associated with the emissivity lookup tables and is
 * later used in radiative transfer calculations.
 *
 * The temperature grid is uniformly sampled between @c TMIN and @c TMAX
 * with @c TBLNS points. The spectral integration is performed on the finest
 * frequency spacing present in the filter function.
 *
 * @param[in]  ctl  Pointer to the control structure defining the number of
 *                  channels and their central wavenumbers.
 * @param[in,out] tbl  Pointer to the lookup-table structure in which the
 *                  source function tables and temperature grid are stored.
 *
 * @note The computation is parallelized over temperature grid points using
 *       OpenMP.
 *
 * @note The source function is normalized by the integral of the filter
 *       function to yield a band-averaged Planck radiance.
 *
 * @author Lars Hoffmann
 */
void init_srcfunc(
  const ctl_t * ctl,
  tbl_t * tbl);

/**
 * @brief Interpolate atmospheric state variables at a given altitude.
 *
 * Computes pressure, temperature, volume mixing ratios, and extinction
 * coefficients at the specified altitude by interpolating between adjacent
 * model levels in the atmospheric profile.
 *
 * @param[in]  ctl  Control structure defining the number of gases (@ref ctl_t::ng)
 *                  and spectral windows (@ref ctl_t::nw).
 * @param[in]  atm  Atmospheric profile providing altitude, pressure, temperature,
 *                  gas mixing ratios, and extinction data.
 * @param[in]  z    Target altitude [km].
 * @param[out] p    Interpolated pressure [hPa].
 * @param[out] t    Interpolated temperature [K].
 * @param[out] q    Interpolated gas volume mixing ratios [ppv], length @ref ctl_t::ng.
 * @param[out] k    Interpolated extinction coefficients [km⁻¹], length @ref ctl_t::nw.
 *
 * @note Pressure is interpolated logarithmically using @ref LOGY,
 *       while other quantities use linear interpolation (@ref LIN).
 *
 * @see ctl_t, atm_t, locate_irr, LIN, LOGY
 * 
 * @author Lars Hoffmann
 */
void intpol_atm(
  const ctl_t * ctl,
  const atm_t * atm,
  const double z,
  double *p,
  double *t,
  double *q,
  double *k);

/**
 * @brief Interpolate emissivities and transmittances using the
 *        Curtis–Godson approximation (CGA).
 *
 * Computes gas transmittance along a line-of-sight segment by
 * interpolating precomputed emissivity values from lookup tables.
 * The interpolation is performed in pressure, temperature, and
 * column density space using bilinear (in p, T) and logarithmic
 * (in p) interpolation.
 *
 * @param[in]  ctl        Control structure defining number of gases
 *                        (@ref ctl_t::ng) and channels (@ref ctl_t::nd).
 * @param[in]  tbl        Emissivity lookup tables (@ref tbl_t)
 *                        containing tabulated pressure, temperature,
 *                        and column density grids.
 * @param[in]  los        Line-of-sight structure providing Curtis–Godson
 *                        mean parameters and column densities.
 * @param[in]  ip         Index of the current LOS point.
 * @param[in,out] tau_path  Path transmittance array [nd][ng];
 *                          updated cumulatively for each gas.
 * @param[out] tau_seg    Total segment transmittance per channel [nd].
 *
 * @details
 * - Uses pretabulated emissivity data (`tbl->eps`) for each gas and channel.
 * - Applies logarithmic interpolation in pressure (@ref LOGX)
 *   and linear interpolation in temperature (@ref LIN).
 * - Enforces emissivity limits in the range [0, 1].
 * - Returns unity transmittance if data are missing or column density ≤ 0.
 *
 * @see ctl_t, tbl_t, los_t, intpol_tbl_eps, LIN, LOGX, locate_irr, locate_reg
 * 
 * @author Lars Hoffmann
 */
void intpol_tbl_cga(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/**
 * @brief Interpolate emissivities and transmittances using the
 *        Emissivity Growth Approximation (EGA).
 *
 * Computes gas transmittance along a line-of-sight segment by
 * interpolating emissivity values from lookup tables based on
 * the Emissivity Growth Approximation (EGA). The interpolation
 * is performed in pressure, temperature, and effective column
 * density space derived from the local LOS properties.
 *
 * @param[in]  ctl        Control structure defining number of gases
 *                        (@ref ctl_t::ng) and channels (@ref ctl_t::nd).
 * @param[in]  tbl        Emissivity lookup tables (@ref tbl_t)
 *                        containing tabulated pressure, temperature,
 *                        and column density grids.
 * @param[in]  los        Line-of-sight structure providing local
 *                        pressure, temperature, and column density data.
 * @param[in]  ip         Index of the current LOS point.
 * @param[in,out] tau_path  Path transmittance array [nd][ng];
 *                          updated cumulatively for each gas.
 * @param[out] tau_seg    Total segment transmittance per channel [nd].
 *
 * @details
 * - Uses pretabulated emissivity data (`tbl->eps`) and performs
 *   bilinear interpolation in pressure and temperature.
 * - Column density interpolation is handled by @ref intpol_tbl_u
 *   according to the emissivity growth relation.
 * - Enforces emissivity limits within [0, 1].
 * - Returns unity transmittance if lookup data are invalid or
 *   column density ≤ 0.
 *
 * @see ctl_t, tbl_t, los_t, intpol_tbl_u, intpol_tbl_eps, LIN, locate_irr, locate_reg
 * 
 * @note Implements the EGA variant of the forward model,
 *       selected when @ref ctl_t::formod = 1.
 *
 * @author Lars Hoffmann
 */
void intpol_tbl_ega(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const los_t * los,
  const int ip,
  double tau_path[ND][NG],
  double tau_seg[ND]);

/**
 * @brief Interpolate gas emissivity as a function of column amount.
 *
 * Computes emissivity \f$\varepsilon(u)\f$ from tabulated data using
 * linear interpolation in \f$\log u\f$–\f$\log\varepsilon\f$ space.
 * The input column amount must be provided as \f$\log(u)\f$.
 *
 * Behavior:
 * - For \f$u < u_{\min}\f$, emissivity scales linearly with column amount.
 * - For \f$u > u_{\max}\f$, emissivity asymptotically approaches unity using
 *   an exponential tail matched at \f$(u_{\max}, \varepsilon_{\max})\f$.
 *
 * The implementation uses `log1p`/`expm1` for numerical stability and
 * minimizes conversions between linear and logarithmic space.
 *
 * @param tbl   Lookup table structure.
 * @param ig    Gas index.
 * @param id    Spectral/channel index.
 * @param ip    Pressure index.
 * @param it    Temperature index.
 * @param logu  Natural logarithm of the column amount [molecules/cm²].
 *
 * @return Emissivity in the range \f$[0,1]\f$.
 *
 * @author Lars Hoffmann
 */
double intpol_tbl_eps(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double logu);

/**
 * @brief Interpolate column amount as a function of emissivity.
 *
 * Computes the column amount \f$u(\varepsilon)\f$ from tabulated data using
 * linear interpolation in \f$\log\varepsilon\f$–\f$\log u\f$ space.
 * The emissivity must be provided as \f$\log(\varepsilon)\f$ together with
 * \f$\log(1-\varepsilon)\f$.
 *
 * Behavior:
 * - For \f$\varepsilon < \varepsilon_{\min}\f$, column amount scales linearly
 *   with emissivity.
 * - For \f$\varepsilon > \varepsilon_{\max}\f$, the exponential tail used in
 *   `intpol_tbl_eps()` is analytically inverted.
 *
 * The implementation operates primarily in log-space and uses `log1p` for
 * numerical stability near \f$\varepsilon \rightarrow 1\f$.
 *
 * @param tbl     Lookup table structure.
 * @param ig      Gas index.
 * @param id      Spectral/channel index.
 * @param ip      Pressure index.
 * @param it      Temperature index.
 * @param logeps  Natural logarithm of emissivity (\f$\log\varepsilon\f$).
 *
 * @return Column amount \f$u\f$.
 *
 * @see intpol_tbl_eps
 *
 * @author Lars Hoffmann
 */
double intpol_tbl_u(
  const tbl_t * tbl,
  const int ig,
  const int id,
  const int ip,
  const int it,
  const double logeps);

/**
 * @brief Converts Julian seconds to calendar date and time components.
 *
 * This function converts Julian seconds to calendar date and time
 * components, including year, month, day, hour, minute, and
 * second. It also calculates the fractional part of the seconds.
 *
 * @param jsec Julian seconds to convert.
 * @param year Pointer to store the year.
 * @param mon Pointer to store the month.
 * @param day Pointer to store the day.
 * @param hour Pointer to store the hour.
 * @param min Pointer to store the minute.
 * @param sec Pointer to store the second.
 * @param remain Pointer to store the fractional part of seconds.
 *
 * The function initializes a time structure `t0` with a fixed
 * starting date and time. It then converts the Julian seconds to a
 * time_t type by adding the seconds to the epoch time. Next, it
 * converts the time_t value to a UTC time structure `t1`. Finally, it
 * extracts the year, month, day, hour, minute, and second components
 * from `t1` and calculates the fractional part of seconds, which is
 * stored in `remain`.
 *
 * @author Lars Hoffmann
 */
void jsec2time(
  const double jsec,
  int *year,
  int *mon,
  int *day,
  int *hour,
  int *min,
  int *sec,
  double *remain);

/*!
 * @brief Compute the Jacobian (kernel) matrix by finite differences.
 *
 * Evaluates the sensitivity of the simulated radiances to each
 * element of the atmospheric state vector by perturbing one
 * parameter at a time and re-running the forward model. The result
 * is the Jacobian matrix \f$ K = \partial y / \partial x \f$,
 * where *y* is the measurement vector and *x* is the state vector.
 *
 * @param[in]  ctl  Control structure defining retrieval configuration
 *                  and model setup.
 * @param[in]  tbl  Emissivity lookup tables used by the forward model.
 * @param[in]  atm  Atmospheric state vector and profile data.
 * @param[in]  obs  Observation geometry and radiance data.
 * @param[out] k    Jacobian matrix [m×n], where *m* is the number of
 *                  measurements and *n* the number of state variables.
 *
 * @details
 * - The undisturbed forward model is first computed to obtain the
 *   reference measurement vector.
 * - Each state vector element is perturbed by an adaptive step `h`
 *   depending on its physical type (pressure, temperature, VMR, etc.).
 * - For each perturbation, the forward model is re-evaluated, and
 *   the corresponding column of *K* is estimated using finite differences.
 * - Parallelized over state vector elements using OpenMP.
 *
 * @note
 * Typical perturbation sizes:
 * - Pressure: 1 % or ≥ 1e–7 hPa  
 * - Temperature: 1 K  
 * - VMR: 1 % or ≥ 1e–15  
 * - Extinction: 1e–4 km⁻¹  
 * - Cloud and surface parameters: 1 K or 1e–2 as appropriate.
 *
 * @see ctl_t, tbl_t, atm_t, obs_t, formod, x2atm, atm2x, obs2y, copy_atm, copy_obs
 * 
 * @warning Computationally intensive; requires one forward model
 *          evaluation per state vector element.
 *
 * @author Lars Hoffmann
 */
void kernel(
  const ctl_t * ctl,
  const tbl_t * tbl,
  atm_t * atm,
  obs_t * obs,
  gsl_matrix * k);

/**
 * @brief Locate index for interpolation on an irregular grid.
 *
 * Finds the lower index `ilo` such that \f$ xx[ilo] \le x < xx[ilo+1] \f$
 * for monotonically increasing or decreasing grids.  
 * Used in interpolation routines for altitude, pressure, temperature,
 * or wavenumber profiles that are not evenly spaced.
 *
 * @param[in] xx  Array of monotonic grid values (increasing or decreasing).
 * @param[in] n   Number of grid points.
 * @param[in] x   Target value to locate within the grid range.
 * @return Index `ilo` of the lower grid point surrounding `x`.
 *
 * @details
 * - Uses a binary search algorithm with \f$ O(\log n) \f$ complexity.  
 * - Handles both increasing and decreasing grids automatically.  
 * - Returns the index of the lower neighbor suitable for use in
 *   interpolation routines such as @ref LIN, @ref LOGX, or @ref LOGY.
 *
 * @see LIN, LOGX, LOGY, locate_reg, locate_tbl
 *
 * @warning Assumes `x` lies within the range of `xx`; no bounds checking
 *          beyond the first and last grid points is performed.
 *
 * @author Lars Hoffmann
 */
int locate_irr(
  const double *xx,
  const int n,
  const double x);

/**
 * @brief Locate index for interpolation on a regular (uniform) grid.
 *
 * Computes the lower index `i` such that \f$ xx[i] \le x < xx[i+1] \f$
 * for evenly spaced grid points. Used for quick index lookup when the
 * grid spacing is constant.
 *
 * @param[in] xx  Array of regularly spaced grid values.
 * @param[in] n   Number of grid points.
 * @param[in] x   Target value to locate within the grid range.
 * @return Index `i` of the lower grid point surrounding `x`.
 *
 * @details
 * - Computes the index directly from the grid spacing using
 *   \f$ i = (x - xx_0) / (xx_1 - xx_0) \f$.
 * - Clamps the result to `[0, n - 2]` to avoid out-of-bounds indices.
 * - Suitable for use with uniform grids such as pressure, temperature,
 *   or wavelength tables.
 *
 * @see locate_irr, locate_tbl, LIN, LOGX, LOGY
 * 
 * @warning Assumes uniform grid spacing; results are invalid for
 *          irregularly spaced arrays.
 *
 * @author Lars Hoffmann
 */
int locate_reg(
  const double *xx,
  const int n,
  const double x);

/**
 * @brief Locate index for interpolation within emissivity table grids.
 *
 * Finds the lower index `ilo` such that \f$ xx[ilo] \le x < xx[ilo+1] \f$
 * in a monotonically increasing single-precision grid.  
 * Used for emissivity and column density interpolation in table-based
 * routines such as @ref intpol_tbl_eps and @ref intpol_tbl_u.
 *
 * @param[in] xx  Monotonic (increasing) single-precision grid array.
 * @param[in] n   Number of grid points.
 * @param[in] x   Target value to locate within the grid range.
 * @return Index `ilo` of the lower grid point surrounding `x`.
 *
 * @details
 * - Implements a binary search with \f$ O(\log n) \f$ complexity.
 * - Optimized for lookup tables stored in `float` to minimize memory use.
 * - Returns an index suitable for linear interpolation with @ref LIN.
 *
 * @see intpol_tbl_eps, intpol_tbl_u, LIN, locate_irr, locate_reg
 *
 * @warning Assumes monotonic grid input (increasing order) and
 *          that `x` lies within the table range.
 *
 * @author Lars Hoffmann
 */
int locate_tbl(
  const float *xx,
  const int n,
  const double x);

/*!
 * @brief Invert a square matrix, optimized for diagonal or symmetric positive-definite matrices.
 *
 * Performs in-place inversion of the matrix \f$\mathbf{A}\f$ using either:
 * - **Fast diagonal inversion**, if the matrix is strictly diagonal.
 * - **Cholesky decomposition**, if the matrix is full and symmetric positive-definite.
 *
 * @param[in,out] a  Square matrix (`gsl_matrix`) to be inverted in place.
 *
 * @details
 * The function first checks whether the input matrix is diagonal by testing
 * all off-diagonal elements. If diagonal, each diagonal element \f$a_{ii}\f$
 * is replaced by its reciprocal \f$1/a_{ii}\f$.
 *
 * For non-diagonal matrices, a Cholesky decomposition is performed:
 * \f[
 * \mathbf{A} = \mathbf{L}\mathbf{L}^T
 * \f]
 * followed by inversion using the Cholesky factors, yielding
 * \f$\mathbf{A}^{-1}\f$.
 *
 * This approach assumes \f$\mathbf{A}\f$ is **symmetric and positive-definite**.
 *
 * @see gsl_linalg_cholesky_decomp, gsl_linalg_cholesky_invert
 *
 * @note
 * - The inversion is performed **in place**; the input matrix is overwritten.
 * - No explicit symmetry or definiteness checks are performed — invalid input
 *   may result in numerical instability or GSL errors.
 * - Diagonal detection assumes exact zeros for off-diagonal elements.
 *
 * @warning
 * For ill-conditioned matrices, consider using singular value decomposition (SVD)
 * or regularization methods instead of direct inversion.
 *
 * @author
 * Lars Hoffmann
 */
void matrix_invert(
  gsl_matrix * a);

/**
 * @brief Compute structured matrix products of the form \f$A^T B A\f$ or \f$A B A^T\f$.
 *
 * Evaluates matrix products commonly used in covariance propagation and
 * optimal estimation, depending on the specified transpose mode:
 * - **transpose = 1** → computes \f$\mathbf{A}^T \mathbf{B} \mathbf{A}\f$
 * - **transpose = 2** → computes \f$\mathbf{A} \mathbf{B} \mathbf{A}^T\f$
 *
 * The vector \f$\mathbf{b}\f$ represents the diagonal elements of
 * \f$\mathbf{B}\f$, i.e. a diagonal weighting or covariance matrix.
 *
 * @param[in]  a          Input matrix \f$\mathbf{A}\f$ (size m×n).
 * @param[in]  b          Vector representing the diagonal of \f$\mathbf{B}\f$ (length m or n).
 * @param[in]  transpose  Operation selector:
 *                        - 1 → compute \f$A^T B A\f$
 *                        - 2 → compute \f$A B A^T\f$
 * @param[out] c          Output matrix to store the resulting product.
 *
 * @details
 * - The function internally forms the scaled matrix
 *   \f$(B^{1/2} A)\f$ or \f$(A B^{1/2})\f$, then multiplies it using
 *   BLAS `dgemm` routines for efficiency:
 *   \f[
 *   A^T B A = (B^{1/2}A)^T (B^{1/2}A), \quad
 *   A B A^T = (A B^{1/2}) (A B^{1/2})^T
 *   \f]
 * - The input matrix \f$\mathbf{A}\f$ is not modified.
 * - This operation is typically used in computing gain matrices,
 *   propagated covariances, or sensitivity matrices in retrieval algorithms.
 *
 * @see gsl_blas_dgemm, matrix_invert
 *
 * @note
 * - Assumes \f$\mathbf{B}\f$ is diagonal (provided as a vector of its diagonal elements).
 * - The output matrix \f$\mathbf{C}\f$ must be pre-allocated to the correct size.
 * - No symmetry enforcement or normalization is applied.
 *
 * @warning
 * - If `transpose` is not 1 or 2, the function performs no operation.
 * - Numerical stability depends on the conditioning of A and the scaling of B.
 *
 * @author
 * Lars Hoffmann
 */
void matrix_product(
  const gsl_matrix * a,
  const gsl_vector * b,
  const int transpose,
  gsl_matrix * c);

/**
 * @brief Convert observation radiances into a measurement vector.
 *
 * Extracts all finite radiance values from the observation structure
 * and stores them sequentially in a GSL vector.  
 * Optionally records detector (`id`) and ray path (`ir`) indices for
 * each measurement element.
 *
 * @param[in]  ctl   Control structure containing observation setup (e.g. number of detectors).
 * @param[in]  obs   Observation data structure containing radiances.
 * @param[out] y     Measurement vector to store radiances (may be NULL).
 * @param[out] ida   Optional array to store detector indices (may be NULL).
 * @param[out] ira   Optional array to store ray path indices (may be NULL).
 * @return Number of valid (finite) radiance values added to the vector.
 *
 * @details
 * - Loops over all detector channels (`nd`) and ray paths (`nr`).
 * - Skips non-finite (`NaN` or `Inf`) radiance values.
 * - Produces a compact measurement vector for use in retrievals and Jacobian computations.
 *
 * @see atm2x, x2atm, kernel, formod
 * 
 * @warning Arrays `ida` and `ira` must be preallocated with sufficient size
 *          to hold all finite radiances if provided.
 *
 * @author Lars Hoffmann
 */
size_t obs2y(
  const ctl_t * ctl,
  const obs_t * obs,
  gsl_vector * y,
  int *ida,
  int *ira);

/**
 * @brief Perform optimal estimation retrieval using Levenberg–Marquardt minimization.
 *
 * This function performs an optimal estimation of atmospheric state variables based on
 * measured observations, a priori information, and forward radiative transfer modeling.
 * The estimation follows the Rodgers (2000) formalism and uses a Levenberg–Marquardt
 * algorithm to iteratively minimize the cost function:
 *
 *   χ² = (x - x_a)^T S_a⁻¹ (x - x_a) + (y - F(x))^T S_ε⁻¹ (y - F(x)),
 *
 * where `x` is the atmospheric state vector, `x_a` is its a priori estimate,
 * `S_a` is the a priori covariance matrix, `y` is the measurement vector,
 * `F(x)` is the forward model, and `S_ε` is the measurement error covariance.
 *
 * The routine updates the atmospheric state until convergence criteria are met
 * or the maximum number of iterations is reached. Optionally, the full retrieval
 * error budget and averaging kernel analysis are computed.
 *
 * @param[out] ret       Retrieval configuration and output container. Determines convergence,
 *                       kernel recomputation frequency, and error analysis options.
 * @param[in]  ctl       Control parameters describing problem setup (grids, species, etc.).
 * @param[in]  tbl       Lookup tables required by the forward model.
 * @param[in]  obs_meas  Measured observations used as input.
 * @param[out] obs_i     Intermediate and final modeled observations corresponding to the retrieved state.
 * @param[in]  atm_apr   A priori atmospheric state used as reference.
 * @param[out] atm_i     Atmospheric state vector to be iteratively retrieved and updated.
 * @param[out] chisq     Final value of the cost function (χ²) upon convergence.
 *
 * @note
 * - Aborts early if the problem dimension is zero (no observations or unknowns).
 * - State updates are constrained to physically meaningful bounds (pressure, temperature, etc.).
 * - Matrix computations are performed using GSL (GNU Scientific Library).
 * - If retrieval error analysis is enabled (`ret->err_ana`), the function produces:
 *   - Retrieval covariance matrix
 *   - Error decomposition (noise, forward model)
 *   - Gain matrix
 *   - Averaging kernel matrix and diagnostic analysis
 *
 * @warning
 * Input structures must be properly initialized. The function allocates several GSL matrices
 * and vectors, all of which are freed before returning. The caller is responsible only for
 * memory outside this function.
 *
 * @see formod(), cost_function(), analyze_avk(), set_cov_apr(), set_cov_meas()
 *
 * @par Reference
 *   Rodgers, C. D. (2000). *Inverse Methods for Atmospheric Sounding: Theory and Practice.*
 *
 * @author Lars Hoffmann
 */
void optimal_estimation(
  ret_t * ret,
  ctl_t * ctl,
  tbl_t * tbl,
  obs_t * obs_meas,
  obs_t * obs_i,
  atm_t * atm_apr,
  atm_t * atm_i,
  double *chisq);

/**
 * @brief Perform line-of-sight (LOS) ray tracing through the atmosphere.
 *
 * Computes the geometric path of a viewing ray from the observer to the
 * atmosphere (and possibly the surface), accounting for spherical geometry,
 * optional refraction, and cloud or surface interactions.  
 * Fills the LOS structure with pressure, temperature, gas concentrations,
 * extinction, and path length at each step.
 *
 * @param[in]  ctl  Control structure containing model and numerical settings.
 * @param[in]  atm  Atmospheric state structure (profiles of p, T, q, k, etc.).
 * @param[in,out] obs  Observation geometry and radiance data; updated tangent point.
 * @param[out] los  Line-of-sight structure to be populated with sampled quantities.
 * @param[in]  ir   Index of the current ray path in the observation set.
 *
 * @details
 * - Integrates along the viewing ray starting at the observer position.
 * - Performs stepwise propagation with step length `ds` determined by
 *   altitude and user-specified controls (`rayds`, `raydz`).
 * - Interpolates atmospheric variables at each step using @ref intpol_atm.
 * - Detects surface intersection or top-of-atmosphere exit and terminates accordingly.
 * - Optionally accounts for **refraction** via the refractive index `n(p, T)`.
 * - Accumulates **column densities** and **Curtis–Godson means** for each gas.
 * - Supports **cloud extinction** and **surface emissivity** interpolation.
 *
 * @note
 * The routine enforces that atmospheric grids include the surface (z = 0 km).
 * Rays starting above the atmosphere are propagated downward until entry.
 *
 * @see intpol_atm, tangent_point, formod_pencil, hydrostatic
 *
 * @warning
 * - Fails if the observer is below the surface or the atmosphere lacks z = 0.
 * - Aborts if the number of LOS points exceeds `NLOS`.
 * - Assumes monotonic altitude ordering in atmospheric data.
 *
 * @author Lars Hoffmann
 */
void raytrace(
  const ctl_t * ctl,
  const atm_t * atm,
  obs_t * obs,
  los_t * los,
  const int ir);

/**
 * @brief Read atmospheric input data from a file.
 *
 * This function reads atmospheric data from the file specified by
 * `filename`, optionally prefixed by the directory `dirname`. The
 * file format (ASCII or binary) is determined by the control structure
 * `ctl`. The data are stored in the atmospheric structure `atm`.
 *
 * Supported file formats are:
 * - ASCII format (`ctl->atmfmt == 1`)
 * - Binary format (`ctl->atmfmt == 2`)
 *
 * The function initializes the atmospheric data container, opens the
 * specified file, reads its contents according to the selected format, and
 * performs sanity checks on the number of data points. It also logs basic
 * statistical information about the loaded atmospheric fields (time, altitude,
 * longitude, latitude, pressure, temperature, and species/emitter mixing ratios).
 *
 * @param dirname
 *        Optional directory path where the atmospheric file resides.
 *        If NULL, only `filename` is used.
 *
 * @param filename
 *        Name of the atmospheric data file to read.
 *
 * @param ctl
 *        Pointer to a control structure specifying input parameters,
 *        file format, number of emitters, spectral windows, and additional
 *        atmospheric configuration.
 *
 * @param atm
 *        Pointer to an atmospheric data structure that will be filled with
 *        the values read from the file. The structure must be allocated
 *        before calling this function.
 *
 * @note The function aborts execution using ERRMSG on critical errors,
 *       such as failure to open the file, unknown file format, or absence
 *       of readable data.
 *
 * @warning Ensure that the `atm` structure has been properly allocated,
 *          and that the control parameters in `ctl` are valid before
 *          calling this function.
 *
 * @author Lars Hoffmann
 */
void read_atm(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Read atmospheric data in ASCII format.
 *
 * This function reads and parses atmospheric input data from an ASCII file
 * specified by its filename and stores the values in the atmospheric
 * structure \c atm. The number and type of fields to read are determined by
 * the control structure \c ctl. Each line of the ASCII file corresponds to a
 * single atmospheric data point.
 *
 * The expected order of fields in each line is:
 *   - Time
 *   - Altitude
 *   - Longitude
 *   - Latitude
 *   - Pressure
 *   - Temperature
 *   - Mixing ratios for each gas/emitter (\c ctl->ng values)
 *   - Extinction coefficients for each spectral window (\c ctl->nw values)
 *
 * Additionally, if cloud or surface layer parameters are enabled in \c ctl,
 * they are read once from the first line only:
 *   - Cloud layer: altitude, thickness, and extinction (\c ctl->ncl values)
 *   - Surface layer: temperature and emissivity (\c ctl->nsf values)
 *
 * @param filename
 *        Path to the ASCII file containing atmospheric data.
 *
 * @param ctl
 *        Pointer to a control structure defining the number of gases,
 *        spectral windows, and whether cloud or surface layer information
 *        should be read.
 *
 * @param atm
 *        Pointer to an initialized atmospheric structure where the parsed
 *        data points will be stored. The function updates \c atm->np to
 *        reflect the number of successfully read data records.
 *
 * @note The function reads until end-of-file is reached. Each successfully
 *       parsed line increments the atmospheric data point counter.
 *
 * @warning The function terminates execution using \c ERRMSG if more than
 *          \c NP data points are encountered, or if the input format deviates
 *          from expectations.
 *
 * @author Lars Hoffmann
 */
void read_atm_asc(
  const char *filename,
  const ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Read atmospheric data in binary format.
 *
 * This function reads atmospheric input data from a binary file specified
 * by its filename and stores the decoded values in the atmospheric structure
 * \c atm. The expected binary format is predefined and must match the
 * configuration provided in the control structure \c ctl. The function reads
 * a header containing metadata describing the dataset, followed by the
 * atmospheric fields and optional cloud and surface layer properties.
 *
 * The binary file layout is expected to follow this structure:
 *   1. **Magic identifier** (4 bytes, ignored except for presence)
 *   2. **Header integers**:
 *      - Number of gas species (\c ng)
 *      - Number of spectral windows (\c nw)
 *      - Number of cloud layer extinction elements (\c ncl)
 *      - Number of surface emissivity elements (\c nsf)
 *      These must match the corresponding values in \c ctl.
 *   3. **Data payload**:
 *      - Number of points (\c np)
 *      - Arrays of size \c np for time, altitude, longitude, latitude,
 *        pressure, and temperature
 *      - For each gas species: mixing ratio array of length \c np
 *      - For each spectral window: extinction coefficient array of length \c np
 *   4. **Optional layered parameters**:
 *      - Cloud layer altitude, thickness, and extinction (\c ctl->ncl values)
 *      - Surface temperature and emissivity (\c ctl->nsf values)
 *
 * @param filename
 *        Path to the binary file containing atmospheric data.
 *
 * @param ctl
 *        Pointer to a control structure specifying the expected dimensions of
 *        atmospheric fields and optional layers. Used for header validation.
 *
 * @param atm
 *        Pointer to an allocated atmospheric data structure that will be filled
 *        with the contents of the binary file. All arrays must be allocated
 *        prior to calling this function.
 *
 * @note This function does **not** allocate memory; it assumes storage for all
 *       atmospheric variables already exists and matches the expected sizes.
 *
 * @warning Execution is terminated via \c ERRMSG if:
 *          - The binary header does not match the control structure.
 *          - The binary file does not contain the expected amount of data.
 *
 * @author Lars Hoffmann
 */
void read_atm_bin(
  const char *filename,
  const ctl_t * ctl,
  atm_t * atm);

/**
 * @brief Read one atmospheric profile from a netCDF file.
 *
 * This routine reads a single profile (record) from a netCDF file into an
 * @ref atm_t structure. The file is expected to follow the JURASSIC netCDF
 * layout produced by @ref write_atm_nc(), using an unlimited @c profile
 * dimension and a fixed @c level dimension.
 *
 * The number of valid vertical levels for the requested profile is read from
 * @c nlev(profile) and stored in @c atm->np. All 2-D profile variables are
 * then read for the hyperslab @c (profile,0:atm->np-1).
 *
 * ## File layout (expected)
 * Dimensions:
 * - @c profile : unlimited (record dimension)
 * - @c level   : fixed (maximum number of vertical levels stored in the file)
 *
 * Variables:
 * - @c nlev(profile) : number of valid levels per profile (required)
 * - Core state (2-D): @c time,z,lon,lat,p,t(profile,level)
 * - Trace gases (2-D): one variable per species named @c ctl->emitter[ig]
 *   with dimensions @c (profile,level)
 * - Extinction windows (2-D): @c ext_win_%d(profile,level)
 * - Clouds (1-D, optional if @c ctl->ncl>0): @c cld_z, cld_dz, cld_k_%d(profile)
 * - Surface (1-D, optional if @c ctl->nsf>0): @c sf_t, sf_eps_%d(profile)
 *
 * @param filename Input netCDF file name.
 * @param ctl      Control structure (defines numbers/names of optional fields,
 *                 e.g., @c ng, @c nw, @c ncl, @c nsf, and @c emitter[]).
 * @param atm      Atmospheric profile structure to fill. On return, @c atm->np
 *                 contains the number of valid levels for the requested profile.
 * @param profile  Record index along the unlimited @c profile dimension.
 *
 * @note This function assumes that arrays in @p atm have static capacity @c NP
 *       and checks that @c atm->np is within @c [1,NP]. Errors are handled via
 *       the @c NC(...) macro.
 *
 * @author Lars Hoffmann
 */
void read_atm_nc(
  const char *filename,
  const ctl_t * ctl,
  atm_t * atm,
  int profile);

/**
 * @brief Read model control parameters from command-line and configuration input.
 *
 * Parses all numerical and string parameters required to initialize a
 * JURASSIC simulation, including atmospheric composition, radiative channels,
 * cloud and surface options, continua, ray-tracing setup, retrieval parameters,
 * and output settings.  
 * Populates the @ref ctl_t structure with all configuration values.
 *
 * @param[in]  argc  Argument count from the command line.
 * @param[in]  argv  Argument vector containing user-specified options.
 * @param[out] ctl   Control structure to be filled with parsed settings.
 *
 * @details
 * - Uses @ref scan_ctl to extract key–value pairs from the command line or control file.
 * - Initializes:
 *   - **Emitters and gases** (`NG`, `EMITTER`),
 *   - **Spectral channels** (`ND`, `NU`, `WINDOW`),
 *   - **Cloud parameters** (`NCL`, `CLNU`),
 *   - **Surface parameters** (`NSF`, `SFNU`, `SFTYPE`, `SFSZA`),
 *   - **Hydrostatic reference height** (`HYDZ`),
 *   - **Continuum flags** (`CTM_CO2`, `CTM_H2O`, `CTM_N2`, `CTM_O2`),
 *   - **Ray-tracing options** (`REFRAC`, `RAYDS`, `RAYDZ`),
 *   - **Field-of-view** (`FOV`),
 *   - **Retrieval limits** (`RETP_ZMIN`, `RETQ_ZMAX`, etc.),
 *   - **Output flags** (`WRITE_BBT`, `WRITE_MATRIX`),
 *   - **External forward model paths** (`RFMBIN`, `RFMHIT`, `RFMXSC`).
 * - Validates array bounds and logical parameter ranges.
 * - Automatically detects major gas species indices (e.g. CO₂, H₂O, N₂, O₂).
 * - Logs the executable name, version, and compilation time at startup.
 *
 * @see ctl_t, scan_ctl, find_emitter, read_shape, formod
 *
 * @warning
 * - Aborts if mandatory keys are missing or exceed defined limits (e.g. `NG > NG_MAX`).
 * - Requires consistent index ordering (e.g. `NCL > 1`, `NSF > 1`).
 * - Undefined or invalid parameters trigger `ERRMSG()` aborts.
 *
 * @author Lars Hoffmann
 */
void read_ctl(
  int argc,
  char *argv[],
  ctl_t * ctl);

/**
 * @brief Read a numerical matrix from an ASCII file.
 *
 * Loads values into a GSL matrix from a text file containing
 * sparse or indexed entries in tabular format.  
 * Each valid line is parsed for row and column indices and the corresponding value.
 *
 * @param[in]  dirname   Directory path containing the matrix file (may be NULL).
 * @param[in]  filename  Name of the matrix file to read.
 * @param[out] matrix    Pointer to the GSL matrix to be filled with values.
 *
 * @details
 * - Opens the specified file and scans it line by line.
 * - Initializes the matrix to zero before filling.
 * - Each line is expected to contain at least 13 formatted fields, where:
 *   - The first integer gives the row index,
 *   - The seventh integer gives the column index,
 *   - The thirteenth floating-point number is the matrix element value.  
 * - All successfully parsed entries are written to the corresponding positions
 *   in the GSL matrix using @c gsl_matrix_set().
 * - Non-matching lines are ignored.
 *
 * @see gsl_matrix, read_ctl, read_atm
 *
 * @warning
 * - Aborts if the file cannot be opened.
 * - Expects 0-based integer indices consistent with matrix dimensions.
 * - Lines not matching the expected 13-field format are skipped silently.
 *
 * @author Lars Hoffmann
 */
void read_matrix(
  const char *dirname,
  const char *filename,
  gsl_matrix * matrix);

/**
 * @brief Read observation data from an input file.
 *
 * This function reads atmospheric observation data from the specified file
 * and stores the results in the provided ::obs_t structure. The file may be
 * in ASCII or binary format, depending on the control settings passed via
 * ::ctl_t. After reading the data, the routine performs basic validation
 * (e.g., verifies that at least one observation entry was loaded) and logs
 * diagnostic statistics such as ranges of times, observer coordinates, view
 * point coordinates, tangent point coordinates, radiance or brightness
 * temperature values, and transmittances.
 *
 * The input file path is constructed from the provided directory and filename.
 * If a directory name is given, the file is assumed to reside within it;
 * otherwise, the filename is used as-is. Depending on the value of
 * `ctl->obsfmt`, the function dispatches either to ::read_obs_asc() for ASCII
 * files or ::read_obs_bin() for binary files.
 *
 * @param[in] dirname  Directory containing the input file, or `NULL` to use only @p filename.
 * @param[in] filename Name of the observation file to read.
 * @param[in] ctl      Pointer to a control structure specifying file format and other options.
 * @param[out] obs     Pointer to an observation structure where the data will be stored.
 *
 * @note The function terminates with an error message if the file cannot be
 *       opened, the observation format is unknown, or no valid data is read.
 *
 * @warning The @p obs structure must be properly allocated before calling this
 *          function. The function assumes its arrays are large enough to store
 *          all values contained in the input file.
 *
 * @see read_obs_asc(), read_obs_bin(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void read_obs(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  obs_t * obs);

/**
 * @brief Read ASCII-formatted observation data from a file.
 *
 * This function reads atmospheric observation data from an ASCII text file
 * specified by its filename and stores it in the provided ::obs_t structure.
 * Each line in the input file is expected to contain numerical values
 * representing a single observation record, including time, observer
 * coordinates, view point coordinates, tangent point coordinates, radiance
 * or brightness temperature values, and transmittances. The number of radiance
 * and transmittance values per record is determined by \c ctl->nd.
 *
 * The function reads the file line by line, tokenizes the data fields, and
 * fills the corresponding observation arrays. The number of successfully
 * read observation entries is stored in \c obs->nr.
 *
 * @param filename
 *        Path to the ASCII file containing the observation data.
 *
 * @param ctl
 *        Control structure containing metadata such as the number of spectral
 *        channels (\c nd).
 *
 * @param obs
 *        Observation structure where the parsed data will be stored.
 *
 * @note This is a C function and assumes that the \c obs structure has been
 *       preallocated with sufficient space for all records and spectral
 *       channels. No memory allocation is performed inside this routine.
 *
 * @warning The function terminates with an error message if the number of
 *          entries exceeds the predefined limit \c NR.
 *
 * @see read_obs(), read_obs_bin(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void read_obs_asc(
  const char *filename,
  const ctl_t * ctl,
  obs_t * obs);

/**
 * @brief Read binary-formatted observation data from a file.
 *
 * This C function reads observation data stored in a compact binary format
 * from a file specified by its filename and initializes the provided ::obs_t
 * structure with the values retrieved. The binary format begins with a header
 * that contains a magic identifier and the expected number of spectral
 * channels. The number of channels in the file must match \c ctl->nd, otherwise
 * the routine aborts with an error.
 *
 * After verifying the header, the function reads the number of ray paths and
 * then sequentially loads arrays corresponding to observation time, observer
 * location, view point location, tangent point location, radiance (or
 * brightness temperature), and transmittance data. The number of ray paths is
 * assigned to \c obs->nr. All arrays must have been allocated prior to calling
 * this function.
 *
 * @param filename
 *        Path to the binary file containing the observation data.
 *
 * @param ctl
 *        Pointer to a control structure specifying the number of spectral
 *        channels (\c nd) and other configuration settings.
 *
 * @param obs
 *        Pointer to an observation structure where the decoded binary data
 *        will be stored.
 *
 * @note This function does not perform any memory allocation. The caller must
 *       ensure that all arrays in \c obs have sufficient capacity for the data
 *       being read.
 *
 * @warning The function terminates with an error message if the binary header
 *          does not match the expected channel count, if more data than allowed
 *          by \c NR is encountered, or if any read operation fails.
 *
 * @see read_obs(), read_obs_asc(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void read_obs_bin(
  const char *filename,
  const ctl_t * ctl,
  obs_t * obs);

/**
 * @brief Read one observation profile from a NetCDF file.
 *
 * This function reads all geometric and spectral observation data for a
 * single profile from a NetCDF file using the variable-per-channel layout.
 *
 * The NetCDF file is expected to contain:
 *  - A dimension \c profile (unlimited)
 *  - A dimension \c ray (number of ray paths)
 *  - A variable \c nray(profile) giving the number of rays for each profile
 *  - Geometry variables with dimensions \c (profile, ray):
 *      - \c time, \c obs_z, \c obs_lon, \c obs_lat
 *      - \c vp_z,  \c vp_lon,  \c vp_lat
 *      - \c tp_z,  \c tp_lon,  \c tp_lat
 *  - Spectral variables with dimensions \c (profile, ray), one per channel:
 *      - \c rad_%.4f  (radiance or brightness temperature)
 *      - \c tau_%.4f  (transmittance)
 *    where the formatted frequency corresponds to \c ctl->nu[id].
 *
 * All data for the selected profile are read into the supplied \c obs_t
 * structure. The number of rays is read from \c nray(profile) and stored
 * in \c obs->nr. If the number of rays exceeds \c NR or is less than one,
 * an error is raised.
 *
 * @param filename  Path to the NetCDF observation file.
 * @param ctl       Pointer to the control structure defining spectral channels.
 * @param obs       Pointer to the observation structure to be filled.
 * @param profile   Zero-based index of the profile to read.
 *
 * @author Lars Hoffmann
 */
void read_obs_nc(
  const char *filename,
  const ctl_t * ctl,
  obs_t * obs,
  const int profile);

/**
 * @brief Read and spectrally convolve an RFM output spectrum.
 *
 * Opens the appropriate RFM ASCII spectrum file for a given tangent or
 * observation altitude and convolves the high-resolution spectrum with a
 * provided instrument filter function.  
 * Returns the integrated (filtered) radiance value.
 *
 * @param[in]  basename  Base filename of the RFM output (e.g. "rad" or "tra").
 * @param[in]  z         Altitude [km] used to select the corresponding RFM file.
 * @param[in]  nu        Wavenumber grid [cm⁻¹] of the filter function.
 * @param[in]  f         Filter transmission values corresponding to @p nu.
 * @param[in]  n         Number of points in @p nu and @p f arrays.
 *
 * @return Filtered radiance value integrated over the instrument bandpass.
 *
 * @details
 * - The routine looks for an RFM output file named
 *   `basename_<altitude_in_meters>.asc` (e.g. `rad_04500.asc`).
 * - If not found, it retries with altitude+1 meter to tolerate rounding.
 * - The file is read using @ref read_rfm_spec() into arrays of wavenumbers
 *   (`nurfm`) and radiances (`rad`).
 * - The input filter function @f$ f(\nu) @f$ is linearly interpolated onto the
 *   RFM wavenumber grid, and the spectrum is convolved as
 *   @f[
 *       R = \frac{\int f(\nu) \, I(\nu) \, d\nu}{\int f(\nu) \, d\nu}
 *   @f]
 * - Linear interpolation is used for both spectral alignment and filter sampling.
 * - Returns the resulting band-averaged radiance in the same units as the input.
 *
 * @see read_rfm_spec, locate_irr, LIN, formod_rfm
 *
 * @warning
 * - Aborts if the corresponding RFM file cannot be found.
 * - Assumes `nu` and `f` arrays are monotonic and have at least two points.
 * - Files must contain RFM ASCII spectra in expected column format.
 *
 * @author Lars Hoffmann
 */
double read_obs_rfm(
  const char *basename,
  const double z,
  const double *nu,
  const double *f,
  const int n);

/**
 * @brief Read retrieval configuration and error parameters.
 *
 * Initializes the retrieval control structure (`ret_t`) by reading all
 * iteration and uncertainty parameters from the command line or an input
 * control file using the `scan_ctl()` interface.
 *
 * @param[in]  argc  Number of command-line arguments.
 * @param[in]  argv  Command-line argument vector.
 * @param[in]  ctl   Pointer to global control structure (`ctl_t`) defining
 *                   the number of emitters, detectors, windows, clouds, etc.
 * @param[out] ret   Pointer to retrieval configuration structure (`ret_t`)
 *                   to be populated with iteration and error settings.
 *
 * @details
 * The function performs the following initialization steps:
 *
 * 1. **Iteration control parameters**
 *    - `KERNEL_RECOMP` — number of iterations between kernel recomputations.  
 *    - `CONV_ITMAX` — maximum number of retrieval iterations.  
 *    - `CONV_DMIN` — minimum normalized step size for convergence.
 *
 * 2. **Error analysis flag**
 *    - `ERR_ANA` — enables or disables retrieval error analysis (0 = off, 1 = on).
 *
 * 3. **Instrument and forward model errors**
 *    - `ERR_FORMOD[id]` — relative (%) forward model uncertainty per detector channel.  
 *    - `ERR_NOISE[id]` — absolute instrument noise per detector channel  
 *      [W/(m²·sr·cm⁻¹)] or [K] depending on `write_bbt`.
 *
 * 4. **Pressure and temperature retrieval uncertainties**
 *    - `ERR_PRESS`, `ERR_PRESS_CZ`, `ERR_PRESS_CH` — pressure error [%] and correlation lengths [km].  
 *    - `ERR_TEMP`, `ERR_TEMP_CZ`, `ERR_TEMP_CH` — temperature error [K] and correlation lengths [km].
 *
 * 5. **Volume mixing ratio (VMR) errors**
 *    - `ERR_Q[ig]`, `ERR_Q_CZ[ig]`, `ERR_Q_CH[ig]` — per gas [%] and correlation lengths [km].
 *
 * 6. **Extinction errors**
 *    - `ERR_K[iw]`, `ERR_K_CZ[iw]`, `ERR_K_CH[iw]` — per spectral window [km⁻¹] and correlation lengths [km].
 *
 * 7. **Cloud retrieval parameters**
 *    - `ERR_CLZ` — cloud top height error [km].  
 *    - `ERR_CLDZ` — cloud depth error [km].  
 *    - `ERR_CLK[icl]` — cloud extinction error per frequency [km⁻¹].
 *
 * 8. **Surface retrieval parameters**
 *    - `ERR_SFT` — surface temperature error [K].  
 *    - `ERR_SFEPS[isf]` — surface emissivity errors (dimensionless).
 *
 * @see scan_ctl, set_cov_apr, set_cov_meas, ret_t, ctl_t
 *
 * @note
 * - Each parameter can be specified either in the control file or on the
 *   command line (the latter overrides file values).
 * - Default values are used when a parameter is not explicitly defined.
 * - Correlation lengths of `-999` indicate uncorrelated (diagonal) treatment.
 *
 * @warning
 * - Input validation is minimal; ensure consistency between `ctl` and `ret` dimensions.
 * - Missing mandatory parameters trigger runtime errors.
 *
 * @author
 * Lars Hoffmann
 */
void read_ret(
  int argc,
  char *argv[],
  const ctl_t * ctl,
  ret_t * ret);

/**
 * @brief Read a Reference Forward Model (RFM) ASCII spectrum.
 *
 * Parses an RFM output file containing high-resolution spectral radiances
 * and fills the provided arrays with wavenumber and radiance values.
 *
 * @param[in]  filename  Name of the RFM ASCII spectrum file (e.g. "rad_04500.asc").
 * @param[out] nu        Array to receive the spectral wavenumber grid [cm⁻¹].
 * @param[out] rad       Array to receive the corresponding radiances.
 * @param[out] npts      Pointer to integer receiving the number of spectral points.
 *
 * @details
 * - Expects the RFM file to begin with a four-line header.
 *   The final header line must contain, in order:
 *   - the number of spectral points (`npts`),
 *   - starting wavenumber `nu0` [cm⁻¹],
 *   - spectral increment `dnu` [cm⁻¹],
 *   - ending wavenumber `nu1` [cm⁻¹].
 * - Radiance data follow as a sequence of floating-point values, separated
 *   by spaces, tabs, or line breaks.
 * - The wavenumber grid is reconstructed using linear interpolation between
 *   `nu0` and `nu1`:
 *   @f[
 *       \nu_i = \nu_0 + i \, \frac{\nu_1 - \nu_0}{N - 1}, \quad i = 0,\dots,N-1
 *   @f]
 * - Uses dynamic line buffering and token-based parsing for efficiency.
 *
 * @see read_obs_rfm, locate_irr, LIN, formod_rfm
 *
 * @warning
 * - Aborts if the file cannot be opened or has an invalid header format.
 * - Aborts if the number of grid points exceeds `RFMNPTS`.
 * - Assumes ASCII format consistent with standard RFM `.asc` output.
 *
 * @note
 * This routine reads only the spectral intensity data — additional
 * metadata (e.g., gas profiles, geometry) must be handled separately.
 *
 * @author Lars Hoffmann
 */
void read_rfm_spec(
  const char *filename,
  double *nu,
  double *rad,
  int *npts);

/**
 * @brief Read a two-column shape function from an ASCII file.
 *
 * Loads tabulated x–y data pairs (e.g., filter transmission or field-of-view
 * weighting function) into the provided arrays.
 *
 * @param[in]  filename  Name of the ASCII file containing the shape function.
 * @param[out] x         Array to receive the abscissa values.
 * @param[out] y         Array to receive the ordinate values.
 * @param[out] n         Pointer to integer receiving the number of data points read.
 *
 * @details
 * - The input file must contain at least two whitespace-separated columns:
 *   - Column 1: abscissa (`x`), typically wavenumber [cm⁻¹] or angular offset [deg].
 *   - Column 2: ordinate (`y`), typically transmission or weighting value.
 * - Comment lines or malformed entries are ignored.
 * - The routine logs the number of data points and their value ranges.
 * - Data are stored directly in the provided arrays for subsequent interpolation
 *   or convolution.
 *
 * @see gsl_stats_minmax, formod_fov, init_srcfunc, read_obs_rfm
 *
 * @warning
 * - Aborts if the file cannot be opened or if fewer than two valid points are read.
 * - Aborts if the number of data points exceeds `NSHAPE`.
 * - Assumes numeric ASCII format and monotonic ordering of `x` values.
 *
 * @note
 * This generic shape reader is used for both spectral filter functions and
 * angular field-of-view profiles.
 *
 * @author Lars Hoffmann
 */
void read_shape(
  const char *filename,
  double *x,
  double *y,
  int *n);

/**
 * @brief Read emissivity lookup tables from disk.
 *
 * Loads emissivity lookup tables for all detector/channel indices @c id and
 * all emitter/gas indices @c ig according to the configured table format
 * @c ctl->tblfmt:
 *  - ASCII (@c ctl->tblfmt == 1)
 *  - compact binary (@c ctl->tblfmt == 2)
 *  - netCDF (@c ctl->tblfmt == 3)
 *
 * The function allocates and initializes a @c tbl_t structure, reads the
 * corresponding lookup tables, and (depending on the table format) also loads
 * filter functions.
 *
 * Missing per-table files (ASCII/binary) or missing netCDF files/variables are
 * not fatal: the reader emits a warning and the corresponding lookup table is
 * treated as empty (i.e., @c tbl->np[id][ig] == 0).
 *
 * Diagnostic information about loaded tables (pressure levels, temperature
 * ranges, absorber amounts, emissivity ranges) may be written to the log.
 *
 * @param[in] ctl  Pointer to the control structure defining lookup table format,
 *                 filenames, emitters, and detector/channel frequencies.
 *
 * @return Pointer to an allocated and initialized @c tbl_t structure containing
 *         the loaded lookup tables.
 *
 * @warning If an unsupported lookup table format is specified via
 *          @c ctl->tblfmt, the function aborts with an error message.
 *
 * @note Memory for the returned structure is dynamically allocated and must be
 *       released by the caller.
 *
 * @author Lars Hoffmann
 */
tbl_t *read_tbl(
  const ctl_t * ctl);

/**
 * @brief Read a single ASCII emissivity lookup table.
 *
 * Reads one ASCII table for detector/channel index @p id and emitter/gas index
 * @p ig from a file named:
 *
 *     <ctl->tblbase>_<ctl->nu[id]>_<ctl->emitter[ig]>.tab
 *
 * The table file contains four columns:
 *  1. pressure [hPa]
 *  2. temperature [K]
 *  3. column density [molecules/cm^2]
 *  4. emissivity [-]
 *
 * The routine infers the pressure/temperature/column-density indices from
 * changes in the values read from the file. Values outside the supported ranges
 * for column density or emissivity are skipped and counted.
 *
 * If the file cannot be opened, a warning is issued and the table is treated
 * as empty (i.e., @c tbl->np[id][ig] == 0).
 *
 * @param[in]     ctl  Pointer to control structure specifying filenames and grids.
 * @param[in,out] tbl  Pointer to the table structure to be filled.
 * @param[in]     id   Detector/channel index.
 * @param[in]     ig   Emitter/gas index.
 *
 * @warning Aborts via @c ERRMSG() if table dimensions exceed @c TBLNP / @c TBLNT / @c TBLNU.
 *
 * @note Implementation detail: the ASCII reader may use an internal sentinel
 *       during parsing (historically @c np == -1 such that the first increment
 *       yields index 0). On return, @c tbl->np[id][ig] always represents the
 *       number of pressure levels (0 for an empty table).
 *
 * @author Lars Hoffmann
 */
void read_tbl_asc(
  const ctl_t * ctl,
  tbl_t * tbl,
  const int id,
  const int ig);

/**
 * @brief Read a single compact binary emissivity lookup table.
 *
 * Reads one binary table for detector/channel index @p id and emitter/gas index
 * @p ig from a file named:
 *
 *     <ctl->tblbase>_<ctl->nu[id]>_<ctl->emitter[ig]>.bin
 *
 * The file contains a length field followed by a packed byte blob created by
 * @c tbl_pack():
 *  - @c size_t nbytes   : number of bytes in the packed blob
 *  - @c uint8_t blob[]  : packed lookup table data
 *
 * The blob is unpacked into @p tbl using @c tbl_unpack().
 *
 * If the file cannot be opened, a warning is issued and the table is treated
 * as empty (i.e., @c tbl->np[id][ig] == 0).
 *
 * @param[in]     ctl  Pointer to control structure specifying filenames and grids.
 * @param[in,out] tbl  Pointer to the table structure to be filled.
 * @param[in]     id   Detector/channel index.
 * @param[in]     ig   Emitter/gas index.
 *
 * @warning Aborts via @c ERRMSG() if the file is malformed or unpacking fails.
 *
 * @see write_tbl_bin()
 * @see tbl_pack()
 * @see tbl_unpack()
 *
 * @author Lars Hoffmann
 */
void read_tbl_bin(
  const ctl_t * ctl,
  tbl_t * tbl,
  const int id,
  const int ig);

/**
 * @brief Read one packed emissivity lookup table from an open NetCDF file.
 *
 * Loads a previously stored emissivity lookup table for detector/channel
 * index @p id and emitter (trace gas) index @p ig from an already opened
 * NetCDF file. The file is expected to have been created by
 * @c write_tbl_nc().
 *
 * The NetCDF variable name is derived from the detector/channel frequency:
 *
 *     tbl_XXXX
 *
 * where @c XXXX is @c ctl->nu[id] formatted to four decimal places.
 *
 * Each variable is stored as a one-dimensional @c NC_UBYTE array containing
 * a packed lookup-table byte stream. The data are unpacked into @p tbl
 * using @c tbl_unpack().
 *
 * Missing variables are not fatal: a warning is issued and the corresponding
 * table is left empty (i.e., @c tbl->np[id][ig] remains zero).
 *
 * The NetCDF file handle @p ncid must refer to an open file and is not
 * opened or closed by this function.
 *
 * @param[in]     ctl   Control structure defining emitters, detectors, and
 *                      detector/channel frequencies.
 * @param[in,out] tbl   Table structure that receives the unpacked data.
 * @param[in]     id    Detector/channel index selecting @c ctl->nu[id].
 * @param[in]     ig    Emitter (trace gas) index selecting @c ctl->emitter[ig].
 * @param[in]     ncid  NetCDF file identifier of an already opened file.
 *
 * @note A temporary buffer is allocated to hold the packed lookup-table data
 *       prior to unpacking.
 *
 * @warning Aborts via @c ERRMSG() if the NetCDF variable exists but cannot be
 *          read or unpacked (e.g., corrupted contents).
 *
 * @see write_tbl_nc()
 * @see tbl_unpack()
 *
 * @author Lars Hoffmann
 */
void read_tbl_nc_channel(
  const ctl_t * ctl,
  tbl_t * tbl,
  int id,
  int ig,
  int ncid);

/**
 * @brief Scan control file or command-line arguments for a configuration variable.
 *
 * Searches for a named variable in the JURASSIC control file or command-line
 * arguments, returning its value as a double. Optionally stores the value as
 * a string and supports array-style parameters (e.g., `EMITTER[0]`, `NU[5]`).
 *
 * @param[in]  argc       Number of command-line arguments.
 * @param[in]  argv       Command-line argument vector.
 * @param[in]  varname    Name of the control variable to read.
 * @param[in]  arridx     Array index (use -1 for scalar variables).
 * @param[in]  defvalue   Default value if variable is not found (can be empty).
 * @param[out] value      Optional pointer to a string buffer receiving the value
 *                        (may be `NULL` if only numeric output is required).
 *
 * @return The variable value converted to `double`.
 *
 * @details
 * - The routine first attempts to open the control file provided as the first
 *   command-line argument (`argv[1]`), unless it starts with '-'.
 * - Variable names may appear as either:
 *   - `VAR` (scalar)
 *   - `VAR[index]` (explicit array index)
 *   - `VAR[*]` (wildcard entry applying to all indices)
 * - The search order is:
 *   1. Control file lines of the form:
 *      @code
 *      VAR[index] = VALUE
 *      VAR[*]     = VALUE
 *      @endcode
 *   2. Command-line arguments:
 *      @code
 *      ./jurassic ctlfile VAR[index] VALUE
 *      @endcode
 * - If no match is found:
 *   - The default value `defvalue` is used (if non-empty).
 *   - Otherwise, the routine aborts with an error.
 * - The variable value is printed to the log at verbosity level 1.
 *
 * @see read_ctl, LOG, ERRMSG
 *
 * @warning
 * - Aborts if the control file cannot be opened (unless skipped with '-').
 * - Aborts if a required variable is missing and no default is provided.
 * - Array bounds are not validated against internal limits; use `ctl_t` checks
 *   to ensure consistency.
 *
 * @note
 * - This utility simplifies control input parsing by supporting both command-line
 *   overrides and configuration files with the same syntax.
 * - String comparisons are case-insensitive.
 *
 * @author Lars Hoffmann
 */
double scan_ctl(
  int argc,
  char *argv[],
  const char *varname,
  const int arridx,
  const char *defvalue,
  char *value);

/**
 * @brief Construct the a priori covariance matrix \f$\mathbf{S_a}\f$ for retrieval parameters.
 *
 * Builds the full a priori covariance matrix based on specified retrieval
 * error assumptions and correlation lengths defined in the `ret_t` structure.
 * Each diagonal element represents the variance of an individual state
 * vector element, while off-diagonal terms encode spatial correlations
 * between parameters of the same type.
 *
 * @param[in]  ret   Retrieval configuration and error parameters (`ret_t`).
 * @param[in]  ctl   Control structure defining retrieval setup and state vector mapping (`ctl_t`).
 * @param[in]  atm   Atmospheric profile structure containing geolocation and altitude information (`atm_t`).
 * @param[in]  iqa   Index array linking state vector elements to physical quantities (e.g. pressure, temperature, gas, extinction, etc.).
 * @param[in]  ipa   Index array linking state vector elements to atmospheric grid points.
 * @param[out] s_a   Output a priori covariance matrix \f$\mathbf{S_a}\f$ (`gsl_matrix`), dimension n×n.
 *
 * @details
 * - The function first converts atmospheric quantities to a state vector (`atm2x`)
 *   and scales them according to their a priori uncertainties:
 *   - Pressure and trace gas errors are relative (% of nominal value).
 *   - Temperature, extinction, cloud, and surface errors are absolute.
 * - The diagonal of \f$\mathbf{S_a}\f$ is filled with the variances
 *   \f$\sigma_i^2\f$ of each element.
 * - Off-diagonal elements are populated according to an exponential
 *   correlation model:
 *   \f[
 *   \rho_{ij} = \exp\left(-\frac{d_{ij}}{L_h} - \frac{|z_i - z_j|}{L_v}\right)
 *   \f]
 *   where:
 *   - \f$d_{ij}\f$ is the great-circle (horizontal) distance between
 *     state vector locations,
 *   - \f$L_h\f$ and \f$L_v\f$ are horizontal and vertical correlation lengths
 *     for the parameter type.
 *
 * @see atm2x, geo2cart, DIST, ret_t, ctl_t
 *
 * @note
 * - Parameters with identical type indices (`iqa[i] == iqa[j]`)
 *   are assumed to share correlation properties.
 * - Correlation lengths are taken from `ret_t`, differing by parameter type:
 *   - Pressure: `err_press_cz`, `err_press_ch`
 *   - Temperature: `err_temp_cz`, `err_temp_ch`
 *   - Volume mixing ratio: `err_q_cz[]`, `err_q_ch[]`
 *   - Extinction: `err_k_cz[]`, `err_k_ch[]`
 * - Cloud and surface parameters are assumed uncorrelated (diagonal only).
 *
 * @warning
 * - A zero or negative variance triggers a runtime error.
 * - The matrix is constructed in full (dense), which may be large for
 *   high-resolution retrieval grids.
 *
 * @author Lars Hoffmann
 */
void set_cov_apr(
  const ret_t * ret,
  const ctl_t * ctl,
  const atm_t * atm,
  const int *iqa,
  const int *ipa,
  gsl_matrix * s_a);

/*!
 * @brief Construct measurement error standard deviations and their inverse.
 *
 * Builds the total measurement uncertainty vector used in the
 * optimal estimation retrieval, accounting for both instrument noise
 * and forward model (systematic) errors.
 *
 * @param[in]  ret          Retrieval configuration and error parameters (`ret_t`).
 * @param[in]  ctl          Control structure defining spectral channels and setup (`ctl_t`).
 * @param[in]  obs          Observation dataset (`obs_t`), containing measured radiances or brightness temperatures.
 * @param[out] sig_noise    Vector of instrument noise standard deviations (`gsl_vector`), length m.
 * @param[out] sig_formod   Vector of forward model error standard deviations (`gsl_vector`), length m.
 * @param[out] sig_eps_inv  Vector of inverse total standard deviations, \f$\sigma_\epsilon^{-1}\f$, used for normalization.
 *
 * @details
 * - The function computes the total measurement uncertainty for each
 *   observation element \f$i\f$ as:
 *   \f[
 *   \sigma_{\epsilon,i}^2 = \sigma_{\text{noise},i}^2 + \sigma_{\text{formod},i}^2
 *   \f]
 *   and stores its reciprocal square root:
 *   \f[
 *   (\sigma_{\epsilon,i}^{-1}) = \frac{1}{\sqrt{\sigma_{\epsilon,i}^2}}
 *   \f]
 *
 * - **Noise error (`sig_noise`)**  
 *   Determined from the instrument noise level defined in `ret->err_noise[id]`
 *   for each spectral channel. The noise term is always included in the fit.
 *
 * - **Forward model error (`sig_formod`)**  
 *   Computed as a fixed percentage (`ret->err_formod[id]`) of the
 *   measured radiance (or brightness temperature) per channel.
 *   This represents uncertainty due to imperfect forward modeling.
 *
 * - The inverse total standard deviation vector (`sig_eps_inv`)
 *   is used to normalize the measurement residuals \f$(y - F(x))\f$
 *   in the cost function.
 *
 * @see obs2y, copy_obs, cost_function, set_cov_apr
 *
 * @note
 * - Only finite observation elements are considered; invalid values are set to `NAN`.
 * - Units correspond to the observation quantity:
 *   - Radiance: [W/(m²·sr·cm⁻¹)]
 *   - Brightness temperature: [K]
 * - The forward model error is always relative, expressed in percent (%).
 *
 * @warning
 * - A zero or negative uncertainty triggers a runtime error.
 * - Assumes `obs` and `ctl` are consistent in dimension and indexing.
 *
 * @author Lars Hoffmann
 */
void set_cov_meas(
  const ret_t * ret,
  const ctl_t * ctl,
  const obs_t * obs,
  gsl_vector * sig_noise,
  gsl_vector * sig_formod,
  gsl_vector * sig_eps_inv);

/**
 * @brief Compute the solar zenith angle for a given time and location.
 *
 * Calculates the apparent solar zenith angle (SZA) [deg] based on
 * the observer’s longitude, latitude, and time since 2000-01-01 T00:00 Z.
 *
 * @param[in] sec  Seconds since 2000-01-01 T00:00 Z.
 * @param[in] lon  Observer longitude [deg].
 * @param[in] lat  Observer latitude [deg].
 *
 * @return Solar zenith angle [deg].
 *
 * @details
 * - Implements a simplified astronomical model based on the Sun’s
 *   apparent ecliptic longitude and the Earth’s mean obliquity.
 * - Uses the following steps:
 *   1. Compute the number of days since J2000 noon epoch.
 *   2. Determine the Sun’s apparent ecliptic longitude and declination.
 *   3. Compute the local hour angle from Greenwich Mean Sidereal Time.
 *   4. Derive the solar zenith angle via spherical trigonometry:
 *      @f[
 *      \cos(\theta) = \sin(\varphi)\sin(\delta)
 *                   + \cos(\varphi)\cos(\delta)\cos(h)
 *      @f]
 *   where @f$\varphi@f$ = latitude, @f$\delta@f$ = declination, @f$h@f$ = hour angle.
 *
 * @see formod_pencil, ctl_t::sfsza
 *
 * @note
 * - Neglects atmospheric refraction and seasonal perturbations.
 * - Accuracy is sufficient for radiative transfer applications (<0.1°).
 * - Longitude positive eastward, latitude positive northward.
 *
 * @author Lars Hoffmann
 */
double sza(
  double sec,
  double lon,
  double lat);

/**
 * @brief Determine the tangent point along a line of sight (LOS).
 *
 * Computes the location of the tangent point — the point of minimum altitude —
 * along the current line of sight, based on the LOS geometry stored in @ref los_t.
 *
 * @param[in]  los     Pointer to the line-of-sight (LOS) structure containing
 *                     altitude, longitude, latitude, and segment length data.
 * @param[out] tpz     Pointer to variable receiving tangent point altitude [km].
 * @param[out] tplon   Pointer to variable receiving tangent point longitude [deg].
 * @param[out] tplat   Pointer to variable receiving tangent point latitude [deg].
 *
 * @details
 * - For limb or occultation geometry, the routine:
 *   1. Identifies the LOS grid point with minimum altitude.
 *   2. Fits a quadratic interpolation polynomial:
 *      @f$ z = a x^2 + b x + c @f$
 *      through the altitudes of the three neighboring LOS points.
 *   3. Solves analytically for the vertex position @f$ x = -b / (2a) @f$,
 *      corresponding to the tangent point.
 *   4. Converts this interpolated position back to geographic coordinates.
 * - For nadir or zenith viewing (minimum altitude at the LOS endpoint),
 *   the tangent point defaults to the last grid point.
 *
 * @see raytrace, geo2cart, cart2geo, los_t
 *
 * @note
 * - The LOS segment lengths (`ds`) must be consistent with the geometric spacing
 *   between altitude points for the interpolation to be accurate.
 * - The quadratic interpolation provides sub-kilometer precision for smooth
 *   limb rays.
 * - Longitude and latitude are returned in degrees.
 *
 * @warning
 * - If the LOS contains fewer than three valid points, or the geometry is strongly
 *   curved, the tangent point estimate may be unreliable.
 *
 * @author Lars Hoffmann
 */
void tangent_point(
  const los_t * los,
  double *tpz,
  double *tplon,
  double *tplat);

/**
 * @brief Free lookup table and all internally allocated memory.
 *
 * Frees all dynamically allocated memory owned by a ::tbl_t object,
 * including the spectral lookup arrays (`logu` and `logeps`) for all
 * detector/emitter/pressure/temperature combinations. The ::tbl_t
 * structure itself is freed at the end.
 *
 * The function is safe to call with a NULL pointer and will return
 * immediately in that case. Partially initialized tables are handled
 * safely.
 *
 * @param[in,out] tbl  Pointer to the lookup table to be freed.
 * @param[in]     ctl  Control structure providing the number of
 *                     detectors and emitters used for looping.
 *
 * @note This function must be used instead of `free(tbl)` since ::tbl_t
 *       contains nested dynamically allocated members.
 *
 * @author Lars Hoffmann
 */
void tbl_free(
  const ctl_t * ctl,
  tbl_t * tbl);

/**
 * @brief Pack a lookup table into a contiguous binary buffer.
 *
 * This function serializes the lookup table data for a given detector
 * (`id`) and emitter (`ig`) into a compact, platform-native binary
 * representation. The packed data can be written to disk (e.g., in a
 * NetCDF variable) and later reconstructed using `tbl_unpack()`.
 *
 * The packed layout in `buf` is, in order:
 *  - `int np`                          : number of pressure grid points
 *  - `double p[np]`                   : pressure grid
 *  - for each pressure index `ip`:
 *      - `int nt`                     : number of temperature grid points
 *      - `double t[nt]`               : temperature grid
 *      - for each temperature index `it`:
 *          - `int nu`                 : number of spectral grid points
 *          - `float  u[nu]`           : spectral values
 *          - `float  eps[nu]`         : associated epsilon values
 *
 * No byte-order conversion or padding is applied; the data are written
 * exactly as laid out in memory.
 *
 * @param[in]  tbl         Table structure containing the lookup data.
 * @param[in]  id          Detector index.
 * @param[in]  ig          Emitter index.
 * @param[out] buf         Destination buffer that will receive the packed
 *                         binary data.
 * @param[out] bytes_used  Number of bytes written to `buf`.
 *
 * @warning The caller must ensure that `buf` is large enough to hold the
 *          packed data for the selected detector/emitter pair. The packed
 *          representation includes the detector filter function (mandatory for
 *          binary/netCDF formats). No bounds checking is performed inside this
 *          function.
 *
 * @see tbl_unpack()
 *
 * @author Lars Hoffmann
 */
void tbl_pack(
  const tbl_t * tbl,
  int id,
  int ig,
  uint8_t * buf,
  size_t *bytes_used);

/**
 * @brief Compute required buffer size (in bytes) for tbl_pack().
 *
 * Returns the exact number of bytes that tbl_pack() will write for the
 * given detector/emitter pair (including the mandatory filter function).
 *
 * @see tbl_pack()
 *
 * @author Lars Hoffmann
 */
size_t tbl_packed_size(
  const tbl_t * tbl,
  int id,
  int ig);

/**
 * @brief Unpack a lookup table from a contiguous binary buffer.
 *
 * This function reconstructs the lookup table data for a given detector
 * (`id`) and emitter (`ig`) from a binary buffer previously produced by
 * `tbl_pack()`. The packed data are read sequentially and copied into
 * the corresponding fields of the `tbl` structure.
 *
 * The expected layout of `buf` is:
 *  - `int np`
 *  - `double p[np]`
 *  - for each pressure index `ip`:
 *      - `int nt`
 *      - `double t[nt]`
 *      - for each temperature index `it`:
 *          - `int nu`
 *          - `float u[nu]`
 *          - `float eps[nu]`
 *
 * Range checks are applied to `np`, `nt`, and `nu` against the compile-time
 * limits `TBLNP`, `TBLNT`, and `TBLNU` to prevent buffer overruns and
 * invalid table sizes.
 *
 * @param[in,out] tbl  Table structure that will receive the unpacked data.
 * @param[in]     id   Detector index.
 * @param[in]     ig   Emitter index.
 * @param[in]     buf  Source buffer containing the packed binary table.
 *
 * @return The number of bytes consumed from `buf`.
 *
 * @warning The buffer must contain a valid packed table created by
 *          `tbl_pack()`. Supplying malformed or truncated data will
 *          result in undefined behavior or an error.
 *
 * @see tbl_pack()
 *
 * @author Lars Hoffmann
 */
size_t tbl_unpack(
  tbl_t * tbl,
  int id,
  int ig,
  const uint8_t * buf);


/**
 * @brief Converts time components to seconds since January 1, 2000, 12:00:00 UTC.
 *
 * This function calculates the number of seconds elapsed since
 * January 1, 2000, 12:00:00 UTC, based on the provided year, month,
 * day, hour, minute, and second. It also includes a fractional part
 * to represent the remaining seconds.
 *
 * @param year The year.
 * @param mon The month (1-12).
 * @param day The day of the month (1-31).
 * @param hour The hour of the day (0-23).
 * @param min The minute (0-59).
 * @param sec The second (0-59).
 * @param remain The fractional part of seconds.
 * @param jsec Pointer to store the calculated number of seconds since January 1, 2000, 12:00:00 UTC.
 *
 * The function calculates the time elapsed since January 1, 2000, 12:00:00 UTC, up to the specified time and includes
 * any fractional seconds indicated by the "remain" parameter.
 *
 * @note The function uses the timegm function, which is similar to mktime but operates in UTC.
 *
 * @author Lars Hoffmann
 */
void time2jsec(
  const int year,
  const int mon,
  const int day,
  const int hour,
  const int min,
  const int sec,
  const double remain,
  double *jsec);

/**
 * @brief Simple wall-clock timer for runtime diagnostics.
 *
 * Provides a lightweight timing utility based on `omp_get_wtime()`
 * to measure wall-clock durations between marked code regions.
 * The function supports up to ten concurrent nested timers.
 *
 * @param[in] name  Name or label of the timed code section.
 * @param[in] file  Source file name (usually `__FILE__` macro).
 * @param[in] func  Function name (usually `__func__` macro).
 * @param[in] line  Source line number (usually `__LINE__` macro).
 * @param[in] mode  Timer operation mode:
 *   - `1`: Start new timer.
 *   - `2`: Write elapsed time since last start (without stopping).
 *   - `3`: Write elapsed time and stop timer (pop one level).
 *
 * @details
 * - Each call with `mode == 1` starts a new timer instance and stores
 *   its start time and corresponding source line.
 * - When `mode == 2` or `mode == 3` is called, the elapsed wall-clock
 *   time (in seconds) is computed using:
 *   @f[
 *     \Delta t = t_{\text{now}} - t_{\text{start}}
 *   @f]
 *   and written to the log via the @ref LOG macro.
 * - Supports nested timers (up to 10 levels). Exceeding this limit
 *   triggers a runtime error via @ref ERRMSG.
 *
 * @see LOG, ERRMSG, omp_get_wtime
 *
 * @note
 * - The timing precision and resolution depend on the OpenMP runtime.
 * - Intended for coarse profiling and diagnostic output; not thread-safe.
 * - Lines reported in log messages indicate the start–stop interval.
 *
 * @warning
 * - Exceeding 10 nested timers results in an error.
 * - Calling `mode == 2` or `3` without a prior start causes an internal error.
 *
 * @author Lars Hoffmann
 */
void timer(
  const char *name,
  const char *file,
  const char *func,
  int line,
  int mode);

/**
 * @brief Write atmospheric data to a file.
 *
 * This function writes the atmospheric dataset stored in `atm` to the file
 * specified by `filename`, optionally prefixed by `dirname`. The output format
 * (ASCII or binary) is selected based on the atmospheric format flag
 * `ctl->atmfmt`. The function creates the output file, delegates the writing
 * process to the appropriate format-specific routine, and logs summary
 * statistics of the written atmospheric data.
 *
 * Supported output formats:
 *   - ASCII format (`ctl->atmfmt == 1`), written by `write_atm_asc()`
 *   - Binary format (`ctl->atmfmt == 2`), written by `write_atm_bin()`
 *
 * The function writes:
 *   - Atmospheric profiles: time, altitude, longitude, latitude, pressure,
 *     temperature
 *   - Gas mixing ratios for each emitter (`ctl->ng`)
 *   - Extinction coefficients for each spectral window (`ctl->nw`)
 *   - Optional cloud layer or surface parameters if enabled in `ctl`
 *
 * After writing, the function logs minimum and maximum values of the written
 * fields for verification and diagnostic purposes.
 *
 * @param dirname
 *        Optional directory in which the output file will be created.
 *        If NULL, only `filename` is used.
 *
 * @param filename
 *        Name of the output file to be created and populated with atmospheric data.
 *
 * @param ctl
 *        Pointer to a control structure defining the output format and the sizes
 *        of gas, spectral, cloud, and surface parameter arrays.
 *
 * @param atm
 *        Pointer to the atmospheric data structure whose contents will be written.
 *        All required fields must be initialized and contain `atm->np` valid data points.
 *
 * @note The function aborts execution using `ERRMSG` if the file cannot be created
 *       or if an unsupported output format is requested.
 *
 * @author Lars Hoffmann
 */
void write_atm(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/**
 * @brief Write atmospheric data to an ASCII file.
 *
 * This function writes the contents of an atmospheric structure \c atm to an
 * ASCII file specified by its filename. A descriptive column header is written
 * first, documenting the meaning, units, and ordering of each data field.
 * Atmospheric data points are then written line by line, with optional cloud
 * and surface layer parameters appended if they are enabled in the control
 * structure \c ctl.
 *
 * The output columns include, in order:
 *   1. Time (seconds since 2000-01-01T00:00Z)
 *   2. Altitude [km]
 *   3. Longitude [deg]
 *   4. Latitude [deg]
 *   5. Pressure [hPa]
 *   6. Temperature [K]
 *   + Gas/emitter mixing ratios for each species (\c ctl->ng) [ppv]
 *   + Extinction values for each spectral window (\c ctl->nw) [km^-1]
 *
 * If cloud layer properties are enabled (\c ctl->ncl > 0), the following are added:
 *   - Cloud layer height [km]
 *   - Cloud layer depth [km]
 *   - Cloud extinction values for each frequency (\c ctl->ncl) [km^-1]
 *
 * If surface layer properties are enabled (\c ctl->nsf > 0), the following are added:
 *   - Surface layer height [km]
 *   - Surface layer pressure [hPa]
 *   - Surface layer temperature [K]
 *   - Surface emissivity values (\c ctl->nsf)
 *
 * @param filename
 *        Path to the ASCII output file.
 *
 * @param ctl
 *        Pointer to a control structure defining the number of gases,
 *        spectral windows, and whether cloud or surface layer information
 *        should be included.
 *
 * @param atm
 *        Pointer to the atmospheric structure containing the data to be written.
 *        The function writes all \c atm->np data points.
 *
 * @note A blank line is inserted each time the time coordinate changes, grouping
 *       data points belonging to different timestamps.
 *
 * @warning The function assumes that all arrays in \c atm are properly allocated
 *          and populated. No validation of data ranges is performed here.
 *
 * @author Lars Hoffmann
 */
void write_atm_asc(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/**
 * @brief Write atmospheric data to a binary file.
 *
 * This function writes the atmospheric dataset contained in \c atm to a binary
 * file specified by its filename. The output format is compact and includes a
 * file header followed by the serialized atmospheric fields. The format is
 * compatible with \c read_atm_bin(), ensuring that files written by this
 * function can be read back without loss of information.
 *
 * The binary file structure written is as follows:
 *   1. **Magic identifier** \c "ATM1" (4 bytes)
 *   2. **Header integers** describing dataset layout:
 *        - Number of gas/emitter species (\c ctl->ng)
 *        - Number of spectral windows (\c ctl->nw)
 *        - Number of cloud extinction values (\c ctl->ncl)
 *        - Number of surface emissivity values (\c ctl->nsf)
 *   3. **Data payload**:
 *        - Number of atmospheric points \c np
 *        - Arrays of length \c np containing:
 *            * Time
 *            * Altitude
 *            * Longitude
 *            * Latitude
 *            * Pressure
 *            * Temperature
 *        - Gas mixing ratios for all emitters (\c ctl->ng × \c np)
 *        - Extinction coefficients for all spectral windows (\c ctl->nw × \c np)
 *   4. **Optional parameters** written only if enabled in \c ctl:
 *        - Cloud layer height, depth, and extinction values (\c ctl->ncl)
 *        - Surface temperature and emissivity values (\c ctl->nsf)
 *
 * @param filename
 *        Path to the binary output file.
 *
 * @param ctl
 *        Pointer to a control structure specifying the number of gases,
 *        spectral windows, and whether cloud or surface layer parameters
 *        must be included.
 *
 * @param atm
 *        Pointer to the atmospheric data structure containing values to be
 *        written. All arrays must be populated and \c atm->np must contain the
 *        number of valid atmospheric records.
 *
 * @note This function performs no range checking or validation of the \c atm
 *       contents. It assumes that the memory layout matches expectations.
 *
 * @warning The binary structure must remain consistent with
 *          \c read_atm_bin(); modifying either implementation requires
 *          updating the other accordingly.
 *
 * @author Lars Hoffmann
 */
void write_atm_bin(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/**
 * @brief Write one atmospheric profile to a netCDF file.
 *
 * This routine writes a single profile from an @ref atm_t structure into a
 * netCDF file using an unlimited record dimension. Multiple profiles can be
 * stored in the same file by calling this function repeatedly with increasing
 * @p profile indices.
 *
 * The function creates the file if it does not exist (netCDF-4/HDF5) and
 * defines missing dimensions/variables on the fly. Existing definitions are
 * reused.
 *
 * ## File layout
 * Dimensions:
 * - @c profile : unlimited (record dimension)
 * - @c level   : fixed (maximum number of vertical levels stored in the file)
 *
 * Variables:
 * - @c nlev(profile) : number of valid levels for each profile (written as
 *   @c atm->np)
 * - Core state (2-D): @c time,z,lon,lat,p,t(profile,level)
 * - Trace gases (2-D): one variable per species named @c ctl->emitter[ig]
 *   with dimensions @c (profile,level)
 * - Extinction windows (2-D): @c ext_win_%d(profile,level)
 * - Clouds (1-D, optional if @c ctl->ncl>0): @c cld_z, cld_dz, cld_k_%d(profile)
 * - Surface (1-D, optional if @c ctl->nsf>0): @c sf_t, sf_eps_%d(profile)
 *
 * Only the first @c atm->np elements along @c level are written for 2-D
 * variables. The effective number of levels is stored in @c nlev(profile).
 * (If a profile is overwritten with fewer levels than previously written, old
 * values above @c nlev may remain unless explicitly filled by the caller.)
 *
 * @param filename Output netCDF file name.
 * @param ctl      Control structure (defines numbers/names of optional fields,
 *                 e.g., @c ng, @c nw, @c ncl, @c nsf, and @c emitter[]).
 * @param atm      Atmospheric profile data to write (uses @c atm->np levels).
 * @param profile  Record index along the unlimited @c profile dimension.
 *
 * @note This function requires netCDF support (@c HAVE_NETCDF) and uses the
 *       netCDF C API. Errors are handled via the @c NC(...) macro.
 *
 * @author Lars Hoffmann
 */
void write_atm_nc(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm,
  int profile);

/**
 * @brief Write atmospheric profile in RFM-compatible format.
 *
 * Exports the current atmospheric state to a file formatted for use
 * with the Reference Forward Model (RFM). The file includes altitude,
 * pressure, temperature, and volume mixing ratio profiles for each
 * active emitter.
 *
 * @param[in] filename  Output file name for the RFM atmosphere file.
 * @param[in] ctl        Pointer to the control structure defining active emitters.
 * @param[in] atm        Pointer to the atmospheric profile to export.
 *
 * @details
 * - Produces a plain-text RFM atmosphere file with the following sections:
 *   @code
 *   NLAYERS
 *   *HGT [km]
 *   <altitude_1>
 *   ...
 *   *PRE [mb]
 *   <pressure_1>
 *   ...
 *   *TEM [K]
 *   <temperature_1>
 *   ...
 *   *<EMITTER> [ppmv]
 *   <mixing_ratio_1>
 *   ...
 *   *END
 *   @endcode
 * - The first line specifies the number of vertical layers (`atm->np`).
 * - Each subsequent block begins with a keyword (e.g., `*HGT`, `*PRE`, etc.)
 *   and lists one value per line.
 * - Mixing ratios are converted from parts per volume (ppv) to parts per million (ppmv).
 *
 * @see read_atm, ctl_t, atm_t
 *
 * @note
 * - Compatible with the RFM “ATM” input file format.
 * - Units:
 *   - Altitude in kilometers [km]
 *   - Pressure in millibars [mb]
 *   - Temperature in Kelvin [K]
 *   - Mixing ratios in parts per million by volume [ppmv]
 *
 * @warning
 * - Existing files with the same name will be overwritten.
 * - The function assumes consistent vertical ordering (surface → top of atmosphere).
 *
 * @author Lars Hoffmann
 */
void write_atm_rfm(
  const char *filename,
  const ctl_t * ctl,
  const atm_t * atm);

/**
 * @brief Write a fully annotated matrix (e.g., Jacobian or gain matrix) to file.
 *
 * Outputs a numerical matrix along with detailed metadata describing
 * the row and column spaces. Depending on configuration, the rows
 * and columns may correspond to measurement or state variables.
 *
 * @param[in] dirname   Output directory path (may be `NULL`).
 * @param[in] filename  Output file name.
 * @param[in] ctl       Pointer to control structure defining model setup and metadata.
 * @param[in] matrix    Pointer to GSL matrix to write (e.g., Jacobian, kernel, or covariance).
 * @param[in] atm       Pointer to atmospheric data structure (used when state-space indexing applies).
 * @param[in] obs       Pointer to observation data structure (used when measurement-space indexing applies).
 * @param[in] rowspace  Selects row labeling: `"y"` = measurement space, otherwise state space.
 * @param[in] colspace  Selects column labeling: `"y"` = measurement space, otherwise state space.
 * @param[in] sort      Determines writing order: `"r"` = row-major, otherwise column-major.
 *
 * @details
 * - This routine writes one matrix element per line, including descriptive metadata:
 *   @code
 *   RowIndex RowMeta... ColIndex ColMeta... MatrixValue
 *   @endcode
 * - The row and column metadata differ depending on space selection:
 *   - **Measurement space (`'y'`)**:
 *     - Channel wavenumber [cm⁻¹]
 *     - Observation time [s since 2000-01-01T00:00Z]
 *     - View point altitude [km], longitude [°], latitude [°]
 *   - **State space**:
 *     - Quantity name (e.g., TEMPERATURE, H2O)
 *     - Time, altitude, longitude, latitude of the profile point
 * - The header clearly documents all output columns for traceability.
 *
 * @see read_matrix, kernel, atm2x, obs2y, idx2name
 *
 * @note
 * - The function respects `ctl->write_matrix` — output is skipped if disabled.
 * - Output is human-readable and can be post-processed using external tools (e.g., Python, MATLAB, GNU Octave).
 * - Typically used for writing Jacobians, gain matrices, or averaging kernels.
 * - Matrix orientation can be changed with `sort` to support row-major or column-major output.
 *
 * @warning
 * - Large matrices may produce very large output files.
 * - Memory allocation is performed for temporary indexing arrays; ensure sufficient resources for large N, M.
 * - The function overwrites existing files without confirmation.
 *
 * @author Lars Hoffmann
 */
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

/**
 * @brief Write observation data to an output file in ASCII or binary format.
 *
 * This C function constructs the full output file path from the provided
 * directory and filename, opens the file for writing, and exports the
 * contents of the ::obs_t structure in either ASCII or binary format,
 * depending on the observation format specified by `ctl->obsfmt`. The actual
 * writing of formatted data is delegated to ::write_obs_asc() or
 * ::write_obs_bin().
 *
 * After writing, the function prints diagnostic information showing ranges of
 * times, observer coordinates, view point coordinates, tangent point
 * coordinates, radiance or brightness temperature values (depending on
 * `ctl->write_bbt`), and transmittances. These diagnostics provide useful
 * verification that the output data is valid and consistent.
 *
 * @param[in] dirname  Optional directory path. If NULL, only @p filename is used.
 * @param[in] filename Name of the output observation file.
 * @param[in] ctl      Control structure specifying output format, spectral
 *                     channel configuration, and brightness-temperature mode.
 * @param[in] obs      Observation structure containing the data to be written.
 *
 * @note This is a C function. The output file is always overwritten if it
 *       already exists.
 *
 * @warning The routine aborts with an error message if the output file cannot
 *          be created, or if `ctl->obsfmt` specifies an unsupported format.
 *
 * @see write_obs_asc(), write_obs_bin(), read_obs(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void write_obs(
  const char *dirname,
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs);

/**
 * @brief Write observation data to an ASCII text file.
 *
 * This C function writes the contents of the ::obs_t observation structure as
 * human-readable ASCII text to a file specified by its filename. It first
 * prints a descriptive header that documents each column of the output format,
 * including observation time, observer and view geometry, tangent point
 * information, and spectral values. The number and meaning of spectral fields
 * depend on \c ctl->nd and whether brightness temperature output is enabled
 * via \c ctl->write_bbt.
 *
 * The function then writes one line of data per ray path, including the base
 * geometric information followed by radiance or brightness temperature values
 * and transmittances for each spectral channel. Blank lines are inserted
 * whenever the time stamp changes, providing visual separation of distinct
 * observation groups.
 *
 * @param filename
 *        Path to the ASCII output file.
 *
 * @param ctl
 *        Control structure specifying the number of spectral
 *        channels (\c nd), wavenumbers (\c nu), and output mode
 *        (\c write_bbt).
 *
 * @param obs
 *        Observation structure containing the data to be written.
 *
 * @note This routine produces plain-text output intended for inspection,
 *       debugging, and compatibility with external processing tools.
 *
 * @warning The caller must ensure that the file specified by \c filename is
 *          writable. Existing files may be overwritten.
 *
 * @see write_obs(), write_obs_bin(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void write_obs_asc(
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs);

/**
 * @brief Write observation data in binary format to a file.
 *
 * This C function serializes the contents of the ::obs_t structure into a
 * compact binary format and writes it to a file specified by its filename.
 * The binary format begins with a header consisting of a magic identifier
 * (\c "OBS1") and the number of spectral channels (\c ctl->nd). This header is
 * used by ::read_obs_bin() to validate compatibility when reading.
 *
 * Following the header, the function writes the number of ray paths and then
 * sequentially outputs arrays of observation metadata, geometric parameters,
 * radiance or brightness temperature values, and transmittances. All values
 * are written in native binary representation using the FWRITE() macro, which
 * performs buffered writes and error checking.
 *
 * @param filename
 *        Path to the binary output file.
 *
 * @param ctl
 *        Control structure specifying the number of spectral
 *        channels (\c nd) and corresponding configuration parameters.
 *
 * @param obs
 *        Observation structure containing the data to be written.
 *
 * @note This routine does not perform any formatting or conversion. The
 *       resulting file is portable only to systems with compatible binary
 *       layouts (integer size, floating-point format, and endianness).
 *
 * @warning The caller must ensure that the file specified by \c filename is
 *          writable. Existing files may be overwritten.
 *
 * @see write_obs(), write_obs_asc(), read_obs_bin(), ctl_t, obs_t
 *
 * @author Lars Hoffmann
 */
void write_obs_bin(
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs);

/**
 * @brief Write one observation profile to a NetCDF file.
 *
 * This function writes all geometric and spectral observation data for a
 * single profile into a NetCDF file using a variable-per-channel layout.
 * If the file does not exist it is created; if it already exists, the
 * required dimensions and variables are created if missing and then
 * extended by writing the specified profile.
 *
 * The NetCDF file will contain:
 *  - A dimension \c profile (unlimited)
 *  - A dimension \c ray (fixed, number of ray paths)
 *  - A variable \c nray(profile) giving the number of rays for each profile
 *  - Geometry variables with dimensions \c (profile, ray):
 *      - \c time, \c obs_z, \c obs_lon, \c obs_lat
 *      - \c vp_z,  \c vp_lon,  \c vp_lat
 *      - \c tp_z,  \c tp_lon,  \c tp_lat
 *  - Spectral variables with dimensions \c (profile, ray), one per channel:
 *      - \c rad_%.4f  (radiance or brightness temperature)
 *      - \c tau_%.4f  (transmittance)
 *    where the formatted frequency corresponds to \c ctl->nu[id].
 *
 * Radiance is written either as physical radiance or as brightness
 * temperature depending on \c ctl->write_bbt. The ray dimension is fixed
 * on first creation and all subsequently written profiles must not exceed
 * this number of rays.
 *
 * @param filename  Path to the NetCDF observation file.
 * @param ctl       Pointer to the control structure defining spectral channels.
 * @param obs       Pointer to the observation structure containing the data
 *                  to be written.
 * @param profile   Zero-based index of the profile to write.
 *
 * @author Lars Hoffmann
 */
void write_obs_nc(
  const char *filename,
  const ctl_t * ctl,
  const obs_t * obs,
  const int profile);

/**
 * @brief Write tabulated shape function data to a text file.
 *
 * Exports a shape function (typically a weighting or field-of-view profile)
 * defined by paired arrays of *x* and *y* values to an ASCII file.
 *
 * @param[in] filename  Output file name.
 * @param[in] x         Pointer to array of x-values (independent variable).
 * @param[in] y         Pointer to array of y-values (dependent variable).
 * @param[in] n         Number of data points to write.
 *
 * @details
 * - Writes a plain-text table with two columns:
 *   @code
 *   # $1 = shape function x-value [-]
 *   # $2 = shape function y-value [-]
 *   @endcode
 * - Each line contains one (*x*, *y*) pair written with high precision.
 * - Typically used to export field-of-view functions, apodization kernels,
 *   or any other normalized shape profiles used by the model.
 *
 * @see read_shape
 *
 * @note
 * - Units are dimensionless unless otherwise defined by the application.
 * - The file is fully compatible with `read_shape()` for re-import.
 * - Output precision is set to 10 significant digits for numerical stability.
 *
 * @warning
 * - Existing files with the same name will be overwritten.
 * - The number of points *n* must be consistent with the size of *x* and *y* arrays.
 *
 * @author Lars Hoffmann
 */
void write_shape(
  const char *filename,
  const double *x,
  const double *y,
  const int n);

/**
 * @brief Write retrieval standard deviation profiles to disk.
 *
 * Extracts the diagonal elements of a covariance matrix (a priori,
 * posterior, or error covariance) to obtain the standard deviations
 * of retrieved quantities and writes them as an atmospheric profile file.
 *
 * @param[in]  quantity  Name of the retrieved quantity (e.g., `"apr"`, `"pos"`, `"err"`),
 *                       used to label the output file.
 * @param[in]  ret       Retrieval configuration structure (`ret_t`),
 *                       providing the working directory for output files.
 * @param[in]  ctl       Global control structure (`ctl_t`) defining retrieval setup and quantities.
 * @param[in]  atm       Reference atmospheric state (`atm_t`) for spatial/geometric metadata.
 * @param[in]  s         Covariance matrix (`gsl_matrix`, n×n) from which standard deviations
 *                       are derived (typically posterior covariance \f$\mathbf{S}\f$).
 *
 * @details
 * This function performs the following operations:
 * 1. Extracts the standard deviation vector \f$\sigma_i = \sqrt{S_{ii}}\f$
 *    from the diagonal of the covariance matrix \f$\mathbf{S}\f$.
 * 2. Copies the reference atmospheric structure (`atm`) into an auxiliary
 *    structure (`atm_aux`) to preserve coordinate and geometric metadata.
 * 3. Converts the standard deviation vector into the atmospheric representation
 *    using `x2atm()`, thereby mapping elements of the state vector to the
 *    corresponding atmospheric quantities.
 * 4. Writes the result to disk as a diagnostic file:
 *    \f[
 *    \texttt{<ret->dir>/atm\_err\_<quantity>.tab}
 *    \f]
 *    using the standard JURASSIC atmospheric file format.
 *
 * @see x2atm, copy_atm, write_atm, set_cov_apr, set_cov_meas
 *
 * @note
 * - The file naming convention follows `atm_err_<quantity>.tab`
 *   (e.g., `atm_err_apr.tab`, `atm_err_pos.tab`).
 * - The output profile includes all state quantities defined in `ctl`.
 * - Only diagonal uncertainties are written; correlations are not stored.
 *
 * @warning
 * - The covariance matrix `s` must be symmetric and positive-definite.
 * - The state vector mapping (`x2atm`) must correspond to the matrix ordering.
 *
 * @author Lars Hoffmann
 */
void write_stddev(
  const char *quantity,
  const ret_t * ret,
  const ctl_t * ctl,
  const atm_t * atm,
  const gsl_matrix * s);

/**
 * @brief Write emissivity lookup tables to disk.
 *
 * Writes emissivity lookup tables stored in @p tbl using the output format
 * specified by @c ctl->tblfmt:
 *  - ASCII (@c ctl->tblfmt == 1)
 *  - compact binary (@c ctl->tblfmt == 2)
 *  - netCDF (@c ctl->tblfmt == 3)
 *
 * The function iterates over all emitters (@p ig) and detectors/channels
 * (@p id) and dispatches the corresponding per-table writer.
 *
 * Lookup tables with no pressure levels (i.e. @c tbl->np[id][ig] <= 0) are
 * skipped intentionally. In this case, a warning is issued and no output
 * artifact (file or NetCDF variable) is written for the affected table.
 * This allows trace gases or channels with negligible contribution to be
 * omitted cleanly from the output.
 *
 * After writing the lookup tables, the associated spectral filter functions
 * are written to separate files (format-dependent).
 *
 * @param[in] ctl  Pointer to the control structure defining the lookup table
 *                 output format, base filenames, emitters, and spectral channels.
 * @param[in] tbl  Pointer to the lookup-table structure containing the
 *                 emissivity data and filter functions to be written.
 *
 * @warning If an unsupported lookup table format is specified via
 *          @c ctl->tblfmt, the function aborts with an error message.
 *
 * @note Missing output artifacts caused by skipped tables are handled
 *       gracefully by the corresponding readers, which treat missing
 *       tables as empty.
 *
 * @author Lars Hoffmann
 */
void write_tbl(
  const ctl_t * ctl,
  const tbl_t * tbl);

/**
 * @brief Write one emissivity lookup table in human-readable ASCII format.
 *
 * Writes the lookup table for detector/channel index @p id and emitter/gas
 * index @p ig to a file named:
 *
 *     <ctl->tblbase>_<ctl->nu[id]>_<ctl->emitter[ig]>.tab
 *
 * The ASCII file contains four columns:
 *  1. pressure [hPa]
 *  2. temperature [K]
 *  3. column density [molecules/cm^2]
 *  4. emissivity [-]
 *
 * A header describing the columns is always written. If the table is empty
 * (i.e. @c tbl->np[id][ig] == 0), the file will contain only the header.
 *
 * @param[in] ctl  Control structure providing filename base, frequencies, and
 *                 emitter names.
 * @param[in] tbl  Lookup table data to be written.
 * @param[in] id   Detector/channel index.
 * @param[in] ig   Emitter/gas index.
 *
 * @author Lars Hoffmann
 */
void write_tbl_asc(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const int id,
  const int ig);

/**
 * @brief Write one emissivity lookup table in compact binary format.
 *
 * Writes the lookup table for detector/channel index @p id and emitter/gas
 * index @p ig to a file named:
 *
 *     <ctl->tblbase>_<ctl->nu[id]>_<ctl->emitter[ig]>.bin
 *
 * The file contains a length field followed by a packed byte blob produced by
 * tbl_pack():
 *  - @c size_t nbytes   : number of bytes in the packed blob
 *  - @c uint8_t blob[]  : packed lookup table data (see tbl_pack())
 *
 * If the table is empty (@c tbl->np[id][ig] == 0), the packed blob still
 * encodes this (it includes @c np = 0), so the file remains a valid
 * representation of an empty table.
 *
 * @param[in] ctl  Control structure providing filename base, frequencies, and
 *                 emitter names.
 * @param[in] tbl  Lookup table data to be packed and written.
 * @param[in] id   Detector/channel index.
 * @param[in] ig   Emitter/gas index.
 *
 * @see tbl_pack()
 *
 * @author Lars Hoffmann
 */
void write_tbl_bin(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const int id,
  const int ig);

/**
 * @brief Write one packed lookup table to a NetCDF file.
 *
 * Writes the lookup table for detector/channel index @p id and emitter/gas
 * index @p ig to the NetCDF file:
 *
 *     <ctl->tblbase>_<ctl->emitter[ig]>.nc
 *
 * The table is stored as a packed byte blob (produced by tbl_pack()) in a
 * variable named @c tbl_XXXX with a corresponding dimension @c len_XXXX,
 * where @c XXXX is the formatted frequency @c ctl->nu[id]. The variable type
 * is @c NC_UBYTE.
 *
 * If the NetCDF file does not exist, it is created and initialized with the
 * global attributes:
 *  - @c format_version = "1"
 *  - @c emitter = ctl->emitter[ig]
 *
 * To prevent accidental overwrites, the function terminates with an error if
 * the target variable already exists in the file.
 *
 * @param[in] ctl  Control structure defining emitter names, frequencies, and
 *                 base filename.
 * @param[in] tbl  Lookup table data structure containing the lookup table to
 *                 be packed and written.
 * @param[in] id   Detector/channel index.
 * @param[in] ig   Emitter/gas index.
 *
 * @see tbl_pack()
 *
 * @author Lars Hoffmann
 */
void write_tbl_nc(
  const ctl_t * ctl,
  const tbl_t * tbl,
  const int id,
  const int ig);

/**
 * @brief Map retrieval state vector back to atmospheric structure.
 *
 * Updates the atmospheric data structure (`atm_t`) from the contents of
 * a retrieval state vector (`x`). This function performs the inverse
 * transformation of `atm2x()`, assigning retrieved quantities such as
 * pressure, temperature, gas volume mixing ratios, extinction, and
 * cloud/surface parameters to the corresponding atmospheric fields.
 *
 * @param[in]  ctl  Pointer to control structure defining retrieval settings,
 *                  vertical range limits, and active retrieval flags.
 * @param[in]  x    Pointer to retrieval state vector containing the updated values.
 * @param[out] atm  Pointer to atmospheric data structure to be updated.
 *
 * @details
 * - Each atmospheric quantity is updated only within its respective retrieval
 *   altitude range (`ctl->ret*_zmin`/`ctl->ret*_zmax`).
 * - For each retrievable parameter, the helper routine `x2atm_help()` is called
 *   to sequentially read the next element from the state vector.
 * - The order of assignments must match that in `atm2x()` to ensure
 *   one-to-one correspondence between state vector indices and atmospheric fields.
 *
 * Quantities mapped:
 * - Pressure (`p[zmin:zmax]`)
 * - Temperature (`t[zmin:zmax]`)
 * - Gas volume mixing ratios (`q[ig][zmin:zmax]`)
 * - Extinction coefficients (`k[iw][zmin:zmax]`)
 * - Cloud parameters (`clz`, `cldz`, `clk`)
 * - Surface parameters (`sft`, `sfeps`)
 *
 * @see atm2x, x2atm_help, ctl_t, atm_t
 *
 * @note
 * - Only parameters marked as retrievable in `ctl` (e.g. `ret_sft`, `ret_clk`)
 *   are modified.
 * - The state vector index (`n`) advances automatically as each element is read.
 * - The helper function `x2atm_help()` abstracts sequential access to vector elements.
 *
 * @warning
 * - The atmospheric profile must be initialized before calling this function.
 * - Retrieval ranges and flags in `ctl` must correspond exactly to those used
 *   in the forward-model configuration.
 * - Mismatch between `atm2x()` and `x2atm()` ordering will cause incorrect mappings.
 *
 * @author Lars Hoffmann
 */
void x2atm(
  const ctl_t * ctl,
  const gsl_vector * x,
  atm_t * atm);

/**
 * @brief Helper function to extract a single value from the retrieval state vector.
 *
 * Retrieves the next element from the state vector `x` and assigns it
 * to the provided scalar variable. This function is used by `x2atm()`
 * to sequentially map the contents of the retrieval vector into the
 * corresponding fields of the atmospheric structure.
 *
 * @param[out] value  Pointer to the scalar variable to be updated.
 * @param[in]  x      Pointer to the retrieval state vector (`gsl_vector`).
 * @param[in,out] n   Pointer to the current index in the state vector.
 *                    The index is incremented after each extraction.
 *
 * @details
 * - Acts as a lightweight iterator over the state vector elements.
 * - Ensures consistent and sequential assignment order between `atm2x()` and `x2atm()`.
 * - Increments the index counter `*n` after reading a value, maintaining
 *   the correct position in the vector for subsequent calls.
 *
 * @see x2atm, atm2x
 *
 * @note
 * - The function performs no range checking; it assumes that the index
 *   is within valid bounds of the state vector length.
 * - Typically used only internally by retrieval mapping routines.
 *
 * @author Lars Hoffmann
 */
void x2atm_help(
  double *value,
  const gsl_vector * x,
  size_t *n);

/**
 * @brief Copy elements from the measurement vector @p y into the observation structure.
 *
 * Decomposes the 1-D measurement vector @p y into its radiance components and
 * writes them into the 2-D observation array @c obs->rad[id][ir], using the
 * same ordering as produced by the corresponding forward model.
 *
 * Only entries for which @c obs->rad[id][ir] is finite are updated.  This allows
 * missing or masked radiances to remain untouched in the observation structure.
 *
 * @param[in]  ctl  Control settings defining the number of detector channels
 *                  (@c ctl->nd) and other retrieval configuration parameters.
 * @param[in]  y    Measurement vector containing radiances in forward-model order.
 * @param[out] obs  Observation structure whose radiance array (@c obs->rad)
 *                  is to be filled with values from @p y.
 *
 * @details
 * The function loops over all ray paths (`obs->nr`) and detector channels
 * (`ctl->nd`).  For each pair (detector @c id, ray @c ir) where an existing value
 * in @c obs->rad[id][ir] is finite, the next element from the measurement vector
 * @p y is inserted.  The counter @c m tracks progression through @p y.
 *
 * This function is the inverse operation of the packing performed when
 * constructing the measurement vector from an @ref obs_t structure.
 *
 * @note The measurement vector @p y must contain as many finite elements as the
 *       number of finite entries in @c obs->rad, in the same scanning order.
 *
 * @see obs_t, ctl_t
 *
 * @author Lars Hoffmann
 */
void y2obs(
  const ctl_t * ctl,
  const gsl_vector * y,
  obs_t * obs);

#endif
