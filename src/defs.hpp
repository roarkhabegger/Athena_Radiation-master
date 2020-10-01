#ifndef DEFS_HPP_
#define DEFS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file defs.hpp.in
//  \brief Template file for defs.hpp.  When the configure.py script is run, a new
//  defs.hpp file will be created (overwriting the last) from this template.  This new
//  file contains Athena++ specific cpp macros and definitions set by configure.

//----------------------------------------------------------------------------------------
// macros which define physics and algorithms

// configure.py dict(definitions) string values:
// problem generator
#define PROBLEM_GENERATOR "CR_Parker"

// coordinate system
#define COORDINATE_SYSTEM "cartesian"

// Riemann solver
#define RIEMANN_SOLVER "hlld"

// configure.py dict(definitions) Boolean values:
// enable shearing box? default=0 (false)
#define SHEARING_BOX 0

// Equation of state
#define EQUATION_OF_STATE "adiabatic"

// use general EOS framework default=0 (false).
#define GENERAL_EOS 0

// use EOS table default=0 (false).
#define EOS_TABLE_ENABLED 0

// non-barotropic equation of state (i.e. P not simply a func of rho)? default=1 (true)
#define NON_BAROTROPIC_EOS 1

// include magnetic fields? default=0 (false)
#define MAGNETIC_FIELDS_ENABLED 1

// include super-time-stepping? default=0 (false)
#define STS_ENABLED 0

// include self gravity? default=0 (false)
#define SELF_GRAVITY_ENABLED 0

// include radiative transfer? default=0 (false)
#define RADIATION_ENABLED 0

// include cosmic ray transport? default=0 (false)
#define CR_ENABLED 1

// include thermal conduction? default=0 (false)
#define TC_ENABLED 1

// enable special or general relativity? default=0 (false)
#define RELATIVISTIC_DYNAMICS 0

// enable general relativity? default=0 (false)
#define GENERAL_RELATIVITY 0

// enable GR frame transformations? default=0 (false)
#define FRAME_TRANSFORMATIONS 0

// use single precision floating-point values (binary32)? default=0 (false; use binary64)
#define SINGLE_PRECISION_ENABLED 0

// use double precision for HDF5 output? default=0 (false; write out binary32)
#define H5_DOUBLE_PRECISION_ENABLED 0


// configure.py dict(definitions) Boolean string macros:
// (these options have the latter (false) option as defaults, unless noted otherwise)
// make use of FFT? (FFT or NO_FFT)
#define NO_FFT

// MPI parallelization (MPI_PARALLEL or NOT_MPI_PARALLEL)
#define NOT_MPI_PARALLEL

// OpenMP parallelization (OPENMP_PARALLEL or NOT_OPENMP_PARALLEL)
#define NOT_OPENMP_PARALLEL

// HDF5 output (HDF5OUTPUT or NO_HDF5OUTPUT)
#define NO_HDF5OUTPUT

// debug build macros (DEBUG or NOT_DEBUG)
#define NOT_DEBUG

#define ALI_LEN 32
// try/throw/catch C++ exception handling (ENABLE_EXCEPTIONS or DISABLE_EXCEPTIONS)
// (enabled by default)
#define ENABLE_EXCEPTIONS

// compiler options
#define COMPILED_WITH "g++"
#define COMPILER_COMMAND "g++"
#define COMPILED_WITH_OPTIONS " -O3 -std=c++11  " // NOLINT

//----------------------------------------------------------------------------------------
// macros associated with numerical algorithm (rarely modified)

#define NHYDRO 5
#define NFIELD 3
#define NWAVE 7
#define NSCALARS 0
#define NGHOST 2
#define NCR 4   // cosmic ray transport module variable
#define NTC 4   // thermal conduction variable
#define MAX_NSTAGE 5     // maximum number of stages per cycle for time-integrator
#define MAX_NREGISTER 3  // maximum number of (u, b) register pairs for time-integrator

//----------------------------------------------------------------------------------------
// general purpose macros (never modified)

// all constants specified to 17 total digits of precision = max_digits10 for "double"
#define PI 3.1415926535897932
#define TWO_PI 6.2831853071795862
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.70710678118654752
#define ONE_3RD 0.33333333333333333
#define TWO_3RD 0.66666666666666667
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+36
#define SQR(x) ( (x)*(x) )
#define SIGN(x) ( ((x) < 0.0) ? -1.0 : 1.0 )

#ifdef ENABLE_EXCEPTIONS
#define ATHENA_ERROR(x) throw std::runtime_error(x.str().c_str())
#else
#define ATHENA_ERROR(x) std::cout << x.str(); std::exit(EXIT_FAILURE)
#endif

#endif // DEFS_HPP_
