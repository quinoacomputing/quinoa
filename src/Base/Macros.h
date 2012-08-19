// -----------------------------------------------------------------------------
// \file    src/Base/Macros.h
// \author  jbakosi
// \date    The Aug 14 9:32:00 2012
// \brief   Compile-time macros
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------

#include <cstdlib>
#include <sys/time.h>
#include <ctime>

#ifdef __INTEL_COMPILER
// disregard Intel compiler's remark #981:
// operands are evaluated in unspecified order
#pragma warning(disable:981)
#endif

// code modifying macros
#define SAVERESTART     // save restartfile
#define MICROMIXING     // compute micromixing, otherwise only velocity field
#define VCIEM		// for micromixing use VCIEM, otherwise IEM
//#define PROJECTION    // use Fox's projection method for VCIEM
//#define WALLFUNCTIONS // use wall-functions instead of elliptic relaxation
//#define ME		// use 2-stage modified Euler, otherwise: Euler-Maruyama

// MKL's random number generator specifics
#define SEED            1
#define BRNG_TABLE      VSL_BRNG_MCG59  // period ~ 1.4e+17
#define BRNG_FEW        VSL_BRNG_MCG31  // period ~ 2.1e+9
#define UNIFORM_METHOD  VSL_METHOD_DUNIFORM_STD
#define GAUSSIAN_METHOD VSL_METHOD_DGAUSSIAN_BOXMULLER
#define GAMMA_METHOD    VSL_METHOD_DGAMMA_GNORM

// labels for preconditioners
#define PC_NONE         0       // no preconditioner
#define PC_JACOBI       1       // Jacobi preconditioner
#define PC_SGS          2       // symmetric Gauss-Seidel preconditioner

// labels for outputing results
#define INST            0       // label for outputing instantaneous fields
#define TAV             1       // label for outputing time-averaged fields

// macros for fine-grained profiling
#define STARTTIME \
struct timeval START_TIME, END_TIME; \
int total_usecs; \
gettimeofday(&START_TIME, (struct timezone*)0);

#define ENDTIME \
gettimeofday(&END_TIME, (struct timezone*)0); \
total_usecs = (END_TIME.tv_sec-START_TIME.tv_sec) * 1000000 + (END_TIME.tv_usec-START_TIME.tv_usec); \
printf("Total time was %d uSec.\n", total_usecs);

// general purpose function-macros
#define SWAP(a,b,tmp) {tmp=a; a=b; b=tmp;}
#define ERR(s) {printf("Error: %s\n",s); exit(-1);}

// Duff's device macros for loop unrolling (not used currently)
#define XXDUFF(N,CMD) XXDUFF16(N,CMD)

#define XXDUFF4(N,CMD) {                        \
    register int __wn = ((N) + 3) / 4;          \
    if ((N) > 0) {                              \
      switch (N % 4)                            \
        {                                       \
              case 0:      do {  CMD            \
              case 3:            CMD            \
              case 2:            CMD            \
              case 1:            CMD            \
              } while (--__wn > 0);             \
        }                                       \
    }                                           \
  }

#define XXDUFF8(N,CMD) {                        \
    register int __wn = ((N) + 7) / 8;          \
    if ((N) > 0) {                              \
      switch (N % 8)                            \
        {                                       \
              case 0:      do {  CMD            \
              case 7:            CMD            \
              case 6:            CMD            \
              case 5:            CMD            \
              case 4:            CMD            \
              case 3:            CMD            \
              case 2:            CMD            \
              case 1:            CMD            \
              } while (--__wn > 0);             \
        }                                       \
    }                                           \
  }

#define XXDUFF16(N,CMD) {                       \
    register int __wn = ((N) + 15) / 16;        \
    if ((N) > 0) {                              \
      switch (N % 16)                           \
        {                                       \
              case 0:      do {  CMD            \
              case 15:           CMD            \
              case 14:           CMD            \
              case 13:           CMD            \
              case 12:           CMD            \
              case 11:           CMD            \
              case 10:           CMD            \
              case 9:            CMD            \
              case 8:            CMD            \
              case 7:            CMD            \
              case 6:            CMD            \
              case 5:            CMD            \
              case 4:            CMD            \
              case 3:            CMD            \
              case 2:            CMD            \
              case 1:            CMD            \
              } while (--__wn > 0);             \
        }                                       \
    }                                           \
  }
