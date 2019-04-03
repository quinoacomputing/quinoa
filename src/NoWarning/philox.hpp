// *****************************************************************************
/*!
  \file      src/NoWarning/philox.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include Random123/philox.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_philox_h
#define nowarning_philox_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wexpansion-to-defined"
#endif

#ifdef __powerpc__
  #define POWERPC
  #undef __powerpc__
  #define __x86_64__
  #define R123_USE_MULHILO64_MULHI_INTRIN 0
  #define R123_USE_GNU_UINT128 1
#endif

#if defined(__PGI) || defined(_CRAYC)
  #undef R123_USE_GNU_UINT128
  #undef R123_USE_MULHILO64_C99
  #define R123_USE_MULHILO64_C99 1
#endif

#include <Random123/philox.h>

#ifdef POWERPC
  #define __powerpc__
  #undef __x86_64__
#endif

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_philox_h
