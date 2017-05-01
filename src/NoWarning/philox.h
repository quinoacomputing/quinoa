// *****************************************************************************
/*!
  \file      src/NoWarning/philox.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include Random123/philox.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_philox_h
#define nowarning_philox_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
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

#include <Random123/philox.h>

#ifdef POWERPC
  #define __powerpc__
  #undef __x86_64__
#endif

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_philox_h
