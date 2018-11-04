// *****************************************************************************
/*!
  \file      src/NoWarning/rngtest.def.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include rngtest.def.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_rngtest_def_h
#define nowarning_rngtest_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wunused-variable"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wunused-variable"
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

#include "../Main/rngtest.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_rngtest_def_h
