// *****************************************************************************
/*!
  \file      src/NoWarning/rngtest.def.h
  \author    J. Bakosi
  \date      Tue 03 May 2016 06:58:06 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include rngtest.def.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_rngtest_def_h
#define nowarning_rngtest_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../Main/rngtest.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_rngtest_def_h
