// *****************************************************************************
/*!
  \file      src/NoWarning/quietcerr.def.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include quietcerr.def.h with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_quietcerr_def_h
#define nowarning_quietcerr_def_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include "../UnitTest/quietcerr.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_quietcerr_def_h
