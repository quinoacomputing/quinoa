// *****************************************************************************
/*!
  \file      src/NoWarning/bool.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/bool.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_bool_h
#define nowarning_bool_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/mpl/bool.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_bool_h
