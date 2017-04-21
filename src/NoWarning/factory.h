// *****************************************************************************
/*!
  \file      src/NoWarning/factory.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/functional/factory.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_factory_h
#define nowarning_factory_h

#include "Macro.h"

#if defined(__clang__)
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/functional/factory.hpp>

#if defined(__clang__)
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_factory_h
