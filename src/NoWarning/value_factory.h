// *****************************************************************************
/*!
  \file      src/NoWarning/value_factory.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include boost/functional/value_factory.hpp with turning off
             specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_value_factory_h
#define nowarning_value_factory_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wsign-conversion"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/functional/value_factory.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_value_factory_h
