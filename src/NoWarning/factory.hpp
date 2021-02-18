// *****************************************************************************
/*!
  \file      src/NoWarning/factory.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include boost/functional/factory.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_factory_h
#define nowarning_factory_h

#include "Macro.hpp"

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
