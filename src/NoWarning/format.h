// *****************************************************************************
/*!
  \file      src/NoWarning/format.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/format.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_format_h
#define nowarning_format_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wimplicit-fallthrough"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wcomma"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 367 )
#endif

#include <boost/format.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_format_h
