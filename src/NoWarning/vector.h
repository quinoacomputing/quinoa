// *****************************************************************************
/*!
  \file      src/NoWarning/vector.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/vector.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_vector_h
#define nowarning_vector_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/mpl/vector.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_vector_h
