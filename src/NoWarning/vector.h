// *****************************************************************************
/*!
  \file      src/NoWarning/vector.h
  \author    J. Bakosi
  \date      Wed 11 May 2016 06:59:08 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/vector.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_vector_h
#define nowarning_vector_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#elif STRICT_GNUC
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/mpl/vector.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif STRICT_GNUC
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_vector_h
