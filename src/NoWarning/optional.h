// *****************************************************************************
/*!
  \file      src/NoWarning/optional.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 02:02:50 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/optional.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_optional_h
#define nowarning_optional_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wundef"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/optional.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_optional_h
