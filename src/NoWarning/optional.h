//******************************************************************************
/*!
  \file      src/NoWarning/optional.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:50:24 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/optional.hpp with turning off specific
             compiler warnings
*/
//******************************************************************************
#ifndef nowarning_optional_h
#define nowarning_optional_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/optional.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_optional_h
