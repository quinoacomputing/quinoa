//******************************************************************************
/*!
  \file      src/NoWarning/replace.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 04:08:53 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/algorithm/string/replace.hpp with turning off
             specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_replace_h
#define nowarning_replace_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/algorithm/string/replace.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_replace_h
