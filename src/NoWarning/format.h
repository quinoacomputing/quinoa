//******************************************************************************
/*!
  \file      src/NoWarning/format.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:59:49 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/format.hpp with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_format_h
#define nowarning_format_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include <boost/format.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_format_h
