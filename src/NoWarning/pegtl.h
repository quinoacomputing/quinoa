// *****************************************************************************
/*!
  \file      src/NoWarning/pegtl.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 02:38:11 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pegtl/pegtl.hh with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_pegtl_h
#define nowarning_pegtl_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <pegtl/pegtl.hh>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_pegtl_h
