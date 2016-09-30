// *****************************************************************************
/*!
  \file      src/NoWarning/pegtl.h
  \author    J. Bakosi
  \date      Fri 30 Sep 2016 12:41:48 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pegtl/pegtl.hh with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_pegtl_h
#define nowarning_pegtl_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <pegtl/pegtl.hh>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_pegtl_h
