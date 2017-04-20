// *****************************************************************************
/*!
  \file      src/NoWarning/pegtl.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pegtl/pegtl.hh with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_pegtl_h
#define nowarning_pegtl_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshadow"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshadow-field-in-constructor"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <pegtl.hh>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_pegtl_h
