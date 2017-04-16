// *****************************************************************************
/*!
  \file      src/NoWarning/tut_result.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 11:29:17 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include tut/tut_result.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_tut_result_h
#define nowarning_tut_result_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <tut/tut_result.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_tut_result_h
