// *****************************************************************************
/*!
  \file      src/NoWarning/replace.h
  \author    J. Bakosi
  \date      Wed 11 May 2016 07:02:13 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/algorithm/string/replace.hpp with turning off
             specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_replace_h
#define nowarning_replace_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
  #pragma clang diagnostic ignored "-Wdocumentation"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

#include <boost/algorithm/string/replace.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_replace_h
