// *****************************************************************************
/*!
  \file      src/NoWarning/threefry.h
  \author    J. Bakosi
  \date      Tue 26 Jul 2016 07:45:56 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include Random123/threefry.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_threefry_h
#define nowarning_threefry_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wconversion"
#endif

#include <Random123/threefry.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_threefry_h
