// *****************************************************************************
/*!
  \file      src/NoWarning/TString.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include <Root>/TString.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TString_h
#define nowarning_TString_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include <TString.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TString_h
