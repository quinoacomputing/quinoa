// *****************************************************************************
/*!
  \file      src/NoWarning/TString.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include Root's TString.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_TString_h
#define nowarning_TString_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include <TString.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TString_h
