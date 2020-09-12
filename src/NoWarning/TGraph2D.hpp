// *****************************************************************************
/*!
  \file      src/NoWarning/TGraph2D.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include Root's TGraph2D.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_TGraph2D_h
#define nowarning_TGraph2D_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include <TGraph2D.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TGraph2D_h
