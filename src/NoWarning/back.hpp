// *****************************************************************************
/*!
  \file      src/NoWarning/back.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include brigand/sequences/back.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_back_h
#define nowarning_back_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wundef"
#endif

#include <brigand/sequences/back.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_back_h
