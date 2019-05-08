// *****************************************************************************
/*!
  \file      src/NoWarning/remove.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include brigand/algorithms/remove.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_remove_h
#define nowarning_remove_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wundef"
#endif

#include <brigand/algorithms/remove.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_remove_h
