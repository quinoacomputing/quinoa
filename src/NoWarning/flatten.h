// *****************************************************************************
/*!
  \file      src/NoWarning/flatten.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include brigand/algorithms/flatten.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_flatten_h
#define nowarning_flatten_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include <brigand/algorithms/flatten.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_flatten_h
