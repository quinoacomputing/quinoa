// *****************************************************************************
/*!
  \file      src/NoWarning/remove.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include brigand/algorithms/remove.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_remove_h
#define nowarning_remove_h

#include "Macro.h"

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
