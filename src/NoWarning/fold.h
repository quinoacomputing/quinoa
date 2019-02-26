// *****************************************************************************
/*!
  \file      src/NoWarning/fold.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Include brigand/algorithms/fold.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_fold_h
#define nowarning_fold_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include <brigand/algorithms/fold.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_fold_h
