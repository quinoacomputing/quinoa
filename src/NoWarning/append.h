// *****************************************************************************
/*!
  \file      src/NoWarning/append.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Include brigand/sequences/append.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_append_h
#define nowarning_append_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include <brigand/sequences/append.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_append_h
