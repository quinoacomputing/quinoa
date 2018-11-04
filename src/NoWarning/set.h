// *****************************************************************************
/*!
  \file      src/NoWarning/set.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include brigand/sequences/set.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_set_h
#define nowarning_set_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wundef"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wpedantic"
#endif

#undef I
#include <brigand/sequences/set.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_set_h
