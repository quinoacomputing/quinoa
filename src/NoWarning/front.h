// *****************************************************************************
/*!
  \file      src/NoWarning/front.h
  \copyright 2016-2018, Triad National Security, LLC.
  \brief     Include brigand/sequences/front.hpp with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_front_h
#define nowarning_front_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wundef"
#endif

#include <brigand/sequences/front.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_front_h
