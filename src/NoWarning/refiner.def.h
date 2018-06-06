// *****************************************************************************
/*!
  \file      src/NoWarning/refiner.def.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include refiner.def.h with turning off specific compiler warnings.
*/
// *****************************************************************************
#ifndef nowarning_refiner_def_h
#define nowarning_refiner_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
#endif

#include "../Inciter/refiner.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_refiner_def_h
