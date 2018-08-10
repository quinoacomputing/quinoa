// *****************************************************************************
/*!
  \file      src/NoWarning/charestatecollector.def.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include charestatecollector.def.h with turning off specific
             compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_charestatecollector_def_h
#define nowarning_charestatecollector_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
  #pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include "../Base/charestatecollector.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_charestatecollector_def_h
