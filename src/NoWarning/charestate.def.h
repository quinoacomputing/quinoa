// *****************************************************************************
/*!
  \file      src/NoWarning/charestate.def.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include charestate.def.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_charestate_def_h
#define nowarning_charestate_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic ignored "-Wunused-variable"
// #elif defined(STRICT_GNUC)
//   #pragma GCC diagnostic push
//   #pragma GCC diagnostic ignored "-Wcast-qual"
//   #pragma GCC diagnostic ignored "-Wshadow"
//   #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//   #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../Base/charestate.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
// #elif defined(STRICT_GNUC)
//   #pragma GCC diagnostic pop
#endif

#endif // nowarning_charestate_def_h
