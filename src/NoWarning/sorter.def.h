// *****************************************************************************
/*!
  \file      src/NoWarning/sorter.def.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include sorter.def.h with turning off specific compiler warnings.
*/
// *****************************************************************************
#ifndef nowarning_sorter_def_h
#define nowarning_sorter_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-variable"
  #pragma clang diagnostic ignored "-Wcast-qual"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-variable"
  #pragma GCC diagnostic ignored "-Wsuggest-attribute=noreturn"
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "../Inciter/sorter.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_sorter_def_h
