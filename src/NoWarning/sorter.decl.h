// *****************************************************************************
/*!
  \file      src/NoWarning/sorter.decl.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include sorter.decl.h with turning off specific compiler warnings.
*/
// *****************************************************************************
#ifndef nowarning_sorter_decl_h
#define nowarning_sorter_decl_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wunused-parameter"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../Inciter/sorter.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_sorter_decl_h
