// *****************************************************************************
/*!
  \file      src/NoWarning/partitioner.decl.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include partitioner.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_partitioner_decl_h
#define nowarning_partitioner_decl_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../Inciter/partitioner.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_partitioner_decl_h
