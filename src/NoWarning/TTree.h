// *****************************************************************************
/*!
  \file      src/NoWarning/TTree.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include Root's TTree.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TTree_h
#define nowarning_TTree_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 1 )
  #pragma warning( disable: 181 )
  #pragma warning( disable: 522 )
  #pragma warning( disable: 2282 )
#endif

#include <TTree.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_TTree_h
