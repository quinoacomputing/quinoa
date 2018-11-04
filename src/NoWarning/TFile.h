// *****************************************************************************
/*!
  \file      src/NoWarning/TFile.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include Root's TFile.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TFile_h
#define nowarning_TFile_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 522 )
  #pragma warning( disable: 2282 )
#endif

#include <TFile.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_TFile_h
