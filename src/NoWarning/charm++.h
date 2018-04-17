// *****************************************************************************
/*!
  \file      src/NoWarning/charm++.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include charm++.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_charmpp_h
#define nowarning_charmpp_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wsign-compare"
  #pragma clang diagnostic ignored "-Wheader-hygiene"
  #pragma clang diagnostic ignored "-Wunused-variable"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wdouble-promotion"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wnon-virtual-dtor"
  #pragma clang diagnostic ignored "-Wzero-length-array"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wshadow"
  #pragma clang diagnostic ignored "-Wcovered-switch-default"
  #pragma clang diagnostic ignored "-Wswitch-enum"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wmismatched-tags"
  #pragma clang diagnostic ignored "-Wunused-private-field"
  #pragma clang diagnostic ignored "-Wundefined-func-template"
  #pragma clang diagnostic ignored "-Wcomma"
  #pragma clang diagnostic ignored "-Woverloaded-virtual"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wcast-qual"
  #pragma clang diagnostic ignored "-Wshadow-field"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wredundant-decls"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wfloat-equal"
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wpedantic"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wswitch-default"
  #pragma GCC diagnostic ignored "-Wextra"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 181 )
  #pragma warning( disable: 1720 )
  #pragma warning( disable: 1125 )
#endif

#include <charm++.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_charmpp_h
