// *****************************************************************************
/*!
  \file      src/NoWarning/controller.decl.h
  \copyright 2020 Charmworks, Inc.
             All rights reserved. See the LICENSE file for details.
  \brief     Include controller.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_controller_decl_h
#define nowarning_controller_decl_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wunused-private-field"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wcast-qual"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wheader-hygiene"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wdouble-promotion"
  #pragma clang diagnostic ignored "-Wnon-virtual-dtor"
  #pragma clang diagnostic ignored "-Wshadow"
  #pragma clang diagnostic ignored "-Wshadow-field"
  #pragma clang diagnostic ignored "-Wshadow-field-in-constructor"
  #pragma clang diagnostic ignored "-Wswitch-enum"
  #pragma clang diagnostic ignored "-Wcovered-switch-default"
  #pragma clang diagnostic ignored "-Wzero-length-array"
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wundefined-func-template"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
  #pragma GCC diagnostic ignored "-Wredundant-decls"
  #pragma GCC diagnostic ignored "-Wfloat-equal"
  #pragma GCC diagnostic ignored "-Wextra"
  #pragma GCC diagnostic ignored "-Wdeprecated-copy"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 181 )
  #pragma warning( disable: 1720 )
  #pragma warning( disable: 2282 )
#endif

#include "../Transfer/controller.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_controller_decl_h
