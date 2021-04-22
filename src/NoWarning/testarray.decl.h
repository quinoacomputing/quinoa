// *****************************************************************************
/*!
  \file      src/NoWarning/testarray.decl.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include testarray.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_testarray_decl_h
#define nowarning_testarray_decl_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wundefined-func-template"
  #pragma clang diagnostic ignored "-Wdocumentation"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wsign-compare"
  #pragma clang diagnostic ignored "-Wcast-qual"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wheader-hygiene"
  #pragma clang diagnostic ignored "-Wdouble-promotion"
  #pragma clang diagnostic ignored "-Wunused-private-field"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wshadow"
  #pragma clang diagnostic ignored "-Wshadow-field"
  #pragma clang diagnostic ignored "-Wshadow-field-in-constructor"
  #pragma clang diagnostic ignored "-Wnon-virtual-dtor"
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wcovered-switch-default"
  #pragma clang diagnostic ignored "-Woverloaded-virtual"
  #pragma clang diagnostic ignored "-Wswitch-enum"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../UnitTest/testarray.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_testarray_decl_h
