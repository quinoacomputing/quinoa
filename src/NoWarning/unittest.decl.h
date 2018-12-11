// *****************************************************************************
/*!
  \file      src/NoWarning/unittest.decl.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include unittest.decl.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_unittest_decl_h
#define nowarning_unittest_decl_h

#include "Macro.h"
#include "QuinoaConfig.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wredundant-decls"
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#ifdef ENABLE_INCITER
  #include "../Main/unittestinciter.decl.h"
#else
  #include "../Main/unittest.decl.h"
#endif

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_unittest_decl_h
