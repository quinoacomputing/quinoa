// *****************************************************************************
/*!
  \file      src/NoWarning/charmmesh.def.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include charmmesh.def.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_charmmesh_def_h
#define nowarning_charmmesh_def_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-variable"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wcast-qual"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
//  #pragma GCC diagnostic ignored "-Wredundant-decls"
//  #pragma GCC diagnostic ignored "-Wlong-long"
//  #pragma GCC diagnostic ignored "-Wunused-parameter"
//  #pragma GCC diagnostic ignored "-Wcast-qual"
//  #pragma GCC diagnostic ignored "-Wpedantic"
//  #pragma GCC diagnostic ignored "-Wshadow"
//  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
//  #pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif

#include "../Mesh/charmmesh.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_charmmesh_def_h
