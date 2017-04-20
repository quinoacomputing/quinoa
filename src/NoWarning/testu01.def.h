// *****************************************************************************
/*!
  \file      src/NoWarning/testu01.def.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include testu01.def.h with turning off specific compiler warnings
*/
// *****************************************************************************

// Note the absence of include guards. This is because the Charm++ definition
// file included below may need to be included multiple times in a single
// translation unit: one guarded by CK_TEMPLATES_ONLY and one without, where

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wreorder"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wreorder"
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../RNGTest/testu01.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif
