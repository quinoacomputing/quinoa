// *****************************************************************************
/*!
  \file      src/NoWarning/backward.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Include backward.hpp with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_backward_h
#define nowarning_backward_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include "backward.hpp"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_backward_h
