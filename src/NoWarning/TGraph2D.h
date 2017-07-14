// *****************************************************************************
/*!
  \file      src/NoWarning/TGraph2D.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include <Root>/TGraph2D.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TGraph2D_h
#define nowarning_TGraph2D_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
#endif

#include <TGraph2D.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TGraph2D_h
