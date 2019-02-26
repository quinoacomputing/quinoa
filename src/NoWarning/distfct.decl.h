// *****************************************************************************
/*!
  \file      src/NoWarning/distfct.decl.h
  \copyright 2012-2015, J. Bakosi, 2016-2019, Los Alamos National Security, LLC.
  \brief     Include distfct.decl.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_distfct_decl_h
#define nowarning_distfct_decl_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include "../Inciter/distfct.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_distfct_decl_h
