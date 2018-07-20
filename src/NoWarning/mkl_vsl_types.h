// *****************************************************************************
/*!
  \file      src/NoWarning/mkl_vsl_types.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include mkl_vsl_types.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_mkl_vsl_types_h
#define nowarning_mkl_vsl_types_h

#include "Macro.h"

#ifdef HAS_MKL

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <mkl_vsl_types.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif  // HAS_MKL

#endif // nowarning_mkl_vsl_types_h
