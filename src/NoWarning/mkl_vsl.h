// *****************************************************************************
/*!
  \file      src/NoWarning/mkl_vsl.h
  \author    J. Bakosi
  \date      Sat 15 Apr 2017 11:40:53 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include mkl_vsl.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_mkl_vsl_h
#define nowarning_mkl_vsl_h

#include "Macro.h"

#ifdef HAS_MKL

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <mkl_vsl.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif  // HAS_MKL

#endif // nowarning_mkl_vsl_h
