// *****************************************************************************
/*!
  \file      src/NoWarning/mkl_lapacke.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Include mkl_lapacke.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_mkl_lapacke_h
#define nowarning_mkl_lapacke_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <mkl_lapacke.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_mkl_lapacke_h
