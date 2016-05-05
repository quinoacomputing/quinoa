//******************************************************************************
/*!
  \file      src/NoWarning/mkl_lapacke.h
  \author    J. Bakosi
  \date      Tue 03 May 2016 08:41:37 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include mkl_lapacke.h with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_mkl_lapacke_h
#define nowarning_mkl_lapacke_h

#ifdef HAS_MKL

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <mkl_lapacke.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif  // HAS_MKL

#endif // nowarning_mkl_lapacke_h
