//******************************************************************************
/*!
  \file      src/NoWarning/mkl_vsl.h
  \author    J. Bakosi
  \date      Mon 02 May 2016 08:02:04 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include mkl_vsl.h with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_mkl_vsl_h
#define nowarning_mkl_vsl_h

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
