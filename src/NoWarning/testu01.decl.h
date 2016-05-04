//******************************************************************************
/*!
  \file      src/NoWarning/testu01.decl.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:58:23 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include testu01.decl.h with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_testu01_decl_h
#define nowarning_testu01_decl_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../RNGTest/testu01.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_testu01_decl_h
