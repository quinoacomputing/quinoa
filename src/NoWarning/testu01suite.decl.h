// *****************************************************************************
/*!
  \file      src/NoWarning/testu01suite.decl.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:34:32 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include testu01suite.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_testu01suite_decl_h
#define nowarning_testu01suite_decl_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../RNGTest/testu01suite.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_testu01suite_decl_h
