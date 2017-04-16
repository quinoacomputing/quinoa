// *****************************************************************************
/*!
  \file      src/NoWarning/unittest.decl.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 09:35:31 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include unittest.decl.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_unittest_decl_h
#define nowarning_unittest_decl_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wunused-parameter"
#elif STRICT_GNUC
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wredundant-decls"
#endif

#include "../Main/unittest.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif STRICT_GNUC
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_unittest_decl_h
