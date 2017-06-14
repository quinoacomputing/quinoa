// *****************************************************************************
/*!
  \file      src/NoWarning/TTree.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include <Root>/TTree.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TTree_h
#define nowarning_TTree_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wconversion"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <TTree.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TTree_h
