// *****************************************************************************
/*!
  \file      src/NoWarning/TFile.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include <Root>/TFile.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_TFile_h
#define nowarning_TFile_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wundef"
#endif

#include <TFile.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TFile_h
