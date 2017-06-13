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
#endif

#include <TFile.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_TFile_h
