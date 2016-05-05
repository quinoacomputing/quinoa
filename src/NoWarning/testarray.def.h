//******************************************************************************
/*!
  \file      src/NoWarning/testarray.def.h
  \author    J. Bakosi
  \date      Tue 03 May 2016 08:38:03 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include testarray.def.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_testarray_def_h
#define nowarning_testarray_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include "../UnitTest/testarray.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_testarray_def_h
