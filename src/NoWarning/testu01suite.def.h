//******************************************************************************
/*!
  \file      src/NoWarning/testu01suite.def.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:40:54 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include testu01suite.def.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_testu01suite_def_h
#define nowarning_testu01suite_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "../RNGTest/testu01suite.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_testu01suite_def_h
