//******************************************************************************
/*!
  \file      src/NoWarning/meshconv.def.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:39:49 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include meshconv.def.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_meshconv_def_h
#define nowarning_meshconv_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "../Main/meshconv.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_meshconv_def_h
