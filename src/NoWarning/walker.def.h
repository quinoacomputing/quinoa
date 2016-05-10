// *****************************************************************************
/*!
  \file      src/NoWarning/walker.def.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 11:59:26 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include walker.def.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_walker_def_h
#define nowarning_walker_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include "../Main/walker.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_walker_def_h
