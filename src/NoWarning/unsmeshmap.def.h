// *****************************************************************************
/*!
  \file      src/NoWarning/unsmeshmap.def.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 08:44:28 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include unsmeshmap.def.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_unsmeshmap_def_h
#define nowarning_unsmeshmap_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "../LoadBalance/unsmeshmap.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_unsmeshmap_def_h
