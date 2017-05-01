// *****************************************************************************
/*!
  \file      src/NoWarning/linearmap.def.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include linearmap.def.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_linearmap_def_h
#define nowarning_linearmap_def_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include "../LoadBalance/linearmap.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_linearmap_def_h
