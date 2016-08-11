// *****************************************************************************
/*!
  \file      src/NoWarning/performer.def.h
  \author    J. Bakosi
  \date      Wed 03 Aug 2016 02:06:18 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include performer.def.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_performer_def_h
#define nowarning_performer_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wunused-variable"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#include "../Inciter/performer.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_performer_def_h
