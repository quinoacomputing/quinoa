// *****************************************************************************
/*!
  \file      src/NoWarning/particlewriter.def.h
  \author    F.J. Gonzalez
  \date      Wed 04 May 2016 09:41:18 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include particlewriter.def.h with turning off specific compiler
             warnings.
*/
// *****************************************************************************
#ifndef nowarning_particlewriter_def_h
#define nowarning_particlewriter_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include "../IO/particlewriter.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_particlewriter_def_h
