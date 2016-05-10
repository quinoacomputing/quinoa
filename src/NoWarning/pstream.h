//******************************************************************************
/*!
  \file      src/NoWarning/pstream.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:49:38 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pstreams/pstream.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_pstream_h
#define nowarning_pstream_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif

#include <pstreams/pstream.h>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_pstream_h
