//******************************************************************************
/*!
  \file      src/NoWarning/mpi.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 08:02:07 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include mpi.h with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_mpi_h
#define nowarning_mpi_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wcast-align"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
  #pragma GCC diagnostic ignored "-Wlong-long"
#endif

#include <mpi.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif


#endif // nowarning_mpi_h
