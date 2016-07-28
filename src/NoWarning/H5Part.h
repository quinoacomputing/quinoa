// *****************************************************************************
/*!
  \file      src/NoWarning/H5Part.h
  \author    J. Bakosi
  \date      Thu 28 Jul 2016 08:31:55 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include H5Part.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_H5Part_h
#define nowarning_H5Part_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#define PARALLEL_IO
#include <H5Part.h>
#undef PARALLEL_IO

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif


#endif // nowarning_H5Part_h
