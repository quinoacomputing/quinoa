// *****************************************************************************
/*!
  \file      src/NoWarning/H5Block.h
  \author    J. Bakosi
  \date      Thu 28 Jul 2016 08:31:43 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include H5Block.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_H5Block_h
#define nowarning_H5Block_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#define PARALLEL_IO
#include <H5Block.h>
#undef PARALLEL_IO

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif


#endif // nowarning_H5Block_h
