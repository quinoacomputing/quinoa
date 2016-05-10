// *****************************************************************************
/*!
  \file      src/NoWarning/partitioner.decl.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:22:56 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include partitioner.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_partitioner_decl_h
#define nowarning_partitioner_decl_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "../Inciter/partitioner.decl.h"

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_partitioner_decl_h
