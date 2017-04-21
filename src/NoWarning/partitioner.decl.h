// *****************************************************************************
/*!
  \file      src/NoWarning/partitioner.decl.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include partitioner.decl.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_partitioner_decl_h
#define nowarning_partitioner_decl_h

#include "Macro.h"

#if defined(__clang__)
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "../Inciter/partitioner.decl.h"

#if defined(__clang__)
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_partitioner_decl_h
