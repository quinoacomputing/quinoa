// *****************************************************************************
/*!
  \file      src/NoWarning/hash.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include boost/functional/hash.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_hash_h
#define nowarning_hash_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include <boost/functional/hash.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_hash_h
