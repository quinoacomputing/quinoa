// *****************************************************************************
/*!
  \file      src/NoWarning/sip_hash.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include highwayhash/sip_hash.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_sip_hash_h
#define nowarning_sip_hash_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include <highwayhash/sip_hash.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_sip_hash_h
