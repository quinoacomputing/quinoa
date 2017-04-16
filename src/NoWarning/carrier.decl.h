// *****************************************************************************
/*!
  \file      src/NoWarning/carrier.decl.h
  \author    J. Bakosi
  \date      Tue 16 Aug 2016 09:21:13 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include carrier.decl.h with turning off specific compiler
             warnings.
*/
// *****************************************************************************
#ifndef nowarning_carrier_decl_h
#define nowarning_carrier_decl_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
#elif STRICT_GNUC
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include "../Inciter/carrier.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif STRICT_GNUC
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_carrier_decl_h
