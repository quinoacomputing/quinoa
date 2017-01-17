// *****************************************************************************
/*!
  \file      src/NoWarning/philox.h
  \author    J. Bakosi
  \date      Thu 12 Jan 2017 10:53:05 AM MST
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include Random123/philox.h with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_philox_h
#define nowarning_philox_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi"
#endif

#include <Random123/philox.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_philox_h
