// *****************************************************************************
/*!
  \file      src/NoWarning/ice_and.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 01:57:43 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/type_traits/detail/ice_and.hpp with turning off
             specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_ice_and_h
#define nowarning_ice_and_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-W#pragma-messages"
#endif

#include <boost/type_traits/detail/ice_and.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_ice_and_h
