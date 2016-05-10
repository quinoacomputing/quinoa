//******************************************************************************
/*!
  \file      src/NoWarning/value_factory.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 04:15:53 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/functional/value_factory.hpp with turning off
             specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_value_factory_h
#define nowarning_value_factory_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/functional/value_factory.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_value_factory_h
