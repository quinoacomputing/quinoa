//******************************************************************************
/*!
  \file      src/NoWarning/factory.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 09:50:49 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/functional/factory.hpp with turning off specific
             compiler warnings
*/
//******************************************************************************
#ifndef nowarning_factory_h
#define nowarning_factory_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/functional/factory.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_factory_h
