//******************************************************************************
/*!
  \file      src/NoWarning/bind.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:38:38 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/bind/bind.hpp with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_bind_h
#define nowarning_bind_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/bind/bind.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_bind_h
