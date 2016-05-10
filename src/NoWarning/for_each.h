//******************************************************************************
/*!
  \file      src/NoWarning/for_each.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:52:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/for_each.hpp with turning off specific
             compiler warnings
*/
//******************************************************************************
#ifndef nowarning_for_each_h
#define nowarning_for_each_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <boost/mpl/for_each.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_for_each_h
