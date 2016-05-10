// *****************************************************************************
/*!
  \file      src/NoWarning/vector.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 10:04:38 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/vector.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_vector_h
#define nowarning_vector_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <boost/mpl/vector.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_vector_h
