// *****************************************************************************
/*!
  \file      src/NoWarning/set.h
  \author    J. Bakosi
  \date      Mon 09 May 2016 03:38:38 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/set.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_set_h
#define nowarning_set_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

#include <boost/mpl/set.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_set_h
