//******************************************************************************
/*!
  \file      src/NoWarning/pugixml.h
  \author    J. Bakosi
  \date      Thu 05 May 2016 12:49:38 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pugixml.hpp with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_pugixml_h
#define nowarning_pugixml_h

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Werror=effc++"
#elif defined(__INTEL_COMPILER)
#endif

#include <pugixml.hpp>

#if defined(__clang__)
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
#endif

#endif // nowarning_pugixml_h
