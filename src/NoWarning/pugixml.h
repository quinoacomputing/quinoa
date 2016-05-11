// *****************************************************************************
/*!
  \file      src/NoWarning/pugixml.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 02:50:32 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pugixml.hpp with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_pugixml_h
#define nowarning_pugixml_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <pugixml.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_pugixml_h
