// *****************************************************************************
/*!
  \file      src/NoWarning/pugixml.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include pugixml.hpp with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_pugixml_h
#define nowarning_pugixml_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
#endif

#include <pugixml.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_pugixml_h
