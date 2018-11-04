// *****************************************************************************
/*!
  \file      src/NoWarning/uniform.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include Random123/uniform.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_uniform_h
#define nowarning_uniform_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundef"
#endif

#include <Random123/uniform.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_uniform_h
