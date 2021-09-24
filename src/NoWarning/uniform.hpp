// *****************************************************************************
/*!
  \file      src/NoWarning/uniform.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include Random123/uniform.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_uniform_h
#define nowarning_uniform_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wundef"
  #pragma clang diagnostic ignored "-Wimplicit-int-float-conversion"
#endif

#include <Random123/uniform.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_uniform_h
