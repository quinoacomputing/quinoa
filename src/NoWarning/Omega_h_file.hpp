// *****************************************************************************
/*!
  \file      src/NoWarning/Omega_h_file.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Include Omega_h_file.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_omega_h_file_h
#define nowarning_omega_h_file_h

#include "Macro.hpp"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wextra-semi-stmt"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 3253 )
#endif

#include <Omega_h_file.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_omega_h_file_h
