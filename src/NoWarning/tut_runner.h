// *****************************************************************************
/*!
  \file      src/NoWarning/tut_runner.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Include tut/tut_runner.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_tut_runner_h
#define nowarning_tut_runner_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 1720 )
#endif

#include <tut/tut_runner.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_tut_runner_h
