//******************************************************************************
/*!
  \file      src/NoWarning/tut.h
  \author    J. Bakosi
  \date      Thu 05 May 2016 08:54:37 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include tut/tut.hpp with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_tut_h
#define nowarning_tut_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
  #pragma clang diagnostic ignored "-Wfloat-equal"
  #pragma clang diagnostic ignored "-Wmissing-noreturn"
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wunused-function"
  #pragma clang diagnostic ignored "-Wunused-variable"
  #pragma clang diagnostic ignored "-Wsign-compare"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
  #pragma GCC diagnostic ignored "-Wfloat-equal"
  #pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
  #pragma GCC diagnostic ignored "-Wswitch-default"
#elif defined(__INTEL_COMPILER)
  #pragma warning( push )
  #pragma warning( disable: 3346 )
#endif

#include <tut/tut.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#elif defined(__INTEL_COMPILER)
  #pragma warning( pop )
#endif

#endif // nowarning_tut_h
