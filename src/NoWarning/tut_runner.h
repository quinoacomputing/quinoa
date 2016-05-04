//******************************************************************************
/*!
  \file      src/NoWarning/tut_runner.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 08:55:18 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include tut/tut_runner.hpp with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_tut_runner_h
#define nowarning_tut_runner_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdeprecated"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#endif

#include <tut/tut_runner.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_tut_runner_h
