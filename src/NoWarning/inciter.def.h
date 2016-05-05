//******************************************************************************
/*!
  \file      src/NoWarning/inciter.def.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:46:33 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include inciter.def.h with turning off specific compiler
             warnings
*/
//******************************************************************************
#ifndef nowarning_inciter_def_h
#define nowarning_inciter_def_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wextra-semi"
  #pragma clang diagnostic ignored "-Wmissing-prototypes"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

#include "../Main/inciter.def.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_inciter_def_h
