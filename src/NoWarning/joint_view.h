// *****************************************************************************
/*!
  \file      src/NoWarning/joint_view.h
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     Include boost/mpl/joint_view.hpp with turning off specific compiler
             warnings
*/
// *****************************************************************************
#ifndef nowarning_joint_view_h
#define nowarning_joint_view_h

#include "Macro.h"

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wdisabled-macro-expansion"
#endif

#include <boost/mpl/joint_view.hpp>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_joint_view_h
