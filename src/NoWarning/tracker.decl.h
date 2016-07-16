// *****************************************************************************
/*!
  \file      src/NoWarning/tracker.decl.h
  \author    F.J. Gonzalez
  \date      Wed 04 May 2016 09:41:18 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include tracker.decl.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_tracker_decl_h
#define nowarning_tracker_decl_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wunused-parameter"
  #pragma clang diagnostic ignored "-Wshorten-64-to-32"
  #pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include "../Inciter/tracker.decl.h"

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_tracker_decl_h
