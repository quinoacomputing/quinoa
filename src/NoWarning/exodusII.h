// *****************************************************************************
/*!
  \file      src/NoWarning/exodusII.h
  \author    J. Bakosi
  \date      Fri 30 Sep 2016 12:41:26 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include exodusII.h with turning off specific compiler warnings
*/
// *****************************************************************************
#ifndef nowarning_exodusII_h
#define nowarning_exodusII_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
  #pragma clang diagnostic ignored "-Wold-style-cast"
  #pragma clang diagnostic ignored "-Wsign-conversion"
  #pragma clang diagnostic ignored "-Wdeprecated"
  #pragma clang diagnostic ignored "-Wcast-align"
  #pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#endif

#include <exodusII.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_exodusII_h
