//******************************************************************************
/*!
  \file      src/NoWarning/exodusII.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 08:45:20 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Include exodusII.h with turning off specific compiler warnings
*/
//******************************************************************************
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
#elif defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Weffc++"
#endif

#include <exodusII.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_exodusII_h
