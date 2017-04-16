// *****************************************************************************
/*!
  \file      src/NoWarning/exodusII.h
  \author    J. Bakosi
  \date      Wed 16 Nov 2016 10:18:04 PM MST
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
  #pragma clang diagnostic ignored "-Wextra-semi"
#elif STRICT_GNUC
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wcast-qual"
#endif

#include <exodusII.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif STRICT_GNUC
  #pragma GCC diagnostic pop
#endif

#endif // nowarning_exodusII_h
