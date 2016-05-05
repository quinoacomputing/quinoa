//******************************************************************************
/*!
  \file      src/NoWarning/ne_nemesisI.h
  \author    J. Bakosi
  \date      Mon 02 May 2016 12:30:16 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Include ne_nemesisI.h with turning off specific compiler warnings
*/
//******************************************************************************
#ifndef nowarning_ne_nemesisI_h
#define nowarning_ne_nemesisI_h

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wreserved-id-macro"
#endif

#include <ne_nemesisI.h>

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

#endif // nowarning_ne_nemesisI_h
